#ifndef HISTOGRAM_IO_H
#define HISTOGRAM_IO_H

// ---------------------------------------------------------------------------
// histogram_io.h
//   Shared TFile-based I/O for the histogramming/plotting split executables
//   (clustering_hist/clustering_plot, rpt_v5_hist/rpt_v5_plot, and any future
//   util/ executable adopting the same pattern): write a named set of
//   TH1-derived histograms plus scalar run metadata into one ROOT file after
//   the event loop, then read them back by name in a separate, ntuple-free
//   plotting-stage process.
//
//   Naming convention: every object keeps its own existing ROOT name
//   (GetName()) as the TFile key -- no TDirectory nesting. This works because
//   every histogram booked in this codebase already embeds enough context
//   (score/variable/scenario/slice) in its name to be globally unique within
//   a run (see PlotObj's constructor in plotting_utilities.h, or rpt_v5's
//   Scenario histogram names like "HS_zonly_lo").
//
//   Scalar metadata (event counts, accumulator doubles, the energy label)
//   uses TParameter<Long64_t>/TParameter<Double_t>/TObjString under a
//   "meta_"-prefixed key, kept in the same file alongside the histograms.
// ---------------------------------------------------------------------------

#include <TFile.h>
#include <TH1.h>
#include <TObjString.h>
#include <TParameter.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // HistWriter
  //   RAII wrapper around a RECREATE-mode TFile. Call WriteHist/WriteHists for
  //   every histogram produced by the event loop, WriteScalar for accumulator
  //   scalars, and WriteRunMeta once at the end for provenance. Close() (or
  //   the destructor) finalizes the file.
  // ---------------------------------------------------------------------------
  class HistWriter {
  public:
    explicit HistWriter(const std::string& path)
      : file_(TFile::Open(path.c_str(), "RECREATE")) {
      if (!file_ || file_->IsZombie())
        throw std::runtime_error("HistWriter: failed to open '" + path + "' for writing");
    }

    ~HistWriter() { Close(); }

    HistWriter(const HistWriter&) = delete;
    HistWriter& operator=(const HistWriter&) = delete;

    void Close() {
      if (file_) {
        file_->Close();
        delete file_;
        file_ = nullptr;
      }
    }

    // Writes any TObject (typically a TH1D/TH2D/TProfile2D) under its own
    // GetName() into the open file.
    void WriteHist(const TObject* obj) {
      if (!file_) throw std::runtime_error("HistWriter: file not open");
      if (!obj)   throw std::runtime_error("HistWriter::WriteHist: null object");
      file_->cd();
      obj->Write(obj->GetName(), TObject::kOverwrite);
    }

    // Convenience: write every element of a container of raw pointers or
    // unique_ptrs (anything with -> yielding a TObject*) in one call.
    template <typename Container>
    void WriteHists(const Container& objs) {
      for (const auto& obj : objs) WriteHist(&*obj);
    }

    void WriteScalar(const std::string& name, Long64_t value) {
      if (!file_) throw std::runtime_error("HistWriter: file not open");
      file_->cd();
      TParameter<Long64_t> p(name.c_str(), value);
      p.Write(name.c_str(), TObject::kOverwrite);
    }

    void WriteScalar(const std::string& name, Double_t value) {
      if (!file_) throw std::runtime_error("HistWriter: file not open");
      file_->cd();
      TParameter<Double_t> p(name.c_str(), value);
      p.Write(name.c_str(), TObject::kOverwrite);
    }

    // Run metadata written once, after all histograms/scalars are in place.
    //   The VBS selection is REQUIRED rather than optional, so it cannot be
    //   forgotten. It is provenance the filename can no longer carry: the
    //   SELECTION_TAG in the path only records DEVIATIONS from the current
    //   defaults, so when those defaults changed (2026-08-26, |Deta| >= 3 /
    //   m_jj >= 0 -> no |Deta| / m_jj >= 200) every previously-written
    //   "untagged" file kept a name that now denotes a different selection.
    //   Recording the actual values inside the file makes that detectable --
    //   see HistReader::CheckSelection, which the plot stages call.
    void WriteRunMeta(const std::string& energyLabel, Long64_t nEventsTotal,
                      double vbsDeta, double vbsMjj) {
      if (!file_) throw std::runtime_error("HistWriter: file not open");
      file_->cd();
      TObjString label(energyLabel.c_str());
      label.Write("meta_energy_label", TObject::kOverwrite);
      WriteScalar("meta_n_events_total", nEventsTotal);
      WriteScalar("meta_vbs_deta", vbsDeta);
      WriteScalar("meta_vbs_mjj",  vbsMjj);
    }

  private:
    TFile* file_;
  };

  // ---------------------------------------------------------------------------
  // HistReader
  //   RAII wrapper around a READ-mode TFile. LoadInto reads a saved histogram
  //   by name into an already-booked destination histogram of matching
  //   binning (construct-empty-then-load: the plotting stage re-runs the same
  //   booking logic the histogramming stage used, then fills it from the
  //   file) -- this reuses the existing binning/title/name construction code
  //   instead of duplicating it in a second place.
  // ---------------------------------------------------------------------------
  class HistReader {
  public:
    explicit HistReader(const std::string& path)
      : file_(TFile::Open(path.c_str(), "READ")) {
      if (!file_ || file_->IsZombie())
        throw std::runtime_error("HistReader: failed to open '" + path + "' for reading");
    }

    ~HistReader() { Close(); }

    HistReader(const HistReader&) = delete;
    HistReader& operator=(const HistReader&) = delete;

    void Close() {
      if (file_) {
        file_->Close();
        delete file_;
        file_ = nullptr;
      }
    }

    // Loads the saved histogram named `name` into `dest` in place: clears
    // dest's own contents/errors/stats (but not its identity/binning) via
    // Reset("ICES"), validates the saved histogram's binning matches dest's,
    // then TH1::Add()s the saved contents in. Throws on a missing key or a
    // binning mismatch -- a stale file / current-code constant mismatch must
    // fail loudly, not silently misalign bins.
    template <typename H>
    void LoadInto(H* dest, const std::string& name) {
      if (!file_) throw std::runtime_error("HistReader: file not open");
      if (!dest)  throw std::runtime_error("HistReader::LoadInto: null destination");

      auto* src = file_->Get<H>(name.c_str());
      if (!src)
        throw std::runtime_error("HistReader::LoadInto: key '" + name + "' not found in file");

      if (src->GetNbinsX() != dest->GetNbinsX() ||
          src->GetNbinsY() != dest->GetNbinsY()) {
        throw std::runtime_error(
          "HistReader::LoadInto: binning mismatch for '" + name +
          "' (saved file vs. current booking constants)");
      }

      dest->Reset("ICES");
      dest->Add(src);
    }

    Long64_t ReadScalarLong(const std::string& name) {
      if (!file_) throw std::runtime_error("HistReader: file not open");
      auto* p = file_->Get<TParameter<Long64_t>>(name.c_str());
      if (!p) throw std::runtime_error("HistReader::ReadScalarLong: key '" + name + "' not found");
      return p->GetVal();
    }

    Double_t ReadScalarDouble(const std::string& name) {
      if (!file_) throw std::runtime_error("HistReader: file not open");
      auto* p = file_->Get<TParameter<Double_t>>(name.c_str());
      if (!p) throw std::runtime_error("HistReader::ReadScalarDouble: key '" + name + "' not found");
      return p->GetVal();
    }

    std::string ReadEnergyLabel() {
      if (!file_) throw std::runtime_error("HistReader: file not open");
      auto* s = file_->Get<TObjString>("meta_energy_label");
      if (!s) throw std::runtime_error("HistReader::ReadEnergyLabel: 'meta_energy_label' not found");
      return s->GetString().Data();
    }

    Long64_t ReadNEventsTotal() { return ReadScalarLong("meta_n_events_total"); }

    // -----------------------------------------------------------------------
    // CheckSelection
    //   Verifies the file was produced with the VBS selection this process is
    //   configured for, and aborts if not: the plot stage stamps its own
    //   selection onto every axis label and summary line, so reading a file
    //   made under a different one silently mislabels the output rather than
    //   failing.
    //
    //   A file written before meta_vbs_* existed cannot be checked -- warn
    //   loudly and continue, because refusing would make every pre-existing
    //   histogram file unreadable. That warning is the whole point: those are
    //   exactly the files whose name no longer implies their selection.
    // -----------------------------------------------------------------------
    void CheckSelection(double vbsDeta, double vbsMjj) {
      if (!file_) throw std::runtime_error("HistReader: file not open");
      auto* d = file_->Get<TParameter<Double_t>>("meta_vbs_deta");
      auto* m = file_->Get<TParameter<Double_t>>("meta_vbs_mjj");
      if (!d || !m) {
        std::cerr << "WARNING: this histogram file predates selection metadata, so the\n"
                     "         selection it was produced with CANNOT be verified. The\n"
                     "         defaults changed on 2026-08-26 (|Deta| >= 3, m_jj >= 0\n"
                     "         -> no |Deta|, m_jj >= 200), so an untagged file may have\n"
                     "         been made under the OLD selection. Plots will be labelled\n"
                     "         |Deta| >= " << vbsDeta << ", m_jj >= " << vbsMjj
                  << " GeV regardless. Re-run the\n"
                     "         histogramming stage if the provenance matters.\n";
        return;
      }
      const double tol = 1e-9;
      if (std::abs(d->GetVal() - vbsDeta) > tol || std::abs(m->GetVal() - vbsMjj) > tol) {
        std::ostringstream os;
        os << "HistReader::CheckSelection: this file was produced with |Deta| >= "
           << d->GetVal() << ", m_jj >= " << m->GetVal()
           << " GeV, but this process is configured for |Deta| >= " << vbsDeta
           << ", m_jj >= " << vbsMjj << " GeV. Plots would be mislabelled. "
              "Pass matching --vbs-deta=/--vbs-mjj= flags, or re-run the "
              "histogramming stage.";
        throw std::runtime_error(os.str());
      }
    }

  private:
    TFile* file_;
  };

}  // namespace MyUtl

#endif  // HISTOGRAM_IO_H
