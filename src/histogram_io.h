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
    void WriteRunMeta(const std::string& energyLabel, Long64_t nEventsTotal) {
      if (!file_) throw std::runtime_error("HistWriter: file not open");
      file_->cd();
      TObjString label(energyLabel.c_str());
      label.Write("meta_energy_label", TObject::kOverwrite);
      WriteScalar("meta_n_events_total", nEventsTotal);
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

  private:
    TFile* file_;
  };

}  // namespace MyUtl

#endif  // HISTOGRAM_IO_H
