#ifndef RPT_V5_COMMON_H
#define RPT_V5_COMMON_H

// ---------------------------------------------------------------------------
// rpt_v5_common.h
//   Shared between rpt_v5_hist.cxx (event loop, writes histograms to a ROOT
//   file) and rpt_v5_plot.cxx (reads that file, produces the RpT PDF) --
//   split out of the former single rpt_v5.cxx. Houses the Scenario container,
//   the shared RpT histogram binning/booking, and the save/load helpers that
//   round-trip a Scenario set through histogram_io.h.
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "histogram_io.h"

#include <TH1.h>
#include <TString.h>

#include <string>
#include <vector>

using namespace MyUtl;

// -----------------------------------------------------------------------------
// Scenario container.
// -----------------------------------------------------------------------------
struct Scenario {
  std::string name;
  std::string legend;
  Color_t color;
  TH1D* h_hs;
  TH1D* h_pu;
};

// -----------------------------------------------------------------------------
// Histogram binning + Scenario factory (file scope, not main()-local, so a
// TTreeProcessorMT worker thread can build one full set per thread, and so
// the plotting stage can build an identically-binned empty set to load into).
//
// Non-uniform binning: fine bins in [0, 2.5] for ROC granularity, then one
// wide bin to capture the tail without bloating memory.  250 bins of width
// 0.01 — coarse enough that the ROC + error bars aren't overcrowded.
// -----------------------------------------------------------------------------
inline const std::vector<double> rpt_bins = []() {
  std::vector<double> b;
  b.reserve(252);
  for (int i = 0; i <= 250; ++i) b.push_back(0.01 * i);
  b.push_back(375.0);
  return b;
}();
inline const int rpt_nbin = (int)rpt_bins.size() - 1;  // 251

inline TH1D* makeHist(const char* name, const char* title) {
  return new TH1D(name, title, rpt_nbin, rpt_bins.data());
}

// Scenario set: the no-timing baseline plus the three timing algorithms under
// comparison. zonly is kept (not one of the three) because a ROC has nothing to
// improve *over* without it -- it is the reference every ratio panel divides by.
//
// The earlier waves_misas (event-level HS-timing-purity oracle) and truth
// (perfect vertex t0 ceiling) scenarios were dropped to keep the region-split
// histogram count manageable; both are a two-line re-add here plus their time
// source in rpt_v5_hist.cxx if the ceiling/oracle context is wanted again.
inline std::vector<Scenario> makeScenarios(const std::string& suffix) {
  std::vector<Scenario> s = {
    {"zonly",       "ITk-only",                    C05, nullptr, nullptr},
    {"hgtd",        "HGTD t_{0} (Athena)",         C01, nullptr, nullptr},
    {"trkptz",      "TRKPTZ t_{0}",                C02, nullptr, nullptr},
    {"waves",       "WAVeS t_{0}",                 C03, nullptr, nullptr},
    // The real WAVeS vertex time, but with idealised 30 ps track times. Every
    // event and every jet still enters, so unlike an event-purity oracle this
    // changes no denominator. Paired with the truth row below it isolates one
    // variable: the two differ ONLY in where the vertex time comes from, so
    // the gap between them is the cost of our t0 selection and the gap from
    // WAVeS up to here is the cost of real per-track timing.
    {"waves_smear", "WAVeS t_{0}, 30 ps tracks",   C04, nullptr, nullptr},
    // Reference-study truth t0: vertex time = Gaus(truth HS vertex, 10 ps),
    // track time = Gaus(the track's OWN truth vertex, 30 ps) -- the 30 ps is
    // what carries HS/PU separation, since a pileup track inherits its own
    // vertex's time ~175 ps away. Gated with the same pull as everything else.
    {"truth",       "Truth t_{0} (10#oplus30 ps)",  C06, nullptr, nullptr},
  };
  for (auto& sc : s) {
    sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                       ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
    sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                       ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
  }
  return s;
}

// -----------------------------------------------------------------------------
// saveScenarios / loadScenarios
//   Persist/restore a Scenario set's raw histograms through a
//   HistWriter/HistReader. No derived data on Scenario to exclude (unlike
//   AnalysisObj) -- every field is event-loop-filled.
// -----------------------------------------------------------------------------
inline void saveScenarios(MyUtl::HistWriter& w, const std::vector<Scenario>& sv) {
  for (const auto& s : sv) {
    w.WriteHist(s.h_hs);
    w.WriteHist(s.h_pu);
  }
}

inline void loadScenarios(MyUtl::HistReader& r, std::vector<Scenario>& sv) {
  for (auto& s : sv) {
    r.LoadInto(s.h_hs, s.h_hs->GetName());
    r.LoadInto(s.h_pu, s.h_pu->GetName());
  }
}

#endif  // RPT_V5_COMMON_H
