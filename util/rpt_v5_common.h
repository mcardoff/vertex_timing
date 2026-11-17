#ifndef RPT_V5_COMMON_H
#define RPT_V5_COMMON_H

// ---------------------------------------------------------------------------
// rpt_v5_common.h
//   Shared between rpt_v5_hist.cxx (event loop, writes histograms to a ROOT
//   file) and rpt_v5_plot.cxx (reads that file, produces the RpT PDF) --
//   split out of the former single rpt_v5.cxx. Houses the Scenario container,
//   the shared RpT histogram binning/booking, the save/load helpers that
//   round-trip a Scenario set through histogram_io.h, and the per-sample
//   vertex-t0 inflation table (also used by util/vbs_time_veto.cxx).
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
    // Idealised WAVeS: our clustering, but with BOTH idealisations at once --
    // 30 ps track times, and the cluster selected by smallest |t - t_truth|
    // instead of by score. It is the ceiling of our approach, so the gap up
    // from the WAVeS row is what real times plus imperfect selection cost, and
    // the gap remaining to the truth row is what the CLUSTERING itself costs.
    // Needs truth for the selection, so it is an oracle, not deployable.
    {"waves_ideal", "Idealised WAVeS [oracle]",    C04, nullptr, nullptr},
    {"truth",       "Truth t_{0} (10#oplus30 ps)",  C06, nullptr, nullptr},
    // TZP: the classical selector (Score::TRKPTZ_TZQ, short name TRKPTZ_TZP)
    // with its guarded in-jet time -- see the TZP section in CLAUDE.md's
    // Physics Findings. Appended LAST so the load-bearing indices 0-5 are
    // untouched; fillJets fills sv[6] positionally.
    {"tzp",         "TZP t_{0}",                   C08, nullptr, nullptr},
  };
  for (auto& sc : s) {
    sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                       ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
    sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                       ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
  }
  return s;
}

// ── Vertex-time error inflation, per scenario ────────────────────────────────
// Shared by util/rpt_v5_hist.cxx (the R_pT time gate) and util/vbs_time_veto.cxx
// (the event-level jet-vs-t0 veto), so both use one calibration.
//
// The quoted vertex-time uncertainty understates the true spread, so a nominal
// N-sigma gate behaves like a much tighter one and discards genuine HS tracks.
// rpt_v6 measured this on the Athena vertex time (sigma_trk 27.0 ps, sigma_vtx
// 9.1 ps quoted vs a 50.8 ps observed core: a 1.78x understatement) and showed
// that correcting it turns a ratio that fell BELOW 1 at 0.875 efficiency into a
// sustained ~1.3-1.45.
//
// The factor is NOT shared: each scenario derives its vertex time differently,
// so each has its own calibration. Values below are measured by the
// PRINT_PULL_DIAG block at the end of util/rpt_v5_hist.cxx -- run it, read the
// "sigma ratio" column, and set these to it. Applied as sigma_vtx *= f, i.e. var_vtx
// *= f^2, replacing the previous blanket 2.25 (= 1.5^2) applied to every
// scenario alike.
//
// The factor is also PER SAMPLE, not just per scenario: how badly the quoted
// vertex-time error understates the truth depends on how often a usable vertex
// time exists at all, which varies enormously between samples (VBF 80.1% valid,
// dijet 73.6%, Z+jets only 16.5%). Z+jets needs ~1.6 where VBF needs ~1.4, so
// letting it inherit VBF's numbers ran its gate ~18% tighter than nominal,
// over-trimming genuine HS tracks -- which shows up downstream as a depressed
// maximum reachable HS efficiency (its RpT==0 spike, and so its ROC endpoint,
// is the worst of the three).
//
// MEASURED by rpt_v5_hist's PRINT_PULL_DIAG block (truth-HS tracks, |dt| < 150 ps
// core), read straight off its "ratio" column, one full-statistics grid run per
// sample:
//
//   scenario      vbf     zjets     dijet
//   hgtd         1.48      1.61      1.53
//   trkptz       1.39      1.63      1.45
//   waves        1.38      1.65      1.46
//
// The measurement is NOT circular: the accumulator uses the RAW quoted var_vtx
// and the ungated track list, and the vertex times themselves (cluster
// selection) do not depend on the gate -- so a single pass measures the ratio
// that the next run should apply, with no need to iterate to a fixed point.
//
// The Athena vertex time is the least well calibrated of the three on VBF, as
// expected; on Z+jets all three are comparably bad. Compared to rpt_v6's 1.78
// on the same Athena time, part of the difference is population rather than
// calibration: rpt_v6 measured over every timed HS track, while these are
// measured after the z-association, which already removes the worst outliers.
struct Inflation { double hgtd, trkptz, waves, tzp; };

// Keyed on MyUtl::SAMPLE_NAME. The local default run (no --sample, so an empty
// SAMPLE_NAME) reads the VBF ntuples, so it correctly falls through to the VBF
// row rather than needing an entry of its own.
inline Inflation inflationFor(const std::string& sample) {
  // tzp inflation seeded from the waves value per sample; the PRINT_PULL_DIAG
  // table prints the freshly measured ratio next to it -- update after the
  // first run per sample if they disagree.
  if (sample == "zjets") return {1.61, 1.63, 1.65, 1.65};
  if (sample == "dijet") return {1.53, 1.45, 1.46, 1.46};
  // SEEDED, NOT MEASURED: no rpt_v5 run has measured ttbar's ratios yet. Seeded
  // from the dijet row (ttbar has tracked dijet in every core-fraction study
  // this branch). Check this run's calibration table and replace.
  if (sample == "ttbar") return {1.53, 1.45, 1.46, 1.46};
  return {1.48, 1.39, 1.38, 1.38};  // vbf, and the local default run
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
