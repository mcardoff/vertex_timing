#ifndef RPT_V5_COMMON_H
#define RPT_V5_COMMON_H

// ---------------------------------------------------------------------------
// rpt_v5_common.h
//   Shared between rpt_v5_hist.cxx (event loop, writes histograms to a ROOT
//   file) and rpt_v5_plot.cxx (reads that file, produces the RpT PDF) --
//   split out of the former single rpt_v5.cxx. Houses the Scenario container,
//   the shared RpT histogram binning/booking, the save/load helpers that
//   round-trip a Scenario set through histogram_io.h, the per-sample
//   vertex-t0 inflation table (also used by util/vbs_time_veto.cxx), and the
//   --jvt parsing and file check both stages share.
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "histogram_io.h"
#include "jet_tagging.h"
#include "sample_config.h"   // recordRunSetting / recordDzPara

#include <TH1.h>
#include <TString.h>

#include <stdexcept>
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
// Shared by util/rpt_v5_hist.cxx (the R_pT time gate), util/vbs_time_veto.cxx
// (the event-level jet-vs-t0 veto) and the tdr_* programs (via util/tdr_common.h),
// so all of them use one calibration.
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

// One sample's row, and where its numbers came from:
//   MEASURED  the hgtd/trkptz/waves columns are PRINT_PULL_DIAG ratios from a
//             run on this sample (the table above).
//   SEEDED    no measurement on this sample is in use: the row is a copy of
//             another sample's row, the closest stand-in available. It stays a
//             copy (re-measuring the source moves it too) until it is replaced
//             by a MEASURED row of its own.
// The tzp column is seeded from the same row's waves value in EVERY row, so
// MEASURED describes the other three columns.
//
// Keyed on MyUtl::SAMPLE_NAME: a --sample= registry name (src/sample_config.h),
// or "" for the local default run. EVERY registered sample needs a row here,
// and an unknown name aborts. This used to fall through to the VBF row for any
// sample without one, which silently ran zeejets and the three mu=0 samples on
// VBF's factors. Registering a new sample means adding it below, measured or
// seeded.
struct InflationRow {
  Inflation   f;
  bool        measured;    // false = SEEDED
  std::string provenance;  // for the startup [calibration] echo
};

inline InflationRow inflationRow(const std::string& sample) {
  // MEASURED. Each row's tzp value is a copy of its waves value; the
  // PRINT_PULL_DIAG table prints the freshly measured tzp ratio beside it, so
  // replace the copy when a run's measurement disagrees. The 2026-09-01 grid
  // run measured tzp at 1.59 (zjets) and 1.43 (dijet) against the 1.65 and
  // 1.46 in use: 2-4%, under 1% on the gate width, so the copies were kept.
  static const std::map<std::string, Inflation> measured = {
    {"vbf",   {1.48, 1.39, 1.38, 1.38}},
    {"zjets", {1.61, 1.63, 1.65, 1.65}},
    {"dijet", {1.53, 1.45, 1.46, 1.46}},
  };

  // SEEDED: sample -> {the row it copies, a note for the echo}.
  struct Seed { std::string from, note; };
  static const std::map<std::string, Seed> seeded = {
    // Seeded from dijet, which ttbar has tracked in every core-fraction study.
    // The 2026-09-01 grid run then measured ttbar's own ratios (tzp 1.43
    // against this row's 1.46) and the dijet-seeded row held up in every
    // scenario, so the seed was kept rather than replaced by those numbers.
    {"ttbar",       {"dijet",   "checked against a 2026-09-01 ttbar measurement and kept"}},
    // Z->ee at mu=200, the 601189 half of zjets: zjets's row, not VBF's.
    {"zeejets",     {"zjets",   ""}},
    // mu=0: each seeded from its mu=200 counterpart, the pairing documented in
    // the registry (src/sample_config.h). These are PLACEHOLDERS: the factor
    // tracks how often a usable vertex time exists, which pileup changes, so a
    // pileup-free sample has no reason to share its counterpart's value. Run
    // PRINT_PULL_DIAG on the sample before relying on a gated mu=0 result.
    {"vbf_mu0",     {"vbf",     "mu=200 counterpart; not measured at mu=0"}},
    {"zeejets_mu0", {"zeejets", "mu=200 counterpart; not measured at mu=0"}},
    {"ttbar_mu0",   {"ttbar",   "mu=200 counterpart; not measured at mu=0"}},
  };

  // The local default run (no --sample) reads the local VBF ntuples, so it IS
  // the VBF row, not a fallback to it.
  if (sample.empty()) {
    InflationRow r = inflationRow("vbf");
    r.provenance += " (local default run reads VBF)";
    return r;
  }
  if (auto m = measured.find(sample); m != measured.end())
    return {m->second, true, "MEASURED on " + sample};
  if (auto s = seeded.find(sample); s != seeded.end()) {
    const Seed& seed = s->second;
    InflationRow r = inflationRow(seed.from);
    std::string p = "SEEDED from " + seed.from;
    if (!seed.note.empty()) p += " (" + seed.note + ")";
    if (!r.measured) p += ", itself " + r.provenance;
    r.measured   = false;
    r.provenance = p;
    return r;
  }
  std::cerr << "inflationFor: FATAL: no vertex-time inflation row for sample '"
            << sample << "'. Add one to inflationRow() in util/rpt_v5_common.h, "
               "MEASURED (rpt_v5_hist's PRINT_PULL_DIAG ratios on that sample) or "
               "SEEDED from the sample it most resembles. Refusing to borrow "
               "another sample's factors silently.\n";
  std::exit(1);
}

inline Inflation inflationFor(const std::string& sample) { return inflationRow(sample).f; }

// The startup echo, printed by every consumer once inflationFor() has resolved,
// so a run's own stdout records which row produced its outputs and whether that
// row was measured on this sample or seeded from another.
inline void printInflation(const std::string& sample) {
  const InflationRow r = inflationRow(sample);
  std::printf("[calibration] vertex-time inflation for '%s': "
              "hgtd %.2f  trkptz %.2f  waves %.2f  tzp %.2f  -- %s\n",
              sample.empty() ? "local (vbf)" : sample.c_str(),
              r.f.hgtd, r.f.trkptz, r.f.waves, r.f.tzp, r.provenance.c_str());
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

// -----------------------------------------------------------------------------
// recordRptAssoc
//   --dzpara (the WAVeS clustering input's association) and --rpt-signif (the
//   R_pT track list by z0 significance) each change every rpt_v5 histogram, so a
//   run with one tags its outputs ("dzpara", "rptsignif") and records the
//   setting in them (MyUtl::recordRunSetting). rpt_v5_hist passes what it runs
//   with, rpt_v5_plot the flags it was given: the same tags in the same order,
//   so the plot stage finds the file.
// -----------------------------------------------------------------------------
inline void recordRptAssoc(bool clusterDzPara, bool rptSignif) {
  MyUtl::recordDzPara(clusterDzPara);
  MyUtl::recordRunSetting("rpt_signif", "--rpt-signif", rptSignif ? 1.0 : 0.0, 0.0, "rptsignif");
}

// -----------------------------------------------------------------------------
// resolveRptJvt
//   --jvt=<wp> for both stages (src/jet_tagging.h's resolveJvtWP), so the plot
//   stage builds the same output tag as the hist stage and finds its file: a
//   working point that tags (every one but none) appends its tag (jvtLoose,
//   ...) to SELECTION_TAG. Call after resolveSelection. nullptr after an error
//   message, for an unknown working point or one rpt_v5 cannot evaluate.
// -----------------------------------------------------------------------------
inline const JvtWP* resolveRptJvt(int argc, char** argv, const std::string& prog) {
  const JvtWP* wp = resolveJvtWP(argc, argv, prog);
  if (wp && wp->rptMin >= 0.0) {
    const std::string tag = std::string(wp->tag).substr(1);
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  return wp;
}

// -----------------------------------------------------------------------------
// checkJvtWP
//   The plot stage's check that the file was made with the requested --jvt
//   working point, as HistReader::CheckSelection does for the VBS cuts: the
//   tag in the path no longer says which working point a --hist-file= file
//   has. rpt_v5_hist writes meta_vbs_jvt_wp (the index in JVT_WPS) for every
//   tagged run and has since --jvt existed (630cb9b), so a file without the
//   key was made without pileup-jet tagging: that is fine for --jvt=none and a
//   mismatch otherwise. Throws on a mismatch.
// -----------------------------------------------------------------------------
inline void checkJvtWP(MyUtl::HistReader& r, const JvtWP& wp) {
  const Long64_t want = static_cast<Long64_t>(&wp - JVT_WPS);
  const Long64_t nWP  = static_cast<Long64_t>(sizeof(JVT_WPS) / sizeof(JVT_WPS[0]));
  Long64_t have = 0;   // none: what a file without the key was made with
  r.FindScalarLong("meta_vbs_jvt_wp", have);
  if (have == want) return;
  const std::string haveName = (have >= 0 && have < nWP) ? JVT_WPS[have].name
                                                         : "unknown index " + std::to_string(have);
  throw std::runtime_error("checkJvtWP: this file was produced with --jvt=" + haveName +
                           ", but this process is configured for --jvt=" + wp.name +
                           ". Pass the matching --jvt=, or re-run the histogramming stage.");
}

#endif  // RPT_V5_COMMON_H
