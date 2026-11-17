// rpt_v5_hist.cxx — Histogramming stage of the rpt_v5 split.
//
// Runs the same TTreeProcessorMT event loop and per-thread ThreadState merge
// that the former monolithic rpt_v5.cxx did, then writes the merged
// Hard-Scatter/Pile-Up RpT histograms plus the scalar event-count/floor-counter
// accumulators via histogram_io.h, instead of doing any ROC/console-summary/PDF
// work directly. rpt_v5_plot.cxx reads that file back and does all of that --
// see CLAUDE.MD's "Main Executables" section.
//
// Scenarios: the seven rows of makeScenarios() (util/rpt_v5_common.h) --
// zonly, hgtd, trkptz, waves, waves_ideal, truth, tzp. Their indices are
// load-bearing (fillJets/fillRegion fill sv[0..6] positionally).
//
// No event-level selection besides vertex quality (and Z->ll for Z+jets) --
// every forward jet in the acceptance contributes an independent RpT
// measurement. The VBS knobs (--vbs-deta/--vbs-mjj) gate the regions only.
//
// Jet populations: 30-40 and >40 GeV at 2.4 < |eta| < 3.8 (forward,
// HGTD-covered), the same two slices at |eta| < 2.4 (central baseline), and
// the narrow VBS regions _r1/_r2 (both forward legs 2.4-3.8).
//
// Outputs (sample/shard/selection tags applied by histFilePath):
//   <prefix>rpt_v5_hist.root     histograms + scalars   -- merge with hist_merge
//   <prefix>rpt_v5_regions.root  TTree "regions": one row per WIDE-window
//                                (forward |eta| > 2.4, no upper edge) VBS-region
//                                event, every scenario's per-leg R_pT and gate
//                                inputs -- merge with hadd (see RegionRow)

#include <TChain.h>
#include <TH1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <ROOT/TTreeProcessorMT.hxx>

#include <boost/filesystem.hpp>
#include <algorithm>
#include <atomic>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#include "AtlasStyle.h"
#include "clustering_constants.h"
#include "sample_config.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "rpt_v5_common.h"
#include "histogram_io.h"
#include "idealised_timing.h"
#include "jet_tagging.h"

using namespace MyUtl;

// Jet acceptance for this script (paper values, distinct from rpt_v2).
static constexpr double JET_ETA_MIN = 2.4;
static constexpr double JET_ETA_MAX = 3.8;

// "Central" for the region-2 topology: inside HGTD's inner edge, i.e. a jet the
// detector has no timing coverage for at all. Deliberately the complement of
// JET_ETA_MIN rather than the reference study's tighter |eta| < 1.5 -- the
// physically meaningful boundary here is "HGTD can/can't see it".
static constexpr double CENTRAL_ETA_MAX = JET_ETA_MIN;

// Event-level HS-timing-purity threshold for the "clean timing" WAVeS scenario:
// a jet's RpT is filled only in events where ≥ this fraction of HS pT is timed
// within |pull| < 3σ (calcHSTimingPurity).  Mirrors the Score::WAVES_MISAS
// oracle.
static constexpr float MISAS_PURITY_CUT = 0.75f;

// Reference-study truth-t0 smearing: TRUTH_VTX_SMEAR (10 ps) and
// TRUTH_TRK_SMEAR (30 ps) now live in src/idealised_timing.h, alongside the
// idealised-world generation itself, so the event-display macro can replay
// exactly the draws these histograms were filled with.

// Per-track time-gate half-width in σ.  A 2σ cut over-trims genuine HS tracks
// when the vertex time is slightly mis-estimated, dragging the high-efficiency
// end of the ROC below ITk-only (worst in the >40 GeV slice).  Loosening it
// recovers that region at a small low-efficiency cost.
//
// 3.0 is the reference study's forward working point (jet_pileup_studies,
// 5 Oct 2023, slide 30: "For Forward eta region 2.4 <|eta|< 3.8, best setting
// is : dR<0.2, sigma_z0(pT,eta) cut < 1.4 & t30 <3.0"), chosen there by a scan
// over 2.0/2.5/3.0.  Note the two gates are not literally the same quantity but
// are close in effect: the reference divides by a FLAT 30 ps (every track in
// that study carries a 30 ps smear), so its t30 < 3.0 is a fixed +-90 ps window,
// while ours divides by the per-track pull denominator
// sqrt(infl^2 * var_vtx + sigma_trk^2), which for a typical HGTD track
// (sigma_trk ~26 ps, inflated sigma_vtx ~12 ps) is ~29 ps -> ~+-86 ps.
static constexpr double GATE_SIGMA = 3.0;

// ── Vertex-time error inflation, per scenario ────────────────────────────────
// The quoted vertex-time uncertainty understates the true spread, so a nominal
// N-sigma gate behaves like a much tighter one and discards genuine HS tracks.
// rpt_v6 measured this on the Athena vertex time (sigma_trk 27.0 ps, sigma_vtx
// 9.1 ps quoted vs a 50.8 ps observed core: a 1.78x understatement) and showed
// that correcting it turns a ratio that fell BELOW 1 at 0.875 efficiency into a
// sustained ~1.3-1.45.
//
// The factor is NOT shared: each scenario derives its vertex time differently,
// so each has its own calibration. Values below are measured by the
// PRINT_PULL_DIAG block at the end of this file -- run it, read the "sigma
// ratio" column, and set these to it. Applied as sigma_vtx *= f, i.e. var_vtx
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
// MEASURED by the PRINT_PULL_DIAG block below (truth-HS tracks, |dt| < 150 ps
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
static Inflation inflationFor(const std::string& sample) {
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

// Resolved once in main() before the event loop starts, then only read by the
// worker threads -- write-once-before-fork, so no synchronisation is needed.
static Inflation INFL = {1.48, 1.39, 1.38, 1.38};

// --jvt=none|loose|tight: pileup-jet tagging (src/jet_tagging.h) applied to
// the pT-passing jets BEFORE the VBS pair is formed, as an analysis would --
// the same taggers and working points vbs_region_diag --jvt uses. It gates the
// REGIONS only (the narrow _r1/_r2 histograms and the wide region tree), like
// the --vbs-* knobs: the inclusive forward slices stay the selection-free
// reference measurement. Resolved in main() before any worker starts.
static const JvtWP* JVT_SEL = &JVT_WPS[0];
static bool         JVT_ON  = false;

// Set true to print the measured per-scenario pull widths after the event loop.
static constexpr bool PRINT_PULL_DIAG = true;

static inline double dR(double j_eta, double j_phi, double t_eta, double t_phi) {
  double deta = j_eta - t_eta;
  double dphi = TVector2::Phi_mpi_pi(j_phi - t_phi);
  return std::sqrt(deta * deta + dphi * dphi);
}

// -----------------------------------------------------------------------------
// Track-to-vertex association, reference-study definition.
//
// getNewDzpara now lives in src/clustering_constants.h (MyUtl), shared with the
// clustering's own --dzpara association so the two cannot drift apart. It was
// ported verbatim from util/myJet_ana_fr.C, where the reference study calls it
// sigma_z0(pT, eta): an empirical parameterization of the OBSERVED z0 resolution,
// used in place of the per-track covariance sqrt(var_z0).
//
// The point of the substitution is that a significance cut is self-referential --
// it scales each track's window by that track's own error estimate, so a badly
// measured track earns a WIDER window. The parameterization sets the window from
// (eta, pT) alone, which the track cannot influence. The reference study scanned
// both and reported (5 Oct 2023, slide 30) the parameterization as the best
// forward choice and plain significance as the best central one.
// -----------------------------------------------------------------------------

// Multiple of getNewDzpara a track's |z0 - z_vtx| must stay within. 1.4 is the
// reference study's working point (its comment notes 0.8 discriminates better;
// left at 1.4 to match what the study actually ran).
static constexpr double DZ0_PARA_SCALE = 1.4;

// Forward R_pT track-to-vertex association. The parameterization is the default:
// the reference study scanned both and reported it as the best forward choice,
// and on local VBF it is worth roughly double the WAVeS improvement factor in the
// 30-40 GeV slice (+19% vs +5% at 0.80 HS efficiency, +26% vs +14% at 0.85).
//
// --rpt-signif flips this list back to a plain z0-significance cut at
// RPT_Z_SIGNIF, which is what rpt_v5 ran before 2026-08-13. Kept reachable
// because the choice moves the ITk-only baseline (126.6 -> 139.1 at 0.80) and
// therefore every "HGTD improves rejection by Nx" ratio quoted against it, so
// the two must stay comparable on demand.
//
// This is INDEPENDENT of MyUtl::USE_DZ_PARA (--dzpara), which governs the
// CLUSTERING's track selection. The two were briefly wired together; they are
// separate questions and are kept separate.
static constexpr double RPT_Z_SIGNIF = 2.5;
inline bool RPT_USE_DZ_PARA = true;

// Central (|eta| < 2.4) track-to-vertex association: |z0 - z_vtx| / sigma_z0.
//
// getNewDzpara comes from myJet_ana_fr.C -- a FORWARD-region study -- and does
// not extrapolate inward: it returns sigma_z0 = 31 um at eta = 0 rising to
// 4.4 mm at eta = 3.8, and the central end is several times tighter than ITk
// z0 resolution at 1 GeV. The reference study does not use it centrally either;
// its central plots scan a track significance cut instead.
//
// 5.0 is that scan's best working point, and an independent scan of our own
// (rpt_v6 --etamin=0 --etamax=1.5 --zcut=signif, ITk-only rejection at 0.800
// HS efficiency) reproduces it: para 77, then 90 / 124 / 145 / 156 / 163 for
// s = 2.0 / 2.5 / 3.0 / 4.0 / 5.0 -- monotonic, and every significance value
// beats the parameterization.
static constexpr double CENTRAL_Z_SIGNIF = 5.0;

// dR(track, jet) cone for the RpT numerator. The reference study applies this
// ON TOP of ghost association -- ghost-associated tracks outside the cone are
// dropped -- so it is strictly tighter than ghost association alone.
static constexpr double RPT_TRACK_JET_DR = 0.2;

// -----------------------------------------------------------------------------
// Event-display diagnostics: jet-level WAVeS ("mine") vs HGTD RpT comparison,
// restored from the pre-split rpt_v4.cxx (dropped somewhere in the rpt_v5
// rewrite). Ported with three changes: (a) WAVeS (t_waves/set_waves) stands
// in for v4's TRKPTZ, matching rpt_v5's WAVeS-based scenario set; (b) the new
// full-file-path event-display interface (--file_path) instead of the old
// fragile file-number-string extraction; (c) TTreeProcessorMT thread safety
// via per-thread ThreadState-local top-N vectors, merged after the event loop.
//
// Deliberately does NOT reintroduce v4's extra passBasicCuts()/passJetPtCut()
// (VBS jet-pair topology) gate on the jet-comparison diagnostic -- rpt_v5 as
// a whole intentionally has no such gate (every forward jet is an independent
// RpT measurement; see the file header), and restricting only this
// diagnostic to a smaller jet population than what's actually filled into the
// main RpT histograms would be confusing. hgtd_vtx_valid && waves_ok are
// kept, though: those aren't a topology cut, they guard against comparing
// against a degenerate/no-op time gate (applyTimeGate falls back to no
// gating at all when the vertex time is invalid).
// -----------------------------------------------------------------------------

// Set to true to print event-display commands to stdout after the event loop.
static constexpr bool PRINT_EVENT_DISPLAYS = true;

// How many event-display candidates to keep per region.
static constexpr int N_REGION_DISPLAYS = 12;

// -----------------------------------------------------------------------------
// RegionCase — an R1/R2 event worth drawing.
//
// This REPLACED the older JetCompCase/HurtJet display collection (WAVeS-vs-HGTD
// RpT disagreement, and HS jets "hurt" by the WAVeS time gate). Those ranked on
// criteria unrelated to VBS topology, so none of their events were guaranteed to
// be in a region at all -- and mixing them into the same stdout made it
// impossible to tell which commands were which. The regions are the focus now,
// so the display output is theirs exclusively.
//
// metric semantics differ by region, both signed so the printed line can say
// which direction timing moved things:
//   R1  margin(HS) - margin(PU), WAVeS minus no-timing. Positive => timing
//       widened the correct gap; negative => timing eroded or inverted it.
//   R2  forward-PU RpT, WAVeS minus no-timing. Negative => timing suppressed
//       the fake (the desired direction).
// -----------------------------------------------------------------------------
struct RegionCase {
  std::string file_path;
  Long64_t    entry;
  int    idx_hs, idx_pu;      // reco jet indices; idx_hs = -1 in R2
  double hs_pt, hs_eta;
  double pu_pt, pu_eta;
  double val_zonly, val_waves, metric;
  double t_waves;             // --extra_time annotation
};

static void insertRegionCase(std::vector<RegionCase>& v, RegionCase c,
                             int max_n = N_REGION_DISPLAYS) {
  v.push_back(std::move(c));
  std::sort(v.begin(), v.end(), [](const RegionCase& a, const RegionCase& b) {
    return std::abs(a.metric) > std::abs(b.metric);
  });
  if ((int)v.size() > max_n) v.resize(max_n);
}

// Merge one thread's top-N into the running top-N. Correct for the same reason
// as the old mergeCases: the global top-N is a subset of the union of each
// thread's already-truncated list.
static void mergeRegionCases(std::vector<RegionCase>& dst,
                             std::vector<RegionCase>& src,
                             int max_n = N_REGION_DISPLAYS) {
  for (auto& c : src) dst.push_back(std::move(c));
  std::sort(dst.begin(), dst.end(), [](const RegionCase& a, const RegionCase& b) {
    return std::abs(a.metric) > std::abs(b.metric);
  });
  if ((int)dst.size() > max_n) dst.resize(max_n);
}

// -----------------------------------------------------------------------------
// RegionRow — one WIDE-window VBS-region event, written to the side-file tree
// <prefix>rpt_v5_regions.root (tree "regions").
//
// "Wide" means forward = |eta| > JET_ETA_MIN with NO upper edge (central
// |eta| < CENTRAL_ETA_MAX): the same pair-level R1/R2 the VBS composition
// plots use, so legs past |eta| 3.8 -- and past HGTD/ITk at 4.0, where a jet
// has no tracks and R_pT is 0 under every scenario -- are included. The
// _r1/_r2 HISTOGRAMS above keep the narrow 2.4-3.8 window (their ROCs need a
// timeable HS leg); `core` flags the rows that are also in that narrow
// region, which is a strict subset (same pair, same paper labels, only the
// forward upper edge differs).
//
// One row per event, not per leg: R1 carries both forward legs, R2 its
// forward PU leg plus the central HS leg's kinematics (rpt_hs_* = -1 there --
// the central leg is outside HGTD and never a timed leg).
//
// The per-scenario arrays are index-aligned with makeScenarios() (checked at
// startup against SCEN_NAMES) and are written as NAMED scalar branches
// (rpt_hs_waves, t0_truth, ...), so readers never depend on the index.
//
// It lives in its own file because util/hist_merge.cxx refuses unknown key
// types: shards of this file are merged with hadd (TTrees only), the
// histogram file with hist_merge.
// -----------------------------------------------------------------------------
static constexpr int N_SCEN = 7;
static const char* const SCEN_NAMES[N_SCEN] = {
  "zonly", "hgtd", "trkptz", "waves", "waves_ideal", "truth", "tzp"};

// The real-t0 scenarios whose gate outcomes are attributed (see RegionRow's
// "why" columns): hgtd, trkptz, waves, as indices into SCEN_NAMES.
static constexpr int N_WHY = 3;
static constexpr int WHY_SCEN[N_WHY] = {1, 2, 3};

struct RegionRow {
  std::string file_path;
  Long64_t    entry  = -1;
  int         region = 0;      // 1 = R1, 2 = R2 (wide window)
  bool        core   = false;  // also in the narrow 2.4-3.8 region
  int         idx_hs = -1, idx_pu = -1;   // reco jet indices (R2: idx_hs central)
  double      mjj = -1.0, deta = -1.0;
  double      hs_pt = 0, hs_eta = 0, hs_phi = 0;
  double      pu_pt = 0, pu_eta = 0, pu_phi = 0;
  // z-associated ghost tracks within RPT_TRACK_JET_DR of each leg, and how
  // many of them carry a valid HGTD time -- what timing can act on at all.
  int         hs_ntrk = 0, hs_ntimed = 0, pu_ntrk = 0, pu_ntimed = 0;
  // Jets > MIN_JET_PT (not overlap-removed) in MIN_ABS_ETA_JET..MAX_ABS_ETA_JET.
  // vbs_region_diag's preselection requires >= 1; the regions here do not, so
  // this is what reproduces the composition-plot counts exactly.
  int         n_jets_fwd_acc = 0;
  double      rpt_hs[N_SCEN] = {}, rpt_pu[N_SCEN] = {};
  double      t0[N_SCEN] = {}, sig0[N_SCEN] = {}, infl[N_SCEN] = {};
  bool        ok[N_SCEN] = {};
  double      t_truth = 0.0;       // unsmeared TruthVtx_time[0]
  // --jvt: the legs' own tagger discriminants (-1 when no tagger is applied),
  // and which pT-passing jets the tagger removed before pairing -- what a
  // display needs to show why a visibly larger jet is not a leg.
  double      hs_jvt_rpt = -1, hs_fjvt = -1, pu_jvt_rpt = -1, pu_fjvt = -1;
  std::vector<int> rm_jvt, rm_fjvt;
  // WHY the gate moved pT, per real-t0 case (WHY_SCEN), judged against the
  // truth HS vertex time ONLY -- one counterfactual: would this track's gate
  // outcome change if t0 were t_HS (same sigma_t0 and inflation)? No pileup
  // truth is needed, which the grid samples do not carry. Sums of track pT
  // over the leg's R_pT cone (z-associated ghost tracks within dR < 0.2).
  //   HS leg (R1), HS-vertex tracks the gate REMOVED:
  //     t0far / t0near  kept with t0 = t_HS: the reco t0 caused it, split at
  //                     |t0 - t_HS| = PASS_SIGMA (a wrong cluster vs an offset)
  //     trk             removed even with t0 = t_HS: the track's own time is
  //                     incompatible with the HS time
  //   PU leg, non-HS tracks the gate KEPT:
  //     untimed (all cases), not0 (this case had no t0: no gate at all),
  //     intime (kept even with t0 = t_HS), t0far / t0near (removed with t0 = t_HS)
  double hs_hspt = 0.0, pu_pupt = 0.0, pu_keep_untimed = 0.0;   // cone totals
  double hs_rm_t0far[N_WHY] = {}, hs_rm_t0near[N_WHY] = {}, hs_rm_trk[N_WHY] = {};
  double pu_keep_not0[N_WHY] = {}, pu_keep_intime[N_WHY] = {};
  double pu_keep_t0far[N_WHY] = {}, pu_keep_t0near[N_WHY] = {};
};

// -----------------------------------------------------------------------------
// EventCase — a display candidate for the non-region categories.
//   Separate from RegionCase because these rank on event-level timing quality
//   rather than on a jet pair, and carry no leg indices: no --jet_idx is
//   emitted for them, so the wedges keep their truth colouring.
// -----------------------------------------------------------------------------
struct EventCase {
  std::string file_path;
  Long64_t    entry;
  double      t_show, metric, v1, v2, v3;
};

static void insertEventCase(std::vector<EventCase>& v, EventCase c,
                            int max_n = N_REGION_DISPLAYS) {
  v.push_back(std::move(c));
  std::sort(v.begin(), v.end(), [](const EventCase& a, const EventCase& b) {
    return a.metric > b.metric; });
  if ((int)v.size() > max_n) v.resize(max_n);
}

static void mergeEventCases(std::vector<EventCase>& dst,
                            std::vector<EventCase>& src,
                            int max_n = N_REGION_DISPLAYS) {
  for (auto& c : src) dst.push_back(std::move(c));
  std::sort(dst.begin(), dst.end(), [](const EventCase& a, const EventCase& b) {
    return a.metric > b.metric; });
  if ((int)dst.size() > max_n) dst.resize(max_n);
}

// -----------------------------------------------------------------------------
// ThreadState
//   Everything one worker thread accumulates across whatever task ranges it
//   services: its own copy of the two pT-slice Scenario sets (each worker's
//   histograms are identically named to every other worker's -- harmless,
//   since TH1::AddDirectory(kFALSE) means none of them register into any
//   global directory) plus the event/floor-diagnostic counters. Built once
//   per worker thread (lazily, on first use) and merged into one after the
//   event loop.
// -----------------------------------------------------------------------------
struct ThreadState {
  // Re-seeded per event from event-stable quantities, so the truth-t0 smearing
  // does not depend on how TTreeProcessorMT distributes entries across threads.
  TRandom3 rng{1};

  std::vector<Scenario> scen_lo = makeScenarios("_lo");  // 30–40 GeV, forward
  std::vector<Scenario> scen_hi = makeScenarios("_hi");  // >40 GeV,   forward
  // Same two pT slices at CENTRAL eta, as an ITk-only baseline: |eta| < 2.4 is
  // outside HGTD acceptance, so its tracks have no valid time, the gate is a
  // no-op there and all four scenarios should land on top of each other. That
  // makes these both the reference the forward gain is measured against AND a
  // standing check that the timing machinery stays inert where it has no data.
  std::vector<Scenario> scen_lo_cen = makeScenarios("_lo_cen");  // 30–40 GeV, |eta|<2.4
  std::vector<Scenario> scen_hi_cen = makeScenarios("_hi_cen");  // >40 GeV,   |eta|<2.4
  // VBS-topology regions (see classifyRegion below). pT-inclusive (>MIN_JET_PT)
  // rather than sliced: both are rare topologies and splitting them further
  // would leave the ROCs statistics-limited.
  std::vector<Scenario> scen_r1 = makeScenarios("_r1");  // both VBS legs forward
  std::vector<Scenario> scen_r2 = makeScenarios("_r2");  // fwd PU leg + central HS leg
  long n_total = 0, n_pass_basic = 0, n_hgtd_valid = 0;
  // Z+jets-only breakdown (see event_processing.h EventResult::code doc for the
  // analogous clustering-side counters): n_pass_basic above is vertex-quality
  // only (checked BEFORE lepton selection here); n_pass_lepton_sel further
  // requires the Z->ll selection. No-ops (both equal n_pass_basic) elsewhere.
  long n_pass_lepton_sel = 0;
  // Finer breakdown of the lepton-selection failures (mirrors event_processing.h's
  // EventResult::code -2/-6/-7): among vertex-passing events that fail
  // passLeptonSelection, how many had 0 / exactly 1 / >=2 qualifying
  // (isGoodLepton) leptons but no OS-SF pair.
  long n_rej_no_lepton = 0, n_rej_one_lepton = 0, n_rej_no_ossf_pair = 0;
  double pu_tot_pt = 0, pu_floor_pt = 0, hs_tot_pt = 0, hs_floor_pt = 0;  // >40
  double pu_tot_lo = 0, pu_floor_lo = 0, hs_tot_lo = 0, hs_floor_lo = 0;  // 30-40

  // Event-display candidates, R1/R2 only (see RegionCase doc comment above).
  std::vector<RegionCase> cases_r1, cases_r2, cases_r2_fail;
  std::vector<EventCase>  cases_mis, cases_wwin;
  // Wide-window region events for the side-file tree (see RegionRow), and a
  // standing check that the narrow region really is a subset of the wide one.
  std::vector<RegionRow> region_rows;
  long n_narrow_not_subset = 0;
  // Per-scenario pull-width accumulators (see PRINT_PULL_DIAG).
  double pull_dt2_hgtd = 0, pull_var_hgtd = 0;
  double pull_dt2_trkptz = 0, pull_var_trkptz = 0;
  double pull_dt2_waves = 0, pull_var_waves = 0;
  double pull_dt2_tzp = 0, pull_var_tzp = 0;
  long   pull_n_hgtd = 0, pull_n_trkptz = 0, pull_n_waves = 0, pull_n_tzp = 0;
  // First moment and out-of-core counts, so the diagnostic can separate a
  // systematic OFFSET (which an inflation cannot fix) from genuine spread, and
  // report how much of each scenario's distribution the core window discards.
  double pull_dt_hgtd = 0, pull_dt_trkptz = 0, pull_dt_waves = 0, pull_dt_tzp = 0;
  long   pull_ntail_hgtd = 0, pull_ntail_trkptz = 0, pull_ntail_waves = 0, pull_ntail_tzp = 0;
  // Per-jet effect of the WAVeS time gate vs ITk-only, forward slices.
  // 'up' must stay 0: applyTimeGate returns a SUBSET of its input, so R_pT can
  // only fall. It is counted anyway as a standing assertion on that property.
  long rpt_n[2] = {0, 0}, rpt_up[2] = {0, 0}, rpt_same[2] = {0, 0};
  long rpt_down[2] = {0, 0}, rpt_zeroed[2] = {0, 0};
  double rpt_relloss[2] = {0.0, 0.0};
};

// RpT numerator: sum pT of tracks that are ghost-associated to the jet AND in
// the caller's z₀-selected set AND within RPT_TRACK_JET_DR of the jet axis.
// The dR term is why the jet direction has to be passed in -- it is per
// (track, jet), so unlike the z₀ cut it cannot be folded into assoc_set.
static double computeRpT(BranchPointerWrapper* b,
                         const std::vector<int>& ghost_indices,
                         double j_pt, double j_eta, double j_phi,
                         const std::unordered_set<int>& assoc_set) {
  double sumpt = 0.0;
  for (int idx : ghost_indices) {
    if (!assoc_set.count(idx)) continue;
    if (dR(j_eta, j_phi, b->trackEta[idx], b->trackPhi[idx]) > RPT_TRACK_JET_DR) continue;
    sumpt += b->trackPt[idx];
  }
  return sumpt / j_pt;
}

int main(int argc, char** argv) {
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  // Must precede any histogram construction -- see the identical note in
  // src/clustering_hist.cxx's main(). Every worker thread's ThreadState
  // produces identically-named histograms, harmless only because of this.
  TH1::AddDirectory(kFALSE);

  // SCEN_NAMES names the region tree's per-scenario branches, so it must stay
  // index-aligned with makeScenarios(). Fail loudly rather than mislabel a
  // column if a scenario is ever inserted or reordered.
  {
    auto chk = makeScenarios("_namecheck");
    bool aligned = ((int)chk.size() == N_SCEN);
    for (int k = 0; aligned && k < N_SCEN; ++k) aligned = (chk[k].name == SCEN_NAMES[k]);
    for (auto& s : chk) { delete s.h_hs; delete s.h_pu; }
    if (!aligned) {
      std::cerr << "SCEN_NAMES is out of sync with makeScenarios() -- fix before running.\n";
      return 1;
    }
  }

  // Flushed phase timestamps -- see the identical rationale in
  // src/clustering_hist.cxx and MyUtl::PhaseTimer.
  MyUtl::PhaseTimer phase;

  // --- Sample selection (--sample=vbf|zjets|dijet; default: local VBF ntuple) ---
  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);  // --vbs-deta=<x>; sets SELECTION_TAG
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  MyUtl::FILE_SHARD   = MyUtl::resolveShard(argc, argv);
  MyUtl::OVERLAP_REMOVAL = sample.overlapRemoval;  // Z+jets lepton–jet overlap removal
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  if (MyUtl::SAMPLE_NAME.empty())
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/hists");
  unsigned nThreads = MyUtl::resolveThreads(argc, argv);

  // --dzpara: swap the z0-significance association for the reference study's
  // getNewDzpara parameterization, in BOTH the R_pT track list and the WAVeS
  // clustering input (they share passTrackVertexAssociation). Set before any
  // worker thread starts; read-only afterwards. Default off -- see the trk_all
  // comment for why the choice moves every improvement factor.
  for (int i = 1; i < argc; ++i) {
    std::string a(argv[i]);
    if (a == "--dzpara")     MyUtl::USE_DZ_PARA = true;   // clustering input
    if (a == "--rpt-signif") RPT_USE_DZ_PARA    = false;  // R_pT track list
  }
  for (int i = 1; i < argc; ++i) {
    const std::string a(argv[i]);
    if (a.rfind("--jvt=", 0) != 0) continue;
    const JvtWP* wp = jvtWPByName(a.substr(6));
    if (!wp) { std::cerr << "unknown " << a << " (none|loose|tight)\n"; return 1; }
    JVT_SEL = wp;
  }
  JVT_ON = (JVT_SEL->rptMin >= 0.0);
  if (JVT_ON) {
    // The discriminants need the vertex fit's track assignment, an EXTENDED
    // branch; set before any BranchPointerWrapper binds.
    MyUtl::EXTENDED_BRANCHES = true;
    // A tagged run must not overwrite an untagged one: jvtLoose/jvtTight join
    // SELECTION_TAG, so every output name (hist file and region tree) carries it.
    const std::string tag = std::string(JVT_SEL->tag).substr(1);
    MyUtl::SELECTION_TAG = MyUtl::SELECTION_TAG.empty() ? tag : MyUtl::SELECTION_TAG + "_" + tag;
  }
  std::cout << "[jvt] " << JVT_SEL->name;
  if (JVT_ON)
    std::cout << " before VBS pairing: keep R_pT(vtx fit) >= " << JVT_SEL->rptMin
              << " for |eta| < " << JVT_ETA_MAX << ", pT < " << JVT_PT_MAX
              << "; fJVT <= " << JVT_SEL->fjvtMax << " for " << FJVT_ETA_MIN
              << " <= |eta| < " << FJVT_ETA_MAX << ", pT < " << FJVT_PT_MAX
              << "  (tag: " << MyUtl::SELECTION_TAG << ")";
  std::cout << '\n';
  std::cout << "[assoc] R_pT fwd: "
            << (RPT_USE_DZ_PARA ? "getNewDzpara x " + std::to_string(DZ0_PARA_SCALE)
                                : "z0 significance < " + std::to_string(RPT_Z_SIGNIF))
            << " | R_pT central: z0 significance < " << CENTRAL_Z_SIGNIF
            << " | clustering: "
            << (MyUtl::USE_DZ_PARA ? "getNewDzpara" : "z0 significance") << '\n';

  // Per-sample vertex-time calibration. Resolved here, before any worker thread
  // exists, so the event loop only ever reads it. Echoed so a run's stdout
  // records which calibration produced its histograms -- the PRINT_PULL_DIAG
  // table at the end then prints the freshly measured ratio beside these in its
  // "in use" column, making a stale entry visible in the same log.
  INFL = inflationFor(MyUtl::SAMPLE_NAME);
  std::printf("[calibration] vertex-time inflation for '%s': "
              "hgtd %.2f  trkptz %.2f  waves %.2f  tzp %.2f\n",
              MyUtl::SAMPLE_NAME.empty() ? "local (vbf)" : MyUtl::SAMPLE_NAME.c_str(),
              INFL.hgtd, INFL.trkptz, INFL.waves, INFL.tzp);

  TChain chain("ntuple");
  // setupChain now validates the file list itself (and aborts if empty) with
  // no file I/O. The chain.GetEntries() == 0 check that used to sit here
  // opened every file in the chain -- 0.3-1.2 s each on the AF's /data, so
  // 5-30 minutes of dead time -- to answer a question the file list already
  // answers. See setupChain's note in src/event_processing.h.
  setupChain(chain, sample.ntupleDir.c_str(), MyUtl::FILE_SHARD);
  if (JVT_ON) {
    // With EXTENDED_BRANCHES on, a branch the sample lacks (Track_btagIp_* on
    // every grid sample) would make the reader iterate ZERO entries and exit 0.
    recordAvailableBranches(chain);
    if (!MyUtl::hasBranch("Track_recoVtx_idx")) {
      std::cerr << "--jvt needs Track_recoVtx_idx, which this sample lacks; a silent "
                   "fallback would give every jet a JVT R_pT of 0 and remove them all\n";
      return 1;
    }
  }
  phase.mark("chain built");
  ROOT::EnableImplicitMT(nThreads);

  // --- Per-thread state registry, merged into one after the event loop.
  //     ThreadState (via Scenario's raw TH1D*) is trivially copy-constructible,
  //     which would make ROOT::TThreadedObject silently *shallow*-copy the
  //     TH1D pointers across "cloned" slots -- every worker would then Fill()
  //     the same histograms with no synchronization. Use the same hand-rolled
  //     mutex-guarded-registry + thread_local-pointer-cache pattern as
  //     src/clustering_hist.cxx's per-thread AnalysisObj map instead.
  //
  //     The mutex also serializes the ThreadState construction itself, not
  //     just the push_back -- defensive, matching a real crash found in
  //     clustering_dt.cxx's AnalysisObj construction (SetFillColorAlpha()
  //     racing on ROOT's global TColor table). ThreadState's construction
  //     doesn't call any color-touching ROOT functions today, but the
  //     serialization is a one-time per-thread cost, so there's no reason
  //     to leave that door open. ---
  std::mutex stateRegistryMutex;
  std::vector<std::unique_ptr<ThreadState>> stateRegistry;

  std::atomic<Long64_t> progressCounter{0};

  // Optional --max-events=<N> cap for quick local checks (see resolveMaxEvents
  // in sample_config.h). -1 means unlimited -- every existing invocation
  // without the flag is unaffected. The denominator is deliberately NOT
  // chain.GetEntries(): see the note above and in setupChain. Uncapped, the
  // total is simply unknown and progress prints as a bare count.
  const Long64_t maxEvents     = MyUtl::resolveMaxEvents(argc, argv);
  const Long64_t progressDenom = (maxEvents > 0) ? maxEvents : -1;
  if (maxEvents > 0)
    std::cout << "Restricting to first " << progressDenom << " events (--max-events)\n";

  ROOT::TTreeProcessorMT proc(chain, nThreads);
  proc.Process([&](TTreeReader& reader) {
    // Fresh per invocation: a worker thread can be handed a different
    // TTreeReader across task ranges, so this cannot be thread_local.
    BranchPointerWrapper branch(reader);

    // Lazily build this thread's state once; reused across however many
    // task ranges this worker thread services.
    thread_local ThreadState* tlState = nullptr;
    if (!tlState) {
      std::lock_guard<std::mutex> lock(stateRegistryMutex);
      stateRegistry.push_back(std::make_unique<ThreadState>());
      tlState = stateRegistry.back().get();
    }
    ThreadState& state = *tlState;

    // Redefined fresh per invocation alongside branch (cheap -- no
    // precomputation, just captures branch by reference).
    // Paper jet-label helper (ATL-HGTD-PUB-2022-001 Sec. 3).
    // HS  : dR(reco, truthHS) < 0.3  AND  truthHS_pT > 10 GeV
    auto paperIsHS = [&](double j_eta, double j_phi) {
      for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
        if (branch.truthHSJetPt[t] < 10.0) continue;
        if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.3)
          return true;
      }
      return false;
    };
    // PU: dR > 0.6 from any truth HS jet with pT > 4 GeV (Sec. 3).
    auto paperIsPU = [&](double j_eta, double j_phi) {
      for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
        if (branch.truthHSJetPt[t] < 4.0) continue;
        if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.6)
          return false;
      }
      return true;
    };

    while (reader.Next()) {
      Long64_t n = ++progressCounter;
      // Stop this worker's current task range once the shared counter has
      // already crossed the cap; see the identical pattern/rationale in
      // src/clustering_hist.cxx.
      if (maxEvents > 0 && n > maxEvents) break;
      ++state.n_total;

      if (n % 5000 == 0)
      {
        // '\r' only when stdout is a terminal; under condor it is a file,
        // where overwriting one line yields an unreadable mega-line and no
        // sense of rate. See MyUtl::STDOUT_IS_TTY.
        if (MyUtl::STDOUT_IS_TTY)
          std::cout << "Progress: " << n << (progressDenom > 0
                       ? "/" + std::to_string(progressDenom) : "") << "\r" << std::flush;
        else if (n % 100000 == 0)
          phase.mark("processed " + std::to_string(n) + " events");
      }

      // ── Require only vertex quality (paper Sec. 3: |z_reco − z_truth| < 2 mm).
      if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
      if (std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]) > MAX_VTX_DZ) continue;
      ++state.n_pass_basic;

      // Lepton–jet overlap removal (Z+jets only): flag reco jets within
      // LEPTON_JET_DR of a lepton track. No-op otherwise. SKIP_EVENT mode
      // vetoes the whole event on any overlap; REMOVE_JETS mode instead has
      // fillJets skip the flagged jets via isJetRemoved.
      branch.computeOverlapRemoval();
      // Z→ℓℓ selection (Z+jets only): require an opposite-sign same-flavour
      // lepton pair with pt > LEPTON_MIN_PT; no-op on other samples. Classify
      // the failure the same way event_processing.h's EventResult::code does
      // (0 / 1 / >=2 qualifying leptons) so a low survival rate can be
      // attributed to too-few-leptons vs. a pairing-logic issue.
      if (!branch.passLeptonSelection()) {
        int nGood = branch.countGoodLeptons();
        if      (nGood == 0) ++state.n_rej_no_lepton;
        else if (nGood == 1) ++state.n_rej_one_lepton;
        else                 ++state.n_rej_no_ossf_pair;
        continue;
      }
      ++state.n_pass_lepton_sel;
      if (branch.vetoLeptonOverlap()) continue;

      // ── Track selection: all tracks (no eta cut) associated to the primary
      //    vertex, for the z-only baseline / ITk-only scenario.
      //
      //    Defaults to the reference study's getNewDzpara parameterization; see
      //    RPT_USE_DZ_PARA above, and --rpt-signif for the z0-significance form.
      //
      //    The choice matters for how results are REPORTED, not just internally:
      //    it moves the ITk-only baseline (126.6 significance vs 139.1
      //    parameterized, forward 30-40 GeV at 0.80 HS efficiency on local VBF),
      //    and every "HGTD improves rejection by Nx" number is a ratio to that
      //    baseline. The same WAVeS curve reads as +5% or +19% depending only on
      //    which association the denominator used. Always quote the association
      //    alongside the factor.
      std::vector<int> trk_all;
      const double vtxZ = branch.recoVtxZ[0];
      for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
        double trkPt = branch.trackPt[trk];
        if (trkPt < MIN_TRACK_PT || trkPt > MAX_TRACK_PT) continue;
        if (!branch.trackQuality[trk]) continue;
        if (RPT_USE_DZ_PARA) {
          double dz = std::abs(branch.trackZ0[trk] - vtxZ);
          if (dz / MyUtl::getNewDzpara(branch.trackEta[trk], trkPt) > DZ0_PARA_SCALE) continue;
        } else if (!passTrackVertexAssociation((int)trk, 0, &branch, RPT_Z_SIGNIF)) {
          continue;
        }
        trk_all.push_back((int)trk);
      }

      // Central counterpart of the list above, on a z0-significance cut rather
      // than the forward parameterization (see CENTRAL_Z_SIGNIF). Kept as a
      // separate list rather than switching per track inside one loop: a
      // forward jet's dR < 0.2 cone can reach below |eta| = 2.4, and mixing
      // the two associations inside a single list would silently change the
      // forward numbers this baseline is meant to be compared against.
      std::vector<int> trk_all_cen;
      for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
        double trkPt = branch.trackPt[trk];
        if (trkPt < MIN_TRACK_PT || trkPt > MAX_TRACK_PT) continue;
        if (!branch.trackQuality[trk]) continue;
        double vz = branch.trackVarZ0[trk];
        if (!(vz > 0)) continue;
        double dz = std::abs(branch.trackZ0[trk] - vtxZ);
        if (dz / std::sqrt(vz) > CENTRAL_Z_SIGNIF) continue;
        trk_all_cen.push_back((int)trk);
      }

      // ── HGTD-acceptance tracks only (used for WAVeS clustering). ────────────
      std::vector<int> trk_z = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 2.5);

      // ── WAVeS clustering + selection. ────────────────────────────────────────
      auto clusters = clusterTracksInTime(
          trk_z, &branch, DIST_CUT_CONE,
          /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
          ClusteringMethod::ITERATIVE, /*useZ0=*/false,
          /*sortTracks=*/false, /*calcPurityFlag=*/true);

      // WAVeS selection: highest WAVeS-score cluster; time via in-jet refinement.
      double t_waves = 0.0, var_waves = 0.0;
      bool   waves_ok = false;
      // TRKPTZ selection: the baseline Σ pT e^{-1.5|Δz|} score, over the SAME
      // cluster collection -- so waves-vs-trkptz here isolates the choice of
      // selection score, with clustering held fixed.
      double t_trkptz = 0.0, var_trkptz = 0.0;
      bool   trkptz_ok = false;
      double t_tzp = 0.0, var_tzp = 0.0;
      bool   tzp_ok = false;
      if (!clusters.empty()) {
        auto best   = chooseCluster(clusters, Score::WAVES);
        t_waves     = best.calculateTime(Score::WAVES, &branch);  // in-jet refined
        var_waves   = best.sigmas[0] * best.sigmas[0];
        waves_ok    = true;

        auto bestT  = chooseCluster(clusters, Score::TRKPTZ);
        t_trkptz    = bestT.calculateTime(Score::TRKPTZ, &branch);
        var_trkptz  = bestT.sigmas[0] * bestT.sigmas[0];
        trkptz_ok   = true;

        // TZP: the classical selector, same cluster collection -- so
        // tzp-vs-trkptz isolates the selection + guarded in-jet timing.
        auto bestP  = chooseCluster(clusters, Score::TRKPTZ_TZQ);
        t_tzp       = bestP.calculateTime(Score::TRKPTZ_TZQ, &branch);  // guarded in-jet
        var_tzp     = bestP.sigmas[0] * bestP.sigmas[0];
        tzp_ok      = true;
      }

      // ── HGTD ntuple vertex time. ─────────────────────────────────────────────
      double t_hgtd         = branch.recoVtxTime[0];
      double var_hgtd       = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
      bool   hgtd_vtx_valid = (branch.recoVtxValid[0] == 1);
      if (hgtd_vtx_valid) ++state.n_hgtd_valid;

      // ── Per-track time gate.  pull width ~1.5 → var_vtx ×2.25. ──────────────
      auto applyTimeGate = [&](const std::vector<int>& base,
                                double t_vtx, double var_vtx, bool vtx_valid,
                                double sigma = 2.0, double infl = 1.5) {
        std::vector<int> out;
        out.reserve(base.size());
        for (int idx : base) {
          bool apply = vtx_valid && branch.trackTimeValid[idx] == 1;
          if (!apply) { out.push_back(idx); continue; }
          double dt    = branch.trackTime[idx] - t_vtx;
          double var_t = branch.trackTimeRes[idx] * branch.trackTimeRes[idx];
          double pull  = std::abs(dt) / std::sqrt(infl * infl * var_vtx + var_t);
          if (pull < sigma) out.push_back(idx);
        }
        return out;
      };

      std::vector<int> trk_hgtd   = applyTimeGate(trk_all, t_hgtd,   var_hgtd,   hgtd_vtx_valid, GATE_SIGMA, INFL.hgtd);
      std::vector<int> trk_waves  = applyTimeGate(trk_all, t_waves,  var_waves,  waves_ok,       GATE_SIGMA, INFL.waves);
      std::vector<int> trk_trkptz = applyTimeGate(trk_all, t_trkptz, var_trkptz, trkptz_ok,      GATE_SIGMA, INFL.trkptz);
      std::vector<int> trk_tzp    = applyTimeGate(trk_all, t_tzp,    var_tzp,    tzp_ok,         GATE_SIGMA, INFL.tzp);
      // Central list gets the identical gate. It is a no-op there in practice
      // (no HGTD coverage below |eta| 2.4, so no track carries a valid time),
      // but applying it keeps all four central scenarios defined exactly as
      // their forward counterparts -- which is what makes their agreement a
      // meaningful check rather than a tautology.
      std::vector<int> trk_hgtd_cen   = applyTimeGate(trk_all_cen, t_hgtd,   var_hgtd,   hgtd_vtx_valid, GATE_SIGMA, INFL.hgtd);
      std::vector<int> trk_waves_cen  = applyTimeGate(trk_all_cen, t_waves,  var_waves,  waves_ok,       GATE_SIGMA, INFL.waves);
      std::vector<int> trk_trkptz_cen = applyTimeGate(trk_all_cen, t_trkptz, var_trkptz, trkptz_ok,      GATE_SIGMA, INFL.trkptz);
      std::vector<int> trk_tzp_cen    = applyTimeGate(trk_all_cen, t_tzp,    var_tzp,    tzp_ok,         GATE_SIGMA, INFL.tzp);

      // ── Pull-width measurement, one accumulator set per scenario ───────────
      // Truth-HS tracks only (trackToTruthvtx == 0) with a valid time, in
      // events where that scenario produced a vertex time. Accumulates the
      // observed dt spread inside a generous core window against the QUOTED
      // uncertainty; the ratio of the two is the inflation factor to set above.
      if (PRINT_PULL_DIAG) {
        auto accum = [&](double t_vtx, double var_vtx, bool ok,
                         double& sum_dt2, double& sum_var, long& n,
                         double& sum_dt, long& n_tail) {
          if (!ok) return;
          for (int idx : trk_all) {
            if (branch.trackTimeValid[idx] != 1) continue;
            if (branch.trackToTruthvtx[idx] != 0) continue;   // truth-HS only
            double dt = branch.trackTime[idx] - t_vtx;
            if (std::abs(dt) > 150.0) { ++n_tail; continue; } // core window
            double var_t = branch.trackTimeRes[idx] * branch.trackTimeRes[idx];
            sum_dt2 += dt * dt;
            sum_dt  += dt;
            sum_var += var_vtx + var_t;
            ++n;
          }
        };
        accum(t_hgtd,   var_hgtd,   hgtd_vtx_valid, state.pull_dt2_hgtd,   state.pull_var_hgtd,   state.pull_n_hgtd, state.pull_dt_hgtd, state.pull_ntail_hgtd);
        accum(t_trkptz, var_trkptz, trkptz_ok,      state.pull_dt2_trkptz, state.pull_var_trkptz, state.pull_n_trkptz, state.pull_dt_trkptz, state.pull_ntail_trkptz);
        accum(t_waves,  var_waves,  waves_ok,       state.pull_dt2_waves,  state.pull_var_waves,  state.pull_n_waves, state.pull_dt_waves, state.pull_ntail_waves);
        accum(t_tzp,    var_tzp,    tzp_ok,         state.pull_dt2_tzp,    state.pull_var_tzp,    state.pull_n_tzp,   state.pull_dt_tzp,   state.pull_ntail_tzp);
      }

      // ── Idealised-timing reference scenarios ──────────────────────────────
      // Smeared times are generated ONCE per track here and then used for both
      // the re-clustering below and the gates, so the vertex time and the track
      // times a row is gated against always come from the same world. Building
      // them separately would mean clustering on one draw and gating on
      // another, which is what made the earlier version incoherent.
      //
      // The generation itself lives in src/idealised_timing.h (see its header
      // for the seeding / draw-order contract), shared with the event-display
      // macro so a display of the truth or waves_ideal row shows exactly the
      // times these histograms were filled with.
      const IdealisedTiming ideal = makeIdealisedTiming(&branch, state.rng);
      const auto& smTimes = ideal.times;

      // Re-cluster in the idealised world, then select the cluster closest in
      // time to truth rather than the highest-scoring one (the second
      // idealisation -- waves_ideal never uses the WAVeS score). The cluster
      // structure depends on the track times, so reusing the t0 built from the
      // real times would evaluate good tracks against a vertex time derived
      // from noisy ones. values[0] is the cluster's own weighted mean of the
      // smeared times. Ranked on |dt|, not purity: purity was tried and fails
      // badly (HS jets at R_pT = 0 rise 3.0% -> 11.8%) -- see
      // closestToTruthCluster.
      double t_wsm = 0.0, var_wsm = 0.0;
      bool   wsm_ok = false;
      {
        const auto cl_sm = idealisedClusters(trk_z, &branch, ideal);
        const int  bi    = closestToTruthCluster(cl_sm, branch.truthVtxTime[0]);
        if (bi >= 0) {
          t_wsm   = cl_sm[bi].values[0];
          var_wsm = cl_sm[bi].sigmas[0] * cl_sm[bi].sigmas[0];
          wsm_ok  = true;
        }
      }

      const double t_truth_vtx = ideal.tTruthVtx;

      // Shared gate: idealised track times against whichever vertex time.
      auto smearedGate = [&](const std::vector<int>& base, double t_vtx,
                             double var_vtx, bool vtx_ok, double infl) {
        std::vector<int> out; out.reserve(base.size());
        const double den = std::sqrt(infl * infl * var_vtx
                                   + TRUTH_TRK_SMEAR * TRUTH_TRK_SMEAR);
        for (int idx : base) {
          auto it = smTimes.find(idx);
          if (!vtx_ok || it == smTimes.end() || den <= 0) { out.push_back(idx); continue; }
          if (std::abs(it->second - t_vtx) / den < GATE_SIGMA) out.push_back(idx);
        }
        return out;
      };

      // Truth row: ideal vertex time (10 ps) + the same ideal track times, so
      // it differs from the WAVeS-smeared row in exactly one variable.
      std::vector<int> trk_truth     = smearedGate(trk_all,     t_truth_vtx,
                                                   TRUTH_VTX_SMEAR * TRUTH_VTX_SMEAR, true, 1.0);
      std::vector<int> trk_truth_cen = smearedGate(trk_all_cen, t_truth_vtx,
                                                   TRUTH_VTX_SMEAR * TRUTH_VTX_SMEAR, true, 1.0);

      std::vector<int> trk_wsm     = smearedGate(trk_all,     t_wsm, var_wsm, wsm_ok, INFL.waves);
      std::vector<int> trk_wsm_cen = smearedGate(trk_all_cen, t_wsm, var_wsm, wsm_ok, INFL.waves);

      // Build per-scenario sets once per event for O(1) ghost-index lookup.
      struct TrackSets { std::unordered_set<int> all, hgtd, trkptz, waves, waves_ideal, truth, tzp; };
      TrackSets fwd{ {trk_all.begin(),    trk_all.end()},
                     {trk_hgtd.begin(),   trk_hgtd.end()},
                     {trk_trkptz.begin(), trk_trkptz.end()},
                     {trk_waves.begin(),  trk_waves.end()},
                     {trk_wsm.begin(),    trk_wsm.end()},
                     {trk_truth.begin(),  trk_truth.end()},
                     {trk_tzp.begin(),    trk_tzp.end()} };
      TrackSets cen{ {trk_all_cen.begin(),    trk_all_cen.end()},
                     {trk_hgtd_cen.begin(),   trk_hgtd_cen.end()},
                     {trk_trkptz_cen.begin(), trk_trkptz_cen.end()},
                     {trk_waves_cen.begin(),  trk_waves_cen.end()},
                     {trk_wsm_cen.begin(),    trk_wsm_cen.end()},
                     {trk_truth_cen.begin(),  trk_truth_cen.end()},
                     {trk_tzp_cen.begin(),    trk_tzp_cen.end()} };

      // ── Fill jets into pT slices. ─────────────────────────────────────────────
      // eta_min/eta_max select the acceptance; do_floor is false for the central
      // baseline so its jets cannot contaminate the forward untimed-floor
      // counters, which are reported as a forward quantity.
      auto fillJets = [&](std::vector<Scenario>& sv, double pt_lo, double pt_hi,
                          double eta_min, double eta_max, bool do_floor,
                          const TrackSets& S) {
        for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
          if (branch.isJetRemoved(j)) continue;  // lepton-overlap removed (Z+jets)
          double j_pt  = branch.topoJetPt[j];
          double j_eta = branch.topoJetEta[j];
          double j_phi = branch.topoJetPhi[j];
          if (j_pt <= pt_lo || j_pt >= pt_hi) continue;
          if (std::abs(j_eta) < eta_min || std::abs(j_eta) > eta_max) continue;
          bool isHS = paperIsHS(j_eta, j_phi);
          bool isPU = paperIsPU(j_eta, j_phi);
          if (!isHS && !isPU) continue;
          const auto& ghost = branch.topoJetGhostTrackIdx[j];
          auto fill = [&](Scenario& s, const std::unordered_set<int>& s_set, bool ok = true) {
            if (!ok) return;
            double r = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, s_set);
            if (isHS) s.h_hs->Fill(r);
            else      s.h_pu->Fill(r);
          };
          fill(sv[0], S.all);                         // ITk-only
          fill(sv[1], S.hgtd);                        // HGTD t0 (Athena)
          fill(sv[2], S.trkptz);                      // TRKPTZ t0
          fill(sv[3], S.waves);                       // WAVeS t0
          fill(sv[4], S.waves_ideal);                 // ORACLE: 30 ps tracks + perfect selection
          fill(sv[5], S.truth);                       // truth t0, 10 (+) 30 ps
          fill(sv[6], S.tzp);                         // TZP t0 (classical selector)

          // How the WAVeS gate moved this jet's R_pT relative to ITk-only.
          // Forward only (do_floor marks the forward calls): central is outside
          // HGTD acceptance, so its gate is a no-op and would only dilute this
          // with a wall of "unchanged".
          if (do_floor) {
            double r_z = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, S.all);
            double r_w = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, S.waves);
            int k = isHS ? 0 : 1;
            state.rpt_n[k]++;
            if      (r_w > r_z + 1e-9) state.rpt_up[k]++;
            else if (r_w > r_z - 1e-9) state.rpt_same[k]++;
            else {
              state.rpt_down[k]++;
              if (r_w <= 0.0) state.rpt_zeroed[k]++;
              if (r_z > 0.0)  state.rpt_relloss[k] += (r_z - r_w) / r_z;
            }
          }

          // Untimed floor accounting, per slice (forward only).
          if (!do_floor) continue;
          for (int idx : ghost) {
            if (!S.all.count(idx)) continue;
            double pt = branch.trackPt[idx];
            bool untimed = (branch.trackTimeValid[idx] != 1);
            if (pt_lo >= 40.0) {
              if (isHS) { state.hs_tot_pt += pt; if (untimed) state.hs_floor_pt += pt; }
              else      { state.pu_tot_pt += pt; if (untimed) state.pu_floor_pt += pt; }
            } else {
              if (isHS) { state.hs_tot_lo += pt; if (untimed) state.hs_floor_lo += pt; }
              else      { state.pu_tot_lo += pt; if (untimed) state.pu_floor_lo += pt; }
            }
          }
        }
      };

      fillJets(state.scen_lo, 30.0, 40.0, JET_ETA_MIN, JET_ETA_MAX, true,  fwd);
      fillJets(state.scen_hi, 40.0, 1e9,  JET_ETA_MIN, JET_ETA_MAX, true,  fwd);
      fillJets(state.scen_lo_cen, 30.0, 40.0, 0.0, CENTRAL_ETA_MAX, false, cen);
      fillJets(state.scen_hi_cen, 40.0, 1e9,  0.0, CENTRAL_ETA_MAX, false, cen);

      // Full ntuple file path + local (per-file) entry number for the
      // event-display commands printed after the loop. Same
      // TTreeProcessorMT-safe pattern as src/clustering_hist.cxx: no outer
      // TChain is available inside this lambda, and reader.GetTree() gives the
      // currently-loaded per-file constituent tree directly (a TChain-bound
      // sequential reader would return a chain-global entry number instead).
      std::string filePath   = reader.GetTree()->GetCurrentFile()->GetName();
      Long64_t    localEntry = reader.GetTree()->GetReadEntry();

      // ── Display candidates: track-time misassociation, and WAVeS-vs-TRKPTZ.
      if (PRINT_EVENT_DISPLAYS) {
        const double tTruth = branch.truthVtxTime[0];
        // Count forward HS tracks that are timed, and how many of those carry a
        // time inconsistent with truth -- the same |pull| >= TRUTH_PULL_CUT
        // test calcHSTimingPurity uses, but counting TRACKS rather than
        // pT-weighting, since what makes a display legible is how many wrong
        // bars are visible.
        int nTimedHS = 0, nMisTimed = 0;
        for (int idx : trk_z) {
          if (branch.trackToTruthvtx[idx] != 0 || branch.trackTimeValid[idx] != 1) continue;
          ++nTimedHS;
          const int    pi = branch.trackToParticle[idx];
          const double tTrue = (pi != -1) ? branch.particleT[pi]
                                          : branch.truthVtxTime[branch.trackToTruthvtx[idx]];
          if (std::abs(branch.trackTime[idx] - tTrue) / branch.trackTimeRes[idx] >= TRUTH_PULL_CUT)
            ++nMisTimed;
        }
        const float pur = calcHSTimingPurity(trk_z, &branch);

        // Misassociation: ranked by the NUMBER of mis-timed forward HS tracks.
        // Ranking on purity alone pinned every pick at the 8-track floor (a
        // sparse event reaches a low fraction easily); ranking on
        // (1 - purity) x multiplicity favoured busy events and never went below
        // ~0.2 purity. The count is the quantity that is actually low-purity
        // AND populated, which is what makes the wrong times visible on a slide.
        // The two criteria pull against each other: a low purity FRACTION is
        // easiest to reach in a sparse event, while a high mis-timed COUNT
        // favours busy ones with middling purity. So gate on genuinely low
        // purity, then rank by count within that -- the busiest of the
        // genuinely bad events, rather than the extreme of either alone.
        // The event must also FAIL. Gating on low purity alone selected events
        // where the tracks are badly mis-timed and the algorithm recovers the
        // right vertex time anyway -- which illustrates the opposite of the
        // point. failure_decomposition's timing-misassignment category requires
        // both the mis-timed tracks AND the resulting time being wrong, so the
        // delivered time must miss truth by more than the efficiency window.
        if (pur < 0.25f && waves_ok && std::abs(t_waves - tTruth) > PASS_SIGMA)
          insertEventCase(state.cases_mis,
                          {filePath, localEntry, t_waves, (double)nMisTimed,
                           pur, (double)nMisTimed, t_waves - tTruth});

        // WAVeS lands on truth where the TRKPTZ baseline does not -- the case
        // the score change exists for.
        if (waves_ok && trkptz_ok) {
          const double dW = std::abs(t_waves - tTruth), dT = std::abs(t_trkptz - tTruth);
          if (dW < 40.0 && dT > 120.0)
            insertEventCase(state.cases_wwin,
                            {filePath, localEntry, t_waves, dT - dW, dW, dT, 0.0});
        }
      }

      // ── VBS-topology regions ────────────────────────────────────────────────
      // Two topologies where forward timing is the deciding information, taken
      // off the SAME VBS candidate pair the clustering-side selection uses
      // (calcBestVbsPair: opposite-hemisphere, pT-passing, max m_jj):
      //
      //   R1  both legs forward, one truth-HS + one truth-PU.
      //       Both are in HGTD acceptance, so timing has to say WHICH is the
      //       hard-scatter jet. Fills h_hs from the HS leg, h_pu from the PU leg
      //       -- a self-contained ROC.
      //
      //   R2  one leg forward + truth-PU, the other central + truth-HS.
      //       Only the forward PU leg is filled (into h_pu): the HS leg is
      //       central, outside HGTD acceptance entirely, so it carries no timing
      //       information and would only dilute the signal side. R2's h_hs is
      //       therefore left EMPTY by construction -- rpt_v5_plot pairs
      //       r2.h_pu against r1.h_hs to build its ROC. Changing that here
      //       (e.g. to also fill the central HS leg) would silently make the R2
      //       ROC compare two different detector acceptances.
      //
      // Labels use paperIsHS/paperIsPU, matching the inclusive histograms above,
      // rather than BranchPointerWrapper::isJetTruthHS -- mixing the two
      // definitions across the same plot set would make the regions
      // incomparable to the inclusive result.
      //
      // The VBS topology knobs (--vbs-deta=, --vbs-mjj=) gate the REGIONS ONLY,
      // not the inclusive _lo/_hi histograms above. That asymmetry is
      // deliberate and load-bearing:
      //   - rpt_v5's inclusive measurement is by design selection-free ("every
      //     forward jet in the acceptance contributes an independent RpT
      //     measurement", see the file header), so gating it would silently
      //     redefine the primary result.
      //   - the regions ARE a VBS topology selection, so the knobs that define
      //     that topology have to apply here or they mean nothing.
      // Note calcBestVbsPair only FINDS the max-m_jj opposite-hemisphere pair;
      // it does not apply either cut. passJetPtCut is what normally applies
      // them on the clustering side, and rpt_v5_hist deliberately never calls
      // it -- so without the explicit test below the knobs would be inert here,
      // changing only the output filename via SELECTION_TAG and making two
      // different --vbs-mjj runs produce byte-identical physics.
      //
      // Region membership itself comes from BranchPointerWrapper::
      // classifyVbsRegion, which also owns the paper-HS/PU labelling that used
      // to live in this file's local lambdas.
      //
      // The eta window passed here is JET_ETA_MIN/JET_ETA_MAX, NOT the
      // clustering side's unbounded VBS_FWD_ETA_MAX, and that difference is
      // deliberate: the lines below fill an R_pT histogram FOR each leg, and a
      // leg past |eta| 4.00 is outside HGTD, so it would enter every timed
      // scenario carrying its ITk-only R_pT and dilute the R1 rejection with
      // jets timing never touched. Region membership on the clustering side is
      // a topology statement; here it has to be a timeable-jet statement.
      // (These two windows have never matched -- 2.4/3.8 against 2.38/4.00 --
      // despite an older comment on both sides claiming they did.)
      //
      // The narrow result is kept for the wide-window block below, which must
      // contain it (see RegionRow).
      // Pileup-jet tagging before pairing (--jvt; see JVT_SEL). keepPtr stays
      // nullptr without it, which is exactly the historical classification.
      std::vector<JetDisc> jetDisc;
      std::vector<char>    jvtKeep;
      std::vector<int>     rmJvtIdx, rmFjvtIdx;   // pT-passing jets removed
      if (JVT_ON) {
        jetDisc = computeJetDiscriminants(branch);
        jvtKeep.assign(jetDisc.size(), 1);
        for (int j = 0; j < (int)jetDisc.size(); ++j) {
          const bool rmJ = jvtRemoves(jetDisc[j], *JVT_SEL);
          const bool rmF = !rmJ && fjvtRemoves(jetDisc[j], *JVT_SEL);
          if (!rmJ && !rmF) continue;
          jvtKeep[j] = 0;
          if (branch.topoJetPt[j] > MIN_JET_PT && !branch.isJetRemoved(j))
            (rmJ ? rmJvtIdx : rmFjvtIdx).push_back(j);
        }
      }
      const std::vector<char>* keepPtr = JVT_ON ? &jvtKeep : nullptr;

      VbsRegion narrowRegion = VbsRegion::NONE;
      int narrowHS = -1, narrowPU = -1;
      {
        int fwdHS = -1, fwdPU = -1;
        auto region = branch.classifyVbsRegion(JET_ETA_MIN, JET_ETA_MAX,
                                               CENTRAL_ETA_MAX, &fwdHS, &fwdPU,
                                               nullptr, keepPtr);
        narrowRegion = region;
        narrowHS = fwdHS;
        narrowPU = fwdPU;

        if (region != VbsRegion::NONE) {
          // Fill one jet into a region's scenario set, as HS or PU.
          auto fillRegion = [&](std::vector<Scenario>& sv, int j, bool asHS) {
            double j_pt  = branch.topoJetPt[j];
            double j_eta = branch.topoJetEta[j];
            double j_phi = branch.topoJetPhi[j];
            const auto& ghost = branch.topoJetGhostTrackIdx[j];
            auto put = [&](Scenario& s, const std::unordered_set<int>& s_set) {
              double r = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, s_set);
              (asHS ? s.h_hs : s.h_pu)->Fill(r);
            };
            put(sv[0], fwd.all);
            put(sv[1], fwd.hgtd);
            put(sv[2], fwd.trkptz);
            put(sv[3], fwd.waves);
            put(sv[4], fwd.waves_ideal);
            put(sv[5], fwd.truth);
            put(sv[6], fwd.tzp);
          };

          // Per-jet RpT under the no-timing baseline and under WAVeS, used both
          // to fill and to rank this event as an display candidate below.
          auto rptOf = [&](int j, const std::unordered_set<int>& s_set) {
            return computeRpT(&branch, branch.topoJetGhostTrackIdx[j],
                              branch.topoJetPt[j], branch.topoJetEta[j],
                              branch.topoJetPhi[j], s_set);
          };

          if (region == VbsRegion::R1) {
            fillRegion(state.scen_r1, fwdHS, true);
            fillRegion(state.scen_r1, fwdPU, false);

            // Rank by how much timing changes the HS-vs-PU RpT margin. The
            // sign matters (positive = timing widened the correct gap,
            // negative = timing eroded or inverted it), so rank on |delta| and
            // let the printed line say which -- one list surfaces both the
            // rescues and the regressions rather than needing two.
            double mZ = rptOf(fwdHS, fwd.all)   - rptOf(fwdPU, fwd.all);
            double mW = rptOf(fwdHS, fwd.waves) - rptOf(fwdPU, fwd.waves);
            insertRegionCase(state.cases_r1,
                             {filePath, localEntry, fwdHS, fwdPU,
                              branch.topoJetPt[fwdHS], branch.topoJetEta[fwdHS],
                              branch.topoJetPt[fwdPU], branch.topoJetEta[fwdPU],
                              mZ, mW, mW - mZ, t_waves});
          } else {  // R2 -- forward PU leg only; no forward HS leg exists.
            fillRegion(state.scen_r2, fwdPU, false);

            // Rank by how far timing pushes the fake's RpT down: a forward PU
            // jet with high no-timing RpT is precisely the one that fakes a
            // tagging jet, and the drop is the rejection actually delivered.
            double rZ = rptOf(fwdPU, fwd.all);
            double rW = rptOf(fwdPU, fwd.waves);
            insertRegionCase(state.cases_r2,
                             {filePath, localEntry, -1, fwdPU,
                              0.0, 0.0,
                              branch.topoJetPt[fwdPU], branch.topoJetEta[fwdPU],
                              rZ, rW, rW - rZ, t_waves});
            // R2's SUCCESS list above is one-sided by construction:
            // applyTimeGate returns a subset, so rW <= rZ always and every
            // entry is "suppressed". It therefore cannot show a failure.
            // R2 fails when a fake SURVIVES -- high RpT still standing after
            // the gate -- which has |rW - rZ| ~ 0 and sits at the bottom of
            // that ranking. Rank those by the surviving RpT instead.
            insertRegionCase(state.cases_r2_fail,
                             {filePath, localEntry, -1, fwdPU,
                              0.0, 0.0,
                              branch.topoJetPt[fwdPU], branch.topoJetEta[fwdPU],
                              rZ, rW, rW, t_waves});
          }
        }
      }

      // ── Wide-window VBS regions -> side-file tree (see RegionRow) ─────────
      // Same classifier, same max-m_jj pair, same paper labels as the block
      // above; only the forward upper edge is lifted to VBS_FWD_ETA_MAX, so
      // legs past |eta| 3.8 are kept. Nothing here fills a histogram -- the
      // narrow _r1/_r2 sets above are exactly as before.
      {
        int wHS = -1, wPU = -1;
        BranchPointerWrapper::VbsPair wPair;
        const VbsRegion wRegion = branch.classifyVbsRegion(
            JET_ETA_MIN, VBS_FWD_ETA_MAX, CENTRAL_ETA_MAX, &wHS, &wPU, &wPair, keepPtr);

        // Lifting an upper edge can only add events: a narrow-region event
        // must be the same region with the same legs in the wide window.
        const bool core = (narrowRegion != VbsRegion::NONE);
        if (core && (narrowRegion != wRegion || narrowHS != wHS || narrowPU != wPU))
          ++state.n_narrow_not_subset;

        if (wRegion != VbsRegion::NONE) {
          const bool isR1 = (wRegion == VbsRegion::R1);
          // Index-aligned with makeScenarios() / SCEN_NAMES, and with the
          // sets fillRegion above fills sv[0..6] from.
          const std::unordered_set<int>* sets[N_SCEN] = {
            &fwd.all, &fwd.hgtd, &fwd.trkptz, &fwd.waves,
            &fwd.waves_ideal, &fwd.truth, &fwd.tzp};
          // What each scenario gated against: t0, sigma_t0 BEFORE inflation,
          // whether a gate applied at all, and the inflation used -- exactly
          // the arguments of applyTimeGate / smearedGate above, so a display
          // can reproduce the gate without re-deriving anything.
          const double t0s[N_SCEN]   = {0.0, t_hgtd, t_trkptz, t_waves,
                                        t_wsm, t_truth_vtx, t_tzp};
          const double sig0s[N_SCEN] = {0.0, std::sqrt(var_hgtd), std::sqrt(var_trkptz),
                                        std::sqrt(var_waves), std::sqrt(var_wsm),
                                        TRUTH_VTX_SMEAR, std::sqrt(var_tzp)};
          const bool   oks[N_SCEN]   = {false, hgtd_vtx_valid, trkptz_ok, waves_ok,
                                        wsm_ok, true, tzp_ok};
          const double infls[N_SCEN] = {0.0, INFL.hgtd, INFL.trkptz, INFL.waves,
                                        INFL.waves, 1.0, INFL.tzp};

          // z-associated ghost tracks in the R_pT cone, and the timed subset.
          auto countTracks = [&](int j, int& ntrk, int& ntimed) {
            ntrk = ntimed = 0;
            for (int idx : branch.topoJetGhostTrackIdx[j]) {
              if (!fwd.all.count(idx)) continue;
              if (dR(branch.topoJetEta[j], branch.topoJetPhi[j],
                     branch.trackEta[idx], branch.trackPhi[idx]) > RPT_TRACK_JET_DR) continue;
              ++ntrk;
              if (branch.trackTimeValid[idx] == 1) ++ntimed;
            }
          };
          auto rptOf = [&](int j, const std::unordered_set<int>& s_set) {
            return computeRpT(&branch, branch.topoJetGhostTrackIdx[j],
                              branch.topoJetPt[j], branch.topoJetEta[j],
                              branch.topoJetPhi[j], s_set);
          };

          RegionRow row;
          row.file_path = filePath;
          row.entry     = localEntry;
          row.region    = isR1 ? 1 : 2;
          row.core      = core;
          row.idx_pu    = wPU;
          // R2 has no forward HS leg; its HS leg is the pair's other (central) jet.
          row.idx_hs    = isR1 ? wHS : (wPair.idxI == wPU ? wPair.idxJ : wPair.idxI);
          row.mjj       = wPair.mjj;
          row.deta      = wPair.dEta;
          row.hs_pt  = branch.topoJetPt[row.idx_hs];
          row.hs_eta = branch.topoJetEta[row.idx_hs];
          row.hs_phi = branch.topoJetPhi[row.idx_hs];
          row.pu_pt  = branch.topoJetPt[wPU];
          row.pu_eta = branch.topoJetEta[wPU];
          row.pu_phi = branch.topoJetPhi[wPU];
          countTracks(wPU, row.pu_ntrk, row.pu_ntimed);
          if (isR1) countTracks(wHS, row.hs_ntrk, row.hs_ntimed);
          else      row.hs_ntrk = row.hs_ntimed = -1;  // central leg: never timed
          {
            std::vector<int> passPtIdx;
            int nPt = 0, nPtEta = 0;
            branch.collectPtPassingJets(passPtIdx, nPt, nPtEta);
            if (keepPtr) {        // count only the jets the tagger kept
              nPtEta = 0;
              for (int j : passPtIdx) {
                const double ae = std::abs((double)branch.topoJetEta[j]);
                if ((*keepPtr)[j] && ae > MIN_ABS_ETA_JET && ae < MAX_ABS_ETA_JET) ++nPtEta;
              }
            }
            row.n_jets_fwd_acc = nPtEta;
          }
          for (int k = 0; k < N_SCEN; ++k) {
            row.rpt_pu[k] = rptOf(wPU, *sets[k]);
            row.rpt_hs[k] = isR1 ? rptOf(wHS, *sets[k]) : -1.0;
            row.t0[k]     = t0s[k];
            row.sig0[k]   = sig0s[k];
            row.ok[k]     = oks[k];
            row.infl[k]   = infls[k];
          }
          row.t_truth = branch.truthVtxTime[0];

          // Why the gate moved pT (see RegionRow). Same cone as computeRpT.
          {
            const double tHS = branch.truthVtxTime[0];
            auto cone = [&](int j, auto&& fn) {
              for (int idx : branch.topoJetGhostTrackIdx[j]) {
                if (!fwd.all.count(idx)) continue;
                if (dR(branch.topoJetEta[j], branch.topoJetPhi[j],
                       branch.trackEta[idx], branch.trackPhi[idx]) > RPT_TRACK_JET_DR) continue;
                fn(idx);
              }
            };
            if (isR1)
              cone(wHS, [&](int idx) {
                if (branch.trackToTruthvtx[idx] == 0) row.hs_hspt += branch.trackPt[idx]; });
            cone(wPU, [&](int idx) {
              if (branch.trackToTruthvtx[idx] == 0) return;
              row.pu_pupt += branch.trackPt[idx];
              if (branch.trackTimeValid[idx] != 1) row.pu_keep_untimed += branch.trackPt[idx];
            });
            for (int w = 0; w < N_WHY; ++w) {
              const int k = WHY_SCEN[w];
              const std::unordered_set<int>& kept = *sets[k];
              const double den0  = infls[k] * infls[k] * sig0s[k] * sig0s[k];
              const bool   farT0 = std::abs(t0s[k] - tHS) >= PASS_SIGMA;
              // The counterfactual: this case's gate with t0 moved to t_HS.
              auto keptWithTrueT0 = [&](int idx) {
                const double st = branch.trackTimeRes[idx];
                return std::abs(branch.trackTime[idx] - tHS) / std::sqrt(den0 + st * st) < GATE_SIGMA;
              };
              if (isR1)
                cone(wHS, [&](int idx) {
                  if (branch.trackToTruthvtx[idx] != 0 || kept.count(idx)) return;
                  const double pt = branch.trackPt[idx];   // an HS track the gate removed
                  if (!keptWithTrueT0(idx)) row.hs_rm_trk[w] += pt;
                  else (farT0 ? row.hs_rm_t0far[w] : row.hs_rm_t0near[w]) += pt;
                });
              cone(wPU, [&](int idx) {
                if (branch.trackToTruthvtx[idx] == 0 || !kept.count(idx)) return;
                if (branch.trackTimeValid[idx] != 1) return;   // pu_keep_untimed
                const double pt = branch.trackPt[idx];     // a timed PU track that survived
                if (!oks[k])                  row.pu_keep_not0[w]   += pt;
                else if (keptWithTrueT0(idx)) row.pu_keep_intime[w] += pt;
                else (farT0 ? row.pu_keep_t0far[w] : row.pu_keep_t0near[w]) += pt;
              });
            }
          }
          if (JVT_ON) {
            row.pu_jvt_rpt = jetDisc[wPU].rpt;
            row.pu_fjvt    = jetDisc[wPU].fjvt;
            row.hs_jvt_rpt = jetDisc[row.idx_hs].rpt;
            row.hs_fjvt    = jetDisc[row.idx_hs].fjvt;
            row.rm_jvt     = rmJvtIdx;
            row.rm_fjvt    = rmFjvtIdx;
          }
          state.region_rows.push_back(std::move(row));
        }
      }

    }
  });
  std::cout << "\n";

  // --- Merge per-thread state into one ---
  if (stateRegistry.empty()) {
    std::cerr << "No events processed.  Aborting.\n";
    return 1;
  }
  ThreadState& merged = *stateRegistry.front();
  for (size_t i = 1; i < stateRegistry.size(); ++i) {
    ThreadState& other = *stateRegistry[i];
    for (size_t k = 0; k < merged.scen_lo.size(); ++k) {
      merged.scen_lo[k].h_hs->Add(other.scen_lo[k].h_hs);
      merged.scen_lo[k].h_pu->Add(other.scen_lo[k].h_pu);
    }
    for (size_t k = 0; k < merged.scen_r1.size(); ++k) {
      merged.scen_r1[k].h_hs->Add(other.scen_r1[k].h_hs);
      merged.scen_r1[k].h_pu->Add(other.scen_r1[k].h_pu);
    }
    for (size_t k = 0; k < merged.scen_r2.size(); ++k) {
      merged.scen_r2[k].h_hs->Add(other.scen_r2[k].h_hs);
      merged.scen_r2[k].h_pu->Add(other.scen_r2[k].h_pu);
    }
    for (size_t k = 0; k < merged.scen_hi.size(); ++k) {
      merged.scen_hi[k].h_hs->Add(other.scen_hi[k].h_hs);
      merged.scen_hi[k].h_pu->Add(other.scen_hi[k].h_pu);
    }
    for (int k = 0; k < 2; ++k) {
      merged.rpt_n[k]      += other.rpt_n[k];
      merged.rpt_up[k]     += other.rpt_up[k];
      merged.rpt_same[k]   += other.rpt_same[k];
      merged.rpt_down[k]   += other.rpt_down[k];
      merged.rpt_zeroed[k] += other.rpt_zeroed[k];
      merged.rpt_relloss[k]+= other.rpt_relloss[k];
    }
    for (size_t k = 0; k < merged.scen_lo_cen.size(); ++k) {
      merged.scen_lo_cen[k].h_hs->Add(other.scen_lo_cen[k].h_hs);
      merged.scen_lo_cen[k].h_pu->Add(other.scen_lo_cen[k].h_pu);
    }
    for (size_t k = 0; k < merged.scen_hi_cen.size(); ++k) {
      merged.scen_hi_cen[k].h_hs->Add(other.scen_hi_cen[k].h_hs);
      merged.scen_hi_cen[k].h_pu->Add(other.scen_hi_cen[k].h_pu);
    }
    merged.n_total      += other.n_total;
    merged.n_pass_basic += other.n_pass_basic;
    merged.n_hgtd_valid += other.n_hgtd_valid;
    merged.n_pass_lepton_sel += other.n_pass_lepton_sel;
    merged.n_rej_no_lepton   += other.n_rej_no_lepton;
    merged.n_rej_one_lepton  += other.n_rej_one_lepton;
    merged.n_rej_no_ossf_pair += other.n_rej_no_ossf_pair;
    merged.pu_tot_pt    += other.pu_tot_pt;
    merged.pu_floor_pt  += other.pu_floor_pt;
    merged.hs_tot_pt    += other.hs_tot_pt;
    merged.hs_floor_pt  += other.hs_floor_pt;
    merged.pu_tot_lo    += other.pu_tot_lo;
    merged.pu_floor_lo  += other.pu_floor_lo;
    merged.hs_tot_lo    += other.hs_tot_lo;
    merged.hs_floor_lo  += other.hs_floor_lo;

    // Event-display diagnostic candidates: merge each category's top-N
    // (see mergeCases doc comment near the top of this file).
    merged.pull_dt2_hgtd   += other.pull_dt2_hgtd;   merged.pull_var_hgtd   += other.pull_var_hgtd;   merged.pull_n_hgtd   += other.pull_n_hgtd;  merged.pull_dt_hgtd += other.pull_dt_hgtd;  merged.pull_ntail_hgtd += other.pull_ntail_hgtd;
    merged.pull_dt2_trkptz += other.pull_dt2_trkptz; merged.pull_var_trkptz += other.pull_var_trkptz; merged.pull_n_trkptz += other.pull_n_trkptz;  merged.pull_dt_trkptz += other.pull_dt_trkptz;  merged.pull_ntail_trkptz += other.pull_ntail_trkptz;
    merged.pull_dt2_waves  += other.pull_dt2_waves;  merged.pull_var_waves  += other.pull_var_waves;  merged.pull_n_waves  += other.pull_n_waves;  merged.pull_dt_waves += other.pull_dt_waves;  merged.pull_ntail_waves += other.pull_ntail_waves;
    merged.pull_dt2_tzp    += other.pull_dt2_tzp;    merged.pull_var_tzp    += other.pull_var_tzp;    merged.pull_n_tzp    += other.pull_n_tzp;    merged.pull_dt_tzp   += other.pull_dt_tzp;   merged.pull_ntail_tzp   += other.pull_ntail_tzp;
    mergeRegionCases(merged.cases_r1, other.cases_r1);
    mergeRegionCases(merged.cases_r2, other.cases_r2);
    mergeRegionCases(merged.cases_r2_fail, other.cases_r2_fail);
    mergeEventCases (merged.cases_mis,  other.cases_mis);
    mergeEventCases (merged.cases_wwin, other.cases_wwin);
    merged.region_rows.insert(merged.region_rows.end(),
                              std::make_move_iterator(other.region_rows.begin()),
                              std::make_move_iterator(other.region_rows.end()));
    merged.n_narrow_not_subset += other.n_narrow_not_subset;
  }

  std::cout << "\nFINISHED PROCESSING\n";
  phase.mark("event loop done");

  // --- Save every histogram + scalar accumulator to a ROOT file ---
  const std::string histPath = MyUtl::histFilePath("rpt_v5_hist.root");
  MyUtl::HistWriter writer(histPath);
  saveScenarios(writer, merged.scen_lo);
  saveScenarios(writer, merged.scen_hi);
  saveScenarios(writer, merged.scen_r1);
  saveScenarios(writer, merged.scen_r2);
  saveScenarios(writer, merged.scen_lo_cen);
  saveScenarios(writer, merged.scen_hi_cen);
  writer.WriteScalar("meta_n_total",      static_cast<Long64_t>(merged.n_total));
  writer.WriteScalar("meta_n_pass_basic", static_cast<Long64_t>(merged.n_pass_basic));
  writer.WriteScalar("meta_n_hgtd_valid", static_cast<Long64_t>(merged.n_hgtd_valid));
  writer.WriteScalar("meta_n_pass_lepton_sel", static_cast<Long64_t>(merged.n_pass_lepton_sel));
  writer.WriteScalar("meta_n_rej_no_lepton",    static_cast<Long64_t>(merged.n_rej_no_lepton));
  writer.WriteScalar("meta_n_rej_one_lepton",   static_cast<Long64_t>(merged.n_rej_one_lepton));
  writer.WriteScalar("meta_n_rej_no_ossf_pair", static_cast<Long64_t>(merged.n_rej_no_ossf_pair));
  writer.WriteScalar("meta_pu_tot_pt",    merged.pu_tot_pt);
  writer.WriteScalar("meta_pu_floor_pt",  merged.pu_floor_pt);
  writer.WriteScalar("meta_hs_tot_pt",    merged.hs_tot_pt);
  writer.WriteScalar("meta_hs_floor_pt",  merged.hs_floor_pt);
  writer.WriteScalar("meta_pu_tot_lo",    merged.pu_tot_lo);
  writer.WriteScalar("meta_pu_floor_lo",  merged.pu_floor_lo);
  writer.WriteScalar("meta_hs_tot_lo",    merged.hs_tot_lo);
  writer.WriteScalar("meta_hs_floor_lo",  merged.hs_floor_lo);
  writer.WriteRunMeta(MyUtl::ENERGY_LABEL, merged.n_total, MyUtl::VBS_JET_D_ETA, MyUtl::VBS_JET_MJJ);
  // meta_vbs_ prefix: hist_merge requires it to agree across shards rather
  // than summing it. Written only for a tagged run, so an untagged file's key
  // set is exactly what it was before --jvt existed.
  if (JVT_ON)
    writer.WriteScalar("meta_vbs_jvt_wp", static_cast<Long64_t>(JVT_SEL - JVT_WPS));
  writer.Close();
  std::cout << "Wrote histograms to " << histPath << "\n";
  phase.mark("histograms written");

  // --- Wide-window region tree, side file (see RegionRow) ---
  //     Always written, even empty, so the condor template can always name it.
  //     Sorted so the file does not depend on how TTreeProcessorMT scheduled
  //     entries across threads.
  {
    auto& rows = merged.region_rows;
    std::sort(rows.begin(), rows.end(), [](const RegionRow& a, const RegionRow& b) {
      return a.file_path != b.file_path ? a.file_path < b.file_path : a.entry < b.entry;
    });
    const std::string regPath = MyUtl::histFilePath("rpt_v5_regions.root");
    TFile fout(regPath.c_str(), "RECREATE");
    if (fout.IsZombie()) {
      std::cerr << "Could not open " << regPath << " for writing.\n";
      return 1;
    }
    TTree tree("regions", "rpt_v5 wide-window VBS-region events "
                          "(forward |eta| > 2.4, no upper edge; one row per event)");
    RegionRow r;
    // Run constants, repeated per row so the tree is self-describing and
    // survives hadd (which would keep only one copy of a TParameter).
    double gate   = GATE_SIGMA;
    bool   dzpara = RPT_USE_DZ_PARA;
    double cutMjj = MyUtl::VBS_JET_MJJ, cutDeta = MyUtl::VBS_JET_D_ETA;
    tree.Branch("file_path", &r.file_path);
    tree.Branch("entry",  &r.entry,  "entry/L");
    tree.Branch("region", &r.region, "region/I");
    tree.Branch("core",   &r.core,   "core/O");
    tree.Branch("idx_hs", &r.idx_hs, "idx_hs/I");
    tree.Branch("idx_pu", &r.idx_pu, "idx_pu/I");
    tree.Branch("mjj",    &r.mjj,    "mjj/D");
    tree.Branch("deta",   &r.deta,   "deta/D");
    tree.Branch("hs_pt",  &r.hs_pt,  "hs_pt/D");
    tree.Branch("hs_eta", &r.hs_eta, "hs_eta/D");
    tree.Branch("hs_phi", &r.hs_phi, "hs_phi/D");
    tree.Branch("pu_pt",  &r.pu_pt,  "pu_pt/D");
    tree.Branch("pu_eta", &r.pu_eta, "pu_eta/D");
    tree.Branch("pu_phi", &r.pu_phi, "pu_phi/D");
    tree.Branch("hs_ntrk",   &r.hs_ntrk,   "hs_ntrk/I");
    tree.Branch("hs_ntimed", &r.hs_ntimed, "hs_ntimed/I");
    tree.Branch("pu_ntrk",   &r.pu_ntrk,   "pu_ntrk/I");
    tree.Branch("pu_ntimed", &r.pu_ntimed, "pu_ntimed/I");
    tree.Branch("n_jets_fwd_acc", &r.n_jets_fwd_acc, "n_jets_fwd_acc/I");
    for (int k = 0; k < N_SCEN; ++k) {
      const std::string s = SCEN_NAMES[k];
      tree.Branch(("rpt_hs_" + s).c_str(), &r.rpt_hs[k], ("rpt_hs_" + s + "/D").c_str());
      tree.Branch(("rpt_pu_" + s).c_str(), &r.rpt_pu[k], ("rpt_pu_" + s + "/D").c_str());
      tree.Branch(("t0_"     + s).c_str(), &r.t0[k],     ("t0_"     + s + "/D").c_str());
      tree.Branch(("sig0_"   + s).c_str(), &r.sig0[k],   ("sig0_"   + s + "/D").c_str());
      tree.Branch(("infl_"   + s).c_str(), &r.infl[k],   ("infl_"   + s + "/D").c_str());
      tree.Branch(("ok_"     + s).c_str(), &r.ok[k],     ("ok_"     + s + "/O").c_str());
    }
    tree.Branch("t_truth",    &r.t_truth, "t_truth/D");
    tree.Branch("gate_sigma", &gate,      "gate_sigma/D");
    tree.Branch("rpt_dzpara", &dzpara,    "rpt_dzpara/O");
    tree.Branch("vbs_mjj_cut",  &cutMjj,  "vbs_mjj_cut/D");
    tree.Branch("vbs_deta_cut", &cutDeta, "vbs_deta_cut/D");
    // --jvt: 0 none, 1 loose, 2 tight (JVT_WPS order), and what it did.
    int jvtWp = (int)(JVT_SEL - JVT_WPS);
    tree.Branch("jvt_wp",     &jvtWp,        "jvt_wp/I");
    tree.Branch("hs_jvt_rpt", &r.hs_jvt_rpt, "hs_jvt_rpt/D");
    tree.Branch("hs_fjvt",    &r.hs_fjvt,    "hs_fjvt/D");
    tree.Branch("pu_jvt_rpt", &r.pu_jvt_rpt, "pu_jvt_rpt/D");
    tree.Branch("pu_fjvt",    &r.pu_fjvt,    "pu_fjvt/D");
    tree.Branch("rm_jvt",     &r.rm_jvt);
    tree.Branch("rm_fjvt",    &r.rm_fjvt);
    tree.Branch("hs_hspt",         &r.hs_hspt,         "hs_hspt/D");
    tree.Branch("pu_pupt",         &r.pu_pupt,         "pu_pupt/D");
    tree.Branch("pu_keep_untimed", &r.pu_keep_untimed, "pu_keep_untimed/D");
    for (int w = 0; w < N_WHY; ++w) {
      const std::string s = SCEN_NAMES[WHY_SCEN[w]];
      auto B = [&](const char* stem, double* v) {
        const std::string n = std::string(stem) + "_" + s;
        tree.Branch(n.c_str(), v, (n + "/D").c_str());
      };
      B("hs_rm_t0far", &r.hs_rm_t0far[w]);     B("hs_rm_t0near", &r.hs_rm_t0near[w]);
      B("hs_rm_trk",   &r.hs_rm_trk[w]);
      B("pu_keep_not0", &r.pu_keep_not0[w]);   B("pu_keep_intime", &r.pu_keep_intime[w]);
      B("pu_keep_t0far", &r.pu_keep_t0far[w]); B("pu_keep_t0near", &r.pu_keep_t0near[w]);
    }
    long nR1 = 0, nR2 = 0, nR1core = 0, nR2core = 0;
    for (const auto& row : rows) {
      r = row;
      tree.Fill();
      (row.region == 1 ? nR1 : nR2)++;
      if (row.core) (row.region == 1 ? nR1core : nR2core)++;
    }
    tree.Write();
    fout.Close();
    std::printf("\n=== WIDE-WINDOW VBS REGIONS (forward |eta| > %.1f, no upper edge) ===\n",
                JET_ETA_MIN);
    std::printf("  R1 events : %8ld   (of which in the narrow 2.4-3.8 region: %ld)\n", nR1, nR1core);
    std::printf("  R2 events : %8ld   (of which in the narrow 2.4-3.8 region: %ld)\n", nR2, nR2core);
    std::printf("  narrow events NOT contained in the wide region: %ld  (must be 0)\n",
                merged.n_narrow_not_subset);
    std::cout << "Wrote region tree to " << regPath << "\n";
  }
  phase.mark("region tree written");

  // --- Z+jets event-selection breakdown (no-op elsewhere: n_pass_lepton_sel
  //     == n_pass_basic when OVERLAP_REMOVAL is unset). Printed directly here
  //     rather than added to rpt_v5_plot's console summary, since that reads
  //     mandatory scalars from the hist file and would break on older files
  //     that predate meta_n_pass_lepton_sel. ---
  std::cout << "\n=== Z+JETS EVENT SELECTION (of " << merged.n_pass_basic
            << " passing vertex quality) ===\n";
  std::cout << "  Pass Z->ll lepton selection : " << merged.n_pass_lepton_sel
            << " (" << std::fixed << std::setprecision(1)
            << (merged.n_pass_basic > 0
                  ? 100.0 * merged.n_pass_lepton_sel / merged.n_pass_basic : 0.0)
            << "%)\n";
  std::cout << "    Rejected, 0 good leptons  : " << merged.n_rej_no_lepton    << '\n';
  std::cout << "    Rejected, 1 good lepton   : " << merged.n_rej_one_lepton   << '\n';
  std::cout << "    Rejected, no OS-SF pair   : " << merged.n_rej_no_ossf_pair << '\n';

  // --- Per-scenario vertex-time calibration -----------------------------------
  //     sigma ratio = observed core spread of (t_trk - t_vtx) for truth-HS
  //     tracks, divided by the QUOTED sqrt(var_vtx + var_trk). A value above 1
  //     means the quoted error is understated by that factor, so a nominal
  //     GATE_SIGMA cut behaves like GATE_SIGMA/ratio. Set INFL_* to the ratio.
  if (PRINT_PULL_DIAG) {
    std::printf("\n=== VERTEX-TIME CALIBRATION (truth-HS tracks, |dt| < 150 ps core) ===\n");
    std::printf("  %-8s %9s %9s %10s %10s %8s %8s %7s\n",
                "scenario", "n core", "mean dt", "core sig", "quoted", "ratio", "tail>150", "in use");
    auto row = [](const char* nm, double dt2, double dt1, double var, long n,
                  long ntail, double inuse) {
      if (n < 100) { std::printf("  %-8s %9ld   (too few)\n", nm, n); return; }
      double mean = dt1 / n;
      // Width ABOUT THE MEAN: a systematic offset is a bias, not resolution, and
      // scaling sigma_vtx cannot correct it -- so it must not be folded into the
      // inflation the way an RMS about zero would fold it.
      double sig  = std::sqrt(std::max(0.0, dt2 / n - mean * mean));
      double quo  = std::sqrt(var / n);
      double tail = 100.0 * ntail / double(n + ntail);
      std::printf("  %-8s %9ld %7.1fps %8.1fps %8.1fps %8.2f %7.1f%% %7.2f\n",
                  nm, n, mean, sig, quo, quo > 0 ? sig / quo : 0.0, tail, inuse);
    };
    row("hgtd",   merged.pull_dt2_hgtd,   merged.pull_dt_hgtd,   merged.pull_var_hgtd,
        merged.pull_n_hgtd,   merged.pull_ntail_hgtd,   INFL.hgtd);
    row("trkptz", merged.pull_dt2_trkptz, merged.pull_dt_trkptz, merged.pull_var_trkptz,
        merged.pull_n_trkptz, merged.pull_ntail_trkptz, INFL.trkptz);
    row("waves",  merged.pull_dt2_waves,  merged.pull_dt_waves,  merged.pull_var_waves,
        merged.pull_n_waves,  merged.pull_ntail_waves,  INFL.waves);
    row("tzp",    merged.pull_dt2_tzp,    merged.pull_dt_tzp,    merged.pull_var_tzp,
        merged.pull_n_tzp,    merged.pull_ntail_tzp,    INFL.tzp);
    std::printf("  mean dt : systematic offset -- an inflation CANNOT correct this.\n");
    std::printf("  tail    : fraction outside the core window, i.e. how much\n");
    std::printf("            structure the core-width calibration does not see.\n");
    std::printf("  (ratio != in use -> update this sample's row in inflationFor() and rebuild)\n");
  }

  {
    std::printf("\n=== PER-JET EFFECT OF THE WAVeS TIME GATE (forward, vs ITk-only) ===\n");
    std::printf("  %-4s %9s %7s %9s %9s %9s %12s\n",
                "jet", "n", "raised", "unchanged", "lowered", "->zero", "mean loss*");
    const char* nm[2] = {"HS", "PU"};
    for (int k = 0; k < 2; ++k) {
      long n = merged.rpt_n[k];
      if (n == 0) continue;
      std::printf("  %-4s %9ld %7ld %8.1f%% %8.1f%% %8.1f%% %11.1f%%\n",
                  nm[k], n, merged.rpt_up[k],
                  100.0 * merged.rpt_same[k]   / n,
                  100.0 * merged.rpt_down[k]   / n,
                  100.0 * merged.rpt_zeroed[k] / n,
                  merged.rpt_down[k] ? 100.0 * merged.rpt_relloss[k] / merged.rpt_down[k] : 0.0);
    }
    std::printf("  *mean fractional R_pT loss, averaged over the LOWERED jets only.\n");
    std::printf("  'raised' must be 0 by construction (the gate returns a subset).\n");
  }

  // --- Event displays, R1/R2 ONLY -------------------------------------------
  //     The older WAVeS-vs-HGTD jet comparison and "timing-hurt HS jets"
  //     listings were removed rather than kept alongside these: they ranked on
  //     criteria unrelated to VBS topology, so their events were not
  //     necessarily in any region, and interleaving both sets on stdout made
  //     it impossible to tell which commands belonged to which study. The
  //     regions are the focus, so the display output is theirs alone.
  if (PRINT_EVENT_DISPLAYS) {
    auto printEvents = [](const char* title, const char* desc,
                          const std::vector<EventCase>& cases, const char* fmt) {
      std::cout << "\n=== " << title << " ===\n  " << desc << "\n\n";
      if (cases.empty()) { std::cout << "  (none found)\n"; return; }
      for (const auto& c : cases) {
        std::printf(fmt, c.v1, c.v2, c.v3);
        std::printf("  cd python && python3 event_display.py --file_path \"%s\""
                    " --event_num %lld --extra_time %.2f\n\n",
                    c.file_path.c_str(), c.entry, c.t_show);
      }
    };
    printEvents("TRACK-TIME MISASSOCIATION",
                "Forward HS tracks whose HGTD times disagree with truth: the "
                "times themselves are wrong, not merely imprecise.",
                merged.cases_mis,
                "  HS timing purity=%.2f  --  %.0f mis-timed HS tracks  --  delivered time off truth by %+.1f ps\n");
    printEvents("WAVeS BEATS TRKPTZ",
                "WAVeS lands on the truth time while the TRKPTZ baseline does not.",
                merged.cases_wwin,
                "  |dt|: WAVeS=%.1f ps  TRKPTZ=%.1f ps%.0s\n");

    auto printRegion = [](const char* title, const char* metric_desc,
                          const std::vector<RegionCase>& cases, bool isR1) {
      std::cout << "\n=== " << title << " ===\n";
      std::cout << "  " << metric_desc << "\n\n";
      if (cases.empty()) { std::cout << "  (none found)\n"; return; }
      for (const auto& c : cases) {
        if (isR1) {
          std::printf("  HS leg pT=%.1f eta=%+.2f | PU leg pT=%.1f eta=%+.2f"
                      "  margin: %.3f -> %.3f  (%+.3f, timing %s)\n",
                      c.hs_pt, c.hs_eta, c.pu_pt, c.pu_eta,
                      c.val_zonly, c.val_waves, c.metric,
                      c.metric > 0 ? "HELPED" : "HURT");
        } else {
          std::printf("  fwd PU leg pT=%.1f eta=%+.2f  RpT: %.3f -> %.3f"
                      "  (%+.3f, timing %s)\n",
                      c.pu_pt, c.pu_eta, c.val_zonly, c.val_waves, c.metric,
                      c.metric < 0 ? "SUPPRESSED the fake" : "left it untouched");
        }
        // --jet_idx highlights the leg the region is about: the HS leg in R1
        // (which one is the hard scatter?), the fake in R2 (can we kill it?).
        int    hi_idx   = isR1 ? c.idx_hs : c.idx_pu;
        const char* lbl = isR1 ? "HS" : "PU";
        std::printf("  cd python && python3 event_display.py --file_path \"%s\""
                    " --event_num %lld --extra_time %.2f --jet_idx %d --jet_label %s\n\n",
                    c.file_path.c_str(), c.entry, c.t_waves, hi_idx, lbl);
      }
    };

    printRegion("VBS REGION R1 - both candidate legs forward (HS vs PU)",
                "Ranked by |change in HS-minus-PU RpT margin| between ITk-only and WAVeS.",
                merged.cases_r1, true);
    printRegion("VBS REGION R2 FAILURES - fakes that SURVIVED the gate",
                "Ranked by the forward-PU RpT still standing after WAVeS. "
                "The success list below cannot show these: rW <= rZ always, so "
                "every entry there is 'suppressed' and failures rank last.",
                merged.cases_r2_fail, false);

    printRegion("VBS REGION R2 - forward PU leg + central HS leg",
                "Ranked by |change in forward-PU RpT| between ITk-only and WAVeS.",
                merged.cases_r2, false);
  }

  return 0;
}
