// rpt_v5_hist.cxx — Histogramming stage of the rpt_v5 split.
//
// Runs the same TTreeProcessorMT event loop and per-thread ThreadState merge
// that the former monolithic rpt_v5.cxx did, then writes the merged
// Hard-Scatter/Pile-Up RpT histograms (5 scenarios x 2 pT slices) plus the
// scalar event-count/floor-counter accumulators to
// <OUTPUT_DIR>/hists/rpt_v5_hist.root via histogram_io.h, instead of doing
// any ROC/console-summary/PDF work directly. rpt_v5_plot.cxx reads that file
// back and does all of that -- see CLAUDE.md's "Main Executables" section.
//
// Scenarios:
//   1. zonly        — z-significance tracks only, no time gate   ("ITk-only")
//   2. hgtd         — ntuple RecoVtx_time / RecoVtx_timeRes gate  ("HGTD t_{0}")
//   3. waves        — WAVeS-selected cluster time gate            ("WAVeS t_{0}")
//   4. waves_misas  — WAVeS time gate, events gated on HS timing  ("WAVeS t_{0} + clean timing")
//                     purity ≥ MISAS_PURITY_CUT (event filter)
//   5. truth        — reco track times gated vs truth vertex t_{0} ("Truth t_{0}")
//
// No event-level selection besides the MISAS filter on scenario 4 — every forward
// jet in the acceptance contributes an independent RpT measurement.
//
// Two jet pT windows:
//   Slice A: 30 < pT < 40 GeV
//   Slice B: pT > 40 GeV
// Jet eta acceptance: 2.4 < |eta| < 3.8 (forward, HGTD-covered), plus a
// |eta| < 2.4 central baseline filled in the same run and pT slices.
//
// Output: <OUTPUT_DIR>/hists/rpt_v5_hist.root

#include <TChain.h>
#include <TH1.h>
#include <TStyle.h>
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

// Reference-study truth-t0 smearing (util/myJet_ana_fr.C).
static constexpr double TRUTH_VTX_SMEAR = 10.0;  // ps, on the HS vertex time
static constexpr double TRUTH_TRK_SMEAR = 30.0;  // ps, on each track's own vertex time

// Per-track time-gate half-width in σ.  A 2σ cut over-trims genuine HS tracks
// when the vertex time is slightly mis-estimated, dragging the high-efficiency
// end of the ROC below ITk-only (worst in the >40 GeV slice).  Loosening it
// recovers that region at a small low-efficiency cost.
static constexpr double GATE_SIGMA = 2.5;

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
struct Inflation { double hgtd, trkptz, waves; };

// Keyed on MyUtl::SAMPLE_NAME. The local default run (no --sample, so an empty
// SAMPLE_NAME) reads the VBF ntuples, so it correctly falls through to the VBF
// row rather than needing an entry of its own.
static Inflation inflationFor(const std::string& sample) {
  if (sample == "zjets") return {1.61, 1.63, 1.65};
  if (sample == "dijet") return {1.53, 1.45, 1.46};
  return {1.48, 1.39, 1.38};  // vbf, and the local default run
}

// Resolved once in main() before the event loop starts, then only read by the
// worker threads -- write-once-before-fork, so no synchronisation is needed.
static Inflation INFL = {1.48, 1.39, 1.38};

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
// Ported verbatim (coefficients and all) from a colleague's RpT study,
// util/myJet_ana_fr.C::getNewDzpara. Returns the expected |z0 - z_vtx| scale in
// mm for a track of given |eta| and pT: a 6th-order polynomial in |eta| whose
// coefficients are selected from five pT bins.
//
// This REPLACES the z0-significance cut (|z0-z_vtx|/sqrt(var_z0) < N) that
// rpt_v5 used previously. The two are not equivalent: the significance form
// leans on the per-track covariance, which the z0_pull_diag study found
// underestimates the true spread by ~15% (hence Z0_VAR_INFLATION elsewhere in
// this codebase), while this is an empirical parameterization of the observed
// resolution. Using it makes our RpT numerator directly comparable to the
// reference study's.
//
// Note the pT bins are half-open on the low side and the <=1.5 bin catches
// everything below, matching the original's if-chain exactly (no else-if, so a
// value landing on a boundary takes the LAST matching branch -- preserved here
// deliberately rather than "cleaned up", since changing it would silently shift
// which coefficients boundary-pT tracks get).
// -----------------------------------------------------------------------------
static double getNewDzpara(double eta, double pt) {
  eta = std::abs(eta);
  double p[7] = {0, 0, 0, 0, 0, 0, 0};
  auto set = [&p](const double (&src)[7]) { std::copy(src, src + 7, p); };

  if (pt <= 1.5) {
    static const double c[7] = { 0.0314036, 0.790955, -2.65987, 3.62073, -2.18228, 0.614866, -0.0634521 };
    set(c);
  }
  if (pt > 1.5 && pt <= 2.5) {
    static const double c[7] = { 0.0229273, 0.540101, -1.80727, 2.45187, -1.47382, 0.414345, -0.0426769 };
    set(c);
  }
  if (pt > 2.5 && pt <= 5.0) {
    static const double c[7] = { 0.0163773, 0.345112, -1.14474, 1.54382, -0.923523, 0.258617, -0.0265446 };
    set(c);
  }
  if (pt > 5.0 && pt <= 10.0) {
    static const double c[7] = { 0.010919, 0.179329, -0.581971, 0.773186, -0.45679, 0.126608, -0.012875 };
    set(c);
  }
  if (pt > 10.0) {
    static const double c[7] = { 0.00835945, 0.0957783, -0.299255, 0.38722, -0.22351, 0.0607521, -0.00606524 };
    set(c);
  }

  double d = p[0];
  double e = 1.0;
  for (int k = 1; k < 7; ++k) { e *= eta; d += p[k] * e; }
  return std::abs(d);
}

// Multiple of getNewDzpara a track's |z0 - z_vtx| must stay within. 1.4 is the
// reference study's working point (its comment notes 0.8 discriminates better;
// left at 1.4 to match what the study actually ran).
static constexpr double DZ0_PARA_SCALE = 1.4;

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
  std::vector<RegionCase> cases_r1, cases_r2;
  // Per-scenario pull-width accumulators (see PRINT_PULL_DIAG).
  double pull_dt2_hgtd = 0, pull_var_hgtd = 0;
  double pull_dt2_trkptz = 0, pull_var_trkptz = 0;
  double pull_dt2_waves = 0, pull_var_waves = 0;
  long   pull_n_hgtd = 0, pull_n_trkptz = 0, pull_n_waves = 0;
  // First moment and out-of-core counts, so the diagnostic can separate a
  // systematic OFFSET (which an inflation cannot fix) from genuine spread, and
  // report how much of each scenario's distribution the core window discards.
  double pull_dt_hgtd = 0, pull_dt_trkptz = 0, pull_dt_waves = 0;
  long   pull_ntail_hgtd = 0, pull_ntail_trkptz = 0, pull_ntail_waves = 0;
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

  // --- Sample selection (--sample=vbf|zjets|dijet; default: local VBF ntuple) ---
  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);  // --vbs-deta=<x>; sets SELECTION_TAG
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  MyUtl::OVERLAP_REMOVAL = sample.overlapRemoval;  // Z+jets lepton–jet overlap removal
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  if (MyUtl::SAMPLE_NAME.empty())
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/hists");
  unsigned nThreads = MyUtl::resolveThreads(argc, argv);

  // Per-sample vertex-time calibration. Resolved here, before any worker thread
  // exists, so the event loop only ever reads it. Echoed so a run's stdout
  // records which calibration produced its histograms -- the PRINT_PULL_DIAG
  // table at the end then prints the freshly measured ratio beside these in its
  // "in use" column, making a stale entry visible in the same log.
  INFL = inflationFor(MyUtl::SAMPLE_NAME);
  std::printf("[calibration] vertex-time inflation for '%s': "
              "hgtd %.2f  trkptz %.2f  waves %.2f\n",
              MyUtl::SAMPLE_NAME.empty() ? "local (vbf)" : MyUtl::SAMPLE_NAME.c_str(),
              INFL.hgtd, INFL.trkptz, INFL.waves);

  TChain chain("ntuple");
  setupChain(chain, sample.ntupleDir.c_str());
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found.  Aborting.\n";
    return 1;
  }
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

  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  // Optional --max-events=<N> cap for quick local checks (see resolveMaxEvents
  // in sample_config.h). -1 means unlimited -- every existing invocation
  // without the flag is unaffected.
  const Long64_t maxEvents     = MyUtl::resolveMaxEvents(argc, argv);
  const Long64_t progressDenom = (maxEvents > 0) ? std::min(maxEvents, N_EVENT) : N_EVENT;
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
        std::cout << "Progress: " << n << "/" << progressDenom << "\r" << std::flush;

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
      //    vertex by the reference study's parameterized |Δz₀| cut, for the
      //    z-only baseline / ITk-only scenario. See getNewDzpara above for why
      //    this replaced the previous z₀-significance cut. ──────────────────────
      std::vector<int> trk_all;
      const double vtxZ = branch.recoVtxZ[0];
      for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
        double trkPt = branch.trackPt[trk];
        if (trkPt < MIN_TRACK_PT || trkPt > MAX_TRACK_PT) continue;
        if (!branch.trackQuality[trk]) continue;
        double dz = std::abs(branch.trackZ0[trk] - vtxZ);
        if (dz / getNewDzpara(branch.trackEta[trk], trkPt) > DZ0_PARA_SCALE) continue;
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
      if (!clusters.empty()) {
        auto best   = chooseCluster(clusters, Score::WAVES);
        t_waves     = best.calculateTime(Score::WAVES, &branch);  // in-jet refined
        var_waves   = best.sigmas[0] * best.sigmas[0];
        waves_ok    = true;

        auto bestT  = chooseCluster(clusters, Score::TRKPTZ);
        t_trkptz    = bestT.calculateTime(Score::TRKPTZ, &branch);
        var_trkptz  = bestT.sigmas[0] * bestT.sigmas[0];
        trkptz_ok   = true;
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
      // Central list gets the identical gate. It is a no-op there in practice
      // (no HGTD coverage below |eta| 2.4, so no track carries a valid time),
      // but applying it keeps all four central scenarios defined exactly as
      // their forward counterparts -- which is what makes their agreement a
      // meaningful check rather than a tautology.
      std::vector<int> trk_hgtd_cen   = applyTimeGate(trk_all_cen, t_hgtd,   var_hgtd,   hgtd_vtx_valid, GATE_SIGMA, INFL.hgtd);
      std::vector<int> trk_waves_cen  = applyTimeGate(trk_all_cen, t_waves,  var_waves,  waves_ok,       GATE_SIGMA, INFL.waves);
      std::vector<int> trk_trkptz_cen = applyTimeGate(trk_all_cen, t_trkptz, var_trkptz, trkptz_ok,      GATE_SIGMA, INFL.trkptz);

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
      }

      // ── Truth-t0 reference scenario ────────────────────────────────────────
      // Smear once per TRACK (not per list), so a track appearing in both the
      // forward and central association lists carries the same time in each.
      {
        UInt_t sd = (UInt_t)std::llround(std::abs(branch.truthVtxZ[0])    * 1e4)
                  ^ ((UInt_t)std::llround(std::abs(branch.truthVtxTime[0]) * 1e3) << 11)
                  ^ ((UInt_t)branch.trackZ0.GetSize() << 23);
        state.rng.SetSeed(sd ? sd : 1u);
      }
      const double t_truth_vtx = state.rng.Gaus(branch.truthVtxTime[0], TRUTH_VTX_SMEAR);
      const size_t nTrkAll = branch.trackZ0.GetSize();
      std::vector<float> t_truth_trk(nTrkAll, 0.f);
      std::vector<char>  t_truth_ok (nTrkAll, 0);
      for (size_t i = 0; i < nTrkAll; ++i) {
        int tv = branch.trackToTruthvtx[i];
        if (tv < 0 || tv >= (int)branch.truthVtxTime.GetSize()) continue;
        t_truth_trk[i] = (float)state.rng.Gaus(branch.truthVtxTime[tv], TRUTH_TRK_SMEAR);
        t_truth_ok[i]  = 1;
      }
      // Same pull form as applyTimeGate, with infl = 1: these sigmas are exact
      // by construction, so there is nothing to inflate. A track with no truth
      // vertex link is kept unconditionally, as an untimed track would be.
      auto truthGate = [&](const std::vector<int>& base) {
        std::vector<int> out; out.reserve(base.size());
        const double den = std::sqrt(TRUTH_VTX_SMEAR * TRUTH_VTX_SMEAR
                                   + TRUTH_TRK_SMEAR * TRUTH_TRK_SMEAR);
        for (int idx : base) {
          // Only where HGTD actually measured. Idealising the TIME of a track
          // the detector never timed would invent coverage it does not have --
          // which shows up immediately in the central baseline, where |eta| <
          // 2.4 is outside acceptance and every timing row must therefore lie
          // exactly on ITk-only. The reference macro omits this check because
          // it only ever looks at forward jets, where the distinction is moot.
          if (branch.trackTimeValid[idx] != 1) { out.push_back(idx); continue; }
          if (!t_truth_ok[idx]) { out.push_back(idx); continue; }
          if (std::abs(t_truth_trk[idx] - t_truth_vtx) / den < GATE_SIGMA) out.push_back(idx);
        }
        return out;
      };
      std::vector<int> trk_truth     = truthGate(trk_all);
      std::vector<int> trk_truth_cen = truthGate(trk_all_cen);

      // Same idealised track times, gated against the REAL WAVeS vertex time
      // instead of the smeared truth one, so this row and the truth row differ
      // in exactly one variable: where the vertex time comes from. Both respect
      // HGTD coverage (see truthGate), so an untimed track is untouched here
      // just as it is in every real scenario. The vertex sigma is still the
      // real cluster's and carries its inflation; the track sigma is exactly
      // 30 ps by construction.
      auto smearGate = [&](const std::vector<int>& base, double t_vtx,
                           double var_vtx, bool vtx_ok, double infl) {
        std::vector<int> out; out.reserve(base.size());
        const double den = std::sqrt(infl * infl * var_vtx
                                   + TRUTH_TRK_SMEAR * TRUTH_TRK_SMEAR);
        for (int idx : base) {
          if (branch.trackTimeValid[idx] != 1) { out.push_back(idx); continue; }
          if (!vtx_ok || !t_truth_ok[idx] || den <= 0) { out.push_back(idx); continue; }
          if (std::abs(t_truth_trk[idx] - t_vtx) / den < GATE_SIGMA) out.push_back(idx);
        }
        return out;
      };
      std::vector<int> trk_wsm     = smearGate(trk_all,     t_waves, var_waves, waves_ok, INFL.waves);
      std::vector<int> trk_wsm_cen = smearGate(trk_all_cen, t_waves, var_waves, waves_ok, INFL.waves);

      // Build per-scenario sets once per event for O(1) ghost-index lookup.
      struct TrackSets { std::unordered_set<int> all, hgtd, trkptz, waves, waves_smear, truth; };
      TrackSets fwd{ {trk_all.begin(),    trk_all.end()},
                     {trk_hgtd.begin(),   trk_hgtd.end()},
                     {trk_trkptz.begin(), trk_trkptz.end()},
                     {trk_waves.begin(),  trk_waves.end()},
                     {trk_wsm.begin(),    trk_wsm.end()},
                     {trk_truth.begin(),  trk_truth.end()} };
      TrackSets cen{ {trk_all_cen.begin(),    trk_all_cen.end()},
                     {trk_hgtd_cen.begin(),   trk_hgtd_cen.end()},
                     {trk_trkptz_cen.begin(), trk_trkptz_cen.end()},
                     {trk_waves_cen.begin(),  trk_waves_cen.end()},
                     {trk_wsm_cen.begin(),    trk_wsm_cen.end()},
                     {trk_truth_cen.begin(),  trk_truth_cen.end()} };

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
          fill(sv[4], S.waves_smear);                 // WAVeS t0, idealised 30 ps tracks
          fill(sv[5], S.truth);                       // truth t0, 10 (+) 30 ps

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
      // classifyVbsRegion -- shared with the clustering side's VBF_R1/VBF_R2
      // scores specifically so the two analyses cannot drift into meaning
      // different event sets by the same name. It also owns the paper-HS/PU
      // labelling that used to live in this file's local lambdas.
      {
        int fwdHS = -1, fwdPU = -1;
        auto region = branch.classifyVbsRegion(JET_ETA_MIN, JET_ETA_MAX,
                                               CENTRAL_ETA_MAX, &fwdHS, &fwdPU);

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
            put(sv[4], fwd.waves_smear);
            put(sv[5], fwd.truth);
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
          }
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
    mergeRegionCases(merged.cases_r1, other.cases_r1);
    mergeRegionCases(merged.cases_r2, other.cases_r2);
  }

  std::cout << "\nFINISHED PROCESSING\n";

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
  writer.WriteRunMeta(MyUtl::ENERGY_LABEL, merged.n_total);
  writer.Close();
  std::cout << "Wrote histograms to " << histPath << "\n";

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
                      c.metric < 0 ? "SUPPRESSED the fake" : "raised it");
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
    printRegion("VBS REGION R2 - forward PU leg + central HS leg",
                "Ranked by |change in forward-PU RpT| between ITk-only and WAVeS.",
                merged.cases_r2, false);
  }

  return 0;
}
