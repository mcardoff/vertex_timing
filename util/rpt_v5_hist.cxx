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
// Jet eta acceptance: 2.4 < |eta| < 3.8.
//
// Output: <OUTPUT_DIR>/hists/rpt_v5_hist.root

#include <TChain.h>
#include <TH1.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
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

// Per-track time-gate half-width in σ.  A 2σ cut over-trims genuine HS tracks
// when the vertex time is slightly mis-estimated, dragging the high-efficiency
// end of the ROC below ITk-only (worst in the >40 GeV slice).  Loosening it
// recovers that region at a small low-efficiency cost.
static constexpr double GATE_SIGMA = 2.5;

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

struct JetCompCase {
  std::string file_path;
  Long64_t    entry;
  int    jet_idx;
  double j_pt, j_eta;
  double rpt_mine, rpt_hgtd;
  double t_mine;    // WAVeS cluster time, used as the --extra_time annotation
  double delta;     // rpt_mine - rpt_hgtd  (positive -> WAVeS better)
  bool   isHS;
};

struct HurtJet {
  std::string file_path;
  Long64_t    entry;
  double      j_pt, j_eta, j_phi;
  double      rpt_z, rpt_mine;
  int         n_lost;
  double      t_mine;  // WAVeS cluster time, for the --extra_time annotation.
                        // Always meaningful whenever a HurtJet is recorded --
                        // see the collection site's comment for why.
};

// Keep the top max_n cases by |delta|.
static void insertCase(std::vector<JetCompCase>& v, JetCompCase c, int max_n = 5) {
  v.push_back(std::move(c));
  std::sort(v.begin(), v.end(),
            [](const JetCompCase& a, const JetCompCase& b) {
              return std::abs(a.delta) > std::abs(b.delta);
            });
  if ((int)v.size() > max_n) v.resize(max_n);
}

// Merge one thread's top-N candidates into the running merged top-N. Correct
// because the true global top-N is always a subset of the union of each
// thread's own (already-truncated) top-N lists.
static void mergeCases(std::vector<JetCompCase>& dst, std::vector<JetCompCase>& src, int max_n = 5) {
  for (auto& c : src) dst.push_back(std::move(c));
  std::sort(dst.begin(), dst.end(),
            [](const JetCompCase& a, const JetCompCase& b) {
              return std::abs(a.delta) > std::abs(b.delta);
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
  std::vector<Scenario> scen_lo = makeScenarios("_lo");  // 30–40 GeV
  std::vector<Scenario> scen_hi = makeScenarios("_hi");  // >40 GeV
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

  // Event-display diagnostic candidates (see JetCompCase/HurtJet doc comment above).
  std::vector<JetCompCase> cases_mine_lo, cases_hgtd_lo;   // 30-40 GeV
  std::vector<JetCompCase> cases_mine_hi, cases_hgtd_hi;   // >40 GeV
  std::vector<JetCompCase> cases_pu_mine_corrects_lo, cases_pu_mine_corrects_hi;
  std::vector<JetCompCase> cases_pu_mine_worse_lo, cases_pu_mine_worse_hi;
  std::vector<HurtJet>     hurt_events;
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
                                double sigma = 2.0) {
        std::vector<int> out;
        out.reserve(base.size());
        for (int idx : base) {
          bool apply = vtx_valid && branch.trackTimeValid[idx] == 1;
          if (!apply) { out.push_back(idx); continue; }
          double dt    = branch.trackTime[idx] - t_vtx;
          double var_t = branch.trackTimeRes[idx] * branch.trackTimeRes[idx];
          double pull  = std::abs(dt) / std::sqrt(2.25 * var_vtx + var_t);
          if (pull < sigma) out.push_back(idx);
        }
        return out;
      };

      std::vector<int> trk_hgtd   = applyTimeGate(trk_all, t_hgtd,   var_hgtd,   hgtd_vtx_valid, GATE_SIGMA);
      std::vector<int> trk_waves  = applyTimeGate(trk_all, t_waves,  var_waves,  waves_ok,       GATE_SIGMA);
      std::vector<int> trk_trkptz = applyTimeGate(trk_all, t_trkptz, var_trkptz, trkptz_ok,      GATE_SIGMA);

      // Build per-scenario sets once per event for O(1) ghost-index lookup.
      std::unordered_set<int> set_all   (trk_all.begin(),    trk_all.end());
      std::unordered_set<int> set_hgtd  (trk_hgtd.begin(),   trk_hgtd.end());
      std::unordered_set<int> set_trkptz(trk_trkptz.begin(), trk_trkptz.end());
      std::unordered_set<int> set_waves (trk_waves.begin(),  trk_waves.end());

      // ── Fill jets into pT slices. ─────────────────────────────────────────────
      auto fillJets = [&](std::vector<Scenario>& sv, double pt_lo, double pt_hi) {
        for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
          if (branch.isJetRemoved(j)) continue;  // lepton-overlap removed (Z+jets)
          double j_pt  = branch.topoJetPt[j];
          double j_eta = branch.topoJetEta[j];
          double j_phi = branch.topoJetPhi[j];
          if (j_pt <= pt_lo || j_pt >= pt_hi) continue;
          if (std::abs(j_eta) < JET_ETA_MIN || std::abs(j_eta) > JET_ETA_MAX) continue;
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
          fill(sv[0], set_all);                       // ITk-only
          fill(sv[1], set_hgtd);                      // HGTD t0 (Athena)
          fill(sv[2], set_trkptz);                    // TRKPTZ t0
          fill(sv[3], set_waves);                     // WAVeS t0

          // Untimed floor accounting, per slice.
          for (int idx : ghost) {
            if (!set_all.count(idx)) continue;
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

      fillJets(state.scen_lo, 30.0, 40.0);
      fillJets(state.scen_hi, 40.0, 1e9);

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
      {
        std::vector<int> passPtIdx;
        int nPt = 0, nPtEta = 0;
        branch.collectPtPassingJets(passPtIdx, nPt, nPtEta);
        auto pair = branch.calcBestVbsPair(passPtIdx);

        const bool passVbsTopology = pair.valid()
                                  && pair.mjj  >= VBS_JET_MJJ
                                  && pair.dEta >= VBS_JET_D_ETA;

        if (passVbsTopology) {
          auto isFwd = [&](int j) {
            double e = std::abs((double)branch.topoJetEta[j]);
            return e > JET_ETA_MIN && e < JET_ETA_MAX;
          };
          auto isCentral = [&](int j) {
            return std::abs((double)branch.topoJetEta[j]) < CENTRAL_ETA_MAX;
          };
          auto label = [&](int j, bool& hs, bool& pu) {
            double e = branch.topoJetEta[j], p = branch.topoJetPhi[j];
            hs = paperIsHS(e, p);
            pu = paperIsPU(e, p);
          };

          const int a = pair.idxI, b = pair.idxJ;
          bool aHS, aPU, bHS, bPU;
          label(a, aHS, aPU);
          label(b, bHS, bPU);

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
            put(sv[0], set_all);
            put(sv[1], set_hgtd);
            put(sv[2], set_trkptz);
            put(sv[3], set_waves);
          };

          // R1: both forward, exactly one HS and one PU leg.
          if (isFwd(a) && isFwd(b)) {
            if (aHS && bPU) { fillRegion(state.scen_r1, a, true);
                              fillRegion(state.scen_r1, b, false); }
            else if (bHS && aPU) { fillRegion(state.scen_r1, b, true);
                                   fillRegion(state.scen_r1, a, false); }
          }

          // R2: forward PU leg + central HS leg. Forward PU leg only.
          if (isFwd(a) && aPU && isCentral(b) && bHS) fillRegion(state.scen_r2, a, false);
          if (isFwd(b) && bPU && isCentral(a) && aHS) fillRegion(state.scen_r2, b, false);
        }
      }

      // ── Event-display diagnostics (see JetCompCase/HurtJet doc comment
      //    near the top of this file). Full ntuple file path + local
      //    (per-file) entry number -- same TTreeProcessorMT-safe pattern as
      //    src/clustering_hist.cxx (no outer TChain available inside this
      //    lambda; reader.GetTree() gives the currently-loaded per-file
      //    constituent tree directly). Collection always runs (cheap --
      //    a handful of comparisons per event); only the final print is
      //    gated behind PRINT_EVENT_DISPLAYS.
      std::string filePath   = reader.GetTree()->GetCurrentFile()->GetName();
      Long64_t    localEntry = reader.GetTree()->GetReadEntry();

      // WAVeS ("mine") vs HGTD jet-level comparison.
      if (hgtd_vtx_valid && waves_ok) {
        for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
          if (branch.isJetRemoved(j)) continue;  // lepton-overlap removed (Z+jets)
          double j_pt  = branch.topoJetPt[j];
          double j_eta = branch.topoJetEta[j];
          double j_phi = branch.topoJetPhi[j];
          bool in_lo = (j_pt > 30.0 && j_pt < 40.0);
          bool in_hi = (j_pt > 40.0);
          if (!in_lo && !in_hi) continue;
          if (std::abs(j_eta) < JET_ETA_MIN || std::abs(j_eta) > JET_ETA_MAX) continue;
          bool isHS = paperIsHS(j_eta, j_phi);
          bool isPU = paperIsPU(j_eta, j_phi);
          if (!isHS && !isPU) continue;

          const auto& ghost = branch.topoJetGhostTrackIdx[j];
          double rpt_hgtd_j = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, set_hgtd);
          double rpt_mine_j = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, set_waves);
          double delta = rpt_mine_j - rpt_hgtd_j;
          if (std::abs(delta) < 0.05) continue;  // skip trivial differences

          JetCompCase c{filePath, localEntry, j, j_pt, j_eta,
                        rpt_mine_j, rpt_hgtd_j, t_waves, delta, isHS};
          if (in_lo) {
            if (delta > 0) insertCase(state.cases_mine_lo, c);
            else           insertCase(state.cases_hgtd_lo, c);
          } else {
            if (delta > 0) insertCase(state.cases_mine_hi, c);
            else           insertCase(state.cases_hgtd_hi, c);
          }
          // PU jets where HGTD assigns significant RpT (mistag) but WAVeS suppresses.
          if (isPU && delta < 0 && rpt_hgtd_j > 0.1) {
            if (in_lo) insertCase(state.cases_pu_mine_corrects_lo, c);
            else       insertCase(state.cases_pu_mine_corrects_hi, c);
          }
          // PU jets where WAVeS gives higher RpT than HGTD (makes PU look more HS-like).
          if (isPU && delta > 0 && rpt_mine_j > 0.1) {
            if (in_lo) insertCase(state.cases_pu_mine_worse_lo, c);
            else       insertCase(state.cases_pu_mine_worse_hi, c);
          }
        }
      }

      // Hurt-HS diagnostic: HS jets (30-40 GeV) where the WAVeS time gate
      // removed ghost-associated tracks relative to the z-only baseline,
      // lowering RpT. Not gated on waves_ok explicitly -- when it's false,
      // set_waves falls back to set_all (applyTimeGate's no-op path), so
      // rpt_mine>=rpt_z always holds and the n_lost/rpt_mine<rpt_z checks
      // below naturally skip these events, matching v4's original structure.
      // (This also means t_waves below is always meaningful whenever a
      // HurtJet is actually recorded -- waves_ok must have been true, or the
      // event would already have been skipped.)
      // Capped per-thread (not globally, since each worker races
      // independently); merged lists are concatenated, not re-sorted,
      // matching v4's original first-come-first-served collection order.
      if (state.hurt_events.size() < 25) {
        for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
          if (branch.isJetRemoved(j)) continue;
          double j_pt  = branch.topoJetPt[j];
          double j_eta = branch.topoJetEta[j];
          double j_phi = branch.topoJetPhi[j];
          if (j_pt <= 30.0 || j_pt >= 40.0) continue;
          if (std::abs(j_eta) < JET_ETA_MIN || std::abs(j_eta) > JET_ETA_MAX) continue;
          if (!paperIsHS(j_eta, j_phi)) continue;

          const auto& ghost = branch.topoJetGhostTrackIdx[j];
          int n_lost = 0;
          for (int idx : ghost)
            if (set_all.count(idx) && !set_waves.count(idx)) ++n_lost;
          if (n_lost == 0) continue;

          double rpt_z    = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, set_all);
          double rpt_mine = computeRpT(&branch, ghost, j_pt, j_eta, j_phi, set_waves);
          if (rpt_mine >= rpt_z) continue;  // shouldn't happen, but skip no-op cases

          state.hurt_events.push_back({filePath, localEntry, j_pt, j_eta, j_phi,
                                        rpt_z, rpt_mine, n_lost, t_waves});
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
    mergeCases(merged.cases_mine_lo, other.cases_mine_lo);
    mergeCases(merged.cases_hgtd_lo, other.cases_hgtd_lo);
    mergeCases(merged.cases_mine_hi, other.cases_mine_hi);
    mergeCases(merged.cases_hgtd_hi, other.cases_hgtd_hi);
    mergeCases(merged.cases_pu_mine_corrects_lo, other.cases_pu_mine_corrects_lo);
    mergeCases(merged.cases_pu_mine_corrects_hi, other.cases_pu_mine_corrects_hi);
    mergeCases(merged.cases_pu_mine_worse_lo, other.cases_pu_mine_worse_lo);
    mergeCases(merged.cases_pu_mine_worse_hi, other.cases_pu_mine_worse_hi);
    for (auto& h : other.hurt_events)
      if (merged.hurt_events.size() < 25) merged.hurt_events.push_back(h);
  }

  std::cout << "\nFINISHED PROCESSING\n";

  // --- Save every histogram + scalar accumulator to a ROOT file ---
  const std::string histPath = MyUtl::histFilePath("rpt_v5_hist.root");
  MyUtl::HistWriter writer(histPath);
  saveScenarios(writer, merged.scen_lo);
  saveScenarios(writer, merged.scen_hi);
  saveScenarios(writer, merged.scen_r1);
  saveScenarios(writer, merged.scen_r2);
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

  // --- Event-display diagnostics: WAVeS ("mine") vs HGTD jet-level
  //     comparison, restored from rpt_v4.cxx (see doc comment near the top
  //     of this file). Gated behind PRINT_EVENT_DISPLAYS -- flip that flag
  //     to true and rerun to print ready-to-run event_display.py commands. ---
  if (PRINT_EVENT_DISPLAYS) {
    std::cout << "\n=== TIMING-HURT HS JETS (30-40 GeV, WAVeS scenario) ===\n";
    std::cout << "  HS jets where the WAVeS time gate removed >=1 ghost-associated track, lowering RpT.\n\n";
    for (auto& h : merged.hurt_events) {
      std::printf("  jet pT=%.1f  eta=%.2f  phi=%.2f  RpT: %.3f->%.3f  tracks_lost=%d\n",
                  h.j_pt, h.j_eta, h.j_phi, h.rpt_z, h.rpt_mine, h.n_lost);
      std::printf("  cd python && python3 event_display.py --file_path \"%s\" --event_num %lld --extra_time %.2f\n\n",
                  h.file_path.c_str(), h.entry, h.t_mine);
    }
    if (merged.hurt_events.empty())
      std::cout << "  (none found)\n";

    auto printCases = [](const char* title, const std::vector<JetCompCase>& cases) {
      std::cout << "\n" << title << ":\n";
      if (cases.empty()) { std::cout << "  (none found)\n"; return; }
      for (const auto& c : cases) {
        std::printf("  jet pT=%.1f GeV  eta=%.2f  %s  RpT: mine=%.3f  hgtd=%.3f  delta=%.3f\n",
                    c.j_pt, c.j_eta, c.isHS ? "HS" : "PU",
                    c.rpt_mine, c.rpt_hgtd, c.delta);
        std::printf("  cd python && python3 event_display.py --file_path \"%s\" --event_num %lld"
                    " --extra_time %.2f --jet_idx %d --jet_label %s"
                    " --rpt_hgtd %.3f --rpt_mine %.3f\n\n",
                    c.file_path.c_str(), c.entry, c.t_mine, c.jet_idx,
                    c.isHS ? "HS" : "PU", c.rpt_hgtd, c.rpt_mine);
      }
    };

    std::cout << "\n=== WAVeS vs HGTD: JET-LEVEL COMPARISON ===\n";
    printCases("CASE 1 - WAVeS improves (30-40 GeV, WAVeS RpT > HGTD RpT)", merged.cases_mine_lo);
    printCases("CASE 2 - WAVeS improves (>40 GeV, WAVeS RpT > HGTD RpT)",   merged.cases_mine_hi);
    printCases("CASE 3 - HGTD better  (30-40 GeV, HGTD RpT > WAVeS RpT)",   merged.cases_hgtd_lo);
    printCases("CASE 4 - HGTD better  (>40 GeV, HGTD RpT > WAVeS RpT)",     merged.cases_hgtd_hi);

    std::cout << "\n=== PU MISTAG CORRECTION: HGTD mistags PU as HS, WAVeS corrects ===\n";
    std::cout << "  (PU jets with rpt_hgtd > 0.1 where WAVeS gives lower RpT)\n";
    printCases("CASE 5 - WAVeS corrects PU mistag (30-40 GeV)", merged.cases_pu_mine_corrects_lo);
    printCases("CASE 6 - WAVeS corrects PU mistag (>40 GeV)",   merged.cases_pu_mine_corrects_hi);

    std::cout << "\n=== PU MISTAG WORSENING: WAVeS pushes PU jet closer to HS ===\n";
    std::cout << "  (PU jets with rpt_mine > 0.1 where WAVeS gives higher RpT than HGTD)\n";
    printCases("CASE 7 - WAVeS worsens PU mistag (30-40 GeV)", merged.cases_pu_mine_worse_lo);
    printCases("CASE 8 - WAVeS worsens PU mistag (>40 GeV)",   merged.cases_pu_mine_worse_hi);
  }

  return 0;
}
