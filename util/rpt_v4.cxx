// rpt_v4.cxx — Jet-level RpT analysis, advisor-requested pT split: 30–40 GeV and >40 GeV.
//
// Three scenarios (no oracle gates):
//   1. zonly  — z-significance tracks only, no time gate
//   2. hgtd   — ntuple RecoVtx_time / RecoVtx_timeRes (per-track gate)
//   3. mine   — TRKPTZ-selected cluster time (per-track gate)
//
// No event-level selection — every event is processed, every forward jet
// in the acceptance contributes an independent RpT measurement.
//
// Two jet pT windows (filled from the same per-event track lists):
//   Slice A: 30 < pT < 40 GeV
//   Slice B: pT > 40 GeV
// Jet eta acceptance: 2.4 < |eta| < 3.8.
//
// Output: figs/rpt_plots/rpt_v4.pdf

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <boost/filesystem.hpp>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "clustering_constants.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

#define debug false

using namespace MyUtl;

// Jet acceptance for this script (paper values, distinct from rpt_v2).
static constexpr double JET_ETA_MIN = 2.4;
static constexpr double JET_ETA_MAX = 3.8;

// -----------------------------------------------------------------------------
// ROC: HS efficiency vs PU rejection (1 / mistag).
// -----------------------------------------------------------------------------
TGraph* generate_roc(TH1D* PU_hist, TH1D* HS_hist) {
  int bin = PU_hist->GetNbinsX();
  std::vector<float> vx, vy;
  for (int i = 1; i <= bin; ++i) {
    double HS_eff    = HS_hist->Integral(i, bin + 1) / HS_hist->Integral();
    double PU_mistag = PU_hist->Integral(i, bin + 1) / PU_hist->Integral();
    if (std::abs(PU_mistag) > 1e-6 && HS_eff < 0.99) {
      vx.push_back(HS_eff);
      vy.push_back(1.0 / PU_mistag);
    }
  }
  return new TGraph((int)vx.size(), vx.data(), vy.data());
}

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

static inline double dR(double j_eta, double j_phi, double t_eta, double t_phi) {
  double deta = j_eta - t_eta;
  double dphi = TVector2::Phi_mpi_pi(j_phi - t_phi);
  return std::sqrt(deta * deta + dphi * dphi);
}

// RpT using ghost association (paper definition): sum pT of tracks that are
// both ghost-associated to the jet AND in the caller's z₀-selected set.
static double computeRpT(BranchPointerWrapper* b,
                         const std::vector<int>& ghost_indices,
                         double j_pt,
                         const std::unordered_set<int>& assoc_set) {
  double sumpt = 0.0;
  for (int idx : ghost_indices)
    if (assoc_set.count(idx)) sumpt += b->trackPt[idx];
  return sumpt / j_pt;
}

struct HurtJet {
  std::string file_num;
  long        entry;
  double      j_pt, j_eta, j_phi;
  double      rpt_z, rpt_mine;
  int         n_lost;
};

struct JetCompCase {
  std::string file_num;
  long   entry;
  int    jet_idx;
  double j_pt, j_eta;
  double rpt_mine, rpt_hgtd;
  double t_trkptz;
  double delta;    // rpt_mine − rpt_hgtd  (positive → mine better)
  bool   isHS;
};

// Keep the top N cases by |delta|.
static void insertCase(std::vector<JetCompCase>& v, JetCompCase c, int max_n = 5) {
  v.push_back(std::move(c));
  std::sort(v.begin(), v.end(),
            [](const JetCompCase& a, const JetCompCase& b){
              return std::abs(a.delta) > std::abs(b.delta);
            });
  if ((int)v.size() > max_n) v.resize(max_n);
}

int main() {
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found.  Aborting.\n";
    return 1;
  }

  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  // Paper jet-label helper (ATL-HGTD-PUB-2022-001 Sec. 3).
  // HS  : dR(reco, truthHS) < 0.3  AND  truthHS_pT > 10 GeV
  // PU  : dR > 0.6 from ALL truth jets (HS + ITPU + OOTPU) with pT > 4 GeV
  auto paperIsHS = [&](double j_eta, double j_phi) {
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      if (branch.truthHSJetPt[t] < 10.0) continue;
      if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.3)
        return true;
    }
    return false;
  };
  // PU: dR > 0.6 from any truth HS jet with pT > 4 GeV (Sec. 3).
  // Only HS truth jets enter — ITPU/OOTPU are irrelevant for this cut.
  auto paperIsPU = [&](double j_eta, double j_phi) {
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      if (branch.truthHSJetPt[t] < 4.0) continue;
      if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.6)
        return false;
    }
    return true;
  };

  // Non-uniform binning: 0.01-wide bins in [0, 2.5] for ROC granularity,
  // then one wide bin [2.5, 375] to capture the tail without bloating memory.
  static std::vector<double> rpt_bins = []() {
    std::vector<double> b;
    b.reserve(202);
    for (int i = 0; i <= 200; ++i) b.push_back(0.0125 * i);
    b.push_back(375.0);
    return b;
  }();
  const int nbin = (int)rpt_bins.size() - 1;  // 201

  auto makeHist = [&](const char* name, const char* title) {
    return new TH1D(name, title, nbin, rpt_bins.data());
  };

  auto makeScenarios = [&](const std::string& suffix) {
    std::vector<Scenario> s = {
      {"zonly",      "z-only (no timing)",            C05, nullptr, nullptr},
      {"hgtd",       "HGTD ntuple t_{0}",             C01, nullptr, nullptr},
      {"mine",       "My algo (TRKPTZ)",                    C02, nullptr, nullptr},
      {"mine_pure",  "My algo (TRKPTZ) + pure cluster",  C03, nullptr, nullptr},
      {"mine_misas",  "My algo (TRKPTZ) + clean timing",  C04, nullptr, nullptr},
      {"mine_truth",  "Truth vertex t_{0}",               C06, nullptr, nullptr},
      {"mine_hsclust", "HS + in-time + top-5 PU clustering", C07, nullptr, nullptr},
    };
    for (auto& sc : s) {
      sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                         ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
      sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                         ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
    }
    return s;
  };

  std::vector<Scenario> scen_lo = makeScenarios("_lo");  // 30–40 GeV
  std::vector<Scenario> scen_hi = makeScenarios("_hi");  // >40 GeV

  // "Flagged vs Correcting" scenarios (advisor request):
  //   Correcting (mine_misas): clean-timing gate — apply timing when right, SKIP jet when wrong
  //   Flagging   (flagged):    clean-timing gate — apply timing when right, z-only when wrong
  //   The comparison reveals whether including wrong-timing jets at z-only performance
  //   helps or hurts relative to excluding them entirely.
  auto makeFlaggedScenarios = [&](const std::string& suffix) {
    std::vector<Scenario> s = {
      {"zonly_f",    "z-only (no timing)",                              C05, nullptr, nullptr},
      {"correcting", "My algo + correcting showers",                  C04, nullptr, nullptr},
      {"flagging",   "My algo + tagging showers",                      C02, nullptr, nullptr},
    };
    for (auto& sc : s) {
      sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                         ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
      sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                         ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
    }
    return s;
  };
  std::vector<Scenario> flag_lo = makeFlaggedScenarios("_flo");  // 30–40 GeV
  std::vector<Scenario> flag_hi = makeFlaggedScenarios("_fhi");  // >40 GeV
  // Append baseline TRKPTZ as a reference (reuses scen_lo/hi histograms; no extra fill needed).
  flag_lo.push_back({"mine_ref", "My algo (TRKPTZ)",  C02, nullptr, nullptr});
  flag_hi.push_back({"mine_ref", "My algo (TRKPTZ)",  C02, nullptr, nullptr});

  // Jet-level mine-vs-HGTD comparison cases (top 5 by |Δ RpT| each).
  std::vector<JetCompCase> cases_mine_lo, cases_hgtd_lo;  // 30–50 GeV
  std::vector<JetCompCase> cases_mine_hi, cases_hgtd_hi;  // 50–70 GeV
  // PU jets where HGTD gives high RpT (mistag) but mine suppresses (delta < 0).
  std::vector<JetCompCase> cases_pu_mine_corrects_lo, cases_pu_mine_corrects_hi;
  // PU jets where mine gives higher RpT than HGTD (mine makes PU look more HS-like, delta > 0).
  std::vector<JetCompCase> cases_pu_mine_worse_lo, cases_pu_mine_worse_hi;

  long n_total = 0, n_pass_basic = 0, n_hgtd_valid = 0;
  std::vector<HurtJet> hurt_events;
  hurt_events.reserve(30);

  // Per-jet improvement-vs-HGTD categorization.
  //   "Improved": HS jet with rpt_mine > rpt_hgtd, OR PU jet with rpt_mine < rpt_hgtd.
  //   Bucket by cause:
  //     no_hgtd_time : HGTD vertex time invalid (mine has any time to offer)
  //     mine_t_close : both valid, |t_mine − t_truth| < |t_hgtd − t_truth|
  //     other        : both valid, mine's t no closer (mostly tied/HGTD slightly closer)
  struct ImproveStats {
    long hs_total = 0, hs_no_hgtd = 0, hs_mine_t = 0, hs_other = 0;
    long pu_total = 0, pu_no_hgtd = 0, pu_mine_t = 0, pu_other = 0;
    long hs_seen = 0, pu_seen = 0;
  };
  ImproveStats stats_lo, stats_hi;

  while (reader.Next()) {
    ++n_total;

    // ── Require only vertex quality (paper Sec. 3: |z_reco − z_truth| < 2 mm).
    //    No jet-count requirement — every event contributes its jets independently.
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    if (std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]) > MAX_VTX_DZ) continue;
    ++n_pass_basic;

    // ── Track selection: all tracks by z-significance (no eta cut) for the
    //    z-only baseline, matching the paper's ITk-only scenario. ────────────
    std::vector<int> trk_all;
    for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
      double trkPt = branch.trackPt[trk];
      if (trkPt < MIN_TRACK_PT || trkPt > MAX_TRACK_PT) continue;
      if (!branch.trackQuality[trk]) continue;
      if (passTrackVertexAssociation((int)trk, 0, &branch, 2.5))
        trk_all.push_back((int)trk);
    }

    // ── HGTD-acceptance tracks only (used for TRKPTZ clustering). ───────────
    std::vector<int> trk_z = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 2.5);

    // ── TRKPTZ clustering (mine scenario). ───────────────────────────────────
    auto clusters = clusterTracksInTime(
        trk_z, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
        /*sortTracks=*/false, /*calcPurityFlag=*/true);

    // TRKPTZ selection: highest-score cluster (no oracle gate).
    double t_trkptz = 0.0, var_trkptz = 0.0;
    float  trkptz_purity = 0.0f;
    bool   trkptz_ok = false;
    if (!clusters.empty()) {
      auto best      = chooseCluster(clusters, Score::TRKPTZ);
      t_trkptz       = best.values[0];
      var_trkptz     = best.sigmas[0] * best.sigmas[0];
      trkptz_purity  = best.purity;
      trkptz_ok      = true;
    }

    // ── HGTD ntuple vertex time. ─────────────────────────────────────────────
    double t_hgtd         = branch.recoVtxTime[0];
    double var_hgtd       = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
    bool   hgtd_vtx_valid = (branch.recoVtxValid[0] == 1);
    if (hgtd_vtx_valid) ++n_hgtd_valid;

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

    std::vector<int> trk_hgtd   = applyTimeGate(trk_all, t_hgtd,   var_hgtd,   hgtd_vtx_valid);
    std::vector<int> trk_trkptz = applyTimeGate(trk_all, t_trkptz, var_trkptz, trkptz_ok);

    // ── Truth-vertex-t₀ cross-check: gate with the perfect HS vertex time.
    //    var_vtx → 0 (no vertex-time uncertainty); applyTimeGate still folds in
    //    the per-track resolution (var_t), so separation stays finite — a sanity
    //    bound that must NOT yield infinite PU rejection.
    double t_truth = branch.truthVtxTime[0];
    std::vector<int> trk_truth = applyTimeGate(trk_all, t_truth, 0.0, true);

    // ── "HS + in-time + top-K PU" clustering ceiling.  Pool consists of:
    //      • truth-matched HS tracks (valid HGTD time),
    //      • PU tracks whose measured time is within ±PASS_SIGMA of truth HS,
    //      • the top-K highest-pT PU tracks (the score is most sensitive to
    //        these, regardless of timing).
    //    Models realistic cluster purity in the presence of the high-pT PU
    //    tracks that dominate the TRKPTZ sum.  No oracle event selection.
    constexpr int N_TOP_PU = 5;
    const double t_truth_hs = branch.truthVtxTime[0];
    std::unordered_set<int> hs_only_set;
    for (int idx : trk_z) {
      if (branch.trackTimeValid[idx] != 1) continue;
      if (branch.trackToTruthvtx[idx] == 0) {
        hs_only_set.insert(idx);
      } else if (std::abs(branch.trackTime[idx] - t_truth_hs) < PASS_SIGMA) {
        hs_only_set.insert(idx);
      }
    }
    // Collect remaining PU tracks and rank by pT; add top-K.
    std::vector<std::pair<double,int>> pu_by_pt;
    for (int idx : trk_z) {
      if (branch.trackTimeValid[idx] != 1)  continue;
      if (branch.trackToTruthvtx[idx] == 0) continue;
      if (hs_only_set.count(idx))           continue;
      pu_by_pt.emplace_back(branch.trackPt[idx], idx);
    }
    int n_add = std::min((int)pu_by_pt.size(), N_TOP_PU);
    std::partial_sort(pu_by_pt.begin(), pu_by_pt.begin() + n_add, pu_by_pt.end(),
                      [](const auto& a, const auto& b){ return a.first > b.first; });
    for (int i = 0; i < n_add; ++i) hs_only_set.insert(pu_by_pt[i].second);
    std::vector<int> trk_hs_only(hs_only_set.begin(), hs_only_set.end());
    auto clusters_hs = clusterTracksInTime(
        trk_hs_only, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
        /*sortTracks=*/false, /*calcPurityFlag=*/false);
    double t_hsclust = 0.0, var_hsclust = 0.0;
    bool   hsclust_ok = false;
    if (!clusters_hs.empty()) {
      auto best    = chooseCluster(clusters_hs, Score::TRKPTZ);
      t_hsclust    = best.values[0];
      var_hsclust  = best.sigmas[0] * best.sigmas[0];
      hsclust_ok   = true;
    }
    std::vector<int> trk_hsclust = applyTimeGate(trk_all, t_hsclust, var_hsclust, hsclust_ok);

    bool  gate_pure  = trkptz_ok && (trkptz_purity > 0.50f);
    bool  gate_misas = (calcHSTimingPurity(trk_z, &branch) >= 0.80f);

    // Build per-scenario sets once per event for O(1) ghost-index lookup.
    std::unordered_set<int> set_all     (trk_all.begin(),     trk_all.end());
    std::unordered_set<int> set_hgtd    (trk_hgtd.begin(),    trk_hgtd.end());
    std::unordered_set<int> set_trkptz  (trk_trkptz.begin(),  trk_trkptz.end());
    std::unordered_set<int> set_truth   (trk_truth.begin(),   trk_truth.end());
    std::unordered_set<int> set_hsclust (trk_hsclust.begin(), trk_hsclust.end());

    // "Flagging" set: remove tracks whose HGTD time is badly assigned.
    // A track is excluded only when it has a valid time AND its truth-pull
    // |t_reco − t_truth_particle| / trackTimeRes ≥ 3.  Tracks with no valid
    // time (trackTimeValid==0) are kept so z-only information is preserved.
    // "Correcting showers" set: fix bad HGTD tracks by replacing their time with
    //   smeared truth particle time (30 ps Gaussian), then gate all timed tracks
    //   against the truth HS vertex time.  Good tracks use their actual HGTD time.
    // "Tagging showers" set: drop tracks whose HGTD time fails the truth-particle
    //   pull cut; keep all others (including untimed tracks).
    constexpr float FIX_RES = 0.030f;
    float t_hs = branch.truthVtxTime[0];
    std::unordered_set<int> set_corrected, set_flagged;
    for (int idx : trk_all) {
      bool keep_corrected = true, keep_flagged = true;
      if (branch.trackTimeValid[idx] == 1) {
        int part_idx = branch.trackToParticle[idx];
        float tres = branch.trackTimeRes[idx];
        if (tres > 0.0f && part_idx >= 0 &&
            part_idx < (int)branch.particleT.GetSize()) {
          float t_part  = branch.particleT[part_idx];
          float t_track = branch.trackTime[idx];
          bool good = std::abs(t_track - t_part) / tres < 3.0f;
          // Correcting: fix bad tracks with smeared truth time, then gate vs truth vtx.
          float t_eff   = good ? t_track : t_part + (float)gRandom->Gaus(0.0, FIX_RES);
          float sig_eff = good ? tres : FIX_RES;
          if (std::abs(t_eff - t_hs) / sig_eff >= 3.0f) keep_corrected = false;
          // Tagging: drop bad tracks entirely; good tracks gate vs truth vtx.
          if (!good) keep_flagged = false;
          else if (std::abs(t_track - t_hs) / tres >= 3.0f) keep_flagged = false;
        }
      }
      if (keep_corrected) set_corrected.insert(idx);
      if (keep_flagged)   set_flagged.insert(idx);
    }

    // ── Fill jets into pT slices. ─────────────────────────────────────────────
    auto fillJets = [&](std::vector<Scenario>& sv, double pt_lo, double pt_hi,
                        ImproveStats& st) {
      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
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
          double r = computeRpT(&branch, ghost, j_pt, s_set);
          if (isHS) s.h_hs->Fill(r);
          else      s.h_pu->Fill(r);
        };
        fill(sv[0], set_all);
        fill(sv[1], set_hgtd);
        fill(sv[2], set_trkptz);
        fill(sv[3], set_trkptz, gate_pure);
        fill(sv[4], set_trkptz, gate_misas);
        fill(sv[5], set_truth);
        fill(sv[6], set_hsclust);

        // ── Per-jet improvement-vs-HGTD categorization ────────────────────────
        double rpt_mine = computeRpT(&branch, ghost, j_pt, set_trkptz);
        double rpt_hgtd = computeRpT(&branch, ghost, j_pt, set_hgtd);
        bool improved   = isHS ? (rpt_mine > rpt_hgtd) : (rpt_mine < rpt_hgtd);
        (isHS ? st.hs_seen : st.pu_seen)++;
        if (improved) {
          bool no_hgtd = !hgtd_vtx_valid;
          bool mine_t_close = !no_hgtd && trkptz_ok &&
            (std::abs(t_trkptz - t_truth_hs) < std::abs(t_hgtd - t_truth_hs));
          if (isHS) {
            ++st.hs_total;
            if (no_hgtd)           ++st.hs_no_hgtd;
            else if (mine_t_close) ++st.hs_mine_t;
            else                   ++st.hs_other;
          } else {
            ++st.pu_total;
            if (no_hgtd)           ++st.pu_no_hgtd;
            else if (mine_t_close) ++st.pu_mine_t;
            else                   ++st.pu_other;
          }
        }
      }
    };

    fillJets(scen_lo, 30.0, 40.0, stats_lo);
    fillJets(scen_hi, 40.0, 1e9,  stats_hi);

    // ── Flagged vs Correcting scenarios. ─────────────────────────────────────
    // Correcting vs track-level timing comparison.
    //   sv[1] correcting: event-level gate (gate_misas ≥ 80%) + track-level time
    //          gate on passing events.  Jets in bad-timing events are skipped.
    //   sv[2] flagging: NO event gate.  Remove only tracks whose HGTD time is
    //          truth-misassigned (|t_reco−t_truth_part|/σ_t ≥ 3); keep all others
    //          (including untimed tracks).  Shows the ceiling of per-track truth tagging.
    auto fillFlagged = [&](std::vector<Scenario>& sv, double pt_lo, double pt_hi) {
      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
        double j_pt  = branch.topoJetPt[j];
        double j_eta = branch.topoJetEta[j];
        double j_phi = branch.topoJetPhi[j];
        if (j_pt <= pt_lo || j_pt >= pt_hi) continue;
        if (std::abs(j_eta) < JET_ETA_MIN || std::abs(j_eta) > JET_ETA_MAX) continue;
        bool isHS = paperIsHS(j_eta, j_phi);
        bool isPU = paperIsPU(j_eta, j_phi);
        if (!isHS && !isPU) continue;
        const auto& ghost = branch.topoJetGhostTrackIdx[j];
        auto fillRpT = [&](Scenario& s, const std::unordered_set<int>& s_set) {
          double r = computeRpT(&branch, ghost, j_pt, s_set);
          if (isHS) s.h_hs->Fill(r);
          else      s.h_pu->Fill(r);
        };
        fillRpT(sv[0], set_all);                            // z-only baseline: always fill
        if (gate_misas) fillRpT(sv[1], set_corrected);    // correcting: event gate + fixed tracks
        fillRpT(sv[2], set_flagged);                       // tagging: drop truth-misassigned tracks
      }
    };
    fillFlagged(flag_lo, 30.0, 40.0);
    fillFlagged(flag_hi, 40.0, 1e9);

    // ── Mine vs HGTD jet-level comparison (requires both times valid + full cuts). ──
    if (hgtd_vtx_valid && trkptz_ok &&
        branch.passBasicCuts() && branch.passJetPtCut()) {
      std::string fnum = "??????";
      {
        std::string fname = chain.GetCurrentFile()->GetName();
        auto us = fname.rfind('_');
        if (us != std::string::npos) {
          auto dot = fname.find('.', us + 1);
          fnum = fname.substr(us + 1, dot == std::string::npos ? std::string::npos : dot - us - 1);
        }
      }
      long entry = chain.GetTree()->GetReadEntry();

      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
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
        double rpt_hgtd_j = computeRpT(&branch, ghost, j_pt, set_hgtd);
        double rpt_mine_j = computeRpT(&branch, ghost, j_pt, set_trkptz);
        double delta = rpt_mine_j - rpt_hgtd_j;
        if (std::abs(delta) < 0.05) continue;  // skip trivial differences

        JetCompCase c{fnum, entry, j, j_pt, j_eta, rpt_mine_j, rpt_hgtd_j, t_trkptz, delta, isHS};
        if (in_lo) {
          if (delta > 0) insertCase(cases_mine_lo, c);
          else           insertCase(cases_hgtd_lo, c);
        } else {
          if (delta > 0) insertCase(cases_mine_hi, c);
          else           insertCase(cases_hgtd_hi, c);
        }
        // PU jets where HGTD assigns significant RpT (mistag) but mine suppresses.
        if (isPU && delta < 0 && rpt_hgtd_j > 0.1) {
          if (in_lo) insertCase(cases_pu_mine_corrects_lo, c);
          else       insertCase(cases_pu_mine_corrects_hi, c);
        }
        // PU jets where mine gives higher RpT than HGTD (mine makes PU look more HS-like).
        if (isPU && delta > 0 && rpt_mine_j > 0.1) {
          if (in_lo) insertCase(cases_pu_mine_worse_lo, c);
          else       insertCase(cases_pu_mine_worse_hi, c);
        }
      }
    }

    // ── Hurt-HS diagnostic: find HS jets where the mine time gate removed
    //    tracks within ΔR < 0.3, reducing RpT relative to z-only. ──────────
    if (hurt_events.size() < 25) {
      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
        double j_pt  = branch.topoJetPt[j];
        double j_eta = branch.topoJetEta[j];
        double j_phi = branch.topoJetPhi[j];
        if (j_pt <= 30.0 || j_pt >= 40.0) continue;
        if (std::abs(j_eta) < JET_ETA_MIN || std::abs(j_eta) > JET_ETA_MAX) continue;
        if (!paperIsHS(j_eta, j_phi)) continue;

        // Count ghost-associated tracks lost to the TRKPTZ time gate.
        const auto& ghost = branch.topoJetGhostTrackIdx[j];
        int n_lost = 0;
        for (int idx : ghost)
          if (set_all.count(idx) && !set_trkptz.count(idx)) ++n_lost;
        if (n_lost == 0) continue;

        double rpt_z    = computeRpT(&branch, ghost, j_pt, set_all);
        double rpt_mine = computeRpT(&branch, ghost, j_pt, set_trkptz);
        if (rpt_mine >= rpt_z) continue;  // shouldn't happen, but skip no-op cases

        // Extract 6-digit file number from path, e.g. ".../ntuple_000009.root".
        std::string fname = chain.GetCurrentFile()->GetName();
        std::string fnum  = "??????";
        auto us = fname.rfind('_');
        if (us != std::string::npos) {
          // Find the first '.' after the underscore to handle "000009.SuperNtuple.root"
          auto dot = fname.find('.', us + 1);
          fnum = fname.substr(us + 1, dot == std::string::npos ? std::string::npos : dot - us - 1);
        }
        long entry = chain.GetTree()->GetReadEntry();

        hurt_events.push_back({fnum, entry, j_pt, j_eta, j_phi,
                                rpt_z, rpt_mine, n_lost});
      }
    }
  }

  // ── Output paths ─────────────────────────────────────────────────────────────
  boost::filesystem::create_directories("../figs/rpt_plots");
  const TString out_pdf = "../figs/rpt_plots/rpt_v4.pdf";
  const std::string out_dir_mine = "./eventdisplays/pu_fix_mine_over_hgtd";
  const std::string out_dir_hgtd = "./eventdisplays/pu_worse_hgtd_over_mine";

  TCanvas* canvas = new TCanvas("canvas", "RpT v3", 800, 700);
  canvas->Print(out_pdf + "[");

  auto drawLabels = [](const char* extra = nullptr) {
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82);
    if (extra) {
      TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
      t.DrawLatex(0.18, 0.76, extra);
    }
  };
  const char* lbl_lo = "30 < p_{T}^{jet} < 40 GeV, 2.4 < |#eta| < 3.8";
  const char* lbl_hi = "p_{T}^{jet} > 40 GeV, 2.4 < |#eta| < 3.8";

  auto styleScen = [](std::vector<Scenario>& sv) {
    for (auto& s : sv) {
      if (!s.h_hs || !s.h_pu) continue;  // skip reference entries with no histograms
      for (auto* h : {s.h_hs, s.h_pu}) {
        h->GetXaxis()->SetRangeUser(0.0, 1.5);
        h->GetXaxis()->SetNdivisions(515);
        h->SetLineWidth(2);
      }
      s.h_hs->SetLineColor(s.color);
      s.h_pu->SetLineColor(s.color);
      s.h_pu->SetLineStyle(2);
    }
  };
  styleScen(scen_lo);
  styleScen(scen_hi);
  styleScen(flag_lo);
  styleScen(flag_hi);

  const double roc_xmin = 0.8, roc_xmax = 1.0;
  std::vector<TGraph*> rocs_lo, rocs_hi;
  for (auto& s : scen_lo) rocs_lo.push_back(generate_roc(s.h_pu, s.h_hs));
  for (auto& s : scen_hi) rocs_hi.push_back(generate_roc(s.h_pu, s.h_hs));

  std::vector<TGraph*> rocs_flo, rocs_fhi;
  // Generate only for the 3 scenarios that have histograms; mine_ref (index 3) is appended below.
  for (size_t i = 0; i < 3; ++i) rocs_flo.push_back(generate_roc(flag_lo[i].h_pu, flag_lo[i].h_hs));
  for (size_t i = 0; i < 3; ++i) rocs_fhi.push_back(generate_roc(flag_hi[i].h_pu, flag_hi[i].h_hs));
  // Baseline TRKPTZ: reuse already-computed ROC from main comparison.
  rocs_flo.push_back(generate_roc(scen_lo[2].h_pu, scen_lo[2].h_hs));
  rocs_fhi.push_back(generate_roc(scen_hi[2].h_pu, scen_hi[2].h_hs));

  auto styleRoc = [&](TGraph* g, Color_t col) {
    g->SetTitle("R_{pT} Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->SetLineWidth(2);
    g->GetXaxis()->SetLimits(roc_xmin, roc_xmax);
    g->GetXaxis()->SetNdivisions(810);
    g->SetMinimum(1.0);
  };
  for (size_t i = 0; i < rocs_lo.size(); ++i) {
    styleRoc(rocs_lo[i], scen_lo[i].color);
    styleRoc(rocs_hi[i], scen_hi[i].color);
  }
  for (auto* g : rocs_lo) g->SetMaximum(300);
  for (auto* g : rocs_hi) g->SetMaximum(500);

  // ROC page with ratio-to-zonly panel.
  auto drawRocWithRatio = [&](std::vector<TGraph*>& gs,
                               std::vector<Scenario>& sc,
                               double ymax, double ratio_ymax,
                               const char* extra_label) {
    canvas->Clear();
    canvas->SetLogy(false);

    TPad* pad_main = new TPad("pad_main", "", 0.0, 0.30, 1.0, 1.0);
    pad_main->SetBottomMargin(0.02);
    pad_main->SetTopMargin(0.07);
    pad_main->SetLeftMargin(0.16);
    pad_main->SetRightMargin(0.05);
    pad_main->Draw();
    pad_main->cd();

    gs[0]->SetMaximum(ymax);
    gs[0]->GetXaxis()->SetLabelSize(0);
    gs[0]->GetXaxis()->SetTitleSize(0);
    gs[0]->Draw("ALP");
    for (size_t i = 1; i < gs.size(); ++i) gs[i]->Draw("LP SAME");

    TLegend* L = new TLegend(0.58, 0.60, 0.93, 0.88);
    StyleLegend(L);
    for (size_t i = 0; i < gs.size(); ++i) L->AddEntry(gs[i], sc[i].legend.c_str());
    L->Draw();

    {
      TLatex t; t.SetNDC(); t.SetTextColor(1);
      t.SetTextFont(72); t.SetTextSize(0.05); t.DrawLatex(0.18, 0.88, "ATLAS");
      t.SetTextFont(42); t.SetTextSize(0.05); t.DrawLatex(0.28, 0.88, "Simulation Internal");
      t.SetTextSize(0.044); t.DrawLatex(0.18, 0.80, "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.");
      if (extra_label) { t.SetTextSize(0.044); t.DrawLatex(0.18, 0.72, extra_label); }
    }

    canvas->cd();

    TPad* pad_ratio = new TPad("pad_ratio", "", 0.0, 0.0, 1.0, 0.30);
    pad_ratio->SetTopMargin(0.02);
    pad_ratio->SetBottomMargin(0.38);
    pad_ratio->SetLeftMargin(0.16);
    pad_ratio->SetRightMargin(0.05);
    pad_ratio->Draw();
    pad_ratio->cd();

    TGraph* ref = gs[0];
    std::vector<TGraph*> ratios;
    for (size_t i = 0; i < gs.size(); ++i) {
      TGraph* g = gs[i];
      int n = g->GetN();
      std::vector<double> rx, ry;
      for (int j = 0; j < n; ++j) {
        double x = g->GetX()[j];
        if (x < roc_xmin || x > roc_xmax) continue;
        double ref_y = ref->Eval(x);
        if (ref_y <= 0) continue;
        rx.push_back(x);
        ry.push_back(g->GetY()[j] / ref_y);
      }
      TGraph* r = new TGraph(rx.size(), rx.data(), ry.data());
      r->SetLineColor(g->GetLineColor());
      r->SetLineWidth(2);
      ratios.push_back(r);
    }

    TH1F* frame = pad_ratio->DrawFrame(roc_xmin, 0.5, roc_xmax, ratio_ymax,
        ";Hard Scatter Efficiency;Ratio to z-only");
    frame->GetXaxis()->SetNdivisions(810);
    frame->GetXaxis()->SetLabelSize(0.13);
    frame->GetXaxis()->SetTitleSize(0.13);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetLabelSize(0.11);
    frame->GetYaxis()->SetTitleSize(0.11);
    frame->GetYaxis()->SetTitleOffset(0.65);
    frame->GetYaxis()->SetNdivisions(504);

    TLine* line1 = new TLine(roc_xmin, 1.0, roc_xmax, 1.0);
    line1->SetLineColor(kGray + 1);
    line1->SetLineStyle(2);
    line1->Draw();
    for (auto* r : ratios) r->Draw("L SAME");

    canvas->cd();
    canvas->Print(out_pdf);
  };

  // (1) ROC — 30–50 GeV.
  drawRocWithRatio(rocs_lo, scen_lo, 300.0, 4.0, lbl_lo);

  // (2) ROC — 50–70 GeV.
  drawRocWithRatio(rocs_hi, scen_hi, 500.0, 5.0, lbl_hi);

  // (3–5) Per-scenario HS vs PU, 30–50 GeV slice, log-y.
  canvas->Clear();
  canvas->SetLogy(true);
  for (auto& s : scen_lo) {
    TLegend* L = new TLegend(0.60, 0.75, 0.92, 0.88);
    StyleLegend(L);
    L->AddEntry(s.h_hs, "Hard Scatter");
    L->AddEntry(s.h_pu, "Pile-Up");
    double ymax = 50.0 * std::max(s.h_hs->GetMaximum(), s.h_pu->GetMaximum());
    s.h_pu->SetMaximum(ymax);
    s.h_pu->SetTitle((std::string("R_{pT}: ") + s.legend).c_str());
    s.h_pu->Draw("HIST");
    s.h_hs->Draw("HIST SAME");
    drawLabels(lbl_lo);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // (6) All-scenario HS overlay — 30–50 GeV.
  {
    TLegend* L = new TLegend(0.55, 0.68, 0.92, 0.88);
    StyleLegend(L);
    double ymax = 0;
    for (auto& s : scen_lo) ymax = std::max(ymax, s.h_hs->GetMaximum());
    scen_lo[0].h_hs->SetMaximum(50.0 * ymax);
    scen_lo[0].h_hs->SetTitle("Hard Scatter R_{pT} — all scenarios");
    for (size_t i = 0; i < scen_lo.size(); ++i) {
      scen_lo[i].h_hs->SetLineStyle(1);
      scen_lo[i].h_hs->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(scen_lo[i].h_hs, scen_lo[i].legend.c_str());
    }
    drawLabels(lbl_lo);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // (7) All-scenario PU overlay — 30–50 GeV.
  {
    TLegend* L = new TLegend(0.55, 0.68, 0.92, 0.88);
    StyleLegend(L);
    double ymax = 0;
    for (auto& s : scen_lo) ymax = std::max(ymax, s.h_pu->GetMaximum());
    scen_lo[0].h_pu->SetMaximum(50.0 * ymax);
    scen_lo[0].h_pu->SetTitle("Pile-Up R_{pT} — all scenarios");
    for (size_t i = 0; i < scen_lo.size(); ++i) {
      scen_lo[i].h_pu->SetLineStyle(1);
      scen_lo[i].h_pu->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(scen_lo[i].h_pu, scen_lo[i].legend.c_str());
    }
    drawLabels(lbl_lo);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // ── Flagging vs Correcting pages ─────────────────────────────────────────────
  // Style the flagged ROCs: correcting = solid, flagging = dashed.
  for (size_t i = 0; i < rocs_flo.size(); ++i) {
    styleRoc(rocs_flo[i], flag_lo[i].color);
    styleRoc(rocs_fhi[i], flag_hi[i].color);
    rocs_flo[i]->SetMaximum(300.0);
    rocs_fhi[i]->SetMaximum(500.0);
  }
  // correcting = solid, flagging = dashed, baseline TRKPTZ = dotted
  rocs_flo[1]->SetLineStyle(1);  rocs_fhi[1]->SetLineStyle(1);
  rocs_flo[2]->SetLineStyle(2);  rocs_fhi[2]->SetLineStyle(2);
  rocs_flo[3]->SetLineStyle(1);  rocs_fhi[3]->SetLineStyle(1);

  // Flagging vs Correcting ROC — 30–50 GeV.
  drawRocWithRatio(rocs_flo, flag_lo, 300.0, 4.0, lbl_lo);

  // Flagging vs Correcting ROC — 50–70 GeV.
  drawRocWithRatio(rocs_fhi, flag_hi, 500.0, 5.0, lbl_hi);

  canvas->Print(out_pdf + "]");

  // ── Console summary ───────────────────────────────────────────────────────────
  std::cout << "\n=== EVENT COUNTS ===\n";
  std::cout << "Total events        : " << n_total      << '\n';
  std::cout << "Pass basic cuts     : " << n_pass_basic << '\n';
  std::cout << "HGTD vtx valid      : " << n_hgtd_valid << " / " << n_pass_basic
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * n_hgtd_valid / n_pass_basic) << "%)\n";

  std::cout << "\n=== ENTRIES PER SCENARIO (30-40 GeV) ===\n";
  for (auto& s : scen_lo)
    std::cout << "  " << s.name
              << "  HS: " << (long)s.h_hs->Integral()
              << "  PU: " << (long)s.h_pu->Integral() << '\n';

  std::cout << "\n=== ENTRIES PER SCENARIO (>40 GeV) ===\n";
  for (auto& s : scen_hi)
    std::cout << "  " << s.name
              << "  HS: " << (long)s.h_hs->Integral()
              << "  PU: " << (long)s.h_pu->Integral() << '\n';

  // Per-jet improvement-vs-HGTD breakdown.
  auto printImproveBlock = [](const char* hdr, const ImproveStats& s) {
    auto pct = [](long n, long d) -> double {
      return d > 0 ? 100.0 * (double)n / (double)d : 0.0;
    };
    std::cout << "\n=== MINE vs HGTD IMPROVEMENT BREAKDOWN " << hdr << " ===\n";
    std::cout << "  HS jets seen           : " << s.hs_seen  << "\n";
    std::cout << "  HS improved over HGTD  : " << s.hs_total
              << " (" << std::fixed << std::setprecision(1)
              << pct(s.hs_total, s.hs_seen) << "% of HS jets)\n";
    if (s.hs_total > 0) {
      std::cout << "    1) HGTD has no time  : " << s.hs_no_hgtd
                << "  (" << pct(s.hs_no_hgtd, s.hs_total) << "% of improved)\n";
      std::cout << "    2) Mine t closer     : " << s.hs_mine_t
                << "  (" << pct(s.hs_mine_t, s.hs_total) << "% of improved)\n";
      std::cout << "    3) Neither           : " << s.hs_other
                << "  (" << pct(s.hs_other, s.hs_total) << "% of improved)\n";
    }
    std::cout << "  PU jets seen           : " << s.pu_seen  << "\n";
    std::cout << "  PU improved over HGTD  : " << s.pu_total
              << " (" << pct(s.pu_total, s.pu_seen) << "% of PU jets)\n";
    if (s.pu_total > 0) {
      std::cout << "    1) HGTD has no time  : " << s.pu_no_hgtd
                << "  (" << pct(s.pu_no_hgtd, s.pu_total) << "% of improved)\n";
      std::cout << "    2) Mine t closer     : " << s.pu_mine_t
                << "  (" << pct(s.pu_mine_t, s.pu_total) << "% of improved)\n";
      std::cout << "    3) Neither           : " << s.pu_other
                << "  (" << pct(s.pu_other, s.pu_total) << "% of improved)\n";
    }
  };
  printImproveBlock("(30-40 GeV)", stats_lo);
  printImproveBlock("(>40 GeV)",   stats_hi);

  std::cout << "\nWrote " << out_pdf << "\n";

  // ── Hurt-HS events: HS jets where the mine time gate removed tracks ────────
  std::cout << "\n=== TIMING-HURT HS JETS (30–50 GeV, mine scenario) ===\n";
  std::cout << "  HS jets where time gate removed ≥1 track within ΔR<0.3, lowering RpT.\n\n";
  for (auto& h : hurt_events) {
    std::printf("  jet pT=%.1f  eta=%.2f  phi=%.2f  RpT: %.3f→%.3f  tracks_lost=%d\n",
                h.j_pt, h.j_eta, h.j_phi, h.rpt_z, h.rpt_mine, h.n_lost);
    std::printf("  python event_display.py --file_num %s --event_num %ld --extra_time 0.00\n\n",
                h.file_num.c_str(), h.entry);
  }
  if (hurt_events.empty())
    std::cout << "  (none found)\n";

  // ── Mine vs HGTD comparison event display commands ───────────────────────────
  // Prints a summary to stdout and writes a run_displays.sh into out_dir (if set).
  auto printCases = [&](const char* title,
                        const std::vector<JetCompCase>& cases,
                        const std::string& out_dir = "") {
    std::cout << "\n" << title << ":\n";
    if (cases.empty()) { std::cout << "  (none found)\n"; return; }
    for (const auto& c : cases) {
      std::printf("  jet pT=%.1f GeV  eta=%.2f  %s  RpT: mine=%.3f  hgtd=%.3f  Δ=%.3f\n",
                  c.j_pt, c.j_eta, c.isHS ? "HS" : "PU",
                  c.rpt_mine, c.rpt_hgtd, c.delta);
      std::printf("  python3 event_display.py --file_num %s --event_num %ld"
                  " --extra_time %.2f --jet_idx %d"
                  " --jet_label %s --rpt_hgtd %.3f --rpt_mine %.3f"
                  " --output_dir %s\n\n",
                  c.file_num.c_str(), c.entry, c.t_trkptz, c.jet_idx,
                  c.isHS ? "HS" : "PU", c.rpt_hgtd, c.rpt_mine,
                  out_dir.empty() ? "event_displays" : out_dir.c_str());
    }
  };

  std::cout << "\n=== MINE vs HGTD: JET-LEVEL COMPARISON ===\n";
  printCases("CASE 1 — Mine improves (30–40 GeV, mine RpT > HGTD RpT)", cases_mine_lo);
  printCases("CASE 2 — Mine improves (>40 GeV, mine RpT > HGTD RpT)", cases_mine_hi);
  printCases("CASE 3 — HGTD better  (30–40 GeV, HGTD RpT > mine RpT)", cases_hgtd_lo);
  printCases("CASE 4 — HGTD better  (>40 GeV, HGTD RpT > mine RpT)", cases_hgtd_hi);

  std::cout << "\n=== PU MISTAG CORRECTION: HGTD mistags PU as HS, mine corrects ===\n";
  std::cout << "  (PU jets with rpt_hgtd > 0.1 where mine gives lower RpT)\n";
  printCases("CASE 5 — Mine corrects PU mistag (30–40 GeV)", cases_pu_mine_corrects_lo, out_dir_mine);
  printCases("CASE 6 — Mine corrects PU mistag (>40 GeV)",   cases_pu_mine_corrects_hi, out_dir_mine);

  std::cout << "\n=== PU MISTAG WORSENING: Mine pushes PU jet closer to HS ===\n";
  std::cout << "  (PU jets with rpt_mine > 0.1 where mine gives higher RpT than HGTD)\n";
  printCases("CASE 7 — Mine worsens PU mistag (30–40 GeV)", cases_pu_mine_worse_lo, out_dir_hgtd);
  printCases("CASE 8 — Mine worsens PU mistag (>40 GeV)",   cases_pu_mine_worse_hi, out_dir_hgtd);

  return 0;
}
