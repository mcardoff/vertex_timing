// -----------------------------------------------------------------------------
// tdr_common.h -- shared pieces of the TDR-structured programs (util/tdr_pujet.cxx,
// util/tdr_vbf.cxx), which present the HGTD timing results the way ATLAS-TDR-031
// does: Sec. 3.3.3 (pileup-jet suppression with R_pT and timing) and Sec. 3.4.1
// (VBF H->inv against Z+jets).
//
// Deliberately SEPARATE from rpt_v5 and vbs_time_veto: their selections and
// conventions differ, and this keeps the TDR presentation free to follow the TDR.
// The ingredients are COPIED from them, so the numbers stay comparable:
//   R_pT tracks   util/rpt_v5_hist.cxx: ghost-associated, 1 < pT < 30 GeV, quality,
//                 within dR < 0.2 of the axis, associated to PV0 by the forward
//                 parameterisation (getNewDzpara x 1.4) or, for the central list,
//                 a z0 significance < 5
//   t0s           util/vbs_time_veto.cxx (rpt_v5's, verbatim): HGTD = Athena's
//                 RecoVtx_time; TRKPTZ and WAVeS from one time clustering of the
//                 associated tracks; truth = t_HS (reference only)
//   t0 gate       rpt_v5's applyTimeGate: keep a timed track when
//                 |t_trk - t0| / sqrt(infl^2 var_t0 + var_trk) < gate; untimed tracks
//                 always kept. The TDR's gate is 2 sigma (rpt_v5 uses 3).
// TDR timing taggers (Sec. 3.3.3):
//   t0          R_pT recomputed after the gate
//   self        jets with >= 2 timed R_pT tracks: the timed ones clustered in time
//               (2 sigma, the jet-time clustering), R_pT per cluster counting the
//               cluster's tracks plus the untimed ones, the jet's R_pT the maximum;
//               < 2 timed tracks: unchanged
//   combined    the t0 version where that t0 exists, self-tagging where it does not
//               (the TDR's "when no t0 is found"); only Athena's t0 can be absent
// -----------------------------------------------------------------------------
#pragma once

#include <TVector2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "clustering_constants.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "rpt_v5_common.h"

namespace tdr {

using namespace MyUtl;

// ── Conventions (TDR Sec. 3.3.3 / 3.4.1) ─────────────────────────────────────
constexpr double JET_PT_MIN      = 30.0;   // jets entering the taggers
constexpr double FWD_ETA_MIN     = 2.4;    // TDR forward window 2.4 < |eta| < 4.0
constexpr double FWD_ETA_MAX     = 4.0;
constexpr double GATE_SIGMA_TDR  = 2.0;    // the TDR's t0 gate
constexpr double SELF_DIST_CUT   = 2.0;    // self-tagging time clustering (the jet-time clustering's)

// ── R_pT ingredients, copied from util/rpt_v5_hist.cxx ───────────────────────
constexpr double DZ0_PARA_SCALE   = 1.4;   // forward association: |dz| / getNewDzpara < 1.4
constexpr double CENTRAL_Z_SIGNIF = 5.0;   // central association: |dz| / sigma_z0 < 5
constexpr double RPT_TRACK_JET_DR = 0.2;   // R_pT cone
constexpr double T0_ASSOC_NSIGMA  = 2.5;   // t0 clustering input (rpt_v5 / vbs_time_veto)

// t0 sources, in this order everywhere
constexpr int N_T0 = 4;
enum T0Idx { T0_HGTD = 0, T0_TRKPTZ = 1, T0_WAVES = 2, T0_TRUTH = 3 };
inline const char* const T0_NAMES[N_T0] = {"hgtd", "trkptz", "waves", "truth"};

struct T0Set {
  std::array<double, N_T0> t{}, sig{}, infl{};
  std::array<bool, N_T0>   ok{};
};

inline double dR(double e1, double p1, double e2, double p2) {
  return std::hypot(e1 - e2, TVector2::Phi_mpi_pi(p1 - p2));
}

// The event's t0s -- util/vbs_time_veto.cxx's block, verbatim, without TZP.
inline T0Set computeT0s(BranchPointerWrapper& b, const Inflation& INFL) {
  T0Set s;
  s.t[T0_HGTD]   = b.recoVtxTime[0];
  s.sig[T0_HGTD] = b.recoVtxTimeRes[0];
  s.ok[T0_HGTD]  = (b.recoVtxValid[0] == 1);
  const std::vector<int> trk_z = getAssociatedTracks(&b, MIN_TRACK_PT, MAX_TRACK_PT, T0_ASSOC_NSIGMA);
  const auto clusters = clusterTracksInTime(
      trk_z, &b, DIST_CUT_CONE,
      /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);
  if (!clusters.empty()) {
    auto bT = chooseCluster(clusters, Score::TRKPTZ);
    s.t[T0_TRKPTZ] = bT.calculateTime(Score::TRKPTZ, &b);
    s.sig[T0_TRKPTZ] = bT.sigmas[0];  s.ok[T0_TRKPTZ] = true;
    auto bW = chooseCluster(clusters, Score::WAVES);
    s.t[T0_WAVES] = bW.calculateTime(Score::WAVES, &b);       // in-jet refined
    s.sig[T0_WAVES] = bW.sigmas[0];  s.ok[T0_WAVES] = true;
  }
  s.t[T0_TRUTH] = b.truthVtxTime[0];  s.sig[T0_TRUTH] = 0.0;  s.ok[T0_TRUTH] = true;
  s.infl = {INFL.hgtd, INFL.trkptz, INFL.waves, 1.0};
  return s;
}

// PV0-associated track lists: forward (rpt_v5's dz parameterisation) and central
// (z0 significance) -- the same split rpt_v5 makes, chosen per jet by its |eta|.
struct TrackLists {
  std::unordered_set<int> fwd, cen;
};
inline TrackLists rptTrackLists(BranchPointerWrapper& b) {
  TrackLists L;
  const double vtxZ = b.recoVtxZ[0];
  for (size_t trk = 0; trk < b.trackZ0.GetSize(); ++trk) {
    const double pt = b.trackPt[trk];
    if (pt < MIN_TRACK_PT || pt > MAX_TRACK_PT || !b.trackQuality[trk]) continue;
    const double dz = std::abs(b.trackZ0[trk] - vtxZ);
    if (dz / MyUtl::getNewDzpara(b.trackEta[trk], pt) <= DZ0_PARA_SCALE) L.fwd.insert((int)trk);
    const double vz = b.trackVarZ0[trk];
    if (vz > 0 && dz / std::sqrt(vz) <= CENTRAL_Z_SIGNIF) L.cen.insert((int)trk);
  }
  return L;
}

// The jet's R_pT tracks: ghost-associated, in the PV0 list, within the R_pT cone.
inline std::vector<int> jetRptTracks(BranchPointerWrapper& b, int j, const std::unordered_set<int>& assoc) {
  std::vector<int> out;
  const double je = b.topoJetEta[j], jp = b.topoJetPhi[j];
  for (int idx : b.topoJetGhostTrackIdx[j]) {
    if (!assoc.count(idx)) continue;
    if (dR(je, jp, b.trackEta[idx], b.trackPhi[idx]) > RPT_TRACK_JET_DR) continue;
    out.push_back(idx);
  }
  return out;
}

inline double sumPt(BranchPointerWrapper& b, const std::vector<int>& trk) {
  double s = 0.0;
  for (int i : trk) s += b.trackPt[i];
  return s;
}

// rpt_v5's applyTimeGate
inline std::vector<int> timeGate(BranchPointerWrapper& b, const std::vector<int>& base,
                                 double t0, double sig0, bool ok, double infl, double gate) {
  std::vector<int> out;
  out.reserve(base.size());
  for (int idx : base) {
    if (!ok || b.trackTimeValid[idx] != 1) { out.push_back(idx); continue; }
    const double dt = b.trackTime[idx] - t0;
    const double vt = b.trackTimeRes[idx] * b.trackTimeRes[idx];
    if (std::abs(dt) / std::sqrt(infl * infl * sig0 * sig0 + vt) < gate) out.push_back(idx);
  }
  return out;
}

// Self-tagging R_pT numerator: max over in-time clusters of the timed R_pT tracks of
// (cluster sum pT + untimed sum pT); the full sum when fewer than 2 tracks are timed.
inline double selfTagSumPt(BranchPointerWrapper& b, const std::vector<int>& trk, int* nTimedOut = nullptr) {
  std::vector<int> timed;
  double untimed = 0.0;
  for (int i : trk) {
    if (b.trackTimeValid[i] == 1 && b.trackTimeRes[i] > 0.0) timed.push_back(i);
    else untimed += b.trackPt[i];
  }
  if (nTimedOut) *nTimedOut = (int)timed.size();
  if (timed.size() < 2) return sumPt(b, trk);
  static const std::unordered_map<int, double> kNone;
  std::vector<Cluster> cl = makeSimpleClusters(timed, &b, /*useSmearedTimes=*/false, kNone, kNone,
                                               /*checkTimeValid=*/true, /*usez0=*/false);
  doIterativeClustering(&cl, SELF_DIST_CUT);
  double best = 0.0;
  for (const Cluster& c : cl) {
    double s = 0.0;
    for (int i : c.trackIndices) s += b.trackPt[i];
    best = std::max(best, s);
  }
  return best + untimed;
}

// Every tagger's R_pT for one jet.
struct JetRpt {
  int    ntrk = 0, ntimed = 0;          // R_pT tracks, and how many carry an HGTD time
  double itk = 0.0;                     // ITk only (no timing)
  std::array<double, N_T0> t0{};        // t0-gated, per t0 source
  double self = 0.0;                    // self-tagging
  double comb = 0.0;                    // Athena t0 where it exists, self-tagging otherwise
};
inline JetRpt jetRpt(BranchPointerWrapper& b, int j, const TrackLists& L, const T0Set& T, double gate) {
  JetRpt r;
  const double pt = b.topoJetPt[j];
  const bool fwd = std::abs((double)b.topoJetEta[j]) >= FWD_ETA_MIN;
  const std::vector<int> trk = jetRptTracks(b, j, fwd ? L.fwd : L.cen);
  r.ntrk = (int)trk.size();
  r.itk = sumPt(b, trk) / pt;
  for (int k = 0; k < N_T0; ++k)
    r.t0[k] = sumPt(b, timeGate(b, trk, T.t[k], T.sig[k], T.ok[k], T.infl[k], gate)) / pt;
  r.self = selfTagSumPt(b, trk, &r.ntimed) / pt;
  r.comb = T.ok[T0_HGTD] ? r.t0[T0_HGTD] : r.self;
  return r;
}

}  // namespace tdr
