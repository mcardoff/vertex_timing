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
//                 associated to PV0 by the forward parameterisation (getNewDzpara x 1.4,
//                 = 2.4 sigma of the HS z0 resolution: the 2022 PUB note's "2.5 times")
//                 or, for the central list, a z0 significance < 5.
//                 Matching (--rpt-match): ghost association only by DEFAULT -- the TDR's
//                 (its ref. [12]) and the 2022 PUB note's definition; "dr02" adds rpt_v5's
//                 dR < 0.2 to the axis (from the PUB note's macro, not its text), tagged
//                 dr02_ and kept for the presentation's backup. The cone strengthens ITk
//                 only and takes the work timing would do (notes/tdr_presentation.md).
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
// Track-time cleaning (ON by default since 2026-11-03, the user's call: "fold this into
// our presentation"): the ntuples' track times are NOT hit-cleaned, and a third of the
// timed forward tracks carry a misassigned time. A timed track counts as untimed unless
// it passes the TDR's two reconstruction-level cuts (Sec. 3.2.1): its last ITk hit on the
// HGTD (Track_hitInLastLayer), and >= 2 HGTD hits at 3.5 < |eta| < 3.9. The third
// TDR step (outlier removal among a track's hits) looks already applied upstream:
// Track_hasValidTime is 0 for 14.5% of 2-hit tracks, ~1% of 3-hit, 0.1% of 4-hit and
// never for 1-hit ones, with Track_time = 0 (inferred from the data, not the Athena
// code). Modes (--clean-times=all|t0|off):
//   all  (default, untagged)  every track time used here: the TRKPTZ/WAVeS t0 inputs,
//                             the gate and self-tagging -- the TDR's own reconstruction
//   t0   (tag cleant0_)       the t0 inputs only; the jets' tracks keep the ntuple times
//   off  (tag rawtimes_)      the ntuple as is (what rpt_v5 and every other study use)
// Athena's t0 (RecoVtx_time) comes from upstream and cannot be recomputed. Needs
// Track_hitInLastLayer, which the skims drop: read the raw ntuples (--ntuple-dir=).
// -----------------------------------------------------------------------------
#pragma once

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TVector2.h>

#include <memory>

#include <algorithm>
#include <array>
#include <limits>
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

// t0 sources, in this order everywhere. The last three (2026-11-13, todo vt-23j) are the WAVeS t0 with its failures
// repaired by truth, exactly as clustering_hist's rows of the talk's "perfect timing / perfect selection" slide
// (src/failure_fixes.h's applyFixes on the same clustering input):
//   wfix_sel   Fix::SELECT             the hard scatter's cluster chosen (WFIX_SELECT, "perfect selection")
//   wfix_time  Fix::MISTIMED           the HS tracks >= 3 sigma from t_HS dropped, re-clustered (WFIX_MISTIMED, "perfect timing")
//   wfix_both  Fix::MISTIMED | SELECT  both (WFIX_C3)
// They use WAVeS's inflation; no t0 (ok = false) when no cluster qualifies.
// The next five (2026-11-14, todo vt-5t0) are two reco-level SELECTION levers against the > 50 GeV capture (a hard
// forward pileup jet's collision wins WAVeS's choice, notes/gp_talk.md): WAVeS's choice re-made on the same clusters
// with a modified score, then WAVeS's in-jet time and sigma for the new pick -- so they differ from WAVeS by the
// selection alone (wavesLever below):
//   waves_a0 / waves_a025 / waves_a05   pT_jet^alpha in place of WAVeS's pT_jet, alpha = 0 / 0.25 / 0.5 (WAVeS: 1)
//   waves_pv0   WAVeS x the share of the cluster's pT the vertex fit gave PV0 (Track_recoVtx_idx)
//   waves_npu   WAVeS x (1 - the share it gave ANOTHER reco vertex): tracks outside every fit count as PV0's
// WAVeS's inflation. No t0 when no cluster qualifies, and for the vertex-fit pair when the input lacks Track_recoVtx_idx.
// The last four (2026-11-14) size the room above 50 GeV, where perfect selection + timing (wfix_both = WFIX_C3) reach
// ITk x1.58 at 85% and the truth t0 x2.88 (results/wfix_downstream.md):
//   wfix_c4 / wfix_c5 / wfix_c6   the rest of clustering_hist's cumulative fix chain on WAVeS: + misclustering (the
//                                 chosen cluster timed from its HS tracks alone), + HS track efficiency (untimed or
//                                 cleaned HS tracks in HGTD given Gaus(t_HS, 30 ps)), + HS time resolution (timed HS
//                                 tracks at Gaus(t_HS, 1 ps)); WAVeS's inflation; misassigned HS times untimed in the
//                                 gate, as wfix_time / wfix_both (fixesTimes)
//   truth_wsig                    t_HS exactly, but gated with WAVeS's sigma_t0 and inflation: the truth row's gate is
//                                 narrower (sigma_t0 = 0), so truth - truth_wsig is the part of its gain that no real
//                                 t0 can have; ok wherever WAVeS has a t0 (else sigma_t0 = 0, as truth)
constexpr int N_T0 = 16;
enum T0Idx { T0_HGTD = 0, T0_TRKPTZ = 1, T0_WAVES = 2, T0_TRUTH = 3, T0_WFIX_SEL = 4, T0_WFIX_TIME = 5, T0_WFIX_BOTH = 6,
             T0_WA0 = 7, T0_WA025 = 8, T0_WA05 = 9, T0_WPV0 = 10, T0_WNPU = 11,
             T0_WFIX_C4 = 12, T0_WFIX_C5 = 13, T0_WFIX_C6 = 14, T0_TRUTH_WSIG = 15 };
inline const char* const T0_NAMES[N_T0] = {"hgtd", "trkptz", "waves", "truth", "wfix_sel", "wfix_time", "wfix_both",
                                           "waves_a0", "waves_a025", "waves_a05", "waves_pv0", "waves_npu",
                                           "wfix_c4", "wfix_c5", "wfix_c6", "truth_wsig"};

struct T0Set {
  std::array<double, N_T0> t{}, sig{}, infl{};
  std::array<bool, N_T0>   ok{};
};

inline double dR(double e1, double p1, double e2, double p2) {
  return std::hypot(e1 - e2, TVector2::Phi_mpi_pi(p1 - p2));
}

// ── Track-time cleaning ───────────────────────────────────────────────────────
enum class CleanMode { OFF, T0, ALL };
inline const char* cleanModeName(CleanMode m) { return m == CleanMode::OFF ? "off" : m == CleanMode::T0 ? "t0" : "all"; }
struct TimeCleaning {
  CleanMode mode = CleanMode::OFF;
  std::unique_ptr<TTreeReaderValue<std::vector<bool>>> lastLayer;   // Track_hitInLastLayer
  bool t0() const { return mode != CleanMode::OFF; }                // clean the t0 inputs?
  bool jets() const { return mode == CleanMode::ALL; }              // clean the jets' tracks?
};
// --clean-times=all|t0|off (default all). Call before the output path is built: it tags the file.
inline CleanMode cleaningMode(int argc, char** argv) {
  CleanMode m = CleanMode::ALL;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--clean-times", 0) != 0) continue;
    const std::string v = a.size() > 14 && a[13] == '=' ? a.substr(14) : "all";
    if (v == "all") m = CleanMode::ALL;
    else if (v == "t0") m = CleanMode::T0;
    else if (v == "off") m = CleanMode::OFF;
    else { std::cerr << "[tdr] --clean-times=" << v << ": expected all, t0 or off\n"; std::exit(1); }
  }
  const char* tag = m == CleanMode::OFF ? "rawtimes" : m == CleanMode::T0 ? "cleant0" : nullptr;
  if (tag) SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  return m;
}
// bind the branch the cleaning needs (before the first reader.Next())
inline TimeCleaning setupCleaning(TTreeReader& reader, CleanMode m) {
  TimeCleaning C;
  C.mode = m;
  if (m == CleanMode::OFF) return C;
  if (!hasBranch("Track_hitInLastLayer")) {
    std::cerr << "[tdr] track-time cleaning needs Track_hitInLastLayer, which this input lacks"
                 " (the skims drop it: read the raw ntuples with --ntuple-dir=, or run --clean-times=off)\n";
    std::exit(1);
  }
  C.lastLayer = std::make_unique<TTreeReaderValue<std::vector<bool>>>(reader, "Track_hitInLastLayer");
  return C;
}
// does this track carry a usable time? (clean: apply the TDR cuts)
inline bool timed(BranchPointerWrapper& b, int idx, const TimeCleaning& C, bool clean) {
  if (b.trackTimeValid[idx] != 1) return false;
  if (!clean) return true;
  return MyUtl::tdrTimeSurvivesCleaning((**C.lastLayer)[idx], b.trackEta[idx], b.trackHgtdHits[idx]);   // shared with clustering_hist
}

// What the clustering offered, for the failure breakdown (util/tdr_t0, 2026-11-09 Global Performance talk): the
// number of time clusters, the cluster time (full cluster, TRKPTZ's) closest to t_HS, and the full-cluster time of
// the cluster WAVeS chose (before its in-jet re-timing). NaN without clusters.
struct T0Extra {
  int nClusters = 0;
  double bestDt = std::numeric_limits<double>::quiet_NaN();
  double wavesFullDt = std::numeric_limits<double>::quiet_NaN();
  // the hard scatter's cluster, by its tracks' truth origin: the one holding the most hard-scatter pT. bestDt alone
  // overstates "a good cluster existed": a pileup cluster lands within 60 ps of t_HS by chance (~21% per cluster)
  double hsDt = std::numeric_limits<double>::quiet_NaN();
  double hsPurity = 0.0;   // its pT fraction from the hard scatter (0: no cluster has a hard-scatter track)
  // the hard scatter's own timed tracks in the t0 input (truth vertex 0, a cleaned time), for "are its track times
  // wrong?" and "where is its timing weight?" (the user, 2026-11-09):
  int nHsTimed = 0, nHsMistimed = 0;
  double hsPtTimed = 0.0;
  double hsPtMistimed = 0.0;    // |t - t_HS| >= 3 sigma_trk: the time is wrong
  double hsPtMisassigned = 0.0; // the TDR's definition: <= 50% of its HGTD hits from its own particle
  double hsPtInJet = 0.0;       // within dR < 0.4 of a forward jet (pT > 30, 2.38 < |eta| < 4.0): what WAVeS up-weights
};

// ── Selection levers (todo vt-5t0) ───────────────────────────────────────────
// WAVeS's score (Cluster::updateScores) generalised:
//   S(C) = [ SUM_trk pT_trk * pT_jet^alpha / max(dR, WAVES_DR_FLOOR) ] * exp(-1.5 |dz_C|) * f_vtx(C)
// with the nearest qualifying forward jet per track, and TRKPTZ's score in place of the bracket when the event has no
// qualifying forward jet -- updateScores' loop and arithmetic exactly, so alpha = 1 with f_vtx = 1 IS Score::WAVES.
enum class VtxFactor { NONE, PV0, NOT_OTHER };
struct WavesLever { int k; double alpha; VtxFactor vf; };
inline constexpr std::array<WavesLever, 5> WAVES_LEVERS = {{
    {T0_WA0, 0.0, VtxFactor::NONE}, {T0_WA025, 0.25, VtxFactor::NONE}, {T0_WA05, 0.5, VtxFactor::NONE},
    {T0_WPV0, 1.0, VtxFactor::PV0}, {T0_WNPU, 1.0, VtxFactor::NOT_OTHER}}};
// The index of the cluster the lever picks: the first maximum, as chooseCluster. When every cluster scores 0 (a vertex
// factor of 0 everywhere: none of the input reached PV0) the factor says nothing, and WAVeS's own pick (iWaves) stands.
inline int wavesLeverPick(const std::vector<Cluster>& clusters, BranchPointerWrapper& b, const WavesLever& L,
                          const std::vector<double>& zTerm, int iWaves) {
  int best = 0;
  double bestS = -1.0;
  for (int c = 0; c < (int)clusters.size(); ++c) {
    const Cluster& cl = clusters[c];
    double sum = 0.0;
    for (int idx : cl.trackIndices) {
      const double trkEta = b.trackEta[idx], trkPhi = b.trackPhi[idx], trkPt = b.trackPt[idx];
      double minDR = 1e6, nearJetPt = 0.0;
      for (int j = 0; j < (int)b.topoJetEta.GetSize(); ++j) {
        if (b.isJetRemoved(j)) continue;
        const double jEta = b.topoJetEta[j], jPt = b.topoJetPt[j];
        if (jPt < MIN_JET_PT) continue;
        if (std::abs(jEta) < MIN_ABS_ETA_JET || std::abs(jEta) > MAX_ABS_ETA_JET) continue;
        const double dr = std::hypot(jEta - trkEta, TVector2::Phi_mpi_pi(b.topoJetPhi[j] - trkPhi));
        if (dr < minDR) { minDR = dr; nearJetPt = jPt; }
      }
      if (nearJetPt <= 0.0) continue;
      sum += trkPt * (L.alpha == 1.0 ? nearJetPt : std::pow(nearJetPt, L.alpha)) / std::max(minDR, WAVES_DR_FLOOR);
    }
    double S = sum > 0.0 ? sum * zTerm[c] : cl.scores.at(Score::TRKPTZ.id);
    if (L.vf != VtxFactor::NONE) {
      double all = 0.0, pv0 = 0.0, other = 0.0;
      for (int idx : cl.trackIndices) {
        const double pt = b.trackPt[idx];
        const int v = b.recoVtxOf(idx);
        all += pt;
        if (v == 0) pv0 += pt; else if (v > 0) other += pt;   // v < 0: in no vertex fit
      }
      S *= L.vf == VtxFactor::PV0 ? pv0 / all : (all - other) / all;
    }
    if (S > bestS) { bestS = S; best = c; }
  }
  return (L.vf != VtxFactor::NONE && !(bestS > 0.0)) ? iWaves : best;
}

// The event's t0s -- util/vbs_time_veto.cxx's block, verbatim, without TZP.
inline T0Set computeT0s(BranchPointerWrapper& b, const Inflation& INFL, const TimeCleaning& C, T0Extra* x = nullptr) {
  T0Set s;
  s.t[T0_HGTD]   = b.recoVtxTime[0];
  s.sig[T0_HGTD] = b.recoVtxTimeRes[0];
  s.ok[T0_HGTD]  = (b.recoVtxValid[0] == 1);
  std::vector<int> trk_z = getAssociatedTracks(&b, MIN_TRACK_PT, MAX_TRACK_PT, T0_ASSOC_NSIGMA);
  const std::vector<int> trk_pool = trk_z;   // before cleaning: the HS track-efficiency fix's pool (applyFixes' teffPool)
  if (C.t0())   // the clustering drops untimed tracks itself; drop the timed ones cleaning un-times
    trk_z.erase(std::remove_if(trk_z.begin(), trk_z.end(),
                               [&](int i) { return b.trackTimeValid[i] == 1 && !timed(b, i, C, true); }), trk_z.end());
  const auto clusters = clusterTracksInTime(
      trk_z, &b, DIST_CUT_CONE,
      /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);
  FixCache cache;   // the fixes' re-clusterings, shared by the rows with the same input bits
  if (!clusters.empty()) {
    auto bT = chooseCluster(clusters, Score::TRKPTZ);
    s.t[T0_TRKPTZ] = bT.calculateTime(Score::TRKPTZ, &b);
    s.sig[T0_TRKPTZ] = bT.sigmas[0];  s.ok[T0_TRKPTZ] = true;
    auto bW = chooseCluster(clusters, Score::WAVES);
    if (x) { auto full = bW; x->wavesFullDt = full.calculateTime(Score::TRKPTZ, &b) - b.truthVtxTime[0]; }
    s.t[T0_WAVES] = bW.calculateTime(Score::WAVES, &b);       // in-jet refined
    s.sig[T0_WAVES] = bW.sigmas[0];  s.ok[T0_WAVES] = true;
    // the truth-fixed WAVeS t0s (the cache lets wfix_both reuse wfix_time's re-clustering)
    for (const auto& [k, mask] : {std::pair<int, int>{T0_WFIX_SEL, Fix::SELECT}, {T0_WFIX_TIME, Fix::MISTIMED},
                                  {T0_WFIX_BOTH, Fix::MISTIMED | Fix::SELECT}}) {
      const auto c = applyFixes(mask, trk_z, trk_z, clusters, &b, cache);
      if (c) { s.t[k] = c->values[0];  s.sig[k] = c->sigmas[0];  s.ok[k] = true; }
    }
    // the selection levers (todo vt-5t0): the same clusters re-scored, WAVeS's in-jet time for the pick
    int iW = 0;
    for (int c = 1; c < (int)clusters.size(); ++c)
      if (clusters[c].scores.at(Score::WAVES.id) > clusters[iW].scores.at(Score::WAVES.id)) iW = c;
    std::vector<double> zTerm(clusters.size());
    for (size_t c = 0; c < clusters.size(); ++c) {
      const float dz = std::get<1>(Cluster(clusters[c]).calcFeatures(&b));   // updateScores' rawDeltaZ
      zTerm[c] = std::exp(-1.5 * std::abs(dz));
    }
    const bool haveVtx = (bool)b.trackRecoVtxIdx;   // bound only under EXTENDED_BRANCHES, and when the input has it
    for (const WavesLever& L : WAVES_LEVERS) {
      if (L.vf != VtxFactor::NONE && !haveVtx) continue;
      const Cluster& p = clusters[wavesLeverPick(clusters, b, L, zTerm, iW)];
      s.t[L.k] = p.calculateTime(Score::WAVES, &b);   // in-jet refined, as WAVeS
      s.sig[L.k] = p.sigmas[0];  s.ok[L.k] = true;
    }
  }
  // the rest of the cumulative fix chain (clustering_hist's WFIX_C4..C6). Outside the block above: the track-efficiency
  // fix can make a t0 in an event with no timed track, as in clustering_hist
  for (const auto& [k, mask] : {std::pair<int, int>{T0_WFIX_C4, Fix::MISTIMED | Fix::SELECT | Fix::MISCL},
                                {T0_WFIX_C5, Fix::MISTIMED | Fix::SELECT | Fix::MISCL | Fix::TEFF},
                                {T0_WFIX_C6, Fix::MISTIMED | Fix::SELECT | Fix::MISCL | Fix::TEFF | Fix::RESOL}}) {
    const auto c = applyFixes(mask, trk_z, trk_pool, clusters, &b, cache);
    if (c) { s.t[k] = c->values[0];  s.sig[k] = c->sigmas[0];  s.ok[k] = true; }
  }
  // the truth time with a real t0's gate width
  s.t[T0_TRUTH_WSIG] = b.truthVtxTime[0];
  s.sig[T0_TRUTH_WSIG] = s.ok[T0_WAVES] ? s.sig[T0_WAVES] : 0.0;
  s.ok[T0_TRUTH_WSIG] = true;
  if (x) {
    for (int idx : trk_z) {
      if (b.trackToTruthvtx[idx] != 0 || b.trackTimeValid[idx] != 1) continue;   // cleaned-away times already dropped
      const double pt = b.trackPt[idx];
      ++x->nHsTimed;  x->hsPtTimed += pt;
      if (b.trackTimeRes[idx] > 0 && std::abs(b.trackTime[idx] - b.truthVtxTime[0]) / b.trackTimeRes[idx] >= 3.0) {
        ++x->nHsMistimed;  x->hsPtMistimed += pt;
      }
      if (b.trackPrimHits[idx] * 2 <= b.trackHgtdHits[idx]) x->hsPtMisassigned += pt;
      for (int j = 0; j < (int)b.topoJetPt.GetSize(); ++j) {
        if (b.isJetRemoved(j) || b.topoJetPt[j] < MIN_JET_PT) continue;
        const double ae = std::abs((double)b.topoJetEta[j]);
        if (ae < MIN_ABS_ETA_JET || ae > MAX_ABS_ETA_JET) continue;
        if (dR(b.trackEta[idx], b.trackPhi[idx], b.topoJetEta[j], b.topoJetPhi[j]) < 0.4) { x->hsPtInJet += pt; break; }
      }
    }
    x->nClusters = (int)clusters.size();
    double best = std::numeric_limits<double>::infinity();
    double hsMax = 0.0;
    for (const auto& c0 : clusters) {
      auto c = c0;
      const double d = c.calculateTime(Score::TRKPTZ, &b) - b.truthVtxTime[0];
      if (std::abs(d) < best) { best = std::abs(d); x->bestDt = d; }
      double sum = 0.0;
      for (int t : c.trackIndices) sum += b.trackPt[t];
      if (c.purity * sum > hsMax) { hsMax = c.purity * sum; x->hsDt = d; x->hsPurity = c.purity; }
    }
  }
  s.t[T0_TRUTH] = b.truthVtxTime[0];  s.sig[T0_TRUTH] = 0.0;  s.ok[T0_TRUTH] = true;
  s.infl = {INFL.hgtd, INFL.trkptz, INFL.waves, 1.0, INFL.waves, INFL.waves, INFL.waves,
            INFL.waves, INFL.waves, INFL.waves, INFL.waves, INFL.waves,
            INFL.waves, INFL.waves, INFL.waves, INFL.waves};
  return s;
}

// The WAVeS t0 with one jet left out (todo vt-9zl, 2026-11-13). Above 50 GeV a forward pileup jet whose collision
// sits close in z to the PV pulls the event's WAVeS t0 onto its own collision, so the gate keeps it (notes/gp_talk.md,
// "Why the timing gain drops above 50 GeV"). Testing a jet against a t0 it did not help make removes that:
//   LeaveOut::JET_AND_TRACKS  the jet is absent from WAVeS's jet list (score, in-jet time) AND the tracks within
//                             dR < 0.4 of it leave the t0 input, which is re-clustered: the t0 as if the jet were not there
//   LeaveOut::JET             the jet is absent from WAVeS's jet list only; its tracks stay in the clustering
// computeT0s's WAVeS block otherwise (same input, clustering, choice, in-jet time and sigma). ok = false when nothing
// is left to cluster: the gate then keeps every track, as for any event without a t0. nInJet: the chosen cluster's
// tracks within dR < 0.4 of ANOTHER qualifying forward jet (WAVeS's jet definition), i.e. whether the t0 still has
// a jet behind it -- for a guard that falls back to the event t0 when it does not.
enum class LeaveOut { JET_AND_TRACKS, JET };
struct WavesT0 { double t = 0.0, sig = 0.0; bool ok = false; int nInJet = 0; };
inline WavesT0 wavesT0LeaveOut(BranchPointerWrapper& b, int j, LeaveOut mode, const TimeCleaning& C) {
  std::vector<int> trk_z = getAssociatedTracks(&b, MIN_TRACK_PT, MAX_TRACK_PT, T0_ASSOC_NSIGMA);
  if (C.t0())
    trk_z.erase(std::remove_if(trk_z.begin(), trk_z.end(),
                               [&](int i) { return b.trackTimeValid[i] == 1 && !timed(b, i, C, true); }), trk_z.end());
  if (mode == LeaveOut::JET_AND_TRACKS) {
    const double je = b.topoJetEta[j], jp = b.topoJetPhi[j];
    trk_z.erase(std::remove_if(trk_z.begin(), trk_z.end(),
                               [&](int i) { return dR(b.trackEta[i], b.trackPhi[i], je, jp) < 0.4; }), trk_z.end());
  }
  b.leaveOutJet = j;   // WAVeS's score and in-jet time skip it; reset before returning
  WavesT0 w;
  const auto clusters = clusterTracksInTime(
      trk_z, &b, DIST_CUT_CONE,
      /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);
  if (!clusters.empty()) {
    auto bW = chooseCluster(clusters, Score::WAVES);
    w.t = bW.calculateTime(Score::WAVES, &b);
    w.sig = bW.sigmas[0];  w.ok = true;
    for (int idx : bW.trackIndices)
      for (int k = 0; k < (int)b.topoJetPt.GetSize(); ++k) {
        if (b.isJetRemoved(k) || b.topoJetPt[k] < MIN_JET_PT) continue;   // skips j too
        const double ak = std::abs((double)b.topoJetEta[k]);
        if (ak < MIN_ABS_ETA_JET || ak > MAX_ABS_ETA_JET) continue;
        if (dR(b.trackEta[idx], b.trackPhi[idx], b.topoJetEta[k], b.topoJetPhi[k]) < 0.4) { ++w.nInJet; break; }
      }
  }
  b.leaveOutJet = -1;
  return w;
}

// PV0-associated track lists: forward (rpt_v5's dz parameterisation) and central
// (z0 significance) -- the same split rpt_v5 makes, chosen per jet by its |eta|.
struct TrackLists {
  std::unordered_set<int> fwd, cen;
};
// The R_pT tracks' upper pT limit. It defaults to the time clustering's MAX_TRACK_PT (30 GeV), which rpt_v5 copied
// and this header copied from it -- never chosen for R_pT (the PUB note's macro, util/myJet_ana_fr.C, uses 45). The user,
// 2026-11-12: "why are we restricting track pT to be between 1 and 30 GeV???". --rpt-trk-ptmax=X|none moves it for
// the R_pT tracks ONLY (the t0 clustering input keeps MAX_TRACK_PT) and tags the file (rptptmaxnone / rptptmax45p0).
inline double RPT_TRK_PT_MAX = MAX_TRACK_PT;
inline void rptTrkPtMaxMode(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--rpt-trk-ptmax=", 0) != 0) continue;
    const std::string v = a.substr(16);
    RPT_TRK_PT_MAX = (v == "none") ? std::numeric_limits<double>::infinity() : std::stod(v);
    if (!(RPT_TRK_PT_MAX > MIN_TRACK_PT)) { std::cerr << "[tdr] --rpt-trk-ptmax=" << v << ": must exceed " << MIN_TRACK_PT << "\n"; std::exit(1); }
  }
  if (RPT_TRK_PT_MAX == MAX_TRACK_PT) return;   // the default: untagged
  std::string tag = "rptptmaxnone";
  if (std::isfinite(RPT_TRK_PT_MAX)) {
    char buf[32]; std::snprintf(buf, sizeof buf, "rptptmax%.1f", RPT_TRK_PT_MAX); tag = buf;
    std::replace(tag.begin(), tag.end(), '.', 'p');
  }
  SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
}
inline TrackLists rptTrackLists(BranchPointerWrapper& b) {
  TrackLists L;
  const double vtxZ = b.recoVtxZ[0];
  for (size_t trk = 0; trk < b.trackZ0.GetSize(); ++trk) {
    const double pt = b.trackPt[trk];
    if (pt < MIN_TRACK_PT || pt > RPT_TRK_PT_MAX || !b.trackQuality[trk]) continue;
    const double dz = std::abs(b.trackZ0[trk] - vtxZ);
    if (dz / MyUtl::getNewDzpara(b.trackEta[trk], pt) <= DZ0_PARA_SCALE) L.fwd.insert((int)trk);
    const double vz = b.trackVarZ0[trk];
    if (vz > 0 && dz / std::sqrt(vz) <= CENTRAL_Z_SIGNIF) L.cen.insert((int)trk);
  }
  return L;
}

enum class RptMatch { GHOST, DR02 };
inline const char* rptMatchName(RptMatch m) { return m == RptMatch::GHOST ? "ghost" : "dr02"; }
// --rpt-match=ghost|dr02 (default ghost). Call before the output path is built: dr02 tags the file.
inline RptMatch rptMatchMode(int argc, char** argv) {
  RptMatch m = RptMatch::GHOST;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--rpt-match=", 0) != 0) continue;
    const std::string v = a.substr(12);
    if (v == "ghost") m = RptMatch::GHOST;
    else if (v == "dr02") m = RptMatch::DR02;
    else { std::cerr << "[tdr] --rpt-match=" << v << ": expected ghost or dr02\n"; std::exit(1); }
  }
  if (m == RptMatch::DR02) SELECTION_TAG = SELECTION_TAG.empty() ? "dr02" : SELECTION_TAG + "_dr02";
  return m;
}

// The jet's R_pT tracks: ghost-associated, in the PV0 list, and (dr02) within the cone.
inline std::vector<int> jetRptTracks(BranchPointerWrapper& b, int j, const std::unordered_set<int>& assoc,
                                     RptMatch match) {
  std::vector<int> out;
  const double je = b.topoJetEta[j], jp = b.topoJetPhi[j];
  for (int idx : b.topoJetGhostTrackIdx[j]) {
    if (!assoc.count(idx)) continue;
    if (match == RptMatch::DR02 && dR(je, jp, b.trackEta[idx], b.trackPhi[idx]) > RPT_TRACK_JET_DR) continue;
    out.push_back(idx);
  }
  return out;
}

inline double sumPt(BranchPointerWrapper& b, const std::vector<int>& trk) {
  double s = 0.0;
  for (int i : trk) s += b.trackPt[i];
  return s;
}

// A hard-scatter track whose HGTD time is misassigned: >= TRUTH_PULL_CUT sigma from t_HS, the definition Fix::MISTIMED
// drops from the t0 input (src/failure_fixes.h). The "perfect timing" t0 sources (wfix_time, wfix_both) treat such a
// track as UNTIMED everywhere downstream too -- in the R_pT / JVF gate and in tdr_fwdpair's fixed core jet time --
// so perfect timing is one consistent world. Fixing only the t0 made the gate reject these tracks on hard-scatter
// jets that the WAVeS t0, built partly from them, had kept (2026-11-13: 99% of the HS jets whose R_pT fell).
inline bool misassignedHs(BranchPointerWrapper& b, int idx) {
  return b.trackToTruthvtx[idx] == 0 && b.trackTimeValid[idx] == 1 && b.trackTimeRes[idx] > 0 &&
         std::abs(b.trackTime[idx] - b.truthVtxTime[0]) / b.trackTimeRes[idx] >= TRUTH_PULL_CUT;
}
inline bool fixesTimes(int k) {
  return k == T0_WFIX_TIME || k == T0_WFIX_BOTH || k == T0_WFIX_C4 || k == T0_WFIX_C5 || k == T0_WFIX_C6;
}

// rpt_v5's applyTimeGate. untimeMisHs: a misassigned hard-scatter track passes as untimed (the perfect-timing sources).
inline std::vector<int> timeGate(BranchPointerWrapper& b, const std::vector<int>& base,
                                 double t0, double sig0, bool ok, double infl, double gate,
                                 const TimeCleaning& C, bool untimeMisHs = false) {
  std::vector<int> out;
  out.reserve(base.size());
  for (int idx : base) {
    if (!ok || !timed(b, idx, C, C.jets()) || (untimeMisHs && misassignedHs(b, idx))) { out.push_back(idx); continue; }
    const double dt = b.trackTime[idx] - t0;
    const double vt = b.trackTimeRes[idx] * b.trackTimeRes[idx];
    if (std::abs(dt) / std::sqrt(infl * infl * sig0 * sig0 + vt) < gate) out.push_back(idx);
  }
  return out;
}

// Self-tagging R_pT numerator: max over in-time clusters of the timed R_pT tracks of
// (cluster sum pT + untimed sum pT); the full sum when fewer than 2 tracks are timed.
inline double selfTagSumPt(BranchPointerWrapper& b, const std::vector<int>& trk, const TimeCleaning& C,
                          int* nTimedOut = nullptr) {
  std::vector<int> tt;
  double untimed = 0.0;
  for (int i : trk) {
    if (timed(b, i, C, C.jets()) && b.trackTimeRes[i] > 0.0) tt.push_back(i);
    else untimed += b.trackPt[i];
  }
  if (nTimedOut) *nTimedOut = (int)tt.size();
  if (tt.size() < 2) return sumPt(b, trk);
  static const std::unordered_map<int, double> kNone;
  std::vector<Cluster> cl = makeSimpleClusters(tt, &b, /*useSmearedTimes=*/false, kNone, kNone,
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
inline JetRpt jetRpt(BranchPointerWrapper& b, int j, const TrackLists& L, const T0Set& T, double gate,
                     const TimeCleaning& C, RptMatch match) {
  JetRpt r;
  const double pt = b.topoJetPt[j];
  const bool fwd = std::abs((double)b.topoJetEta[j]) >= FWD_ETA_MIN;
  const std::vector<int> trk = jetRptTracks(b, j, fwd ? L.fwd : L.cen, match);
  r.ntrk = (int)trk.size();
  r.itk = sumPt(b, trk) / pt;
  for (int k = 0; k < N_T0; ++k)
    r.t0[k] = sumPt(b, timeGate(b, trk, T.t[k], T.sig[k], T.ok[k], T.infl[k], gate, C, fixesTimes(k))) / pt;
  r.self = selfTagSumPt(b, trk, C, &r.ntimed) / pt;
  r.comb = T.ok[T0_HGTD] ? r.t0[T0_HGTD] : r.self;
  return r;
}

// JVT's R_pT (src/jet_tagging.h's: the pT of the jet's ghost tracks the vertex fit gave PV0, over the jet pT), plain
// and with the t0 gate on those same tracks, per t0 source: a timed PV0 track counts only within `gate` sigma of the
// t0; untimed tracks, and every track of an event without that t0, count (timeGate's rule). The time-gated JVT of
// the 2026-11-09 mini deck: at <mu> = 200 a forward pileup jet usually passes JVT on ONE PV0 track, often from another
// collision at the same z, which the gate can see. Needs EXTENDED_BRANCHES (Track_recoVtx_idx).
struct JvtRpt {
  double plain = 0.0;
  std::array<double, N_T0> t0{};
};
inline JvtRpt jvtRpt(BranchPointerWrapper& b, int j, const T0Set& T, double gate, const TimeCleaning& C) {
  std::vector<int> pv;
  for (int idx : b.topoJetGhostTrackIdx[j]) if (b.recoVtxOf(idx) == 0) pv.push_back(idx);
  JvtRpt r;
  const double pt = b.topoJetPt[j];
  r.plain = sumPt(b, pv) / pt;
  for (int k = 0; k < N_T0; ++k)
    r.t0[k] = sumPt(b, timeGate(b, pv, T.t[k], T.sig[k], T.ok[k], T.infl[k], gate, C, fixesTimes(k))) / pt;
  return r;
}

}  // namespace tdr
