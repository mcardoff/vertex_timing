#pragma once
// ---------------------------------------------------------------------------
// failure_fixes.h -- the failure-mode fixes (MyUtl::Fix, 2026-11-11)
//
//   Each fix repairs ONE step of the WAVeS t0 with truth, on the same events,
//   so a fix row's gain over the WAVES row is what that failure mode costs.
//   The steps, in the order applyFixes runs them:
//     1. inputs     MISTIMED drops the HS tracks >= TRUTH_PULL_CUT sigma from
//                   t_HS; TEFF gives every untimed (or cleaned) HS track in
//                   HGTD a time Gaus(t_HS, FIX_TEFF_RES); RESOL puts every
//                   timed HS track at Gaus(t_HS, IDEAL_TRACK_RES). Any of the
//                   three re-clusters (iterative, DIST_CUT_CONE, as the main
//                   collection); otherwise the main collection is used as is.
//     2. selection  SELECT takes the hard scatter's cluster when a selector
//                   could realistically find it (user, 2026-11-11): >= 1 HS
//                   track and >= FIX_SELECT_MIN_PT_SHARE of the pT of all the
//                   tracks the clustering saw; of those, the most HS pT. Else
//                   WAVeS's choice. (Was: most HS pT and > 50% HS-pure.)
//     3. time       MISCL times the chosen cluster from its HS tracks alone;
//                   otherwise WAVeS's in-jet time, or with RETIME the full-
//                   cluster time when the in-jet subset is the less pure of
//                   the two (HS pT fraction).
//   Every choice is by provenance, never by the time itself: a fix that took
//   the time closer to t_HS (an earlier RETIME, the better of the in-jet and
//   full times) narrowed the WAVeS core from 17.4 to 12.1 ps on local VBF --
//   a best-of-two artefact, not a resolution.
//   Only HS truth is used, so the fixes are valid on the grid samples, whose
//   pileup tracks carry no truth link (unlike the IDEAL_RES/EFF scenarios,
//   which redraw pileup times at random there).
//
//   The draws are seeded per event from event-stable quantities, so a fix's
//   times do not depend on how TTreeProcessorMT distributes entries, and two
//   rows sharing the TEFF fill see the same filled times.
// ---------------------------------------------------------------------------

#include <TRandom3.h>
#include <TVector2.h>

#include <cmath>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"

namespace MyUtl {

  using TimeOverrideMap = std::unordered_map<int, std::pair<double,double>>;

  // A re-clustered collection and the time overrides it was built with,
  // cached per event by the input bits (MISTIMED | TEFF | RESOL), so the rows
  // that share their inputs (C2..C4 with the MISTIMED row) cluster once.
  struct FixInputs {
    std::vector<Cluster> clusters;
    TimeOverrideMap      overrides;
  };
  using FixCache = std::unordered_map<int, FixInputs>;

  inline UInt_t fixSeed(BranchPointerWrapper* b) {
    UInt_t sd = (UInt_t)std::llround(std::abs(b->truthVtxZ[0])    * 1e4)
              ^ ((UInt_t)std::llround(std::abs(b->truthVtxTime[0]) * 1e3) << 11)
              ^ ((UInt_t)b->trackZ0.GetSize() << 23);
    sd ^= 0x9E3779B9u;   // not idealisedTimingSeed's stream
    return sd ? sd : 1u;
  }

  inline const FixInputs& fixInputs(int inputBits, const std::vector<int>& tracks,
                                    const std::vector<int>& teffPool,
                                    BranchPointerWrapper* b, FixCache& cache) {
    auto it = cache.find(inputBits);
    if (it != cache.end()) return it->second;
    FixInputs out;
    const double tHS = b->truthVtxTime[0];
    std::vector<int> in;
    in.reserve(tracks.size());
    for (int i : tracks) {
      if ((inputBits & Fix::MISTIMED) && b->trackToTruthvtx[i] == 0 && b->trackTimeValid[i] == 1 &&
          std::abs(b->trackTime[i] - tHS) / b->trackTimeRes[i] >= TRUTH_PULL_CUT) continue;
      in.push_back(i);
    }
    TRandom3 rng(fixSeed(b));
    if (inputBits & Fix::TEFF) {
      // the HS tracks in HGTD with no usable time: untimed, or un-timed by the
      // cleaning (teffPool is the clustering list before it). A track MISTIMED
      // dropped is not refilled: that would correct, not drop, its time.
      std::unordered_set<int> timedIn;
      for (int i : tracks) if (b->trackTimeValid[i] == 1) timedIn.insert(i);
      std::unordered_set<int> inSet(in.begin(), in.end());
      for (int i : teffPool) {
        if (b->trackToTruthvtx[i] != 0 || timedIn.count(i)) continue;
        const double ae = std::abs((double)b->trackEta[i]);
        if (ae <= MIN_HGTD_ETA || ae >= MAX_HGTD_ETA) continue;
        out.overrides[i] = {rng.Gaus(tHS, FIX_TEFF_RES), FIX_TEFF_RES};
        if (!inSet.count(i)) { in.push_back(i); inSet.insert(i); }
      }
    }
    if (inputBits & Fix::RESOL)
      for (int i : in)
        if (b->trackToTruthvtx[i] == 0 && (out.overrides.count(i) || b->trackTimeValid[i] == 1))
          out.overrides[i] = {rng.Gaus(tHS, IDEAL_TRACK_RES), IDEAL_TRACK_RES};
    out.clusters = clusterTracksInTime(in, b, DIST_CUT_CONE, /*useSmearedTimes=*/false,
                                       /*checkTimeValid=*/true, /*smearRes=*/-1.0,
                                       ClusteringMethod::ITERATIVE, /*useZ0=*/false,
                                       /*sortTracks=*/false, /*calcPurityFlag=*/true,
                                       out.overrides.empty() ? nullptr : &out.overrides);
    return cache.emplace(inputBits, std::move(out)).first->second;
  }

  // WAVeS's in-jet time (Cluster::calculateTime, TimeSource::IN_JET with no
  // minimum subset) with the resolutions of any overridden tracks.
  // With purityOut, also the in-jet subset's HS pT fraction (the full cluster's
  // when it falls back to it).
  inline double fixInJetTime(const Cluster& c, BranchPointerWrapper* b, const TimeOverrideMap* ov,
                             double* purityOut = nullptr) {
    std::vector<std::pair<double,double>> jets;
    const int nJets = (int)b->topoJetPt.GetSize();
    for (int j = 0; j < nJets; ++j) {
      if (b->isJetRemoved(j)) continue;
      const double jEta = b->topoJetEta[j];
      if (b->topoJetPt[j] < MIN_JET_PT) continue;
      if (std::abs(jEta) < MIN_ABS_ETA_JET || std::abs(jEta) > MAX_ABS_ETA_JET) continue;
      jets.push_back({jEta, b->topoJetPhi[j]});
    }
    double sumW = 0.0, sumWT = 0.0, ptAll = 0.0, ptHS = 0.0;
    for (size_t k = 0; k < c.trackIndices.size(); ++k) {
      const int idx = c.trackIndices[k];
      bool inJet = false;
      for (auto& [je, jp] : jets)
        if (std::hypot(b->trackEta[idx] - je, TVector2::Phi_mpi_pi(b->trackPhi[idx] - jp)) < 0.4) { inJet = true; break; }
      if (!inJet) continue;
      double s = b->trackTimeRes[idx];
      if (ov) { auto it = ov->find(idx); if (it != ov->end()) s = it->second.second; }
      if (s <= 0.0) continue;
      sumW += 1.0 / (s * s);  sumWT += c.allTimes[k] / (s * s);
      ptAll += b->trackPt[idx];  if (b->trackToTruthvtx[idx] == 0) ptHS += b->trackPt[idx];
    }
    if (purityOut) *purityOut = sumW > 0.0 && ptAll > 0.0 ? ptHS / ptAll : c.purity;
    return sumW > 0.0 ? sumWT / sumW : c.values[0];
  }

  // The fixed WAVeS cluster for one row: values[0] / sigmas[0] carry its time
  // and uncertainty, purity its HS pT fraction. nullopt when no cluster
  // qualifies (the event then counts as a failure, as for the WAVES row).
  //   tracks       the clustering list (after cleaning)
  //   teffPool     the clustering list before cleaning (== tracks without it)
  //   mainClusters the main collection, purity computed
  inline std::optional<Cluster> applyFixes(int mask, const std::vector<int>& tracks,
                                           const std::vector<int>& teffPool,
                                           const std::vector<Cluster>& mainClusters,
                                           BranchPointerWrapper* b, FixCache& cache) {
    const int inputBits = mask & (Fix::MISTIMED | Fix::TEFF | Fix::RESOL);
    const std::vector<Cluster>* col = &mainClusters;
    const TimeOverrideMap* ov = nullptr;
    if (inputBits) {
      const FixInputs& fi = fixInputs(inputBits, tracks, teffPool, b, cache);
      col = &fi.clusters;
      if (!fi.overrides.empty()) ov = &fi.overrides;
    }
    std::vector<const Cluster*> qual;
    for (const auto& c : *col) if (c.nConstituents >= MIN_CLUSTER_TRACKS) qual.push_back(&c);
    if (qual.empty()) return std::nullopt;

    const Cluster* ch = nullptr;
    if (mask & Fix::SELECT) {
      double ptAll = 0.0;                    // every track the clustering saw
      for (const auto& c : *col) for (int t : c.trackIndices) ptAll += b->trackPt[t];
      double best = 0.0; const Cluster* h = nullptr;
      for (const Cluster* c : qual) {
        double pt = 0.0, ptHS = 0.0;
        for (int t : c->trackIndices) { pt += b->trackPt[t]; if (b->trackToTruthvtx[t] == 0) ptHS += b->trackPt[t]; }
        if (ptHS <= 0.0 || pt < FIX_SELECT_MIN_PT_SHARE * ptAll) continue;   // no HS track, or not selectable
        if (ptHS > best) { best = ptHS; h = c; }
      }
      ch = h;
    }
    if (!ch) {                  // WAVeS's choice (chooseCluster: the first maximum)
      ch = qual[0];
      for (const Cluster* c : qual)
        if (c->scores.at(Score::WAVES.id) > ch->scores.at(Score::WAVES.id)) ch = c;
    }

    if (mask & Fix::MISCL) {    // the chosen cluster's HS tracks alone, as one cluster
      std::vector<int> core;
      for (int t : ch->trackIndices) if (b->trackToTruthvtx[t] == 0) core.push_back(t);
      if (!core.empty()) {
        auto cc = clusterTracksInTime(core, b, /*distCut=*/1e9, false, true, -1.0,
                                      ClusteringMethod::ITERATIVE, false, false, true, ov);
        if (cc.size() == 1) return cc[0];
      }
    }
    Cluster out = *ch;
    const double tIn = fixInJetTime(out, b, ov, (mask & Fix::RETIME) ? &out.purity : nullptr);
    out.values[0] = tIn;
    if ((mask & Fix::RETIME) && out.purity < ch->purity) out.values[0] = ch->values[0];   // in-jet subset less pure
    out.purity = ch->purity;
    return out;
  }

}  // namespace MyUtl
