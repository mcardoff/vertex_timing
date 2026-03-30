#ifndef EVENT_PROCESSING_H
#define EVENT_PROCESSING_H

// ---------------------------------------------------------------------------
// event_processing.h
//   Implements the per-event analysis pipeline: ntuple chaining, track
//   pre-selection, clustering, cluster selection, and histogram filling.
//   All functions live inside the MyUtl namespace.  Sections:
//     1. setupChain (directory overload) — add all ROOT files in a directory
//     2. setupChain (single-file overload) — add one file by run number
//     3. passTrackVertexAssociation — |z₀ - vtx_z| / σ significance cut
//     4. pileupRemoval              — reject tracks near a non-HS vertex
//     5. filterTracksInJets         — keep only tracks within ΔR of a jet
//     6. getAssociatedTracks        — full HGTD-acceptance + pTV selection
//     7. EventCounts                — per-event count and histogram-fill values
//     8. Histogram fill helpers     — fillTotals / fillPasses / fillDiffs / fillPurities
//     9. selectClusters             — choose best cluster for every active score
//    10. processEventData           — main per-event analysis orchestrator
// ---------------------------------------------------------------------------

#include "clustering_includes.h"
#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "plotting_utilities.h"

using boost::filesystem::directory_iterator;

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // 1. setupChain  [directory overload]
  //   Iterates over all files in ntupleDir and adds each to the TChain.
  //   Non-file entries (sub-directories) are skipped.  Exits with an error
  //   message if the directory contains no ROOT files.  Used for the primary
  //   VBF H→Invisible sample spread across many per-run ROOT files.
  // ---------------------------------------------------------------------------
  void setupChain(
    TChain &chain, const char* ntupleDir
  ) {
    // VBF H->Invisible sample
    for (const auto& entry : directory_iterator(ntupleDir)) {
      if (entry.is_directory()) continue;
      if(DEBUG) {std::cout << "Adding file: " << entry.path() << '\n';}
      chain.Add(entry.path().c_str());
      // break;
    }

    if (chain.GetEntries() == 0) {
      std::cerr << "No ROOT files found in directory\n";
      return;
    }
  }

  // ---------------------------------------------------------------------------
  // 2. setupChain  [single-file overload]
  //   Adds a single SuperNtuple file to chain, selected by its numeric run
  //   identifier.  Convenience wrapper for quick single-file checks.
  // ---------------------------------------------------------------------------
  void setupChain(
    TChain &chain, const std::string& number
  ) {
    chain.Add(TString::Format("../ntuple-hgtd/user.mcardiff.45809429.Output._%s.SuperNtuple.root", number.c_str()));
  }

  // ---------------------------------------------------------------------------
  // 3. passTrackVertexAssociation
  //   Returns true if the z₀ distance between track trackIdx and vertex
  //   vertexIdx is within significanceCut standard deviations:
  //     |z0_trk - z_vtx| / sqrt(var_z0_trk) < significanceCut
  //   The vertex z variance is assumed to be zero (reco-vertex z is treated
  //   as exact).  Called by getAssociatedTracks and for the tighter MAX_NSIGMA
  //   filter applied after counting statistics.
  // ---------------------------------------------------------------------------
  bool passTrackVertexAssociation(
    int trackIdx, int vertexIdx,
    BranchPointerWrapper *branch,
    double significanceCut
  ) {
    double
      trkZ0    = branch->trackZ0[trackIdx],
      trkZVar = branch->trackVarZ0[trackIdx],
      vxZ      = branch->recoVtxZ[vertexIdx],
      vxZVar  = 0.0;

    double nsigmaPrim = std::abs(trkZ0 - vxZ) / std::sqrt(trkZVar + vxZVar);
    return nsigmaPrim < significanceCut;
  }

  // ---------------------------------------------------------------------------
  // 3b. Cluster::calculateTime — out-of-line implementation
  //   Defined here (after passTrackVertexAssociation) to avoid circular includes:
  //   clustering_structs.h cannot include event_processing.h.
  // ---------------------------------------------------------------------------
  inline double Cluster::calculateTime(Score score, BranchPointerWrapper* branch, double idealRes) const {
    if (score == Score::HGTD)
      return branch->recoVtxTime[0];

    // Z_REFINED: precision-weighted mean of tracks passing z₀ TVA at TVA_CUT_Z_REFINED σ
    if (score == Score::Z_REFINED) {
      double sumW = 0.0, sumWT = 0.0;
      for (size_t i = 0; i < trackIndices.size(); ++i) {
        int idx = trackIndices[i];
        if (!passTrackVertexAssociation(idx, 0, branch, TVA_CUT_Z_REFINED)) continue;
        double t = allTimes[i];
        double s = (idealRes > 0.0) ? idealRes : (double)branch->trackTimeRes[idx];
        if (s <= 0.0) continue;
        double w = 1.0 / (s * s);
        sumW += w; sumWT += w * t;
      }
      if (sumW > 0.0) return sumWT / sumW;
      return values[0];  // fallback: no tracks survived z₀ TVA
    }

    // ZT_REFINED uses its own re-clustered collection (z₀-filtered tracks re-clustered
    // at DIST_CUT_T_REFINED); values[0] is already set correctly — fall through.

    return values[0];
  }

  // ---------------------------------------------------------------------------
  // 4. pileupRemoval
  //   Removes tracks that are significantly associated with a non-HS reco
  //   vertex.  For each track the nearest reco vertex (stored in
  //   Track_nearestVtx_idx) is identified; if that vertex is not vertex 0
  //   (the HS vertex) and the track-to-nearest-vertex significance is below
  //   significanceCut, the track is dropped.  Returns the surviving subset.
  //   Note: z₀sinθ and its uncertainty are corrected back to plain z₀ using
  //   the track's θ and θ variance.
  // ---------------------------------------------------------------------------
  std::vector<int> pileupRemoval(
    const std::vector<int>& tracks,
    BranchPointerWrapper *branch,
    double significanceCut
  ) {
    std::vector<int> output;
    for (const auto& trk: tracks) {
      int near = (int)branch->trackNearIdx[trk];
      double nearZ0 = branch->trackNearZ0sin[trk] / std::sin(branch->trackTheta[trk]);
      double nearVarZ0 = (pow(branch->trackNearZ0sinUnc[trk],2)-pow(nearZ0*std::cos(branch->trackTheta[trk])*branch->trackVarTheta[trk],2))/std::sin(branch->trackTheta[trk])*std::sin(branch->trackTheta[trk]);

      double sigNear = std::abs(nearZ0 - branch->recoVtxZ[near]) / std::sqrt(nearVarZ0);

      if (near != 0 && sigNear < significanceCut)
	continue;

      output.push_back(trk);
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 5. filterTracksInJets
  //   Retains only tracks that lie within minDRCut of at least one reco jet
  //   with pT > MIN_JET_PT.  ΔR is computed using the standard
  //   sqrt(Δη² + Δφ²) metric.  Used to build the FILTJET collection, which
  //   restricts clustering to tracks geometrically associated with jets.
  // ---------------------------------------------------------------------------
  std::vector<int> filterTracksInJets(
    const std::vector<int>& tracks,
    BranchPointerWrapper *branch,
    double minDRCut
  ) {
    std::vector<int> output;
    output.reserve(tracks.size()); // upper bound; shrink_to_fit not worth the cost
    const int N_JETS = branch->topoJetEta.GetSize();
    for (const auto& trk: tracks) {
      double
	trkEta = branch->trackEta[trk],
	trkPhi = branch->trackPhi[trk];
      bool inCone = false;
      for (int jetIdx = 0; jetIdx < N_JETS; ++jetIdx) {
	if (branch->topoJetPt[jetIdx] < MIN_JET_PT) continue;
	double
	  deta = branch->topoJetEta[jetIdx] - trkEta,
	  dphi = TVector2::Phi_mpi_pi(branch->topoJetPhi[jetIdx] - trkPhi);
	if (std::hypot(deta, dphi) < minDRCut) {
	  inCone = true;
	  break; // found a qualifying jet — no need to check the rest
	}
      }
      if (inCone)
	output.push_back(trk);
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 6. getAssociatedTracks
  //   Scans the full track array and returns indices of tracks that pass:
  //     • HGTD η acceptance: MIN_HGTD_ETA < |η| < MAX_HGTD_ETA
  //     • pT window: minTrkPt < pT < maxTrkPt
  //     • Track quality flag
  //     • passTrackVertexAssociation at significance_cut
  //   In processEventData this is first called at 3σ to collect statistics,
  //   then the list is filtered down to MAX_NSIGMA for the clustering step.
  // ---------------------------------------------------------------------------
  std::vector<int> getAssociatedTracks(
      BranchPointerWrapper *branch,
      double minTrkPt, double maxTrkPt,
      double significanceCut
    ) {
    std::vector<int> goodTracks;

    for (size_t trk = 0; trk < branch->trackZ0.GetSize(); ++trk) {
      double
	trkEta = branch->trackEta[trk],
	trkPt  = branch->trackPt[trk],
	trkQuality = branch->trackQuality[trk];

      if (std::abs(trkEta) < MIN_HGTD_ETA or
	  std::abs(trkEta) > MAX_HGTD_ETA)
        continue;

      if (trkPt < minTrkPt or trkPt > maxTrkPt)
	continue;

      if (not trkQuality)
	continue;

      if (passTrackVertexAssociation(trk, 0, branch, significanceCut))
	goodTracks.push_back(trk);
    }

    return goodTracks;
  }

  // ---------------------------------------------------------------------------
  // 5b. calcHSTimingPurity
  //   Computes the fraction of forward HS track pT whose HGTD time is consistent
  //   with the truth vertex (|pull| < HS_TIMING_QUALITY_CUT).
  //   Used as an event-level gate for Score::TEST_MISAS: only events where
  //   >= 75% of HS pT has a good timing measurement enter its denominator.
  //   Returns 1.0 (vacuously clean) when no HS tracks with valid HGTD time exist.
  //
  //   Truth time resolution:
  //     Primary tracks   (trackToParticle != -1): use particleT (most precise).
  //     Secondary tracks (no particle link, vtxIdx != -1): use truthVtxTime.
  //     Unlinked fakes   (vtxIdx == -1): excluded from denominator.
  // ---------------------------------------------------------------------------
  float calcHSTimingPurity(
    const std::vector<int>& tracks,
    BranchPointerWrapper* branch
  ) {
    double num = 0.0, denom = 0.0;
    for (int idx : tracks) {
      if (branch->trackToTruthvtx[idx] != 0) continue;  // not HS
      if (branch->trackTimeValid[idx]  != 1) continue;  // no HGTD measurement

      int    partIdx   = branch->trackToParticle[idx];
      double truthTime = (partIdx != -1)
                         ? branch->particleT[partIdx]
                         : branch->truthVtxTime[branch->trackToTruthvtx[idx]];

      double pull = std::abs(branch->trackTime[idx] - truthTime)
                    / branch->trackTimeRes[idx];
      double pT   = branch->trackPt[idx];
      denom += pT;
      if (pull < HS_TIMING_QUALITY_CUT) num += pT;
    }
    return (denom > 0.0) ? static_cast<float>(num / denom) : 0.0f;
  }

  // ---------------------------------------------------------------------------
  // 5b2. calcClusterHSTimingPurity
  //   Cluster-scoped variant of calcHSTimingPurity.  Iterates only over the
  //   tracks in `cluster.trackIndices`, returning the fraction of HS pT that
  //   passes the |pull| < HS_TIMING_QUALITY_CUT check.  Returns 0.0 if no
  //   qualifying HS tracks are found (pure PU cluster — excluded from CTIME
  //   denominator since there is no HS timing to evaluate).
  // ---------------------------------------------------------------------------
  float calcClusterHSTimingPurity(
    const Cluster& cluster,
    BranchPointerWrapper* branch
  ) {
    double num = 0.0, denom = 0.0;
    for (int idx : cluster.trackIndices) {
      if (branch->trackToTruthvtx[idx] != 0) continue;
      if (branch->trackTimeValid[idx]  != 1) continue;

      int    partIdx   = branch->trackToParticle[idx];
      double truthTime = (partIdx != -1)
                         ? branch->particleT[partIdx]
                         : branch->truthVtxTime[branch->trackToTruthvtx[idx]];

      double pull = std::abs(branch->trackTime[idx] - truthTime)
                    / branch->trackTimeRes[idx];
      double pT   = branch->trackPt[idx];
      denom += pT;
      if (pull < HS_TIMING_QUALITY_CUT) num += pT;
    }
    return (denom > 0.0) ? static_cast<float>(num / denom) : 0.0f;
  }

  // ---------------------------------------------------------------------------
  // 5c. filterHSTracks
  //   Returns the subset of tracks that truth-link to the hard-scatter vertex
  //   (Track_truthVtx_idx == 0) and have a valid HGTD time measurement.
  //   Removes all pileup-origin and unlinked (fake) tracks from the pool.
  //
  //   Note: HS tracks can still carry a misassigned HGTD hit; this filter does
  //   NOT apply a pull check.  TEST_HS therefore isolates the effect of PU
  //   contamination and misclustering, leaving HS-track misassignment in place.
  // ---------------------------------------------------------------------------
  std::vector<int> filterHSTracks(
    const std::vector<int>& tracks,
    BranchPointerWrapper* branch
  ) {
    std::vector<int> output;
    output.reserve(tracks.size());
    for (int idx : tracks) {
      if (branch->trackTimeValid[idx] != 1) continue;   // no HGTD measurement
      if (branch->trackToTruthvtx[idx]  != 0) continue; // not from HS vertex
      output.push_back(idx);
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 7. EventCounts
  //   Holds all per-event counts and pre-folded histogram fill values
  //   derived from a 3σ track selection.  Computed once at the top of
  //   processEventData and passed to the histogram fill helpers below.
  // ---------------------------------------------------------------------------
  struct EventCounts {
    // Raw forward track/jet counts
    int nForwardJet      = 0;
    int nForwardTrack    = 0;
    int nForwardTrackHS  = 0;
    int nForwardTrackPU  = 0;
    double puRatio       = 0.0;
    // Pre-folded x-values for efficiency / purity histograms
    int    effFillValFjet;
    int    effFillValTrack;
    int    effFillValHSTrack;
    int    effFillValPUTrack;
    double effFillValPURatio;
    double effFillValVtxDz;

    EventCounts(BranchPointerWrapper* branch,
                const std::vector<int>& tracks,
                bool checkValidTimes) {
      branch->countForwardJets(nForwardJet);
      branch->countForwardTracks(nForwardTrack, nForwardTrackHS, nForwardTrackPU,
                                 tracks, checkValidTimes);
      puRatio          = (double)nForwardTrackPU / (double)nForwardTrack;
      effFillValFjet    = folded(nForwardJet,     (int)FOLD_FJET);
      effFillValTrack   = folded(nForwardTrack,   (int)FOLD_TRACK);
      effFillValHSTrack = folded(nForwardTrackHS, (int)FOLD_HS_TRACK);
      effFillValPUTrack = folded(nForwardTrackPU, (int)FOLD_PU_TRACK);
      effFillValPURatio = folded(puRatio,         FOLD_PU_FRAC);
      effFillValVtxDz   = std::abs(branch->recoVtxZ[0] - branch->truthVtxZ[0]);
    }
  };

  // ---------------------------------------------------------------------------
  // 8. Histogram fill helpers
  //   Each function fills one histogram operation (total/pass/diff/purity)
  //   across all six analysis sub-categories (fjet, vtx_dz, ftrack, pu_frac,
  //   hs_track, pu_track) in a single call, eliminating the repeated pattern
  //   of six individual fill calls at each call site.
  // ---------------------------------------------------------------------------

  void fillTotals(AnalysisObj& a, const EventCounts& ev) {
    a.ptrFjet->    fillTotal(ev.effFillValFjet   );
    a.ptrVtxDz->   fillTotal(ev.effFillValVtxDz  );
    a.ptrFtrack->  fillTotal(ev.effFillValTrack  );
    a.ptrPuFrac->  fillTotal(ev.effFillValPURatio);
    a.ptrHSTrack-> fillTotal(ev.effFillValHSTrack);
    a.ptrPUTrack-> fillTotal(ev.effFillValPUTrack);
  }

  void fillPasses(AnalysisObj& a, const EventCounts& ev) {
    a.ptrFjet->    fillPass(ev.effFillValFjet   );
    a.ptrVtxDz->   fillPass(ev.effFillValVtxDz  );
    a.ptrFtrack->  fillPass(ev.effFillValTrack  );
    a.ptrPuFrac->  fillPass(ev.effFillValPURatio);
    a.ptrHSTrack-> fillPass(ev.effFillValHSTrack);
    a.ptrPUTrack-> fillPass(ev.effFillValPUTrack);
  }

  void fillDiffs(AnalysisObj& a, const EventCounts& ev, double diff) {
    a.ptrFjet->    fillDiff(ev.effFillValFjet,    diff);
    a.ptrVtxDz->   fillDiff(ev.effFillValVtxDz,  diff);
    a.ptrFtrack->  fillDiff(ev.effFillValTrack,   diff);
    a.ptrPuFrac->  fillDiff(ev.effFillValPURatio, diff);
    a.ptrHSTrack-> fillDiff(ev.effFillValHSTrack, diff);
    a.ptrPUTrack-> fillDiff(ev.effFillValPUTrack, diff);
  }

  void fillPurities(AnalysisObj& a, const EventCounts& ev, double purity) {
    a.ptrFjet->    fillPurity(ev.nForwardJet,     purity);
    a.ptrVtxDz->   fillPurity(ev.effFillValVtxDz, purity);
    a.ptrFtrack->  fillPurity(ev.nForwardTrack,   purity);
    a.ptrPuFrac->  fillPurity(ev.puRatio,         purity);
    a.ptrHSTrack-> fillPurity(ev.nForwardTrackHS, purity);
    a.ptrPUTrack-> fillPurity(ev.nForwardTrackPU, purity);
  }

  // ---------------------------------------------------------------------------
  // 9. selectClusters
  //   Chooses the best cluster for every active score and returns a map of
  //   Score.id → Cluster.  Two collection types:
  //     main clusters   — iterative clustering on all tracks; all scores with
  //                       usesOwnCollection=false are chosen here via the
  //                       all-scores chooseCluster overload.
  //     auxCollections  — map built in section E, one entry per dedicated-
  //                       collection score; each selected by TRKPTZ.
  //   Scores handled separately (non-TRKPTZ selection):
  //     CONE_BDT        — main clusters, TMVA BDT selector
  //     HGTD            — simultaneous clustering on real HGTD times
  //     HGTD_SORT       — pT-sorted simultaneous clustering, TMVA BDT
  // ---------------------------------------------------------------------------
  auto selectClusters(
    const std::vector<Cluster>&                         clusters,
    const std::vector<int>&                             tracks,
    BranchPointerWrapper*                               branch,
    const std::map<Score,AnalysisObj>&                  analyses,
    const std::unordered_map<int,std::vector<Cluster>>& auxCollections
  ) -> std::unordered_map<int,Cluster> {
    std::unordered_map<int,Cluster> chosen;

    auto filterClusters = [](const std::vector<Cluster>& col) {
      std::vector<Cluster> out;
      for (const auto& c : col)
        if (c.nConstituents >= MIN_CLUSTER_TRACKS) out.push_back(c);
      return out;
    };

    // Main collection: all scores with usesOwnCollection=false
    auto qualMain = filterClusters(clusters);
    if (!qualMain.empty())
      chosen = chooseCluster(qualMain, branch);

    // Dedicated collections: each selected by TRKPTZ
    for (const auto& [id, col] : auxCollections) {
      auto qual = filterClusters(col);
      if (!qual.empty())
        chosen[id] = chooseCluster(qual, Score::TRKPTZ);
    }

    // CONE_BDT: main clusters, TMVA BDT selector (aliases HGTD_SORT score)
    if (!qualMain.empty() && analyses.count(Score::CONE_BDT)) {
      auto c = chooseHGTDSortCluster(qualMain, branch);
      c.scores[Score::CONE_BDT.id] = c.scores[Score::HGTD_SORT.id];
      chosen[Score::CONE_BDT.id] = std::move(c);
    }

    // HGTD: simultaneous clustering on real HGTD times
    if (branch->recoVtxValid[0] == 1 && analyses.count(Score::HGTD)) {
      auto col = clusterTracksInTime(tracks, branch, 3.0,
                                     false, true, -1,
                                     ClusteringMethod::SIMULTANEOUS, false);
      auto qual = filterClusters(col);
      if (!qual.empty())
        chosen[Score::HGTD.id] = chooseHGTDCluster(qual, branch);
    }

    // HGTD_SORT: pT-sorted simultaneous clustering, TMVA BDT selector
    if (analyses.count(Score::HGTD_SORT)) {
      auto col = clusterTracksInTime(tracks, branch, 3.0,
                                     false, true, -1,
                                     ClusteringMethod::SIMULTANEOUS, false, true);
      auto qual = filterClusters(col);
      if (!qual.empty())
        chosen[Score::HGTD_SORT.id] = chooseHGTDSortCluster(qual, branch);
    }

    return chosen;
  }

  // ---------------------------------------------------------------------------
  // 10. processEventData
  //   Main per-event analysis orchestrator.  Called once per event per timing
  //   scenario (HGTD / IdealRes / IdealEff).  Pipeline:
  //     A) Event selection — passBasicCuts, passJetPtCut.
  //     B) Track selection — 3σ scan for counting statistics; optionally
  //        tighten to MAX_NSIGMA for the clustering step.
  //     C) Per-event counts — forward jet / track counts and folded fill values.
  //     D) Main cone clustering — covers TRKPT, TRKPTZ, PASS, TESTML, TEST_MISCL.
  //     E) FILTJET collection — jet-cone-filtered tracks (if score active).
  //     F) Denominator fills — fillTotal for all scores except TEST_MISCL
  //        (deferred until the selected cluster's purity is known).
  //     G) Cluster selection — selectClusters() builds the Score→Cluster map.
  //     H) Per-score fills — inclusive reso, efficiency pass/total, diff, purity.
  //
  // ---------------------------------------------------------------------------
  // EventResult — returned by processEventData
  //   code        — -1: rejected by selection; 0: normal; 2: MISCL fail;
  //                 3: MISAS fail
  //   time        — TRKPTZ-selected (or TEST_MISCL) cluster time
  //   nFwdHS      — n forward HS tracks (3σ counting step)
  //   trkptzPass  — true if TRKPTZ passed the PASS_SIGMA timing window
  // ---------------------------------------------------------------------------
  struct EventResult {
    int    code       = -1;
    double time       = -1.0;
    int    nFwdHS     =  0;
    bool   trkptzPass = false;
  };

  EventResult processEventData(
    BranchPointerWrapper *branch,
    bool useSmearedTimes,
    bool checkValidTimes,
    std::map<Score,AnalysisObj>& analyses
  ) {
    // ── A. Event selection ──────────────────────────────────────────────────
    if (!branch->passBasicCuts()) return {};
    if (!branch->passJetPtCut())  return {};

    // ── B. Track selection ──────────────────────────────────────────────────
    // Scan once at 3σ so counting statistics include slightly-displaced tracks,
    // then optionally tighten to MAX_NSIGMA for the clustering step.
    std::vector<int> tracks = getAssociatedTracks(branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);

    // ── C. Per-event counts ─────────────────────────────────────────────────
    EventCounts ev(branch, tracks, checkValidTimes);

    if (DEBUG) {
      std::cout << "nForwardJet = "    << ev.nForwardJet     << '\n';
      std::cout << "nForwardTrack = "  << ev.nForwardTrack   << '\n';
      std::cout << "nForwardTrack_HS = " << ev.nForwardTrackHS << '\n';
      std::cout << "nForwardTrack_PU = " << ev.nForwardTrackPU << '\n';
    }

    if (MAX_NSIGMA != 3.0)
      tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
        [&](int trk) {
          return !passTrackVertexAssociation(trk, 0, branch, MAX_NSIGMA);
        }), tracks.end());

    // ── D. Main iterative clustering ─────────────────────────────────────────
    // Purity is only computed when TEST_MISCL is active (it is expensive).

    const bool NEEDS_PURITY = analyses.count(Score::TEST_MISCL) > 0;
    auto clusters = clusterTracksInTime(tracks, branch, DIST_CUT_CONE,
                                        useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                                        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
                                        /*sortTracks=*/false, NEEDS_PURITY);

    // ── E. Dedicated per-score cluster collections ──────────────────────────
    // Each Score with buildsCollection()==true carries its own clustering spec
    // (distCut, method, useZ0, filter).  A single loop over SCORE_REGISTRY
    // builds all active dedicated collections — no separate table needed.
    auto applyFilter = [&](TrackFilterType ft) -> std::vector<int> {
      switch (ft) {
        case TrackFilterType::JET:    return filterTracksInJets(tracks, branch, 0.4);
        case TrackFilterType::Z0_TVA: {
          std::vector<int> out;
          for (int idx : tracks)
            if (passTrackVertexAssociation(idx, 0, branch, TVA_CUT_Z_REFINED))
              out.push_back(idx);
          return out;
        }
        case TrackFilterType::HS_ONLY: return filterHSTracks(tracks, branch);
        default:                       return tracks;
      }
    };

    std::unordered_map<int, std::vector<Cluster>> auxCollections;
    for (const auto& s : SCORE_REGISTRY) {
      if (!s.buildsCollection() || !analyses.count(s) || tracks.empty()) continue;
      auto t = applyFilter(s.filter);
      if (t.empty()) continue;
      auxCollections[s.id] = clusterTracksInTime(t, branch, s.distCut,
                               useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                               s.method, s.useZ0, /*sortTracks=*/false, NEEDS_PURITY);
    }

    // ── F. Fill denominator histograms ──────────────────────────────────────
    // TEST_MISCL denominator is deferred to step H where cluster purity is known.
    for (auto& [score, analysis] : analyses)
      if (!score.requiresPurity)
        fillTotals(analysis, ev);

    // ── G. Select best cluster for each score ───────────────────────────────
    auto chosen = selectClusters(clusters, tracks, branch, analyses, auxCollections);
    if (DEBUG) std::cout << "Chose clusters\n";

    // ── G½. HS timing purity for TEST_MISAS gate ────────────────────────────
    float hsTimingPurity = 1.0f;
    if (analyses.count(Score::TEST_MISAS) || analyses.count(Score::TEST_CTIME))
      hsTimingPurity = calcHSTimingPurity(tracks, branch);

    // ── H. Per-score histogram filling ──────────────────────────────────────
    int    returnCode = 0;
    double returnVal  = -1.;
    bool passesMine = false, passesMiscl = false, passesMisas = false;
    bool misclInDenominator = false, misasInDenominator = false, ctimeInDenominator = false;
    bool perfEvtInDenominator = false, perfCltInDenominator = false;

    for (auto& [score, analysis] : analyses) {
      if (DEBUG) std::cout << "Filling: " << toString(score) << '\n';

      // Skip scores with invalid or missing cluster data
      if (branch->recoVtxValid[0] == 0 && score == Score::HGTD)   continue;
      if (clusters.empty() && !score.usesOwnCollection)           continue;
      if (!chosen.count(score.id) && score != Score::HGTD)         continue;

      Cluster& scored    = chosen.at(score.id);
      double iResH       = useSmearedTimes ? IDEAL_TRACK_RES : -1.0;
      double t           = scored.calculateTime(score, branch, iResH);
      scored.values[0]   = t;   // keep passEfficiency consistent with diff
      double purity      = clusters.empty() ? 0.0 : scored.purity;
      double diff        = t - branch->truthVtxTime[0];

      // Expose time for the caller's event-display collection
      if (score == Score::TRKPTZ || score.requiresPurity) returnVal = t;

      auto fillResoStack = [&](TH1D* sig, TH1D* mix, TH1D* bkg) {
        if      (purity > 0.75)  sig->Fill(diff);
        else if (purity >= 0.50) mix->Fill(diff);
        else                     bkg->Fill(diff);
      };

      // Efficiency check: base timing window, plus optional score threshold.
      // For requiresPurity scores threshold encodes the purity cut, not a cluster
      // score gate, so the hasThreshold branch must be skipped for them.
      bool passes = scored.passEfficiency(branch);
      if (!score.requiresPurity && score.hasThreshold() && passes)
        passes = scored.scores.at(score.id) > score.threshold;

      // TEST_MISCL / TEST_MISAS: restrict both denominator and numerator to pure
      // clusters.  Purity cut comes from score.threshold when set (≥ 0), otherwise
      // falls back to 1.00.  Each score tracks its own in-denominator flag so that
      // returnCode=2 (MISCL) and returnCode=3 (MISAS) stay independent.
      bool inDenominator = !score.requiresPurity;
      if (score.requiresPurity) {
        bool passesGate;
        if (score == Score::TEST_MISAS) {
          // Gate on event-level HS timing purity: 100% of HS pT must have |pull|<3σ
          passesGate = (hsTimingPurity >= 1.0f);
        } else if (score == Score::TEST_CTIME) {
          // Gate on cluster-level HS timing purity: 100% of HS pT must have |pull|<3σ
          float cp = calcClusterHSTimingPurity(scored, branch);
          passesGate = (cp >= 1.0f);
        } else if (score == Score::PERF_EVT) {
          // MISCL ∧ MISAS: pure cluster AND all event-level HS tracks correctly timed
          passesGate = (purity > 0.75f) && (hsTimingPurity >= 1.0f);
        } else if (score == Score::PERF_CLT) {
          // CTIME ∧ MISCL: cluster HS tracks correctly timed AND cluster is pure (>75% HS)
          // Note: CTIME ∧ MISAS is trivially equivalent to MISAS alone because MISAS
          // (event-level 100% timing) strictly implies CTIME (cluster is a subset).
          float cp = calcClusterHSTimingPurity(scored, branch);
          passesGate = (cp >= 1.0f) && (purity > 0.75f);
        } else {
          float purityCut = score.hasThreshold() ? score.threshold : 0.75f;
          passesGate = (purity > purityCut);
        }
        inDenominator = passesGate;
        if (passesGate) {
          if (score == Score::TEST_MISCL) misclInDenominator = true;
          if (score == Score::TEST_MISAS) misasInDenominator = true;
          if (score == Score::TEST_CTIME) ctimeInDenominator = true;
          if (score == Score::PERF_EVT)   perfEvtInDenominator = true;
          if (score == Score::PERF_CLT)   perfCltInDenominator = true;
          fillTotals(analysis, ev);
        } else {
          passes = false;
        }
      }

      // Inclusive resolution split by cluster purity; only fill if in denominator
      analysis.inclusivePurity->Fill(purity);
      if (inDenominator) {
        fillResoStack(analysis.inclusiveResoSig.get(),
                      analysis.inclusiveResoMix.get(),
                      analysis.inclusiveResoBkg.get());
        if (ev.nForwardTrackHS <= 5)
          fillResoStack(analysis.inclusiveResoLowTrackSig.get(),
                        analysis.inclusiveResoLowTrackMix.get(),
                        analysis.inclusiveResoLowTrackBkg.get());
      }

      if (passes) {
        fillPasses(analysis, ev);
        if (score == Score::TRKPTZ)     passesMine  = true;
        if (score == Score::TEST_MISCL) passesMiscl = true;
        if (score == Score::TEST_MISAS) passesMisas = true;
      }

      // 2D timing residual and purity distributions
      // (requiresPurity scores: only fill events that entered the denominator)
      if (!score.requiresPurity ||
          (score == Score::TEST_MISCL && misclInDenominator) ||
          (score == Score::TEST_MISAS && misasInDenominator) ||
          (score == Score::TEST_CTIME && ctimeInDenominator) ||
          (score == Score::PERF_EVT  && perfEvtInDenominator) ||
          (score == Score::PERF_CLT  && perfCltInDenominator)) {
        fillDiffs   (analysis, ev, diff);
        fillPurities(analysis, ev, purity);
      }
    }

    // returnCode == 2: event was in TEST_MISCL denominator but failed timing window
    if (misclInDenominator && !passesMiscl) returnCode = 2;

    // returnCode == 3: event has clean HS timing (hsTimingPurity > 0.75, in TEST_MISAS
    // denominator) but TRKPTZ still fails — the failure is not caused by timing
    // misassignment.  Only meaningful in the HGTD scenario where TEST_MISAS is active.
    if (analyses.count(Score::TEST_MISAS) &&
        misasInDenominator && !passesMine) returnCode = 3;

    return {returnCode, returnVal, ev.nForwardTrackHS, passesMine};
  }
}
#endif // EVENT_PROCESSING_H
