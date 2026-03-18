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
    chain.Add(Form("../ntuple-hgtd/user.mcardiff.45809429.Output._%s.SuperNtuple.root", number.c_str()));
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
    const int nJets = branch->topoJetEta.GetSize();
    for (const auto& trk: tracks) {
      double
	trkEta = branch->trackEta[trk],
	trkPhi = branch->trackPhi[trk];
      bool inCone = false;
      for (int jetIdx = 0; jetIdx < nJets; ++jetIdx) {
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
  // 5b. filterTruthMatchedTracks
  //   Returns the subset of tracks whose HGTD-measured time is consistent with
  //   their true production time, thereby removing tracks with a misassigned
  //   HGTD hit.  Misassignment can affect both HS and PU tracks — any track
  //   whose reconstructed HGTD time is inconsistent with its origin vertex is
  //   removed.  Three categories are handled:
  //
  //     (1) Has truth PARTICLE link (Track_truthPart_idx != -1):
  //           Use particleT (production vertex time) as the reference.
  //           Keep if |trackTime - particleT| / trackTimeRes < pullCut.
  //     (2) No particle link but has truth VERTEX link (Track_truthVtx_idx != -1):
  //           Track is a secondary (hadronic interaction product) from that vertex.
  //           Use truthVtxTime[vtxIdx] as the reference.
  //           Keep if |trackTime - truthVtxTime| / trackTimeRes < pullCut.
  //     (3) No truth link at all: genuine fake track — always remove.
  //
  //   Tracks without a valid HGTD time (trackTimeValid != 1) are also removed
  //   since there is no measured time to check.
  //
  //   Used to build the TEST_MISAS collection: TRKPTZ clustering on a pool free
  //   of any misassigned tracks, isolating the misassignment error source.
  // ---------------------------------------------------------------------------
  std::vector<int> filterTruthMatchedTracks(
    const std::vector<int>& tracks,
    BranchPointerWrapper* branch,
    double pullCut
  ) {
    std::vector<int> output;
    output.reserve(tracks.size());
    for (int idx : tracks) {
      if (branch->trackTimeValid[idx] != 1) continue;  // no HGTD measurement

      // Determine the truth time for this track.
      // Primary tracks: use the truth particle's production vertex time (most precise).
      // Secondary tracks (hadronic interaction products, no direct particle link):
      //   fall back to the parent truth vertex time — they are real tracks from that
      //   vertex and their HGTD time should be consistent with it.
      // Completely unlinked tracks (no vertex link either) are genuine fakes: remove.
      double truthTime = NAN;
      int partIdx = branch->trackToParticle[idx];
      if (partIdx != -1) {
        truthTime = branch->particleT[partIdx];
      } else {
        int vtxIdx = branch->trackToTruthvtx[idx];
        if (vtxIdx == -1) continue;                   // completely fake → remove
        truthTime = branch->truthVtxTime[vtxIdx];
      }

      double pull = std::abs(branch->trackTime[idx] - truthTime)
                    / branch->trackTimeRes[idx];
      if (pull >= pullCut) continue;                   // HGTD time inconsistent → misassigned
      output.push_back(idx);
    }
    return output;
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
  //   Comparing TEST_MISAS with TEST_HS separates the two error sources:
  //     TRKPTZ  → TEST_MISAS : effect of fixing time misassignment (any origin)
  //     TRKPTZ  → TEST_HS    : effect of removing PU tracks from the pool
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
    a.ptr_fjet->     fillTotal(ev.effFillValFjet   );
    a.ptr_vtx_dz->   fillTotal(ev.effFillValVtxDz  );
    a.ptr_ftrack->   fillTotal(ev.effFillValTrack  );
    a.ptr_pu_frac->  fillTotal(ev.effFillValPURatio);
    a.ptr_hs_track-> fillTotal(ev.effFillValHSTrack);
    a.ptr_pu_track-> fillTotal(ev.effFillValPUTrack);
  }

  void fillPasses(AnalysisObj& a, const EventCounts& ev) {
    a.ptr_fjet->     fillPass(ev.effFillValFjet   );
    a.ptr_vtx_dz->   fillPass(ev.effFillValVtxDz  );
    a.ptr_ftrack->   fillPass(ev.effFillValTrack  );
    a.ptr_pu_frac->  fillPass(ev.effFillValPURatio);
    a.ptr_hs_track-> fillPass(ev.effFillValHSTrack);
    a.ptr_pu_track-> fillPass(ev.effFillValPUTrack);
  }

  void fillDiffs(AnalysisObj& a, const EventCounts& ev, double diff) {
    a.ptr_fjet->     fillDiff(ev.effFillValFjet,    diff);
    a.ptr_vtx_dz->   fillDiff(ev.effFillValVtxDz,  diff);
    a.ptr_ftrack->   fillDiff(ev.effFillValTrack,   diff);
    a.ptr_pu_frac->  fillDiff(ev.effFillValPURatio, diff);
    a.ptr_hs_track-> fillDiff(ev.effFillValHSTrack, diff);
    a.ptr_pu_track-> fillDiff(ev.effFillValPUTrack, diff);
  }

  void fillPurities(AnalysisObj& a, const EventCounts& ev, double purity) {
    a.ptr_fjet->     fillPurity(ev.nForwardJet,     purity);
    a.ptr_vtx_dz->   fillPurity(ev.effFillValVtxDz, purity);
    a.ptr_ftrack->   fillPurity(ev.nForwardTrack,   purity);
    a.ptr_pu_frac->  fillPurity(ev.puRatio,         purity);
    a.ptr_hs_track-> fillPurity(ev.nForwardTrackHS, purity);
    a.ptr_pu_track-> fillPurity(ev.nForwardTrackPU, purity);
  }

  // ---------------------------------------------------------------------------
  // 9. selectClusters
  //   Chooses the best cluster for every active score and returns a map of
  //   Score → Cluster.  Handles five distinct cluster collections:
  //     main clusters      — cone clustering on unfiltered tracks; covers TRKPT,
  //                          TRKPTZ, PASS, TESTML, TEST_MISCL via chooseCluster.
  //     FILTJET            — cone clustering on jet-cone-filtered tracks; selected
  //                          by TRKPTZ scoring via the single-score overload.
  //     ITERATIVE          — anti-KT-style iterative clustering on unfiltered
  //                          tracks; selected by TRKPTZ via the single-score
  //                          overload.
  //     HGTD               — simultaneous clustering on unfiltered tracks using
  //                          real HGTD times; selected by chooseHGTDCluster.
  //     HGTD_SORT          — pT-sorted simultaneous clustering; selected by the
  //                          TMVA BDT via chooseHGTDSortCluster.
  // ---------------------------------------------------------------------------
  auto selectClusters(
    const std::vector<Cluster>&        clusters,
    const std::vector<int>&            tracks,
    BranchPointerWrapper*              branch,
    const std::map<Score,AnalysisObj>& analyses,
    const std::vector<Cluster>&        filtjetClusters,
    const std::vector<Cluster>&        iterativeClusters,
    const std::vector<Cluster>&        refinedClusters,
    const std::vector<Cluster>&        testMisasClusters,
    const std::vector<Cluster>&        testHsClusters
  ) -> std::map<Score,Cluster> {
    std::map<Score,Cluster> chosen;

    // Filter: only consider clusters with enough tracks for a reliable time estimate.
    // When no qualifying cluster exists the score is absent from `chosen`, so the
    // event lands in the denominator only (null selection).
    auto qualifies = [](const Cluster& c) {
      return c.nConstituents >= MIN_CLUSTER_TRACKS;
    };
    auto filterClusters = [&](const std::vector<Cluster>& col) {
      std::vector<Cluster> out;
      for (const auto& c : col)
        if (qualifies(c)) out.push_back(c);
      return out;
    };

    // Main cone-clustering collection covers most scores
    auto qualClusters = filterClusters(clusters);
    if (!qualClusters.empty())
      chosen = chooseCluster(qualClusters, branch);

    // FILTJET: separate collection of jet-cone-filtered tracks
    if (!filtjetClusters.empty()) {
      auto qualFilt = filterClusters(filtjetClusters);
      if (!qualFilt.empty())
        chosen[Score::FILTJET] = chooseCluster(qualFilt, Score::TRKPTZ);
    }

    // ITERATIVE: nearest-neighbour iterative clustering, best selected by TRKPTZ
    if (!iterativeClusters.empty()) {
      auto qualIter = filterClusters(iterativeClusters);
      if (!qualIter.empty())
        chosen[Score::ITERATIVE] = chooseCluster(qualIter, Score::TRKPTZ);
    }

    // REFINED: iterative sub-clusters of top-k cone clusters, selected by TRKPTZ
    if (!refinedClusters.empty()) {
      auto qualRefined = filterClusters(refinedClusters);
      if (!qualRefined.empty())
        chosen[Score::REFINED] = chooseCluster(qualRefined, Score::TRKPTZ);
    }

    // TEST_MISAS: cone-clustering on non-misassigned tracks only (HS or PU).
    // Selection by TRKPTZ score; isolates the misassignment error source.
    if (!testMisasClusters.empty()) {
      auto qualMisas = filterClusters(testMisasClusters);
      if (!qualMisas.empty())
        chosen[Score::TEST_MISAS] = chooseCluster(qualMisas, Score::TRKPTZ);
    }

    // TEST_HS: cone-clustering on HS-truth-linked tracks only. Removes PU
    // contamination from the pool; HS tracks may still be misassigned.
    if (!testHsClusters.empty()) {
      auto qualHS = filterClusters(testHsClusters);
      if (!qualHS.empty())
        chosen[Score::TEST_HS] = chooseCluster(qualHS, Score::TRKPTZ);
    }

    // CONE_BDT: same cone clusters as the main collection, evaluated by the TMVA BDT.
    // chooseHGTDSortCluster hardcodes its output to scores[HGTD_SORT]; alias it to
    // CONE_BDT so that the threshold check (scores.at(CONE_BDT) > 0.3) succeeds.
    if (!qualClusters.empty() && analyses.count(Score::CONE_BDT)) {
      auto c = chooseHGTDSortCluster(qualClusters, branch);
      c.scores[Score::CONE_BDT] = c.scores[Score::HGTD_SORT];
      chosen[Score::CONE_BDT] = std::move(c);
    }

    // HGTD: simultaneous clustering using real HGTD track times.
    // Only built when HGTD is actually active (real-HGTD scenario).
    if (branch->recoVtxValid[0] == 1 && analyses.count(Score::HGTD)) {
      auto clustersHGTD = clusterTracksInTime(tracks, branch, 3.0,
                                              /*useSmearedTimes=*/false,
                                              /*checkTimeValid=*/true,
                                              /*smearRes=*/-1,
                                              ClusteringMethod::SIMULTANEOUS,
                                              /*usez0=*/false);
      auto qualHGTD = filterClusters(clustersHGTD);
      if (!qualHGTD.empty())
        chosen[Score::HGTD] = chooseHGTDCluster(qualHGTD, branch);
    }

    // HGTD_SORT: pT-sorted simultaneous clustering evaluated by TMVA BDT
    if (analyses.count(Score::HGTD_SORT)) {
      auto sortedClusters = clusterTracksInTime(tracks, branch, 3.0,
                                                /*useSmearedTimes=*/false,
                                                /*checkTimeValid=*/true,
                                                /*smearRes=*/-1,
                                                ClusteringMethod::SIMULTANEOUS,
                                                /*usez0=*/false,
                                                /*sortTracks=*/true);
      auto qualSorted = filterClusters(sortedClusters);
      if (!qualSorted.empty())
        chosen[Score::HGTD_SORT] = chooseHGTDSortCluster(qualSorted, branch);
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
  //   Return: pair<returnCode, returnVal>
  //     returnCode == 0   normal event
  //     returnCode == -1  rejected by event selection
  //     returnCode == 2   event entered TEST_MISCL denominator (purity > 0.75)
  //                       but failed the timing window (useful for event display)
  //     returnVal         TRKPTZ-selected (or TEST_MISCL) cluster time
  // ---------------------------------------------------------------------------
  std::pair<int,double> processEventData(
    BranchPointerWrapper *branch,
    bool useSmearedTimes,
    bool checkValidTimes,
    bool useZ0,
    std::map<Score,AnalysisObj>& analyses
  ) {
    // ── A. Event selection ──────────────────────────────────────────────────
    if (!branch->passBasicCuts()) return {-1, 1.};
    if (!branch->passJetPtCut())  return {-1, 1.};

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

    // ── D. Main cone clustering ─────────────────────────────────────────────
    // Purity is only computed when TEST_MISCL is active (it is expensive).
    
    const bool NEEDS_PURITY = analyses.count(Score::TEST_MISCL) > 0;
    auto clusters = clusterTracksInTime(tracks, branch, DIST_CUT_CONE,
                                        useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                                        ClusteringMethod::CONE, useZ0,
                                        /*sortTracks=*/false, NEEDS_PURITY);

    // ── E. FILTJET filtered collection ──────────────────────────────────────
    std::vector<Cluster> filtjetClusters;
    if (analyses.count(Score::FILTJET)) {
      auto filtTracks = filterTracksInJets(tracks, branch, 0.4);
      filtjetClusters = clusterTracksInTime(filtTracks, branch, DIST_CUT_CONE,
                                            useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                                            ClusteringMethod::CONE, useZ0);
    }

    // ── E-2. ITERATIVE collection ────────────────────────────────────────────
    std::vector<Cluster> iterativeClusters;
    if (analyses.count(Score::ITERATIVE)) {
      iterativeClusters = clusterTracksInTime(tracks, branch, DIST_CUT_ITER,
                                              useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                                              ClusteringMethod::ITERATIVE, useZ0);
    }

    // ── E-3. REFINED collection ──────────────────────────────────────────────
    // Two-pass timing refinement: cluster membership is determined by the 3σ
    // cone pass above; the TRKPTZ winner's reported time is then recomputed
    // using only the subset of its tracks that are within DIST_CUT_REFINE σ of
    // the centroid.  No re-clustering; only values[0]/sigmas[0] are updated.
    std::vector<Cluster> refinedClusters;
    if (analyses.count(Score::REFINED) && !clusters.empty()) {
      auto bestIt = std::max_element(clusters.begin(), clusters.end(),
          [](const Cluster& a, const Cluster& b) {
            return a.scores.at(Score::TRKPTZ) < b.scores.at(Score::TRKPTZ);
          });
      double iRes = useSmearedTimes ? IDEAL_TRACK_RES : -1.0;
      refinedClusters = { refineClusterTiming(*bestIt, branch, DIST_CUT_REFINE, iRes) };
    }

    // ── E-4. TEST_MISAS collection ───────────────────────────────────────────
    // Cone clustering on non-misassigned tracks only (any origin — HS or PU).
    // Removes tracks whose HGTD time pull vs truth exceeds TRUTH_PULL_CUT.
    // Only built when TEST_MISAS is active (real-HGTD scenario).
    std::vector<Cluster> testMisasClusters;
    if (analyses.count(Score::TEST_MISAS)) {
      auto misasTracks = filterTruthMatchedTracks(tracks, branch, TRUTH_PULL_CUT);
      if (!misasTracks.empty())
        testMisasClusters = clusterTracksInTime(misasTracks, branch, DIST_CUT_CONE,
                                                useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                                                ClusteringMethod::CONE, useZ0,false,
						NEEDS_PURITY
						);
    }

    // ── E-5. TEST_HS collection ──────────────────────────────────────────────
    // Cone clustering on HS-truth-origin tracks only. Removes all PU tracks
    // from the pool; isolates the PU-contamination / misclustering effect.
    // HS tracks may still carry misassigned HGTD times.
    // Only built when TEST_HS is active (real-HGTD scenario).
    std::vector<Cluster> testHsClusters;
    if (analyses.count(Score::TEST_HS)) {
      auto hsTracks = filterHSTracks(tracks, branch);
      if (!hsTracks.empty())
        testHsClusters = clusterTracksInTime(hsTracks, branch, DIST_CUT_CONE,
                                             useSmearedTimes, checkValidTimes, IDEAL_TRACK_RES,
                                             ClusteringMethod::CONE, useZ0);
    }

    // ── F. Fill denominator histograms ──────────────────────────────────────
    // TEST_MISCL denominator is deferred to step H where cluster purity is known.
    for (auto& [score, analysis] : analyses)
      if (!score.requiresPurity)
        fillTotals(analysis, ev);

    // ── G. Select best cluster for each score ───────────────────────────────
    auto chosen = selectClusters(clusters, tracks, branch, analyses,
                                 filtjetClusters, iterativeClusters, refinedClusters,
                                 testMisasClusters, testHsClusters);
    if (DEBUG) std::cout << "Chose clusters\n";

    // ── H. Per-score histogram filling ──────────────────────────────────────
    int    returnCode = 0;
    double returnVal  = -1.;
    bool passesMine = false, passesMiscl = false, passesMisas = false;
    bool misclInDenominator = false, misasInDenominator = false;

    for (auto& [score, analysis] : analyses) {
      if (DEBUG) std::cout << "Filling: " << toString(score) << '\n';

      // Skip scores with invalid or missing cluster data
      if (branch->recoVtxValid[0] == 0 && score == Score::HGTD)  continue;
      // if (clusters.empty() && !score.usesOwnCollection)           continue;
      if (!chosen.count(score) && score != Score::HGTD)           continue;

      Cluster& scored    = chosen[score];
      double t           = (score == Score::HGTD) ? branch->recoVtxTime[0]
                                                   : scored.values.at(0);
      double purity      = clusters.empty() ? 0.0 : scored.purity;
      double diff        = t - branch->truthVtxTime[0];

      // Expose time for the caller's event-display collection
      if (score == Score::TRKPTZ || score.requiresPurity) returnVal = t;

      // Inclusive timing resolution and purity (purity-gated for requiresPurity scores)
      if (!score.requiresPurity || purity > 0.75) {
        analysis.inclusiveReso->Fill(diff);
        analysis.inclusivePurity->Fill(purity);
        if (ev.nForwardTrackHS <= 5)
          analysis.inclusiveResoLowTrack->Fill(diff);
      }

      // Efficiency check: base timing window, plus optional score threshold.
      // For requiresPurity scores threshold encodes the purity cut, not a cluster
      // score gate, so the hasThreshold branch must be skipped for them.
      bool passes = scored.passEfficiency(branch);
      if (!score.requiresPurity && score.hasThreshold() && passes)
        passes = scored.scores.at(score) > score.threshold;

      // TEST_MISCL / TEST_MISAS: restrict both denominator and numerator to pure
      // clusters.  Purity cut comes from score.threshold when set (≥ 0), otherwise
      // falls back to 0.75.  Each score tracks its own in-denominator flag so that
      // returnCode=2 (MISCL) and returnCode=3 (MISAS) stay independent.
      if (score.requiresPurity) {
        float purityCut = score.hasThreshold() ? score.threshold : 0.75f;
        if (purity > purityCut) {
          if (score == Score::TEST_MISCL) misclInDenominator = true;
          if (score == Score::TEST_MISAS) misasInDenominator = true;
          fillTotals(analysis, ev);
        } else {
          passes = false;  // impure cluster: skip numerator entirely
        }
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
          (score == Score::TEST_MISAS && misasInDenominator)) {
        fillDiffs   (analysis, ev, diff);
        fillPurities(analysis, ev, purity);
      }
    }

    // returnCode == 2: event was in TEST_MISCL denominator but failed timing window
    if (misclInDenominator && !passesMiscl) returnCode = 2;

    // returnCode == 3: TRKPTZ passes but TEST_MISAS selected a pure cluster that
    // still fails the timing window.  Only fires when TEST_MISAS entered the
    // denominator (purity > 0.75) — impure clusters and null selections are excluded.
    // Only meaningful in the HGTD scenario where TEST_MISAS is active.
    if (analyses.count(Score::TEST_MISAS) &&
        misasInDenominator && passesMine && !passesMisas) returnCode = 3;

    return {returnCode, returnVal};
  }
}
#endif // EVENT_PROCESSING_H
