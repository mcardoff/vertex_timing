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
  //   with pT > MIN_JETPT.  ΔR is computed using the standard
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
	if (branch->topoJetPt[jetIdx] < MIN_JETPT) continue;
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
    a.ptr_fjet->     fillDiff(ev.nForwardJet,     diff);
    a.ptr_vtx_dz->   fillDiff(ev.effFillValVtxDz, diff);
    a.ptr_ftrack->   fillDiff(ev.nForwardTrack,   diff);
    a.ptr_pu_frac->  fillDiff(ev.puRatio,         diff);
    a.ptr_hs_track-> fillDiff(ev.nForwardTrackHS, diff);
    a.ptr_pu_track-> fillDiff(ev.nForwardTrackPU, diff);
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
  //   Score → Cluster.  Handles four distinct cluster collections:
  //     main clusters  — cone clustering on unfiltered tracks; covers TRKPT,
  //                      TRKPTZ, PASS, TESTML, TEST_MISCL via chooseCluster.
  //     FILTJET        — cone clustering on jet-cone-filtered tracks; selected
  //                      by TRKPTZ scoring via the single-score overload.
  //     HGTD           — simultaneous clustering on unfiltered tracks using
  //                      real HGTD times; selected by chooseHGTDCluster.
  //     HGTD_SORT      — pT-sorted simultaneous clustering; selected by the
  //                      TMVA BDT via chooseHGTDSortCluster.
  // ---------------------------------------------------------------------------
  auto selectClusters(
    const std::vector<Cluster>&        clusters,
    const std::vector<int>&            tracks,
    BranchPointerWrapper*              branch,
    const std::map<Score,AnalysisObj>& analyses,
    const std::vector<Cluster>&        filtjetClusters
  ) -> std::map<Score,Cluster> {
    std::map<Score,Cluster> chosen;

    // Main cone-clustering collection covers most scores
    if (!clusters.empty())
      chosen = chooseCluster(clusters, branch);

    // FILTJET: separate collection of jet-cone-filtered tracks
    if (!filtjetClusters.empty())
      chosen[Score::FILTJET] = chooseCluster(filtjetClusters, Score::TRKPTZ);

    // HGTD: simultaneous clustering using real HGTD track times.
    // Only built when HGTD is actually active (real-HGTD scenario).
    if (branch->recoVtxValid[0] == 1 && analyses.count(Score::HGTD)) {
      auto clustersHGTD = clusterTracksInTime(tracks, branch, 3.0,
                                              /*useSmearedTimes=*/false,
                                              /*checkTimeValid=*/true,
                                              /*smearRes=*/-1,
                                              /*useCone=*/false,
                                              /*usez0=*/false);
      chosen[Score::HGTD] = chooseHGTDCluster(clustersHGTD, branch);
    }

    // HGTD_SORT: pT-sorted simultaneous clustering evaluated by TMVA BDT
    if (analyses.count(Score::HGTD_SORT)) {
      auto sortedClusters = clusterTracksInTime(tracks, branch, 3.0,
                                                /*useSmearedTimes=*/false,
                                                /*checkTimeValid=*/true,
                                                /*smearRes=*/-1,
                                                /*useCone=*/false,
                                                /*usez0=*/false,
                                                /*sortTracks=*/true);
      if (!sortedClusters.empty())
        chosen[Score::HGTD_SORT] = chooseHGTDSortCluster(sortedClusters, branch);
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
    auto clusters = clusterTracksInTime(tracks, branch, 3.0,
                                        useSmearedTimes, checkValidTimes, 10.0,
                                        /*useCone=*/true, useZ0,
                                        /*sortTracks=*/false, NEEDS_PURITY);

    // ── E. FILTJET filtered collection ──────────────────────────────────────
    std::vector<Cluster> filtjetClusters;
    if (analyses.count(Score::FILTJET)) {
      auto filtTracks = filterTracksInJets(tracks, branch, 0.4);
      filtjetClusters = clusterTracksInTime(filtTracks, branch, 3.0,
                                            useSmearedTimes, checkValidTimes, 20.0,
                                            /*useCone=*/true, useZ0);
    }

    // ── F. Fill denominator histograms ──────────────────────────────────────
    // TEST_MISCL denominator is deferred to step H where cluster purity is known.
    for (auto& [score, analysis] : analyses)
      if (!score.requiresPurity)
        fillTotals(analysis, ev);

    // ── G. Select best cluster for each score ───────────────────────────────
    auto chosen = selectClusters(clusters, tracks, branch, analyses, filtjetClusters);
    if (DEBUG) std::cout << "Chose clusters\n";

    // ── H. Per-score histogram filling ──────────────────────────────────────
    int    returnCode = 0;
    double returnVal  = -1.;
    bool passesMine = false, passesMiscl = false, misclInDenominator = false;

    for (auto& [score, analysis] : analyses) {
      if (DEBUG) std::cout << "Filling: " << toString(score) << '\n';

      // Skip scores with invalid or missing cluster data
      if (branch->recoVtxValid[0] == 0 && score == Score::HGTD)  continue;
      if (clusters.empty() && !score.usesOwnCollection)           continue;
      if (!chosen.count(score) && score != Score::HGTD)           continue;

      Cluster& scored    = chosen[score];
      double t           = (score == Score::HGTD) ? branch->recoVtxTime[0]
                                                   : scored.values.at(0);
      double purity      = clusters.empty() ? 0.0 : scored.purity;
      double diff        = t - branch->truthVtxTime[0];

      // Expose time for the caller's event-display collection
      if (score == Score::TRKPTZ || score.requiresPurity) returnVal = t;

      // Inclusive timing resolution and purity (purity-gated for TEST_MISCL)
      if (!score.requiresPurity || purity > 0.75) {
        analysis.inclusiveReso->Fill(diff);
        analysis.inclusivePurity->Fill(purity);
      }

      // Efficiency check: base timing window, plus optional score threshold
      bool passes = scored.passEfficiency(branch);
      if (score.hasThreshold() && passes)
        passes = scored.scores.at(score) > score.threshold;

      // TEST_MISCL: restrict both denominator and numerator to pure clusters
      if (score.requiresPurity) {
        if (purity > 0.75) {
          misclInDenominator = true;
          fillTotals(analysis, ev);
        } else {
          passes = false;  // impure cluster: skip numerator entirely
        }
      }

      if (passes) {
        fillPasses(analysis, ev);
        if (score == Score::TRKPTZ)   passesMine   = true;
        if (score.requiresPurity)     passesMiscl  = true;
      }

      // 2D timing residual and purity distributions
      // (TEST_MISCL: only fill events that entered the denominator)
      if (!score.requiresPurity || misclInDenominator) {
        fillDiffs   (analysis, ev, diff);
        fillPurities(analysis, ev, purity);
      }
    }

    // returnCode == 2: event was in TEST_MISCL denominator but failed timing window
    if (misclInDenominator && !passesMiscl) returnCode = 2;

    return {returnCode, returnVal};
  }
}
#endif // EVENT_PROCESSING_H
