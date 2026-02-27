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
//     5. filterCaloTracks           — reject tracks inconsistent with calo time
//     6. filterTracksInJets         — keep only tracks within ΔR of a jet
//     7. getAssociatedTracks        — full HGTD-acceptance + pTV selection
//     8. processEventData           — main per-event analysis function
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
  // 5. filterCaloTracks
  //   Retains only tracks whose HGTD time is within significanceCut standard
  //   deviations of a calorimeter-smeared event time.  The calorimeter time
  //   is drawn from Gaus(truthVtxTime[0], caloRes) once per call.  Tracks
  //   where |t_trk - t_calo| / sqrt(σ_trk² + caloRes²) > significanceCut
  //   are removed.  Used to build the FILT60 and FILT90 filtered-track
  //   collections.
  // ---------------------------------------------------------------------------
  std::vector<int> filterCaloTracks(
    const std::vector<int>& tracks,
    BranchPointerWrapper *branch,
    double caloRes,
    double signficanceCut
  ) {
    std::vector<int> output;
    double caloTime = gRandom->Gaus(branch->truthVtxTime[0], caloRes);
    for (const auto& trk: tracks) {
      double trkTime = branch->trackTime[trk];
      double trkTimeRes = branch->trackTimeRes[trk];
      double caloDiff = std::abs(trkTime-caloTime);
      double caloSigma = std::hypot(trkTimeRes, caloRes);
      double nsigma = caloDiff/caloSigma;
      if (nsigma > signficanceCut)
	continue;
      
      output.push_back(trk);
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 6. filterTracksInJets
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
    for (const auto& trk: tracks) {
      double
	trkEta = branch->trackEta[trk],
	trkPhi = branch->trackPhi[trk];
      double minDR = 1e6;
      for (size_t jetIdx=0; jetIdx < branch->topoJetEta.GetSize(); ++jetIdx) {
	if (branch->topoJetPt[jetIdx] < MIN_JETPT) continue;
	double
	  jetEta = branch->topoJetEta[jetIdx],
	  jetPhi = branch->topoJetPhi[jetIdx];
	double
	  deta = jetEta-trkEta,
	  dphi = TVector2::Phi_mpi_pi(jetPhi - trkPhi);
	double thisDR = std::hypot(deta, dphi);
	if (thisDR < minDR)
	  minDR = thisDR;
      }
      if (minDR < minDRCut)
	output.push_back(trk);
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 7. getAssociatedTracks
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
  // 8. processEventData
  //   Main per-event analysis function.  Called once per event per timing
  //   scenario (HGTD / IdealRes / IdealEff).  Pipeline:
  //     a) Basic event selection (passBasicCuts, passJetPtCut).
  //     b) Track selection at 3σ and counting statistics (nForwardJet,
  //        nForwardTrack, nForwardTrackHS, nForwardTrackPU, puRatio).
  //     c) Tighten track association to MAX_NSIGMA if different from 3σ.
  //     d) Cluster the main track list (cone algorithm).
  //     e) Optionally build FILT60, FILT90, FILTJET filtered collections
  //        (only when the corresponding scores are present in analyses).
  //     f) Run HGTD simultaneous clustering on the unfiltered list (for
  //        the HGTD score).
  //     g) Choose the best cluster for every active score.
  //     h) Fill denominator (fillTotal) and numerator (fillPass) histograms
  //        for each score; apply special purity gate for TEST_MISCL.
  //     i) Fill timing residual (fillDiff) and purity (fillPurity) 2D hists.
  //
  //   Return value: pair<returnCode, returnVal>
  //     returnCode == 0  normal event
  //     returnCode == -1 event rejected by selection
  //     returnCode == 2  event passed TEST_MISCL denominator (purity > 0.75)
  //                      but failed the timing window — useful for event
  //                      display collection
  //     returnVal        TRKPTZ-selected cluster time (or TEST_MISCL time)
  // ---------------------------------------------------------------------------
  std::pair<int,double> processEventData(
    BranchPointerWrapper *branch,
    bool useSmearedTimes,
    bool checkValidTimes,
    bool useZ0, // whether or not to use 2D Clustering
    std::map<Score,AnalysisObj>& analyses
  ) {
    int returnCode = 0;
    double returnVal = -1.;
    // check if vertex selection is correct & number of jets
    if (not branch->passBasicCuts()) return std::make_pair(-1,1.);

    // check if there is one forward jet with pt > 30 GeV
    if (not branch->passJetPtCut()) return std::make_pair(-1, 1.);

    int nForwardJet=0;
    branch->countForwardJets(nForwardJet);
    
    // Scan the track array once at the looser cut (3σ) for counting stats,
    // then tighten to MAX_NSIGMA by filtering the already-selected subset.
    std::vector<int> tracks = getAssociatedTracks(branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);

    int nForwardTrack=0, nForwardTrackHS=0, nForwardTrackPU=0;
    branch->countForwardTracks(nForwardTrack,nForwardTrackHS,nForwardTrackPU,tracks, checkValidTimes);

    // if (not branch->pass_forward_hs_tracks(nForwardTrack_HS)) return std::make_pair(-1,1.);

    double puRatio = (double)nForwardTrackPU / (double)nForwardTrack;
    double recoZ = branch->recoVtxZ[0];

    auto effFillValFjet    = folded(nForwardJet    , (int)FOLD_FJET) ;
    auto effFillValTrack   = folded(nForwardTrack  , (int)FOLD_TRACK);
    auto effFillValHSTrack = folded(nForwardTrackHS, (int)FOLD_HS_TRACK);
    auto effFillValPUTrack = folded(nForwardTrackPU, (int)FOLD_PU_TRACK);
    auto effFillValPURatio = puRatio >= 1.0 ? 0.99 : puRatio;
    double effFillValVtxDz = std::abs(branch->recoVtxZ[0] - branch->truthVtxZ[0]);

    // Tighten association: keep only tracks that also pass MAX_NSIGMA.
    // Filtering the 3σ subset avoids a second full scan of the track array.
    if (MAX_NSIGMA != 3.0)
      tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                                  [&](int trk) {
                                    return !passTrackVertexAssociation(
                                        trk, 0, branch, MAX_NSIGMA);
                                  }),
                   tracks.end());    

    if (DEBUG) {
      std::cout << "nForwardJet = " << nForwardJet << '\n';
      std::cout << "nForwardTrack = " << nForwardTrack << '\n';
      std::cout << "nForwardTrack_HS = " << nForwardTrackHS << '\n';
      std::cout << "nForwardTrack_PU = " << nForwardTrackPU << '\n';
    }

    // Only pay for purity calculation when TEST_MISCL is actually in this map
    const bool NEEDS_PURITY = analyses.count(Score::TEST_MISCL) > 0;

    std::vector<Cluster> clusters =
      clusterTracksInTime(
        tracks, branch, 3.0,
	useSmearedTimes, checkValidTimes, 10.0,
	true, useZ0, NEEDS_PURITY);

    std::vector<Cluster> filt60Clusters, filt90Clusters, filtjetClusters;

    if (analyses.count(Score::FILT60)) {
      std::vector<int> filteredTracks = filterCaloTracks(tracks, branch, 60, 2.0);
      filt60Clusters =
	clusterTracksInTime(
          filteredTracks, branch, 3.0,
          useSmearedTimes, checkValidTimes, 20.0,
	  true, useZ0);
    }

    if (analyses.count(Score::FILT90)) {
      std::vector<int> filteredTracks = filterCaloTracks(tracks, branch, 90, 2.0);
      filt90Clusters =
	clusterTracksInTime(
          filteredTracks, branch, 3.0,
          useSmearedTimes, checkValidTimes, 20.0,
	  true, useZ0);
    }

    if (analyses.count(Score::FILTJET)) {
      std::vector<int> filteredTracks = filterTracksInJets(tracks, branch, 0.4);
      filtjetClusters =
	clusterTracksInTime(
          filteredTracks, branch, 3.0,
          useSmearedTimes, checkValidTimes, 20.0,
	  true, useZ0);
    }

    if (DEBUG) std::cout << "LEFT CLUSTERING\n";
    for (auto& [score,analysis]: analyses) {
      // TEST_MISCL denominator is restricted to events where the TRKPTZ-selected
      // cluster has purity > 0.75 — fillTotal is deferred to the per-score loop
      // below where the chosen cluster (and its purity) is already available.
      if (score == Score::TEST_MISCL) continue;
      analysis["fjet"]->     fillTotal(effFillValFjet   );
      analysis["vtx_dz"]->   fillTotal(effFillValVtxDz  );
      analysis["ftrack"]->   fillTotal(effFillValTrack  );
      analysis["pu_frac"]->  fillTotal(effFillValPURatio);
      analysis["hs_track"]-> fillTotal(effFillValHSTrack);
      analysis["pu_track"]-> fillTotal(effFillValPUTrack);
    }
    
    std::map<Score,Cluster> chosen;
    if (clusters.size() != 0)
      chosen = chooseCluster(clusters, branch);    
    if (DEBUG) std::cout << "Chose My Clusters\n";

    if (not filt60Clusters.empty())
      chosen[FILT60] = chooseCluster(filt60Clusters, TRKPTZ);

    if (not filt90Clusters.empty())
      chosen[FILT90] = chooseCluster(filt90Clusters, TRKPTZ);

    if (not filtjetClusters.empty())
      chosen[FILTJET] = chooseCluster(filtjetClusters, TRKPTZ);

    // run HGTD Clustering (simultaneous)
    std::vector<Cluster> clustersHGTD =
      clusterTracksInTime(tracks, branch, 3.0, false, true , -1, false, false);

    if (branch->recoVtxValid[0] == 1 and analyses.count(Score::HGTD))
      chosen[HGTD] = chooseHGTDCluster(clustersHGTD, branch);

    if (DEBUG) std::cout << "Chose clusters\n";

    bool hasPassingCluster = false, passesMine = false, passesMiscl = false, misclInDenominator = false;

    for (auto& [score, analysis] : analyses) {
      if (DEBUG)
        std::cout << "Filling Scores: " << toString(score) << '\n';
      if (branch->recoVtxValid[0] == 0 and score == Score::HGTD)
	continue;

      if (clusters.size() == 0 and score != Score::HGTD)
	continue;

      if (!chosen.count(score) and score != Score::HGTD)
	continue;

      if (DEBUG)
        std::cout << "Attempting to access chosen[" << toString(score) << "]\n";
      Cluster scored = chosen[score];
      if (DEBUG)
        std::cout << "Attempting to access scores values, size: " << scored.values.size() << '\n';
      double scoreBasedTime = score == Score::HGTD ? branch->recoVtxTime[0] : scored.values.at(0);
      if (DEBUG)
        std::cout << "Attempting to access purity\n";
      double clusterPurity = clusters.size() != 0 ? scored.purity : 0;
      if (DEBUG)
        std::cout << "Purity: " << clusterPurity << '\n';
      double diff = scoreBasedTime - branch->truthVtxTime[0];
      if (DEBUG)
        std::cout << "Diff: " << diff << '\n';

      if (score == Score::TRKPTZ)
	returnVal = scoreBasedTime;
      if (score == Score::TEST_MISCL)
        returnVal = scoreBasedTime;
      
      // For TEST_MISCL, only fill inclusive histos when the cluster is pure (purity > 0.75)
      if (score != Score::TEST_MISCL || scored.purity > 0.75) {
        analysis.inclusiveReso->Fill(diff);
        analysis.inclusivePurity->Fill(clusterPurity);
      }

      // Check if this score passes the efficiency test
      // For TESTML:     require passEfficiency AND ML score > 0.5
      // For TEST_MISCL: uses TRKPTZ scoring; both the denominator (fillTotal) and
      //                 numerator (fillPass) are restricted to events where the
      //                 selected cluster has purity > 0.75.
      bool passesEfficiency = scored.passEfficiency(branch);
      if (score == Score::TESTML && passesEfficiency) {
        passesEfficiency = scored.scores.at(Score::TESTML) > 0.3;
      }
      if (score == Score::TEST_MISCL) {
        // Only count this event in the denominator if the selected cluster is pure
        if (scored.purity > 0.75) {
          misclInDenominator = true;
          analysis["fjet"]->     fillTotal(effFillValFjet   );
          analysis["vtx_dz"]->   fillTotal(effFillValVtxDz  );
          analysis["ftrack"]->   fillTotal(effFillValTrack  );
          analysis["pu_frac"]->  fillTotal(effFillValPURatio);
          analysis["hs_track"]-> fillTotal(effFillValHSTrack);
          analysis["pu_track"]-> fillTotal(effFillValPUTrack);
        } else {
          passesEfficiency = false;  // impure cluster: skip fillPass entirely
        }
      }

      if (passesEfficiency) {
	analysis["fjet"]->     fillPass(effFillValFjet   );
	analysis["vtx_dz"]->   fillPass(effFillValVtxDz  );
	analysis["ftrack"]->   fillPass(effFillValTrack  );
	analysis["pu_frac"]->  fillPass(effFillValPURatio);
	analysis["hs_track"]-> fillPass(effFillValHSTrack);
	analysis["pu_track"]-> fillPass(effFillValPUTrack);
	if (score == Score::PASS) {
	  hasPassingCluster = true;
	}
	if (score == Score::TRKPTZ)
	  passesMine = true;
	if (score == Score::TEST_MISCL)
	  passesMiscl = true;
      }

      // fill diff hists — for TEST_MISCL only fill when in the denominator
      if (score != Score::TEST_MISCL || misclInDenominator) {
        analysis["fjet"]->     fillDiff(nForwardJet    , diff);
        analysis["vtx_dz"]->   fillDiff(effFillValVtxDz, diff);
        analysis["ftrack"]->   fillDiff(nForwardTrack  , diff);
        analysis["pu_frac"]->  fillDiff(puRatio        , diff);
        analysis["hs_track"]-> fillDiff(nForwardTrackHS, diff);
        analysis["pu_track"]-> fillDiff(nForwardTrackPU, diff);
      }

      // fill purities — same gate for TEST_MISCL
      if (score != Score::TEST_MISCL || misclInDenominator) {
        analysis["fjet"]->     fillPurity(nForwardJet    , clusterPurity);
        analysis["vtx_dz"]->   fillPurity(effFillValVtxDz, clusterPurity);
        analysis["ftrack"]->   fillPurity(nForwardTrack  , clusterPurity);
        analysis["pu_frac"]->  fillPurity(puRatio        , clusterPurity);
        analysis["hs_track"]-> fillPurity(nForwardTrackHS, clusterPurity);
        analysis["pu_track"]-> fillPurity(nForwardTrackPU, clusterPurity);
      }
    }
    // Flag code 2 only for events that entered the TEST_MISCL denominator
    // (pure cluster, purity > 0.75) but did NOT pass the timing window.
    if (misclInDenominator && !passesMiscl) returnCode = 2;
    
    return std::make_pair(returnCode, returnVal);
  }
}
#endif // EVENT_PROCESSING_H
