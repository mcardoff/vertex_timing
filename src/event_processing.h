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
//     4. filterTracksInJets         — keep only tracks within ΔR of a jet
//     5. getAssociatedTracks        — full HGTD-acceptance + pTV selection
//     6. selectClusters             — choose best cluster for every active score
//     7. processEventData           — main per-event analysis orchestrator
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

    double nsigmaPrim = std::abs(trkZ0 - vxZ) / std::sqrt(Z0_VAR_INFLATION * trkZVar + vxZVar);
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

    // WAVES: precision-weighted mean using only constituent tracks that fall within
    // dR < 0.4 of at least one forward reco jet (pT > MIN_JET_PT, |η| in HGTD range).
    // The cluster is still selected by the WAVeS score; this in-jet timing refinement
    // purifies the time by dropping absorbed PU tracks, worth ~+1% efficiency.
    // WAVES_MISCL / WAVES_MISAS share this path so the oracle rows report the same
    // time as the WAVES row they gate.
    if (score == Score::WAVES || score == Score::WAVES_MISCL || score == Score::WAVES_MISAS) {
      // Collect (eta, phi) of qualifying forward reco jets — no truth matching
      std::vector<std::pair<double,double>> hsJets;
      const int nJets = (int)branch->topoJetPt.GetSize();
      for (int j = 0; j < nJets; ++j) {
        double jEta = branch->topoJetEta[j];
        if (branch->topoJetPt[j] < MIN_JET_PT) continue;
        if (std::abs(jEta) < MIN_ABS_ETA_JET || std::abs(jEta) > MAX_ABS_ETA_JET) continue;
        hsJets.push_back({jEta, branch->topoJetPhi[j]});
      }
      double sumW = 0.0, sumWT = 0.0;
      for (size_t i = 0; i < trackIndices.size(); ++i) {
        int    idx     = trackIndices[i];
        double trkEta  = branch->trackEta[idx];
        double trkPhi  = branch->trackPhi[idx];
        bool   inJet   = false;
        for (auto& [jeta, jphi] : hsJets) {
          double deta = trkEta - jeta;
          double dphi = TVector2::Phi_mpi_pi(trkPhi - jphi);
          if (std::hypot(deta, dphi) < 0.4) { inJet = true; break; }
        }
        if (!inJet) continue;
        double t = allTimes[i];
        double s = (idealRes > 0.0) ? idealRes
                                    : static_cast<double>(branch->trackTimeRes[idx]);
        if (s <= 0.0) continue;
        double w = 1.0 / (s * s);
        sumW += w;  sumWT += w * t;
      }
      if (sumW > 0.0) return sumWT / sumW;
      return values[0];  // fallback: no tracks survived jet filter
    }

    return values[0];
  }

  // ---------------------------------------------------------------------------
  // 3c. Cluster::calculatePurity — out-of-line implementation
  //   For most scores: returns the pre-computed this->purity (HS ΣpT fraction
  //   of the full cluster).  For WAVES: re-evaluates purity using only the
  //   constituent tracks within dR < 0.4 of a forward reco jet, so the reported
  //   purity matches the tracks actually used for the in-jet timing.
  //   WAVES_MISCL / WAVES_MISAS intentionally use the default full-cluster purity
  //   so their purity gates stay directly comparable to TEST_MISCL / TEST_MISAS.
  // ---------------------------------------------------------------------------
  inline double Cluster::calculatePurity(Score score, BranchPointerWrapper* branch) const {
    if (score != Score::WAVES)
      return this->purity;

    // Collect (eta, phi) of qualifying forward reco jets — no truth matching.
    std::vector<std::pair<double,double>> hsJets;
    const int nJets = (int)branch->topoJetPt.GetSize();
    for (int j = 0; j < nJets; ++j) {
      double jEta = branch->topoJetEta[j];
      if (branch->topoJetPt[j] < MIN_JET_PT) continue;
      if (std::abs(jEta) < MIN_ABS_ETA_JET || std::abs(jEta) > MAX_ABS_ETA_JET) continue;
      hsJets.push_back({jEta, branch->topoJetPhi[j]});
    }

    double sumPt = 0.0, hsPt = 0.0;
    for (int idx : trackIndices) {
      double trkEta = branch->trackEta[idx];
      double trkPhi = branch->trackPhi[idx];
      bool inJet = false;
      for (auto& [jeta, jphi] : hsJets) {
        double deta = trkEta - jeta;
        double dphi = TVector2::Phi_mpi_pi(trkPhi - jphi);
        if (std::hypot(deta, dphi) < 0.4) { inJet = true; break; }
      }
      if (!inJet) continue;
      double pt = branch->trackPt[idx];
      sumPt += pt;
      // HS tracks have truthVtx index 0 — mirrors calcPurity() convention
      if (branch->trackToTruthvtx[idx] == 0) hsPt += pt;
    }
    if (sumPt > 0.0) return hsPt / sumPt;
    return this->purity;  // fallback: no in-jet tracks found
  }

  // ---------------------------------------------------------------------------
  // 4. filterTracksInJets
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
  //   Returns 0.0 when no HS tracks with valid HGTD time exist, so such events
  //   are excluded from the TEST_MISAS denominator rather than counted as clean.
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
      if (pull < TRUTH_PULL_CUT) num += pT;
    }
    return (denom > 0.0) ? static_cast<float>(num / denom) : 0.0f;
  }

  // ---------------------------------------------------------------------------
  // 5b2. filterHSTracks
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
  // 7. selectClusters
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

    // HGTD: simultaneous clustering on real HGTD times.
    // HGTD selects by timing proximity to the reco vertex, not by TRKPTZ score, so
    // MIN_CLUSTER_TRACKS does not apply — any non-empty cluster is valid.
    if (branch->recoVtxValid[0] == 1 && analyses.count(Score::HGTD)) {
      auto col = clusterTracksInTime(tracks, branch, 3.0,
                                     false, true, -1,
                                     ClusteringMethod::SIMULTANEOUS, false, false, true);
      if (!col.empty())
        chosen[Score::HGTD.id] = chooseHGTDCluster(col, branch);
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
  // 8. processEventData
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
  //   code           — -1: rejected by selection; 0: normal; 2: MISCL fail;
  //                    3: MISAS fail
  //   time           — TRKPTZ-selected (or TEST_MISCL) cluster time
  //   nFwdHS         — n forward HS tracks (3σ counting step)
  //   trkptzPass     — true if TRKPTZ passed the PASS_SIGMA timing window
  //   tRefinedPass   — true if T_REFINED passed the PASS_SIGMA timing window
  //   misclInDenom   — true if event entered the TEST_MISCL denominator (cluster purity > 75%)
  //   misasInDenom   — true if event entered the TEST_MISAS denominator (hsTimingPurity ≥ 95%)
  //   misclPass      — true if event was in TEST_MISCL denominator AND passed
  //   misasPass      — true if event was in TEST_MISAS denominator AND passed
  // ---------------------------------------------------------------------------
  struct EventResult {
    int    code          = -1;
    double time          = -1.0;
    int    nFwdHS        =  0;
    bool   trkptzPass    = false;
    bool   tRefinedPass  = false;
    bool   misclInDenom  = false;
    bool   misasInDenom  = false;
    bool   misclPass     = false;
    bool   misasPass     = false;
    // Combined cluster quality for the TRKPTZ-selected cluster:
    // (1 - clusPuFrac)² * clamp(avgNHGTD/2, 0, 1) * sigmaTFactor
    double clusQuality   =  0.0;
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

    const bool NEEDS_PURITY = analyses.count(Score::TEST_MISCL) > 0 ||
                              analyses.count(Score::WAVES_MISCL) > 0;
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
        analysis.fillTotals(ev);

    // ── G. Select best cluster for each score ───────────────────────────────
    auto chosen = selectClusters(clusters, tracks, branch, analyses, auxCollections);
    if (DEBUG) std::cout << "Chose clusters\n";

    // ── G½. HS timing purity for TEST_MISAS gate ────────────────────────────
    float hsTimingPurity = 1.0f;
    if (analyses.count(Score::TEST_MISAS) || analyses.count(Score::WAVES_MISAS))
      hsTimingPurity = calcHSTimingPurity(tracks, branch);

    // ── H. Per-score histogram filling ──────────────────────────────────────
    int    returnCode = 0;
    double returnVal  = -1.;
    double trkptzClusQuality = 0.0;  // quality of TRKPTZ-selected cluster; exported in EventResult
    bool passesMine = false, passesTRefined = false;
    bool passesMiscl = false, passesMisas = false;
    bool misclInDenominator = false, misasInDenominator = false;
    bool perfEvtInDenominator = false;

    for (auto& [score, analysis] : analyses) {
      if (DEBUG) std::cout << "Filling: " << score.toString() << '\n';

      // Skip scores with invalid or missing cluster data
      if (branch->recoVtxValid[0] == 0 && score == Score::HGTD)   continue;
      if (clusters.empty() && !score.usesOwnCollection)           continue;
      if (!chosen.count(score.id))                                continue;

      Cluster& scored    = chosen.at(score.id);
      double iResH       = useSmearedTimes ? IDEAL_TRACK_RES : -1.0;
      double t           = scored.calculateTime(score, branch, iResH);
      scored.values[0]   = t;   // keep passEfficiency consistent with diff
      double purity      = clusters.empty() ? 0.0 : scored.calculatePurity(score, branch);
      double diff        = t - branch->truthVtxTime[0];

      // Avg nHGTD hits per track in the selected cluster (computed once, used in
      // multiple fill sites below: inclusive histograms, PlotObj total/pass/diff/purity)
      double avgNHGTD = 0.0;
      if (!scored.trackIndices.empty()) {
        for (int trk : scored.trackIndices)
          avgNHGTD += branch->trackHgtdHits[trk];
        avgNHGTD /= static_cast<double>(scored.trackIndices.size());
      }

      // Cluster-level PU fraction by track count (consistent with event-level puRatio).
      double clusPuFrac = 0.0;
      if (!scored.trackIndices.empty()) {
        int nPU = 0;
        for (int trk : scored.trackIndices)
          if (branch->trackToTruthvtx[trk] != 0) nPU++;
        clusPuFrac = static_cast<double>(nPU) / static_cast<double>(scored.trackIndices.size());
      }

      // Cluster timing uncertainty — hoisted so all fill sites (total/pass/diff/purity)
      // can reach it without re-accessing the sigmas vector, and so it can feed clusQuality
      double clusSigmaT = scored.sigmas.at(0);

      // σ_t factor: linear roll-off between FLOOR and CEIL.
      // Well-determined clusters (σ_t ≤ floor) get factor 1; poorly-determined ones
      // (σ_t ≥ ceil) get 0. This pushes high-Q outliers — clusters that look clean
      // by composition but report large self-uncertainty — toward LOW quality.
      //
      // Note: nConstituents was tested as a separate fourth factor; it produced identical
      // results because σ_t already encodes √N (averaging power) plus per-track resolution
      // variance and track-time disagreement. nConstituents is fully subsumed by σ_t.
      double sigmaTFactor = std::clamp(
        (CLUS_SIGMA_T_FACTOR_CEIL - clusSigmaT) /
          (CLUS_SIGMA_T_FACTOR_CEIL - CLUS_SIGMA_T_FACTOR_FLOOR),
        0.0, 1.0);

      // Combined cluster quality: (1 - PU frac)² × clamp(avgNHGTD/2, 0, 1) × σ_t factor
      // Ranges 0→1; captures PU contamination, nhit, and timing certainty in one scalar.
      // Quadratic puFrac because the linear form let high-puFrac clusters (~0.7-0.8) survive
      // into the MID tier — these are the ~30 ps "middle background" clusters. σ_t doesn't
      // catch them either because PU tracks at the same z have correlated (wrong) times.
      double puFracFactor = (1.0 - clusPuFrac) * (1.0 - clusPuFrac);
      double clusQuality = puFracFactor * std::min(avgNHGTD / 2.0, 1.0) * sigmaTFactor;
      if (score == Score::TRKPTZ) trkptzClusQuality = clusQuality;

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
        if (score == Score::TEST_MISAS || score == Score::WAVES_MISAS) {
          // Gate on event-level HS timing purity: 100% of HS pT must have |pull|<3σ
          passesGate = (hsTimingPurity >= 0.95f);
        } else if (score == Score::PERF_EVT) {
          // MISCL ∧ MISAS: pure cluster AND all event-level HS tracks correctly timed
          passesGate = (purity > 0.75f) && (hsTimingPurity >= 0.95f);
        } else {
          float purityCut = score.hasThreshold() ? score.threshold : 0.75f;
          passesGate = (purity > purityCut);
        }
        inDenominator = passesGate;
        if (passesGate) {
          if (score == Score::TEST_MISCL) misclInDenominator = true;
          if (score == Score::TEST_MISAS) misasInDenominator = true;
          if (score == Score::PERF_EVT)   perfEvtInDenominator = true;
          analysis.fillTotals(ev);
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

        // PlotObj denominators for cluster-level variables (filled here rather than at
        // the pre-loop fillTotals call since these values require the selected cluster)
        analysis.ptrNhit->fillTotal(avgNHGTD);
        analysis.ptrClusPuFrac->fillTotal(clusPuFrac);
        analysis.ptrClusQuality->fillTotal(clusQuality);
        if (clusSigmaT > 0.0) analysis.ptrClusSigmaT->fillTotal(clusSigmaT);
        if (avgNHGTD < 1.5)
          fillResoStack(analysis.inclusiveResoNhit1Sig.get(),
                        analysis.inclusiveResoNhit1Mix.get(),
                        analysis.inclusiveResoNhit1Bkg.get());
        else if (avgNHGTD < 2.5)
          fillResoStack(analysis.inclusiveResoNhit2Sig.get(),
                        analysis.inclusiveResoNhit2Mix.get(),
                        analysis.inclusiveResoNhit2Bkg.get());
        else
          fillResoStack(analysis.inclusiveResoNhit3pSig.get(),
                        analysis.inclusiveResoNhit3pMix.get(),
                        analysis.inclusiveResoNhit3pBkg.get());

        // Cluster quality-binned inclusive reso (two tiers: HIGH Q≥0.5, LOW Q<0.5)
        if (clusQuality >= CLUS_QUALITY_SPLIT) {
          fillResoStack(analysis.inclusiveResoClusQHighSig.get(),
                        analysis.inclusiveResoClusQHighMix.get(),
                        analysis.inclusiveResoClusQHighBkg.get());
        } else {
          fillResoStack(analysis.inclusiveResoClusQLowSig.get(),
                        analysis.inclusiveResoClusQLowMix.get(),
                        analysis.inclusiveResoClusQLowBkg.get());
        }

        // Pull distribution: Δt / σ_cluster — skipped when σ is unphysical (≤ 0)
        if (clusSigmaT > 0.0) {
          double sigmaT = clusSigmaT;
          double pull = diff / sigmaT;
          auto fillPullStack = [&](TH1D* sig, TH1D* mix, TH1D* bkg) {
            if      (purity > 0.75)  sig->Fill(pull);
            else if (purity >= 0.50) mix->Fill(pull);
            else                     bkg->Fill(pull);
          };
          fillPullStack(analysis.inclusivePullSig.get(),
                        analysis.inclusivePullMix.get(),
                        analysis.inclusivePullBkg.get());
          // Per-quality-tier pull — two tiers matching the reso split
          if (clusQuality >= CLUS_QUALITY_SPLIT)
            fillPullStack(analysis.inclusivePullClusQHighSig.get(),
                          analysis.inclusivePullClusQHighMix.get(),
                          analysis.inclusivePullClusQHighBkg.get());
          else
            fillPullStack(analysis.inclusivePullClusQLowSig.get(),
                          analysis.inclusivePullClusQLowMix.get(),
                          analysis.inclusivePullClusQLowBkg.get());
          if (ev.nForwardTrackHS <= 5)
            fillPullStack(analysis.inclusivePullLowTrackSig.get(),
                          analysis.inclusivePullLowTrackMix.get(),
                          analysis.inclusivePullLowTrackBkg.get());
        }
      }

      if (passes) {
        analysis.fillPasses(ev);
        analysis.ptrNhit->fillPass(avgNHGTD);
        analysis.ptrClusPuFrac->fillPass(clusPuFrac);
        analysis.ptrClusQuality->fillPass(clusQuality);
        if (clusSigmaT > 0.0) analysis.ptrClusSigmaT->fillPass(clusSigmaT);
        if (score == Score::TRKPTZ)     passesMine    = true;
        if (score == Score::T_REFINED)  passesTRefined = true;
        if (score == Score::TEST_MISCL) passesMiscl   = true;
        if (score == Score::TEST_MISAS) passesMisas   = true;
      }

      // 2D timing residual and purity distributions
      // (requiresPurity scores: only fill events that entered the denominator;
      //  inDenominator is always true for non-gated scores)
      if (inDenominator) {
        analysis.fillDiffs   (ev, diff);
        analysis.fillPurities(ev, purity);
        analysis.ptrNhit->fillDiff       (avgNHGTD,   diff);
        analysis.ptrNhit->fillPurity     (avgNHGTD,   purity);
        analysis.ptrClusPuFrac->fillDiff  (clusPuFrac, diff);
        analysis.ptrClusPuFrac->fillPurity(clusPuFrac, purity);
        analysis.ptrClusQuality->fillDiff  (clusQuality, diff);
        analysis.ptrClusQuality->fillPurity(clusQuality, purity);
        if (clusSigmaT > 0.0) {
          analysis.ptrClusSigmaT->fillDiff  (clusSigmaT, diff);
          analysis.ptrClusSigmaT->fillPurity(clusSigmaT, purity);
        }
        // 2D heatmaps: mean |Δt| and σ(Δt) as a function of (clusPuFrac, avgNHGTD)
        analysis.prof2dPuFracVsNhit     ->Fill(clusPuFrac, avgNHGTD, std::abs(diff));
        analysis.prof2dPuFracVsNhitSigma->Fill(clusPuFrac, avgNHGTD, diff);

        // In-time PU diagnostic: find the dominant PU truth vertex contributing tracks
        // to this cluster (mode of trackToTruthvtx among PU tracks), then fill
        // Δt(cluster − HS truth) vs Δz(dominant PU vtx − HS vtx).
        // The expected slope σ_t / σ_z ≈ 3.5 ps/mm confirms in-time PU origin.
        if (clusPuFrac > 0.0) {
          std::unordered_map<int, int> puVtxCount;
          for (int trk : scored.trackIndices) {
            int vtx = branch->trackToTruthvtx[trk];
            if (vtx > 0) puVtxCount[vtx]++;   // only PU vertices (HS is index 0)
          }
          if (!puVtxCount.empty()) {
            int dominantPUVtx = -1, maxCount = 0;
            for (const auto& [vtx, n] : puVtxCount) {
              if (n > maxCount) { maxCount = n; dominantPUVtx = vtx; }
            }
            if (dominantPUVtx > 0) {
              double dzPU = branch->truthVtxZ[dominantPUVtx] - branch->truthVtxZ[0];
              analysis.dtClusterVsDzPU->Fill(dzPU, diff);
              if (clusQuality >= CLUS_QUALITY_SPLIT)
                analysis.dtClusterVsDzPUHigh->Fill(dzPU, diff);
              else
                analysis.dtClusterVsDzPULow->Fill(dzPU, diff);
            }
          }
        }
      }
    }

    // returnCode == 2: event was in TEST_MISCL denominator but failed timing window
    if (misclInDenominator && !passesMiscl) returnCode = 2;

    // returnCode == 3: event has clean HS timing (hsTimingPurity > 0.75, in TEST_MISAS
    // denominator) but TRKPTZ still fails — the failure is not caused by timing
    // misassignment.  Only meaningful in the HGTD scenario where TEST_MISAS is active.
    if (analyses.count(Score::TEST_MISAS) &&
        misasInDenominator && !passesMine) returnCode = 3;

    return {returnCode, returnVal, ev.nForwardTrackHS,
            passesMine, passesTRefined,
            misclInDenominator, misasInDenominator,
            passesMiscl, passesMisas,
            trkptzClusQuality};
  }
}
#endif // EVENT_PROCESSING_H
