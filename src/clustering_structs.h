#ifndef STRUCTS_H
#define STRUCTS_H

// ---------------------------------------------------------------------------
// clustering_structs.h
//   Core data structures shared across the entire analysis pipeline.
//   Two structs are defined here:
//
//   BranchPointerWrapper
//     Owns all TTreeReaderArray / TTreeReaderValue handles for a single
//     TTreeReader session.  Also hosts all per-event selection predicates
//     and kinematic helper methods so that event-loop code in
//     event_processing.h stays concise.
//
//   Cluster
//     Represents one time-space cluster produced by the clustering
//     algorithms in clustering_functions.h.  Stores the weighted-mean
//     time (and optionally z₀), per-constituent track indices, all raw
//     times (for spread diagnostics), and a score map keyed on Score.
//     Also hosts the ML feature extraction and score-update logic so
//     that the cluster is self-contained.
// ---------------------------------------------------------------------------

#include "clustering_includes.h"
#include "clustering_constants.h"
#include <cstdlib>

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // BranchPointerWrapper
  //   Aggregates every TTreeReaderArray / TTreeReaderValue needed by the
  //   analysis into a single object that is passed by pointer throughout the
  //   event loop.  The constructor binds each member to its branch name so
  //   that branch-name strings appear only once in the codebase.
  //
  //   Selection methods are grouped as follows:
  //     Basic event cuts         — passBasicCuts, passJetPtCut
  //     Track / jet counters     — countForwardJets, countForwardTracks
  //     Testing-only HS filter   — passForwardHsTracks
  //     Track scoring helpers    — calcJetptDRScore, calcTrkptDRScore
  //     VBF signal region cuts   — see dedicated sub-section below
  // ---------------------------------------------------------------------------
  struct BranchPointerWrapper {
    TTreeReader& reader;

    TTreeReaderValue<float> weight;

    TTreeReaderArray<float> trackZ0;
    TTreeReaderArray<float> trackD0;
    TTreeReaderArray<float> trackPt;
    TTreeReaderArray<float> trackQP;
    TTreeReaderArray<float> trackEta;
    TTreeReaderArray<float> trackPhi;
    TTreeReaderArray<float> trackTheta;
    TTreeReaderArray<float> trackVarZ0;
    TTreeReaderArray<float> trackVarD0;
    TTreeReaderArray<float> trackVarQp;
    TTreeReaderArray<float> trackVarTheta;
    TTreeReaderArray<float> trackVarPhi;
    TTreeReaderArray<float> trackTime;
    TTreeReaderArray<float> trackTimeRes;
    TTreeReaderArray<int>   trackTimeValid;
    TTreeReaderArray<int>   trackToTruthvtx;
    TTreeReaderArray<int>   trackToParticle;
    TTreeReaderArray<bool>  trackQuality;
    TTreeReaderArray<int>   trackHgtdHits;
    TTreeReaderArray<int>   trackPrimHits;

    TTreeReaderArray<float> truthVtxZ;
    TTreeReaderArray<float> truthVtxTime;
    TTreeReaderArray<bool>  truthVtxIshs;

    TTreeReaderArray<float> recoVtxZ;
    TTreeReaderArray<float> recoVtxTime;
    TTreeReaderArray<float> recoVtxTimeRes;
    TTreeReaderArray<int>   recoVtxValid;

    TTreeReaderArray<float> topoJetPt;
    TTreeReaderArray<float> topoJetEta;
    TTreeReaderArray<float> topoJetPhi;

    TTreeReaderArray<float> truthHSJetPt;
    TTreeReaderArray<float> truthHSJetEta;
    TTreeReaderArray<float> truthHSJetPhi;

    TTreeReaderArray<float> truthITPUJetPt;
    TTreeReaderArray<float> truthITPUJetEta;
    TTreeReaderArray<float> truthITPUJetPhi;

    TTreeReaderArray<float> truthOOTPUJetPt;
    TTreeReaderArray<float> truthOOTPUJetEta;
    TTreeReaderArray<float> truthOOTPUJetPhi;

    TTreeReaderArray<std::vector<int>> topoJetTruthHSIdx;
    TTreeReaderArray<std::vector<int>> topoJetGhostTrackIdx;

    TTreeReaderArray<float> particleT;

  BranchPointerWrapper(TTreeReader& r)
    : reader (r),
      weight (r, "weight"),
      trackD0 (r, "Track_d0"), trackZ0 (r, "Track_z0"), trackPt (r, "Track_pt"),
      trackEta (r, "Track_eta"), trackQP (r, "Track_qOverP"),
      trackTheta (r, "Track_theta"), trackPhi (r, "Track_phi"),
      trackVarD0 (r, "Track_var_d0"), trackVarZ0 (r, "Track_var_z0"),
      trackVarQp (r, "Track_var_qOverP"), trackVarTheta (r, "Track_var_theta"),
      trackVarPhi (r, "Track_var_phi0"),
      trackTime (r, "Track_time"), trackTimeRes (r, "Track_timeRes"),
      trackTimeValid (r, "Track_hasValidTime"), trackQuality (r, "Track_quality"),
      trackToTruthvtx (r, "Track_truthVtx_idx"),
      trackToParticle (r, "Track_truthPart_idx"),
      trackHgtdHits (r, "Track_nHGTDHits"),
      trackPrimHits (r, "Track_nHGTDPrimaryHits"),
      truthVtxZ (r, "TruthVtx_z"), truthVtxTime (r, "TruthVtx_time"),
      truthVtxIshs (r, "TruthVtx_isHS"), recoVtxZ (r, "RecoVtx_z"),
      recoVtxTime (r, "RecoVtx_time"), recoVtxTimeRes (r, "RecoVtx_timeRes"),
      recoVtxValid (r, "RecoVtx_hasValidTime"),
      topoJetPt (r, "AntiKt4EMTopoJets_pt"),
      topoJetEta (r, "AntiKt4EMTopoJets_eta"),
      topoJetPhi (r, "AntiKt4EMTopoJets_phi"),
      truthHSJetPt (r, "TruthHSJet_pt"), truthHSJetEta (r, "TruthHSJet_eta"),
      truthHSJetPhi (r, "TruthHSJet_phi"),
      truthITPUJetPt (r, "TruthITPUJet_pt"), truthITPUJetEta (r, "TruthITPUJet_eta"),
      truthITPUJetPhi (r, "TruthITPUJet_phi"),
      truthOOTPUJetPt (r, "TruthOOTPUJet_pt"), truthOOTPUJetEta (r, "TruthOOTPUJet_eta"),
      truthOOTPUJetPhi (r, "TruthOOTPUJet_phi"),
      topoJetTruthHSIdx (r, "AntiKt4EMTopoJets_truthHSJet_idx"),
      topoJetGhostTrackIdx (r, "AntiKt4EMTopoJets_ghostTrack_idx"),
      particleT (r, "TruthPart_prodVtx_time")
    {}

    // -----------------------------------------------------------------------
    // passBasicCuts
    //   Returns false (skip event) if:
    //     • fewer than MIN_JETS reco jets are present (fast pre-filter before
    //       the full HS-matching check in passJetPtCut), or
    //     • the reconstructed HS vertex z is more than MAX_VTX_DZ from the
    //       truth HS vertex z (guards against wrong-vertex events).
    // -----------------------------------------------------------------------
    bool passBasicCuts() {
      if (this->topoJetPt.GetSize() < MIN_JETS) {
	if (DEBUG) std::cout << "Skipping low jet event\n";
        return false;
      }
    
      // check reco HS vertex is with 2mm of truth HS vertex
      if(std::abs(this->truthVtxZ[0] - this->recoVtxZ[0]) > MAX_VTX_DZ) {
	if(DEBUG) std::cout << "Skipping event due to incorrect HS vertex\n";
	return false;
      }
      
      return true;
    }

    // -----------------------------------------------------------------------
    // calcBestVbsDeltaEta
    //   Among the given pT-passing jet indices, finds the pair in opposite
    //   η hemispheres (η_i · η_j < 0) with the largest invariant mass m_jj
    //   — the standard VBS-candidate-jet definition — and returns their |Δη|.
    //   Returns -1 if no opposite-hemisphere pair exists.
    // -----------------------------------------------------------------------
    double calcBestVbsDeltaEta(const std::vector<int>& passPtIdx) const {
      double bestMjj = -1.0, bestDEta = -1.0;
      for (size_t a = 0; a < passPtIdx.size(); ++a) {
        for (size_t b = a + 1; b < passPtIdx.size(); ++b) {
          int i = passPtIdx[a], j = passPtIdx[b];
          float etaI = this->topoJetEta[i], etaJ = this->topoJetEta[j];
          if (etaI * etaJ >= 0) continue;  // require opposite hemispheres

          TLorentzVector ji, jj;
          ji.SetPtEtaPhiM(this->topoJetPt[i], etaI, this->topoJetPhi[i], 0.0);
          jj.SetPtEtaPhiM(this->topoJetPt[j], etaJ, this->topoJetPhi[j], 0.0);
          double mjj = (ji + jj).M();

          if (mjj > bestMjj) {
            bestMjj  = mjj;
            bestDEta = std::abs(etaI - etaJ);
          }
        }
      }
      return bestDEta;
    }

    // -----------------------------------------------------------------------
    // passJetPtCut
    //   Requires at least MIN_PASSPT_JETS reco jets above MIN_JET_PT,
    //   at least MIN_PASSETA_JETS of which are in the forward HGTD acceptance,
    //   and that the VBS candidate pair — the opposite-hemisphere,
    //   pT-passing jet pair with the largest m_jj — has an η separation ≥
    //   VBS_JET_D_ETA. No truth-HS matching required.
    // -----------------------------------------------------------------------
    bool passJetPtCut() {
      int passptcount = 0, passptetacount = 0;
      std::vector<int> passPtIdx;  // indices of pT-passing jets

      for (int jetIdx = 0; jetIdx < (int)this->topoJetPt.GetSize(); ++jetIdx) {
        float eta = std::abs(this->topoJetEta[jetIdx]);
        float pt  = this->topoJetPt[jetIdx];
        if (DEBUG) std::cout << "reco jet pt, eta: " << pt << ", " << eta << '\n';
        bool passpt    = pt > MIN_JET_PT;
        bool passpteta = passpt && eta > MIN_ABS_ETA_JET && eta < MAX_ABS_ETA_JET;
        passptcount    += passpt    ? 1 : 0;
        passptetacount += passpteta ? 1 : 0;
        if (passpt) passPtIdx.push_back(jetIdx);
      }

      bool passesPt   = passptcount   >= MIN_PASSPT_JETS;
      bool passesEta  = passptetacount >= MIN_PASSETA_JETS;
      double bestDEta = calcBestVbsDeltaEta(passPtIdx);
      bool passesDEta = bestDEta >= VBS_JET_D_ETA;
      return passesPt && passesEta && passesDEta;
    }

    // -----------------------------------------------------------------------
    // passForwardHsTracks  [testing only]
    //   Gates on the number of forward hard-scatter tracks being within the
    //   window [MIN_NHS_TRACK, MAX_NHS_TRACK].  Used during development to
    //   restrict studies to events with a controlled HS track multiplicity;
    //   not applied in the primary analysis.
    // -----------------------------------------------------------------------
    bool passForwardHsTracks(int& nForwardHSTrack) {
      return
	MAX_NHS_TRACK >= nForwardHSTrack and
	nForwardHSTrack >= MIN_NHS_TRACK;
    }

    // -----------------------------------------------------------------------
    // countForwardJets
    //   Counts reco jets that fall in the HGTD η acceptance
    //   (MIN_ABS_ETA_JET < |η| < MAX_ABS_ETA_JET) and exceed MIN_JET_PT.
    //   The result is written into nForwardJet and used as the x-axis
    //   variable in several efficiency/resolution plots.
    // -----------------------------------------------------------------------
    void countForwardJets(int& nForwardJet) {
    nForwardJet = 0;
    for(int jetIdx = 0; jetIdx < this->topoJetEta.GetSize(); ++jetIdx) {
      float
	jetEta = std::abs(this->topoJetEta[jetIdx]),
	jetPt = this->topoJetPt[jetIdx];
      if (jetEta > MIN_ABS_ETA_JET and jetEta < MAX_ABS_ETA_JET and jetPt > MIN_JET_PT)
	nForwardJet++;
    }
  }

    // -----------------------------------------------------------------------
    // countTruthHSJets
    //   Counts truth hard-scatter jets in the WHOLE event (no η restriction)
    //   above MIN_JET_PT.  Used as the binning x-axis variable in place of the
    //   forward reco-jet count.
    // -----------------------------------------------------------------------
    void countTruthHSJets(int& nTruthHSJet) {
    nTruthHSJet = 0;
    for (int jetIdx = 0; jetIdx < (int)this->truthHSJetPt.GetSize(); ++jetIdx) {
      if (this->truthHSJetPt[jetIdx] > MIN_JET_PT)
	nTruthHSJet++;
    }
  }

    // -----------------------------------------------------------------------
    // countForwardTracks
    //   Splits a pre-selected track list into forward HS and PU components.
    //   A track is counted if it passes the HGTD η window, pT bounds, track
    //   quality flag, and (when checkTimeValid is true) has a valid HGTD time.
    //   Writes totals into nFTrack, nFTrackHS, nFTrackPU; used as x-axis
    //   variables in efficiency/resolution plots.
    // -----------------------------------------------------------------------
    void countForwardTracks(
      int& nFTrack, int& nFTrackHS, int& nFTrackPU,
      const std::vector<int>& tracks, bool checkTimeValid
    ) {
      nFTrack = 0; nFTrackHS = 0; nFTrackPU = 0;
      for(auto trkIdx: tracks) {
	double eta = std::abs(this->trackEta[trkIdx]);
	double pt = this->trackPt[trkIdx];
	bool quality = this->trackQuality[trkIdx] == true;
	bool hasValidTime = checkTimeValid ? (this->trackTimeValid[trkIdx] == 1) : true;
	// already know these pass association
	if (eta > MIN_ABS_ETA_TRACK and eta < MAX_ABS_ETA_TRACK and
	    pt > MIN_TRACK_PT_COUNT and pt < MAX_TRACK_PT and
	    quality and hasValidTime) {
	  nFTrack++;
	  // Fetch truth vertex index only for tracks that pass the gate
	  int truthVtx = this->trackToTruthvtx[trkIdx];
	  if (truthVtx != -1 and this->truthVtxIshs[truthVtx])
	    nFTrackHS++;
	  else
	    nFTrackPU++;
	}
      }
    }

    // -----------------------------------------------------------------------
    // calcJetptDRScore
    //   Returns jet_pT * exp(−ΔR) where ΔR is the angular distance to the
    //   nearest reco jet.  Higher values indicate the track is close to a
    //   high-pT jet; used as a per-track weight in FILTJET-style scoring.
    // -----------------------------------------------------------------------
    double calcJetptDRScore(int trkIdx) {
      double
	trkEta = this->trackEta[trkIdx],
	trkPhi = this->trackPhi[trkIdx];

      double minDR = 1e6;
      int minIdx = -1;
      for (int jetIdx=0; jetIdx < this->topoJetEta.GetSize(); ++jetIdx) {
	double
	  jetEta = this->topoJetEta[jetIdx],
	  jetPhi = this->topoJetPhi[jetIdx];
	double
	  deta = jetEta-trkEta,
	  dphi = TVector2::Phi_mpi_pi(jetPhi - trkPhi);
	double thisDR = std::hypot(deta, dphi);
	if (thisDR < minDR) {
	  minDR = thisDR;
	  minIdx = jetIdx;
	}
      }
      double returnScore = this->topoJetPt[minIdx]*std::exp(-minDR);
      return returnScore;
    }

    // -----------------------------------------------------------------------
    // calcTrkptDRScore
    //   Returns track_pT * exp(−ΔR) where ΔR is the angular distance to the
    //   nearest reco jet.  Analogous to calcJetptDRScore but uses the track's
    //   own pT rather than the jet pT; used as the per-track contribution to
    //   the TRKPTZ cluster score.
    // -----------------------------------------------------------------------
    double calcTrkptDRScore(int trkIdx) {
      double
	trkEta = this->trackEta[trkIdx],
	trkPhi = this->trackPhi[trkIdx];


      double minDR = 1e6;
      for (int jetIdx = 0; this->topoJetEta.GetSize() > jetIdx; ++jetIdx) {
        double
	  jetEta = this->topoJetEta[jetIdx],
	  jetPhi = this->topoJetPhi[jetIdx];

	double
	  deta = jetEta-trkEta,
	  dphi = TVector2::Phi_mpi_pi(jetPhi - trkPhi);

	double thisDR = std::sqrt(deta*deta + dphi*dphi);
	if (thisDR < minDR) {
	  minDR = thisDR;
	}
      }
      double returnScore = this->trackPt[trkIdx]*std::exp(-minDR);
      return returnScore;
    }

  };

  // ---------------------------------------------------------------------------
  // EventCounts
  //   Holds all per-event counts and pre-folded histogram fill values
  //   derived from a 3σ track selection.  Computed once at the top of
  //   processEventData and passed to the histogram fill helpers below.
  // ---------------------------------------------------------------------------
  struct EventCounts {
    // Raw forward track/jet counts
    int nForwardJet      = 0;
    int nTruthHSJet      = 0;   // truth HS jets in the whole event (> MIN_JET_PT)
    int nForwardTrack    = 0;
    int nForwardTrackHS  = 0;
    int nForwardTrackPU  = 0;
    double puRatio       = 0.0;
    // Pre-folded x-values for efficiency / purity histograms
    int    effFillValFjet;
    int    effFillValHSJet;
    int    effFillValTrack;
    int    effFillValHSTrack;
    int    effFillValPUTrack;
    double effFillValPURatio;

    EventCounts(BranchPointerWrapper* branch,
                const std::vector<int>& tracks,
                bool checkValidTimes) {
      branch->countForwardJets(nForwardJet);
      branch->countTruthHSJets(nTruthHSJet);
      branch->countForwardTracks(nForwardTrack, nForwardTrackHS, nForwardTrackPU,
                                 tracks, checkValidTimes);
      puRatio          = (double)nForwardTrackPU / (double)nForwardTrack;
      effFillValFjet    = folded(nForwardJet,     (int)FOLD_FJET);
      effFillValHSJet   = folded(nTruthHSJet,     (int)FOLD_FJET);
      effFillValTrack   = folded(nForwardTrack,   (int)FOLD_TRACK);
      effFillValHSTrack = folded(nForwardTrackHS, (int)FOLD_HS_TRACK);
      effFillValPUTrack = folded(nForwardTrackPU, (int)FOLD_PU_TRACK);
      effFillValPURatio = folded(puRatio,         FOLD_PU_FRAC);
    }
  };

  // ---------------------------------------------------------------------------
  // Cluster
  //   Represents a single time(-z₀) cluster produced by doSimultaneousClustering
  //   or doConeClustering.  The weighted-mean position is stored in values[],
  //   with per-dimension uncertainties in sigmas[].  allTimes holds the raw
  //   track times for spread diagnostics.  trackIndices lists the original
  //   track indices from the branch arrays so that truth-matching and feature
  //   extraction can look up any per-track quantity.
  //
  //   Member methods:
  //     operator== / !=  — equality by (values[0], sigmas[0], nConstituents)
  //     calcPurity        — fraction of cluster pT belonging to truth HS vertex
  //     updateScores      — computes TRKPTZ, WAVES, JET_T_REFINED, WAVES_MISCL/MISAS, TEST_MISAS scores
  //     timeSpread        — std-dev of raw track times within the cluster
  //     zSpread           — std-dev of track z₀ values within the cluster
  //     passEfficiency    — true if cluster time is within 3·PASS_SIGMA of truth
  //     calcFeatures      — extracts and normalises the 10 ML input features
  // ---------------------------------------------------------------------------
  struct Cluster {
    std::vector<double> values;
    std::vector<double> sigmas;
    std::vector<double> allTimes;
    std::vector<int> trackIndices;
    std::unordered_map<int,double> scores;
    double purity = 0.0;
    bool wasMerged = false;
    int nConstituents=1;

    // -----------------------------------------------------------------------
    // operator== / operator!=
    //   Two clusters are considered identical if they share the same
    //   weighted-mean time (values[0]), time uncertainty (sigmas[0]), and
    //   constituent count.  This is sufficient for the duplicate-removal
    //   checks inside chooseCluster.
    // -----------------------------------------------------------------------
    bool operator==(const Cluster& other) {
      bool sameValues = values.at(0) == other.values.at(0);
      bool sameSigmas = sigmas.at(0) == other.sigmas.at(0);
      bool sameConsts = nConstituents == other.nConstituents;
      // this SHOULD be sufficient
      return sameConsts and sameSigmas and sameValues;
    }

    bool operator!=(const Cluster& other) {
      return !(*this == other);
    }

    // -----------------------------------------------------------------------
    // calcPurity
    //   Computes the pT-weighted fraction of cluster tracks that originate
    //   from the hard-scatter vertex (truth vertex index 0).  Stores the
    //   result in this->purity.  Called only when a purity-gated score is
    //   active (calcPurityFlag == true) to avoid unnecessary truth-matching work.
    // -----------------------------------------------------------------------
    void calcPurity(BranchPointerWrapper *branch) {
      double num = 0.0, denom = 0.0;
      for (auto trk: this->trackIndices) {
	if (branch->trackToTruthvtx[trk] == 0) num += branch->trackPt[trk];
	denom += branch->trackPt[trk];
      }
      this->purity = num/denom;
    }

    // -----------------------------------------------------------------------
    // updateScores
    //   Populates the score map for all derived scoring algorithms after
    //   clustering is complete.  Called once per cluster in
    //   clusterTracksInTime.
    //
    //   TRKPTZ   — TRKPT weighted by exp(−1.5 · |Δz|), where Δz is the
    //              distance from the cluster z₀ to the reco primary vertex.
    //              When z₀ is available as values[1] it is used directly;
    //              otherwise a precision-weighted average over track z₀ is
    //              computed on the fly.
    // -----------------------------------------------------------------------
    void updateScores(
      BranchPointerWrapper *branch
    ) {
      // Call calcFeatures once — it returns the normalised feature vector and
      // the raw deltaZ (cluster z − reco vertex z).  Reusing the returned
      // deltaZ avoids a second precision-weighted z-average pass over tracks.
      auto [features, rawDeltaZ, rawDeltaZResunits] = this->calcFeatures(branch);

      // TRKPTZ: score = TRKPT × exp(−1.5 · |Δz|)
      if (this->values.size() > 1) {
        // Clustering was done with z₀ as a second dimension: use that value.
        double dz = std::abs(this->values.at(1) - branch->recoVtxZ[0]);
        this->scores[Score::TRKPTZ.id] =
          this->scores.at(Score::TRKPT.id) * std::exp(-1.5 * dz);
      } else {
        // Common case (usez0=false): reuse deltaZ already computed by calcFeatures.
        this->scores[Score::TRKPTZ.id] =
          this->scores.at(Score::TRKPT.id) * std::exp(-1.5 * std::abs(rawDeltaZ));
      }

      // WAVES: WAVeS-style score — Σ_i pT_i × pT_jet(i) / max(ΔR_i, DR_FLOOR)
      // multiplied by exp(−1.5|Δz_cluster|), where Δz is the pT-weighted cluster z centroid
      // minus the reco vertex z.  The cluster-level z-term is more effective than per-track
      // z-pull weighting because it averages over track-by-track z noise.
      // Linear pT (not squared) so a couple of high-pT time-misassigned tracks can't
      // outvote a larger time-coherent cluster.
      // Falls back to TRKPTZ if no qualifying forward jets exist.
      {
        double wavesSum = 0.0;
        for (int idx : this->trackIndices) {
          double trkEta  = branch->trackEta[idx];
          double trkPhi  = branch->trackPhi[idx];
          double trkPt   = branch->trackPt[idx];

          double minDR     = 1e6;
          double nearJetPt = 0.0;
          for (int j = 0; j < (int)branch->topoJetEta.GetSize(); ++j) {
            double jEta = branch->topoJetEta[j];
            double jPt  = branch->topoJetPt[j];
            if (jPt < MIN_JET_PT) continue;
            if (std::abs(jEta) < MIN_ABS_ETA_JET || std::abs(jEta) > MAX_ABS_ETA_JET) continue;
            double deta = jEta - trkEta;
            double dphi = TVector2::Phi_mpi_pi(branch->topoJetPhi[j] - trkPhi);
            double dr   = std::hypot(deta, dphi);
            if (dr < minDR) { minDR = dr; nearJetPt = jPt; }
          }
          if (nearJetPt <= 0.0) continue;
          wavesSum += trkPt * nearJetPt
                      / std::max(minDR, WAVES_DR_FLOOR);
        }
        if (wavesSum > 0.0)
          this->scores[Score::WAVES.id] =
              wavesSum * std::exp(-1.5 * std::abs(rawDeltaZ));
        else
          this->scores[Score::WAVES.id] = this->scores.at(Score::TRKPTZ.id);
      }

      // JET_T_REFINED: dedicated collection (jet-filtered tracks at 2σ iterative);
      // cluster selected by TRKPTZ via the aux-collection path in selectClusters().
      this->scores[Score::JET_T_REFINED.id] = this->scores.at(Score::TRKPTZ.id);

      // WAVeS oracle variants: selected by the WAVeS score; denominator gates
      // (cluster purity / HS timing purity) applied at fill time in event_processing.h.
      this->scores[Score::WAVES_MISCL.id] = this->scores.at(Score::WAVES.id);
      this->scores[Score::WAVES_MISAS.id] = this->scores.at(Score::WAVES.id);

      // TEST_MISAS uses TRKPTZ as its selection score; the purity gate is applied
      // at efficiency-check time in event_processing.h (both pass and total fills).
      this->scores[Score::TEST_MISAS.id] = this->scores.at(Score::TRKPTZ.id);
    }
    
    // -----------------------------------------------------------------------
    // timeSpread
    //   Returns the population standard deviation of the raw track times
    //   stored in allTimes.  Used as a diagnostic variable in error-analysis
    //   plots to characterise how well the cluster's tracks agree in time.
    // -----------------------------------------------------------------------
    double timeSpread() {
      // calculate stdev of times
      auto avg = std::accumulate(allTimes.begin(),allTimes.end(),0.0)/this->allTimes.size();
      auto ssd = 0.0;
      for (auto t: allTimes) {
	ssd += (avg - t)*(avg - t);
      }

      return std::sqrt(ssd/this->allTimes.size());
    }

    // -----------------------------------------------------------------------
    // zSpread
    //   Returns the population standard deviation of track z₀ values for
    //   the cluster constituents.  Analogous to timeSpread but in the
    //   longitudinal direction; used in error-analysis diagnostics.
    // -----------------------------------------------------------------------
    double zSpread(BranchPointerWrapper *bpw) {
      // calculate stdev of z0 values
      auto avg = 0.0;
      for (auto trk: trackIndices)
	avg += bpw->trackZ0[trk];
      avg *= 1./this->trackIndices.size();
      auto ssd = 0.0;
      for (auto t: trackIndices) {
	ssd += (avg - bpw->trackZ0[t])*(avg - bpw->trackZ0[t]);
      }
      return std::sqrt(ssd/this->trackIndices.size());
    }

    // -----------------------------------------------------------------------
    // passEfficiency
    //   Returns true if the cluster's weighted-mean time is within
    //   PASS_SIGMA of the truth hard-scatter vertex time.
    //   This is the primary efficiency criterion used when filling the
    //   pass/total histograms in event_processing.h.
    // -----------------------------------------------------------------------
    bool passEfficiency(BranchPointerWrapper *branch) const {
      if (DEBUG) std::cout << "Choosing pass score\n";
      if (this->values.size() == 0)
	return false;

      double diff = std::abs(this->values.at(0)-branch->truthVtxTime[0]);
      if (diff > PASS_SIGMA)
	return false;

      return true;
    }

    // -----------------------------------------------------------------------
    // calculateTime
    //   Returns the time to use for diff/resolution plots and passEfficiency
    //   for the given score.  Implemented out-of-line in event_processing.h
    //   (after passTrackVertexAssociation is defined) to avoid circular includes.
    // -----------------------------------------------------------------------
    double calculateTime(
        Score score,
	BranchPointerWrapper* branch,
	double idealRes = -1.0
    ) const;

    // -----------------------------------------------------------------------
    // calculatePurity
    //   Returns the purity to use for resolution and efficiency plots for the
    //   given score.  For most scores this is the pre-computed this->purity
    //   (fraction of cluster ΣpT from HS tracks).  For WAVES, purity
    //   is re-evaluated using only the constituent tracks that survive the
    //   dR < 0.4 HS-jet filter, mirroring the time computation in calculateTime.
    //   Implemented out-of-line in event_processing.h alongside calculateTime.
    // -----------------------------------------------------------------------
    double calculatePurity(
        Score score,
        BranchPointerWrapper* branch
    ) const;

    // -----------------------------------------------------------------------
    // calcFeatures
    //   Extracts the 8 features used as DNN input and applies stored
    //   normalisation (means/stds from normalization_params.json) to return a
    //   unit-normalised feature vector.  Features (in training order):
    //     0  delta_z              — cluster z relative to reco primary vertex
    //     1  delta_z_resunits     — delta_z in units of cluster_z_sigma
    //     2  cluster_z_sigma      — precision-weighted z₀ uncertainty
    //     3  cluster_d0           — precision-weighted transverse impact parameter
    //     4  cluster_d0_sigma     — d0 uncertainty
    //     5  cluster_qOverP       — precision-weighted charge-over-momentum
    //     6  cluster_qOverP_sigma — q/p uncertainty
    //     7  cluster_sumpt        — scalar sum of constituent track pT
    //
    //   Returns {normalised feature vector, raw deltaZ}.
    //   The raw deltaZ is returned alongside so updateScores can compute
    //   TRKPTZ without repeating the z weighted-average loop.
    // -----------------------------------------------------------------------
    std::tuple<std::vector<float>, float, float> calcFeatures(BranchPointerWrapper *branch) {
      float znum  = 0.f, zden  = 0.f;
      float dnum  = 0.f, dden  = 0.f;
      float qpnum = 0.f, qpden = 0.f;
      float sumpt = 0.f;

      for (auto trk: trackIndices) {
        float trkVarZ = branch->trackVarZ0[trk];
        float trkVarD = branch->trackVarD0[trk];
        float trkVarQ = branch->trackVarQp[trk];

        znum  += branch->trackZ0[trk] / trkVarZ;  zden  += 1.f / trkVarZ;
        dnum  += branch->trackD0[trk] / trkVarD;  dden  += 1.f / trkVarD;
        qpnum += branch->trackQP[trk] / trkVarQ;  qpden += 1.f / trkVarQ;
        sumpt += branch->trackPt[trk];
      }

      float clusterZ           = znum / zden;
      float clusterZSigma      = 1.0f / std::sqrt(zden);
      float clusterD0          = dnum / dden;
      float clusterD0Sigma     = 1.0f / std::sqrt(dden);
      float clusterQOverP      = qpnum / qpden;
      float clusterQOverPSigma = 1.0f / std::sqrt(qpden);

      float deltaZ         = clusterZ - branch->recoVtxZ[0];
      float deltaZResunits = deltaZ / clusterZSigma;

      // Normalization parameters from share/models/normalization_params.json
      static const float MEANS[8] = {
        0.005349026450717473,  // delta_z mean
        0.006244652963375047,  // delta_z_resunits mean
        0.47186453018432917,  // cluster_z_sigma mean
        -0.00025930208416971584,  // cluster_d0 mean
        0.041210234810337115,  // cluster_d0_sigma mean
        2.497333134513254e-07,  // cluster_qOverP mean
        2.399697181725376e-06,  // cluster_qOverP_sigma mean
        22.605882118982347  // cluster_sumpt mean
      };
      static const float STDS[8] = {
        1.2256897153543274,  // delta_z std
        1.8559666108975377,  // delta_z_resunits std
        0.3927622369867232,  // cluster_z_sigma std
        0.15349002263169267,  // cluster_d0 std
        0.025843144355708243,  // cluster_d0_sigma std
        2.2377624188670666e-05,  // cluster_qOverP std
        1.5046303393371763e-06,  // cluster_qOverP_sigma std
        23.948572592936106  // cluster_sumpt std
      };

      std::vector<float> features = {
        (deltaZ             - MEANS[0]) / STDS[0],
        (deltaZResunits     - MEANS[1]) / STDS[1],
        (clusterZSigma      - MEANS[2]) / STDS[2],
        (clusterD0          - MEANS[3]) / STDS[3],
        (clusterD0Sigma     - MEANS[4]) / STDS[4],
        (clusterQOverP      - MEANS[5]) / STDS[5],
        (clusterQOverPSigma - MEANS[6]) / STDS[6],
        (sumpt              - MEANS[7]) / STDS[7],
      };

      return {features, deltaZ, deltaZResunits};
    }
  };
  
}
#endif // STRUCTS_H
