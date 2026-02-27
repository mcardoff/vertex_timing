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
#include "ml_model.h"
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
    TTreeReaderArray<float> trackTime;
    TTreeReaderArray<float> trackTimeRes;
    TTreeReaderArray<int>   trackTimeValid;
    TTreeReaderArray<int>   trackToTruthvtx;
    TTreeReaderArray<int>   trackToParticle;
    TTreeReaderArray<bool>  trackQuality;
    TTreeReaderArray<int>   trackHgtdHits;
    TTreeReaderArray<int>   trackPrimHits;
    TTreeReaderArray<float> trackNearIdx;
    TTreeReaderArray<float> trackNearZ0sin;
    TTreeReaderArray<float> trackNearZ0sinUnc;
    TTreeReaderArray<float> trackNearSig;

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

    TTreeReaderArray<float> particleT;

  BranchPointerWrapper(TTreeReader& r)
    : reader (r),
      weight (r, "weight"),
      trackD0 (r, "Track_d0"),
      trackZ0 (r, "Track_z0"),
      trackPt (r, "Track_pt"),
      trackEta (r, "Track_eta"),
      trackQP (r, "Track_qOverP"),
      trackTheta (r, "Track_theta"),
      trackPhi (r, "Track_phi"),
      trackVarD0 (r, "Track_var_d0"),
      trackVarZ0 (r, "Track_var_z0"),
      trackVarQp (r, "Track_var_qOverP"),
      trackVarTheta (r, "Track_var_theta"),
      trackTime (r, "Track_time"),
      trackTimeRes (r, "Track_timeRes"),
      trackTimeValid (r, "Track_hasValidTime"),
      trackQuality (r, "Track_quality"),
      trackToTruthvtx (r, "Track_truthVtx_idx"),
      trackToParticle (r, "Track_truthPart_idx"),
      trackHgtdHits (r, "Track_nHGTDHits"),
      trackPrimHits (r, "Track_nHGTDPrimaryHits"),
      trackNearIdx (r, "Track_nearestVtx_idx"),
      trackNearZ0sin (r, "Track_nearestVtx_z0SinTheta"),
      trackNearZ0sinUnc (r, "Track_nearestVtx_z0SinThetaUncertainty"),
      trackNearSig (r, "Track_nearestVtx_sig"),
      truthVtxZ (r, "TruthVtx_z"),
      truthVtxTime (r, "TruthVtx_time"),
      truthVtxIshs (r, "TruthVtx_isHS"),
      recoVtxZ (r, "RecoVtx_z"),
      recoVtxTime (r, "RecoVtx_time"),
      recoVtxTimeRes (r, "RecoVtx_timeRes"),
      recoVtxValid (r, "RecoVtx_hasValidTime"),
      topoJetPt (r, "AntiKt4EMTopoJets_pt"),
      topoJetEta (r, "AntiKt4EMTopoJets_eta"),
      topoJetPhi (r, "AntiKt4EMTopoJets_phi"),
      truthHSJetPt (r, "TruthHSJet_pt"),
      truthHSJetEta (r, "TruthHSJet_eta"),
      particleT (r, "TruthPart_prodVtx_time")
    {}

    // -----------------------------------------------------------------------
    // passBasicCuts
    //   Returns false (skip event) if:
    //     • fewer than MIN_JETS truth-HS jets or reco jets are present, or
    //     • the reconstructed HS vertex z is more than MAX_VTX_DZ from the
    //       truth HS vertex z (guards against wrong-vertex events).
    // -----------------------------------------------------------------------
    bool passBasicCuts() {
      if (this->truthHSJetPt.GetSize() < MIN_JETS or
	  this->topoJetPt.GetSize() < MIN_JETS) {
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
    // passJetPtCut
    //   Requires the truth-HS jet collection to contain at least
    //   MIN_PASSPT_JETS jets above MIN_JETPT, at least MIN_PASSETA_JETS of
    //   which are in the forward HGTD acceptance, and that the two leading
    //   jets have an η separation ≥ VBS_JET_D_ETA (VBS topology selection).
    // -----------------------------------------------------------------------
    bool passJetPtCut() {
      int passptcount = 0, passptetacount = 0;
      for(int jetIdx = 0; jetIdx < this->truthHSJetEta.GetSize(); ++jetIdx) {
	float
	  eta = std::abs(truthHSJetEta[jetIdx]),
	  pt = truthHSJetPt[jetIdx];
	if (DEBUG) std::cout << "pt, eta: " << pt << ", " << eta << '\n';
	bool passpt = pt > MIN_JETPT;
	bool passpteta = passpt and eta > MIN_ABS_ETA_JET and eta < MAX_ABS_ETA_JET;
	passptcount += passpt ? 1 : 0;
	passptetacount += passpteta ? 1 : 0;
      }
      
      bool passesPt   = passptcount >= MIN_PASSPT_JETS;
      bool passesEta  = passptetacount >= MIN_PASSETA_JETS;
      bool passesDEta = std::abs(truthHSJetEta[0]-truthHSJetEta[1]) >= VBS_JET_D_ETA;
      return passesPt and passesEta and passesDEta;
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
    //   (MIN_ABS_ETA_JET < |η| < MAX_ABS_ETA_JET) and exceed MIN_JETPT.
    //   The result is written into nForwardJet and used as the x-axis
    //   variable in several efficiency/resolution plots.
    // -----------------------------------------------------------------------
    void countForwardJets(int& nForwardJet) {
    nForwardJet = 0;
    for(int jetIdx = 0; jetIdx < this->topoJetEta.GetSize(); ++jetIdx) {
      float
	jetEta = std::abs(this->topoJetEta[jetIdx]),
	jetPt = this->topoJetPt[jetIdx];
      if (jetEta > MIN_ABS_ETA_JET and jetEta < MAX_ABS_ETA_JET and jetPt > MIN_JETPT)
	nForwardJet++;
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
	int  truthVtx = this->trackToTruthvtx[trkIdx];
	bool quality = this->trackQuality[trkIdx] == true;
	bool hasValidTime = checkTimeValid ? (this->trackTimeValid[trkIdx] == 1) : true;
	// already know these pass association
	if (eta > MIN_ABS_ETA_TRACK and eta < MAX_ABS_ETA_TRACK and
	    pt > MIN_TRACK_PT_COUNT and pt < MAX_TRACK_PT and
	    quality and hasValidTime) {
	  nFTrack++;
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

    // -----------------------------------------------------------------------
    // VBF H→Invisible signal region selection
    //   Implements the kinematic cuts from JHEP08(2022)104.
    //   Jets are assumed to be pT-sorted in descending order (leading jet
    //   at index 0, sub-leading at index 1) — verified to be 100% true in
    //   this sample via check_jet_sorting.cxx.
    //   Some thresholds are loosened relative to the published analysis to
    //   retain more events for timing studies; see individual comments.
    //
    //   Individual predicates (each returns bool unless noted):
    //     countJetsAbovePt    — count jets above a pT threshold
    //     passJetMultiplicity — ≥ 2 jets with pT > 25 GeV
    //     passLeadingJetPt    — leading > 60 GeV, sub-leading > 40 GeV
    //     passJetDeltaPhi     — |Δφ_jj| < 2 (not back-to-back)
    //     passOppositeHemispheres — η_j1 · η_j2 < 0
    //     passLargeDeltaEta   — |Δη_jj| > 3.0
    //     calcDijetMass       — computes m_jj (GeV); used by passLargeDijetMass
    //     passLargeDijetMass  — m_jj > 300 GeV
    //     calcAllJetsPtSum    — scalar sum of all jet pT (GeV)
    //     passAllJetsPtSum    — scalar sum > 140 GeV
    //     passForwardJet      — ≥ 1 leading jet in HGTD η acceptance
    //     calcCentrality      — centrality C_i for jet i relative to j1/j2
    //     calcRelativeMass    — relative mass m_rel_i = min(m_j1i,m_j2i)/m_jj
    //     passFSRCompatibility — C_i < 0.6 and m_rel < 0.05 for j3, j4
    //     passVBFSignalRegion — AND of all above (FSR cut currently disabled)
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // countJetsAbovePt
    //   Returns the number of reco jets with pT above ptThreshold (GeV).
    //   Helper for passJetMultiplicity.
    // -----------------------------------------------------------------------
    int countJetsAbovePt(double ptThreshold = 25.0) const {
      int count = 0;
      for (int jetIdx = 0; jetIdx < this->topoJetPt.GetSize(); ++jetIdx) {
        if (this->topoJetPt[jetIdx] > ptThreshold) {
          count++;
        }
      }
      return count;
    }

    // -----------------------------------------------------------------------
    // passJetMultiplicity
    //   Requires at least 2 reco jets with pT > 25 GeV.
    // -----------------------------------------------------------------------
    bool passJetMultiplicity() const {
      int nJets = countJetsAbovePt(25.0);
      return (nJets >= 2);
    }

    // -----------------------------------------------------------------------
    // passLeadingJetPt
    //   Requires leading jet pT > 60 GeV and sub-leading jet pT > 40 GeV.
    //   Original ATLAS cuts (JHEP08(2022)104): leading > 80 GeV,
    //   sub-leading > 50 GeV — loosened here to retain more events for
    //   timing studies.
    // -----------------------------------------------------------------------
    bool passLeadingJetPt() const {
      if (this->topoJetPt.GetSize() < 2) return false;

      // Jets are assumed to be pT-sorted
      double leadingPt = this->topoJetPt[0];
      double subleadingPt = this->topoJetPt[1];

      return (leadingPt > 60.0 && subleadingPt > 40.0);
    }

    // -----------------------------------------------------------------------
    // passJetDeltaPhi
    //   Requires |Δφ_jj| < 2 between the two leading jets, rejecting
    //   back-to-back topologies inconsistent with VBF kinematics.
    // -----------------------------------------------------------------------
    bool passJetDeltaPhi() const {
      if (this->topoJetPt.GetSize() < 2) return false;

      double phi1 = this->topoJetPhi[0];
      double phi2 = this->topoJetPhi[1];
      double deltaPhi = std::abs(TVector2::Phi_mpi_pi(phi1 - phi2));

      return (deltaPhi < 2.0);
    }

    // -----------------------------------------------------------------------
    // passOppositeHemispheres
    //   Requires the two leading jets to be in opposite η hemispheres
    //   (η_j1 · η_j2 < 0), a hallmark of VBF topology.
    // -----------------------------------------------------------------------
    bool passOppositeHemispheres() const {
      if (this->topoJetEta.GetSize() < 2) return false;

      double eta1 = this->topoJetEta[0];
      double eta2 = this->topoJetEta[1];

      return (eta1 * eta2 < 0);
    }

    // -----------------------------------------------------------------------
    // passLargeDeltaEta
    //   Requires |Δη_jj| > minDeltaEta (default 3.0) between the two
    //   leading jets.  Original ATLAS cut is Δη_jj > 3.8; loosened here
    //   for timing studies.
    // -----------------------------------------------------------------------
    bool passLargeDeltaEta(double minDeltaEta = 3.0) const {
      if (this->topoJetEta.GetSize() < 2) return false;

      double eta1 = this->topoJetEta[0];
      double eta2 = this->topoJetEta[1];
      double deltaEta = std::abs(eta1 - eta2);

      return (deltaEta > minDeltaEta);
    }

    // -----------------------------------------------------------------------
    // calcDijetMass
    //   Builds TLorentzVectors for the two leading jets (assumed massless)
    //   and returns their invariant mass in GeV.  Returns -1 if fewer than
    //   2 jets are present.
    // -----------------------------------------------------------------------
    double calcDijetMass() const {
      if (this->topoJetPt.GetSize() < 2) return -1.0;

      // Build TLorentzVectors for the two leading jets
      TLorentzVector jet1, jet2;

      // Assume massless jets for simplicity
      jet1.SetPtEtaPhiM(this->topoJetPt[0], this->topoJetEta[0],
                        this->topoJetPhi[0], 0.0);
      jet2.SetPtEtaPhiM(this->topoJetPt[1], this->topoJetEta[1],
                        this->topoJetPhi[1], 0.0);

      TLorentzVector dijet = jet1 + jet2;
      return dijet.M(); // Return in GeV (jets already in GeV)
    }

    // -----------------------------------------------------------------------
    // passLargeDijetMass
    //   Requires dijet invariant mass m_jj > minMass GeV (default 300 GeV,
    //   i.e. 0.3 TeV).  Original ATLAS cut is m_jj > 0.8 TeV; loosened
    //   here for timing studies.
    // -----------------------------------------------------------------------
    bool passLargeDijetMass(double minMass = 300.0) const {
      double mjj = calcDijetMass();
      return (mjj > minMass);
    }

    // -----------------------------------------------------------------------
    // calcAllJetsPtSum
    //   Returns the scalar sum of pT over all reco jets (GeV).
    //   Used by passAllJetsPtSum to apply the H_T-like requirement.
    // -----------------------------------------------------------------------
    double calcAllJetsPtSum() const {
      double ptSum = 0.0;
      for (int jetIdx = 0; jetIdx < this->topoJetPt.GetSize(); ++jetIdx) {
        ptSum += this->topoJetPt[jetIdx];
      }
      return ptSum;
    }

    // -----------------------------------------------------------------------
    // passAllJetsPtSum
    //   Requires the scalar sum of all jet pT to exceed minPtSum GeV
    //   (default 140 GeV), applying the H_T-like requirement.
    // -----------------------------------------------------------------------
    bool passAllJetsPtSum(double minPtSum = 140.0) const {
      double ptSum = calcAllJetsPtSum();
      return (ptSum > minPtSum);
    }

    // -----------------------------------------------------------------------
    // passForwardJet
    //   Requires at least one of the two leading jets to fall within the
    //   HGTD η acceptance (MIN_HGTD_ETA < |η| < MAX_HGTD_ETA), ensuring
    //   the event can be studied with HGTD timing information.
    // -----------------------------------------------------------------------
    bool passForwardJet() const {
      if (this->topoJetEta.GetSize() < 2) return false;

      double eta1 = std::abs(this->topoJetEta[0]);
      double eta2 = std::abs(this->topoJetEta[1]);

      bool jet1Forward = (eta1 > MIN_HGTD_ETA && eta1 < MAX_HGTD_ETA);
      bool jet2Forward = (eta2 > MIN_HGTD_ETA && eta2 < MAX_HGTD_ETA);

      return (jet1Forward || jet2Forward);
    }

    // -----------------------------------------------------------------------
    // calcCentrality
    //   Computes the rapidity centrality C_i for jet at jetIdx relative to
    //   the two leading jets:
    //     C_i = exp(−4 · (η_i − η_centre)² / Δη_jj²)
    //   where η_centre = (η_j1 + η_j2) / 2.
    //   C_i → 1 for a jet sitting between j1 and j2; C_i → 0 for a jet
    //   outside the rapidity gap.  Used to identify FSR jets in
    //   passFSRCompatibility.
    // -----------------------------------------------------------------------
    double calcCentrality(int jetIdx) const {
      if (this->topoJetEta.GetSize() < 2) return -1.0;
      if (jetIdx >= this->topoJetEta.GetSize()) return -1.0;

      double etaJ1 = this->topoJetEta[0];
      double etaJ2 = this->topoJetEta[1];
      double etaI = this->topoJetEta[jetIdx];

      // Central position between the two leading jets
      double etaCenter = (etaJ1 + etaJ2) / 2.0;
      double deltaEtaJj = etaJ1 - etaJ2;

      // Calculate centrality
      double numerator = etaI - etaCenter;
      double cI = std::exp(-4.0 / (deltaEtaJj * deltaEtaJj) * numerator * numerator);

      return cI;
    }

    // -----------------------------------------------------------------------
    // calcRelativeMass
    //   Returns the relative invariant mass for jet at jetIdx:
    //     m_rel_i = min(m_j1i, m_j2i) / m_jj
    //   where m_j1i (m_j2i) is the invariant mass of jet i with the leading
    //   (sub-leading) VBF jet.  Small m_rel indicates the extra jet is
    //   consistent with FSR off one of the VBF jets.
    // -----------------------------------------------------------------------
    double calcRelativeMass(int jetIdx) const {
      if (this->topoJetPt.GetSize() < 2) return -1.0;
      if (jetIdx >= this->topoJetPt.GetSize()) return -1.0;

      // Build TLorentzVectors (assume massless jets)
      TLorentzVector jet1, jet2, jetI;
      jet1.SetPtEtaPhiM(this->topoJetPt[0], this->topoJetEta[0], this->topoJetPhi[0], 0.0);
      jet2.SetPtEtaPhiM(this->topoJetPt[1], this->topoJetEta[1], this->topoJetPhi[1], 0.0);
      jetI.SetPtEtaPhiM(this->topoJetPt[jetIdx], this->topoJetEta[jetIdx], this->topoJetPhi[jetIdx], 0.0);

      // Calculate invariant masses
      double mJj = (jet1 + jet2).M();
      double mJ1i = (jet1 + jetI).M();
      double mJ2i = (jet2 + jetI).M();

      // Return relative mass
      double mRel = std::min(mJ1i, mJ2i) / mJj;
      return mRel;
    }

    // -----------------------------------------------------------------------
    // passFSRCompatibility
    //   For the third and fourth reco jets (if present), checks that each
    //   satisfies C_i < 0.6 and m_rel_i < 0.05, consistent with being FSR
    //   rather than additional hard-scatter activity.  Events with only 2
    //   jets pass automatically.
    //   Note: this cut is currently commented out in passVBFSignalRegion.
    // -----------------------------------------------------------------------
    bool passFSRCompatibility() const {
      int nJets = this->topoJetPt.GetSize();

      // If only 2 jets, automatically pass (no FSR jets to check)
      if (nJets <= 2) return true;

      // Check third jet (index 2)
      if (nJets >= 3) {
        double cJ3 = calcCentrality(2);
        double mRelJ3 = calcRelativeMass(2);
        if (cJ3 >= 0.6 || mRelJ3 >= 0.05) return false;
      }

      // Check fourth jet (index 3) if it exists
      if (nJets >= 4) {
        double cJ4 = calcCentrality(3);
        double mRelJ4 = calcRelativeMass(3);
        if (cJ4 >= 0.6 || mRelJ4 >= 0.05) return false;
      }

      return true;
    }

    // -----------------------------------------------------------------------
    // passVBFSignalRegion
    //   AND of all individual VBF predicates above.  Excludes requirements
    //   that are not available in this ntuple: lepton veto, b-tagging, JVT,
    //   E_T^miss, and soft-term requirements.  passFSRCompatibility is
    //   evaluated but currently excluded from the final AND.
    // -----------------------------------------------------------------------
    bool passVBFSignalRegion() const {
      bool passMultiplicity = passJetMultiplicity();
      bool passLeadPt = passLeadingJetPt();
      bool passDeltaPhi = passJetDeltaPhi();
      bool passOppositeEta = passOppositeHemispheres();
      bool passDeltaEta = passLargeDeltaEta();
      bool passMjj = passLargeDijetMass();
      bool passPtSum = passAllJetsPtSum();
      bool passForward = passForwardJet();
      bool passFSR = passFSRCompatibility();

      if (DEBUG) {
        std::cout << "VBF Selection:" << std::endl;
        std::cout << "  Jet multiplicity (2-4 jets, pT>25): " << passMultiplicity << std::endl;
        std::cout << "  Leading jet pT (>60,>40 GeV): " << passLeadPt << std::endl;
        std::cout << "  Delta phi_jj (<2): " << passDeltaPhi << std::endl;
        std::cout << "  Opposite hemispheres: " << passOppositeEta << std::endl;
        std::cout << "  Delta eta_jj (>3.0): " << passDeltaEta << std::endl;
        std::cout << "  Dijet mass (>600 GeV): " << passMjj << std::endl;
        std::cout << "  All jets pT sum (>140 GeV): " << passPtSum << std::endl;
        std::cout << "  At least one forward jet (2.38<|eta|<4.0): " << passForward << std::endl;
        std::cout << "  FSR compatibility (C_i<0.6, m_rel<0.05 for j3,j4): " << passFSR << std::endl;
      }

      return passMultiplicity && passLeadPt && passDeltaPhi &&
             passOppositeEta && passDeltaEta && passMjj && passPtSum && passForward; // && passFSR;
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
//     updateScores      — computes TRKPTZ, CALO*, TESTML, TEST_MISCL scores
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
    std::map<Score,double> scores;
    double purity = 0.0;
    bool maxPtCluster = false;
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
    //   result in this->purity.  Called only when TEST_MISCL is active
    //   (calcPurityFlag == true) to avoid unnecessary truth-matching work.
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
    //   TESTML   — Output of the DNN (MLModel::predict) on the 8 normalised
    //              features returned by calcFeatures.
    //   TEST_MISCL — Copies TRKPTZ so that the same selection is used; the
    //              purity gate is applied later in event_processing.h.
    //   CALO90 / CALO60 — Also copy TRKPTZ (calo-time proximity filtering
    //              is applied in chooseCluster, not here).
    // -----------------------------------------------------------------------
    void updateScores(BranchPointerWrapper *branch, MLModel *mlModel) {
      // Assign TRKPTZ, CALO90, CALO60, TESTML Scores

      if (this->values.size() > 1) {
	double dz = std::abs(this->values.at(1)-branch->recoVtxZ[0]);
	auto oldscore = std::pow(this->scores.at(Score::TRKPT),0.9);
        this->scores[Score::TRKPTZ] = oldscore*exp(-1.5*dz);
      } else {
	double znum=0., zden=0.;
	for (auto trk: this->trackIndices) {
	  auto
	    trkZ = branch->trackZ0[trk],
	    trkVarZ = branch->trackVarZ0[trk];
	  znum += trkZ/(trkVarZ);
	  zden += 1/(trkVarZ);
	}

	double z = znum/zden, zsigma = 1/std::sqrt(zden);
	double dz = std::abs(z-branch->recoVtxZ[0]);
        this->scores[Score::TRKPTZ] = this->scores.at(Score::TRKPT)*exp(-1.5*dz);
      }

      // ML score for TESTML
      std::vector<float> features = this->calcFeatures(branch);
      float mlScore = mlModel->predict(features);
      this->scores[TESTML] = mlScore;

      // TEST_MISCL uses TRKPTZ as its selection score; the purity gate is applied
      // at efficiency-check time in event_processing.h (both pass and total fills).
      this->scores[TEST_MISCL] = this->scores.at(Score::TRKPTZ);

      this->scores[Score::CALO90] = this->scores.at(Score::TRKPTZ);
      this->scores[Score::CALO60] = this->scores.at(Score::TRKPTZ);
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
    //   3 · PASS_SIGMA of the truth hard-scatter vertex time.
    //   This is the primary efficiency criterion used when filling the
    //   pass/total histograms in event_processing.h.
    // -----------------------------------------------------------------------
    bool passEfficiency(BranchPointerWrapper *branch) {
      if (DEBUG) std::cout << "Choosing pass score\n";
      if (this->values.size() == 0)
	return false;

      double diff = std::abs(this->values.at(0)-branch->truthVtxTime[0]);
      if (diff > 3*PASS_SIGMA)
	return false;

      return true;
      // int nHSTrack = 0;
      // for (auto idx: this->trackIndices) {
      // 	if (idx == -1 or branch->trackToTruthvtx[idx] == 0)
      // 	  nHSTrack++;
      // }

      // return nHSTrack > 0;
    }

    // -----------------------------------------------------------------------
    // calcFeatures
    //   Extracts the 10 features used as DNN input and applies the stored
    //   normalisation (means/stds from training) to return a unit-normalised
    //   feature vector.  Features (in training order):
    //     0  delta_z              — cluster z relative to reco primary vertex
    //     1  delta_z_resunits     — delta_z in units of cluster_z_sigma
    //     2  cluster_z_sigma      — precision-weighted z₀ uncertainty
    //     3  cluster_d0           — precision-weighted transverse impact parameter
    //     4  cluster_d0_sigma     — d0 uncertainty
    //     5  cluster_qOverP       — precision-weighted charge-over-momentum
    //     6  cluster_qOverP_sigma — q/p uncertainty
    //     7  cluster_sumpt        — scalar sum of constituent track pT
    //     8  cluster_time_sigma   — precision-weighted HGTD time uncertainty
    //                               (= sigmas[0], analogous to cluster_z_sigma)
    //     9  cluster_n_tracks     — number of constituent tracks
    //   Normalisation parameters are stored as static arrays so they are
    //   allocated only once regardless of how many clusters are evaluated.
    //
    //     10 cluster_dR           — min ΔR between the pT-weighted cluster
    //                               centroid (η,φ) and the nearest topo jet
    //                               passing MIN_JETPT
    //   NOTE: MEANS[10], STDS[10] are placeholders — update from
    //   normalization_params.json after retraining with the new feature.
    // -----------------------------------------------------------------------
    std::vector<float> calcFeatures(BranchPointerWrapper *branch) {
      float znum  = 0., zden  = 0.;
      float dnum  = 0., dden  = 0.;
      float qpnum = 0., qpden = 0.;
      float sumpt = 0.;
      for (auto trk: trackIndices) {
        auto
            trkZ = branch->trackZ0[trk], trkVarZ = branch->trackVarZ0[trk],
	    trkD = branch->trackD0[trk], trkVarD = branch->trackVarD0[trk],
	    trkQ = branch->trackQP[trk], trkVarQ = branch->trackVarQp[trk],
	    trkPt = branch->trackPt[trk];

        znum += trkZ / (trkVarZ);
        zden += 1 / (trkVarZ);

        dnum += trkD / (trkVarD);
        dden += 1 / (trkVarD);

	qpnum += trkQ / (trkVarQ);
        qpden += 1 / (trkVarQ);

	sumpt += trkPt;
      }

      // Calculate cluster properties
      float clusterZ = znum / zden;
      float clusterZSigma = 1.0f / std::sqrt(zden);
      float clusterD0 = dnum / dden;
      float clusterD0Sigma = 1.0f / std::sqrt(dden);
      float clusterQOverP = qpnum / qpden;
      float clusterQOverPSigma = 1.0f / std::sqrt(qpden);

      // Feature 8: precision-weighted time uncertainty
      // sigmas[0] is set by makeSimpleClusters to 1/sqrt(sum(1/sigma_t^2))
      float clusterTimeSigma = static_cast<float>(this->sigmas.at(0));

      // Feature 9: constituent track multiplicity
      float clusterNTracks   = static_cast<float>(this->trackIndices.size());

      // Feature 10: min ΔR between pT-weighted cluster centroid and nearest jet
      // Mirrors clusterMinDR() in export_training_data.cxx exactly.
      float clusterDR = -1.0f;
      {
        double sumPt = 0., etaSum = 0., phiSum = 0.;
        for (int trk : this->trackIndices) {
          double pt  = branch->trackPt[trk];
          double eta = branch->trackEta[trk];
          double phi = branch->trackPhi[trk];
          sumPt  += pt;
          etaSum += pt * eta;
          phiSum += pt * phi;
        }
        if (sumPt > 0.) {
          double cEta  = etaSum / sumPt;
          double cPhi  = phiSum / sumPt;
          double minDR = 1e9;
          for (int j = 0; j < (int)branch->topoJetPt.GetSize(); ++j) {
            if (branch->topoJetPt[j] < MIN_JETPT) continue;
            double deta = branch->topoJetEta[j] - cEta;
            double dphi = TVector2::Phi_mpi_pi(branch->topoJetPhi[j] - cPhi);
            double dR   = std::sqrt(deta*deta + dphi*dphi);
            if (dR < minDR) minDR = dR;
          }
          if (minDR < 1e8) clusterDR = static_cast<float>(minDR);
        }
      }

      // Reference vertex (primary vertex)
      float refVtxZ = branch->recoVtxZ[0];

      // Calculate delta_z and delta_z in resolution units
      float deltaZ = clusterZ - refVtxZ;
      float deltaZResunits = deltaZ / clusterZSigma;

      // Normalization parameters (static: allocated once, never reconstructed)
      static const float MEANS[11] = {
        0.000181542369520147,  // delta_z mean
        0.0022696158408701665,  // delta_z_resunits mean
        1.0721281960578108,  // cluster_z_sigma mean
        -0.0025243291487635593,  // cluster_d0 mean
        0.08665513919653857,  // cluster_d0_sigma mean
        1.9620804663321413e-07,  // cluster_qOverP mean
        4.9269252726781045e-06,  // cluster_qOverP_sigma mean
        14.154102551101227,  // cluster_sumpt mean
        17.277764049475678,  // cluster_time_sigma mean
        5.649361497551156,  // cluster_n_tracks mean
        1.7420735517923134  // cluster_dR mean
      };
      static const float STDS[11] = {
        2.7448773941025766,  // delta_z std
        1.7700623344292996,  // delta_z_resunits std
        1.1223853949974232,  // cluster_z_sigma std
        0.42576499860027983,  // cluster_d0 std
        0.11750566161560266,  // cluster_d0_sigma std
        4.238450693042457e-05,  // cluster_qOverP std
        4.359020547359499e-06,  // cluster_qOverP_sigma std
        21.517411638686664,  // cluster_sumpt std
        9.664553181701693,  // cluster_time_sigma std
        6.224281332216253,  // cluster_n_tracks std
        1.0543692647024774  // cluster_dR std
      };

      // Return normalized features in training order
      std::vector<float> features = {
        (deltaZ             - MEANS[0])  / STDS[0],   // Feature 0
        (deltaZResunits     - MEANS[1])  / STDS[1],   // Feature 1
        (clusterZSigma      - MEANS[2])  / STDS[2],   // Feature 2
        (clusterD0          - MEANS[3])  / STDS[3],   // Feature 3
        (clusterD0Sigma     - MEANS[4])  / STDS[4],   // Feature 4
        (clusterQOverP      - MEANS[5])  / STDS[5],   // Feature 5
        (clusterQOverPSigma - MEANS[6])  / STDS[6],   // Feature 6
        (sumpt              - MEANS[7])  / STDS[7],   // Feature 7
        (clusterTimeSigma   - MEANS[8])  / STDS[8],   // Feature 8
        (clusterNTracks     - MEANS[9])  / STDS[9],   // Feature 9
        (clusterDR          - MEANS[10]) / STDS[10],  // Feature 10
      };

      return features;
    }
  };
  
}
#endif // STRUCTS_H
