#ifndef STRUCTS_H
#define STRUCTS_H

#include "clustering_includes.h"
#include "clustering_constants.h"
#include "ml_model.h"
#include <cstdlib>

namespace MyUtl {

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

    bool passForwardHsTracks(int& nForwardHSTrack) {
      return
	MAX_NHS_TRACK >= nForwardHSTrack and
	nForwardHSTrack >= MIN_NHS_TRACK;
    }

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

    // ============================================================================
    // VBF H->Invisible Signal Region Selection Functions
    // Based on JHEP08(2022)104 - https://doi.org/10.1007/JHEP08(2022)104
    //
    // NOTE: These functions assume jets are pT-sorted in descending order
    //       (verified with check_jet_sorting.cxx - 100% of events are sorted)
    //       Leading jet = topoJetPt[0], Subleading jet = topoJetPt[1]
    // ============================================================================

    // Count jets passing pT threshold of 25 GeV
    int countJetsAbovePt(double ptThreshold = 25.0) const {
      int count = 0;
      for (int jetIdx = 0; jetIdx < this->topoJetPt.GetSize(); ++jetIdx) {
        if (this->topoJetPt[jetIdx] > ptThreshold) {
          count++;
        }
      }
      return count;
    }

    // Check if event has at least 2 jets with pT > 25 GeV
    bool passJetMultiplicity() const {
      int nJets = countJetsAbovePt(25.0);
      return (nJets >= 2);
    }

    // Check leading jet pT > 60 GeV and subleading jet pT > 40 GeV
    // Original ATLAS cuts: leading > 80 GeV, subleading > 50 GeV (loosened for timing studies)
    bool passLeadingJetPt() const {
      if (this->topoJetPt.GetSize() < 2) return false;

      // Jets are assumed to be pT-sorted
      double leadingPt = this->topoJetPt[0];
      double subleadingPt = this->topoJetPt[1];

      return (leadingPt > 60.0 && subleadingPt > 40.0);
    }

    // Check if two leading jets are not back-to-back: Δφ_jj < 2
    bool passJetDeltaPhi() const {
      if (this->topoJetPt.GetSize() < 2) return false;

      double phi1 = this->topoJetPhi[0];
      double phi2 = this->topoJetPhi[1];
      double deltaPhi = std::abs(TVector2::Phi_mpi_pi(phi1 - phi2));

      return (deltaPhi < 2.0);
    }

    // Check VBF topology: opposite hemispheres (η_j1 * η_j2 < 0)
    bool passOppositeHemispheres() const {
      if (this->topoJetEta.GetSize() < 2) return false;

      double eta1 = this->topoJetEta[0];
      double eta2 = this->topoJetEta[1];

      return (eta1 * eta2 < 0);
    }

    // Check VBF topology: large pseudorapidity separation (Δη_jj > 3.0)
    // Original ATLAS cut: Δη_jj > 3.8 (loosened for timing studies)
    bool passLargeDeltaEta(double minDeltaEta = 3.0) const {
      if (this->topoJetEta.GetSize() < 2) return false;

      double eta1 = this->topoJetEta[0];
      double eta2 = this->topoJetEta[1];
      double deltaEta = std::abs(eta1 - eta2);

      return (deltaEta > minDeltaEta);
    }

    // Calculate invariant mass of two leading jets
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

    // Check VBF topology: large invariant mass (m_jj > 0.3 TeV)
    // Original ATLAS cut: m_jj > 0.8 TeV (loosened for timing studies)
    bool passLargeDijetMass(double minMass = 300.0) const {
      double mjj = calcDijetMass();
      return (mjj > minMass);
    }

    // Calculate scalar sum of all jet pT (for pT^all-jets requirement)
    double calcAllJetsPtSum() const {
      double ptSum = 0.0;
      for (int jetIdx = 0; jetIdx < this->topoJetPt.GetSize(); ++jetIdx) {
        ptSum += this->topoJetPt[jetIdx];
      }
      return ptSum;
    }

    // Check if sum of all jet pT > 140 GeV
    bool passAllJetsPtSum(double minPtSum = 140.0) const {
      double ptSum = calcAllJetsPtSum();
      return (ptSum > minPtSum);
    }

    // Check if at least one of the two leading jets is in HGTD acceptance
    bool passForwardJet() const {
      if (this->topoJetEta.GetSize() < 2) return false;

      double eta1 = std::abs(this->topoJetEta[0]);
      double eta2 = std::abs(this->topoJetEta[1]);

      bool jet1Forward = (eta1 > MIN_HGTD_ETA && eta1 < MAX_HGTD_ETA);
      bool jet2Forward = (eta2 > MIN_HGTD_ETA && eta2 < MAX_HGTD_ETA);

      return (jet1Forward || jet2Forward);
    }

    // Calculate centrality C_i for jet i (for i = j3, j4)
    // C_i quantifies how central jet i is relative to the two leading jets in rapidity
    double calcCentrality(int jetIdx) const {
      if (this->topoJetEta.GetSize() < 2) return -1.0;
      if (jetIdx >= this->topoJetEta.GetSize()) return -1.0;

      double eta_j1 = this->topoJetEta[0];
      double eta_j2 = this->topoJetEta[1];
      double eta_i = this->topoJetEta[jetIdx];

      // Central position between the two leading jets
      double eta_center = (eta_j1 + eta_j2) / 2.0;
      double delta_eta_jj = eta_j1 - eta_j2;

      // Calculate centrality
      double numerator = eta_i - eta_center;
      double C_i = std::exp(-4.0 / (delta_eta_jj * delta_eta_jj) * numerator * numerator);

      return C_i;
    }

    // Calculate relative mass m_rel for jet i
    // m_rel_i = min{m_j1i, m_j2i} / m_jj
    double calcRelativeMass(int jetIdx) const {
      if (this->topoJetPt.GetSize() < 2) return -1.0;
      if (jetIdx >= this->topoJetPt.GetSize()) return -1.0;

      // Build TLorentzVectors (assume massless jets)
      TLorentzVector jet1, jet2, jetI;
      jet1.SetPtEtaPhiM(this->topoJetPt[0], this->topoJetEta[0], this->topoJetPhi[0], 0.0);
      jet2.SetPtEtaPhiM(this->topoJetPt[1], this->topoJetEta[1], this->topoJetPhi[1], 0.0);
      jetI.SetPtEtaPhiM(this->topoJetPt[jetIdx], this->topoJetEta[jetIdx], this->topoJetPhi[jetIdx], 0.0);

      // Calculate invariant masses
      double m_jj = (jet1 + jet2).M();
      double m_j1i = (jet1 + jetI).M();
      double m_j2i = (jet2 + jetI).M();

      // Return relative mass
      double m_rel = std::min(m_j1i, m_j2i) / m_jj;
      return m_rel;
    }

    // Check FSR compatibility for third and fourth jets (if they exist)
    // Requires C_i < 0.6 and m_rel_i < 0.05
    bool passFSRCompatibility() const {
      int nJets = this->topoJetPt.GetSize();

      // If only 2 jets, automatically pass (no FSR jets to check)
      if (nJets <= 2) return true;

      // Check third jet (index 2)
      if (nJets >= 3) {
        double C_j3 = calcCentrality(2);
        double m_rel_j3 = calcRelativeMass(2);
        if (C_j3 >= 0.6 || m_rel_j3 >= 0.05) return false;
      }

      // Check fourth jet (index 3) if it exists
      if (nJets >= 4) {
        double C_j4 = calcCentrality(3);
        double m_rel_j4 = calcRelativeMass(3);
        if (C_j4 >= 0.6 || m_rel_j4 >= 0.05) return false;
      }

      return true;
    }

    // Combined VBF signal region selection
    // Returns true if event passes all VBF H->invisible cuts
    // (excluding lepton veto, b-tagging, JVT, EmissT, and soft term requirements)
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

    void calcPurity(BranchPointerWrapper *branch) {
      double num = 0.0, denom = 0.0;
      for (auto trk: this->trackIndices) {
	if (branch->trackToTruthvtx[trk] == 0) num += branch->trackPt[trk];
	denom += branch->trackPt[trk];
      }
      this->purity = num/denom;
    }

    void updateScores(BranchPointerWrapper *branch, MLModel *ml_model) {
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
      float mlScore = ml_model->predict(features);
      this->scores[TESTML] = mlScore;

      // TEST_MISCL uses TRKPTZ as its selection score; the purity gate is applied
      // at efficiency-check time in event_processing.h (both pass and total fills).
      this->scores[TEST_MISCL] = this->scores.at(Score::TRKPTZ);

      this->scores[Score::CALO90] = this->scores.at(Score::TRKPTZ);
      this->scores[Score::CALO60] = this->scores.at(Score::TRKPTZ);
    }
    
    double timeSpread() {
      // calculate stdev of times
      auto avg = std::accumulate(allTimes.begin(),allTimes.end(),0.0)/this->allTimes.size();
      auto ssd = 0.0;
      for (auto t: allTimes) {
	ssd += (avg - t)*(avg - t);
      }

      return std::sqrt(ssd/this->allTimes.size());
    }

    double zSpread(BranchPointerWrapper *bpw) {
      // calculate stdev of times
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

    std::vector<float> calcFeatures(BranchPointerWrapper *branch) {
      float znum  = 0., zden  = 0.;
      float dnum  = 0., dden  = 0.;
      float qpnum = 0., qpden = 0.;
      float sumpt = 0.;
      for (auto trk: trackIndices) {
        auto
            trk_z = branch->trackZ0[trk], trk_var_z = branch->trackVarZ0[trk],
	    trk_d = branch->trackD0[trk], trk_var_d = branch->trackVarD0[trk],
	    trk_q = branch->trackQP[trk], trk_var_q = branch->trackVarQp[trk],
	    trk_pt = branch->trackPt[trk];

        znum += trk_z / (trk_var_z);
        zden += 1 / (trk_var_z);

        dnum += trk_d / (trk_var_d);
        dden += 1 / (trk_var_d);

	qpnum += trk_q / (trk_var_q);
        qpden += 1 / (trk_var_q);

	sumpt += trk_pt;
      }

      // Calculate cluster properties
      float cluster_z = znum / zden;
      float cluster_z_sigma = 1.0f / std::sqrt(zden);
      float cluster_d0 = dnum / dden;
      float cluster_d0_sigma = 1.0f / std::sqrt(dden);
      float cluster_qOverP = qpnum / qpden;
      float cluster_qOverP_sigma = 1.0f / std::sqrt(qpden);

      // Reference vertex (primary vertex)
      float ref_vtx_z = branch->recoVtxZ[0];

      // Calculate delta_z and delta_z in resolution units
      float delta_z = cluster_z - ref_vtx_z;
      float delta_z_resunits = delta_z / cluster_z_sigma;

      // Normalization parameters (static: allocated once, never reconstructed)
      static const float means[8] = {
        0.6658103458145465f,      // delta_z mean
        1.4062413922431898f,      // delta_z_resunits mean
        0.4384938254278939f,      // cluster_z_sigma mean
        -0.0006795810095683315f,  // cluster_d0 mean
        0.03825910434563131f,     // cluster_d0_sigma mean
        2.5463849668322525e-07f,  // cluster_qOverP mean
        2.249843783028648e-06f,   // cluster_qOverP_sigma mean
        36.363541536390116f,      // cluster_sumpt mean
      };
      static const float stds[8] = {
        1.0862063628195677f,      // delta_z std
        1.0457992632101616f,      // delta_z_resunits std
        0.5121707999068104f,      // cluster_z_sigma std
        0.16671231081127164f,     // cluster_d0 std
        0.034856748414000265f,    // cluster_d0_sigma std
        2.33339973420475e-05f,    // cluster_qOverP std
        2.0502623770486916e-06f,  // cluster_qOverP_sigma std
        33.386926190214176f,      // cluster_sumpt std
      };

      // Return normalized features in training order
      std::vector<float> features = {
        (delta_z              - means[0]) / stds[0],  // Feature 0
        (delta_z_resunits     - means[1]) / stds[1],  // Feature 1
        (cluster_z_sigma      - means[2]) / stds[2],  // Feature 2
        (cluster_d0           - means[3]) / stds[3],  // Feature 3
        (cluster_d0_sigma     - means[4]) / stds[4],  // Feature 4
        (cluster_qOverP       - means[5]) / stds[5],  // Feature 5
        (cluster_qOverP_sigma - means[6]) / stds[6],  // Feature 6
        (sumpt                - means[7]) / stds[7],  // Feature 7
      };

      return features;
    }
  };
  
}
#endif // STRUCTS_H
