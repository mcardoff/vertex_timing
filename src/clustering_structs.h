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
      for(int jetIdx = 0; jetIdx < this->topoJetEta.GetSize(); ++jetIdx) {
	float
	  eta = std::abs(topoJetEta[jetIdx]),
	  pt = truthHSJetPt[jetIdx];
	if (DEBUG) std::cout << "pt, eta: " << pt << ", " << eta << '\n';
	bool passpt = pt > MIN_JETPT;
	bool passpteta = passpt and eta > MIN_ABS_ETA_JET and eta < MAX_ABS_ETA_JET;
	passptcount += passpt ? 1 : 0;
	passptetacount += passpteta ? 1 : 0;
      }
      bool aPasses = passptcount >= MIN_PASSPT_JETS;
      bool bPasses = passptetacount >= MIN_PASSETA_JETS;
      return aPasses and bPasses;
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
	    pt > MIN_TRACK_PT and pt < MAX_TRACK_PT and
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

      // have to load model, calculate parameters, and update value in map
      std::vector<float> features = this->calcFeatures(branch);
      this->scores[TESTML] = ml_model->predict(features);
      
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

      // Return features in training order
      std::vector<float> features = {
          delta_z,              // Feature 0
          delta_z_resunits,     // Feature 1
          cluster_z_sigma,      // Feature 2
          cluster_d0,           // Feature 3
          cluster_d0_sigma,     // Feature 4
          cluster_qOverP,       // Feature 5
          cluster_qOverP_sigma, // Feature 6
          sumpt                 // Feature 7
      };

      const std::vector<float> means = {
	// Your scaler.mean_ values in order
	0.6658103458145465,     // delta_z mean
	1.4062413922431898,     // delta_z_resunits mean
	0.4384938254278939,     // cluster_z_sigma mean
	-0.0006795810095683315, // cluster_d0 mean
	0.03825910434563131,    // cluster_d0_sigma mean
	2.5463849668322525e-07, // cluster_qOverP mean
	2.249843783028648e-06,  // cluster_qOverP_sigma mean
	36.363541536390116      // cluster_sumpt mean
      };
      

      const std::vector<float> stds = {
          // Your scaler.scale_ values in order
          1.0862063628195677,     // delta_z std
          1.0457992632101616,     // delta_z_resunits std
          0.5121707999068104,     // cluster_z_sigma std
          0.16671231081127164,    // cluster_d0 std
          0.034856748414000265,   // cluster_d0_sigma std
          2.33339973420475e-05,   // cluster_qOverP std
          2.0502623770486916e-06, // cluster_qOverP_sigma std
          33.386926190214176      // cluster_sumpt std
      };
      
    
    // Normalize features
    for (size_t i = 0; i < features.size(); i++) {
      features[i] = (features[i] - means[i]) / stds[i];
    }

      return features;
    }
  };
  
}
#endif // STRUCTS_H
