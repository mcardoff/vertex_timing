#ifndef STRUCTS_H
#define STRUCTS_H

#include "clustering_includes.h"
#include "clustering_constants.h"
#include <functional>

namespace myutl {

  struct BranchPointerWrapper {
    TTreeReader& reader;

    TTreeReaderValue<float> weight;

    TTreeReaderArray<float> track_z0;
    TTreeReaderArray<float> track_d0;
    TTreeReaderArray<float> track_pt;
    TTreeReaderArray<float> track_qp;
    TTreeReaderArray<float> track_eta;
    TTreeReaderArray<float> track_phi;
    TTreeReaderArray<float> track_theta;
    TTreeReaderArray<float> track_var_z0;
    TTreeReaderArray<float> track_var_d0;
    TTreeReaderArray<float> track_var_qp;
    TTreeReaderArray<float> track_var_theta;
    TTreeReaderArray<float> track_time;
    TTreeReaderArray<float> track_time_res;
    TTreeReaderArray<int>   track_time_valid;
    TTreeReaderArray<int>   track_to_truthvtx;
    TTreeReaderArray<int>   track_to_particle;
    TTreeReaderArray<bool>  track_quality;
    TTreeReaderArray<int>   track_hgtd_hits;
    TTreeReaderArray<int>   track_prim_hits;
    TTreeReaderArray<float> track_near_idx;
    TTreeReaderArray<float> track_near_z0sin;
    TTreeReaderArray<float> track_near_z0sin_unc;
    TTreeReaderArray<float> track_near_sig;

    TTreeReaderArray<float> truth_vtx_z;
    TTreeReaderArray<float> truth_vtx_time;
    TTreeReaderArray<bool>  truth_vtx_ishs;

    TTreeReaderArray<float> reco_vtx_z;
    TTreeReaderArray<float> reco_vtx_time;
    TTreeReaderArray<float> reco_vtx_timeRes;
    TTreeReaderArray<int>   reco_vtx_valid;

    TTreeReaderArray<float> topojet_pt;
    TTreeReaderArray<float> topojet_eta;
    TTreeReaderArray<float> topojet_phi;

    TTreeReaderArray<float> truthhsjet_pt;
    TTreeReaderArray<float> truthhsjet_eta;

    TTreeReaderArray<float> particle_t;

  BranchPointerWrapper(TTreeReader& r)
    : reader (r),
      weight (r, "weight"),
      track_d0 (r, "Track_d0"),
      track_z0 (r, "Track_z0"),
      track_pt (r, "Track_pt"),
      track_eta (r, "Track_eta"),
      track_qp (r, "Track_qOverP"),
      track_theta (r, "Track_theta"),
      track_phi (r, "Track_phi"),
      track_var_d0 (r, "Track_var_d0"),
      track_var_z0 (r, "Track_var_z0"),
      track_var_qp (r, "Track_var_qOverP"),
      track_var_theta (r, "Track_var_theta"),
      track_time (r, "Track_time"),
      track_time_res (r, "Track_timeRes"),
      track_time_valid (r, "Track_hasValidTime"),
      track_quality (r, "Track_quality"),
      track_to_truthvtx (r, "Track_truthVtx_idx"),
      track_to_particle (r, "Track_truthPart_idx"),
      track_hgtd_hits (r, "Track_nHGTDHits"),
      track_prim_hits (r, "Track_nHGTDPrimaryHits"),
      track_near_idx (r, "Track_nearestVtx_idx"),
      track_near_z0sin (r, "Track_nearestVtx_z0SinTheta"),
      track_near_z0sin_unc (r, "Track_nearestVtx_z0SinThetaUncertainty"),
      track_near_sig (r, "Track_nearestVtx_sig"),
      truth_vtx_z (r, "TruthVtx_z"),
      truth_vtx_time (r, "TruthVtx_time"),
      truth_vtx_ishs (r, "TruthVtx_isHS"),
      reco_vtx_z (r, "RecoVtx_z"),
      reco_vtx_time (r, "RecoVtx_time"),
      reco_vtx_timeRes (r, "RecoVtx_timeRes"),
      reco_vtx_valid (r, "RecoVtx_hasValidTime"),
      topojet_pt (r, "AntiKt4EMTopoJets_pt"),
      topojet_eta (r, "AntiKt4EMTopoJets_eta"),
      topojet_phi (r, "AntiKt4EMTopoJets_phi"),
      truthhsjet_pt (r, "TruthHSJet_pt"),
      truthhsjet_eta (r, "TruthHSJet_eta"),
      particle_t (r, "TruthPart_prodVtx_time")
    {}

    bool pass_basic_cuts() {
      if (this->truthhsjet_pt.GetSize() < min_jets or
	  this->topojet_pt.GetSize() < min_jets) {
	if (debug) std::cout << "Skipping low jet event" << std::endl;
        return false;
      }
    
      // check reco HS vertex is with 2mm of truth HS vertex
      if(std::abs(truth_vtx_z[0] - reco_vtx_z[0]) > max_vtx_dz) {
	if(debug) std::cout << "Skipping event due to incorrect HS vertex" << std::endl;
	return false;
      }
      
      return true;
    }

    bool pass_jet_pt_cut() { 
      int passptcount = 0, passptetacount = 0;
      for(int jet_idx = 0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
	float
	  eta = topojet_eta[jet_idx],
	  pt = truthhsjet_pt[jet_idx];
	if (debug) std::cout << "pt, eta: " << pt << ", " << eta << std::endl;
	bool passpt = pt > min_jetpt;
	bool passpteta = pt > min_jetpt and std::abs(eta) > min_abs_eta_jet;
	passptcount += passpt ? 1 : 0;
	passptetacount += passpteta ? 1 : 0;
      }
      bool pass_a = passptcount >= min_passpt_jets;
      bool pass_b = passptetacount >= min_passeta_jets;
      return pass_a and pass_b;
    }

    bool pass_forward_hs_tracks(int& nForwardHSTrack) {
      return
	max_nhs_track >= nForwardHSTrack and
	nForwardHSTrack >= min_nhs_track;
    }

    void count_forward_jets(int& nForwardJet) {
    nForwardJet = 0;
    for(int jet_idx = 0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
      float
	jet_eta = this->topojet_eta[jet_idx],
	jet_pt = this->topojet_pt[jet_idx];
      if (std::abs(jet_eta) > min_abs_eta_jet and
	  jet_pt > min_jetpt)
	nForwardJet++;
    }
  }
    
    void count_forward_tracks(int& nFTrack, int& nFTrack_HS, int& nFTrack_PU,
			      std::vector<int> tracks, bool check_time_valid) {
      nFTrack = 0; nFTrack_HS = 0; nFTrack_PU = 0;
      for(auto trk_idx: tracks) {
	double eta = this->track_eta[trk_idx];
	double pt = this->track_pt[trk_idx];
	int  truthvtx = this->track_to_truthvtx[trk_idx];
	bool quality = this->track_quality[trk_idx] == true;
	bool check_query = (this->track_time_valid[trk_idx] == 1);
	// already know these pass association
	if (std::abs(eta) > min_abs_eta_track and
	    pt > min_track_pt and pt < max_track_pt and
	    quality and check_query) {
	  nFTrack++;
	  if (truthvtx != -1 and this->truth_vtx_ishs[truthvtx])
	    nFTrack_HS++;
	  else
	    nFTrack_PU++;
	}
      }
    }

    double calc_jetpt_dr_score(int trk_idx) {
      double
	trk_eta = this->track_eta[trk_idx],
	trk_phi = this->track_phi[trk_idx];

      double min_dR = 1e6;
      int min_idx = -1;
      for (int jet_idx=0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
	double
	  jet_eta = this->topojet_eta[jet_idx],
	  jet_phi = this->topojet_phi[jet_idx];
	double
	  deta = jet_eta-trk_eta,
	  dphi = TVector2::Phi_mpi_pi(jet_phi - trk_phi);
	double this_dR = std::hypot(deta, dphi);
	if (this_dR < min_dR) {
	  min_dR = this_dR;
	  min_idx = jet_idx;
	}
      }
      double returnScore = this->topojet_pt[min_idx]*std::exp(-min_dR);
      return returnScore;
    }

    double calc_trkpt_dr_score(int trk_idx) {
      double
	trk_eta = this->track_eta[trk_idx],
	trk_phi = this->track_phi[trk_idx];


      double min_dR = 1e6;
      for (int jet_idx=0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
	double
	  jet_eta = this->topojet_eta[jet_idx],
	  jet_phi = this->topojet_phi[jet_idx];

	double
	  deta = jet_eta-trk_eta,
	  dphi = TVector2::Phi_mpi_pi(jet_phi - trk_phi);

	double this_dR = std::sqrt(deta*deta + dphi*dphi);
	if (this_dR < min_dR) {
	  min_dR = this_dR;
	}
      }
      double returnScore = this->track_pt[trk_idx]*std::exp(-min_dR);
      return returnScore;
    }
  };

  struct Cluster {
    std::vector<double> values;
    std::vector<double> sigmas;
    std::vector<double> all_times;
    std::vector<int> track_indices;
    std::map<ScoreType,double> scores;
    double purity = 0.0;
    bool max_pt_cluster = false;
    bool was_merged = false;
    int n_constituents=1;

    bool operator==(const Cluster& other) {
      bool same_values = values.at(0) == other.values.at(0);
      bool same_sigmas = sigmas.at(0) == other.sigmas.at(0);
      bool same_consts = n_constituents == other.n_constituents;
      // this SHOULD be sufficient
      return same_consts and same_sigmas and same_values;
    }

    bool operator!=(const Cluster& other) {
      return !(*this == other);
    }

    void calcPurity(BranchPointerWrapper *branch) {
      double num = 0.0, denom = 0.0;
      for (auto trk: this->track_indices) {
	if (branch->track_to_truthvtx[trk] == 0) num += branch->track_pt[trk];
	denom += branch->track_pt[trk];
      }
      this->purity = num/denom;
    }

    double spread_t() {
      // calculate stdev of times
      auto avg = std::accumulate(all_times.begin(),all_times.end(),0.0)/this->all_times.size();
      auto ssd = 0.0;
      for (auto t: all_times) {
	ssd += (avg - t)*(avg - t);
      }

      return std::sqrt(ssd/this->all_times.size());
    }

    double spread_z(BranchPointerWrapper *bpw) {
      // calculate stdev of times
      auto avg = 0.0;
      for (auto trk: track_indices)
	avg += bpw->track_z0[trk];
      avg *= 1./this->track_indices.size();
      auto ssd = 0.0;
      for (auto t: track_indices) {
	ssd += (avg - bpw->track_z0[t])*(avg - bpw->track_z0[t]);
      }
      return std::sqrt(ssd/this->track_indices.size());
    }

    bool passEfficiency(BranchPointerWrapper *branch) {
      if (debug) std::cout << "Choosing pass score" << std::endl;
      if (this->values.size() == 0)
	return false;

      double diff = std::abs(this->values.at(0)-branch->truth_vtx_time[0]);
      if (diff > 60)
	return false;

      int nHSTrack = 0;
      for (auto idx: this->track_indices) {
	if (branch->track_to_truthvtx[idx] == 0)
	  nHSTrack++;
      }

      return nHSTrack > 0;
    }
  };
  
}
#endif // STRUCTS_H
