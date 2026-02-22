#include <RtypesCore.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector2.h>
#include <TDirectory.h>
#include <TLatex.h>
#include <TColor.h>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <boost/filesystem.hpp>

// Clustering infrastructure
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"

#define debug false

// cut variables
auto min_jets          =  2;    // min number of jets
auto min_jetpt         =  30.0; // self explanatory
auto min_abs_eta       =  2.0;  // min eta for a jet to be considered "forward"
auto min_abs_track_eta =  2.4;  // min eta for a track to be considered "forward"
auto min_track_pt      =  1.0;  // track_pt > 1.0 GeV
auto max_vtx_dz        =  2.0;  // max error for reco HS vertex z
auto max_nsigma = 3.0;          // how close a track can be to PV
auto t_cut = 60;

// Dictionary generation for standalone compilation
// (not needed for ROOT CINT macros)

TGraph* generate_roc(TH1D* PU_hist, TH1D* HS_hist) {
  int bin = PU_hist->GetNbinsX();
  std::vector<float> vector_x, vector_y;
  for (int i = 2; i <= bin; ++i) {
    double HS_eff = HS_hist->Integral(i, bin+1) / HS_hist->Integral();
    double PU_mistag = PU_hist->Integral(i, bin+1) / PU_hist->Integral();
    
    if (std::abs(PU_mistag) > 1e-6 && (HS_eff < 0.99)) {
      vector_x.push_back(HS_eff);
      vector_y.push_back(1.0 / PU_mistag);
    }
  }

  TGraph* graph = new TGraph(vector_x.size(),&vector_x[0],&vector_y[0]);

  return graph;
}

void setup_chain(TChain &chain) {
  // VBF H->Invisible sample
  for (const auto& entry : boost::filesystem::directory_iterator("../../ntuple-hgtd")) {
    if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
    chain.Add(entry.path().c_str());
    // break; // add 1
  }
  
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found in directory: " << std::endl;
    return;
  }
}

// Add before main()
struct ExampleEvent {
  TString file_num;
  Long64_t event_num;
  double extra_time;
  double metric; // For sorting
};

int main() {
  gStyle->SetOptStat(0);
  // gInterpreter not available in standalone mode
  // colors
  auto c01 = TColor::GetColor("#3f90da");
  auto c02 = TColor::GetColor("#ffa90e");
  auto c03 = TColor::GetColor("#bd1f01");
  auto c04 = TColor::GetColor("#94a4a2");
  auto c05 = TColor::GetColor("#832db6");
  auto c06 = TColor::GetColor("#a96b59"); 
  auto c07 = TColor::GetColor("#94a4a2");
  auto c08 = TColor::GetColor("#e76300");
  auto c09 = TColor::GetColor("#b9ac70");
  auto c10 = TColor::GetColor("#717581");
  auto c11 = TColor::GetColor("#92dadd");

// Add at the start of main()
  std::vector<ExampleEvent> examples_rescued_by_time;  // Case 1
  std::vector<ExampleEvent> examples_not_fixed_by_time; // Case 2
  double rpt_cut = 0.5; // Standard threshold for HS/PU classification  

  TChain chain ("ntuple");
  setup_chain(chain);
  TTreeReader reader(&chain);
  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);

  // Initialize branch wrapper for clustering infrastructure
  MyUtl::BranchPointerWrapper branch(reader);

  // jet variables (still need these separately for RpT calculation)
  TTreeReaderArray<float> topojet_pt
    (reader, "AntiKt4EMTopoJets_pt");
  TTreeReaderArray<float> truthhsjet_pt
    (reader, "TruthHSJet_pt");
  TTreeReaderArray<float> topojet_eta
    (reader, "AntiKt4EMTopoJets_eta");
  TTreeReaderArray<float> truthhsjet_eta
    (reader, "TruthHSJet_eta");
  TTreeReaderArray<float> topojet_phi
    (reader, "AntiKt4EMTopoJets_phi");
  TTreeReaderArray<float> truthhsjet_phi
    (reader, "TruthHSJet_phi");

  TTreeReaderArray<std::vector<int>> topojet_hsjet_indices
    (reader, "AntiKt4EMTopoJets_truthHSJet_idx");
  TTreeReaderArray<std::vector<int>> topojet_itpujet_indices
    (reader, "AntiKt4EMTopoJets_truthITPUJet_idx");
  TTreeReaderArray<std::vector<int>> topojet_otpujet_indices
    (reader, "AntiKt4EMTopoJets_truthOOTPUJet_idx");

  // Truth vertex isHS flag
  TTreeReaderArray<bool>  truth_vtx_isHS
    (reader, "TruthVtx_isHS");

  double rpt_min = 0.0, rpt_max = 375.0;
  double diff_min = -0.1, diff_max = 100.0;
  const int bin = (rpt_max-rpt_min)/0.03;
  // Histograms
  TH1D *PU_RpT_z_hist = new TH1D("PU_RpT_z_hist",
				 "Pile-Up RpT: Z sig Tracks;RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *HS_RpT_z_hist = new TH1D("HS_RpT_z_hist",
				 "Hard Scatter RpT: Z sig Tracks;RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *PU_RpT_t_hist = new TH1D("PU_RpT_t_hist",
				 "Pile-Up RpT: Z+T sig Tracks (Cone Clustering);RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *HS_RpT_t_hist = new TH1D("HS_RpT_t_hist",
				 "Hard Scatter RpT: Z+T sig Tracks (Cone Clustering);RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *PU_RpT_hgtd_hist = new TH1D("PU_RpT_hgtd_hist",
				 "Pile-Up RpT: Z+T sig Tracks (HGTD);RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *HS_RpT_hgtd_hist = new TH1D("HS_RpT_hgtd_hist",
				 "Hard Scatter RpT: Z+T sig Tracks (HGTD);RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *PU_RpT_truth_hist = new TH1D("PU_RpT_truth_hist",
				     "Pile-Up RpT: Z+T,Truth Time;RpT;Entries",
				     bin, rpt_min, rpt_max);

  TH1D *HS_RpT_truth_hist = new TH1D("HS_RpT_truth_hist",
				     "Hard Scatter RpT: Z+T,Truth Time;RpT;Entries",
				     bin, rpt_min, rpt_max);

  TH1D *PU_idl_z_hist = new TH1D("PU_idl_z_hist",
				 "Pile-Up RpT: Z Ideal Timing;RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *HS_idl_z_hist = new TH1D("HS_idl_z_hist",
				 "Hard Scatter RpT: Z Ideal Timing;RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *PU_idl_t_hist = new TH1D("PU_idl_t_hist",
				 "Pile-Up RpT: Z+T Ideal Timing;RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *HS_idl_t_hist = new TH1D("HS_idl_t_hist",
				 "Hard Scatter RpT: Z+T Ideal Timing;RpT;Entries",
				 bin, rpt_min, rpt_max);

  TH1D *RpT_diff_015 = new TH1D("RpT_diff_015",
				"Pile-Up R_{pT}(z)-R_{pT}(z,t);RpT(z)-RpT(z,t);Entries",
				bin, diff_min, diff_max);

  TH1D *RpT_diff_020 = new TH1D("RpT_diff_020",
				"Pile-Up R_{pT}(z)-R_{pT}(z,t);RpT(z)-RpT(z,t);Entries",
				bin, diff_min, diff_max);

  TH1D *RpT_diff_all = new TH1D("RpT_diff_all",
				"Pile-Up R_{pT}(z)-R_{pT}(z,t);RpT(z)-RpT(z,t);Entries",
				bin, diff_min, diff_max);

  std::vector<std::string> fail_hs, fail_pu;

  int njets = 0, nevent = 0, nhsjets = 0;
  int n_events_clustering_attempted = 0;
  int n_events_clustering_succeeded = 0;

  while (reader.Next()) {
    if (topojet_pt.GetSize() < 1) {
      if(debug)
	std::cout << "Skipping low jet event" << std::endl;
      continue;  // Skip if no jets are present
    }

    // commented bc we look at every event
    // if(reco_vtx_valid[0] == 0) { 
    //   if(debug)
    // 	std::cout << "Skipping event where reco vertex has no valid time" << std::endl;
    //   continue;  // Skip if reco vertex has no valid time
    // }
    
    // check reco HS vertex is within 2mm of truth HS vertex
    if(std::abs(branch.truthVtxZ[0] - branch.recoVtxZ[0]) > max_vtx_dz ) {
      if(debug)
	std::cout << "Skipping event due to incorrect HS vertex" << std::endl;
      continue;
    }

    if (branch.truthHSJetPt.GetSize() < 1 or
	branch.topoJetPt.GetSize() < 1) {
      continue;
    }

    if (not branch.passJetPtCut()) continue;

    // // check if there is one forward jet with pt > 30 GeV
    // bool pass_jet_pt_cut = false;
    // for(int idx = 0; idx < topojet_eta.GetSize(); ++idx) {
    //   float eta = topojet_eta[idx], pt = topojet_pt[idx];
    //   bool forward = std::abs(eta) - min_abs_eta > 1e-6;
    //   bool truth_matched = topojet_hsjet_indices[idx].size() > 1; // check this one specifically 
    //   // if (topojet_hsjet_indices.GetSize() == 0)
    // 	// std::cout << topojet_hsjet_indices.GetSize() << std::endl;
    //   if (!pass_jet_pt_cut && pt > 30 && forward && truth_matched) {
    // 	// check if the pt > 30
    // 	pass_jet_pt_cut = true;
    // 	break;
    //   }
    // }

    // if(!pass_jet_pt_cut)
    //   continue;

    // njets += topojet_eta.GetSize();
    // nhsjets += truthhsjet_eta.GetSize();
    // nevent++;
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset();

    // Get truth vertex information
    double pri_vtx_z = branch.recoVtxZ[0];
    double tru_vtx_t = branch.truthVtxTime[0];

    // Perform clustering to get vertex time
    std::vector<int> hgtd_track_indices;
    for(int track_idx = 0; track_idx < branch.trackEta.GetSize(); track_idx++) {
      double this_eta = branch.trackEta[track_idx];
      bool hasValidTime = branch.trackTimeValid[track_idx] == 1;
      bool quality = branch.trackQuality[track_idx];

      // HGTD acceptance: 2.38 < |η| < 4.0
      bool in_hgtd_acceptance = std::abs(this_eta) > 2.38 && std::abs(this_eta) < 4.0;

      if (in_hgtd_acceptance && hasValidTime && quality) {
        hgtd_track_indices.push_back(track_idx);
      }
    }

    // Run clustering on HGTD tracks
    // Create simple clusters (one per track)
    std::map<int, double> emptyMap; // empty maps for non-smeared times
    std::vector<MyUtl::Cluster> all_clusters = MyUtl::clusterTracksInTime(hgtd_track_indices, &branch, 3.0, false, true, 30.0, true, false);

    
    // Get best cluster time (highest TRKPT scoring cluster)
    double pri_vtx_t;
    double pri_vtx_t_var;
    bool clustering_succeeded = false;

    if (!all_clusters.empty()) {
      // Find cluster with highest TRKPT score
      double best_score = -10000000.0;
      int best_idx = -1;
      for (size_t i = 0; i < all_clusters.size(); i++) {
        double score = -std::abs(all_clusters[i].values.at(0)-tru_vtx_t);
        if (score > best_score) {
          best_score = score;
          best_idx = i;
        }
      }

      if (best_idx >= 0) {
        const auto& best_cluster = all_clusters[best_idx];
        pri_vtx_t = best_cluster.values[0];  // Time is first value
        pri_vtx_t_var = best_cluster.sigmas[0] * best_cluster.sigmas[0];
        clustering_succeeded = true;
      } else {
        // Fallback
        pri_vtx_t = branch.recoVtxTime[0];
        pri_vtx_t_var = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
      }
    } else {
      // Fallback to ntuple reco vertex time if clustering fails
      pri_vtx_t = branch.recoVtxTime[0];
      pri_vtx_t_var = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
    }

    // find connected tracks to the hard scatter vertex
    std::vector<int> connected_indices_z, connected_indices_t, connected_indices_hgtd, connected_indices_truth_t;

    double del_t = std::abs(pri_vtx_t - tru_vtx_t);

    // HGTD ntuple vertex time (for comparison with clustering)
    double hgtd_vtx_t = branch.recoVtxTime[0];
    double hgtd_vtx_t_var = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
    double del_t_hgtd = std::abs(hgtd_vtx_t - tru_vtx_t);

    // Track clustering statistics
    if (!hgtd_track_indices.empty()) {
      n_events_clustering_attempted++;
    }

    // if (del_t > t_cut) clustering_succeeded = false;

    if (clustering_succeeded) {
      n_events_clustering_succeeded++;
    }

    for(int track_idx = 0; track_idx < branch.trackEta.GetSize(); track_idx++) {
      double this_pt = branch.trackPt[track_idx];

      if (this_pt < 1.0)
        continue;

      if (std::abs(branch.trackEta[track_idx]) < 2.38)
        continue;

      if (std::abs(branch.trackEta[track_idx]) > 4.0)
        continue;

      if (not(branch.trackQuality[track_idx] == true))
	continue;

      double this_z0 = branch.trackZ0[track_idx];
      double this_z0_var = branch.trackVarZ0[track_idx];

      // z cut
      float dz = pri_vtx_z - this_z0;
      float z_nsigma = std::abs(dz / std::sqrt(this_z0_var));
      bool z_sig_cut = z_nsigma < max_nsigma;

      if (z_sig_cut)
	connected_indices_z.push_back(track_idx);
    }

    // filter connected_indices_z by clustering time cut (for Cone Clustering case)
    for (auto track_idx: connected_indices_z) {
      bool apply_t_cut = branch.trackTimeValid[track_idx] == 1 && clustering_succeeded;

      if (!apply_t_cut) {
	connected_indices_t.push_back(track_idx);
	continue;
      }

      // time cut using clustering vertex time
      double this_t = branch.trackTime[track_idx];
      double this_t_var = branch.trackTimeRes[track_idx]*branch.trackTimeRes[track_idx];
      float dt = this_t - pri_vtx_t;
      float t_nsigma = std::abs(dt / std::sqrt(this_t_var));
      bool t_sig_cut = t_nsigma < max_nsigma;

      if (t_sig_cut)
	connected_indices_t.push_back(track_idx);
    }

    // filter connected_indices_z by HGTD ntuple time cut (for HGTD comparison)
    for (auto track_idx: connected_indices_z) {
      bool apply_t_cut = branch.trackTimeValid[track_idx] == 1 && branch.recoVtxValid[0] == 1;

      if (!apply_t_cut) {
	connected_indices_hgtd.push_back(track_idx);        
	continue;
      }

      // time cut using HGTD ntuple vertex time
      double this_t = branch.trackTime[track_idx];
      double this_t_var = branch.trackTimeRes[track_idx]*branch.trackTimeRes[track_idx];
      float dt = this_t - hgtd_vtx_t;
      float t_nsigma = std::abs(dt / std::sqrt(this_t_var));
      bool t_sig_cut = t_nsigma < max_nsigma;

      if (t_sig_cut)
	connected_indices_hgtd.push_back(track_idx);
    }

    // filter connected_indices_z by truth time cut (for truth comparison)
    for (auto track_idx: connected_indices_z) {
      bool apply_t_cut = branch.trackTimeValid[track_idx] == 1;

      if (!apply_t_cut) {
	connected_indices_truth_t.push_back(track_idx);
	continue;
      }

      // time cut using truth vertex time
      double this_t = branch.trackTime[track_idx];
      double this_t_var = branch.trackTimeRes[track_idx]*branch.trackTimeRes[track_idx];
      float dt = this_t - tru_vtx_t;
      float t_nsigma = std::abs(dt / std::sqrt(this_t_var));
      bool t_sig_cut = t_nsigma < max_nsigma;

      if (t_sig_cut)
	connected_indices_truth_t.push_back(track_idx);
    }

    // calculate jet rpts
    for (int jet_idx = 0; jet_idx < topojet_pt.GetSize(); ++jet_idx) {
      double jet_pt  = topojet_pt [jet_idx];
      double jet_eta = topojet_eta[jet_idx];
      double jet_phi = topojet_phi[jet_idx];

      if (30 < jet_pt && jet_pt < 50)
	continue;

      if (std::abs(jet_eta) > 4.0 || std::abs(jet_eta) < 2.4)
	continue;
      
      // find out if it is a hard scatter jet
      bool jet_isHS = topojet_hsjet_indices[jet_idx].size() > 0;
      // for (int tj_idx = 0; tj_idx < truthhsjet_eta.GetSize(); ++tj_idx) {
      // 	double truthjet_pt  = truthhsjet_pt [tj_idx];
      // 	double truthjet_eta = truthhsjet_eta[tj_idx];
      // 	double truthjet_phi = truthhsjet_phi[tj_idx];

      // 	if (truthjet_pt < 10)
      // 	  continue;

      // 	double deta = truthjet_eta - jet_eta;
      // 	double dphi = TVector2::Phi_mpi_pi(truthjet_phi - jet_phi); // handles phi wrapping!
      // 	double dR_truth_reco = std::sqrt(deta*deta + dphi*dphi);

      // 	if (dR_truth_reco < 0.3)
      // 	  jet_isHS = true;
      // }

      // if (jet_isHS != (topojet_hsjet_indices[jet_idx].size() > 0))
      // 	std::cout << "!!!!!!!!!!!!!!COMPUTED HS: " <<
      // 	  "computed: " << jet_isHS << " " <<
      // 	  "ntuple: " << (topojet_hsjet_indices[jet_idx].size() > 0) << " " <<
      // 	  "vector size: " << topojet_hsjet_indices[jet_idx].size() << " " <<
      // 	  "1st element: " << topojet_hsjet_indices[jet_idx][0] <<
      // 	  std::endl;

      // check if it is a pileup jet
      bool jet_isPU = topojet_itpujet_indices[jet_idx].size()+topojet_otpujet_indices[jet_idx].size() > 0;
      // if (!jet_isHS) {
      // 	for (int tj_idx = 0; tj_idx < truthhsjet_eta.GetSize(); ++tj_idx) {
      // 	  double truthjet_pt  = truthhsjet_pt [tj_idx];
      // 	  double truthjet_eta = truthhsjet_eta[tj_idx];
      // 	  double truthjet_phi = truthhsjet_phi[tj_idx];

      // 	  if (truthjet_pt < 4)
      // 	    continue;

      // 	  double deta = truthjet_eta - jet_eta;
      // 	  double dphi = TVector2::Phi_mpi_pi(truthjet_phi - jet_phi); // handles phi wrapping!
      // 	  double dR_truth_reco = std::sqrt(deta*deta + dphi*dphi);

      // 	  if (dR_truth_reco < 0.6) {
      // 	    jet_isPU = false;
      // 	    break;
      //     }
      // 	}
      // }

      // check if it matches status of ntuple indices
      // if (jet_isPU != (topojet_itpujet_indices[jet_idx].size()+topojet_otpujet_indices[jet_idx].size() > 0))
      // 	std::cout << "COMPUTED PU MISMATCH " <<
      // 	  "computed: " << jet_isPU << " " <<
      // 	  "ntuple: " << (topojet_itpujet_indices[jet_idx].size()+topojet_otpujet_indices[jet_idx].size() > 0) << " " <<
      // 	  "vector sizes: " << topojet_itpujet_indices[jet_idx].size() << " " << topojet_otpujet_indices[jet_idx].size() << " " <<
      // 	  "truthHSjet size: " << topojet_hsjet_indices[jet_idx].size() << 
      // 	  std::endl;

      double z_track_pt_sum = 0; // numerator of RpT
      for (auto track_idx: connected_indices_z) {
	double this_track_pt  = branch.trackPt [track_idx];
	double this_track_eta = branch.trackEta[track_idx];
	double this_track_phi = branch.trackPhi[track_idx];

	double deta = (jet_eta - this_track_eta);
	double dphi = TVector2::Phi_mpi_pi(jet_phi - this_track_phi);
	double dR_track_jet = std::sqrt(deta*deta+dphi*dphi);

	if(dR_track_jet > 0.3) continue;
	// all cuts applied
	z_track_pt_sum += this_track_pt;
	// maybe add flag for HS?
      }

      double jet_RpT_z = z_track_pt_sum / jet_pt;

      if (rpt_min > jet_RpT_z || jet_RpT_z > rpt_max)
	std::cout << "Z RPT OOB: " << jet_RpT_z << std::endl;

      if (jet_isHS) {
	HS_RpT_z_hist->Fill(jet_RpT_z);
	if (del_t < t_cut && branch.recoVtxValid[0] == 1)
	  HS_idl_z_hist->Fill(jet_RpT_z);
      } else if (jet_isPU) {
	PU_RpT_z_hist->Fill(jet_RpT_z);
	if (del_t < t_cut && branch.recoVtxValid[0] == 1)
	  PU_idl_z_hist->Fill(jet_RpT_z);
      }

      double t_track_pt_sum = 0; // numerator of RpT
      for (auto track_idx: connected_indices_t) {
	double this_track_pt  = branch.trackPt [track_idx];
	double this_track_eta = branch.trackEta[track_idx];
	double this_track_phi = branch.trackPhi[track_idx];

	double deta = (jet_eta - this_track_eta);
	double dphi = TVector2::Phi_mpi_pi(jet_phi - this_track_phi); // handles phi wrapping!
	double dR_track_jet = std::sqrt(deta*deta+dphi*dphi);

	if(dR_track_jet > 0.3) continue;
	t_track_pt_sum += this_track_pt;
      }

      double jet_RpT_t = t_track_pt_sum / jet_pt;
      if (rpt_min > jet_RpT_t || jet_RpT_t > rpt_max)
	std::cout << "T RPT OOB: " << jet_RpT_t << std::endl;
      if (diff_min > (jet_RpT_z-jet_RpT_t) || (jet_RpT_z-jet_RpT_t) > diff_max)
	std::cout << "DIFF OOB: " << jet_RpT_z-jet_RpT_t << std::endl;

      if (jet_isHS) {
	HS_RpT_t_hist->Fill(jet_RpT_t);
	if (del_t < t_cut && branch.recoVtxValid[0] == 1)
	  HS_idl_t_hist->Fill(jet_RpT_t);
	// else
	  // std::cout << Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0", filename.substr(40,6).c_str(),this_evnt) << std::endl;
	// if (jet_RpT_t < 0.15)
	  // RpT_diff_020->Fill(jet_RpT_z-jet_RpT_t);
	// if (jet_RpT_t < 0.02) // low rpt, but tagged hard-scatter
	  // fail_pu.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0", filename.substr(40,6).c_str(),this_evnt));
      } else if (jet_isPU) {
	PU_RpT_t_hist->Fill(jet_RpT_t);
	if (del_t < t_cut && branch.recoVtxValid[0] == 1) // only check dt if theres a valid time
	  PU_idl_t_hist->Fill(jet_RpT_t);

	if (std::abs(jet_RpT_z-jet_RpT_t) > 1e-6) {
	  RpT_diff_all->Fill(jet_RpT_z-jet_RpT_t);
	  if (jet_RpT_t > 0.15)
	    RpT_diff_015->Fill(jet_RpT_z-jet_RpT_t);
	  if (jet_RpT_t > 0.2)
	    RpT_diff_020->Fill(jet_RpT_z-jet_RpT_t);
	}
	// if (jet_RpT_t > 0.5) // high rpt, but tagged pile-up
	  // fail_pu.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0", filename.substr(40,6).c_str(),this_evnt));
      }

      double hgtd_track_pt_sum = 0; // numerator of RpT (HGTD ntuple)
      for (auto track_idx: connected_indices_hgtd) {
	double this_track_pt  = branch.trackPt [track_idx];
	double this_track_eta = branch.trackEta[track_idx];
	double this_track_phi = branch.trackPhi[track_idx];

	double deta = (jet_eta - this_track_eta);
	double dphi = TVector2::Phi_mpi_pi(jet_phi - this_track_phi); // handles phi wrapping!
	double dR_track_jet = std::sqrt(deta*deta+dphi*dphi);

	if(dR_track_jet > 0.3) continue;
	hgtd_track_pt_sum += this_track_pt;
      }

      double jet_RpT_hgtd = hgtd_track_pt_sum / jet_pt;
      if (rpt_min > jet_RpT_hgtd || jet_RpT_hgtd > rpt_max)
	std::cout << "HGTD RPT OOB: " << jet_RpT_hgtd << std::endl;

      if (jet_isHS) {
	HS_RpT_hgtd_hist->Fill(jet_RpT_hgtd);
      } else if (jet_isPU) {
	PU_RpT_hgtd_hist->Fill(jet_RpT_hgtd);
      }

      double truth_t_track_pt_sum = 0; // numerator of RpT
      for (auto track_idx: connected_indices_truth_t) {
	double this_track_pt  = branch.trackPt [track_idx];
	double this_track_eta = branch.trackEta[track_idx];
	double this_track_phi = branch.trackPhi[track_idx];

	double deta = (jet_eta - this_track_eta);
	double dphi = TVector2::Phi_mpi_pi(jet_phi - this_track_phi); // handles phi wrapping!
	double dR_track_jet = std::sqrt(deta*deta+dphi*dphi);

	if(dR_track_jet > 0.3) continue;
	truth_t_track_pt_sum += this_track_pt;
      }

      double jet_RpT_truth_t = truth_t_track_pt_sum / jet_pt;
      if (rpt_min > jet_RpT_truth_t || jet_RpT_truth_t > rpt_max)
	std::cout << "TRUTH T RPT OOB: " << jet_RpT_truth_t << std::endl;

      if (jet_isHS) {
	HS_RpT_truth_hist->Fill(jet_RpT_truth_t);
      } else if (jet_isPU) {
	PU_RpT_truth_hist->Fill(jet_RpT_truth_t);
      }

      // Inside the jet loop, after RpT calculations
      if (jet_isHS) {
	// Extract file number from path (e.g., ntuple_1.root -> 1)
	TString fileName = branch.reader.GetTree()->GetCurrentFile()->GetName();
	TString file_num = fileName(49, 6);
	file_num.ReplaceAll("ntuple_", "");
	file_num.ReplaceAll(".root", "");

    
	// Case 2: RPT without time fails AND adding time does NOT fix it
	if (jet_RpT_z < rpt_cut && jet_RpT_t < rpt_cut) {
	  examples_not_fixed_by_time.push_back({
	      file_num, this_evnt, pri_vtx_t, (rpt_cut - jet_RpT_t)
	    });
	}
      }

      if (jet_isPU) {
        TString fileName = branch.reader.GetTree()->GetCurrentFile()->GetName();
	TString file_num = fileName(49, 6);
	file_num.ReplaceAll("ntuple_", "");
	file_num.ReplaceAll(".root", "");
	// Case 1: RPT without time fails (< cut) but timing saves it (> cut)
	if (jet_RpT_z > rpt_cut && jet_RpT_t < rpt_cut) {
	  examples_rescued_by_time.push_back({
	      file_num, this_evnt, pri_vtx_t, (jet_RpT_t - jet_RpT_z)
	    });
	}
      }
      
    }
  }

  // Generate ROC(?)
  
  TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);

  for (auto hist: {PU_RpT_t_hist, PU_RpT_hgtd_hist, PU_RpT_z_hist, HS_RpT_t_hist, HS_RpT_hgtd_hist, HS_RpT_z_hist, PU_idl_t_hist, PU_idl_z_hist, HS_idl_t_hist, HS_idl_z_hist,PU_RpT_truth_hist,HS_RpT_truth_hist}) {
    hist->GetXaxis()->SetRangeUser(0.0, 1.5);
    hist->GetXaxis()->SetNdivisions(515); // min-max/0.1
  }

  PU_RpT_t_hist->SetLineColor(c01);
  PU_RpT_t_hist->SetLineWidth(2);
  legend->AddEntry(PU_RpT_t_hist,"Pile-Up (z+t Cone Clustering)");
  HS_RpT_t_hist->SetLineColor(c02);
  HS_RpT_t_hist->SetLineWidth(2);
  legend->AddEntry(HS_RpT_t_hist,"Hard Scatter (z+t Cone Clustering)");

  PU_RpT_hgtd_hist->SetLineColor(c03);
  PU_RpT_hgtd_hist->SetLineWidth(2);
  legend->AddEntry(PU_RpT_hgtd_hist,"Pile-Up (z+t HGTD)");
  HS_RpT_hgtd_hist->SetLineColor(c04);
  HS_RpT_hgtd_hist->SetLineWidth(2);
  legend->AddEntry(HS_RpT_hgtd_hist,"Hard Scatter (z+t HGTD)");

  PU_RpT_z_hist->SetLineColor(c05);
  PU_RpT_z_hist->SetLineWidth(2);
  legend->AddEntry(PU_RpT_z_hist,"Pile-Up (z)");
  HS_RpT_z_hist->SetLineColor(c06);
  HS_RpT_z_hist->SetLineWidth(2);
  legend->AddEntry(HS_RpT_z_hist,"Hard Scatter (z)");

  PU_idl_t_hist->SetLineColor(c07);
  PU_idl_t_hist->SetLineWidth(2);
  legend->AddEntry(PU_idl_t_hist,"Pile-Up (z+t) (Ideal Timing)");
  HS_idl_t_hist->SetLineColor(c08);
  HS_idl_t_hist->SetLineWidth(2);
  legend->AddEntry(HS_idl_t_hist,"Hard Scatter (z+t) (Ideal Timing)");

  PU_idl_z_hist->SetLineColor(c09);
  PU_idl_z_hist->SetLineWidth(2);
  legend->AddEntry(PU_idl_z_hist,"Pile-Up (z) (Ideal Timing)");
  HS_idl_z_hist->SetLineColor(c10);
  HS_idl_z_hist->SetLineWidth(2);
  legend->AddEntry(HS_idl_z_hist,"Hard Scatter (z) (Ideal Timing)");

  PU_RpT_truth_hist->SetLineColor(c11);
  PU_RpT_truth_hist->SetLineWidth(2);
  HS_RpT_truth_hist->SetLineColor(c01);
  HS_RpT_truth_hist->SetLineWidth(2);
  
  canvas->SetLogy(true);
  canvas->Print("../figs/rpt_plots.pdf[");

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);  // Align left-top

  const char* cuttext = "#geq 1 Truth Matched Forward Jet";

  double old_pu_t_max = PU_RpT_t_hist->GetMaximum(),
    old_hs_t_max = HS_RpT_t_hist->GetMaximum(),
    old_pu_z_max = PU_RpT_z_hist->GetMaximum(),
    old_hs_z_max = HS_RpT_z_hist->GetMaximum(),
    old_hs_t_idl_max = HS_idl_t_hist->GetMaximum(),
    old_pu_t_idl_max = PU_idl_t_hist->GetMaximum(),
    old_hs_z_idl_max = HS_idl_z_hist->GetMaximum(),
    old_pu_z_idl_max = PU_idl_z_hist->GetMaximum();

  // PU/HS Comparison z
  TLegend *z_legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  z_legend->AddEntry(PU_RpT_z_hist, "Pile-Up");
  z_legend->AddEntry(HS_RpT_z_hist, "Hard-Scatter");
  
  PU_RpT_z_hist->SetMaximum(5*std::max({old_pu_z_max,old_hs_z_max}));
  PU_RpT_z_hist->Draw("HIST");
  HS_RpT_z_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  PU_RpT_z_hist->SetTitle("RpT (z)");
  z_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");
  
  // PU/HS Comparison z+time
  TLegend *zt_legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  zt_legend->AddEntry(PU_RpT_t_hist, "Pile-Up");
  zt_legend->AddEntry(HS_RpT_t_hist, "Hard-Scatter");
  
  PU_RpT_t_hist->SetMaximum(5*std::max({old_pu_t_max,old_hs_t_max}));
  PU_RpT_t_hist->Draw("HIST");
  HS_RpT_t_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  PU_RpT_t_hist->SetTitle("RpT (z+time)");
  zt_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  // PU/HS Comparison z (IDEAL TIMING)
  TLegend *idealz_legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  idealz_legend->AddEntry(PU_idl_z_hist, "Pile-Up");
  idealz_legend->AddEntry(HS_idl_z_hist, "Hard-Scatter");
  
  PU_idl_z_hist->SetMaximum(5*std::max({old_pu_z_idl_max,old_hs_z_idl_max}));
  PU_idl_z_hist->Draw("HIST");
  HS_idl_z_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  PU_idl_z_hist->SetTitle("RpT (z) Ideal Timing");
  idealz_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");
  
  // PU/HS Comparison z+t (IDEAL TIMING)
  TLegend *idealtime_legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  idealtime_legend->AddEntry(PU_idl_t_hist, "Pile-Up");
  idealtime_legend->AddEntry(HS_idl_t_hist, "Hard-Scatter");
  
  PU_idl_t_hist->SetMaximum(5*std::max({old_pu_t_idl_max,old_hs_t_idl_max}));
  PU_idl_t_hist->Draw("HIST");
  HS_idl_t_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  PU_idl_t_hist->SetTitle("RpT (z+time) Ideal Timing");
  idealtime_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  // PU/HS Comparison z+t (truth time)
  TLegend *truthtime_legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  truthtime_legend->AddEntry(PU_RpT_truth_hist, "Pile-Up");
  truthtime_legend->AddEntry(HS_RpT_truth_hist, "Hard-Scatter");
  
  PU_RpT_truth_hist->SetMaximum(5*std::max({PU_RpT_truth_hist->GetMaximum(),HS_RpT_truth_hist->GetMaximum()}));
  PU_RpT_truth_hist->Draw("HIST");
  HS_RpT_truth_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  PU_RpT_truth_hist->SetTitle("RpT (z+time) Truth Vertex Time");
  truthtime_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  // z/z+time/ideal time Comparison (PU)
  PU_RpT_z_hist->SetMaximum(5*std::max({old_pu_z_max,old_pu_t_max,old_pu_z_idl_max,old_pu_t_idl_max}));
  PU_RpT_z_hist->Draw("HIST");
  PU_RpT_t_hist->Draw("HIST SAME");
  PU_idl_z_hist->Draw("HIST SAME");
  PU_idl_t_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  PU_RpT_z_hist->SetTitle("Pile-Up Comparison");
  legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  // z/z+time Comparison (HS)
  HS_RpT_z_hist->SetMaximum(5*std::max({old_hs_z_max,old_hs_t_max,old_hs_z_idl_max,old_hs_t_idl_max}));
  HS_RpT_z_hist->Draw("HIST");
  HS_RpT_t_hist->Draw("HIST SAME");
  HS_idl_z_hist->Draw("HIST SAME");
  HS_idl_t_hist->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, cuttext);
  HS_RpT_z_hist->SetTitle("Hard Scatter Comparison");
  legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  // Difference Plot
  canvas->SetLogy(false);
  TLegend *diff_legend = new TLegend(0.6,0.7,0.9,0.9);
  
  diff_legend->AddEntry(RpT_diff_020,"Pile-Up, R_{pT}>0.20");
  diff_legend->AddEntry(RpT_diff_015,"Pile-Up, R_{pT}>0.15");
  diff_legend->AddEntry(RpT_diff_all,"Pile-Up, All R_{pT}");

  RpT_diff_015->Scale(1/RpT_diff_015->Integral());
  RpT_diff_020->Scale(1/RpT_diff_020->Integral());
  RpT_diff_all->Scale(1/RpT_diff_all->Integral());
  double max_val = 1.1*std::max({RpT_diff_020->GetMaximum(), RpT_diff_015->GetMaximum(), RpT_diff_all->GetMaximum()});
    
  RpT_diff_020->SetTitle("R_{pT}(z)-R_{pT}(z,t)");
  RpT_diff_020->SetLineColor(c01);
  RpT_diff_020->SetLineWidth(2);
  RpT_diff_020->GetXaxis()->SetRangeUser(-0.01, 0.5);
  RpT_diff_020->SetMaximum(max_val);

  RpT_diff_015->SetTitle("R_{pT}(z)-R_{pT}(z,t)");
  RpT_diff_015->SetLineColor(c02);
  RpT_diff_015->SetLineWidth(2);
  RpT_diff_015->GetXaxis()->SetRangeUser(-0.01, 0.5);
  RpT_diff_015->SetMaximum(max_val);

  RpT_diff_all->SetTitle("R_{pT}(z)-R_{pT}(z,t)");
  RpT_diff_all->SetLineColor(c03);
  RpT_diff_all->SetLineWidth(2);
  RpT_diff_all->GetXaxis()->SetRangeUser(-0.01, 0.5);
  RpT_diff_all->SetMaximum(max_val);

  RpT_diff_020->Draw("HIST");
  RpT_diff_015->Draw("HIST SAME");
  RpT_diff_all->Draw("HIST SAME");
  diff_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  // roc curve(s)
  TGraph* z_graph = generate_roc(PU_RpT_z_hist, HS_RpT_z_hist);
  TGraph* t_graph = generate_roc(PU_RpT_t_hist, HS_RpT_t_hist);
  TGraph* hgtd_graph = generate_roc(PU_RpT_hgtd_hist, HS_RpT_hgtd_hist);
  TGraph* z_ideal = generate_roc(PU_idl_z_hist, HS_idl_z_hist);
  TGraph* t_ideal = generate_roc(PU_idl_t_hist, HS_idl_t_hist);
  // TGraph* t_truth = generate_roc(PU_RpT_truth_hist, HS_RpT_truth_hist);

  std::cout << "\n=== ROC CURVE STATISTICS ===" << std::endl;
  std::cout << "z_graph points       : " << z_graph->GetN() << std::endl;
  std::cout << "t_graph points       : " << t_graph->GetN() << std::endl;
  std::cout << "hgtd_graph points    : " << hgtd_graph->GetN() << std::endl;
  std::cout << "t_ideal points       : " << t_ideal->GetN() << std::endl;

  TLegend *roc_legend = new TLegend(0.6, 0.6, 0.9, 0.9);

  roc_legend->AddEntry(z_graph,"Z Info Only");
  roc_legend->AddEntry(t_graph,"Add MY t_{0}");
  roc_legend->AddEntry(hgtd_graph,"HGTD t_{0}");
  // roc_legend->AddEntry(z_ideal,"z ideal timing");
  roc_legend->AddEntry(t_ideal,"Ideal t_{0}");
  // roc_legend->AddEntry(t_truth,"z+t, truth vtx time");
  roc_legend->SetTextSize(.05);  

  auto roc_xmin = 0.8, roc_xmax = 1.0;

  z_graph->SetTitle("RpT Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
  z_graph->SetLineColor(c01);
  z_graph->SetLineWidth(2);
  z_graph->GetXaxis()->SetRangeUser(roc_xmin, roc_xmax);
  z_graph->GetXaxis()->SetNdivisions(810);
  z_graph->SetMaximum(180);

  t_graph->SetTitle("RpT Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
  t_graph->SetLineColor(c02);
  t_graph->SetLineWidth(2);
  t_graph->GetXaxis()->SetRangeUser(roc_xmin, roc_xmax);
  t_graph->GetXaxis()->SetNdivisions(810);
  t_graph->SetMaximum(180);

  hgtd_graph->SetTitle("RpT Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
  hgtd_graph->SetLineColor(c03);
  hgtd_graph->SetLineWidth(2);
  hgtd_graph->GetXaxis()->SetRangeUser(roc_xmin, roc_xmax);
  hgtd_graph->GetXaxis()->SetNdivisions(810);
  hgtd_graph->SetMaximum(180);

  // z_ideal->SetTitle("RpT Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
  // z_ideal->SetLineColor(c03);
  // z_ideal->SetLineWidth(2);
  // z_ideal->GetXaxis()->SetRangeUser(roc_xmin, roc_xmax);
  // z_ideal->GetXaxis()->SetNdivisions(810);
  // z_ideal->SetMaximum(180);

  t_ideal->SetTitle("RpT Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
  t_ideal->SetLineColor(c04);
  t_ideal->SetLineWidth(2);
  t_ideal->GetXaxis()->SetRangeUser(roc_xmin, roc_xmax);
  t_ideal->GetXaxis()->SetNdivisions(810);
  t_ideal->SetMaximum(180);

  // t_truth->SetTitle("RpT Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
  // t_truth->SetLineColor(c05);
  // t_truth->SetLineWidth(2);
  // t_truth->GetXaxis()->SetRangeUser(roc_xmin, roc_xmax);
  // t_truth->GetXaxis()->SetNdivisions(810);
  // t_truth->SetMaximum(120);
  
  z_graph->Draw("ALP");
  t_graph->Draw("LP SAME");
  hgtd_graph->Draw("LP SAME");
  // z_ideal->Draw("LP SAME");
  t_ideal->Draw("LP SAME");
  // t_truth->Draw("LP SAME");
  roc_legend->Draw();
  canvas->Print("../figs/rpt_plots.pdf");

  canvas->Print("../figs/rpt_plots.pdf]");

  std::cout << "\n=== CLUSTERING STATISTICS ===" << std::endl;
  std::cout << "Events with HGTD tracks           : " << n_events_clustering_attempted << std::endl;
  std::cout << "Events with successful clustering : " << n_events_clustering_succeeded << std::endl;
  if (n_events_clustering_attempted > 0) {
    std::cout << "Clustering success rate           : "
              << (100.0 * n_events_clustering_succeeded / n_events_clustering_attempted)
              << "%" << std::endl;
  }

  std::cout << "\n=== HISTOGRAM STATISTICS ===" << std::endl;
  std::cout << "Number of Pile-Up Jets (z)       : " << PU_RpT_z_hist->Integral() << std::endl;
  std::cout << "Number of HS Jets (z)            : " << HS_RpT_z_hist->Integral() << std::endl;
  std::cout << "Number of Pile-Up Jets (z+t Cone): " << PU_RpT_t_hist->Integral() << std::endl;
  std::cout << "Number of HS Jets (z+t Cone)     : " << HS_RpT_t_hist->Integral() << std::endl;
  std::cout << "Number of Pile-Up Jets (z+t HGTD): " << PU_RpT_hgtd_hist->Integral() << std::endl;
  std::cout << "Number of HS Jets (z+t HGTD)     : " << HS_RpT_hgtd_hist->Integral() << std::endl;

  // Check if histograms are different
  double diff_count = 0;
  for (int i = 1; i <= PU_RpT_t_hist->GetNbinsX(); i++) {
    if (std::abs(PU_RpT_t_hist->GetBinContent(i) - PU_RpT_hgtd_hist->GetBinContent(i)) > 0.01) {
      diff_count++;
    }
  }
  std::cout << "\nDifferent bins (PU Cone vs HGTD) : " << diff_count << " / " << PU_RpT_t_hist->GetNbinsX() << std::endl;

  diff_count = 0;
  for (int i = 1; i <= HS_RpT_t_hist->GetNbinsX(); i++) {
    if (std::abs(HS_RpT_t_hist->GetBinContent(i) - HS_RpT_hgtd_hist->GetBinContent(i)) > 0.01) {
      diff_count++;
    }
  }
  std::cout << "Different bins (HS Cone vs HGTD)  : " << diff_count << " / " << HS_RpT_t_hist->GetNbinsX() << std::endl;
  // std::cout << "Number of total jets        : " << njets << std::endl;
  // std::cout << "Number of total hs jets     : " << nhsjets << std::endl;
  // std::cout << "Number of total event       : " << nevent << std::endl;
  // std::cout << "Number of avg jet/event     : " << njets/nevent << std::endl;

  // std::cout << "z_x = [";
  // for (auto x: z_vector_x)
  //   std::cout << x << ", ";
  // std::cout << "]" << std::endl;

  // std::cout << "z_y = [";
  // for (auto x: z_vector_y)
  //   std::cout << x << ", ";
  // std::cout << "]" << std::endl;

  // std::cout << "t_x = [";
  // for (auto x: t_vector_x)
  //   std::cout << x << ", ";
  // std::cout << "]" << std::endl;

  // std::cout << "t_y = [";
  // for (auto x: t_vector_y)
  //   std::cout << x << ", ";
  // std::cout << "]" << std::endl;

  // std::cout << "------------------------------" << std::endl << "Failing HS"
  // << std::endl; for(int i = 0; i < 20; ++i) {
  //   int idx = rand() % fail_hs.size();
  //   std:: cout << fail_hs[idx] << std::endl;
  // }
  // std::cout << "Total number of fail HS: " << fail_hs.size() << std::endl;
  // std::cout << "------------------------------" << std::endl << "Failing PU"
  // << std::endl; for(auto cmd: fail_pu)
  //   std:: cout << cmd << std::endl;
  // std::cout << "Total number of fail PU: " << fail_pu.size() << std::endl;

  // Add at the very end of main()
  std::cout << "\n========================================" << std::endl;
  std::cout << " EXAMPLE EVENTS FOR EVENT DISPLAY" << std::endl;
  std::cout << "========================================\n" << std::endl;

  auto print_examples = [](const char* title, std::vector<ExampleEvent>& events, int max_print = 10) {
    std::cout << title << ":\n";
    if (events.empty()) {
      std::cout << " (No examples found)\n";
      return;
    }
    // Sort by metric descending to show most impactful cases first
    std::sort(events.begin(), events.end(), [](const ExampleEvent& a, const ExampleEvent& b) {
      return a.metric > b.metric;
    });

    int n_print = std::min(max_print, (int)events.size());
    for (int i = 0; i < n_print; i++) {
      printf("python3 event_display.py --file_num %s --event_num %lld --extra_time %.2f\n", 
	     events[i].file_num.Data(), events[i].event_num, events[i].extra_time);
    }
    std::cout << std::endl;
  };

  print_examples("CASE 1: Rescued by Timing (RpT_z < cut, RpT_t > cut)", examples_rescued_by_time);
  print_examples("CASE 2: Still Fails with Timing (RpT_z < cut, RpT_t < cut)", examples_not_fixed_by_time);  

  return 0;
}
