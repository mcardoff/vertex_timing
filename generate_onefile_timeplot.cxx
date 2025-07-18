#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <THStack.h>
#include <TVector2.h>
#include <cmath>
#include <iostream>
#include <boost/filesystem.hpp>

// Explicitly load the dictionary for ROOT to recognize nested vectors
#ifdef __CLING__
#pragma link C++ class std::vector<std::vector<int>>+;
#endif

#define debug false

using std::endl;

TTree *setup_tree(std::string Number) {
  TFile *file = TFile::Open(Form("../ntuple/user.scheong.42871997.Output._%s.SuperNtuple.root", Number.c_str()));
  TTree *tree = (TTree*)file->Get("ntuple");
  return tree;
}

void generate_onefile_timeplot(std::string Number, int event_num) {
  gStyle->SetOptStat(0);
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  // colors
  auto c1 = TColor::GetColor("#3f90da");
  auto c2 = TColor::GetColor("#ffa90e");
  auto c3 = TColor::GetColor("#bd1f01");
  auto c4 = TColor::GetColor("#94a4a2");
  auto c5 = TColor::GetColor("#832db6");

  TCanvas *canvas = new TCanvas("canvas", "Track Time Histogram", 800, 600);
  TTree * tree = setup_tree(Number);
  
  // Jet variables
  std::vector<float> *jet_pt = nullptr;
  std::vector<float> *jet_phi = nullptr;
  std::vector<float> *jet_eta = nullptr;
  std::vector<std::vector<int>> *jet_track_indices = nullptr;

  // Vertex variables
  std::vector<float> *truth_vtx_z = nullptr;
  std::vector<float> *truth_vtx_time = nullptr;
  std::vector<float> *reco_vtx_z = nullptr;
  std::vector<float> *reco_vtx_time = nullptr;
  std::vector<float> *reco_vtx_time_res = nullptr;
  std::vector<int>   *reco_vtx_valid = nullptr;

  // Track variables
  std::vector<float> *track_z0 = nullptr;
  std::vector<float> *track_pt = nullptr;
  std::vector<float> *track_phi = nullptr;
  std::vector<float> *track_eta = nullptr;
  std::vector<float> *track_time = nullptr;
  std::vector<float> *track_var_z0 = nullptr;
  std::vector<int>   *track_truthVtx_idx = nullptr;
  std::vector<int>   *track_time_valid = nullptr;

  // Set branch addresses
  tree->SetBranchAddress("AntiKt4EMTopoJets_pt", &jet_pt);
  tree->SetBranchAddress("AntiKt4EMTopoJets_phi", &jet_phi);
  tree->SetBranchAddress("AntiKt4EMTopoJets_eta", &jet_eta);
  tree->SetBranchAddress("AntiKt4EMTopoJets_track_idx", &jet_track_indices);
  
  tree->SetBranchAddress("TruthVtx_z", &truth_vtx_z);
  tree->SetBranchAddress("TruthVtx_time", &truth_vtx_time);
  
  tree->SetBranchAddress("RecoVtx_z", &reco_vtx_z);
  tree->SetBranchAddress("RecoVtx_time", &reco_vtx_time);
  tree->SetBranchAddress("RecoVtx_timeRes", &reco_vtx_time_res);
  tree->SetBranchAddress("RecoVtx_hasValidTime", &reco_vtx_valid);

  tree->SetBranchAddress("Track_z0", &track_z0);
  tree->SetBranchAddress("Track_pt", &track_pt);
  tree->SetBranchAddress("Track_phi", &track_phi);
  tree->SetBranchAddress("Track_eta", &track_eta);
  tree->SetBranchAddress("Track_time", &track_time);
  tree->SetBranchAddress("Track_var_z0", &track_var_z0);
  tree->SetBranchAddress("Track_truthVtx_idx", &track_truthVtx_idx);
  tree->SetBranchAddress("Track_hasValidTime", &track_time_valid);

  int numpages = 1;

  // cut variables
  int    min_jets     = 2;    // min number of jets
  double min_jetpt    = 30.0; // self explanatory
  double min_abs_eta  = 2.0;  // min eta for a vtx to be considered "forward"
  double min_abs_eta_track = 2.4;  // min eta for a track to be considered "forward"
  double min_track_pt = 1.0;  // track_pt > 1.0 GeV
  double max_vtx_dz   = 2.0;  // max error for reco HS vertex z
  double max_nsigma   = 3.0;  // how close a track can be to PV

  tree->GetEntry(event_num);

  float truth_vertex_time = (*truth_vtx_time)[0];
  float reco_vertex_time = (*reco_vtx_time)[0];
  float reco_vertex_time_res = (*reco_vtx_time_res)[0];

  // auto filename = Form("figs/trackhist_%s_%d.pdf",Number.c_str(),event_num);
  auto filename = "EventDisplayTest.pdf";
  
  // select valid times
  std::vector<float> selected_times;
  std::vector<float> weight_pt, weight_dR, weight_jR, weight_tR, weight_sig;
  std::vector<bool> hardScatter_Mask;
  int nvalidtime = 0;
  for (int idx = 0; idx < track_z0->size(); ++idx) {
    if (track_time_valid->at(idx) == 0)
      continue; // skip tracks with invalid time

    if (track_pt->at(idx) < min_track_pt)
      continue;

    float nsigma = (track_z0->at(idx) - reco_vtx_z->at(0))/std::sqrt(track_var_z0->at(idx));
    if(std::abs(nsigma) > max_nsigma)
      continue; // skip tracks too far from vertex

    float min_dR = 10000.0, jptmin_dR = 0;
    for (int jdx = 0; jdx < jet_eta->size(); jdx++) {
      double j_eta = jet_eta->at(jdx),
	j_phi = jet_phi->at(jdx),
	j_pt = jet_pt->at(jdx);

      double deta = track_eta->at(idx) - j_eta;                                               
      double dphi = TVector2::Phi_mpi_pi(track_phi->at(idx) - j_phi); // handles phi wrapping!
      double dR_track_jet = std::sqrt(deta*deta + dphi*dphi);

      if (dR_track_jet < min_dR) {
	min_dR = dR_track_jet;
	jptmin_dR = j_pt*exp(-min_dR);
      }
    }


    std::cout << "(idx, t, pt) " << idx << ", " << track_time->at(idx) << ", " << track_pt->at(idx) << std::endl;
    hardScatter_Mask.push_back(track_truthVtx_idx->at(idx) == 0);
    
    selected_times.push_back(track_time->at(idx));
    weight_pt.push_back(track_pt->at(idx));
    weight_dR.push_back(exp(-min_dR));
    weight_jR.push_back(jptmin_dR);
    weight_tR.push_back(track_pt->at(idx)*exp(-min_dR));
    weight_sig.push_back(exp(-std::abs(nsigma)));

    // if((*track_eta)[idx] > min_abs_eta_track)
    nvalidtime++; // nvalidtime is now number of forward tracks     
  }

  std::cout << nvalidtime << std::endl;

  if (!selected_times.empty()) {
    float min_time = std::min({*std::min_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
    float max_time = std::max({*std::max_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
      
    float range = max_time - min_time;
    float extended_min_time = min_time - 0.05 * range;
    float extended_max_time = max_time + 0.05 * range;
    
    TH1F *hist_weight_one = new TH1F("track_time_hist_event",
				     Form("Track Time Distribution Event %d;Time (ps);A.U.", event_num),
				     50, extended_min_time, extended_max_time);
    TH1F *hist_HS         = (TH1F*)hist_weight_one->Clone("track_time_hist_HS");
    TH1F *hist_PU         = (TH1F*)hist_weight_one->Clone("track_time_hist_PU");
    TH1F *hist_weight_pt  = (TH1F*)hist_weight_one->Clone("track_time_hist_weight_pt");
    TH1F *hist_weight_dR  = (TH1F*)hist_weight_one->Clone("track_time_hist_weight_exp(-minDR)");
    TH1F *hist_weight_jR  = (TH1F*)hist_weight_one->Clone("track_time_hist_weight_jetptexp(-minDR)");
    TH1F *hist_weight_tR  = (TH1F*)hist_weight_one->Clone("track_time_hist_weight_trackptexp(-minDR)");
    TH1F *hist_weight_sig = (TH1F*)hist_weight_one->Clone("track_time_hist_weight_exp(-sig)");
    hist_weight_pt ->SetTitle("Track t Dist. Track p_{T} Weighted");
    hist_weight_dR ->SetTitle("Track t Dist. exp(-#DeltaR_{min}) Weighted");
    hist_weight_jR ->SetTitle("Track t Dist. Jet p_{T}*exp(-#DeltaR_{min}) Weighted");
    hist_weight_tR ->SetTitle("Track t Dist. Track p_{T}*exp(-#DeltaR_{min}) Weighted");
    hist_weight_sig->SetTitle("Track t Dist. exp(-|s|) Weighted");
    
    for (int idx = 0; idx < selected_times.size(); idx++) {
      float time = selected_times[idx];
      float pt_wght = weight_pt[idx], dR_wght=weight_dR[idx], jR_wght=weight_jR[idx], tR_wght=weight_tR[idx], sig_wght=weight_sig[idx];
      hist_weight_one->Fill(time);
      hist_weight_pt->Fill(time,pt_wght);
      hist_weight_dR->Fill(time,dR_wght);
      hist_weight_jR->Fill(time,jR_wght);
      hist_weight_tR->Fill(time,tR_wght);
      hist_weight_sig->Fill(time,sig_wght);

      if (hardScatter_Mask[idx])
	hist_HS->Fill(time);
      else
	hist_PU->Fill(time);
    }

    canvas->Print(Form("%s[",filename), "pdf");
    THStack *stack =  new THStack("stackPUHS", "PU/HS Stack;Time (ps);A.U.");

    hist_HS->SetLineColorAlpha(kBlue, 0.35);
    hist_HS->SetFillColorAlpha(kBlue, 0.35);
    hist_PU->SetLineColorAlpha(kRed,  0.35);
    hist_PU->SetFillColorAlpha(kRed,  0.35);
    
    stack->Add(hist_HS);
    stack->Add(hist_PU);
    stack->SetMaximum(1.25*stack->GetMaximum());
    stack->Draw("HIST");

    TLine *truthHSLine = new TLine(truth_vertex_time, 0, truth_vertex_time, stack->GetMaximum());
    truthHSLine->SetLineColor(kBlue);
    truthHSLine->SetLineWidth(2);
    truthHSLine->Draw("SAME");

    TLine *recoHSLine = new TLine(reco_vertex_time, 0, reco_vertex_time, stack->GetMaximum());
    recoHSLine->SetLineColor(kGreen);
    recoHSLine->SetLineWidth(2);
    recoHSLine->Draw("SAME");

    TLegend *stacklegend = new TLegend(0.65,0.8,0.9,0.9);
    stacklegend->AddEntry(hist_HS,"Hard-Scatter Time Histogram");
    stacklegend->AddEntry(hist_PU,"Pile-Up Time Histogram");
    stacklegend->AddEntry(truthHSLine,"Truth HS Time");
    stacklegend->AddEntry(recoHSLine,"Reco HS Time");
    stacklegend->Draw();


    canvas->Print(filename, "pdf");
      
    for (auto hist: {hist_weight_one,hist_weight_sig,hist_weight_pt,hist_weight_dR,hist_weight_jR,hist_weight_tR}) {
      
      hist->SetFillColorAlpha(kBlue, 0.35);
      hist->Draw("HIST F");

      TLegend *legend = new TLegend(0.65,0.8,0.9,0.9);
      
      TLine *truthHSLine = new TLine(truth_vertex_time, 0, truth_vertex_time, hist->GetMaximum());
      truthHSLine->SetLineColor(kBlue);
      truthHSLine->SetLineWidth(2);
      truthHSLine->Draw("SAME");

      TLine *recoHSLine = new TLine(reco_vertex_time, 0, reco_vertex_time, hist->GetMaximum());
      recoHSLine->SetLineColor(kGreen);
      recoHSLine->SetLineWidth(2);
      recoHSLine->Draw("SAME");

      legend->AddEntry(hist,"Time Histograms");
      legend->AddEntry(truthHSLine,"Truth HS Time");
      legend->AddEntry(recoHSLine,"Reco HS Time");
      legend->Draw();

      TLatex *latex = new TLatex();
      latex->SetNDC();  // Use Normalized Device Coordinates (0 to 1)
      latex->SetTextSize(0.03);
      latex->SetTextAlign(13);  // Align text left
      latex->DrawLatex(0.15, 0.88,
		       Form("Truth Vertex Time: %.2f ps, %d entries",
			    truth_vertex_time, nvalidtime));
      latex->DrawLatex(0.15, 0.85,
		       Form("#Deltat: %.2f ps, #Deltat/res: %.2f", truth_vertex_time-reco_vertex_time, (truth_vertex_time-reco_vertex_time)/reco_vertex_time_res));
      
      hist->SetMaximum(1.25*hist->GetMaximum());

      canvas->Print(filename, "pdf");
      
      delete hist;
      delete latex;
      delete truthHSLine;
      delete recoHSLine;
    }
    canvas->Print(Form("%s]",filename), "pdf");
  }

  std::cout << "python event_display_VBF_R25.py --file_num " << Number
            << " --event_num " << event_num << " --vtxID 0" << endl;
}
