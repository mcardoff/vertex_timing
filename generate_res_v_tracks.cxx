#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TStyle.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <iostream>
#include <boost/filesystem.hpp>

void generate_res_v_tracks() {
  bool debug = false;
  gStyle->SetOptStat(0);
  // colors
  auto c1 = TColor::GetColor("#3f90da");
  auto c2 = TColor::GetColor("#ffa90e");
  auto c3 = TColor::GetColor("#bd1f01");
  auto c4 = TColor::GetColor("#94a4a2");
  auto c5 = TColor::GetColor("#832db6");
  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  TChain chain("ntuple");

  // VBF H->Invisible sample
  for (const auto& entry : boost::filesystem::directory_iterator("./ntuple")) {
    if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
    chain.Add(entry.path().c_str());
  }
  
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found in directory: " << std::endl;
    return;
  }

  TTreeReader reader(&chain);
  // jet variables
  TTreeReaderArray<float> jet_pt
    // (reader, "TruthHSJet_pt");
    (reader, "AntiKt4EMTopoJets_pt");
  TTreeReaderArray<float> jet_eta    
    // (reader, "TruthHSJet_eta");
    (reader, "AntiKt4EMTopoJets_eta");

  // vertex variables
  /// truth vertex
  TTreeReaderArray<float> truth_vtx_z
    (reader, "TruthVtx_z");
  TTreeReaderArray<float> truth_vtx_time
    (reader, "TruthVtx_time");

  /// reco vertex
  TTreeReaderArray<float> reco_vtx_z
    (reader, "RecoVtx_z");
  TTreeReaderArray<float> reco_vtx_time
    (reader, "RecoVtx_time");
  TTreeReaderArray<float> reco_vtx_timeRes
    (reader, "RecoVtx_timeRes");
  TTreeReaderArray<int> reco_vtx_valid
    (reader, "RecoVtx_hasValidTime");

  // track variables
  TTreeReaderArray<float> track_z0
    (reader, "Track_z0");
  TTreeReaderArray<float> track_pt
    (reader, "Track_pt");
  TTreeReaderArray<float> track_eta
    (reader, "Track_eta");
  TTreeReaderArray<float> track_time
    (reader, "Track_time");
  TTreeReaderArray<float> track_time_res
    (reader, "Track_timeRes");
  TTreeReaderArray<float> track_var_z0
    (reader, "Track_var_z0");
  TTreeReaderArray<int> track_to_truthvtx
    (reader, "Track_truthVtx_idx");
  TTreeReaderArray<int> track_to_particle
    (reader, "Track_truthPart_idx");
  TTreeReaderArray<int> track_time_valid
    (reader, "Track_hasValidTime");

  int track_min = 0, track_max = 60;
  double diff_min = -1000.0, diff_max = 1000.0;
  double res_min = -110, res_max = 110;
  int bins_x = (track_max-track_min)/4, bins_y=(diff_max-diff_min)/8;
  
  // Histograms
  /// Store vertex dt = reco_vtx_time[0] - truth_vtx_time[0]
  TH2F *hist1 = new TH2F(
			 "hist1",
			 "RecoVtx t_{0} - TruthVtx t_{0} vs Forward Tracks;n Forward Track;#Delta t[ps]",
			 bins_x,               // bin x
			 track_min, track_max, // bounds x
			 bins_y,               // bin y
			 diff_min, diff_max    // bounds y
			 );

  TH2F *hist2 = new TH2F(
			 "hist2",
			 "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t} vs Forward Tracks;n Forward Track;#Delta t[ps]",
			 bins_x,               // bin x
			 track_min, track_max, // bounds x
			 bins_y,               // bin y
			 res_min, res_max      // bounds y
			 );

  // cut variables
  int    min_jets    = 2;         // min number of jets
  bool   check_valid_time = true; // check whether or not recovtx has a valid time
  double min_jetpt   = 30.0;      // self explanatory
  double min_abs_jet_eta = 2.0;   // min eta for a jet to be considered "forward"
  double min_abs_track_eta = 2.3; // min eta for a track to be considered "forward"
  double min_track_pt = 1.0;      // track_pt > 1.0 GeV
  double max_vtx_dz  = 1.0;       // max error for reco HS vertex z
  double max_nsigma  = 2.0;       // how close a track can be to PV

  while (reader.Next()) {
    if (jet_pt.GetSize() < min_jets ) {
      if(debug)
	std::cout << "Skipping low jet event" << std::endl;
      continue;  // Skip if no jets are present
    }

    if(reco_vtx_valid[0] == 0 && check_valid_time) {
      if(debug)
	std::cout << "Skipping event where reco vertex has no valid time" << std::endl;
      continue;  // Skip if reco vertex has no valid time
    }
    
    // check reco HS vertex is with 2mm of truth HS vertex
    if(std::abs(truth_vtx_z[0] - reco_vtx_z[0]) > max_vtx_dz ) {
      if(debug)
	std::cout << "Skipping event due to incorrect HS vertex" << std::endl;
      continue;
    }
    
    // check if jet pt > 30 GeV
    if (jet_pt.GetSize() != 0) {
      float min_pt = *std::min_element(jet_pt.begin(), jet_pt.end());
      if(min_pt < min_jetpt) {
	if(debug)
	  std::cout << "Skipping event due to low jet pt" << std::endl;
	continue;
      }
    }

    if (reco_vtx_valid[0] != 1) {
      continue;
    }
    
    int nForwardTrack = 0;
    int nForwardJet = 0;
    for(int i = 0; i < track_eta.GetSize(); ++i) {
      auto eta = track_eta[i];
      auto valid = track_time_valid[i];
      auto pt = track_pt[i];
      // if (std::abs(eta) > min_abs_track_eta && valid == 1 && pt > min_track_pt)
      auto nsigma = (track_z0[i] - reco_vtx_z[0])/std::sqrt(track_var_z0[i]);
      if (eta > min_abs_track_eta && pt > min_track_pt && std::abs(nsigma) < max_nsigma)
	nForwardTrack++;
    }
    
    float diff = reco_vtx_time[0] - truth_vtx_time[0];
    float res = diff/reco_vtx_timeRes[0];
    if (nForwardTrack > track_max || nForwardTrack < track_min) 
      std::cout << "!!!!!n_ftrack: " << nForwardTrack << std::endl;
    if (diff > diff_max || diff < diff_min)
      std::cout << "!!!!!diff: " << diff << std::endl;
    if (res > res_max || res < res_min)
      std::cout << "!!!!!res: " << res << std::endl;
    // std::cout << "--------------------" << std::endl;
    // std::cout << "reso: " << diff << std::endl;
    // std::cout << "--------------------" << std::endl;
    hist1->Fill(nForwardTrack, diff); // increments value
    hist2->Fill(nForwardTrack, res); // increments value
  }

  // Fit Gaussian to each slice along Y axis
  TF1 *f1 = new TF1("dgaus",
		    "[1] * exp(-0.5 * ((x - [0]) / [3])^2) + [2] * exp(-0.5 * ((x - [0]) / [4])^2)",
		    diff_min, diff_max);
  f1->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
  TF1 *fit_1 = new TF1("fit_", "dgaus", diff_min, diff_max);
  fit_1->SetParameters(0, 0.8*hist1->GetMaximum(), 1e-2, 26.0, 175.0);
  fit_1->FixParameter(4, 175);
  fit_1->SetParLimits(3,0,1.E6);
  hist1->FitSlicesY(fit_1);
  TH1D *hist1_1 = (TH1D*)gDirectory->Get("hist1_0"); // Mean values
  TH1D *hist1_2 = (TH1D*)gDirectory->Get("hist1_3"); // Sigma values

  // Draw histogram
  hist1->Draw("COLZ");
  canvas->SaveAs("figs/track2dhist/2dhist_tracks.pdf(");
  canvas->SaveAs("figs/track2dhist/diff_tracks.pdf");
    
  hist1_1->Draw();
  canvas->SaveAs("figs/track2dhist/2dhist_tracks.pdf");
  canvas->SaveAs("figs/track2dhist/diff_trackfit_mean.pdf");

  // hist1_2->SetMaximum(30);
  // hist1_2->SetMinimum(0);
  hist1_2->Draw();
  canvas->SaveAs("figs/track2dhist/2dhist_tracks.pdf)");
  canvas->SaveAs("figs/track2dhist/diff_trackfit_sigma.pdf");

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);  // Align left-top

  TH1D *hSlice2 = hist1->ProjectionY("hSlice2", 1, 1);
  hSlice2->Scale(1/hSlice2->Integral());
  
  TH1D *hSlice3 = hist1->ProjectionY("hSlice3", 2, 2);
  hSlice3->Scale(1/hSlice3->Integral());
  
  TH1D *hSlice4 = hist1->ProjectionY("hSlice4", 3, 3);
  hSlice4->Scale(1/hSlice4->Integral());
  
  TH1D *hSlice5 = hist1->ProjectionY("hSlice5", 4, 4);
  hSlice5->Scale(1/hSlice5->Integral());
  
  TH1D *hSlice6 = hist1->ProjectionY("hSlice6", 5, 5);
  hSlice6->Scale(1/hSlice6->Integral());
  
  TH1D *hSlice7 = hist1->ProjectionY("hSlice7", 6, 6);
  hSlice7->Scale(1/hSlice7->Integral());
  
  TH1D *hSlice8 = hist1->ProjectionY("hSlice8", 7, 7);
  hSlice8->Scale(1/hSlice8->Integral());
  
  TH1D *hSlice9 = hist1->ProjectionY("hSlice9", 8, 8);
  hSlice9->Scale(1/hSlice9->Integral());

  float max_val = 1.2*std::max({
      hSlice2->GetMaximum(),
      hSlice3->GetMaximum(),
      hSlice4->GetMaximum(),
      hSlice5->GetMaximum(),
      hSlice6->GetMaximum(),
      hSlice7->GetMaximum(),
      hSlice8->GetMaximum(),
      hSlice9->GetMaximum(),
    });

  hSlice2->SetMaximum(max_val);
  hSlice2->SetMinimum(1e-3);
  hSlice3->SetMaximum(max_val);
  hSlice3->SetMinimum(1e-3);
  hSlice4->SetMaximum(max_val);
  hSlice4->SetMinimum(1e-3);
  hSlice5->SetMaximum(max_val);
  hSlice5->SetMinimum(1e-3);
  hSlice6->SetMaximum(max_val);
  hSlice6->SetMinimum(1e-3);
  hSlice7->SetMaximum(max_val);
  hSlice7->SetMinimum(1e-3);
  hSlice8->SetMaximum(max_val);
  hSlice8->SetMinimum(1e-3);
  hSlice9->SetMaximum(max_val);
  hSlice9->SetMinimum(1e-3);
  
  canvas->SetLogy(true);
  hSlice2->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf(");
  hSlice3->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf");
  hSlice4->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf");
  hSlice5->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf");
  hSlice6->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf");
  hSlice7->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf");
  hSlice8->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf");
  hSlice9->Draw("HIST");
  canvas->SaveAs("figs/track2dhist/slices.pdf)");
  canvas->SetLogy(false);

  // hist2->FitSlicesY();
  // TH1D *hist2_1 = (TH1D*)gDirectory->Get("hist2_1"); // Mean values
  // TH1D *hist2_2 = (TH1D*)gDirectory->Get("hist2_2"); // Sigma values

  // // Draw histogram
  // hist2->Draw("COLZ");
  // canvas->SaveAs("figs/2dhist_tracks.pdf");
  // canvas->SaveAs("figs/res_tracks.pdf");
    
  // hist2_1->Draw();
  // canvas->SaveAs("figs/2dhist_tracks.pdf");
  // canvas->SaveAs("figs/res_trackfit_mean.pdf");
    
  // hist2_2->Draw();
  // canvas->SaveAs("figs/2dhist_tracks.pdf)");
  // canvas->SaveAs("figs/res_trackfit_sigma.pdf");

  hist1_2->Print("all");
}
