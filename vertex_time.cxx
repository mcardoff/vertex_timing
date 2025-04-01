#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TStyle.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <boost/filesystem.hpp>

#define debug false

// cut variables
auto min_jets          =  2;    // min number of jets
auto min_jetpt         =  30.0; // self explanatory
auto min_abs_eta       =  2.0;  // min eta for a vtx to be considered "forward"
auto min_abs_track_eta =  2.4;  // min eta for a track to be considered "forward"
auto min_track_pt      =  1.0;  // track_pt > 1.0 GeV
auto max_vtx_dz        =  1.0e6;  // max error for reco HS vertex z
auto max_nsigma        =  3.0;  // how close a track can be to PV

void setup_chain(TChain &chain) {
  // VBF H->Invisible sample
  for (const auto& entry : boost::filesystem::directory_iterator("./ntuple")) {
    if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
    chain.Add(entry.path().c_str());
    // break; // add 1
  }
  
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found in directory: " << std::endl;
    return;
  }
}

void draw_hist_with_fits(TCanvas* canvas, bool is_dgaus, const char* fname,
			 TH1F *hist1, TH1F *hist2, TH1F *hist3, TH1F *hist4,
			 TF1  *fit_1, TF1  *fit_2, TF1  *fit_3, TF1  *fit_4) {
  canvas->SetLogy(true);
  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);  // Align left-top

  hist1->Draw("HIST");
  hist2->Draw("HIST SAME");
  hist3->Draw("HIST SAME");
  hist4->Draw("HIST SAME");
  latex.DrawLatexNDC(0.15, 0.85, Form("N_{jets}\\geq %d, p_{T} > %0.2f GeV", min_jets, min_jetpt));
    
  TLegend *legend = new TLegend(0.65, 0.65, 0.9, 0.9);
  if(is_dgaus){
    legend->AddEntry(hist1, Form("=0 Forward Jet #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_1->GetParameter(3), fit_1->GetChisquare()/fit_1->GetNDF()), "l");
    legend->AddEntry(hist2, Form("=1 Forward Jet #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_2->GetParameter(3), fit_2->GetChisquare()/fit_2->GetNDF()), "l");
    legend->AddEntry(hist3, Form("=2 Forward Jet #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_3->GetParameter(3), fit_3->GetChisquare()/fit_3->GetNDF()), "l");
    legend->AddEntry(hist4, Form(">2 Forward Jet #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_4->GetParameter(3), fit_4->GetChisquare()/fit_4->GetNDF()), "l");
  } else {
    legend->AddEntry(hist1, Form("=0 Forward Jet, #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_1->GetParameter(2), fit_1->GetChisquare()/fit_1->GetNDF()), "l");
    legend->AddEntry(hist2, Form("=1 Forward Jet, #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_2->GetParameter(2), fit_2->GetChisquare()/fit_2->GetNDF()), "l");
    legend->AddEntry(hist3, Form("=2 Forward Jet, #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_3->GetParameter(2), fit_3->GetChisquare()/fit_3->GetNDF()), "l");
    legend->AddEntry(hist4, Form(">2 Forward Jet, #sigma = %.2f, #chi/Ndf=%.2f",
				 fit_4->GetParameter(2), fit_4->GetChisquare()/fit_4->GetNDF()), "l");
  }

  legend->Draw();
  canvas->Print(Form("%s(",fname), "pdf");

  // individual plots with fits
  hist1->Draw("HIST");
  fit_1->Draw("SAME");
  latex.DrawLatexNDC(0.15, 0.85, "N(Forward Jet) = 0");
  legend->Draw();
  canvas->Print(fname, "pdf");

  hist2->Draw("HIST");
  fit_2->Draw("SAME");
  latex.DrawLatexNDC(0.15, 0.85, "N(Forward Jet) = 1");
  legend->Draw();
  canvas->Print(fname, "pdf");

  hist3->Draw("HIST");
  fit_3->Draw("SAME");
  latex.DrawLatexNDC(0.15, 0.85, "N(Forward Jet) = 2");
  legend->Draw();
  canvas->Print(fname, "pdf");

  hist4->Draw("HIST");
  fit_4->Draw("SAME");
  latex.DrawLatexNDC(0.15, 0.85, "N(Forward Jet) > 2");
  legend->Draw();
  canvas->Print(Form("%s)",fname), "pdf");
  canvas->SetLogy(false);
}

TF1* dgaus_fit(TH1F *hist, double fit_bound, Int_t color) {
  TF1 *fit_1 = new TF1(Form("fit_%s",hist->GetName()), "dgaus", -fit_bound, fit_bound);
  fit_1->SetParameters(0, hist->GetMaximum(), 1e-2, 26.0, 175.0);
  fit_1->FixParameter(4, 175);
  fit_1->SetParLimits(3,0,1.E6);
  hist->Fit(fit_1, "R");
  // set color while we're here
  hist->SetLineColor(color);
  fit_1->SetLineColor(color);
  return fit_1;
}

TF1* gaus_fit(TH1F *hist, double fit_bound, Int_t color) {
  TF1 *fit_1 = new TF1(Form("fit_%s",hist->GetName()), "gaus", -fit_bound, fit_bound);
  hist->Fit(fit_1, "R");
  // set color while we're here
  hist->SetLineColor(color);
  fit_1->SetLineColor(color);
  return fit_1;
}

void vertex_time() {
  // gStyle->SetOptStat(0);
  // colors
  auto c1 = TColor::GetColor("#3f90da");
  auto c2 = TColor::GetColor("#ffa90e");
  auto c3 = TColor::GetColor("#bd1f01");
  auto c4 = TColor::GetColor("#94a4a2");
  auto c5 = TColor::GetColor("#832db6");

  TChain chain ("ntuple");
  setup_chain(chain);
  TTreeReader reader(&chain);
  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);

  // jet variables
  TTreeReaderArray<float> jet_pt
    // (reader, "TruthHSJet_pt");
    (reader, "AntiKt4EMTopoJets_pt");
  TTreeReaderArray<float> jet_eta
    // (reader, "TruthHSJet_eta");
    (reader, "AntiKt4EMTopoJets_eta");
  TTreeReaderArray<std::vector<int>> jet_track_indices
    (reader, "AntiKt4EMTopoJets_track_idx");

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
  TTreeReaderArray<int>   reco_vtx_valid
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
  TTreeReaderArray<int>   track_to_truthvtx
    (reader, "Track_truthVtx_idx");
  TTreeReaderArray<int>   track_to_particle
    (reader, "Track_truthPart_idx");
  TTreeReaderArray<int>   track_time_valid
    (reader, "Track_hasValidTime");

  // particle vars
  TTreeReaderArray<float> prod_vtx_z
    (reader, "TruthPart_prodVtx_z");

  int bins = 150;
  double diff_bound = 2000.0;
  double res_bound = 100.0;

  // Histograms
  /// Store vertex dt = reco_vtx_time[0] - truth_vtx_time[0]
  TH1F *hist1 = new TH1F("hist1", "RecoVtx t_{0};#Delta t;Entries", bins, -diff_bound, diff_bound);

  while (reader.Next()) {
    // if (jet_pt.GetSize() < min_jets ) {
      // if(debug)
	// std::cout << "Skipping low jet event" << std::endl;
      // continue;  // Skip if no jets are present
    // }

    if(reco_vtx_valid[0] == 0) {
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
    
    // // check if jet pt > 30 GeV
    // float min_pt = *std::min_element(jet_pt.begin(), jet_pt.end());
    // if(min_pt < min_jetpt) {
    //   if(debug)
    // 	std::cout << "Skipping event due to low jet pt" << std::endl;
    //   continue;
    // }

    float time = reco_vtx_time[0];
    hist1->Fill(time);
  }

  // Normalize Histograms, set line width
  for(auto hist: {hist1}) {
    hist->Scale(1/hist->Integral());
    hist->SetLineWidth(2);
  }

  // NON DIVIDED BY RESOS
  hist1->SetMaximum(1.1*hist1->GetMaximum());
  hist1->Draw("HIST F");
  hist1->SetLineWidth(2);
  hist1->SetLineColor(c1);
  hist1->SetFillColor(c1);
  canvas->SaveAs("figs/vtx_time_dist.pdf");
}
