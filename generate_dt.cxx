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
auto min_abs_eta       =  2.0;  // min eta for a jet to be considered "forward"
auto min_abs_track_eta =  2.4;  // min eta for a track to be considered "forward"
auto min_track_pt      =  1.0;  // track_pt > 1.0 GeV
auto max_vtx_dz        =  2.0;  // max error for reco HS vertex z
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

void set_maximum(TH1F *hist1, TH1F *hist2, TH1F *hist3, TH1F *hist4) {
  // assume they're normalized
  float max_val = 1.2*std::max({hist1->GetMaximum(), hist2->GetMaximum(), hist3->GetMaximum(), hist4->GetMaximum()});
  hist1->SetMaximum(max_val);
  hist2->SetMaximum(max_val);
  hist3->SetMaximum(max_val);
  hist4->SetMaximum(max_val);
  
  float min_val = 1e-3;
  hist1->SetMinimum(min_val);
  hist2->SetMinimum(min_val);
  hist3->SetMinimum(min_val);
  hist4->SetMinimum(min_val);
  return;
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
    legend->AddEntry(hist1, Form("=0 Forward Jet #sigma = %.2f, #chi^{2}=%.2f",
				 fit_1->GetParameter(3), fit_1->GetChisquare()), "l");
    legend->AddEntry(hist2, Form("=1 Forward Jet #sigma = %.2f, #chi^{2}=%.2f",
				 fit_2->GetParameter(3), fit_2->GetChisquare()), "l");
    legend->AddEntry(hist3, Form("=2 Forward Jet #sigma = %.2f, #chi^{2}=%.2f",
				 fit_3->GetParameter(3), fit_3->GetChisquare()), "l");
    legend->AddEntry(hist4, Form(">2 Forward Jet #sigma = %.2f, #chi^{2}=%.2f",
				 fit_4->GetParameter(3), fit_4->GetChisquare()), "l");
  } else {
    legend->AddEntry(hist1, Form("=0 Forward Jet, #sigma = %.2f, #chi^{2}=%.2f",
				 fit_1->GetParameter(2), fit_1->GetChisquare()), "l");
    legend->AddEntry(hist2, Form("=1 Forward Jet, #sigma = %.2f, #chi^{2}=%.2f",
				 fit_2->GetParameter(2), fit_2->GetChisquare()), "l");
    legend->AddEntry(hist3, Form("=2 Forward Jet, #sigma = %.2f, #chi^{2}=%.2f",
				 fit_3->GetParameter(2), fit_3->GetChisquare()), "l");
    legend->AddEntry(hist4, Form(">2 Forward Jet, #sigma = %.2f, #chi^{2}=%.2f",
				 fit_4->GetParameter(2), fit_4->GetChisquare()), "l");
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
  fit_1->SetParameters(hist->GetMaximum(),hist->GetMean(), 20);
  hist->Fit(fit_1, "R");
  // set color while we're here
  hist->SetLineColor(color);
  fit_1->SetLineColor(color);
  return fit_1;
}

void generate_dt() {
  gStyle->SetOptStat(0);
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
  TTreeReaderArray<float> pflowjet_pt
    (reader, "AntiKt4EMPFlowJets_pt");
  TTreeReaderArray<float> topojet_pt
    // (reader, "TruthHSJet_pt");
    (reader, "AntiKt4EMTopoJets_pt");
  TTreeReaderArray<float> pflowjet_eta
    (reader, "AntiKt4EMPFlowJets_eta");
  TTreeReaderArray<float> topojet_eta
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

  int bins = 250;
  double diff_bound = 1000.0;
  double res_bound = 110.0;

  // Histograms
  /// Store vertex dt = reco_vtx_time[0] - truth_vtx_time[0]
  TH1F *hist1 = new TH1F("hist1", "RecoVtx t_{0} - TruthVtx t_{0};#Delta t (ps);Entries", bins, -diff_bound, diff_bound);
  TH1F *hist2 = new TH1F("hist2", "RecoVtx t_{0} - TruthVtx t_{0};#Delta t (ps);Entries", bins, -diff_bound, diff_bound);
  TH1F *hist3 = new TH1F("hist3", "RecoVtx t_{0} - TruthVtx t_{0};#Delta t (ps);Entries", bins, -diff_bound, diff_bound);
  TH1F *hist4 = new TH1F("hist4", "RecoVtx t_{0} - TruthVtx t_{0};#Delta t (ps);Entries", bins, -diff_bound, diff_bound);

  /// Store dt/res
  TH1F *normhist1 = new TH1F("normhist1", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t};res(t);Entries", bins, -res_bound, res_bound);
  TH1F *normhist2 = new TH1F("normhist2", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t};res(t);Entries", bins, -res_bound, res_bound);
  TH1F *normhist3 = new TH1F("normhist3", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t};res(t);Entries", bins, -res_bound, res_bound);
  TH1F *normhist4 = new TH1F("normhist4", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t};res(t);Entries", bins, -res_bound, res_bound);

  TH1I *nForJetHist = new TH1I("forjethist", "number of fjets>2;nJet;Entries",10,-0.5,9.5);

  while (reader.Next()) {
    if (topojet_pt.GetSize() < min_jets ) {
      if(debug)
	std::cout << "Skipping low jet event" << std::endl;
      continue;  // Skip if no jets are present
    }

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

    // check if jet pt > 30 GeV, count forward jets
    float min_topopt  = *std::min_element(topojet_pt.begin() , topojet_pt.end() );
    if(min_topopt < min_jetpt) {	
      if(debug)
	std::cout << "Skipping event due to low jet pt" << std::endl;
      continue;
    }

    int nForwardJet = 0;
    for(auto eta: topojet_eta) {
      if (std::abs(eta) - min_abs_eta > 1e-6)
	nForwardJet++;
    }

    if (debug)
      std::cout << "nForwardJet = " << nForwardJet << std::endl;

    float diff = reco_vtx_time[0] - truth_vtx_time[0];
    float reso = diff / (reco_vtx_timeRes[0]);
    if (std::abs(diff) > diff_bound)
      std::cout << "DIFF OUT OF BOUNDS: " << diff << std::endl;

    if (std::abs(reso) > res_bound)
      std::cout << "RESO OUT OF BOUNDS: " << reso << std::endl;
    
    if (nForwardJet == 0)        { // exactly 0 forward jets
      hist1->Fill(diff);
      normhist1->Fill(reso);
    } else if (nForwardJet == 1) { // exactly 1 forward jet
      hist2->Fill(diff);
      normhist2->Fill(reso);
    } else if (nForwardJet == 2) { // exactly 2 forward jet
      hist3->Fill(diff);
      normhist3->Fill(reso);
    } else { // more
      nForJetHist->Fill(nForwardJet);
      hist4->Fill(diff);
      normhist4->Fill(reso);
    } 
  }

  // Normalize Histograms, set line width
  for(auto hist: {hist1, hist2, hist3, hist4, normhist1, normhist2, normhist3, normhist4}) {
    hist->Scale(1/hist->Integral());
    hist->SetLineWidth(2);
  }

  set_maximum(hist1, hist2, hist3, hist4);
  set_maximum(normhist1, normhist2, normhist3, normhist4);

  double fit_bound = diff_bound;
  double normfit_bound = res_bound;

  TF1 *f1 = new TF1("dgaus",
		    "[1] * exp(-0.5 * ((x - [0]) / [3])^2) + [2] * exp(-0.5 * ((x - [0]) / [4])^2)",
		    -fit_bound, fit_bound);
  f1->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
  
  TF1 *fit_1 = dgaus_fit(hist1, fit_bound, c1);
  TF1 *fit_2 = dgaus_fit(hist2, fit_bound, c2);
  TF1 *fit_3 = dgaus_fit(hist3, fit_bound, c3);
  TF1 *fit_4 = dgaus_fit(hist4, fit_bound, c5);

  TF1 *normfit_1 = gaus_fit(normhist1, normfit_bound, c1);
  TF1 *normfit_2 = gaus_fit(normhist2, normfit_bound, c2);
  TF1 *normfit_3 = gaus_fit(normhist3, normfit_bound, c3);
  TF1 *normfit_4 = gaus_fit(normhist4, normfit_bound, c5);

  draw_hist_with_fits(canvas, true, "figs/dt_plots.pdf",
		      hist1, hist2, hist3, hist4,
		      fit_1, fit_2, fit_3, fit_4);

  draw_hist_with_fits(canvas, false, "figs/resdt_plots.pdf",
		      normhist1, normhist2, normhist3, normhist4,
		      normfit_1, normfit_2, normfit_3, normfit_4);

  nForJetHist->SetMaximum(1.1*nForJetHist->GetMaximum());
  nForJetHist->Draw("HIST F");
  nForJetHist->SetLineWidth(2);
  nForJetHist->SetLineColor(c1);
  nForJetHist->SetFillColor(c1);
  canvas->SaveAs("figs/nforjet.pdf");
  
  if(debug) {
    std::cout << "hist1 has " << hist1->GetEntries() << " Entries" << std::endl;
    std::cout << "    fit 1 has " << 100*fit_1->GetParError(3)/fit_1->GetParameter(3) << "% error on sigma" << std::endl;
    std::cout << "    fit 2 has " << 100*fit_2->GetParError(3)/fit_2->GetParameter(3) << "% error on sigma" << std::endl;
    std::cout << "    fit 3 has " << 100*fit_3->GetParError(3)/fit_3->GetParameter(3) << "% error on sigma" << std::endl;
    std::cout << "    fit 4 has " << 100*fit_4->GetParError(3)/fit_4->GetParameter(3) << "% error on sigma" << std::endl;
    std::cout << "normfit 1 has " << 100*normfit_1->GetParError(2)/normfit_1->GetParameter(2) << "% error on sigma" << std::endl;
    std::cout << "normfit 2 has " << 100*normfit_2->GetParError(2)/normfit_2->GetParameter(2) << "% error on sigma" << std::endl;
    std::cout << "normfit 3 has " << 100*normfit_3->GetParError(2)/normfit_3->GetParameter(2) << "% error on sigma" << std::endl;
    std::cout << "normfit 4 has " << 100*normfit_4->GetParError(2)/normfit_4->GetParameter(2) << "% error on sigma" << std::endl;
  }
}
