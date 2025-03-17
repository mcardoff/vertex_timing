void generate_dt() {
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
  for (const auto& entry : std::filesystem::directory_iterator("./ntuple")) {
    if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
    chain.Add(entry.path().c_str());
    // break; // add 1
  }
  
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found in directory: " << std::endl;
    return;
  }

  TTreeReader reader(&chain);
  // jet variables
  TTreeReaderArray<float> jet_pt (reader, "AntiKt4EMTopoJets_pt");
  TTreeReaderArray<float> jet_eta(reader, "AntiKt4EMTopoJets_eta");
  TTreeReaderArray<std::vector<int>> jet_track_indices(reader, "AntiKt4EMTopoJets_track_idx");

  TTreeReaderArray<float> truth_hsjet_pt (reader, "TruthHSJet_pt");

  // vertex variables
  /// truth vertex
  TTreeReaderArray<float> truth_vtx_z   (reader, "TruthVtx_z");
  TTreeReaderArray<float> truth_vtx_time(reader, "TruthVtx_time");

  /// reco vertex
  TTreeReaderArray<float> reco_vtx_z      (reader, "RecoVtx_z");
  TTreeReaderArray<float> reco_vtx_time   (reader, "RecoVtx_time");
  TTreeReaderArray<float> reco_vtx_timeRes(reader, "RecoVtx_timeRes");
  TTreeReaderArray<int>   reco_vtx_valid  (reader, "RecoVtx_hasValidTime");

  // track variables
  TTreeReaderArray<float> track_z0(reader, "Track_z0");
  TTreeReaderArray<float> track_pt(reader, "Track_pt");
  TTreeReaderArray<float> track_eta(reader, "Track_eta");
  TTreeReaderArray<float> track_time(reader, "Track_time");
  TTreeReaderArray<float> track_time_res(reader, "Track_timeRes");
  TTreeReaderArray<float> track_var_z0(reader, "Track_var_z0");
  TTreeReaderArray<int>   track_to_truthvtx(reader, "Track_truthVtx_idx");
  TTreeReaderArray<int>   track_to_particle(reader, "Track_truthPart_idx");
  TTreeReaderArray<int>   track_time_valid(reader, "Track_hasValidTime");

  // particle vars
  TTreeReaderArray<float> prod_vtx_z  (reader, "TruthPart_prodVtx_z");

  int bins = 50;

  // Histograms
  double non_norm_bound = 200.0;
  /// Store vertex dt = reco_vtx_time[0] - truth_vtx_time[0]
  TH1F *hist1 = new TH1F("hist1", "RecoVtx t_{0} - TruthVtx t_{0};#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *hist2 = new TH1F("hist2", "RecoVtx t_{0} - TruthVtx t_{0},nForwardJet=1;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *hist3 = new TH1F("hist3", "RecoVtx t_{0} - TruthVtx t_{0},nForwardJet=2;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *hist4 = new TH1F("hist4", "RecoVtx t_{0} - TruthVtx t_{0},nForwardJet>2;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);

  /// Store dt/res
  TH1F *hist_norm1 = new TH1F("hist_norm1", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t};#Delta t;Entries", bins, -20, 20);
  TH1F *hist_norm2 = new TH1F("hist_norm2", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t},nForwardJet=1;#Delta t;Entries", bins, -20, 20);
  TH1F *hist_norm3 = new TH1F("hist_norm3", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t},nForwardJet=2;#Delta t;Entries", bins, -20, 20);
  TH1F *hist_norm4 = new TH1F("hist_norm4", "(RecoVtx t_{0} - TruthVtx t_{0})/#sigma_{t},nForwardJet>2;#Delta t;Entries", bins, -20, 20);

  // cut variables
  int    min_jets    = 2;    // min number of jets
  double min_jetpt   = 30.0; // self explanatory
  double min_abs_eta = 2.0;  // min eta for a vtx to be considered "forward"
  double min_abs_track_eta = 2.4;  // min eta for a track to be considered "forward"
  double min_track_pt = 1.0;  // track_pt > 1.0 GeV
  double max_vtx_dz  = 2.0;  // max error for reco HS vertex z
  double max_nsigma  = 2.0;  // how close a track can be to PV

  while (reader.Next()) {
    if (jet_pt.GetSize() < min_jets ) {
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
    
    // check if jet pt > 30 GeV
    float min_pt = *std::min_element(jet_pt.begin(), jet_pt.end());
    if(min_pt < min_jetpt) {
      if(debug)
	std::cout << "Skipping event due to low jet pt" << std::endl;
      continue;
    }

    int nForwardJet = 0;
    for(auto eta : jet_eta) {
      if (std::abs(eta) > min_abs_eta)
	nForwardJet++;
    }

    // new forward jet classification
    // int nForwardJet = 0;
    // for(auto jet_idx: jet_track_indices) {
    //   bool isForward = false;
    //   for(auto idx: jet_idx) {
    // 	if(std::abs(track_eta[idx]) > min_abs_eta) {
    // 	  // this is a forward jet
    // 	  isForward = true;
    // 	  break;
    // 	}
    //   }
    //   if(isForward)
    // 	nForwardJet++;
    // }

    if (debug)
      std::cout << "nForwardJet = " << nForwardJet << std::endl;

    if (reco_vtx_valid[0] == 1) {
      float diff = reco_vtx_time[0] - truth_vtx_time[0];
      float reso = diff / (reco_vtx_timeRes[0]);
    
      if (nForwardJet == 0)        { // exactly 0 forward jets
	hist1->Fill(diff);
	hist_norm1->Fill(reso);
      } else if (nForwardJet == 1) { // exactly 1 forward jet
	hist2->Fill(diff);
	hist_norm2->Fill(reso);
      } else if (nForwardJet == 2) { // exactly 2 forward jet
	hist3->Fill(diff);
	hist_norm3->Fill(reso);
      } else { // more
	hist4->Fill(diff);
	hist_norm4->Fill(reso);
      } 
    }
  }

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);  // Align left-top

  hist_norm1->Scale(1/hist_norm1->Integral());
  hist_norm2->Scale(1/hist_norm2->Integral());
  hist_norm3->Scale(1/hist_norm3->Integral());
  hist_norm4->Scale(1/hist_norm4->Integral());
  
  hist_norm1->SetMaximum(1.2*std::max({hist_norm1->GetMaximum(), hist_norm2->GetMaximum(), hist_norm3->GetMaximum(), hist_norm4->GetMaximum()}));
  hist_norm2->SetMaximum(hist_norm1->GetMaximum());
  hist_norm3->SetMaximum(hist_norm1->GetMaximum());
  hist_norm4->SetMaximum(hist_norm1->GetMaximum());

  canvas->SetLogy(true);
  
  hist_norm1->SetLineColor(c1);
  hist_norm1->SetLineWidth(2);
  hist_norm1->Draw("HIST");
  latex.DrawLatexNDC(0.15, 0.85, Form("N_{jets}\\geq %d, p_{T} > %0.2f GeV", min_jets, min_jetpt));

  hist_norm2->SetLineColor(c2);
  hist_norm2->SetLineWidth(2);
  hist_norm2->Draw("HIST SAME");

  hist_norm3->SetLineColor(c3);
  hist_norm3->SetLineWidth(2);
  hist_norm3->Draw("HIST SAME");

  hist_norm4->SetLineColor(c5);
  hist_norm4->SetLineWidth(2);
  hist_norm4->Draw("HIST SAME");

  double normfit_bound = 20;

  TF1 *normfit_1 = new TF1("normfit_1", "gaus", -normfit_bound, normfit_bound);
  hist_norm1->Fit(normfit_1, "R");
  TF1 *normfit_2 = new TF1("normfit_2", "gaus", -normfit_bound, normfit_bound);
  hist_norm2->Fit(normfit_2, "R");
  TF1 *normfit_3 = new TF1("normfit_3", "gaus", -normfit_bound, normfit_bound);
  hist_norm3->Fit(normfit_3, "R");
  TF1 *normfit_4 = new TF1("normfit_4", "gaus", -normfit_bound, normfit_bound);
  hist_norm4->Fit(normfit_4, "R");
  
  TLegend *normlegend = new TLegend(0.65, 0.65, 0.9, 0.9);
  normlegend->AddEntry(hist_norm1, Form("=0 Forward Jet, #sigma = %.2f", normfit_1->GetParameter(2)), "l");
  normlegend->AddEntry(hist_norm2, Form("=1 Forward Jet, #sigma = %.2f", normfit_2->GetParameter(2)), "l");
  normlegend->AddEntry(hist_norm3, Form("=2 Forward Jet, #sigma = %.2f", normfit_3->GetParameter(2)), "l");
  normlegend->AddEntry(hist_norm4, Form(">2 Forward Jet, #sigma = %.2f", normfit_4->GetParameter(2)), "l");
  normlegend->Draw();
  
  canvas->Print("dtplots.pdf(", "pdf");

  // individual plots with fits
  hist_norm1->Draw("HIST");
  normfit_1->SetLineColor(c1);
  normfit_1->Draw("SAME");
  normlegend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  hist_norm2->Draw("HIST");
  normfit_2->SetLineColor(c2);
  normfit_2->Draw("SAME");
  normlegend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  hist_norm3->Draw("HIST");
  normfit_3->SetLineColor(c3);
  normfit_3->Draw("SAME");
  normlegend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  hist_norm4->Draw("HIST");
  normfit_4->SetLineColor(c5);
  normfit_4->Draw("SAME");
  normlegend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  // NON DIVIDED BY RESOS
  hist1->Scale(1/hist1->Integral());
  hist2->Scale(1/hist2->Integral());
  hist3->Scale(1/hist3->Integral());
  hist4->Scale(1/hist4->Integral());

  hist1->SetMaximum(1.2*std::max({hist1->GetMaximum(), hist2->GetMaximum(), hist3->GetMaximum(), hist4->GetMaximum()}));
  hist2->SetMaximum(hist1->GetMaximum());
  hist3->SetMaximum(hist1->GetMaximum());
  hist4->SetMaximum(hist1->GetMaximum());
  
  hist1->SetLineColor(c1);
  hist1->SetLineWidth(2);
  hist1->Draw("HIST");
  latex.DrawLatexNDC(0.15, 0.85, Form("N_{jets}\\geq %d, p_{T} > %0.2f GeV", min_jets, min_jetpt));

  hist2->SetLineColor(c2);
  hist2->SetLineWidth(2);
  hist2->Draw("HIST SAME");

  hist3->SetLineColor(c3);
  hist3->SetLineWidth(2);
  hist3->Draw("HIST SAME");
  
  hist4->SetLineColor(c5);
  hist4->SetLineWidth(2);
  hist4->Draw("HIST SAME");

  double fit_bound = non_norm_bound;

  TF1 *f1 = new TF1("dgaus",
		    "[1] * exp(-0.5 * ((x - [0]) / [3])^2) + [2] * exp(-0.5 * ((x - [0]) / [4])^2)",
		    -fit_bound, fit_bound);
  f1->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");

  TF1 *fit_1 = new TF1("fit_1", "dgaus", -fit_bound, fit_bound);
  fit_1->SetParameters(0, hist1->GetMaximum(), 1e-2, 26.0, 175.0);
  fit_1->FixParameter(4, 175);
  hist1->Fit(fit_1, "R");
  TF1 *fit_2 = new TF1("fit_2", "dgaus", -fit_bound, fit_bound);
  fit_2->SetParameters(0, hist2->GetMaximum(), 1e-2, 26.0, 175.0);
  fit_2->FixParameter(4, 175);
  hist2->Fit(fit_2, "R");
  TF1 *fit_3 = new TF1("fit_3", "dgaus", -fit_bound, fit_bound);
  fit_3->SetParameters(0, hist3->GetMaximum(), 1e-2, 26.0, 175.0);
  fit_3->FixParameter(4, 175);
  hist3->Fit(fit_3, "R");
  TF1 *fit_4 = new TF1("fit_4", "dgaus", -fit_bound, fit_bound);
  fit_4->SetParameters(0, hist4->GetMaximum(), 1e-2, 26.0, 175.0);
  fit_4->FixParameter(4, 175);
  hist4->Fit(fit_4, "R");
    
  TLegend *legend = new TLegend(0.65, 0.65, 0.9, 0.9);
  legend->AddEntry(hist1, Form("=0 Forward Jet #sigma = %.2f#pm%.2f(%.2f%%)",
			       fit_1->GetParameter(3),
			       fit_1->GetParError(3),
			       100*fit_1->GetParError(3)/fit_1->GetParameter(3)), "l");
  legend->AddEntry(hist2, Form("=1 Forward Jet #sigma = %.2f#pm%.2f(%.2f%%)",
			       fit_2->GetParameter(3),
			       fit_2->GetParError(3),
			       100*fit_2->GetParError(3)/fit_2->GetParameter(3)), "l");
  legend->AddEntry(hist3, Form("=2 Forward Jet #sigma = %.2f#pm%.2f(%.2f%%)",
			       fit_3->GetParameter(3),
			       fit_3->GetParError(3),
			       100*fit_3->GetParError(3)/fit_3->GetParameter(3)), "l");
  legend->AddEntry(hist4, Form(">3 Forward Jet #sigma = %.2f#pm%.2f(%.2f%%)",
			       fit_4->GetParameter(3),
			       fit_4->GetParError(3),
			       100*fit_4->GetParError(3)/fit_4->GetParameter(3)), "l");
  
  legend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  // individual plots with fits
  hist1->Draw("HIST");
  fit_1->SetLineColor(c1);
  fit_1->Draw("SAME");
  legend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  hist2->Draw("HIST");
  fit_2->SetLineColor(c2);
  fit_2->Draw("SAME");
  legend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  hist3->Draw("HIST");
  fit_3->SetLineColor(c3);
  fit_3->Draw("SAME");
  legend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  hist4->Draw("HIST");
  fit_4->SetLineColor(c5);
  fit_4->Draw("SAME");
  legend->Draw();
  canvas->Print("dtplots.pdf)", "pdf");

  canvas->SetLogy(false);

  std::cout << "hist1 has " << hist1->GetEntries() << " Entries" << std::endl;

  if(debug) {
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
