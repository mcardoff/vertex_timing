void generate_dt_tracks() {
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
  TTreeReaderArray<float> jet_pt (reader, "AntiKt4EMPFlowJets_pt");
  TTreeReaderArray<float> jet_eta(reader, "AntiKt4EMPFlowJets_eta");
  TTreeReaderArray<std::vector<int>> jet_track_indices(reader, "AntiKt4EMTopoJets_track_idx");

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
  TTreeReaderArray<int>   track_time_valid(reader, "Track_hasValidTime");

  int bins = 50;

  // Histograms
  double non_norm_bound = 200.0;
  /// Track dt
  TH1F *track_hist1 = new TH1F("trackhist1", "Track time - TruthVtx t_{0};#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *track_hist2 = new TH1F("trackhist2", "Track time - TruthVtx t_{0},nForwardJet=1;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *track_hist3 = new TH1F("trackhist3", "Track time - TruthVtx t_{0},nForwardJet=2;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *track_hist4 = new TH1F("trackhist4", "Track time - TruthVtx t_{0},nForwardJet>2;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);

  /// Track dt/res
  TH1F *track_hist_norm1 = new TH1F("trackhist_norm1", "(Track time - TruthVtx t_{0})/#sigma_{t};#Delta t;Entries", bins, -20, 20);
  TH1F *track_hist_norm2 = new TH1F("trackhist_norm2", "(Track time - TruthVtx t_{0})/#sigma_{t},nForwardJet=1;#Delta t;Entries", bins, -20, 20);
  TH1F *track_hist_norm3 = new TH1F("trackhist_norm3", "(Track time - TruthVtx t_{0})/#sigma_{t},nForwardJet=2;#Delta t;Entries", bins, -20, 20);
  TH1F *track_hist_norm4 = new TH1F("trackhist_norm4", "(Track time - TruthVtx t_{0})/#sigma_{t},nForwardJet>2;#Delta t;Entries", bins, -20, 20);
  
  /// Store vertex dt = reco_vtx_time[0] - truth_vtx_time[0]
  TH1F *hist1 = new TH1F("hist1", "Track time - TruthVtx t_{0};#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *hist2 = new TH1F("hist2", "Track time - TruthVtx t_{0},nForwardJet=1;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *hist3 = new TH1F("hist3", "Track time - TruthVtx t_{0},nForwardJet=2;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);
  TH1F *hist4 = new TH1F("hist4", "Track time - TruthVtx t_{0},nForwardJet>2;#Delta t;Entries", bins, -non_norm_bound, non_norm_bound);

  /// Store dt/res
  TH1F *hist_norm1 = new TH1F("hist_norm1", "(Track time - TruthVtx t_{0})/#sigma_{t};#Delta t;Entries", bins, -20, 20);
  TH1F *hist_norm2 = new TH1F("hist_norm2", "(Track time - TruthVtx t_{0})/#sigma_{t},nForwardJet=1;#Delta t;Entries", bins, -20, 20);
  TH1F *hist_norm3 = new TH1F("hist_norm3", "(Track time - TruthVtx t_{0})/#sigma_{t},nForwardJet=2;#Delta t;Entries", bins, -20, 20);
  TH1F *hist_norm4 = new TH1F("hist_norm4", "(Track time - TruthVtx t_{0})/#sigma_{t},nForwardJet>2;#Delta t;Entries", bins, -20, 20);

  /// Various flavors of ntracks
  //// Number of tracks with valid time in the HS vertex
  TH1I *nTracks1 = new TH1I("nTracks1", "nTracks with valid time in HS Vertex;nTracks;Entries", 80,-0.5,79.5);
  TH1I *nTracks2 = new TH1I("nTracks2", "nTracks with valid time in HS Vertex;nTracks;Entries", 80,-0.5,79.5);
  TH1I *nTracks3 = new TH1I("nTracks3", "nTracks with valid time in HS Vertex;nTracks;Entries", 80,-0.5,79.5);
  TH1I *nTracks4 = new TH1I("nTracks4", "nTracks with valid time in HS Vertex;nTracks;Entries", 80,-0.5,79.5);

  //// ntracks with valid time, including their HS classifications
  ///// Associated vertex (if it exists) is HS
  TH1I *nHSTracks1 = new TH1I("nHSTracks1", "nTracks Forward Jets = 0;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nHSTracks2 = new TH1I("nHSTracks2", "nTracks Forward Jets = 1;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nHSTracks3 = new TH1I("nHSTracks3", "nTracks Forward Jets = 2;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nHSTracks4 = new TH1I("nHSTracks4", "nTracks Forward Jets #geq 2;nTracks;Entries", 50,-0.5,49.5);

  ///// Associated vertex (if it exists) is NOT HS
  TH1I *nNonHSTracks1 = new TH1I("nNonHSTracks1", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks2 = new TH1I("nNonHSTracks2", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks3 = new TH1I("nNonHSTracks3", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks4 = new TH1I("nNonHSTracks4", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);

  ///// Catches tracks with no known vertex association
  TH1I *nUnkHSTracks1 = new TH1I("nUnkHSTracks1", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nUnkHSTracks2 = new TH1I("nUnkHSTracks2", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nUnkHSTracks3 = new TH1I("nUnkHSTracks3", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nUnkHSTracks4 = new TH1I("nUnkHSTracks4", "nTracks with valid time, tagged NOT HS;nTracks;Entries", 50,-0.5,49.5);

  //// forward tracks
  TH1I *nForTracks1 = new TH1I("nForTracks1", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);
  TH1I *nForTracks2 = new TH1I("nForTracks2", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);
  TH1I *nForTracks3 = new TH1I("nForTracks3", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);
  TH1I *nForTracks4 = new TH1I("nForTracks4", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);

  // cut variables
  int    min_jets    = 1;    // min number of jets
  double max_nsigma  = 3.0;  // how close a track can be to PV
  double min_jetpt   = 30.0; // self explanatory
  double min_abs_eta = 2.0;  // min eta for an obj to be considered "forward"
  double max_vtx_dz  = 2.0;  // max error for reco HS vertex z

  while (reader.Next()) {
    if (jet_pt.GetSize() < min_jets ) {
      if(debug)
	std::cout << "Skipping low jet event" << std::endl;
      continue;  // Skip if no jets are present
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

    // int nForwardJet = 0;
    // for(auto eta : jet_eta) {
    //   if (std::abs(eta) > min_abs_eta)
    // 	nForwardJet++;
    // }

    // new forward jet classification
    int nForwardJet = 0;
    for(auto jet_idx: jet_track_indices) {
      bool isForward = false;
      for(auto idx: jet_idx) {
	if(std::abs(track_eta[idx]) > min_abs_eta) {
	  // this is a forward jet
	  isForward = true;
	  break;
	}
      }
      if(isForward)
	nForwardJet++;
    }

    if (debug)
      std::cout << "nForwardJet = " << nForwardJet << std::endl;

    // ntracks calculation
    // counter for number of tracks with valid time
    int nvalidtime    = 0;
    // counter for number of forward tracks with valid time
    int nvalidftracks = 0;
    // counter for number of tracks with valid time (tagged as HS)
    int nTruthHS      = 0;
    // counter for number of tracks with valid time (tagged as assoc with a non-HS vertex)
    int nTruthNonHS   = 0;
    // counter for number of tracks with valid time (Not linked to a vertex)
    int nTruthUnk     = 0;
    // track loops
    for (int idx = 0; idx < track_z0.GetSize(); ++idx) {
      if(debug)
	std::cout << "track: " << idx << " hasvalidtime==" << track_time_valid[idx] << std::endl;

      if (!(track_time_valid[idx] == 1 && track_pt[idx] > 1.0))
	continue;
      
      float nsigma = (track_z0[idx] - reco_vtx_z[0])/std::sqrt(track_var_z0[idx]);

      if(std::abs(nsigma) > max_nsigma) {
	if(debug)
	  std::cout << "Skipping Track due to High n sigma" << std::endl;
	continue;
      }
      
      nvalidtime++;

      if(track_eta[idx] > min_abs_eta)
	nvalidftracks++;

      if(track_to_truthvtx[idx] == 0) {
	// associated to truth hs
	nTruthHS++;
      } else if (track_to_truthvtx[idx] == -1) {
	// unknown associated vertex (not used for vertex fitting)
	nTruthUnk++;
      } else {
	// associated to other primary vtx
	nTruthNonHS++;
      }

      float diff = track_time[idx] - truth_vtx_time[0];
      float reso = diff / track_time_res[idx];
      // exactly 0 forward jets
      if(nForwardJet == 0){
	track_hist1->Fill(diff);
	track_hist_norm1->Fill(reso);
      }
      // exactly 1 forward jet
      if(nForwardJet == 1){
	track_hist2->Fill(diff);
	track_hist_norm2->Fill(reso);
      }
      // exactly 2 forward jet
      if(nForwardJet == 2) {
	track_hist3->Fill(diff);
	track_hist_norm3->Fill(reso);
      }
      // > 2 forward jet
      if(nForwardJet > 2) {
	track_hist4->Fill(diff);
	track_hist_norm4->Fill(reso);
      }
    }

    // exactly 0 forward jets
    if(nForwardJet == 0){
      nTracks1->Fill(nvalidtime);
      nForTracks1->Fill(nvalidftracks);
      nHSTracks1->Fill(nTruthHS);
      nNonHSTracks1->Fill(nTruthNonHS);
      nUnkHSTracks1->Fill(nTruthUnk);
    }
    // exactly 1 forward jet
    if(nForwardJet == 1){
      nTracks2->Fill(nvalidtime);
      nForTracks2->Fill(nvalidftracks);
      nHSTracks2->Fill(nTruthHS);
      nNonHSTracks2->Fill(nTruthNonHS);
      nUnkHSTracks2->Fill(nTruthUnk);
    }
    // exactly 2 forward jet
    if(nForwardJet == 2) {
      nTracks3->Fill(nvalidtime);
      nForTracks3->Fill(nvalidftracks);
      nHSTracks3->Fill(nTruthHS);
      nNonHSTracks3->Fill(nTruthNonHS);
      nUnkHSTracks3->Fill(nTruthUnk);
    }
    // > 2 forward jet
    if(nForwardJet > 2) {
      nTracks4->Fill(nvalidtime);
      nForTracks4->Fill(nvalidftracks);
      nHSTracks4->Fill(nTruthHS);
      nNonHSTracks4->Fill(nTruthNonHS);
      nUnkHSTracks4->Fill(nTruthUnk);
    }
      
    if (reco_vtx_valid[0] == 1) {
    
      float diff = reco_vtx_time[0] - truth_vtx_time[0];
      float reso = diff / reco_vtx_timeRes[0];
    
      if(nForwardJet == 0) {
	hist1->Fill(diff);
	hist_norm1->Fill(reso);
      } else if (nForwardJet == 1) {
	hist2->Fill(diff);
	hist_norm2->Fill(reso);
      } else if(nForwardJet == 2) {
	hist3->Fill(diff);
	hist_norm3->Fill(reso);
      } else if(nForwardJet > 2) {
	hist4->Fill(diff);
	hist_norm4->Fill(reso);
      }
      // exactly 0 forward jets
      // exactly 1 forward jet
      // exactly 2 forward jet
      // > 2 forward jet
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

  // canvas->SetLogy(true);
  
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

  double normfit_bound = 5;

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
  legend->AddEntry(hist1, Form("=0 Forward Jet #sigma = %.2f", fit_1->GetParameter(3)), "l");
  legend->AddEntry(hist2, Form("=1 Forward Jet #sigma = %.2f", fit_2->GetParameter(3)), "l");
  legend->AddEntry(hist3, Form("=2 Forward Jet #sigma = %.2f", fit_3->GetParameter(3)), "l");
  legend->AddEntry(hist4, Form(">2 Forward Jet #sigma = %.2f", fit_4->GetParameter(3)), "l");
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
  canvas->Print("dtplots.pdf", "pdf");

  canvas->SetLogy(false);

  // ntracks plot
  nTracks1->SetMaximum(1.1*std::max({nTracks1->GetMaximum(), nTracks2->GetMaximum(), nTracks3->GetMaximum(), nTracks4->GetMaximum()}));
  
  nTracks1->SetLineColor(c1);
  nTracks1->SetLineWidth(2);
  nTracks1->Draw("HIST");

  nTracks2->SetLineColor(c2);
  nTracks2->SetLineWidth(2);
  nTracks2->Draw("HIST SAME");

  nTracks3->SetLineColor(c3);
  nTracks3->SetLineWidth(2);
  nTracks3->Draw("HIST SAME");
  
  nTracks4->SetLineColor(c5);
  nTracks4->SetLineWidth(2);
  nTracks4->Draw("HIST SAME");

  TLegend *nTracksLegend = new TLegend(0.65, 0.65, 0.9, 0.9);
  nTracksLegend->AddEntry(nTracks1, Form("=0 Forward Jet"), "l");
  nTracksLegend->AddEntry(nTracks2, Form("=1 Forward Jet"), "l");
  nTracksLegend->AddEntry(nTracks3, Form("=2 Forward Jet"), "l");
  nTracksLegend->AddEntry(nTracks4, Form(">2 Forward Jet"), "l");
  nTracksLegend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  // nfortracks plot
  nForTracks1->SetMaximum(1.1*std::max({nForTracks1->GetMaximum(), nForTracks2->GetMaximum(), nForTracks3->GetMaximum(), nForTracks4->GetMaximum()}));
  
  nForTracks1->SetLineColor(c1);
  nForTracks1->SetLineWidth(2);
  nForTracks1->Draw("HIST");

  nForTracks2->SetLineColor(c2);
  nForTracks2->SetLineWidth(2);
  nForTracks2->Draw("HIST SAME");

  nForTracks3->SetLineColor(c3);
  nForTracks3->SetLineWidth(2);
  nForTracks3->Draw("HIST SAME");
  
  nForTracks4->SetLineColor(c5);
  nForTracks4->SetLineWidth(2);
  nForTracks4->Draw("HIST SAME");

  TLegend *nForTracksLegend = new TLegend(0.65, 0.65, 0.9, 0.9);
  nForTracksLegend->AddEntry(nForTracks1, Form("=0 Forward Jet"), "l");
  nForTracksLegend->AddEntry(nForTracks2, Form("=1 Forward Jet"), "l");
  nForTracksLegend->AddEntry(nForTracks3, Form("=2 Forward Jet"), "l");
  nForTracksLegend->AddEntry(nForTracks4, Form(">2 Forward Jet"), "l");
  nForTracksLegend->Draw();
  canvas->Print("dtplots.pdf", "pdf");

  nHSTracks1->SetMaximum(1.1*std::max({nHSTracks1->GetMaximum(),nNonHSTracks1->GetMaximum(),nUnkHSTracks1->GetMaximum()}));

  nHSTracks1->SetLineColor(c1);
  nHSTracks1->SetLineWidth(2);
  nHSTracks1->Draw("HIST");

  nNonHSTracks1->SetLineColor(c2);
  nNonHSTracks1->SetLineWidth(2);
  nNonHSTracks1->Draw("HIST SAME");

  nUnkHSTracks1->SetLineColor(c3);
  nUnkHSTracks1->SetLineWidth(2);
  nUnkHSTracks1->Draw("HIST SAME");

  TLegend *nTracksClass1 = new TLegend(0.65, 0.65, 0.9, 0.9);
  nTracksClass1->AddEntry(nHSTracks1, Form("Truth HS"), "l");
  nTracksClass1->AddEntry(nNonHSTracks1, Form("Truth Non-HS"), "l");
  nTracksClass1->AddEntry(nUnkHSTracks1, Form("Truth Unknown"), "l");
  nTracksClass1->Draw();

  canvas->Print("dtplots.pdf", "pdf");

  nHSTracks2->SetMaximum(1.1*std::max({nHSTracks2->GetMaximum(),nNonHSTracks2->GetMaximum(),nUnkHSTracks2->GetMaximum()}));

  nHSTracks2->SetLineColor(c1);
  nHSTracks2->SetLineWidth(2);
  nHSTracks2->Draw("HIST");

  nNonHSTracks2->SetLineColor(c2);
  nNonHSTracks2->SetLineWidth(2);
  nNonHSTracks2->Draw("HIST SAME");

  nUnkHSTracks2->SetLineColor(c3);
  nUnkHSTracks2->SetLineWidth(2);
  nUnkHSTracks2->Draw("HIST SAME");

  TLegend *nTracksClass2 = new TLegend(0.65, 0.65, 0.9, 0.9);
  nTracksClass2->AddEntry(nHSTracks2, Form("Truth HS"), "l");
  nTracksClass2->AddEntry(nNonHSTracks2, Form("Truth Non-HS"), "l");
  nTracksClass2->AddEntry(nUnkHSTracks2, Form("Truth Unknown"), "l");
  nTracksClass2->Draw();

  canvas->Print("dtplots.pdf", "pdf");

  nHSTracks3->SetMaximum(1.1*std::max({nHSTracks3->GetMaximum(),nNonHSTracks3->GetMaximum(),nUnkHSTracks3->GetMaximum()}));

  nHSTracks3->SetLineColor(c1);
  nHSTracks3->SetLineWidth(2);
  nHSTracks3->Draw("HIST");

  nNonHSTracks3->SetLineColor(c2);
  nNonHSTracks3->SetLineWidth(2);
  nNonHSTracks3->Draw("HIST SAME");

  nUnkHSTracks3->SetLineColor(c3);
  nUnkHSTracks3->SetLineWidth(2);
  nUnkHSTracks3->Draw("HIST SAME");

  TLegend *nTracksClass3 = new TLegend(0.65, 0.65, 0.9, 0.9);
  nTracksClass3->AddEntry(nHSTracks3, Form("Truth HS"), "l");
  nTracksClass3->AddEntry(nNonHSTracks3, Form("Truth Non-HS"), "l");
  nTracksClass3->AddEntry(nUnkHSTracks3, Form("Truth Unknown"), "l");
  nTracksClass3->Draw();

  canvas->Print("dtplots.pdf", "pdf");

  nHSTracks4->SetMaximum(1.1*std::max({nHSTracks4->GetMaximum(),nNonHSTracks4->GetMaximum(),nUnkHSTracks4->GetMaximum()}));

  nHSTracks4->SetLineColor(c1);
  nHSTracks4->SetLineWidth(2);
  nHSTracks4->Draw("HIST");

  nNonHSTracks4->SetLineColor(c2);
  nNonHSTracks4->SetLineWidth(2);
  nNonHSTracks4->Draw("HIST SAME");

  nUnkHSTracks4->SetLineColor(c3);
  nUnkHSTracks4->SetLineWidth(2);
  nUnkHSTracks4->Draw("HIST SAME");

  TLegend *nTracksClass4 = new TLegend(0.65, 0.65, 0.9, 0.9);
  nTracksClass4->AddEntry(nHSTracks4, Form("Truth HS"), "l");
  nTracksClass4->AddEntry(nNonHSTracks4, Form("Truth Non-HS"), "l");
  nTracksClass4->AddEntry(nUnkHSTracks4, Form("Truth Unknown"), "l");
  nTracksClass4->Draw();

  canvas->Print("dtplots.pdf)", "pdf");

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
