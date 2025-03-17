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

  int track_min = 0, track_max = 1500;
  double res_min = -800.0, res_max = 800.0;
  int bins_x = (track_max-track_min)/30, bins_y=50;
  

  // Histograms
  /// Store vertex dt = reco_vtx_time[0] - truth_vtx_time[0]
  TH2F *hist1 = new TH2F(
			 "hist1",
			 "RecoVtx t_{0} - TruthVtx t_{0} vs Forward Tracks;n Forward Track;#Delta t[ps]",
			 bins_x,               // bin x
			 track_min, track_max, // bounds x
			 bins_y,               // bin y
			 res_min, res_max      // bounds y
			 );

  // cut variables
  int    min_jets    = 2;    // min number of jets
  double min_jetpt   = 30.0; // self explanatory
  double min_abs_eta = 2.0;  // min eta for a vtx to be considered "forward"
  double min_abs_track_eta = 2.4;  // min eta for a track to be considered "forward"
  double min_track_pt = 1.0;  // track_pt > 1.0 GeV
  double max_vtx_dz  = 2.0;   // max error for reco HS vertex z
  double max_nsigma  = 2.0;   // how close a track can be to PV

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

    int nForwardTrack = 0;
    for(int i = 0; i < track_eta.GetSize(); ++i) {
      auto eta = track_eta[i];
      auto valid = track_time_valid[i];
      if (std::abs(eta) > min_abs_track_eta && valid == 1)
	nForwardTrack++;
    }

    if (reco_vtx_valid[0] == 1) {
      float diff = reco_vtx_time[0] - truth_vtx_time[0];
      if (nForwardTrack > track_max || nForwardTrack < track_min) 
	std::cout << "!!!!!n_ftrack: " << nForwardTrack << std::endl;
      if (diff > res_max || diff < res_min)
	std::cout << "!!!!!res: " << diff << std::endl;
      // std::cout << "--------------------" << std::endl;
      // std::cout << "reso: " << diff << std::endl;
      // std::cout << "--------------------" << std::endl;
      hist1->Fill(nForwardTrack, diff); // increments value
      
    }
  }

  // Fit Gaussian to each slice along Y axis
  TH1D *h_mean = (TH1D*)hist1->ProjectionX("h_mean");
  TH1D *h_sigma = (TH1D*)hist1->ProjectionX("h_sigma");
    
  hist1->FitSlicesY();
  TH1D *hist1_1 = (TH1D*)gDirectory->Get("hist1_1"); // Mean values
  TH1D *hist1_2 = (TH1D*)gDirectory->Get("hist1_2"); // Sigma values

  // Draw histogram
  hist1->Draw("COLZ");
  canvas->SaveAs("tracks_res.pdf");
    
  hist1_1->Draw();
  canvas->SaveAs("fit_mean.pdf");
    
  hist1_2->Draw();
  canvas->SaveAs("fit_sigma.pdf");
}
