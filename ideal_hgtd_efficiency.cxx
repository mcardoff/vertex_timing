#include "common_includes.h"

using namespace myutl;

void ideal_hgtd_efficiency() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  double diff_min = -1000.0, diff_max = 1000.0;
  double diff_width = 2.0;

  double purity_min = 0, purity_max = 1;
  double purity_width = 0.05;

  double fjet_min = 1, fjet_max = 31, fold_fjet = 5;
  double fjet_width = 1.0;

  double track_min = 0, track_max = 100, fold_track = 20;
  double track_width = 1.0;

  double pu_track_min = track_min, pu_track_max = track_max, fold_hs_track = fold_track;
  double pu_track_width = track_width;

  double hs_track_min = track_min, hs_track_max = track_max, fold_pu_track = fold_track;
  double hs_track_width = track_width;

  double pu_frac_min = 0, pu_frac_max = 1 , fold_pu_frac = pu_frac_max;
  double pu_frac_width = 0.05;

  double z_min = -200, z_max = 200, fold_z = 100;
  double z_width = 10.0;

  // All Needed Histogramming information
  const char* time_type = "Ideal Eff. HGTD";
  std::map<ScoreType,TH1D*> inclusive_resos;
  std::map<ScoreType,TH1D*> inclusive_purity;
  for (ScoreType score: enum_vec) {
    inclusive_resos[score] = new TH1D(Form("inclusivereso_%s",toString(score)),
				      Form("Inclusive %s t_{0} - TruthVtx t_{0} (%s);#Delta t[ps];Entries",
					   toString(score), time_type),
				      (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    inclusive_purity[score] = new TH1D(Form("inclusivepurity_%s",toString(score)), Form("Inclusive Purity (%s);Purity;Entries",time_type),
				       (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);
  }

  PlotObj fjet     ("n Forward Jets", time_type,
		    "fjet", "figs/idealeff_nfjet.pdf",
		    fjet_min     , fjet_max    , fjet_width    ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_fjet    , fold_fjet+fjet_width        );
  
  PlotObj ftrack   ("n Forward Tracks", time_type,
		    "track", "figs/idealeff_ntrack.pdf",
		    track_min    , track_max   , track_width   ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_track   , fold_track+track_width      );
  
  PlotObj pu_frac  ("Pile Up Fraction", time_type,
		    "pu_frac", "figs/idealeff_pufrac.pdf",
		    pu_frac_min  , pu_frac_max , pu_frac_width ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_pu_frac , pu_frac_max                 );
  
  PlotObj hs_track ("n Forward HS Tracks", time_type,
		    "hs_track", "figs/idealeff_nhstrack.pdf",
		    hs_track_min , hs_track_max, hs_track_width,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_hs_track, fold_hs_track+hs_track_width);
  
  PlotObj pu_track ("n Forward PU Tracks", time_type,
		    "pu_track", "figs/idealeff_nputrack.pdf",
		    pu_track_min , pu_track_max, pu_track_width,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_pu_track, fold_pu_track+pu_track_width);
  
  PlotObj recovtx_z("(Reco Vtx z) (mm)", time_type,
		    "recovtx_z", "figs/idealeff_vtx_z.pdf",
		    z_min        , z_max       , z_width       ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    z_max        , z_max                       );

  TH2D* hs_pu_inclusive = new TH2D("hs_vs_pu","N HardScatter vs. N Pile-Up;n HS Tracks;n PU Tracks",
				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  bool smeared_times = true, valid_times = false;

  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed

    auto readnum = chain.GetReadEntry()+1;

    // if (readnum > 1000)
    //   continue;
    if (progress and readnum % 1000 == 0) std::cout << "Progress: " << readnum << "/" << chain.GetEntries() << "\n";

    if (not branch.pass_basic_cuts()) // check n_jet cut and vertex dz cut
      continue;

    // check if there is one forward jet with pt > 30 GeV
    int nForwardJet=0;
    branch.count_forward_jets(nForwardJet);
    
    if (not branch.pass_jet_pt_cut()) continue;

    std::vector<int> tracks = getAssociatedTracks(&branch, min_track_pt);

    int nForwardTrack=0, nForwardTrack_HS=0, nForwardTrack_PU=0;
    branch.count_forward_tracks(nForwardTrack,nForwardTrack_HS,nForwardTrack_PU,tracks);

    hs_pu_inclusive->Fill(nForwardTrack_HS,nForwardTrack_PU);
    
    double pu_ratio = (double)nForwardTrack_PU / (double)nForwardTrack;
    double reco_z = branch.reco_vtx_z[0];

    auto eff_fill_val_fjet     = folded(nForwardJet     , (int)fold_fjet) ;
    auto eff_fill_val_track    = folded(nForwardTrack   , (int)fold_track);
    auto eff_fill_val_hs_track = folded(nForwardTrack_HS, (int)fold_hs_track);
    auto eff_fill_val_pu_track = folded(nForwardTrack_PU, (int)fold_pu_track);
    auto eff_fill_val_z        = reco_z;
    auto eff_fill_val_pu_ratio = pu_ratio;

    if (debug) {
      std::cout << "nForwardJet = " << nForwardJet << std::endl;
      std::cout << "nForwardTrack = " << nForwardTrack << std::endl;
      std::cout << "nForwardTrack_HS = " << nForwardTrack_HS << std::endl;
      std::cout << "nForwardTrack_PU = " << nForwardTrack_PU << std::endl;
      std::cout << "Vertex_z = " << reco_z << std::endl;
    }
    
    std::vector<Cluster> cluster =
      clusterTracksInTime(tracks, &branch, 3.0, 10.0, smeared_times, valid_times, true);

    fjet.eff_total->     Fill(eff_fill_val_fjet    );
    ftrack.eff_total->   Fill(eff_fill_val_track   );
    pu_frac.eff_total->  Fill(eff_fill_val_pu_ratio);
    hs_track.eff_total-> Fill(eff_fill_val_hs_track);
    pu_track.eff_total-> Fill(eff_fill_val_pu_track);
    recovtx_z.eff_total->Fill(eff_fill_val_z       );
    
    // std::cout << "Acessing cluster size" << std::endl;
    if (cluster.size() == 0)
      continue;

    std::map<ScoreType,Cluster> chosen = chooseCluster(cluster, &branch);

    // run HGTD Clustering (simultaneous)
    std::vector<Cluster> hgtd_clusters =
      clusterTracksInTime(tracks, &branch, 3.0, -1, false, true, false);

    if (branch.reco_vtx_valid[0] == 1)
      chosen[ScoreType::HGTD] = chooseHGTDCluster(hgtd_clusters, &branch);

    for (ScoreType score : enum_vec) {
      
      if (branch.reco_vtx_valid[0] == 0 and score == HGTD)
	continue;
      
      if (fjet.eff_pass.count(score) != 1)
	continue;
      
      if (score == ScoreType::INVALID)
	continue;
      
      Cluster scored = chosen[score];
      double score_based_time = scored.values.at(0);
      double cluster_purity = scored.purity;
      double diff = score_based_time - branch.truth_vtx_time[0];
      inclusive_resos.at(score)->Fill(diff);
      inclusive_purity.at(score)->Fill(cluster_purity);
      
      if (passEfficiency(diff, scored, &branch)) {
	fjet.eff_pass.at(score)->     Fill(eff_fill_val_fjet    );
	ftrack.eff_pass.at(score)->   Fill(eff_fill_val_track   );
	pu_frac.eff_pass.at(score)->  Fill(eff_fill_val_pu_ratio);
	hs_track.eff_pass.at(score)-> Fill(eff_fill_val_hs_track);
	pu_track.eff_pass.at(score)-> Fill(eff_fill_val_pu_track);
	recovtx_z.eff_pass.at(score)->Fill(eff_fill_val_z       );
      }
      
      // fill diff hists
      fjet.hist.at(score)->     Fill(nForwardJet     , diff);
      ftrack.hist.at(score)->   Fill(nForwardTrack   , diff);
      pu_frac.hist.at(score)->  Fill(pu_ratio        , diff);
      hs_track.hist.at(score)-> Fill(nForwardTrack_HS, diff);
      pu_track.hist.at(score)-> Fill(nForwardTrack_PU, diff);
      recovtx_z.hist.at(score)->Fill(reco_z          , diff);

      // fill purities
      fjet.purity.at(score)->     Fill(nForwardJet     , cluster_purity);
      ftrack.purity.at(score)->   Fill(nForwardTrack   , cluster_purity);
      pu_frac.purity.at(score)->  Fill(pu_ratio        , cluster_purity);
      hs_track.purity.at(score)-> Fill(nForwardTrack_HS, cluster_purity);
      pu_track.purity.at(score)-> Fill(nForwardTrack_PU, cluster_purity);
      recovtx_z.purity.at(score)->Fill(reco_z          , cluster_purity);
    }
  }
  
  // add all objects to plotobjs before drawing everything
  std::vector<PlotObj*> plots = {&fjet, &ftrack, &hs_track, &pu_track, &pu_frac, &recovtx_z};
  for (auto& plot: plots)
    plot->plot_postprocessing();

  std::cout << "FINISHED CREATING " << std::endl;
  
  for (auto& plot: plots)
    plot->plot_logic(canvas);  

  plot_inclusive("figs/idealeff_resos_logscale.pdf",
		 true, true, -400, 400,
		 canvas, inclusive_resos);

  plot_inclusive("figs/idealeff_resos_linscale.pdf",
		 false, true, -150, 150,
		 canvas, inclusive_resos);

  plot_purity("figs/idealeff_purity.pdf", false,
	      canvas, fjet.algolegend.get(), inclusive_purity);
  
  auto hs_pu_fname = "figs/idealeff_hs_v_pu.pdf";
  canvas->Print(Form("%s[",hs_pu_fname));
  hs_pu_inclusive->GetXaxis()->SetRangeUser(0, 35);
  hs_pu_inclusive->GetYaxis()->SetRangeUser(0, 50);
  hs_pu_inclusive->Draw("COLZ");
  canvas->Print(Form("%s",hs_pu_fname));
  canvas->Print(Form("%s]",hs_pu_fname));
}
