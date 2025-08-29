#include "common_includes.h"

using namespace myutl;

void ideal_hgtd_efficiency() {
  const char* time_type = "Ideal Eff. HGTD";
  const char* file_prefix = "idealeff_nonear";
  bool smeared_times = true, valid_times = false;

  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  // All Needed Histogramming information
  std::map<ScoreType,TH1D*> inclusive_resos;
  std::map<ScoreType,TH1D*> inclusive_purity;
  for (ScoreType score: enum_vec) {
    inclusive_resos[score] = new TH1D(Form("inclusivereso_%s",toString(score)),
				      Form("Inclusive %s t_{0} - TruthVtx t_{0} (%s);#Delta t[ps];Entries", toString(score), time_type),
				      (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    inclusive_purity[score] = new TH1D(Form("inclusivepurity_%s",toString(score)),
				       Form("Inclusive Purity (%s);Purity;Entries",time_type),
				       (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);
  }

  PlotObj fjet     ("n Forward Jets",
		    time_type, Form("figs/%s_nfjet.pdf", file_prefix),
		    fjet_min     , fjet_max    , fjet_width    ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_fjet    , fold_fjet+fjet_width        );
  
  PlotObj ftrack   ("n Forward Tracks", time_type, Form("figs/%s_ntrack.pdf", file_prefix),
		    track_min    , track_max   , track_width   ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_track   , fold_track+track_width      );
  
  PlotObj pu_frac  ("Pile Up Fraction", time_type, Form("figs/%s_pufrac.pdf",file_prefix),
		    pu_frac_min  , pu_frac_max , pu_frac_width ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_pu_frac , pu_frac_max                 );
  
  PlotObj hs_track ("n Forward HS Tracks", time_type, Form("figs/%s_nhstrack.pdf",file_prefix),
		    hs_track_min , hs_track_max, hs_track_width,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_hs_track, fold_hs_track+hs_track_width);
  
  PlotObj pu_track ("n Forward PU Tracks", time_type, Form("figs/%s_nputrack.pdf",file_prefix),
		    pu_track_min , pu_track_max, pu_track_width,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_pu_track, fold_pu_track+pu_track_width);
  
  PlotObj recovtx_z("(Reco Vtx z) (mm)", time_type, Form("figs/%s_vtx_z.pdf",file_prefix),
		    z_min        , z_max       , z_width       ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    z_max        , z_max                       );

  TH2D* hs_pu_inclusive = new TH2D("hs_vs_pu","N HardScatter vs. N Pile-Up;n HS Tracks;n PU Tracks",
				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readnum = chain.GetReadEntry()+1;
    if (progress and readnum % 1000 == 0) std::cout << "Progress: " << readnum << "/" << chain.GetEntries() << "\n";

    process_event_data(&branch, smeared_times, valid_times, inclusive_resos, inclusive_purity,
		       fjet, ftrack, pu_frac, hs_track, pu_track, recovtx_z, hs_pu_inclusive);
  }
  
  // add all objects to plotobjs before drawing everything
  std::vector<PlotObj*> plots = {&fjet, &ftrack, &hs_track, &pu_track, &pu_frac, &recovtx_z};
  for (auto& plot: plots)
    plot->plot_postprocessing();
  std::cout << "FINISHED CREATING " << std::endl;

  for (auto& plot: plots)
    plot->plot_logic(canvas);
  std::cout << "FINISHED PRINTING BULK " << std::endl;
  
  plot_inclusive(Form("figs/%s_resos_logscale.pdf", file_prefix),
		 true, true, -400, 400, canvas, inclusive_resos);

  plot_inclusive(Form("figs/%s_resos_linscale.pdf", file_prefix),
		 false, true, -150, 150, canvas, inclusive_resos);

  plot_purity(Form("figs/%s_purity.pdf", file_prefix),
	      true, canvas, fjet.algolegend.get(), inclusive_purity);

  auto hs_pu_fname = Form("figs/%s_hs_v_pu.pdf", file_prefix);
  canvas->Print(Form("%s[",hs_pu_fname));
  hs_pu_inclusive->GetXaxis()->SetRangeUser(0, 50);
  hs_pu_inclusive->GetYaxis()->SetRangeUser(0, 50);
  hs_pu_inclusive->Draw("COLZ");
  canvas->Print(Form("%s",hs_pu_fname));
  canvas->Print(Form("%s]",hs_pu_fname));

  std::cout << "ALL DONE :3" << std::endl;
}
