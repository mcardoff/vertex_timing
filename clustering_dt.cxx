#include "common_includes.h"

using namespace myutl;

void clustering_dt() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  // All Needed Histogramming information
  std::map<ScoreType,TH1D*> hgtdtimes_resos, idealres_resos, idealeff_resos, hgtdtimes_resos_pur, idealres_resos_pur, idealeff_resos_pur;
  std::map<ScoreType,TH1D*> hgtdtimes_purity, idealres_purity, idealeff_purity, hgtdtimes_purity_pur, idealres_purity_pur, idealeff_purity_pur;
  for (ScoreType score: enum_vec) {
    hgtdtimes_resos[score] = new TH1D(Form("hgtdtimes_reso_%s",toString(score)),
				      Form("Inclusive %s t_{0} - TruthVtx t_{0} (HGTD Times);#Delta t[ps];Entries", toString(score)),
				      (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    hgtdtimes_purity[score] = new TH1D(Form("hgtdtimes_purity_%s",toString(score)),
				       Form("Inclusive Purity (HGTD Times);Purity;Entries"),
				       (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);

    idealres_resos[score] = new TH1D(Form("idealres_reso_%s",toString(score)),
				     Form("Inclusive %s t_{0} - TruthVtx t_{0} (Ideal Res. HGTD);#Delta t[ps];Entries", toString(score)),
				     (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    idealres_purity[score] = new TH1D(Form("idealres_purity_%s",toString(score)),
				      Form("Inclusive Purity (Ideal Res. HGTD);Purity;Entries"),
				      (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);

    idealeff_resos[score] = new TH1D(Form("idealeff_reso_%s",toString(score)),
				     Form("Inclusive %s t_{0} - TruthVtx t_{0} (Ideal Eff. HGTD);#Delta t[ps];Entries", toString(score)),
				     (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    idealeff_purity[score] = new TH1D(Form("idealeff_purity_%s",toString(score)),
				      Form("Inclusive Purity (Ideal Eff. HGTD);Purity;Entries"),
				      (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);
    // PUR
    hgtdtimes_resos_pur[score] = new TH1D(Form("hgtdtimes_reso_pur_%s",toString(score)),
					  Form("Inclusive %s t_{0} - TruthVtx t_{0} (HGTD Times);#Delta t[ps];Entries", toString(score)),
					  (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    hgtdtimes_purity_pur[score] = new TH1D(Form("hgtdtimes_purity_pur_%s",toString(score)),
					   Form("Inclusive Purity (HGTD Times);Purity;Entries"),
					   (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);

    idealres_resos_pur[score] = new TH1D(Form("idealres_reso_pur_%s",toString(score)),
					 Form("Inclusive %s t_{0} - TruthVtx t_{0} (Ideal Res. HGTD);#Delta t[ps];Entries", toString(score)),
					 (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    idealres_purity_pur[score] = new TH1D(Form("idealres_purity_pur_%s",toString(score)),
					  Form("Inclusive Purity (Ideal Res. HGTD);Purity;Entries"),
					  (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);

    idealeff_resos_pur[score] = new TH1D(Form("idealeff_reso_pur_%s",toString(score)),
					 Form("Inclusive %s t_{0} - TruthVtx t_{0} (Ideal Eff. HGTD);#Delta t[ps];Entries", toString(score)),
					 (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
    idealeff_purity_pur[score] = new TH1D(Form("idealeff_purity_pur_%s",toString(score)),
					  Form("Inclusive Purity (Ideal Eff. HGTD);Purity;Entries"),
					  (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);
  }

  const char *nopur_file_prefix = "nonear", *pur_file_prefix = "near";
  // HGTD Times
  PlotObj fjet_hgtdtimes     ("n Forward Jets", "HGTD Times", 
			      Form("figs/hgtdtimes_nfjet.pdf"),
			      fjet_min  , fjet_max  , fjet_width  ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      fold_fjet, fold_fjet+fjet_width     );
  
  PlotObj ftrack_hgtdtimes   ("n Forward Tracks", "HGTD Times", 
			      Form("figs/hgtdtimes_ntrack.pdf"),
			      track_min , track_max , track_width ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      fold_track, fold_track+track_width  );
  
  PlotObj pu_frac_hgtdtimes  ("Pile Up Fraction", "HGTD Times", 
			      Form("figs/hgtdtimes_pufrac.pdf"),
			      pu_frac_min, pu_frac_max, pu_frac_width,
			      diff_min   , diff_max   , diff_width   ,
			      purity_min , purity_max , purity_width ,
			      fold_pu_frac, pu_frac_max              );
  
  PlotObj hs_track_hgtdtimes ("n Forward HS Tracks", "HGTD Times", 
			      Form("figs/hgtdtimes_nhstrack.pdf"),
			      hs_track_min, hs_track_max, hs_track_width ,
			      diff_min    , diff_max    , diff_width     ,
			      purity_min  , purity_max  , purity_width   ,
			      fold_hs_track, fold_hs_track+hs_track_width);
  
  PlotObj pu_track_hgtdtimes ("n Forward PU Tracks", "HGTD Times", 
			      Form("figs/hgtdtimes_nputrack.pdf"),
			      pu_track_min, pu_track_max, pu_track_width ,
			      diff_min    , diff_max    , diff_width     ,
			      purity_min  , purity_max  , purity_width   ,
			      fold_pu_track, fold_pu_track+pu_track_width);
  
  PlotObj recovtx_z_hgtdtimes("(Reco Vtx z) (mm)", "HGTD Times",
			      Form("figs/hgtdtimes_vtx_z.pdf"),
			      z_min     , z_max     , z_width     ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      z_max     , z_max                   );

  // Ideal HGTD Resolution
  // PlotObj fjet_idealres     ("n Forward Jets", "Ideal Res. HGTD", 
  // 			     Form("figs/idealres_nfjet.pdf"),
  // 			     fjet_min     , fjet_max    , fjet_width    ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_fjet    , fold_fjet+fjet_width        );
  
  // PlotObj ftrack_idealres   ("n Forward Tracks", "Ideal Res. HGTD", 
  // 			     Form("figs/idealres_ntrack.pdf"),
  // 			     track_min    , track_max   , track_width   ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_track   , fold_track+track_width      );
  
  // PlotObj pu_frac_idealres  ("Pile Up Fraction", "Ideal Res. HGTD", 
  // 			     Form("figs/idealres_pufrac.pdf"),
  // 			     pu_frac_min  , pu_frac_max , pu_frac_width ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_pu_frac , pu_frac_max                 );
  
  // PlotObj hs_track_idealres ("n Forward HS Tracks", "Ideal Res. HGTD",
  // 			     Form("figs/idealres_nhstrack.pdf"),
  // 			     hs_track_min , hs_track_max, hs_track_width,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_hs_track, fold_hs_track+hs_track_width);
  
  // PlotObj pu_track_idealres ("n Forward PU Tracks", "Ideal Res. HGTD", 
  // 			     Form("figs/idealres_nputrack.pdf"),
  // 			     pu_track_min , pu_track_max, pu_track_width,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_pu_track, fold_pu_track+pu_track_width);
  
  // PlotObj recovtx_z_idealres("(Reco Vtx z) (mm)", "Ideal Res. HGTD", 
  // 			     Form("figs/idealres_vtx_z.pdf"),
  // 			     z_min        , z_max       , z_width       ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     z_max        , z_max                       );

  // // Ideal HGTD Resolution
  // PlotObj fjet_idealeff     ("n Forward Jets", "Ideal Eff. HGTD", 
  // 			     Form("figs/idealeff_nfjet.pdf"),
  // 			     fjet_min     , fjet_max    , fjet_width    ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_fjet    , fold_fjet+fjet_width        );
  
  // PlotObj ftrack_idealeff   ("n Forward Tracks", "Ideal Eff. HGTD", 
  // 			     Form("figs/idealeff_ntrack.pdf"),
  // 			     track_min    , track_max   , track_width   ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_track   , fold_track+track_width      );
  
  // PlotObj pu_frac_idealeff  ("Pile Up Fraction", "Ideal Eff. HGTD", 
  // 			     Form("figs/idealeff_pufrac.pdf"),
  // 			     pu_frac_min  , pu_frac_max , pu_frac_width ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_pu_frac , pu_frac_max                 );
  
  // PlotObj hs_track_idealeff ("n Forward HS Tracks", "Ideal Eff. HGTD", 
  // 			     Form("figs/idealeff_nhstrack.pdf"),
  // 			     hs_track_min , hs_track_max, hs_track_width,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_hs_track, fold_hs_track+hs_track_width);
  
  // PlotObj pu_track_idealeff ("n Forward PU Tracks", "Ideal Eff. HGTD",
  // 			     Form("figs/idealeff_nputrack.pdf"),
  // 			     pu_track_min , pu_track_max, pu_track_width,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     fold_pu_track, fold_pu_track+pu_track_width);
  
  // PlotObj recovtx_z_idealeff("(Reco Vtx z) (mm)", "Ideal Eff. HGTD",
  // 			     Form("figs/idealeff_vtx_z.pdf"),
  // 			     z_min        , z_max       , z_width       ,
  // 			     diff_min     , diff_max    , diff_width    ,
  // 			     purity_min   , purity_max  , purity_width  ,
  // 			     z_max        , z_max                       );

  // // HGTD Times, PUR
  // PlotObj fjet_hgtdtimes_pur     ("n Forward Jets", "HGTD Times, PU Removal", 
  // 				  Form("figs/hgtdtimes_pur_nfjet.pdf"),
  // 				  fjet_min     , fjet_max    , fjet_width    ,
  // 				  diff_min     , diff_max    , diff_width    ,
  // 				  purity_min   , purity_max  , purity_width  ,
  // 				  fold_fjet    , fold_fjet+fjet_width        );
  
  // PlotObj ftrack_hgtdtimes_pur   ("n Forward Tracks", "HGTD Times, PU Removal", 
  // 				  Form("figs/hgtdtimes_pur_ntrack.pdf"),
  // 				  track_min    , track_max   , track_width   , 
  // 				  diff_min     , diff_max    , diff_width    ,
  // 				  purity_min   , purity_max  , purity_width  , 
  // 				  fold_track   , fold_track+track_width      );
  
  // PlotObj pu_frac_hgtdtimes_pur  ("Pile Up Fraction", "HGTD Times, PU Removal", 
  // 				  Form("figs/hgtdtimes_pur_pufrac.pdf"),
  // 				  pu_frac_min  , pu_frac_max , pu_frac_width , 
  // 				  diff_min     , diff_max    , diff_width    ,
  // 				  purity_min   , purity_max  , purity_width  , 
  // 				  fold_pu_frac , pu_frac_max                 );
  
  // PlotObj hs_track_hgtdtimes_pur ("n Forward HS Tracks", "HGTD Times, PU Removal", 
  // 				  Form("figs/hgtdtimes_pur_nhstrack.pdf"),
  // 				  hs_track_min , hs_track_max, hs_track_width,
  // 				  diff_min     , diff_max    , diff_width    ,
  // 				  purity_min   , purity_max  , purity_width  ,
  // 				  fold_hs_track, fold_hs_track+hs_track_width);
  
  // PlotObj pu_track_hgtdtimes_pur ("n Forward PU Tracks", "HGTD Times, PU Removal", 
  // 				  Form("figs/hgtdtimes_pur_nputrack.pdf"),
  // 				  pu_track_min , pu_track_max, pu_track_width,
  // 				  diff_min     , diff_max    , diff_width    ,
  // 				  purity_min   , purity_max  , purity_width  ,
  // 				  fold_pu_track, fold_pu_track+pu_track_width);
  
  // PlotObj recovtx_z_hgtdtimes_pur("(Reco Vtx z) (mm)", "HGTD Times, PU Removal", 
  // 				  Form("figs/hgtdtimes_pur_vtx_z.pdf"),
  // 				  z_min        , z_max       , z_width       ,
  // 				  diff_min     , diff_max    , diff_width    ,
  // 				  purity_min   , purity_max  , purity_width  ,
  // 				  z_max        , z_max                       );

  // // Ideal HGTD Resolution, PUR
  // PlotObj fjet_idealres_pur     ("n Forward Jets", "Ideal Res. HGTD, PU Removal", 
  // 				 Form("figs/idealres_pur_%s_nfjet.pdf", pur_file_prefix),
  // 				 fjet_min     , fjet_max    , fjet_width    ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_fjet    , fold_fjet+fjet_width        );
  
  // PlotObj ftrack_idealres_pur   ("n Forward Tracks", "Ideal Res. HGTD, PU Removal", 
  // 				 Form("figs/idealres_pur_%s_ntrack.pdf", pur_file_prefix),
  // 				 track_min    , track_max   , track_width   ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_track   , fold_track+track_width      );
  
  // PlotObj pu_frac_idealres_pur  ("Pile Up Fraction", "Ideal Res. HGTD, PU Removal", 
  // 				 Form("figs/idealres_pur_%s_pufrac.pdf", pur_file_prefix),
  // 				 pu_frac_min  , pu_frac_max , pu_frac_width ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_pu_frac , pu_frac_max                 );
  
  // PlotObj hs_track_idealres_pur ("n Forward HS Tracks", "Ideal Res. HGTD, PU Removal", 
  // 				 Form("figs/idealres_pur_%s_nhstrack.pdf", pur_file_prefix),
  // 				 hs_track_min , hs_track_max, hs_track_width,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_hs_track, fold_hs_track+hs_track_width);
  
  // PlotObj pu_track_idealres_pur ("n Forward PU Tracks", "Ideal Res. HGTD, PU Removal", 
  // 				 Form("figs/idealres_pur_%s_nputrack.pdf", pur_file_prefix),
  // 				 pu_track_min , pu_track_max, pu_track_width,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_pu_track, fold_pu_track+pu_track_width);
  
  // PlotObj recovtx_z_idealres_pur("(Reco Vtx z) (mm)", "Ideal Res. HGTD, PU Removal", 
  // 				 Form("figs/idealres_pur_%s_vtx_z.pdf",pur_file_prefix),
  // 				 z_min        , z_max       , z_width       ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 z_max        , z_max                       );

  // // Ideal HGTD Resolution+Efficiency, PUR
  // PlotObj fjet_idealeff_pur     ("n Forward Jets", "Ideal Eff. HGTD, PU Removal", 
  // 				 Form("figs/idealeff_pur_%s_nfjet.pdf",  pur_file_prefix),
  // 				 fjet_min     , fjet_max    , fjet_width    ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_fjet    , fold_fjet+fjet_width        );
  
  // PlotObj ftrack_idealeff_pur   ("n Forward Tracks", "Ideal Eff. HGTD, PU Removal", 
  // 				 Form("figs/idealeff_pur_%s_ntrack.pdf",  pur_file_prefix),
  // 				 track_min    , track_max   , track_width   ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_track   , fold_track+track_width      );
  
  // PlotObj pu_frac_idealeff_pur  ("Pile Up Fraction", "Ideal Eff. HGTD, PU Removal", 
  // 				 Form("figs/idealeff_pur_%s_pufrac.pdf", pur_file_prefix),
  // 				 pu_frac_min  , pu_frac_max , pu_frac_width ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_pu_frac , pu_frac_max                 );
  
  // PlotObj hs_track_idealeff_pur ("n Forward HS Tracks", "Ideal Eff. HGTD, PU Removal", 
  // 				 Form("figs/idealeff_pur_%s_nhstrack.pdf", pur_file_prefix),
  // 				 hs_track_min , hs_track_max, hs_track_width,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_hs_track, fold_hs_track+hs_track_width);
  
  // PlotObj pu_track_idealeff_pur ("n Forward PU Tracks", "Ideal Eff. HGTD, PU Removal", 
  // 				 Form("figs/idealeff_pur_%s_nputrack.pdf", pur_file_prefix),
  // 				 pu_track_min , pu_track_max, pu_track_width,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 fold_pu_track, fold_pu_track+pu_track_width);
  
  // PlotObj recovtx_z_idealeff_pur("(Reco Vtx z) (mm)", "Ideal Eff. HGTD, PU Removal", 
  // 				 Form("figs/idealeff_pur_%s_vtx_z.pdf", pur_file_prefix),
  // 				 z_min        , z_max       , z_width       ,
  // 				 diff_min     , diff_max    , diff_width    ,
  // 				 purity_min   , purity_max  , purity_width  ,
  // 				 z_max        , z_max                       );

  TH2D* hs_pu_hgtdtimes = new TH2D("hs_vs_pu_hgtdtimes","N HardScatter vs. N Pile-Up (HGTD Times);n HS Tracks;n PU Tracks",
				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  // TH2D* hs_pu_idealres = new TH2D("hs_vs_pu_idealres","N HardScatter vs. N Pile-Up (Ideal Res. HGTD);n HS Tracks;n PU Tracks",
  // 				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
  // 				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  // TH2D* hs_pu_idealeff = new TH2D("hs_vs_pu_idealeff","N HardScatter vs. N Pile-Up (Ideal Eff. HGTD);n HS Tracks;n PU Tracks",
  // 				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
  // 				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  // TH2D* hs_pu_hgtdtimes_pur = new TH2D("hs_vs_pu_hgtdtimes_pur","N HardScatter vs. N Pile-Up (HGTD Times, PU Removal);n HS Tracks;n PU Tracks",
  // 				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
  // 				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  // TH2D* hs_pu_idealres_pur = new TH2D("hs_vs_pu_idealres_pur","N HardScatter vs. N Pile-Up (Ideal Res. HGTD, PU Removal);n HS Tracks;n PU Tracks",
  // 				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
  // 				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  // TH2D* hs_pu_idealeff_pur = new TH2D("hs_vs_pu_idealeff_pur","N HardScatter vs. N Pile-Up (Ideal Eff. HGTD, PU Removal);n HS Tracks;n PU Tracks",
  // 				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
  // 				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);
  std::vector<TString> degraded_events, improved_events;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop" << std::endl;
  bool progress = false;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readnum = chain.GetReadEntry()+1;
    if (progress) std::cout << "Progress: " << readnum << "/" << chain.GetEntries() << "\n";

    if (readnum > chain.GetEntries()/4) break;
    auto pass_hgtd = process_event_data(&branch, false, true, false, hgtdtimes_resos, hgtdtimes_purity,
					fjet_hgtdtimes, ftrack_hgtdtimes, pu_frac_hgtdtimes, hs_track_hgtdtimes, pu_track_hgtdtimes, recovtx_z_hgtdtimes, hs_pu_hgtdtimes);

    // auto pass_idealres = process_event_data(&branch, true, true, false, idealres_resos, idealres_purity,
    // 		       fjet_idealres, ftrack_idealres, pu_frac_idealres, hs_track_idealres, pu_track_idealres, recovtx_z_idealres, hs_pu_idealres);
    
    // auto pass_idealeff = process_event_data(&branch, true, false, false, idealeff_resos, idealeff_purity,
		       // fjet_idealeff, ftrack_idealeff, pu_frac_idealeff, hs_track_idealeff, pu_track_idealeff, recovtx_z_idealeff, hs_pu_idealeff);
    
    //     if (pass_hgtd == 1 and pass_idealres == 0 and pass_idealeff == 0) {
    //       TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    //       TString file =  fileName(46,6);
    //       std::cout <<
    // 	Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0",
    // 	     file.Data(), this_evnt) << std::endl;
    // }
    // auto pass_pur = process_event_data(&branch, false, true, true, hgtdtimes_resos_pur, hgtdtimes_purity_pur,
    // 		       fjet_hgtdtimes_pur, ftrack_hgtdtimes_pur, pu_frac_hgtdtimes_pur,
    // 		       hs_track_hgtdtimes_pur, pu_track_hgtdtimes_pur, 
    // 		       recovtx_z_hgtdtimes_pur, hs_pu_hgtdtimes_pur);

    // std::cout << "Pass hgtd?" << pass_hgtd << std::endl;
    // std::cout << "Pass hgtd w/ PUR?" << pass_hgtd << std::endl;
    
    if (pass_hgtd.first == 0) {
      TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
      TString file =  fileName(46,6);
      improved_events.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0",
				     file.Data(), this_evnt));
    }
    // process_event_data(&branch, true, true, true, idealres_resos_pur, idealres_purity_pur,
    // 		       fjet_idealres_pur, ftrack_idealres_pur, pu_frac_idealres_pur,
    // 		       hs_track_idealres_pur, pu_track_idealres_pur, 
    // 		       recovtx_z_idealres_pur, hs_pu_idealres_pur);
    
    // process_event_data(&branch, true, false, true, idealeff_resos_pur, idealeff_purity_pur,
    // 		       fjet_idealeff_pur, ftrack_idealeff_pur, pu_frac_idealeff_pur,
    // 		       hs_track_idealeff_pur, pu_track_idealeff_pur, 
    // 		       recovtx_z_idealeff_pur, hs_pu_idealeff_pur);
  }
  
  // add all objects to plotobjs before drawing everything
  std::vector<PlotObj*> plots = {&fjet_hgtdtimes, &ftrack_hgtdtimes, &hs_track_hgtdtimes, &pu_track_hgtdtimes, &pu_frac_hgtdtimes, &recovtx_z_hgtdtimes,
				 &fjet_idealeff, &ftrack_idealeff, &hs_track_idealeff, &pu_track_idealeff, &pu_frac_idealeff, &recovtx_z_idealeff,
				 &fjet_idealres, &ftrack_idealres, &hs_track_idealres, &pu_track_idealres, &pu_frac_idealres, &recovtx_z_idealres,
				 &fjet_hgtdtimes_pur, &ftrack_hgtdtimes_pur, &hs_track_hgtdtimes_pur, &pu_track_hgtdtimes_pur, 
				 &pu_frac_hgtdtimes_pur, &recovtx_z_hgtdtimes_pur};
  for (auto& plot: plots)
    plot->plotPostprocessing();
  std::cout << "FINISHED CREATING " << std::endl;

  for (auto& plot: plots)
    plot->plotLogic(canvas);
  std::cout << "FINISHED PRINTING BULK " << std::endl;

  // plot_money(Form("figs/moneyplot_hstrack.pdf"), canvas,
  // 	     {&hs_track_hgtdtimes, &hs_track_hgtdtimes_pur, &hs_track_idealres, &hs_track_idealeff});

  // plot_money(Form("figs/moneyplot_pufrac.pdf"), canvas,
  // 	     {&pu_frac_hgtdtimes, &pu_frac_hgtdtimes_pur, &pu_frac_idealres, &pu_frac_idealeff});

  // plot_money(Form("figs/moneyplot_track.pdf"), canvas,
  // 	     {&ftrack_hgtdtimes, &ftrack_hgtdtimes_pur, &ftrack_idealres, &ftrack_idealeff});

  // plot_money(Form("figs/moneyplot_fjet.pdf"), canvas,
  // 	     {&fjet_hgtdtimes, &fjet_idealres, &fjet_idealeff, &fjet_hgtdtimes_pur});

  // comparison of w/wo pur
  plot_money(Form("figs/purcomp_hstrack.pdf"), canvas,
	     {&hs_track_hgtdtimes, &hs_track_hgtdtimes_pur});

  plot_money(Form("figs/purcomp_putrack.pdf"), canvas,
	     {&pu_track_hgtdtimes, &pu_track_hgtdtimes_pur});

  plot_money(Form("figs/purcomp_pufrac.pdf"), canvas,
	     {&pu_frac_hgtdtimes, &pu_frac_hgtdtimes_pur});

  plot_money(Form("figs/purcomp_track.pdf"), canvas,
	     {&ftrack_hgtdtimes, &ftrack_hgtdtimes_pur});

  plot_money(Form("figs/purcomp_fjet.pdf"), canvas,
	     {&fjet_hgtdtimes, &fjet_hgtdtimes_pur});
  
  plot_inclusive(Form("figs/hgtdtimes_resos_logscale.pdf"),
		 true, true, -400, 400, canvas, hgtdtimes_resos);

  plot_inclusive(Form("figs/hgtdtimes_resos_linscale.pdf"),
		 false, true, -150, 150, canvas, hgtdtimes_resos);

  plot_inclusive(Form("figs/idealres_resos_logscale.pdf"),
		 true, true, -400, 400, canvas, idealres_resos);

  plot_inclusive(Form("figs/idealres_resos_linscale.pdf"),
		 false, true, -150, 150, canvas, idealres_resos);

  plot_inclusive(Form("figs/idealeff_resos_logscale.pdf"),
		 true, true, -400, 400, canvas, idealeff_resos);

  plot_inclusive(Form("figs/idealeff_resos_linscale.pdf"),
		 false, true, -150, 150, canvas, idealeff_resos);

  plot_purity(Form("figs/hgtdtimes_purity.pdf"),
	      true, canvas, fjet_hgtdtimes.algolegend.get(), hgtdtimes_purity);

  plot_purity(Form("figs/idealres_purity.pdf"),
	      true, canvas, fjet_idealres.algolegend.get(), idealres_purity);

  plot_purity(Form("figs/idealeff_purity.pdf"),
	      true, canvas, fjet_idealeff.algolegend.get(), idealeff_purity);

  plot_inclusive(Form("figs/hgtdtimes_pur_resos_logscale.pdf"),
		 true, true, -400, 400, canvas, hgtdtimes_resos_pur);

  plot_inclusive(Form("figs/hgtdtimes_pur_resos_linscale.pdf"),
		 false, true, -150, 150, canvas, hgtdtimes_resos_pur);

  // plot_inclusive(Form("figs/idealres_pur_resos_logscale.pdf"),
  // 		 true, true, -400, 400, canvas, idealres_resos_pur);

  // plot_inclusive(Form("figs/idealres_pur_resos_linscale.pdf"),
  // 		 false, true, -150, 150, canvas, idealres_resos_pur);

  // plot_inclusive(Form("figs/idealeff_pur_resos_logscale.pdf"),
  // 		 true, true, -400, 400, canvas, idealeff_resos_pur);

  // plot_inclusive(Form("figs/idealeff_pur_resos_linscale.pdf"),
  // 		 false, true, -150, 150, canvas, idealeff_resos_pur);

  plot_purity(Form("figs/hgtdtimes_pur_purity.pdf"),
	      true, canvas, fjet_hgtdtimes_pur.algolegend.get(), hgtdtimes_purity_pur);

  // plot_purity(Form("figs/idealres_%s_purity.pdf", pur_file_prefix),
  // 	      true, canvas, fjet_idealres_pur.algolegend.get(), idealres_purity_pur);

  // plot_purity(Form("figs/idealeff_%s_purity.pdf", pur_file_prefix),
  // 	      true, canvas, fjet_idealeff_pur.algolegend.get(), idealeff_purity_pur);

  auto hs_pu_fname = Form("figs/%s_hs_v_pu.pdf", nopur_file_prefix);
  canvas->Print(Form("%s[",hs_pu_fname));
  hs_pu_hgtdtimes->GetXaxis()->SetRangeUser(0, 50);
  hs_pu_hgtdtimes->GetYaxis()->SetRangeUser(0, 50);
  hs_pu_hgtdtimes->Draw("COLZ");
  canvas->Print(Form("%s",hs_pu_fname));

  hs_pu_idealres->GetXaxis()->SetRangeUser(0, 50);
  hs_pu_idealres->GetYaxis()->SetRangeUser(0, 50);
  hs_pu_idealres->Draw("COLZ");
  canvas->Print(Form("%s",hs_pu_fname));

  hs_pu_idealeff->GetXaxis()->SetRangeUser(0, 50);
  hs_pu_idealeff->GetYaxis()->SetRangeUser(0, 50);
  hs_pu_idealeff->Draw("COLZ");
  canvas->Print(Form("%s",hs_pu_fname));
  
  canvas->Print(Form("%s]",hs_pu_fname));

  std::cout << "ALL DONE :3" << std::endl;

  std::cout << "DEGRADED EVENTS: " << std::endl;
  for (auto s: degraded_events)
    std::cout << s << std::endl;

  std::cout << "IMPROVED EVENTS: " << std::endl;
  for (auto s: improved_events)
    std::cout << s << std::endl;
}
