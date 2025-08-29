#include "common_includes.h"

using namespace myutl;

void clustering_dt() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../../ntuple-hgtd/");
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

  PlotObj fjet_idealres     ("n Forward Jets", "Ideal Res. HGTD", 
			      Form("figs/idealres_nfjet.pdf"),
			      fjet_min  , fjet_max  , fjet_width  ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      fold_fjet, fold_fjet+fjet_width     );
  
  PlotObj ftrack_idealres   ("n Forward Tracks", "Ideal Res. HGTD", 
			      Form("figs/idealres_ntrack.pdf"),
			      track_min , track_max , track_width ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      fold_track, fold_track+track_width  );
  
  PlotObj pu_frac_idealres  ("Pile Up Fraction", "Ideal Res. HGTD", 
			      Form("figs/idealres_pufrac.pdf"),
			      pu_frac_min, pu_frac_max, pu_frac_width,
			      diff_min   , diff_max   , diff_width   ,
			      purity_min , purity_max , purity_width ,
			      fold_pu_frac, pu_frac_max              );
  
  PlotObj hs_track_idealres ("n Forward HS Tracks", "Ideal Res. HGTD", 
			      Form("figs/idealres_nhstrack.pdf"),
			      hs_track_min, hs_track_max, hs_track_width ,
			      diff_min    , diff_max    , diff_width     ,
			      purity_min  , purity_max  , purity_width   ,
			      fold_hs_track, fold_hs_track+hs_track_width);
  
  PlotObj pu_track_idealres ("n Forward PU Tracks", "Ideal Res. HGTD", 
			      Form("figs/idealres_nputrack.pdf"),
			      pu_track_min, pu_track_max, pu_track_width ,
			      diff_min    , diff_max    , diff_width     ,
			      purity_min  , purity_max  , purity_width   ,
			      fold_pu_track, fold_pu_track+pu_track_width);
  
  PlotObj recovtx_z_idealres("(Reco Vtx z) (mm)", "Ideal Res. HGTD",
			      Form("figs/idealres_vtx_z.pdf"),
			      z_min     , z_max     , z_width     ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      z_max     , z_max                   );

  PlotObj fjet_idealeff     ("n Forward Jets", "Ideal Eff. HGTD", 
			      Form("figs/idealeff_nfjet.pdf"),
			      fjet_min  , fjet_max  , fjet_width  ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      fold_fjet, fold_fjet+fjet_width     );
  
  PlotObj ftrack_idealeff   ("n Forward Tracks", "Ideal Eff. HGTD", 
			      Form("figs/idealeff_ntrack.pdf"),
			      track_min , track_max , track_width ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      fold_track, fold_track+track_width  );
  
  PlotObj pu_frac_idealeff  ("Pile Up Fraction", "Ideal Eff. HGTD", 
			      Form("figs/idealeff_pufrac.pdf"),
			      pu_frac_min, pu_frac_max, pu_frac_width,
			      diff_min   , diff_max   , diff_width   ,
			      purity_min , purity_max , purity_width ,
			      fold_pu_frac, pu_frac_max              );
  
  PlotObj hs_track_idealeff ("n Forward HS Tracks", "Ideal Eff. HGTD", 
			      Form("figs/idealeff_nhstrack.pdf"),
			      hs_track_min, hs_track_max, hs_track_width ,
			      diff_min    , diff_max    , diff_width     ,
			      purity_min  , purity_max  , purity_width   ,
			      fold_hs_track, fold_hs_track+hs_track_width);
  
  PlotObj pu_track_idealeff ("n Forward PU Tracks", "Ideal Eff. HGTD", 
			      Form("figs/idealeff_nputrack.pdf"),
			      pu_track_min, pu_track_max, pu_track_width ,
			      diff_min    , diff_max    , diff_width     ,
			      purity_min  , purity_max  , purity_width   ,
			      fold_pu_track, fold_pu_track+pu_track_width);
  
  PlotObj recovtx_z_idealeff("(Reco Vtx z) (mm)", "Ideal Eff. HGTD",
			      Form("figs/idealres_vtx_z.pdf"),
			      z_min     , z_max     , z_width     ,
			      diff_min  , diff_max  , diff_width  ,
			      purity_min, purity_max, purity_width,
			      z_max     , z_max                   );
  
  TH2D* hs_pu_hgtdtimes = new TH2D("hs_vs_pu_hgtdtimes","N HardScatter vs. N Pile-Up (HGTD Times);n HS Tracks;n PU Tracks",
				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);
  std::vector<TString> smear_degrades, smear_improves, both_fail, both_pass;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readnum = chain.GetReadEntry()+1;
    if (progress && readnum % 1000 == 0) std::cout << "Progress: " << readnum << "/" << chain.GetEntries() << "\n";

    if (smear_degrades.size() >= 50) break;
    if (readnum > 10e3) break;
    // if (hs_track_hgtdtimes.hist.at(ScoreType::PASS)->Integral() < 10e3) { // break;
    auto pass_hgtdtimes = process_event_data(&branch, false, true, false, true, hgtdtimes_resos, hgtdtimes_purity,
					     fjet_hgtdtimes, ftrack_hgtdtimes, pu_frac_hgtdtimes,
					     hs_track_hgtdtimes, pu_track_hgtdtimes, recovtx_z_hgtdtimes, hs_pu_hgtdtimes);
    // }
    // if (hs_track_idealres.hist.at(ScoreType::PASS)->Integral() < 10e3) // break;
    auto pass_idealres = process_event_data(&branch, true, true, false, false, hgtdtimes_resos, hgtdtimes_purity,
					    fjet_idealres, ftrack_idealres, pu_frac_idealres,
					    hs_track_idealres, pu_track_idealres, recovtx_z_idealres, hs_pu_hgtdtimes);
    // if (hs_track_idealeff.hist.at(ScoreType::PASS)->Integral() < 10e3) // break;
    auto pass_idealeff = process_event_data(&branch, true, false, false, false, hgtdtimes_resos, hgtdtimes_purity,
					    fjet_idealeff, ftrack_idealeff, pu_frac_idealeff,
					    hs_track_idealeff, pu_track_idealeff, recovtx_z_idealeff, hs_pu_hgtdtimes);
    // if (hs_track_hgtdtimes.hist.at(ScoreType::PASS)->Integral() >= 10e3 and
    // 	hs_track_idealres.hist.at(ScoreType::PASS)->Integral() >= 10e3 and
    // 	hs_track_idealeff.hist.at(ScoreType::PASS)->Integral() >= 10e3)
    //   break;
    // if (pass_hgtdtimes.first == 0 and pass_idealres.first == 1) {
    //   TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    //   TString file =  fileName(46,6);
    //   smear_improves.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0 --extra_time %.2f",
    // 				    file.Data(), this_evnt, pass_hgtdtimes.second));
    // }

    if (pass_hgtdtimes.first == 1) {
      TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
      TString file =  fileName(46,6);
      smear_degrades.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0 --extra_time %.2f",
				    file.Data(), this_evnt, pass_hgtdtimes.second));
    }

    // if (pass_hgtdtimes.first == 0 and pass_idealres.first == 0) {
    //   TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    //   TString file =  fileName(46,6);
    //   both_fail.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0 --extra_time %.2f",
    // 				    file.Data(), this_evnt, pass_hgtdtimes.second));
    // }

    // if (pass_hgtdtimes.first == 1 and pass_idealres.first == 1) {
    //   TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    //   TString file =  fileName(46,6);
    //   both_pass.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0 --extra_time %.2f",
    // 				    file.Data(), this_evnt, pass_hgtdtimes.second));
    // }
      
  }
  
  // add all objects to plotobjs before drawing everything
  std::vector<PlotObj*> plots = {&fjet_hgtdtimes, &ftrack_hgtdtimes, &hs_track_hgtdtimes, &pu_track_hgtdtimes, &pu_frac_hgtdtimes,
				 &fjet_idealres, &ftrack_idealres, &hs_track_idealres, &pu_track_idealres, &pu_frac_idealres,
				 &fjet_idealeff, &ftrack_idealeff, &hs_track_idealeff, &pu_track_idealeff, &pu_frac_idealeff,};
				 // &fjet_hgtdtimes_usez, &ftrack_hgtdtimes_usez, &hs_track_hgtdtimes_usez, &pu_track_hgtdtimes_usez, &pu_frac_hgtdtimes_usez};

  for (auto& plot: plots)
    plot->plotPostprocessing();
  std::cout << "FINISHED CREATING " << std::endl;

  std::cout << "--- EFFICIENCIES FOR TRACKPT SCORE ---" << std::endl;
  hs_track_hgtdtimes.printEfficiencyStats(ScoreType::TRKPT);
  // hs_track_idealres.printEfficiencyStats(ScoreType::TRKPT);
  // hs_track_idealeff.printEfficiencyStats(ScoreType::TRKPT);

  std::cout << "--- EFFICIENCIES FOR TRACKPTZ SCORE ---" << std::endl;
  hs_track_hgtdtimes.printEfficiencyStats(ScoreType::TRKPTZ);
  // hs_track_idealres.printEfficiencyStats(ScoreType::TRKPTZ);
  // hs_track_idealeff.printEfficiencyStats(ScoreType::TRKPTZ);

  std::cout << "--- EFFICIENCIES FOR IDEAL SCORE ---" << std::endl;
  pu_frac_hgtdtimes.printEfficiencyStats(ScoreType::PASS);
  pu_frac_idealres.printEfficiencyStats(ScoreType::PASS);
  pu_frac_idealeff.printEfficiencyStats(ScoreType::PASS);
  

  for (auto& plot: plots)
    plot->plotLogic(canvas);
  std::cout << "FINISHED PRINTING BULK " << std::endl;

  plot_money(Form("figs/moneyplot_hstrack.pdf"), canvas,
	     prep_for_money(
			    {&hs_track_hgtdtimes, &hs_track_idealres, &hs_track_idealeff}, 
			    {{HGTD, TRKPTZ}, {TRKPTZ}, {TRKPTZ, PASS}}
			    ));

  plot_money(Form("figs/moneyplot_nftrack.pdf"), canvas,
	     prep_for_money(
			    {&ftrack_hgtdtimes, &ftrack_idealres, &ftrack_idealeff}, 
			    {{HGTD, TRKPTZ}, {TRKPTZ}, {TRKPTZ, PASS}}
			    ));

  plot_money(Form("figs/moneyplot_putrack.pdf"), canvas,
	     prep_for_money(
			    {&pu_track_hgtdtimes, &pu_track_idealres, &pu_track_idealeff}, 
			    {{HGTD, TRKPTZ}, {TRKPTZ}, {TRKPTZ, PASS}}
			    ));

  plot_money(Form("figs/moneyplot_pufrac.pdf"), canvas,
	     prep_for_money(
			    {&pu_frac_hgtdtimes, &pu_frac_idealres, &pu_frac_idealeff}, 
			    {{HGTD, TRKPTZ}, {TRKPTZ}, {TRKPTZ, PASS}}
			    ));

  plot_money(Form("figs/moneyplot_fjet.pdf"), canvas,
	     prep_for_money(
			    {&fjet_hgtdtimes, &fjet_idealres, &fjet_idealeff}, 
			    {{HGTD, TRKPTZ}, {TRKPTZ}, {TRKPTZ, PASS}}
			    ));

  std::cout << "ALL DONE :3" << std::endl;

  std::cout << smear_degrades.size() << " PASS EVENTS IN BAD RES REGION: " << std::endl;
  for (auto s: smear_degrades)
    std::cout << s << std::endl;
  
  std::cout << smear_improves.size() << " IMPROVED EVENTS: " << std::endl;
  for (auto s: smear_improves)
    std::cout << s << std::endl;
}
