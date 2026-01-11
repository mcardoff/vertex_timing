#include <RtypesCore.h>
#include <TROOT.h>
#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"

using namespace MyUtl;

auto main() -> int {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  ROOT::EnableImplicitMT(); // uses all cores

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  // HGTD Times
  std::map<Score, AnalysisObj> mapHGTD;
  mapHGTD.emplace(HGTD,   AnalysisObj("HGTD Times", HGTD   ));
  mapHGTD.emplace(TRKPTZ, AnalysisObj("HGTD Times", TRKPTZ ));
  // mapHGTD.emplace(CALO60, AnalysisObj("HGTD Times", CALO60 ));
  // mapHGTD.emplace(CALO90, AnalysisObj("HGTD Times", CALO90 ));
  // mapHGTD.emplace(JUST60, AnalysisObj("HGTD Times", JUST60 ));
  // mapHGTD.emplace(JUST90, AnalysisObj("HGTD Times", JUST90 ));
  // mapHGTD.emplace(FILT60, AnalysisObj("HGTD Times", FILT60 ));
  // mapHGTD.emplace(FILT90, AnalysisObj("HGTD Times", FILT90 ));
  mapHGTD.emplace(FILTJET,AnalysisObj("HGTD Times", FILTJET));
  mapHGTD.emplace(PASS  , AnalysisObj("HGTD Times", PASS   ));

  std::map<Score, AnalysisObj> mapIdealRes;
  mapIdealRes.emplace(TRKPTZ , AnalysisObj("Ideal Res. HGTD", TRKPTZ ));
  // mapIdealRes.emplace(CALO60 , AnalysisObj("Ideal Res. HGTD", CALO60 ));
  // mapIdealRes.emplace(CALO90 , AnalysisObj("Ideal Res. HGTD", CALO90 ));
  // mapIdealRes.emplace(JUST60 , AnalysisObj("Ideal Res. HGTD", JUST60 ));
  // mapIdealRes.emplace(JUST90 , AnalysisObj("Ideal Res. HGTD", JUST90 ));
  // mapIdealRes.emplace(FILT60 , AnalysisObj("Ideal Res. HGTD", FILT60 ));
  // mapIdealRes.emplace(FILT90 , AnalysisObj("Ideal Res. HGTD", FILT90 ));
  mapIdealRes.emplace(FILTJET, AnalysisObj("Ideal Res. HGTD", FILTJET));
  mapIdealRes.emplace(PASS   , AnalysisObj("Ideal Res. HGTD", PASS   ));

  std::map<Score, AnalysisObj> mapIdealEff;
  mapIdealEff.emplace(TRKPTZ , AnalysisObj("Ideal Res.+Eff. HGTD", TRKPTZ ));
  // mapIdealEff.emplace(CALO60 , AnalysisObj("Ideal Res.+Eff. HGTD", CALO60 ));
  // mapIdealEff.emplace(CALO90 , AnalysisObj("Ideal Res.+Eff. HGTD", CALO90 ));
  // mapIdealEff.emplace(JUST60 , AnalysisObj("Ideal Res.+Eff. HGTD", JUST60 ));
  // mapIdealEff.emplace(JUST90 , AnalysisObj("Ideal Res.+Eff. HGTD", JUST90 ));
  // mapIdealEff.emplace(FILT60 , AnalysisObj("Ideal Res.+Eff. HGTD", FILT60 ));
  // mapIdealEff.emplace(FILT90 , AnalysisObj("Ideal Res.+Eff. HGTD", FILT90 ));
  mapIdealEff.emplace(FILTJET, AnalysisObj("Ideal Res.+Eff. HGTD", FILTJET));
  mapIdealEff.emplace(PASS   , AnalysisObj("Ideal Res.+Eff. HGTD", PASS   ));

  std::vector<TString> hasPassingClusterButFailsHGTD, hasPassingClusterButFailsIdealRes, hasPassingClusterButFailsIdealEff;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop\n";
  bool progress = true;
  TString fileName, file, eventdisplay="python event_display.py --file_num %s --event_num %lld --extra_time %.2f";
  Long64_t nEvent = chain.GetEntries();
  Long64_t evtMax = 10e3;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t thisEvnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readNum = chain.GetReadEntry()+1;
    if (progress && readNum % 1000 == 0)
      std::cout << "Progress: " << readNum << "/" << nEvent << "\n";

    // bool
    //   a = mapHGTD.at(Score::CALO60).get("pu_frac")->effPass->Integral() < evtMax,
    //   b = mapIdealRes.at(Score::CALO60).get("pu_frac")->effPass->Integral() < evtMax,
    //   c = mapIdealEff.at(Score::CALO60).get("pu_frac")->effPass->Integral() < evtMax;
    
    // if ((!a) and (!b) and (!c)) break;

    auto passHGTD     = processEventData(&branch, false, true , false, mapHGTD);
    auto passIdealRes = processEventData(&branch, true , true , false, mapIdealRes);
    auto passIdealEff = processEventData(&branch, true , false, false, mapIdealEff);

    fileName = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    file = fileName(49,6);

    // Events where my algorithm passes when HGTD fails
    if (passHGTD.first == 2)
     hasPassingClusterButFailsHGTD.push_back(Form(eventdisplay, file.Data(), thisEvnt, passHGTD.second));
      
    // Events where Ideal res passes where HGTD Times Fail
    if (passIdealRes.first == 2)
     hasPassingClusterButFailsIdealRes.push_back(Form(eventdisplay, file.Data(), thisEvnt, passIdealRes.second));

    // Events where Ideal Eff+Res passes where Ideal Res Fails:
    if (passIdealEff.first == 2)
      hasPassingClusterButFailsIdealEff.push_back(Form(eventdisplay, file.Data(), thisEvnt, passIdealEff.second));
  }

  auto maps = {&mapHGTD, &mapIdealRes, &mapIdealEff};
  // auto maps = {&mapHGTD};

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis.postProcessing();

      
  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis["hs_track"]->printEfficiencyStats();

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis.fullPlotting(canvas);

  std::cout << "FINISHED PROCESSING\n";

  for (const auto *const KEY: {"pu_frac", "ftrack", "pu_track", "fjet", "hs_track", "vtx_dz"}) {

    // Comparison between my algo and HGTD algo
    std::cout << "HGTD TRKPTZ" << std::endl;
    moneyPlot(Form("../figs/hgtd_trkptz_%s.pdf", KEY), KEY, canvas,
	      { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ) });

    // HGTD Only
    // std::cout << "HGTD" << std::endl;
    // moneyPlot(Form("../figs/hgtd_%s.pdf", KEY), KEY, canvas,
	      // { &mapHGTD.at(HGTD)});

    // Comparing those same two with calo exclusion
    // std::cout << "CALO EXCLUSION COMP" << std::endl;
    // moneyPlot(Form("../figs/calo_comp_%s.pdf", KEY), KEY, canvas,
    // 	      { &mapHGTD.at(HGTD)  , &mapHGTD.at(TRKPTZ),
    // 		&mapHGTD.at(CALO60), &mapHGTD.at(CALO90), });

    // Filtering out various types of track
    // std::cout << "FILTER COMP" << std::endl;
    // moneyPlot(Form("../figs/filter_comp_%s.pdf", KEY), KEY, canvas,
    // 	      { &mapHGTD.at(HGTD)  ,
    // 		&mapHGTD.at(TRKPTZ),
    // 		&mapHGTD.at(FILT60),
    // 		&mapHGTD.at(FILT90),
    // 		&mapHGTD.at(FILTJET) });

    // Compare calo exclusion with only calo time and filtering (60ps)
    // std::cout << "60ps CASES" << std::endl;
    // moneyPlot(Form("../figs/calo_comp_60_%s.pdf", KEY), KEY, canvas,
    // 	      { &mapHGTD.at(HGTD)  ,
    // 		&mapHGTD.at(TRKPTZ),
    // 	        &mapHGTD.at(CALO60),
    // 		&mapHGTD.at(JUST60),
    // 		&mapHGTD.at(FILT60), });

    // Compare calo exclusion with only calo time and filtering (90ps)
    // std::cout << "90ps CASES" << std::endl;
    // moneyPlot(Form("../figs/calo_comp_90_%s.pdf", KEY), KEY, canvas,
    // 	      { &mapHGTD.at(HGTD)  ,
    // 		&mapHGTD.at(TRKPTZ),
    // 	        &mapHGTD.at(CALO90),
    // 		&mapHGTD.at(JUST90),
    // 		&mapHGTD.at(FILT90), });

    // moneyPlot(Form("../figs/filt_pass_comp_%s.pdf", KEY), KEY, canvas,
    // 	      { &mapHGTD.at(HGTD)  , &mapHGTD.at(TRKPTZ) ,
    // 		&mapHGTD.at(FILT60), &mapHGTD.at(FILTJET),
    // 		&mapHGTD.at(PASS) });

    moneyPlot(Form("../figs/ideal_comp_%s.pdf", KEY), KEY, canvas,
	      { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ), //&mapIdealEff.at(PASS),
	      });

    moneyPlot(Form("../figs/pass_comp_%s.pdf", KEY), KEY, canvas,
	      { &mapHGTD.at(HGTD), &mapHGTD.at(PASS), &mapIdealRes.at(PASS), &mapIdealEff.at(PASS),
	      });

    // moneyPlot(Form("../figs/filt_jet_comp_%s.pdf", KEY), KEY, canvas,
    // 	      { &mapHGTD.at(HGTD)    ,
    // 		&mapHGTD.at(FILTJET) , &mapIdealRes.at(FILTJET), &mapIdealEff.at(FILTJET),
    // 	        &mapIdealEff.at(PASS),
    // 	      });
  }

  inclusivePlot(Form("../figs/inclusivereso_logscale.pdf"), true, false,
		-400, 400, canvas,
	        {
		  &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ), &mapIdealEff.at(PASS),
		});

  inclusivePlot(Form("../figs/inclusivereso_linscale.pdf"), false, false,
		-200, 200, canvas,
	        {
		  &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ), &mapIdealEff.at(PASS),
		});
  
  std::cout << "FINISHED PLOT PRINTING\n";

  // std::cout << (double)mapHGTD.at(TRKPTZ)["vtx_dz"]->effTotal->Integral(2,200)/(double)mapHGTD.at(TRKPTZ)["vtx_dz"]->effTotal->Integral() << std::endl;

  // std::cout << hasPassingClusterButFailsHGTD.size() << " events have a passing cluster but TRKPTZ fails with HGTD times\n";
  // for (const auto& str: hasPassingClusterButFailsHGTD) { std::cout << str << '\n'; }
  
  // std::cout << hasPassingClusterButFailsIdealRes.size() << " events have a passing cluster but TRKPTZ fails with Ideal Res. times\n";
  // for (const auto& str: hasPassingClusterButFailsIdealRes) { std::cout << str << '\n'; }
  
  // std::cout << hasPassingClusterButFailsIdealEff.size() << " events have a passing cluster but TRKPTZ fails with Ideal Res.+Eff. times\n";
  // for (const auto& str: hasPassingClusterButFailsIdealEff) { std::cout << str << '\n'; }
  
  return 0;
}
