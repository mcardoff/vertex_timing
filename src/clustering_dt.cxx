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
  std::map<Score, AnalysisObj> hgtdtimesMap;
  hgtdtimesMap.emplace(Score::HGTD  , AnalysisObj("hgtdtimes", "HGTD Times", Score::HGTD)  );
  hgtdtimesMap.emplace(Score::TRKPTZ, AnalysisObj("hgtdtimes", "HGTD Times", Score::TRKPTZ));
  hgtdtimesMap.emplace(Score::CALO60, AnalysisObj("hgtdtimes", "HGTD Times", Score::CALO60));
  hgtdtimesMap.emplace(Score::CALO90, AnalysisObj("hgtdtimes", "HGTD Times", Score::CALO90));
  hgtdtimesMap.emplace(Score::JUST60, AnalysisObj("hgtdtimes", "HGTD Times", Score::JUST60));
  hgtdtimesMap.emplace(Score::JUST90, AnalysisObj("hgtdtimes", "HGTD Times", Score::JUST90));
  // hgtdtimes_map.emplace(PASS  , AnalysisObj("hgtdtimes", "HGTD Times", PASS)  );

  std::map<Score, AnalysisObj> idealresMap;
  idealresMap.emplace(Score::TRKPTZ, AnalysisObj("idealres", "Ideal Res. HGTD", Score::TRKPTZ));
  idealresMap.emplace(Score::CALO60, AnalysisObj("idealres", "Ideal Res. HGTD", Score::CALO60));
  idealresMap.emplace(Score::CALO90, AnalysisObj("idealres", "Ideal Res. HGTD", Score::CALO90));
  // idealresMap.emplace(PASS  , AnalysisObj("idealres", "Ideal Res. HGTD", PASS)  );

  std::map<Score, AnalysisObj> idealeffMap;
  idealeffMap.emplace(Score::TRKPTZ, AnalysisObj("idealeff", "Ideal Eff.+Res. HGTD", Score::TRKPTZ));
  idealeffMap.emplace(Score::CALO90, AnalysisObj("idealeff", "Ideal Eff.+Res. HGTD", Score::CALO90));
  idealeffMap.emplace(Score::CALO60, AnalysisObj("idealeff", "Ideal Eff.+Res. HGTD", Score::CALO60));
  // idealeffMap.emplace(PASS  , AnalysisObj("idealeff", "Ideal Eff.+Res. HGTD", PASS)  );
    
  std::vector<TString> passMineFailHgtd, passIdealresFailHgtd, passIdealeffFailIdealres;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop\n";
  bool progress = true;
  TString fileName, file, eventdisplay="python event_display.py --file_num %s --event_num %lld --extra_time %.2f";
  Long64_t nEvent = chain.GetEntries();
  Long64_t evtMax = nEvent+1;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t thisEvnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readNum = chain.GetReadEntry()+1;
    if (progress && readNum % 1000 == 0)
      std::cout << "Progress: " << readNum << "/" << nEvent << "\n";

    // bool
    //   a = hgtdtimesMap.at(Score::CALO60).get("pu_frac")->effPass->Integral() < evtMax,
    //   b = idealresMap.at(Score::CALO60).get("pu_frac")->effPass->Integral() < evtMax,
    //   c = idealeffMap.at(Score::CALO60).get("pu_frac")->effPass->Integral() < evtMax;
    
    // if ((!a) and (!b) and (!c)) break;

    auto passHgtd = processEventData(&branch, false, true, false, hgtdtimesMap);
    auto passIdealres = processEventData(&branch, true, true, false, idealresMap);
    auto passIdealeff = processEventData(&branch, true, false, false, idealeffMap);

    // fileName = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    // file = fileName(49,6);

    // // Events where my algorithm passes when HGTD fails
    // if (pass_hgtd.first == 2)
    //  passMineFailHgtd.push_back(Form(eventdisplay, file.Data(), this_evnt, pass_hgtd.second));
      
    // // Events where Ideal res passes where HGTD Times Fail
    // if (pass_hgtd.first == 0 and pass_idealres.first >= 1)
    //  passIdealresFailHgtd.push_back(Form(eventdisplay, file.Data(), this_evnt, pass_idealres.second));

    // // Events where Ideal Eff+Res passes where Ideal Res Fails:
    // if (pass_hgtd.first == 0 and pass_idealres.first == 0 and pass_idealeff.first >= 1)
    //   passIdealeffFailIdealres.push_back(Form(eventdisplay, file.Data(), this_evnt, pass_idealeff.second));
  }

  auto maps = {&hgtdtimesMap, &idealresMap, &idealeffMap};

  for (auto &map: maps) {
    for (auto &[k, analysis] : *map) {
      analysis.postProcessing();
    }
  }
  // for (auto &map: maps)
  //   for (auto &[k, analysis] : *map)
  //     analysis["hs_track"]->printEfficiencyStats();

  std::cout << "FINISHED PROCESSING\n";

  for (const auto *const KEY: {"fjet", "ftrack", "pu_frac", "hs_track", "pu_track"}) {
    // moneyPlot(Form("../figs/fullcomparison_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealresMap.at(TRKPTZ),
    // 	       &idealeffMap.at(TRKPTZ),
    // 	       &idealeffMap.at(PASS),
    // 	       &hgtdtimesMap.at(CALO90),
    // 	      });

    // moneyPlot(Form("../figs/hgtd_trkptz_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	      });

    // moneyPlot(Form("../figs/idealres_calocomp60_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealresMap.at(TRKPTZ),
    // 	       &hgtdtimesMap.at(CALO60),
    // 	       &idealresMap.at(CALO60),
    // 	      });

    // moneyPlot(Form("../figs/idealres_calocomp90_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealresMap.at(TRKPTZ),
    // 	       &hgtdtimesMap.at(CALO90),
    // 	       &idealresMap.at(CALO90),
    // 	      });

    moneyPlot(Form("../figs/idealeff_calocomp60_%s.pdf", KEY), KEY, canvas,
	      {&hgtdtimesMap.at(HGTD),
	       &hgtdtimesMap.at(TRKPTZ),
	       &idealeffMap.at(TRKPTZ),
	       &hgtdtimesMap.at(CALO60),
	       &idealeffMap.at(CALO60),
	      });

    // moneyPlot(Form("../figs/idealeff_calocomp90_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealeffMap.at(TRKPTZ),
    // 	       &hgtdtimesMap.at(CALO90),
    // 	       &idealeffMap.at(CALO90),
    // 	      });

    // moneyPlot(Form("../figs/idealres_calocomp_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealresMap.at(TRKPTZ),
    // 	       &idealresMap.at(CALO60),
    // 	       &idealresMap.at(CALO90),
    // 	      });

    // moneyPlot(Form("../figs/idealeff_calocomp_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealeffMap.at(TRKPTZ),
    // 	       &idealeffMap.at(CALO60),
    // 	       &idealeffMap.at(CALO90),
    // 	      });
    // std::cout << "MONEY\n";
    // moneyPlot(Form("../figs/justcalo_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(Score::HGTD),
    // 	       &hgtdtimesMap.at(Score::JUST60),
    // 	       &hgtdtimesMap.at(Score::JUST90),
    // 	      });

    // moneyPlot(Form("../figs/idealres_comp_%s.pdf", KEY), KEY, canvas,
    // 	      {&hgtdtimesMap.at(HGTD),
    // 	       &hgtdtimesMap.at(TRKPTZ),
    // 	       &idealresMap.at(TRKPTZ),
    // 	       // &idealresMap.at(PASS),
    // 	      });

  }

  // inclusivePlot(Form("../figs/inclusivereso_logscale.pdf"), true, true,
  // 		-400, 400, canvas,
  // 		{&hgtdtimesMap.at(HGTD), &hgtdtimesMap.at(TRKPTZ),
  // 		 &idealresMap.at(TRKPTZ),
  // 		 &idealeffMap.at(TRKPTZ), &idealeffMap.at(PASS)});

  // inclusivePlot(Form("../figs/inclusivereso_linscale.pdf"), false, true,
  // 		-200, 200, canvas,
  // 		{&hgtdtimesMap.at(HGTD), &hgtdtimesMap.at(TRKPTZ),
  // 		 &idealresMap.at(TRKPTZ),
  // 		 &idealeffMap.at(TRKPTZ), &idealeffMap.at(PASS)});
  
  std::cout << "FINISHED PLOT PRINTING\n";

  std::cout << passMineFailHgtd.size() << " events where i'm better than HGTD\n";
  for (const auto& str: passMineFailHgtd) { std::cout << str << '\n'; }
  
  std::cout << passIdealresFailHgtd.size() << " events where ideal res. times better than hgtd times\n";
  for (const auto& str: passIdealresFailHgtd) { std::cout << str << '\n'; }
  
  std::cout << passIdealeffFailIdealres.size() << " events where ideal res.+eff. times better than ideal res hgtd times\n";
  for (const auto& str: passIdealeffFailIdealres) { std::cout << str << '\n'; }
  
  return 0;
}
