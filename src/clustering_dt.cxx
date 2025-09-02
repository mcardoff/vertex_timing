#include <TROOT.h>
#include "event_processing.h"
#include "plotting_utilities.h"

using namespace myutl;

int main() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  ROOT::EnableImplicitMT(); // uses all cores

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  // HGTD Times
  std::map<ScoreType, AnalysisObj> hgtdtimes_map;
  hgtdtimes_map.emplace(HGTD  , AnalysisObj("hgtdtimes", "HGTD Times", HGTD)  );
  hgtdtimes_map.emplace(TRKPTZ, AnalysisObj("hgtdtimes", "HGTD Times", TRKPTZ));
  hgtdtimes_map.emplace(PASS  , AnalysisObj("hgtdtimes", "HGTD Times", PASS)  );

  std::map<ScoreType, AnalysisObj> idealres_map;
  idealres_map.emplace(TRKPTZ, AnalysisObj("idealres", "Ideal Res. HGTD", TRKPTZ));
  idealres_map.emplace(PASS  , AnalysisObj("idealres", "Ideal Res. HGTD", PASS)  );

  std::map<ScoreType, AnalysisObj> idealeff_map;
  idealeff_map.emplace(TRKPTZ, AnalysisObj("idealeff", "Ideal Eff.+Res. HGTD", TRKPTZ));
  idealeff_map.emplace(PASS  , AnalysisObj("idealeff", "Ideal Eff.+Res. HGTD", PASS)  );
    
  std::vector<TString> pass_mine_fail_hgtd, pass_idealres_fail_hgtd, pass_idealeff_fail_idealres;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  int eMax = 10e6;
  TString fileName, file, eventdisplay="python event_display.py --file_num %s --event_num %lld --extra_time %.2f";
  Long64_t nevent = chain.GetEntries();
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readnum = chain.GetReadEntry()+1;
    if (progress && readnum % 1000 == 0)
      std::cout << "Progress: " << readnum << "/" << nevent << "\n";

    // bool
    //   a = hgtdtimes_map.at(ScoreType::PASS).get("pu_frac")->eff_pass->Integral() < eMax,
    //   b = idealres_map.at(ScoreType::PASS).get("pu_frac")->eff_pass->Integral() < eMax,
    //   c = idealeff_map.at(ScoreType::PASS).get("pu_frac")->eff_pass->Integral() < eMax;
    
    // if ((!a) and (!b) and (!c)) break;

    auto pass_hgtd = process_event_data(&branch, false, true, false, hgtdtimes_map);
    auto pass_idealres = process_event_data(&branch, true, true, false, idealres_map);
    auto pass_idealeff = process_event_data(&branch, true, false, false, idealeff_map);

    // fileName = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    // file = fileName(49,6);

    // // Events where my algorithm passes when HGTD fails
    // if (pass_hgtd.first == 2)
    //   pass_mine_fail_hgtd.push_back(Form(eventdisplay, file.Data(), this_evnt, pass_hgtd.second));
      
    // // Events where Ideal res passes where HGTD Times Fail
    // if (pass_hgtd.first == 0 and pass_idealres.first >= 1)
    //   pass_idealres_fail_hgtd.push_back(Form(eventdisplay, file.Data(), this_evnt, pass_idealres.second));

    // // Events where Ideal Eff+Res passes where Ideal Res Fails:
    // if (pass_hgtd.first == 0 and pass_idealres.first == 0 and pass_idealeff.first >= 1)
    //   pass_idealeff_fail_idealres.push_back(Form(eventdisplay, file.Data(), this_evnt, pass_idealeff.second));
  }

  auto maps = {&hgtdtimes_map, &idealres_map, &idealeff_map};

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis.postProcessing();

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis["hs_track"]->printEfficiencyStats();

  std::cout << "FINISHED PROCESSING" << std::endl;

  for (const auto k: {"fjet", "ftrack", "pu_frac", "hs_track", "pu_track"}) {
    moneyPlot(Form("../figs/fullcomparison_%s.pdf", k), k, canvas,
	      {&hgtdtimes_map.at(HGTD),
	       &hgtdtimes_map.at(TRKPTZ),
	       &idealres_map.at(TRKPTZ),
	       &idealeff_map.at(TRKPTZ),
	       &idealeff_map.at(PASS),
	      });

    moneyPlot(Form("../figs/hgtd_trkptz_%s.pdf", k), k, canvas,
	      {&hgtdtimes_map.at(HGTD),
	       &hgtdtimes_map.at(TRKPTZ),
	       // &hgtdtimes_map.at(PASS),
	      });

    moneyPlot(Form("../figs/idealres_comp_%s.pdf", k), k, canvas,
	      {&hgtdtimes_map.at(HGTD),
	       &hgtdtimes_map.at(TRKPTZ),
	       &idealres_map.at(TRKPTZ),
	       // &idealres_map.at(PASS),
	      });

  }

  inclusivePlot(Form("../figs/inclusivereso_logscale.pdf"), true, true,
		-400, 400, canvas,
		{&hgtdtimes_map.at(HGTD), &hgtdtimes_map.at(TRKPTZ),
		 &idealres_map.at(TRKPTZ),
		 &idealeff_map.at(TRKPTZ), &idealeff_map.at(PASS)});

  inclusivePlot(Form("../figs/inclusivereso_linscale.pdf"), false, true,
		-200, 200, canvas,
		{&hgtdtimes_map.at(HGTD), &hgtdtimes_map.at(TRKPTZ),
		 &idealres_map.at(TRKPTZ),
		 &idealeff_map.at(TRKPTZ), &idealeff_map.at(PASS)});
  
  std::cout << "FINISHED PLOT PRINTING" << std::endl;

  std::cout << pass_mine_fail_hgtd.size() << " events where i'm better than HGTD" << std::endl;
  for (auto s: pass_mine_fail_hgtd)
    std::cout << s << std::endl;
  
  std::cout << pass_idealres_fail_hgtd.size() << " events where ideal res. times better than hgtd times" << std::endl;
  for (auto s: pass_idealres_fail_hgtd)
    std::cout << s << std::endl;
  
  std::cout << pass_idealeff_fail_idealres.size() << " events where ideal res.+eff. times better than ideal res hgtd times" << std::endl;
  for (auto s: pass_idealeff_fail_idealres)
    std::cout << s << std::endl;
  
  return 0;
}
