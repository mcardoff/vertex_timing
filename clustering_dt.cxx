#include "event_processing.h"
#include "plotting_utilities.h"

using namespace myutl;

void clustering_dt() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

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
  idealeff_map.emplace(TRKPTZ, AnalysisObj("idealeff", "Ideal Eff. HGTD", TRKPTZ));
  idealeff_map.emplace(PASS  , AnalysisObj("idealeff", "Ideal Eff. HGTD", PASS)  );
    
  std::vector<TString> pass_fail_hgtd;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readnum = chain.GetReadEntry()+1;
    if (progress && readnum % 1000 == 0) std::cout << "Progress: " << readnum << "/" << chain.GetEntries() << "\n";

    // if (readnum > 10e3) break;
    // bool
    //   a = hgtdtimes_map.at(ScoreType::PASS).get("pu_frac")->eff_pass->Integral() < 10e3,
    //   b = idealres_map.at(ScoreType::PASS).get("pu_frac")->eff_pass->Integral() < 10e3,
    //   c = idealeff_map.at(ScoreType::PASS).get("pu_frac")->eff_pass->Integral() < 10e3;
    

    auto pass_hgtd = process_event_data(&branch, false, true, false, hgtdtimes_map);
    // if (pass_hgtd.first == 1) {
    //   TString fileName   = (branch.reader.GetTree()->GetCurrentFile()->GetName());
    //   TString file =  fileName(46,6);
    //   pass_fail_hgtd.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0 --extra_time %.2f",
    // 				    file.Data(), this_evnt, pass_hgtd.second));
    // }
    auto pass_idealres = process_event_data(&branch, true, true, false, idealres_map);
    auto pass_idealeff = process_event_data(&branch, true, false, false, idealeff_map);
  }

  auto maps = {&hgtdtimes_map, &idealres_map, &idealeff_map};

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis.postProcessing();

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis["hs_track"]->printEfficiencyStats();

  std::cout << "FINISHED PROCESSING" << std::endl;

  for (const auto k: {"fjet", "ftrack", "pu_frac", "hs_track", "pu_track"})
    moneyPlot(Form("figs/moneyplot_neweff_%s.pdf", k), k, canvas,
	      {&hgtdtimes_map.at(HGTD),
	       &hgtdtimes_map.at(TRKPTZ),
	       &idealres_map.at(TRKPTZ),
	       &idealeff_map.at(TRKPTZ),
	       &idealeff_map.at(PASS),
	      });

  std::cout << "FINISHED PLOT PRINTING" << std::endl;

  // std::cout << pass_fail_hgtd.size() << " INTERESTING EVENTS" << std::endl;
  // for (const auto& s: pass_fail_hgtd)
  //   std::cout << s << std::endl;
  
}
