#include <RtypesCore.h>
#include <TROOT.h>
#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

// Output directory for saved plots
static constexpr const char* SAVE_DIR = "../figs";

// Event display command template (filled with file number, event number, extra time)
static constexpr const char* EVTDISPLAY_FMT =
  "python3 event_display.py --file_num %s --event_num %lld --extra_time %.2f";

// Set to true to print event display commands to stdout after the event loop.
static constexpr bool PRINT_EVENT_DISPLAYS = false;

// ---------------------------------------------------------------------------
// Helper: build one analysis map for a given scenario label.
//   The active scores are listed explicitly so adding/removing a score
//   only requires changing this one function.
// ---------------------------------------------------------------------------

enum class Scenario { HGTD, IDEAL_RES, IDEAL_EFF };

auto buildAnalysisMap(
   Scenario scenario
) -> std::map<Score, AnalysisObj> {
  const char* label = [&]() -> const char* {
    switch (scenario) {
      case Scenario::HGTD:      return "HGTD Times";
      case Scenario::IDEAL_RES: return "Ideal Res.";
      case Scenario::IDEAL_EFF: return "Ideal Res.+Eff.";
    }
    return "";
  }();

  std::map<Score, AnalysisObj> m;

  // Scores active in all scenarios
  m.emplace(Score::TRKPTZ,     AnalysisObj(label, Score::TRKPTZ    ));
  // m.emplace(Score::TESTML,     AnalysisObj(label, Score::TESTML    ));
  // m.emplace(Score::PASS,       AnalysisObj(label, Score::PASS      ));
  m.emplace(Score::TEST_MISCL, AnalysisObj(label, Score::TEST_MISCL));
  // m.emplace(Score::TEST_MISAS, AnalysisObj(label, Score::TEST_MISAS ));
  // m.emplace(Score::REFINED,    AnalysisObj(label, Score::REFINED    ));

  // Scores active only in the real-HGTD scenario
  if (scenario == Scenario::HGTD) {
    m.emplace(Score::HGTD,       AnalysisObj(label, Score::HGTD      ));
    // m.emplace(Score::TEST_HS,    AnalysisObj(label, Score::TEST_HS    ));
    // m.emplace(Score::HGTD_SORT, AnalysisObj(label, Score::HGTD_SORT));
    // m.emplace(Score::ITERATIVE,  AnalysisObj(label, Score::ITERATIVE ));
    // m.emplace(Score::CONE_BDT,  AnalysisObj(label, Score::CONE_BDT ));
    // m.emplace(Score::FILTJET,   AnalysisObj(label, Score::FILTJET  ));
  }

  return m;
}

// ---------------------------------------------------------------------------
// Helper: add an event-display command string to a list if condition is met.
// ---------------------------------------------------------------------------

void collectEventDisplay(
  std::vector<TString>& list,
  int returnCode,
  const std::pair<int,double>& result,
  const TString& fileNum,
  Long64_t eventNum
) {
  if (result.first == returnCode)
    list.push_back(Form(EVTDISPLAY_FMT, fileNum.Data(), eventNum, result.second));
}

// ---------------------------------------------------------------------------
// Helper: print all collected event display commands for one scenario.
//   Output is gated on PRINT_EVENT_DISPLAYS; set that flag to true at the
//   top of this file to enable printing after the event loop.
// ---------------------------------------------------------------------------

void printEventDisplays(
  const char* label,
  const std::vector<TString>& list
) {
  if (!PRINT_EVENT_DISPLAYS) return;
  std::cout << "\n--- Event displays: " << label << " ("
	    << list.size() << " events) ---\n";
  for (const auto& cmd : list)
    std::cout << cmd << '\n';
}

// ---------------------------------------------------------------------------
// Helper: generate all per-KEY comparison plots for one variable KEY.
// ---------------------------------------------------------------------------

void makeComparisonPlots(
  const char* key,
  TCanvas* canvas,
  std::map<Score, AnalysisObj>& mapHGTD,
  std::map<Score, AnalysisObj>& mapIdealRes,
  std::map<Score, AnalysisObj>& mapIdealEff
) {
  static const std::string compSubdir = std::string(SAVE_DIR) + "/comparisons";

  moneyPlot(Form("%s/theplot_%s.pdf",compSubdir.c_str(), key), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),           // HGTD Base Performance
	      &mapHGTD.at(Score::TRKPTZ),         // TRKPTZ Base Performance
	      &mapHGTD.at(Score::TEST_MISCL),     // TRKPTZ Fix Misclustering
	      &mapIdealRes.at(Score::TRKPTZ),     // TRKPTZ Fix Misclustering
	      &mapIdealRes.at(Score::TEST_MISCL), // TRKPTZ Fix Misclustering
	    });
  
  // HGTD algo vs TRKPTZ
  moneyPlot(Form("%s/hgtd_trkptz_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ)
	    });

  // HGTD algo vs TRKPTZ vs REFINED (cone + iterative refinement)
  // moneyPlot(Form("%s/refined_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::REFINED)
  // 	    });

  // // HGTD algo vs TRKPTZ vs HGTD_SORT (pT-sorted simultaneous + BDT)
  // moneyPlot(Form("%s/hgtd_sort_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::HGTD_SORT)
  // 	    });

  // HGTD algo vs TRKPTZ vs HGTD_SORT (pT-sorted simultaneous + BDT)
  // moneyPlot(Form("%s/iter_v_sort_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::HGTD_SORT),
  // 		&mapHGTD.at(Score::ITERATIVE)
  // 	    });

  // // HGTD base times: TRKPTZ vs DNN
  // moneyPlot(Form("%s/trkptz_dnn_hgtd_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TESTML)
  // 	    });

  // // TRKPTZ full sample vs TRKPTZ restricted to highly pure clusters
  // moneyPlot(Form("%s/pure_clusters_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_MISCL)
  // 	    });

  // // Misassignment effect: HGTD vs TRKPTZ vs TEST_MISAS (misassigned tracks removed)
  // moneyPlot(Form("%s/test_misas_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_MISAS)
  // 	    });

  // // PU-contamination effect: HGTD vs TRKPTZ vs TEST_HS (HS-origin tracks only)
  // moneyPlot(Form("%s/test_hs_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_HS)
  // 	    });

  // // Full error decomposition: misclustering vs misassignment vs PU-contamination
  // moneyPlot(Form("%s/error_decomp_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_MISCL),
  // 		&mapIdealRes.at(Score::TRKPTZ),
  // 		&mapIdealRes.at(Score::PASS),
  // 		&mapIdealEff.at(Score::PASS)
  // 	    });

  // moneyPlot(Form("%s/idealres_error_decomp_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TEST_MISCL),
  //               &mapIdealRes.at(Score::TEST_MISAS)
  // 	    });

  // Check iterative clustering score
  // moneyPlot(Form("%s/iterative_clust_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::ITERATIVE)
  // 	    });

  // HGTD BDT on cone clusters vs HGTD_SORT (BDT on pT-sorted simultaneous) vs TRKPTZ
  // moneyPlot(Form("%s/cone_bdt_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::CONE_BDT)
  // 	    });

  // pure clusters with ideal resolution
  // moneyPlot(Form("%s/pure_clusters_ires_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TEST_MISCL)
  // 	    });

  // moneyPlot(Form("%s/pure_clusters_ieff_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealEff.at(Score::TRKPTZ),
  //               &mapIdealEff.at(Score::TEST_MISCL)
  // 	    });

  // // Ideal-resolution times: TRKPTZ vs DNN
  // moneyPlot(Form("%s/trkptz_dnn_ires_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TESTML)
  // 	    });

  // // Ideal-efficiency times: TRKPTZ vs DNN
  // moneyPlot(Form("%s/trkptz_dnn_ieff_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealEff.at(Score::TRKPTZ),
  //               &mapIdealEff.at(Score::TESTML)
  // 	    });

  // // Full ideal comparison: HGTD → TRKPTZ → IdealRes → IdealEff
  // moneyPlot(Form("%s/ideal_comp_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealEff.at(Score::TRKPTZ)
  // 	    });

  // Effect of fixing HGTD matching alone
  // moneyPlot(Form("%s/fixed_assoc_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(HGTD),
  //               &mapHGTD.at(TRKPTZ),
  //               &mapIdealRes.at(TRKPTZ),
  //               &mapIdealEff.at(TRKPTZ)
  // 	    });

  // Effect of fixing cluster selection alone
  // moneyPlot(Form("%s/fixed_selection_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(HGTD),
  //               &mapHGTD.at(TRKPTZ),
  //               &mapHGTD.at(PASS)
  // 	    });

  // Effect of fixing everything
  // moneyPlot(Form("%s/fixed_all_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(HGTD),
  //               &mapHGTD.at(TRKPTZ),
  //               &mapIdealEff.at(PASS)
  // 	    });

}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

auto main() -> int {
  SetAtlasStyle();  // ATLAS publication style (fonts, margins, ticks, …)

  // --- Data source ---
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  ROOT::EnableImplicitMT(); // use all CPU cores

  gErrorIgnoreLevel = kFatal;

  // --- Canvas ---
  TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);

  // --- Analysis maps (one per timing scenario) ---
  auto mapHGTD     = buildAnalysisMap(Scenario::HGTD      );
  auto mapIdealRes = buildAnalysisMap(Scenario::IDEAL_RES );
  auto mapIdealEff = buildAnalysisMap(Scenario::IDEAL_EFF );

  auto allMaps = { &mapHGTD, &mapIdealRes, &mapIdealEff };

  // --- Event display candidate lists ---
  std::vector<TString> evtDisplayHGTD, evtDisplayIdealRes, evtDisplayIdealEff;

  // --- Event loop ---
  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  while (reader.Next()) {
    const Long64_t READ_NUM  = chain.GetReadEntry() + 1;
    const Long64_t EVENT_NUM = chain.GetReadEntry() - chain.GetChainOffset();

    // if (READ_NUM > 1000) break;

    if (READ_NUM % 100 == 0)
      std::cout << "Progress: " << READ_NUM << "/" << N_EVENT << "\r" << std::flush;

    // Run the three timing scenarios
    auto resHGTD     = processEventData(&branch, false, true,  false, mapHGTD    );
    auto resIdealRes = processEventData(&branch, true,  true,  false, mapIdealRes);
    auto resIdealEff = processEventData(&branch, true,  false, false, mapIdealEff);

    // Extract file identifier from the full path (characters 49–54)
    TString fileName = branch.reader.GetTree()->GetCurrentFile()->GetName();
    TString fileNum  = fileName(49, 6);

    // Collect events where TRKPTZ passes but TEST_MISAS does not (misassignment effect)
    collectEventDisplay(evtDisplayHGTD,      3, resHGTD,     fileNum, EVENT_NUM);
    collectEventDisplay(evtDisplayIdealRes,  3, resIdealRes, fileNum, EVENT_NUM);
    collectEventDisplay(evtDisplayIdealEff,  3, resIdealEff, fileNum, EVENT_NUM);
  }

  for (auto* m : allMaps)
    for (auto& [k, analysis] : *m)
      analysis.postProcessing();

  // --- Per-analysis-object plots ---
  for (auto* m : allMaps)
    for (auto& [k, analysis] : *m)
      analysis.fullPlotting(canvas);

  for (auto* m : allMaps)
      for (auto& [k, analysis] : *m)
	analysis.printEfficiencyStats("hs_track");

  std::cout << "\nFINISHED PROCESSING\n";

  const auto KEYS = {"pu_frac", "fjet", "hs_track"};

  // --- Comparison plots (per variable KEY) ---
  for (const auto* key : KEYS)
    makeComparisonPlots(key, canvas, mapHGTD, mapIdealRes, mapIdealEff);

  // --- Inclusive resolution plots ---
  const std::initializer_list<AnalysisObj*> RESO_SET = {
    &mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::TRKPTZ), &mapHGTD.at(Score::TEST_MISCL), &mapIdealRes.at(Score::TRKPTZ)
    // &mapHGTD.at(Score::REFINED)
    // &mapHGTD.at(Score::TEST_MISAS), &mapHGTD.at(Score::TEST_HS),
    // &mapIdealRes.at(Score::TRKPTZ),
    // &mapIdealEff.at(Score::TRKPTZ), &mapIdealEff.at(Score::PASS),
  };
  inclusivePlot(Form("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET);
  inclusivePlot(Form("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET);

  // Low-track (nHSTrack <= 5) inclusive plots
  auto lowTrackGetter = [](AnalysisObj* a){ return a->inclusiveResoLowTrack.get(); };
  inclusivePlot(Form("%s/inclusive/inclusivereso_lowtrack_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET, lowTrackGetter);
  inclusivePlot(Form("%s/inclusive/inclusivereso_lowtrack_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET, lowTrackGetter);

  std::cout << "FINISHED PLOT PRINTING\n";

  // --- Print event display commands (toggle PRINT_EVENT_DISPLAYS above) ---
  printEventDisplays("HGTD times: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayHGTD    );
  printEventDisplays("Ideal Res.: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayIdealRes);
  printEventDisplays("Ideal Eff.: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayIdealEff);

  return 0;
}
