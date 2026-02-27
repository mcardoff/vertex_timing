#include <RtypesCore.h>
#include <TROOT.h>
#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"

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
  m.emplace(TRKPTZ, AnalysisObj(label, TRKPTZ));
  m.emplace(TESTML, AnalysisObj(label, TESTML));
  m.emplace(PASS,   AnalysisObj(label, PASS  ));
  m.emplace(TEST_MISCL, AnalysisObj(label, TEST_MISCL));

  // Scores active only in the real-HGTD scenario
  if (scenario == Scenario::HGTD) {
    m.emplace(HGTD,       AnalysisObj(label, HGTD      ));
    // Uncomment to re-enable additional algorithms:
    // m.emplace(CALO60,  AnalysisObj(label, CALO60 ));
    // m.emplace(CALO90,  AnalysisObj(label, CALO90 ));
    // m.emplace(JUST60,  AnalysisObj(label, JUST60 ));
    // m.emplace(JUST90,  AnalysisObj(label, JUST90 ));
    // m.emplace(FILT60,  AnalysisObj(label, FILT60 ));
    // m.emplace(FILT90,  AnalysisObj(label, FILT90 ));
    // m.emplace(FILTJET, AnalysisObj(label, FILTJET));
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
  // HGTD algo vs TRKPTZ
  moneyPlot(Form("%s/hgtd_trkptz_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapHGTD.at(TRKPTZ)
	    });

  // HGTD base times: TRKPTZ vs DNN
  moneyPlot(Form("%s/trkptz_dnn_hgtd_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapHGTD.at(TRKPTZ),
                &mapHGTD.at(TESTML)
	    });

  // TRKPTZ full sample vs TRKPTZ restricted to highly pure clusters
  moneyPlot(Form("%s/pure_clusters_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapHGTD.at(TRKPTZ),
                &mapHGTD.at(TEST_MISCL)
	    });

  // pure clusters with ideal resolution
  moneyPlot(Form("%s/pure_clusters_ires_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapIdealRes.at(TRKPTZ),
                &mapIdealRes.at(TEST_MISCL)
	    });

  moneyPlot(Form("%s/pure_clusters_ieff_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapIdealEff.at(TRKPTZ),
                &mapIdealEff.at(TEST_MISCL)
	    });

  // Ideal-resolution times: TRKPTZ vs DNN
  moneyPlot(Form("%s/trkptz_dnn_ires_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapIdealRes.at(TRKPTZ),
                &mapIdealRes.at(TESTML)
	    });

  // Ideal-efficiency times: TRKPTZ vs DNN
  moneyPlot(Form("%s/trkptz_dnn_ieff_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapIdealEff.at(TRKPTZ),
                &mapIdealEff.at(TESTML)
	    });

  // Full ideal comparison: HGTD → TRKPTZ → IdealRes → IdealEff
  moneyPlot(Form("%s/ideal_comp_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(HGTD),
                &mapHGTD.at(TRKPTZ),
                &mapIdealRes.at(TRKPTZ),
                &mapIdealEff.at(TRKPTZ)
	    });

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
  gStyle->SetOptStat(0);

  // --- Data source ---
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  ROOT::EnableImplicitMT(); // use all CPU cores

  gErrorIgnoreLevel = kFatal;

  // --- Canvas ---
  TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

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

    // Collect events where TEST_MISCL fails the efficiency test
    collectEventDisplay(evtDisplayHGTD,      2, resHGTD,     fileNum, EVENT_NUM);
    collectEventDisplay(evtDisplayIdealRes,  2, resIdealRes, fileNum, EVENT_NUM);
    collectEventDisplay(evtDisplayIdealEff,  2, resIdealEff, fileNum, EVENT_NUM);
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
    &mapHGTD.at(HGTD),   &mapHGTD.at(TRKPTZ),
    &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ), &mapIdealEff.at(PASS),
  };
  inclusivePlot(Form("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET);
  inclusivePlot(Form("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET);

  std::cout << "FINISHED PLOT PRINTING\n";

  // --- Print event display commands (toggle PRINT_EVENT_DISPLAYS above) ---
  printEventDisplays("HGTD times: TEST_MISCL fails efficiency", evtDisplayHGTD    );
  printEventDisplays("Ideal Res.: TEST_MISCL fails efficiency", evtDisplayIdealRes);
  printEventDisplays("Ideal Eff.: TEST_MISCL fails efficiency", evtDisplayIdealEff);

  return 0;
}
