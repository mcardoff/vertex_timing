#include <RtypesCore.h>
#include <TROOT.h>
#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "suppress_stdout.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

// Output directory for saved plots
static const char* SAVE_DIR = "../figs-2sigma";

// Event display command template (filled with file number, event number, extra time)
static const char* EVTDISPLAY_FMT =
  "python3 event_display.py --file_num %s --event_num %lld --extra_time %.2f";

// Set to true to print event display commands to stdout after the event loop.
static constexpr bool PRINT_EVENT_DISPLAYS = true;

// ---------------------------------------------------------------------------
// Helper: build one analysis map for a given scenario label.
//   The active scores are listed explicitly so adding/removing a score
//   only requires changing this one function.
// ---------------------------------------------------------------------------

enum class Scenario { HGTD, IdealRes, IdealEff };

auto buildAnalysisMap(Scenario scenario) -> std::map<Score, AnalysisObj> {
  const char* label = [&]() -> const char* {
    switch (scenario) {
      case Scenario::HGTD:     return "HGTD Times";
      case Scenario::IdealRes: return "Ideal Res. HGTD";
      case Scenario::IdealEff: return "Ideal Res.+Eff. HGTD";
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

// Print all collected event display commands for one scenario.
// Controlled by PRINT_EVENT_DISPLAYS at the top of this file.
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
  const char* KEY,
  TCanvas* canvas,
  std::map<Score, AnalysisObj>& mapHGTD,
  std::map<Score, AnalysisObj>& mapIdealRes,
  std::map<Score, AnalysisObj>& mapIdealEff
) {
  // HGTD algo vs TRKPTZ
  moneyPlot(Form("%s/hgtd_trkptz_%s.pdf",            SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ) });

  // HGTD base times: TRKPTZ vs DNN
  moneyPlot(Form("%s/trkptz_dnn_basetimes_%s.pdf",   SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapHGTD.at(TESTML) });

  // Misclustering study: TRKPTZ full sample vs TRKPTZ restricted to pure-cluster events
  moneyPlot(Form("%s/misclustering_study_%s.pdf",    SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapHGTD.at(TEST_MISCL) });

  // Ideal-resolution times: TRKPTZ vs DNN
  moneyPlot(Form("%s/trkptz_dnn_idealrestimes_%s.pdf", SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapIdealRes.at(TRKPTZ), &mapIdealRes.at(TESTML) });

  // Ideal-efficiency times: TRKPTZ vs DNN
  moneyPlot(Form("%s/trkptz_dnn_idealefftimes_%s.pdf", SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapIdealEff.at(TRKPTZ), &mapIdealEff.at(TESTML) });

  // Effect of fixing HGTD matching alone
  moneyPlot(Form("%s/fixed_assoc_%s.pdf",            SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ),
              &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ) });

  // Effect of fixing cluster selection alone
  moneyPlot(Form("%s/fixed_selection_%s.pdf",        SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapHGTD.at(PASS) });

  // Effect of fixing everything
  moneyPlot(Form("%s/fixed_all_%s.pdf",              SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ), &mapIdealEff.at(PASS) });

  // Full ideal comparison: HGTD → TRKPTZ → IdealRes → IdealEff
  moneyPlot(Form("%s/ideal_comp_%s.pdf",             SAVE_DIR, KEY), KEY, canvas,
            { &mapHGTD.at(HGTD), &mapHGTD.at(TRKPTZ),
              &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ) });
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

auto main() -> int {
  gROOT->SetBatch(true);  // batch mode: no display, and suppresses OBJ teardown messages
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
  auto mapHGTD     = buildAnalysisMap(Scenario::HGTD    );
  auto mapIdealRes = buildAnalysisMap(Scenario::IdealRes );
  auto mapIdealEff = buildAnalysisMap(Scenario::IdealEff );

  auto allMaps = { &mapHGTD, &mapIdealRes, &mapIdealEff };

  // --- Event display candidate lists ---
  std::vector<TString> evtDisplayHGTD, evtDisplayIdealRes, evtDisplayIdealEff;

  // --- Event loop ---
  std::cout << "Starting Event Loop\n";
  const Long64_t nEvent = chain.GetEntries();

  while (reader.Next()) {
    const Long64_t readNum  = chain.GetReadEntry() + 1;
    const Long64_t eventNum = chain.GetReadEntry() - chain.GetChainOffset();

    if (readNum % 100 == 0)
      std::cout << "Progress: " << readNum << "/" << nEvent << "\r" << std::flush;

    // Run the three timing scenarios
    auto resHGTD     = processEventData(&branch, false, true,  false, mapHGTD    );
    auto resIdealRes = processEventData(&branch, true,  true,  false, mapIdealRes);
    auto resIdealEff = processEventData(&branch, true,  false, false, mapIdealEff);

    // Extract file identifier from the full path (characters 49–54)
    TString fileName = branch.reader.GetTree()->GetCurrentFile()->GetName();
    TString fileNum  = fileName(49, 6);

    // Collect events where TEST_MISCL fails the efficiency test
    collectEventDisplay(evtDisplayHGTD,      2, resHGTD,     fileNum, eventNum);
    collectEventDisplay(evtDisplayIdealRes,  2, resIdealRes, fileNum, eventNum);
    collectEventDisplay(evtDisplayIdealEff,  2, resIdealEff, fileNum, eventNum);
  }

  // --- Post-processing and plotting (suppress ROOT's OBJ: printf chatter) ---
  {
    SuppressStdout suppress;
    for (auto* m : allMaps)
      for (auto& [k, analysis] : *m)
        analysis.postProcessing();

    // --- Per-analysis-object plots ---
    for (auto* m : allMaps)
      for (auto& [k, analysis] : *m)
        analysis.fullPlotting(canvas);
  }

  // Shut down the thread pool before any ROOT object cleanup runs.
  // Without this, worker threads still hold references to TEfficiency objects
  // when the main-thread destructor fires, causing spurious "OBJ: ..._clone"
  // messages at kInfo level during teardown.
  ROOT::DisableImplicitMT();

  std::cout << "\nFINISHED PROCESSING\n";

  // --- Comparison plots (per variable KEY) ---
  for (const auto* KEY : {"pu_frac", "ftrack", "pu_track", "fjet", "hs_track", "vtx_dz"})
    makeComparisonPlots(KEY, canvas, mapHGTD, mapIdealRes, mapIdealEff);

  // --- Inclusive resolution plots ---
  const std::initializer_list<AnalysisObj*> resoSet = {
    &mapHGTD.at(HGTD),   &mapHGTD.at(TRKPTZ),
    &mapIdealRes.at(TRKPTZ), &mapIdealEff.at(TRKPTZ), &mapIdealEff.at(PASS),
  };
  inclusivePlot(Form("%s/inclusivereso_logscale.pdf", SAVE_DIR),
                true,  false, -400, 400, canvas, resoSet);
  inclusivePlot(Form("%s/inclusivereso_linscale.pdf", SAVE_DIR),
                false, false, -200, 200, canvas, resoSet);

  std::cout << "FINISHED PLOT PRINTING\n";

  // --- Print event display commands (toggle PRINT_EVENT_DISPLAYS above) ---
  printEventDisplays("HGTD times: TEST_MISCL fails efficiency", evtDisplayHGTD    );
  printEventDisplays("Ideal Res.: TEST_MISCL fails efficiency", evtDisplayIdealRes);
  printEventDisplays("Ideal Eff.: TEST_MISCL fails efficiency", evtDisplayIdealEff);

  return 0;
}
