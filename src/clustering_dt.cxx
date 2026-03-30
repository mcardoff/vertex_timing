#include <random>
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

// Max event-display commands to print per low-multiplicity category.
static constexpr int N_LOW_MULT_DISPLAY = 20;
// HS-track count defining "low multiplicity" for the event-display sample.
static constexpr int LOW_MULT_NHS = 5;

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
  // m.emplace(Score::TEST_ML,    AnalysisObj(label, Score::TEST_ML   ));
  // m.emplace(Score::PASS,       AnalysisObj(label, Score::PASS      ));
  m.emplace(Score::TEST_MISCL, AnalysisObj(label, Score::TEST_MISCL));
  m.emplace(Score::Z_REFINED,  AnalysisObj(label, Score::Z_REFINED ));
  m.emplace(Score::ZT_REFINED, AnalysisObj(label, Score::ZT_REFINED));
  m.emplace(Score::T_REFINED,  AnalysisObj(label, Score::T_REFINED ));
  m.emplace(Score::ZT_ITER,    AnalysisObj(label, Score::ZT_ITER   ));
  m.emplace(Score::TEST_MISAS, AnalysisObj(label, Score::TEST_MISAS));
  m.emplace(Score::TEST_CTIME, AnalysisObj(label, Score::TEST_CTIME));
  m.emplace(Score::PERF_EVT,   AnalysisObj(label, Score::PERF_EVT  ));
  m.emplace(Score::PERF_CLT,   AnalysisObj(label, Score::PERF_CLT  ));

  // Scores active only in the real-HGTD scenario
  if (scenario == Scenario::HGTD) {
    m.emplace(Score::HGTD,       AnalysisObj(label, Score::HGTD      ));
    // m.emplace(Score::TEST_HS,    AnalysisObj(label, Score::TEST_HS    ));
    // m.emplace(Score::HGTD_SORT, AnalysisObj(label, Score::HGTD_SORT));
    // m.emplace(Score::CONE,       AnalysisObj(label, Score::CONE      ));
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
  const EventResult& result,
  const TString& fileNum,
  Long64_t eventNum
) {
  if (result.code == returnCode)
    list.push_back(TString::Format(EVTDISPLAY_FMT, fileNum.Data(), eventNum, result.time));
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

  moneyPlot(TString::Format("%s/theplot_%s.pdf",compSubdir.c_str(), key), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::TEST_MISCL),
	      &mapHGTD.at(Score::TEST_MISAS),
	    });
  
  // Timing oracle comparison: MISAS (event-level) vs CTIME (cluster-level)
  moneyPlot(TString::Format("%s/timing_oracle_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
                &mapHGTD.at(Score::TEST_MISCL),
                &mapHGTD.at(Score::TEST_MISAS),
                // &mapHGTD.at(Score::TEST_CTIME),
            });

  moneyPlot(TString::Format("%s/perfect_timing_oracle_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
                &mapHGTD.at(Score::TEST_MISCL),
                &mapHGTD.at(Score::TEST_MISAS),
                // &mapHGTD.at(Score::TEST_CTIME),
                &mapHGTD.at(Score::PERF_EVT),
                // &mapHGTD.at(Score::PERF_CLT),
            });

  // HGTD algo vs TRKPTZ
  moneyPlot(TString::Format("%s/hgtd_trkptz_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ)
	    });

  // HGTD vs TRKPTZ vs T_REFINED vs Z_REFINED
  moneyPlot(TString::Format("%s/z_refined_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
                &mapHGTD.at(Score::T_REFINED),
                &mapHGTD.at(Score::Z_REFINED),
                &mapHGTD.at(Score::ZT_REFINED),
                &mapHGTD.at(Score::ZT_ITER),
	    });

  moneyPlot(TString::Format("%s/2d_clust_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
		&mapHGTD.at(Score::ZT_ITER),
	    });


  // moneyPlot(TString::Format("%s/z_refined_ideal_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               // &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               // &mapHGTD.at(Score::T_REFINED),
  //               // &mapHGTD.at(Score::Z_REFINED),
  //               // &mapHGTD.at(Score::ZT_REFINED),
  //               &mapIdealRes.at(Score::T_REFINED),
  //               &mapIdealRes.at(Score::Z_REFINED)
  //           });

  // // HGTD algo vs TRKPTZ vs HGTD_SORT (pT-sorted simultaneous + BDT)
  // moneyPlot(TString::Format("%s/hgtd_sort_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::HGTD_SORT)
  // 	    });

  // HGTD algo vs TRKPTZ vs HGTD_SORT (pT-sorted simultaneous + BDT)
  // moneyPlot(TString::Format("%s/iter_v_sort_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::HGTD_SORT),
  // 		&mapHGTD.at(Score::CONE)
  // 	    });

  // // HGTD base times: TRKPTZ vs DNN
  // moneyPlot(TString::Format("%s/trkptz_dnn_hgtd_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_ML)
  // 	    });

  // // TRKPTZ full sample vs TRKPTZ restricted to highly pure clusters
  // moneyPlot(TString::Format("%s/pure_clusters_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_MISCL)
  // 	    });

  // // Misassignment effect: HGTD vs TRKPTZ vs TEST_MISAS (misassigned tracks removed)
  // moneyPlot(TString::Format("%s/test_misas_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_MISAS)
  // 	    });

  // // PU-contamination effect: HGTD vs TRKPTZ vs TEST_HS (HS-origin tracks only)
  // moneyPlot(TString::Format("%s/test_hs_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_HS)
  // 	    });

  // // Full error decomposition: misclustering vs misassignment vs PU-contamination
  // moneyPlot(TString::Format("%s/error_decomp_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::TEST_MISCL),
  // 		&mapIdealRes.at(Score::TRKPTZ),
  // 		&mapIdealRes.at(Score::PASS),
  // 		&mapIdealEff.at(Score::PASS)
  // 	    });

  // moneyPlot(TString::Format("%s/idealres_error_decomp_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TEST_MISCL),
  //               &mapIdealRes.at(Score::TEST_MISAS)
  // 	    });

  // Check iterative clustering score
  // moneyPlot(TString::Format("%s/iterative_clust_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::CONE)
  // 	    });

  // HGTD BDT on cone clusters vs HGTD_SORT (BDT on pT-sorted simultaneous) vs TRKPTZ
  // moneyPlot(TString::Format("%s/cone_bdt_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapHGTD.at(Score::CONE_BDT)
  // 	    });

  // pure clusters with ideal resolution
  // moneyPlot(TString::Format("%s/pure_clusters_ires_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TEST_MISCL)
  // 	    });

  // moneyPlot(TString::Format("%s/pure_clusters_ieff_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealEff.at(Score::TRKPTZ),
  //               &mapIdealEff.at(Score::TEST_MISCL)
  // 	    });

  // // Ideal-resolution times: TRKPTZ vs DNN
  // moneyPlot(TString::Format("%s/trkptz_dnn_ires_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TEST_ML)
  // 	    });

  // // Ideal-efficiency times: TRKPTZ vs DNN
  // moneyPlot(TString::Format("%s/trkptz_dnn_ieff_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapIdealEff.at(Score::TRKPTZ),
  //               &mapIdealEff.at(Score::TEST_ML)
  // 	    });

  // // Full ideal comparison: HGTD → TRKPTZ → IdealRes → IdealEff
  // moneyPlot(TString::Format("%s/ideal_comp_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(Score::HGTD),
  //               &mapHGTD.at(Score::TRKPTZ),
  //               &mapIdealRes.at(Score::TRKPTZ),
  //               &mapIdealEff.at(Score::TRKPTZ)
  // 	    });

  // Effect of fixing HGTD matching alone
  // moneyPlot(TString::Format("%s/fixed_assoc_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(HGTD),
  //               &mapHGTD.at(TRKPTZ),
  //               &mapIdealRes.at(TRKPTZ),
  //               &mapIdealEff.at(TRKPTZ)
  // 	    });

  // Effect of fixing cluster selection alone
  // moneyPlot(TString::Format("%s/fixed_selection_%s.pdf", compSubdir.c_str(), key), key, canvas,
  //           {
  //               &mapHGTD.at(HGTD),
  //               &mapHGTD.at(TRKPTZ),
  //               &mapHGTD.at(PASS)
  // 	    });

  // Effect of fixing everything
  // moneyPlot(TString::Format("%s/fixed_all_%s.pdf", compSubdir.c_str(), key), key, canvas,
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

  auto allMaps = { &mapHGTD}; //, &mapIdealRes, &mapIdealEff };

  // --- Event display candidate lists ---
  std::vector<TString> evtDisplayHGTD, evtDisplayIdealRes, evtDisplayIdealEff;
  // Low-multiplicity (nFwdHS == LOW_MULT_NHS) pass and fail, HGTD scenario
  std::vector<TString> lowMultPass, lowMultFail;

  // --- Event loop ---
  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  while (reader.Next()) {
    const Long64_t READ_NUM  = chain.GetReadEntry() + 1;
    const Long64_t EVENT_NUM = chain.GetReadEntry() - chain.GetChainOffset();

    // if (READ_NUM > 10000) break;

    if (READ_NUM % 100 == 0)
      std::cout << "Progress: " << READ_NUM << "/" << N_EVENT << "\r" << std::flush;

    // Run the three timing scenarios
    auto resHGTD     = processEventData(&branch, false, true,  mapHGTD    );
    // auto resIdealRes = processEventData(&branch, true,  true,  mapIdealRes);
    // auto resIdealEff = processEventData(&branch, true,  false, mapIdealEff);

    // Extract file identifier from the full path (characters 49–54)
    TString fileName = branch.reader.GetTree()->GetCurrentFile()->GetName();
    TString fileNum  = fileName(49, 6);

    // Collect events where TRKPTZ passes but TEST_MISAS does not (misassignment effect)
    collectEventDisplay(evtDisplayHGTD,      3, resHGTD,     fileNum, EVENT_NUM);
    // collectEventDisplay(evtDisplayIdealRes,  3, resIdealRes, fileNum, EVENT_NUM);
    // collectEventDisplay(evtDisplayIdealEff,  3, resIdealEff, fileNum, EVENT_NUM);

    // Low-multiplicity event display collection (HGTD scenario, n=LOW_MULT_NHS HS tracks)
    if (resHGTD.code >= 0 && resHGTD.nFwdHS == LOW_MULT_NHS) {
      TString cmd = TString::Format(EVTDISPLAY_FMT, fileNum.Data(), EVENT_NUM, resHGTD.time);
      (resHGTD.trkptzPass ? lowMultPass : lowMultFail).push_back(cmd);
    }
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
	// analysis.printResolutionStats("hs_track");


  std::cout << "\nFINISHED PROCESSING\n";

  const auto KEYS = {"pu_frac", "fjet", "hs_track"};

  // --- Comparison plots (per variable KEY) ---
  for (const auto* key : KEYS)
    makeComparisonPlots(key, canvas, mapHGTD, mapIdealRes, mapIdealEff);

  // --- Inclusive resolution plots ---
  const std::initializer_list<AnalysisObj*> RESO_SET = {
    &mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::TRKPTZ), &mapHGTD.at(Score::TEST_MISCL),
    &mapHGTD.at(Score::TEST_MISAS), &mapHGTD.at(Score::TEST_CTIME),
  };
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET);

  // Low-track (nHSTrack <= 5) inclusive plots
  auto lowTrackGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoLowTrackSig.get(),
             a->inclusiveResoLowTrackMix.get(),
             a->inclusiveResoLowTrackBkg.get() }; };
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_lowtrack_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET, lowTrackGetter);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_lowtrack_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET, lowTrackGetter);

  std::cout << "FINISHED PLOT PRINTING\n";

  // --- Print event display commands (toggle PRINT_EVENT_DISPLAYS above) ---
  printEventDisplays("HGTD times: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayHGTD    );
  // printEventDisplays("Ideal Res.: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayIdealRes);
  // printEventDisplays("Ideal Eff.: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayIdealEff);

  // --- Low-multiplicity event displays (always printed, random subsample) ---
  if (PRINT_EVENT_DISPLAYS) {
    std::mt19937 rng(std::random_device{}());  // non-deterministic seed — different sample each run
    auto subsample = [&](std::vector<TString>& v) {
      std::shuffle(v.begin(), v.end(), rng);
      if ((int)v.size() > N_LOW_MULT_DISPLAY)
        v.resize(N_LOW_MULT_DISPLAY);
    };
    subsample(lowMultPass);
    subsample(lowMultFail);

    auto printGroup = [](const char* header, const std::vector<TString>& cmds) {
      std::cout << "\n=== " << header << " (" << cmds.size() << " events) ===\n";
      for (const auto& c : cmds) std::cout << c << '\n';
    };
    printGroup(
      TString::Format("Low-mult (nFwdHS==%d, HGTD) — TRKPTZ PASS", LOW_MULT_NHS).Data(),
      lowMultPass);
    printGroup(
      TString::Format("Low-mult (nFwdHS==%d, HGTD) — TRKPTZ FAIL", LOW_MULT_NHS).Data(),
      lowMultFail);
  }

  return 0;
}
