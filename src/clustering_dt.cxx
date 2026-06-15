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

// Max commands for each new diagnostic category.
// (All three are gated on PRINT_EVENT_DISPLAYS; comment out individual blocks below
//  to disable a specific category without disabling the others.)
static constexpr int N_T_REFINED_UNCHANGED = 10;  // TRKPTZ↔T_REFINED status agrees
static constexpr int N_MISCL_PASS          =  5;  // pure cluster, PASSES timing window
static constexpr int N_MISAS_PASS          =  5;  // clean HS timing, PASSES timing window

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
  m.emplace(Score::WAVES,   AnalysisObj(label, Score::WAVES  ));
  m.emplace(Score::JET_T_REFINED, AnalysisObj(label, Score::JET_T_REFINED));
  m.emplace(Score::TEST_ML,    AnalysisObj(label, Score::TEST_ML   ));
  // m.emplace(Score::PASS,       AnalysisObj(label, Score::PASS      ));
  m.emplace(Score::TEST_MISCL, AnalysisObj(label, Score::TEST_MISCL));
  m.emplace(Score::Z_REFINED,  AnalysisObj(label, Score::Z_REFINED ));
  m.emplace(Score::ZT_REFINED, AnalysisObj(label, Score::ZT_REFINED));
  m.emplace(Score::T_REFINED,  AnalysisObj(label, Score::T_REFINED ));
  m.emplace(Score::ZT_ITER,    AnalysisObj(label, Score::ZT_ITER   ));
  m.emplace(Score::TEST_MISAS, AnalysisObj(label, Score::TEST_MISAS));
  m.emplace(Score::PERF_EVT,   AnalysisObj(label, Score::PERF_EVT  ));
  m.emplace(Score::WAVES_MISCL, AnalysisObj(label, Score::WAVES_MISCL));
  m.emplace(Score::WAVES_MISAS, AnalysisObj(label, Score::WAVES_MISAS));

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

  // WAVeS oracle comparison: selection by WAVeS score, gated like the TRKPTZ oracles.
  // Color override: WAVES cyan, perfect timing violet, pure clusters yellow.
  moneyPlot(TString::Format("%s/waves_oracle_%s.pdf",compSubdir.c_str(), key), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::WAVES_MISAS),
	      &mapHGTD.at(Score::WAVES_MISCL),
	    },
	    {C01, C02, C05, C04, C03});
  
  moneyPlot(TString::Format("%s/perfect_timing_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
                &mapHGTD.at(Score::TEST_MISCL),
                &mapHGTD.at(Score::TEST_MISAS),
                &mapHGTD.at(Score::PERF_EVT),
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

  // WAVeS: HGTD | TRKPTZ | WAVES (WAVeS score) | JET_T_REFINED (jet-filtered 2σ clustering)
  moneyPlot(TString::Format("%s/jet_refined_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
                &mapHGTD.at(Score::WAVES),
                &mapHGTD.at(Score::JET_T_REFINED),
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

  // HGTD base times: TRKPTZ vs DNN
  moneyPlot(TString::Format("%s/trkptz_dnn_hgtd_%s.pdf", compSubdir.c_str(), key), key, canvas,
            {
                &mapHGTD.at(Score::HGTD),
                &mapHGTD.at(Score::TRKPTZ),
                &mapHGTD.at(Score::TEST_ML)
	    });

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
  // New diagnostic categories — pass and fail variants
  std::vector<TString> tRefinedUnchanged;  // TRKPTZ & T_REFINED have the same pass/fail
  std::vector<TString> tRefinedChanged;    // TRKPTZ & T_REFINED disagree (T_REFINED flipped outcome)
  std::vector<TString> misclPassEvents;    // in TEST_MISCL denominator and PASSES
  std::vector<TString> misclFailEvents;    // in TEST_MISCL denominator and FAILS  (returnCode==2)
  std::vector<TString> misasPassEvents;    // in TEST_MISAS denominator and PASSES
  std::vector<TString> misasFailEvents;    // in TEST_MISAS denominator and FAILS   (returnCode==3)
  // Cluster quality tiers for event display browsing (TRKPTZ-selected cluster, HGTD scenario)
  std::vector<TString> cqHighEvents;  // quality >= CLUS_QUALITY_SPLIT
  std::vector<TString> cqLowEvents;   // quality < CLUS_QUALITY_SPLIT

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

    if (resHGTD.code >= 0) {
      TString cmd = TString::Format(EVTDISPLAY_FMT, fileNum.Data(), EVENT_NUM, resHGTD.time);

      // Category 1: TRKPTZ and T_REFINED agree (both pass or both fail)
      if (resHGTD.trkptzPass == resHGTD.tRefinedPass)
        tRefinedUnchanged.push_back(cmd);
      // Category 1 fail: T_REFINED flipped the outcome relative to TRKPTZ
      else
        tRefinedChanged.push_back(cmd);

      // Category 2: event is in TEST_MISCL denominator (pure cluster) and PASSES
      if (resHGTD.misclPass)
        misclPassEvents.push_back(cmd);
      // Category 2 fail: event did NOT enter the TEST_MISCL denominator (cluster purity ≤ 75%)
      if (!resHGTD.misclInDenom)
        misclFailEvents.push_back(cmd);

      // Category 3: event is in TEST_MISAS denominator (clean HS timing) and PASSES
      if (resHGTD.misasPass)
        misasPassEvents.push_back(cmd);
      // Category 3 fail: event did NOT enter the TEST_MISAS denominator (hsTimingPurity < 95%)
      if (!resHGTD.misasInDenom)
        misasFailEvents.push_back(cmd);

      // Cluster quality tiers (TRKPTZ-selected cluster)
      if (resHGTD.clusQuality >= CLUS_QUALITY_SPLIT) cqHighEvents.push_back(cmd);
      else                                            cqLowEvents .push_back(cmd);
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

  const auto KEYS = {"pu_frac", "fjet", "ftrack", "hs_track", "nhit", "clus_pu_frac", "clus_sigma_t", "clus_quality"};

  // --- Comparison plots (per variable KEY) ---
  for (const auto* key : KEYS)
    makeComparisonPlots(key, canvas, mapHGTD, mapIdealRes, mapIdealEff);

  // --- Inclusive resolution plots ---
  const std::initializer_list<AnalysisObj*> RESO_SET = {
    &mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::TRKPTZ), &mapHGTD.at(Score::TEST_MISCL),
    &mapHGTD.at(Score::TEST_MISAS), &mapHGTD.at(Score::TEST_ML),
  };
  const char* diffLabel = "#Delta t [ps]";
  auto resoGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoSig.get(),
             a->inclusiveResoMix.get(),
             a->inclusiveResoBkg.get() }; };
  // Dynamic y-axis caps: scan every reso/pull histogram across every score in RESO_SET,
  // find the max stack bin content (sig+mix+bkg), then set the cap to 1.1× that for linear
  // and 5× for log (extra log headroom so the ATLAS label doesn't clip the peak).
  // This guarantees all inclusive plots share an identical y-scale so tail reduction
  // between scores or quality tiers is directly visible.
  auto stackMax = [](TH1D* sig, TH1D* mix, TH1D* bkg) -> double {
    double m = 0.0;
    int nb = sig->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
      double s = sig->GetBinContent(i) + mix->GetBinContent(i) + bkg->GetBinContent(i);
      if (s > m) m = s;
    }
    return m;
  };
  double globalResoMax = 0.0, globalPullMax = 0.0;
  for (auto* a : RESO_SET) {
    globalResoMax = std::max({globalResoMax,
      stackMax(a->inclusiveResoSig      .get(), a->inclusiveResoMix      .get(), a->inclusiveResoBkg      .get()),
      stackMax(a->inclusiveResoNhit1Sig .get(), a->inclusiveResoNhit1Mix .get(), a->inclusiveResoNhit1Bkg .get()),
      stackMax(a->inclusiveResoNhit2Sig .get(), a->inclusiveResoNhit2Mix .get(), a->inclusiveResoNhit2Bkg .get()),
      stackMax(a->inclusiveResoNhit3pSig.get(), a->inclusiveResoNhit3pMix.get(), a->inclusiveResoNhit3pBkg.get()),
      stackMax(a->inclusiveResoClusQHighSig.get(), a->inclusiveResoClusQHighMix.get(), a->inclusiveResoClusQHighBkg.get()),
      stackMax(a->inclusiveResoClusQLowSig .get(), a->inclusiveResoClusQLowMix .get(), a->inclusiveResoClusQLowBkg .get())});
    globalPullMax = std::max({globalPullMax,
      stackMax(a->inclusivePullSig         .get(), a->inclusivePullMix         .get(), a->inclusivePullBkg         .get()),
      stackMax(a->inclusivePullClusQHighSig.get(), a->inclusivePullClusQHighMix.get(), a->inclusivePullClusQHighBkg.get()),
      stackMax(a->inclusivePullClusQLowSig .get(), a->inclusivePullClusQLowMix .get(), a->inclusivePullClusQLowBkg .get())});
  }
  const double INCL_RESO_YMAX_LOG = 5.0 * globalResoMax;
  const double INCL_RESO_YMAX_LIN = 1.1 * globalResoMax;
  const double INCL_PULL_YMAX_LOG = 5.0 * globalPullMax;
  const double INCL_PULL_YMAX_LIN = 1.1 * globalPullMax;
  const double INCL_YMIN_LOG      = 0.5;

  // inclusivePlot(TString::Format("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR),
		// true,  false, -400, 400, canvas, RESO_SET);
  // inclusivePlot(TString::Format("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR),
		// false, false, -200, 200, canvas, RESO_SET);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET, resoGetter, nullptr, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET, resoGetter, nullptr, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);

  // nHGTD hits-binned inclusive resolution
  auto nhit1Getter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoNhit1Sig.get(),
             a->inclusiveResoNhit1Mix.get(),
             a->inclusiveResoNhit1Bkg.get() }; };
  auto nhit2Getter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoNhit2Sig.get(),
             a->inclusiveResoNhit2Mix.get(),
             a->inclusiveResoNhit2Bkg.get() }; };
  auto nhit3pGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoNhit3pSig.get(),
             a->inclusiveResoNhit3pMix.get(),
             a->inclusiveResoNhit3pBkg.get() }; };
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_nhit1_logscale.pdf", SAVE_DIR),
      true, false, -400, 400, canvas, RESO_SET, nhit1Getter, nullptr, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_nhit1_linscale.pdf", SAVE_DIR),
      false, false, -200, 200, canvas, RESO_SET, nhit1Getter, nullptr, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_nhit2_logscale.pdf", SAVE_DIR),
      true, false, -400, 400, canvas, RESO_SET, nhit2Getter, nullptr, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_nhit2_linscale.pdf", SAVE_DIR),
      false, false, -200, 200, canvas, RESO_SET, nhit2Getter, nullptr, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_nhit3p_logscale.pdf", SAVE_DIR),
      true, false, -400, 400, canvas, RESO_SET, nhit3pGetter, nullptr, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_nhit3p_linscale.pdf", SAVE_DIR),
      false, false, -200, 200, canvas, RESO_SET, nhit3pGetter, nullptr, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);

  // Cluster quality-binned inclusive resolution (two tiers: Q≥0.5 high, Q<0.5 low)
  auto cqHighGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoClusQHighSig.get(),
             a->inclusiveResoClusQHighMix.get(),
             a->inclusiveResoClusQHighBkg.get() }; };
  auto cqLowGetter  = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusiveResoClusQLowSig.get(),
             a->inclusiveResoClusQLowMix.get(),
             a->inclusiveResoClusQLowBkg.get() }; };
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_cqhigh_logscale.pdf", SAVE_DIR),
      true, false, -400, 400, canvas, RESO_SET, cqHighGetter, nullptr, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_cqhigh_linscale.pdf", SAVE_DIR),
      false, false, -200, 200, canvas, RESO_SET, cqHighGetter, nullptr, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_cqlow_logscale.pdf", SAVE_DIR),
      true, false, -400, 400, canvas, RESO_SET, cqLowGetter,  nullptr, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_cqlow_linscale.pdf", SAVE_DIR),
      false, false, -200, 200, canvas, RESO_SET, cqLowGetter,  nullptr, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);

  // Low-track (nHSTrack <= 5) inclusive plots
  // auto lowTrackGetter = [](AnalysisObj* a) -> ResoTriple {
  //   return { a->inclusiveResoLowTrackSig.get(),
  //            a->inclusiveResoLowTrackMix.get(),
  //            a->inclusiveResoLowTrackBkg.get() }; };
  // inclusivePlot(TString::Format("%s/inclusive/inclusivereso_lowtrack_logscale.pdf", SAVE_DIR),
  // 		true,  false, -400, 400, canvas, RESO_SET, lowTrackGetter);
  // inclusivePlot(TString::Format("%s/inclusive/inclusivereso_lowtrack_linscale.pdf", SAVE_DIR),
  // 		false, false, -200, 200, canvas, RESO_SET, lowTrackGetter);

  // --- Inclusive pull plots: Δt / σ_cluster ---
  // σ_pull ≈ 1 confirms the cluster timing uncertainty is correctly calibrated;
  // σ_pull > 1 means the uncertainty is underestimated.
  auto pullGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusivePullSig.get(),
             a->inclusivePullMix.get(),
             a->inclusivePullBkg.get() }; };
  // auto pullLowTrackGetter = [](AnalysisObj* a) -> ResoTriple {
  //   return { a->inclusivePullLowTrackSig.get(),
  //            a->inclusivePullLowTrackMix.get(),
  //            a->inclusivePullLowTrackBkg.get() }; };
  const char* pullLabel = "(t_{0}-t_{truth})/#sigma_{t}";
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_logscale.pdf", SAVE_DIR),
		true,  false, -10, 10, canvas, RESO_SET,
		pullGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_linscale.pdf", SAVE_DIR),
		false, false, -10, 10, canvas, RESO_SET,
		pullGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LIN, 0.0);

  // --- Per-quality-tier pull plots (two tiers) ---
  auto pullCqHighGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusivePullClusQHighSig.get(),
             a->inclusivePullClusQHighMix.get(),
             a->inclusivePullClusQHighBkg.get() }; };
  auto pullCqLowGetter  = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusivePullClusQLowSig.get(),
             a->inclusivePullClusQLowMix.get(),
             a->inclusivePullClusQLowBkg.get() }; };
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_cqhigh_logscale.pdf", SAVE_DIR),
		true,  false, -10, 10, canvas, RESO_SET, pullCqHighGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_cqhigh_linscale.pdf", SAVE_DIR),
		false, false, -10, 10, canvas, RESO_SET, pullCqHighGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LIN, 0.0);
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_cqlow_logscale.pdf", SAVE_DIR),
		true,  false, -10, 10, canvas, RESO_SET, pullCqLowGetter,  &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_cqlow_linscale.pdf", SAVE_DIR),
		false, false, -10, 10, canvas, RESO_SET, pullCqLowGetter,  &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LIN, 0.0);
  // inclusivePlot(TString::Format("%s/inclusive/inclusivepull_lowtrack_logscale.pdf", SAVE_DIR),
		// true,  false, -10, 10, canvas, RESO_SET,
		// pullLowTrackGetter, &FIT_PULLGAUS, pullLabel);
  // inclusivePlot(TString::Format("%s/inclusive/inclusivepull_lowtrack_linscale.pdf", SAVE_DIR),
		// false, false, -10, 10, canvas, RESO_SET,
		// pullLowTrackGetter, &FIT_PULLGAUS, pullLabel);

  // --- Normalized shape comparison: TRKPTZ vs oracle scores ---
  const std::vector<AnalysisObj*> SHAPE_SET = {
    &mapHGTD.at(Score::TRKPTZ),
    &mapHGTD.at(Score::TEST_MISCL),
    &mapHGTD.at(Score::TEST_MISAS)
  };
  shapeComparisonPlot(
    TString::Format("%s/inclusive/shape_comparison_linscale.pdf", SAVE_DIR),
    false, -200, 200, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlot(
    TString::Format("%s/inclusive/shape_comparison_logscale.pdf", SAVE_DIR),
    true, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlotPair(
    TString::Format("%s/inclusive/shape_comparison.pdf", SAVE_DIR),
    -200, 200, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);


  // --- 2D heatmap: mean |Δt| as a function of cluster PU fraction × avg nHGTD hits ---
  // One page per score in RESO_SET.  Uses a dedicated canvas with extra right margin for
  // the colour-bar palette.  kBird palette: blue (low |Δt|) → yellow/red (high |Δt|).
  {
    TString fname2d = TString::Format("%s/inclusive/cluspufrac_vs_nhit_meandt.pdf", SAVE_DIR);
    TCanvas c2d("c2d_pufrac_nhit", "", 800, 650);
    gStyle->SetPalette(kBird);

    c2d.Print(fname2d + "[");
    for (auto* aObj : RESO_SET) {
      c2d.Clear();
      c2d.SetRightMargin(0.17);
      c2d.SetLeftMargin(0.15);
      c2d.SetBottomMargin(0.15);
      c2d.SetTopMargin(0.08);
      TProfile2D* h2 = aObj->prof2dPuFracVsNhit.get();
      h2->SetMinimum(0);
      h2->SetContour(50);
      h2->GetXaxis()->SetTitleOffset(1.2);
      h2->GetYaxis()->SetTitleOffset(1.4);
      h2->Draw("COLZ");
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.032);
      lat.DrawLatex(0.13, 0.93, TString::Format("%s  [%s]",
                    aObj->score.toString(), aObj->timetypeIDer.c_str()));
      c2d.Print(fname2d);
    }
    c2d.Print(fname2d + "]");
  }

  // --- 2D heatmap: σ(Δt) as a function of cluster PU fraction × avg nHGTD hits ---
  // TProfile2D "S" option stores per-bin std-dev in the error slot; we transfer it
  // to a plain TH2D for COLZ rendering (profile Draw("COLZ") shows the mean, not σ).
  {
    TString fnameSigma = TString::Format("%s/inclusive/cluspufrac_vs_nhit_sigmadt.pdf", SAVE_DIR);
    TCanvas c2ds("c2ds_pufrac_nhit", "", 800, 650);
    gStyle->SetPalette(kBird);

    c2ds.Print(fnameSigma + "[");
    for (auto* aObj : RESO_SET) {
      c2ds.Clear();
      c2ds.SetRightMargin(0.17);
      c2ds.SetLeftMargin(0.15);
      c2ds.SetBottomMargin(0.15);
      c2ds.SetTopMargin(0.08);

      TProfile2D* prof = aObj->prof2dPuFracVsNhitSigma.get();
      // Build a TH2D whose bin contents are the per-bin σ(Δt) from the profile
      TH2D hSigma("hSigma_tmp",
                  TString::Format("%s;Cluster PU Fraction;Avg. nHGTD Hits / Track;#sigma(#Delta t) [ps]",
                                  prof->GetTitle()),
                  prof->GetNbinsX(), prof->GetXaxis()->GetXmin(), prof->GetXaxis()->GetXmax(),
                  prof->GetNbinsY(), prof->GetYaxis()->GetXmin(), prof->GetYaxis()->GetXmax());
      for (int iX = 1; iX <= prof->GetNbinsX(); ++iX)
        for (int iY = 1; iY <= prof->GetNbinsY(); ++iY) {
          double n = prof->GetBinEntries(prof->GetBin(iX, iY));
          if (n > 1)
            hSigma.SetBinContent(iX, iY, prof->GetBinError(iX, iY));
        }
      hSigma.SetMinimum(0);
      hSigma.SetContour(50);
      hSigma.GetXaxis()->SetTitleOffset(1.2);
      hSigma.GetYaxis()->SetTitleOffset(1.4);
      hSigma.Draw("COLZ");
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.032);
      lat.DrawLatex(0.13, 0.93, TString::Format("%s  [%s]",
                    aObj->score.toString(), aObj->timetypeIDer.c_str()));
      c2ds.Print(fnameSigma);
    }
    c2ds.Print(fnameSigma + "]");
    gStyle->SetPalette(kBird);
  }

  // --- In-time PU diagnostic: Δt(cluster − HS truth) vs Δz(dominant PU vtx − HS vtx) ---
  // Three versions: inclusive, HIGH-quality (Q≥0.5), LOW-quality (Q<0.5).
  // Tests whether the 30 ps middle background is in-time pile-up. The "cross" shape
  // (vertical spoke at Δz=0, horizontal spoke at Δt=0) is the signature of in-time PU.
  // We overlay a reference slope line σ_t/σ_z ≈ 3.5 ps/mm (visible mostly as a sanity check).
  auto plotDtDz = [&](const TString& fname, auto getter) {
    TCanvas c2dz("c2dz_dt_dz_pu", "", 800, 650);
    gStyle->SetPalette(kBird);
    c2dz.Print(fname + "[");
    for (auto* aObj : RESO_SET) {
      c2dz.Clear();
      c2dz.SetRightMargin(0.17);
      c2dz.SetLeftMargin(0.15);
      c2dz.SetBottomMargin(0.15);
      c2dz.SetTopMargin(0.08);
      c2dz.SetLogz();

      TH2D* h2 = getter(aObj);
      h2->SetMinimum(1);
      h2->SetContour(50);
      h2->GetXaxis()->SetTitleOffset(1.2);
      h2->GetYaxis()->SetTitleOffset(1.4);
      h2->Draw("COLZ");

      const double slope = PILEUP_SMEAR / 50.0;
      TF1 fLine("fSlope", "[0]*x", -60.0, 60.0);
      fLine.SetParameter(0, slope);
      fLine.SetLineColor(kRed);
      fLine.SetLineWidth(2);
      fLine.SetLineStyle(2);
      fLine.DrawCopy("SAME");

      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.032);
      lat.DrawLatex(0.13, 0.93, TString::Format("%s  [%s]",
                    aObj->score.toString(), aObj->timetypeIDer.c_str()));
      lat.SetTextColor(kRed);
      lat.DrawLatex(0.55, 0.18, TString::Format("Expected slope = %.2f ps/mm", slope));
      c2dz.Print(fname);
    }
    c2dz.Print(fname + "]");
    c2dz.SetLogz(false);
  };

  plotDtDz(TString::Format("%s/inclusive/dt_vs_dz_pu.pdf",      SAVE_DIR),
           [](AnalysisObj* a) { return a->dtClusterVsDzPU.get();     });
  plotDtDz(TString::Format("%s/inclusive/dt_vs_dz_pu_high.pdf", SAVE_DIR),
           [](AnalysisObj* a) { return a->dtClusterVsDzPUHigh.get(); });
  plotDtDz(TString::Format("%s/inclusive/dt_vs_dz_pu_low.pdf",  SAVE_DIR),
           [](AnalysisObj* a) { return a->dtClusterVsDzPULow.get();  });


  // --- nHSTrack CDF / survival diagnostic ---
  // Shows what fraction of events fall below (CDF) or above (survival) each
  // nHSTrack threshold.  Answers: "is the nHSTrack < N region worth optimising?"
  {
    gSystem->mkdir(TString::Format("%s/diagnostics", SAVE_DIR), true);
    TH1D* hRaw = mapHGTD.at(Score::TRKPTZ).ptrHSTrack->effTotal.get();
    auto* hCDF = static_cast<TH1D*>(hRaw->GetCumulative());
    double total = hCDF->GetBinContent(hCDF->GetNbinsX());

    // Build CDF and survival TGraphs
    std::vector<double> xs, yCDF, ySurv;
    if (total > 0) {
      for (int b = 1; b <= hCDF->GetNbinsX(); ++b) {
        double cdf = hCDF->GetBinContent(b) / total;
        xs   .push_back(hCDF->GetBinCenter(b));
        yCDF .push_back(cdf);
        ySurv.push_back(1.0 - cdf);
      }
    }
    delete hCDF;

    int n = (int)xs.size();
    TGraph grCDF (n, xs.data(), yCDF .data());
    TGraph grSurv(n, xs.data(), ySurv.data());

    grCDF .SetLineColor(kBlue+1); grCDF .SetLineWidth(2);
    grSurv.SetLineColor(kRed+1);  grSurv.SetLineWidth(2);
    grCDF .SetMarkerColor(kBlue+1); grCDF .SetMarkerStyle(20); grCDF .SetMarkerSize(0.6);
    grSurv.SetMarkerColor(kRed+1);  grSurv.SetMarkerStyle(20); grSurv.SetMarkerSize(0.6);

    canvas->cd();
    grCDF.SetTitle(";N_{HS tracks} (forward);Fraction of events");
    grCDF.GetXaxis()->SetRangeUser(0, 20);
    grCDF.GetYaxis()->SetRangeUser(0, 1.2);
    grCDF.Draw("APL");
    grSurv.Draw("PL SAME");

    // Reference lines and labels at N=5 and N=10
    auto drawRef = [&](int N, double yC, double yS) {
      TLine* vl = new TLine(N, 0, N, yC);
      vl->SetLineStyle(2); vl->SetLineColor(kGray+1); vl->Draw();
      TLatex lat;
      lat.SetTextSize(0.028); lat.SetTextColor(kBlue+1);
      lat.DrawLatex(N + 0.2, yC + 0.01, TString::Format("%.1f%%", yC * 100));
      lat.SetTextColor(kRed+1);
      lat.DrawLatex(N + 0.2, yS - 0.03, TString::Format("%.1f%%", yS * 100));
    };
    // find CDF values at N=5 and N=10
    auto cdfAt = [&](int N) -> double {
      for (int i = 0; i < n; ++i)
        if ((int)std::round(xs[i]) == N) return yCDF[i];
      return -1;
    };
    for (int refN : {5, 10}) {
      double c = cdfAt(refN);
      if (c >= 0) drawRef(refN, c, 1.0 - c);
    }

    TLegend leg(0.60, 0.82, 0.80, 0.92);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
    leg.AddEntry(&grCDF,  "CDF: frac. with #leq N", "lp");
    leg.AddEntry(&grSurv, "Survival: frac. with > N", "lp");
    leg.Draw();

    ATLASLabel(0.18, 0.88, "Simulation Internal");
    canvas->SaveAs(TString::Format("%s/diagnostics/hs_track_cdf.pdf", SAVE_DIR));
    canvas->Clear();
  }

  std::cout << "FINISHED PLOT PRINTING\n";

  // --- Print event display commands (toggle PRINT_EVENT_DISPLAYS above) ---
  printEventDisplays("HGTD times: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayHGTD    );
  // printEventDisplays("Ideal Res.: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayIdealRes);
  // printEventDisplays("Ideal Eff.: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayIdealEff);

  // Shared helpers for all random-subsample blocks below.
  // Each block is guarded by PRINT_EVENT_DISPLAYS and can be independently
  // commented out without affecting the others.
  std::mt19937 rng(std::random_device{}());  // non-deterministic seed
  auto subsampleN = [&](std::vector<TString>& v, int n) {
    std::shuffle(v.begin(), v.end(), rng);
    if ((int)v.size() > n) v.resize(n);
  };
  auto printGroup = [](const char* header, const std::vector<TString>& cmds) {
    std::cout << "\n=== " << header << " (" << cmds.size() << " events) ===\n";
    for (const auto& c : cmds) std::cout << c << '\n';
  };

  // --- Low-multiplicity (comment out this block to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(lowMultPass, N_LOW_MULT_DISPLAY);
    subsampleN(lowMultFail, N_LOW_MULT_DISPLAY);
    printGroup(
      TString::Format("Low-mult (nFwdHS==%d, HGTD) — TRKPTZ PASS", LOW_MULT_NHS).Data(),
      lowMultPass);
    printGroup(
      TString::Format("Low-mult (nFwdHS==%d, HGTD) — TRKPTZ FAIL", LOW_MULT_NHS).Data(),
      lowMultFail);
  }

  // --- TRKPTZ / T_REFINED status unchanged (comment out this block to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(tRefinedUnchanged, N_T_REFINED_UNCHANGED);
    printGroup("T_REFINED unchanged — TRKPTZ↔T_REFINED both PASS or both FAIL",
               tRefinedUnchanged);
  }

  // --- TRKPTZ / T_REFINED status changed (comment out this block to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(tRefinedChanged, N_T_REFINED_UNCHANGED);
    printGroup("T_REFINED changed — TRKPTZ and T_REFINED DISAGREE (outcome flipped)",
               tRefinedChanged);
  }

  // --- TEST_MISCL passes: pure cluster, passes timing window (comment out to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(misclPassEvents, N_MISCL_PASS);
    printGroup("TEST_MISCL: pure cluster (purity>75%), PASSES timing window",
               misclPassEvents);
  }

  // --- TEST_MISCL not in filter: impure cluster, excluded from denominator (comment out to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(misclFailEvents, N_MISCL_PASS);
    printGroup("TEST_MISCL: EXCLUDED from filter (cluster purity ≤ 75%)",
               misclFailEvents);
  }

  // --- TEST_MISAS passes: clean HS timing, passes timing window (comment out to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(misasPassEvents, N_MISAS_PASS);
    printGroup("TEST_MISAS: clean HS timing (hsTimingPurity≥95%), PASSES timing window",
               misasPassEvents);
  }

  // --- TEST_MISAS not in filter: dirty HS timing, excluded from denominator (comment out to disable) ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(misasFailEvents, N_MISAS_PASS);
    printGroup("TEST_MISAS: EXCLUDED from filter (hsTimingPurity < 95%)",
               misasFailEvents);
  }

  // --- Cluster quality tiers — browse events at two levels of cluster quality ---
  if (PRINT_EVENT_DISPLAYS) {
    subsampleN(cqHighEvents, N_LOW_MULT_DISPLAY);
    subsampleN(cqLowEvents,  N_LOW_MULT_DISPLAY);
    printGroup(TString::Format("Cluster quality HIGH (Q>=%.1f) — expected core σ~10 ps",
                               CLUS_QUALITY_SPLIT).Data(), cqHighEvents);
    printGroup(TString::Format("Cluster quality LOW  (Q<%.1f) — expected background σ~175 ps",
                               CLUS_QUALITY_SPLIT).Data(), cqLowEvents);
  }

  return 0;
}
