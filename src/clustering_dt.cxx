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
  // Focus on the timing-misassignment oracle: HGTD, TRKPTZ, WAVeS, and WAVeS
  // with ideal timing (WAVES_MISAS).  Misclustering case removed.
  // Color override: WAVES yellow, ideal timing violet.
  moneyPlot(TString::Format("%s/waves_oracle_%s.pdf",compSubdir.c_str(), key), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::WAVES_MISAS),
	    },
	    {C01, C02, C03, C04});

  // Simple three-way comparison: HGTD vs TRKPTZ vs WAVeS (no oracles).
  moneyPlot(TString::Format("%s/hgtd_trkptz_waves_%s.pdf",compSubdir.c_str(), key), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::WAVES),
	    },
	    {C01, C02, C03});

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
                &mapHGTD.at(Score::WAVES),
		&mapHGTD.at(Score::ZT_ITER),
	    });
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

auto main() -> int {
  SetAtlasStyle();  // ATLAS publication style (fonts, margins, ticks, …)

  // --- Data source ---
  TChain chain("ntuple");
  setupChain(chain, "../../highstats-ntuple/");
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

  const auto KEYS = {"pu_frac", "fjet", "ftrack", "hs_track", "truthjets"};

  // --- Comparison plots (per variable KEY) ---
  for (const auto* key : KEYS)
    makeComparisonPlots(key, canvas, mapHGTD, mapIdealRes, mapIdealEff);

  // --- Inclusive resolution plots ---
  const std::initializer_list<AnalysisObj*> RESO_SET = {
    &mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::TRKPTZ), &mapHGTD.at(Score::WAVES),
    &mapHGTD.at(Score::WAVES_MISAS), &mapHGTD.at(Score::TEST_ML),
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

  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR),
		true,  false, -400, 400, canvas, RESO_SET, resoGetter, &FIT_TRPGAUS, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR),
		false, false, -200, 200, canvas, RESO_SET, resoGetter, &FIT_TRPGAUS, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);


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

  // --- Normalized shape comparison: TRKPTZ vs oracle scores ---
  const std::vector<AnalysisObj*> SHAPE_SET = {
    &mapHGTD.at(Score::TRKPTZ),
    &mapHGTD.at(Score::WAVES),
    &mapHGTD.at(Score::WAVES_MISAS)
  };
  shapeComparisonPlot(
    TString::Format("%s/diagnostics/shape_comparison_linscale.pdf", SAVE_DIR),
    false, -200, 200, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlot(
    TString::Format("%s/diagnostics/shape_comparison_logscale.pdf", SAVE_DIR),
    true, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlotPair(
    TString::Format("%s/diagnostics/shape_comparison.pdf", SAVE_DIR),
    -200, 200, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);

  std::cout << "FINISHED PLOT PRINTING\n";

  // --- Print event display commands (toggle PRINT_EVENT_DISPLAYS above) ---
  printEventDisplays("HGTD times: TRKPTZ passes, TEST_MISAS fails (misassignment)", evtDisplayHGTD    );

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

  return 0;
}
