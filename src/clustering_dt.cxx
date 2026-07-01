#include <atomic>
#include <memory>
#include <mutex>
#include <random>
#include <RtypesCore.h>
#include <TROOT.h>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include "clustering_constants.h"
#include "sample_config.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

// Output directory for saved plots (mirrors the active --sample; see main()).
static const std::string& SAVE_DIR = MyUtl::OUTPUT_DIR;

// Event display command template (filled with file number, event number, extra time)
static constexpr const char* EVTDISPLAY_FMT =
  "cd python && python3 event_display.py --file_num %s --event_num %lld --extra_time %.2f";

// Set to true to print event display commands to stdout after the event loop.
static constexpr bool PRINT_EVENT_DISPLAYS = false;

// Max event-display commands to print per low-multiplicity category.
static constexpr int N_LOW_MULT_DISPLAY = 20;
// HS-track count defining "low multiplicity" for the event-display sample.
static constexpr int LOW_MULT_NHS = 5;

// Max commands for the MISAS diagnostic category.
// (Gated on PRINT_EVENT_DISPLAYS.)
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
  // m.emplace(Score::PASS,       AnalysisObj(label, Score::PASS      ));
  m.emplace(Score::TEST_MISAS, AnalysisObj(label, Score::TEST_MISAS));
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
  std::map<Score, AnalysisObj>& mapHGTD
) {
  const std::string compSubdir = SAVE_DIR + "/comparisons";

  moneyPlot(TString::Format("%s/theplot_%s.pdf",compSubdir.c_str(), key), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
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

}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

auto main(int argc, char** argv) -> int {
  SetAtlasStyle();  // ATLAS publication style (fonts, margins, ticks, …)

  // Must precede any histogram construction: prevents concurrent TH1
  // self-registration into gDirectory from racing across worker threads
  // below. Every worker thread's buildAnalysisMap() produces identically-
  // named histograms, so removing this would silently reintroduce a
  // gDirectory name collision across threads.
  TH1::AddDirectory(kFALSE);

  // --- Sample selection (--sample=vbf|zjets|dijet; default: local VBF ntuple) ---
  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  for (const char* sub : {"comparisons", "inclusive", "fullplots", "diagnostics"})
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/" + sub);
  unsigned nThreads = MyUtl::resolveThreads(argc, argv);

  // --- Data source ---
  TChain chain("ntuple");
  setupChain(chain, sample.ntupleDir.c_str());
  ROOT::EnableImplicitMT(nThreads);

  gErrorIgnoreLevel = kFatal;

  // --- Canvas ---
  TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);

  // --- Per-thread analysis-map registry, merged into one mapHGTD after the
  //     event loop. AnalysisObj/PlotObj are move-only (unique_ptr members),
  //     so this can't use ROOT::TThreadedObject (which clones its model via
  //     copy-ctor or TObject::Clone(), neither of which applies here) --
  //     each worker thread instead registers its own lazily-built map into
  //     this mutex-guarded vector exactly once, then caches a raw pointer
  //     to it in thread_local storage for lock-free access on every event.
  //
  //     The mutex also serializes the buildAnalysisMap() call itself, not
  //     just the push_back: AnalysisObj's constructor calls
  //     SetFillColorAlpha(), which lazily creates+registers a new TColor
  //     into ROOT's *global* color table the first time a given
  //     (color, alpha) pair is used, and that registration is not
  //     thread-safe. Two worker threads racing to build their first
  //     AnalysisObj map concurrently corrupted the heap (confirmed via a
  //     crash report: SIGABRT inside TColor::GetColorTransparent, malloc
  //     "pointer being freed was not allocated"). Serializing the whole
  //     construction means only one thread ever runs it at a time -- by the
  //     next thread's turn, the colors it needs already exist and no race
  //     occurs. This only costs time once per thread (lazy init), not per
  //     event. ---
  std::mutex mapRegistryMutex;
  std::vector<std::unique_ptr<std::map<Score, AnalysisObj>>> mapRegistry;

  // --- Per-thread event-display candidate lists. TString/std::vector are
  //     genuinely deep-copyable, so ROOT::TThreadedObject's clone-via-copy-
  //     ctor is correct here (unlike the AnalysisObj map above). ---
  ROOT::TThreadedObject<std::vector<TString>> tsEvtDisplayHGTD;
  ROOT::TThreadedObject<std::vector<TString>> tsLowMultPass;
  ROOT::TThreadedObject<std::vector<TString>> tsLowMultFail;
  ROOT::TThreadedObject<std::vector<TString>> tsMisasPassEvents;
  ROOT::TThreadedObject<std::vector<TString>> tsMisasFailEvents;

  std::atomic<Long64_t> progressCounter{0};

  // --- Event loop ---
  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  ROOT::TTreeProcessorMT proc(chain, nThreads);
  proc.Process([&](TTreeReader& reader) {
    // Fresh per invocation: a worker thread can be handed a different
    // TTreeReader across task ranges, so this cannot be thread_local.
    BranchPointerWrapper branch(reader);

    // Lazily build this thread's analysis map once; reused across however
    // many task ranges this worker thread services.
    thread_local std::map<Score, AnalysisObj>* tlMap = nullptr;
    if (!tlMap) {
      std::lock_guard<std::mutex> lock(mapRegistryMutex);
      mapRegistry.push_back(std::make_unique<std::map<Score, AnalysisObj>>(buildAnalysisMap(Scenario::HGTD)));
      tlMap = mapRegistry.back().get();
    }

    // Resolve this thread's event-display slot once per invocation (not per
    // event -- TThreadedObject::Get() involves a thread-id lookup).
    auto evtDisplayHGTD  = tsEvtDisplayHGTD.Get();
    auto lowMultPass     = tsLowMultPass.Get();
    auto lowMultFail     = tsLowMultFail.Get();
    auto misasPassEvents = tsMisasPassEvents.Get();
    auto misasFailEvents = tsMisasFailEvents.Get();

    while (reader.Next()) {
      Long64_t n = ++progressCounter;
      if (n % 5000 == 0)
        std::cout << "Progress: " << n << "/" << N_EVENT << "\r" << std::flush;

      // HGTD timing scenario (the only scenario currently active)
      auto resHGTD = processEventData(&branch, false, true, *tlMap);

      // Extract file identifier from the full path (characters 49–54)
      TString fileName = reader.GetTree()->GetCurrentFile()->GetName();
      TString fileNum  = fileName(49, 6);
      // Local (per-file) entry number: reader.GetTree() is the currently-
      // loaded per-file constituent tree here (no outer TChain is available
      // inside this lambda to replicate the sequential-loop's
      // chain.GetReadEntry()-chain.GetChainOffset() computation), and
      // GetReadEntry() on it gives the correct local number directly --
      // validated against a sequential baseline in util/ttreemt_prototype.cxx.
      Long64_t EVENT_NUM = reader.GetTree()->GetReadEntry();

      // Collect events where TRKPTZ passes but TEST_MISAS does not (misassignment effect)
      collectEventDisplay(*evtDisplayHGTD, 3, resHGTD, fileNum, EVENT_NUM);

      // Low-multiplicity event display collection (HGTD scenario, n=LOW_MULT_NHS HS tracks)
      if (resHGTD.code >= 0 && resHGTD.nFwdHS == LOW_MULT_NHS) {
        TString cmd = TString::Format(EVTDISPLAY_FMT, fileNum.Data(), EVENT_NUM, resHGTD.time);
        (resHGTD.trkptzPass ? *lowMultPass : *lowMultFail).push_back(cmd);
      }

      if (resHGTD.code >= 0) {
        TString cmd = TString::Format(EVTDISPLAY_FMT, fileNum.Data(), EVENT_NUM, resHGTD.time);

        // Category: event is in TEST_MISAS denominator (clean HS timing) and PASSES
        if (resHGTD.misasPass)
          misasPassEvents->push_back(cmd);
        // Category fail: event did NOT enter the TEST_MISAS denominator (hsTimingPurity < 95%)
        if (!resHGTD.misasInDenom)
          misasFailEvents->push_back(cmd);
      }
    }
  });
  std::cout << "\n";

  // --- Merge per-thread analysis maps into one ---
  std::map<Score, AnalysisObj> mapHGTD;
  if (mapRegistry.empty()) {
    mapHGTD = buildAnalysisMap(Scenario::HGTD);
  } else {
    mapHGTD = std::move(*mapRegistry.front());
    for (size_t i = 1; i < mapRegistry.size(); ++i)
      for (auto& [score, analysis] : mapHGTD)
        analysis.mergeFrom(mapRegistry[i]->at(score));
  }

  auto allMaps = { &mapHGTD };

  // --- Merge per-thread event-display candidate lists ---
  auto mergeVec = [](std::vector<TString>& dest, ROOT::TThreadedObject<std::vector<TString>>& src) {
    for (unsigned i = 0; i < src.GetNSlots(); ++i) {
      auto v = src.GetAtSlot(i);
      if (v) dest.insert(dest.end(), v->begin(), v->end());
    }
  };
  std::vector<TString> evtDisplayHGTD, lowMultPass, lowMultFail, misasPassEvents, misasFailEvents;
  mergeVec(evtDisplayHGTD,  tsEvtDisplayHGTD);
  mergeVec(lowMultPass,     tsLowMultPass);
  mergeVec(lowMultFail,     tsLowMultFail);
  mergeVec(misasPassEvents, tsMisasPassEvents);
  mergeVec(misasFailEvents, tsMisasFailEvents);

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
    makeComparisonPlots(key, canvas, mapHGTD);

  // --- Inclusive resolution plots ---
  const std::initializer_list<AnalysisObj*> RESO_SET = {
    &mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::TRKPTZ), &mapHGTD.at(Score::WAVES),
    &mapHGTD.at(Score::WAVES_MISAS),
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
      stackMax(a->inclusiveResoNhit3pSig.get(), a->inclusiveResoNhit3pMix.get(), a->inclusiveResoNhit3pBkg.get())});
    globalPullMax = std::max({globalPullMax,
      stackMax(a->inclusivePullSig         .get(), a->inclusivePullMix         .get(), a->inclusivePullBkg         .get())});
  }
  const double INCL_RESO_YMAX_LOG = 5.0 * globalResoMax;
  const double INCL_RESO_YMAX_LIN = 1.1 * globalResoMax;
  const double INCL_PULL_YMAX_LOG = 5.0 * globalPullMax;
  const double INCL_PULL_YMAX_LIN = 1.1 * globalPullMax;
  const double INCL_YMIN_LOG      = 0.5;

  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_logscale.pdf", SAVE_DIR.c_str()),
		true,  false, -400, 400, canvas, RESO_SET, resoGetter, &FIT_TRPGAUS, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivereso_linscale.pdf", SAVE_DIR.c_str()),
		false, false, -200, 200, canvas, RESO_SET, resoGetter, &FIT_TRPGAUS, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);


  auto pullGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusivePullSig.get(),
             a->inclusivePullMix.get(),
             a->inclusivePullBkg.get() }; };
  const char* pullLabel = "(t_{0}-t_{truth})/#sigma_{t}";
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_logscale.pdf", SAVE_DIR.c_str()),
		true,  false, -10, 10, canvas, RESO_SET,
		pullGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(TString::Format("%s/inclusive/inclusivepull_linscale.pdf", SAVE_DIR.c_str()),
		false, false, -10, 10, canvas, RESO_SET,
		pullGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LIN, 0.0);

  // --- Normalized shape comparison: TRKPTZ vs oracle scores ---
  const std::vector<AnalysisObj*> SHAPE_SET = {
    &mapHGTD.at(Score::TRKPTZ),
    &mapHGTD.at(Score::WAVES),
    &mapHGTD.at(Score::WAVES_MISAS)
  };
  shapeComparisonPlot(
    TString::Format("%s/inclusive/shape_comparison_linscale.pdf", SAVE_DIR.c_str()),
    false, -200, 200, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlot(
    TString::Format("%s/inclusive/shape_comparison_logscale.pdf", SAVE_DIR.c_str()),
    true, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlotPair(
    TString::Format("%s/inclusive/shape_comparison.pdf", SAVE_DIR.c_str()),
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
