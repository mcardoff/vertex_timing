#include <atomic>
#include <memory>
#include <mutex>
#include <random>
#include <RtypesCore.h>
#include <TROOT.h>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include "clustering_constants.h"
#include "clustering_dt_common.h"
#include "sample_config.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "histogram_io.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// clustering_hist.cxx
//   Histogramming stage of the clustering_dt split: runs the full
//   TTreeProcessorMT event loop and merge exactly as the former monolithic
//   clustering_dt.cxx did, then writes every histogram to
//   <OUTPUT_DIR>/hists/clustering_hist.root via histogram_io.h instead of
//   producing PDFs directly. clustering_plot.cxx reads that file back and
//   does all plotting -- see CLAUDE.md's "Main Executables" section.
// ---------------------------------------------------------------------------

// Event display command template (filled with full file path, event number, extra time).
// Quoted since --file_path takes an arbitrary absolute path.
static constexpr const char* EVTDISPLAY_FMT =
  "cd python && python3 event_display.py --file_path \"%s\" --event_num %lld --extra_time %.2f";

// Set to true to print event display commands to stdout after the event loop.
static constexpr bool PRINT_EVENT_DISPLAYS = true;

// Max event-display commands to print per low-multiplicity category.
static constexpr int N_LOW_MULT_DISPLAY = 20;
// HS-track count defining "low multiplicity" for the event-display sample.
static constexpr int LOW_MULT_NHS = 5;

// Max commands for the MISAS diagnostic category.
// (Gated on PRINT_EVENT_DISPLAYS.)
static constexpr int N_MISAS_PASS          =  5;  // clean HS timing, PASSES timing window

// ---------------------------------------------------------------------------
// Helper: add an event-display command string to a list if condition is met.
// ---------------------------------------------------------------------------

void collectEventDisplay(
  std::vector<TString>& list,
  int returnCode,
  const EventResult& result,
  const TString& filePath,
  Long64_t eventNum
) {
  if (result.code == returnCode)
    list.push_back(TString::Format(EVTDISPLAY_FMT, filePath.Data(), eventNum, result.time));
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
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  MyUtl::OVERLAP_REMOVAL = sample.overlapRemoval;  // Z+jets lepton–jet overlap removal
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  if (MyUtl::SAMPLE_NAME.empty())
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/hists");
  unsigned nThreads = MyUtl::resolveThreads(argc, argv);

  // --- Data source ---
  TChain chain("ntuple");
  setupChain(chain, sample.ntupleDir.c_str());
  ROOT::EnableImplicitMT(nThreads);

  gErrorIgnoreLevel = kFatal;

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

  // --- Event-selection rejection-cause tally, keyed on EventResult::code (see
  //     its definition in event_processing.h). Plain atomics suffice here --
  //     no lazy per-thread construction needed, unlike AnalysisObj/TColor.
  //     Lets a low post-selection yield (e.g. on Z+jets) be diagnosed from the
  //     printed breakdown instead of guessed at. ---
  std::atomic<Long64_t> nRejNoLepton{0}, nRejOneLepton{0}, nRejNoOSSFPair{0},
                         nRejOverlapVeto{0}, nRejBasicCuts{0}, nRejJetPtCut{0},
                         nRejOther{0}, nAccepted{0};

  std::atomic<Long64_t> progressCounter{0};

  // --- Event loop ---
  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  // Optional --max-events=<N> cap for quick local checks (see resolveMaxEvents
  // in sample_config.h). -1 means unlimited -- every existing invocation
  // without the flag is unaffected. progressDenom is just for the progress
  // print below; the actual stopping condition is the per-event check inside
  // the loop, keyed off the same shared progressCounter every worker thread
  // increments, so all threads converge on the cap together.
  const Long64_t maxEvents    = MyUtl::resolveMaxEvents(argc, argv);
  const Long64_t progressDenom = (maxEvents > 0) ? std::min(maxEvents, N_EVENT) : N_EVENT;
  if (maxEvents > 0)
    std::cout << "Restricting to first " << progressDenom << " events (--max-events)\n";

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
      // Stop this worker's current task range once the shared counter has
      // already crossed the cap -- other in-flight/future task ranges hit
      // the same check within their own next iteration, so overshoot across
      // threads is bounded by ~nThreads events, not the whole remaining tree.
      if (maxEvents > 0 && n > maxEvents) break;
      if (n % 5000 == 0)
        std::cout << "Progress: " << n << "/" << progressDenom << "\r" << std::flush;

      // HGTD timing scenario (the only scenario currently active)
      auto resHGTD = processEventData(&branch, false, true, *tlMap);

      // Tally rejection cause (see EventResult::code in event_processing.h).
      switch (resHGTD.code) {
        case -2: ++nRejNoLepton;    break;
        case -3: ++nRejOverlapVeto; break;
        case -4: ++nRejBasicCuts;   break;
        case -5: ++nRejJetPtCut;    break;
        case -6: ++nRejOneLepton;   break;
        case -7: ++nRejNoOSSFPair;  break;
        default: (resHGTD.code < 0 ? nRejOther : nAccepted)++; break;
      }

      // Full ntuple file path -- passed to event_display.py as-is (--file_path),
      // which works uniformly across every sample (vbf/zjets/dijet/local),
      // unlike the old fixed-offset substring (fileName(49,6)) that only
      // matched the local default VBF ntuple's one absolute path length.
      TString filePath = reader.GetTree()->GetCurrentFile()->GetName();
      // Local (per-file) entry number: reader.GetTree() is the currently-
      // loaded per-file constituent tree here (no outer TChain is available
      // inside this lambda to replicate the sequential-loop's
      // chain.GetReadEntry()-chain.GetChainOffset() computation), and
      // GetReadEntry() on it gives the correct local number directly --
      // validated against a sequential baseline in util/ttreemt_prototype.cxx.
      Long64_t EVENT_NUM = reader.GetTree()->GetReadEntry();

      // Collect events where TRKPTZ passes but TEST_MISAS does not (misassignment effect)
      collectEventDisplay(*evtDisplayHGTD, 3, resHGTD, filePath, EVENT_NUM);

      // Low-multiplicity event display collection (HGTD scenario, n=LOW_MULT_NHS HS tracks)
      if (resHGTD.code >= 0 && resHGTD.nFwdHS == LOW_MULT_NHS) {
        TString cmd = TString::Format(EVTDISPLAY_FMT, filePath.Data(), EVENT_NUM, resHGTD.time);
        (resHGTD.trkptzPass ? *lowMultPass : *lowMultFail).push_back(cmd);
      }

      if (resHGTD.code >= 0) {
        TString cmd = TString::Format(EVTDISPLAY_FMT, filePath.Data(), EVENT_NUM, resHGTD.time);

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

  std::cout << "\nFINISHED PROCESSING\n";

  // Actual events processed -- equals N_EVENT when --max-events was not given;
  // otherwise reflects the (slightly-overshoot-bounded) cap. Read only after
  // proc.Process() has returned, so every worker thread's increments are in.
  const Long64_t nEventsProcessed = progressCounter.load();

  // --- Event-selection rejection-cause breakdown ---
  std::cout << "\n=== EVENT SELECTION BREAKDOWN (of " << nEventsProcessed << " total) ===\n";
  std::cout << "  Accepted                        : " << nAccepted       << '\n';
  std::cout << "  Rejected: Z->ll, 0 good leptons  : " << nRejNoLepton    << '\n';
  std::cout << "  Rejected: Z->ll, 1 good lepton   : " << nRejOneLepton   << '\n';
  std::cout << "  Rejected: Z->ll, no OS-SF pair   : " << nRejNoOSSFPair  << '\n';
  std::cout << "  Rejected: overlap veto (SKIP)    : " << nRejOverlapVeto << '\n';
  std::cout << "  Rejected: vertex/MIN_JETS        : " << nRejBasicCuts   << '\n';
  std::cout << "  Rejected: jet-pT/VBS topology    : " << nRejJetPtCut    << '\n';
  std::cout << "  Rejected: other/undistinguished  : " << nRejOther       << '\n';

  // --- Save every histogram to a ROOT file for the plotting stage ---
  const std::string histPath = MyUtl::histFilePath("clustering_hist.root");
  MyUtl::HistWriter writer(histPath);
  for (auto& [score, analysis] : mapHGTD) analysis.saveTo(writer);
  writer.WriteScalar("meta_n_accepted",         static_cast<Long64_t>(nAccepted));
  writer.WriteScalar("meta_n_rej_no_lepton",    static_cast<Long64_t>(nRejNoLepton));
  writer.WriteScalar("meta_n_rej_one_lepton",   static_cast<Long64_t>(nRejOneLepton));
  writer.WriteScalar("meta_n_rej_no_ossf_pair", static_cast<Long64_t>(nRejNoOSSFPair));
  writer.WriteScalar("meta_n_rej_overlap_veto", static_cast<Long64_t>(nRejOverlapVeto));
  writer.WriteScalar("meta_n_rej_basic_cuts",   static_cast<Long64_t>(nRejBasicCuts));
  writer.WriteScalar("meta_n_rej_jetpt_cut",    static_cast<Long64_t>(nRejJetPtCut));
  writer.WriteScalar("meta_n_rej_other",        static_cast<Long64_t>(nRejOther));
  writer.WriteRunMeta(MyUtl::ENERGY_LABEL, nEventsProcessed);
  writer.Close();
  std::cout << "Wrote histograms to " << histPath << '\n';

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
