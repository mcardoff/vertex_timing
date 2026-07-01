// ---------------------------------------------------------------------------
// ttreemt_prototype.cxx
//
// Exploratory prototype: does ROOT's TTreeProcessorMT let us parallelize the
// per-event clustering/scoring work in processEventData() itself (not just
// ROOT-internal I/O, which EnableImplicitMT() already covers)?
//
// This runs the SAME event selection/clustering as clustering_dt.cxx's HGTD
// scenario twice against the same sample:
//   1. a plain sequential TTreeReader loop (baseline), and
//   2. a TTreeProcessorMT-driven parallel loop,
// then compares per-event EventResult fields for EXACT equality (the HGTD
// scenario has no RNG and no cross-event state, so results should be
// bit-identical, not just statistically similar), and reports wall-clock
// speedup.
//
// This file is intentionally standalone: it does not modify clustering_dt.cxx
// or any shared header, and is not wired into the default analysis. See
// NOTES_ttreeprocessormt.md / CLAUDE.md for the surrounding investigation.
//
// Histogram merging across threads is explicitly out of scope here -- the
// per-thread AnalysisObj maps built below are required by processEventData's
// signature but are otherwise discarded; they only prove that concurrent
// histogram construction doesn't crash/collide (TH1::AddDirectory(kFALSE)
// is required for that -- see below).
// ---------------------------------------------------------------------------

#include <chrono>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>

#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "sample_config.h"

using namespace MyUtl;

namespace {

  // -------------------------------------------------------------------------
  // Trimmed copy of clustering_dt.cxx's buildAnalysisMap(Scenario::HGTD) --
  // duplicated here (not exported from clustering_dt.cxx) so this prototype
  // stays fully self-contained and clustering_dt.cxx remains untouched.
  // -------------------------------------------------------------------------
  auto buildAnalysisMap() -> std::map<Score, AnalysisObj> {
    const char* label = "HGTD Times";
    std::map<Score, AnalysisObj> m;
    m.emplace(Score::TRKPTZ,        AnalysisObj(label, Score::TRKPTZ));
    m.emplace(Score::WAVES,         AnalysisObj(label, Score::WAVES));
    m.emplace(Score::JET_T_REFINED, AnalysisObj(label, Score::JET_T_REFINED));
    m.emplace(Score::TEST_MISAS,    AnalysisObj(label, Score::TEST_MISAS));
    m.emplace(Score::WAVES_MISCL,   AnalysisObj(label, Score::WAVES_MISCL));
    m.emplace(Score::WAVES_MISAS,   AnalysisObj(label, Score::WAVES_MISAS));
    m.emplace(Score::HGTD,          AnalysisObj(label, Score::HGTD));
    return m;
  }

  // -------------------------------------------------------------------------
  // Robust per-event key: current file's basename + local entry number within
  // that file's tree. TTreeReader::GetTree() returns the currently-loaded
  // per-file constituent tree both for the sequential TChain-backed reader
  // and for whatever tree TTreeProcessorMT hands each task, so GetReadEntry()
  // on it gives the local (not chain-global) entry number in both cases --
  // this sidesteps the production loop's chain.GetReadEntry()-chain.GetChainOffset()
  // arithmetic, which depends on an outer TChain object not available inside
  // the MT lambda. Confirming this key resolves consistently per-slot is the
  // core thing this prototype validates.
  // -------------------------------------------------------------------------
  struct EventKey {
    std::string file;
    Long64_t    entry;
    bool operator==(const EventKey& o) const { return entry == o.entry && file == o.file; }
  };
  struct EventKeyHash {
    size_t operator()(const EventKey& k) const {
      return std::hash<std::string>{}(k.file) ^ (std::hash<Long64_t>{}(k.entry) << 1);
    }
  };

  // For the TTreeProcessorMT pass: each task's TTreeReader is bound to a
  // per-task view with no outer TChain accessible, so GetTree()->GetReadEntry()
  // is the only handle available -- and (confirmed empirically below) it does
  // give the local, per-file entry number in that context.
  EventKey makeKey(TTreeReader& reader) {
    TFile* f = reader.GetTree()->GetCurrentFile();
    return { gSystem->BaseName(f->GetName()), reader.GetTree()->GetReadEntry() };
  }

  // For the sequential pass: the reader is bound directly to the top-level
  // TChain, and TTreeReader::GetTree() returns a handle whose GetReadEntry()
  // reflects the chain-GLOBAL entry number, not the local one -- confirmed by
  // a first run of this prototype, where sequential keys carried local-looking
  // entry numbers up to N_EVENT-1 all attributed to the last-iterated file.
  // This mirrors exactly why clustering_dt.cxx's main() computes EVENT_NUM via
  // chain.GetReadEntry() - chain.GetChainOffset() instead of using the reader's
  // tree directly; reproduce that here so both passes key on the same local
  // per-file entry numbering.
  EventKey makeKeySeq(TChain& chain, TTreeReader& reader) {
    TFile* f = reader.GetTree()->GetCurrentFile();
    Long64_t localEntry = chain.GetReadEntry() - chain.GetChainOffset();
    return { gSystem->BaseName(f->GetName()), localEntry };
  }

  using ResultMap = std::unordered_map<EventKey, EventResult, EventKeyHash>;

  bool sameResult(const EventResult& a, const EventResult& b) {
    return a.code         == b.code
        && a.time         == b.time
        && a.nFwdHS       == b.nFwdHS
        && a.trkptzPass   == b.trkptzPass
        && a.misasInDenom == b.misasInDenom
        && a.misasPass    == b.misasPass
        && a.clusQuality  == b.clusQuality;
  }

  void printResult(const char* label, const EventResult& r) {
    std::cout << "  " << label << ": code=" << r.code << " time=" << r.time
              << " nFwdHS=" << r.nFwdHS << " trkptzPass=" << r.trkptzPass
              << " misasInDenom=" << r.misasInDenom << " misasPass=" << r.misasPass
              << " clusQuality=" << r.clusQuality << "\n";
  }

}  // namespace

int main(int argc, char** argv) {
  auto sample = MyUtl::resolveSample(argc, argv);
  unsigned nThreads = MyUtl::resolveThreads(argc, argv);

  // Must precede any histogram construction: prevents concurrent TH1
  // self-registration into gDirectory from racing across worker threads.
  TH1::AddDirectory(kFALSE);
  gErrorIgnoreLevel = kFatal;

  std::cout << "Sample:            " << sample.ntupleDir << "\n";
  std::cout << "Threads requested: " << nThreads << "\n\n";

  // --- Pass 1: sequential baseline (mirrors clustering_dt.cxx's HGTD loop) --
  ResultMap baseline;
  double tSeq = 0.0;
  {
    TChain chain("ntuple");
    setupChain(chain, sample.ntupleDir.c_str());
    TTreeReader reader(&chain);
    BranchPointerWrapper branch(reader);
    auto analyses = buildAnalysisMap();

    auto t0 = std::chrono::steady_clock::now();
    Long64_t n = 0;
    while (reader.Next()) {
      auto res = processEventData(&branch, /*useSmearedTimes=*/false,
                                   /*checkValidTimes=*/true, analyses);
      baseline.emplace(makeKeySeq(chain, reader), res);
      ++n;
    }
    auto t1 = std::chrono::steady_clock::now();
    tSeq = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Sequential pass: " << n << " events in " << tSeq << " s\n";
  }

  // --- Pass 2: TTreeProcessorMT --------------------------------------------
  ResultMap parallel;
  double tPar = 0.0;
  {
    ROOT::EnableImplicitMT(nThreads);

    TChain chain("ntuple");
    setupChain(chain, sample.ntupleDir.c_str());

    ROOT::TThreadedObject<std::vector<std::pair<EventKey, EventResult>>> tsResults;
    ROOT::TTreeProcessorMT proc(chain, nThreads);

    // Guards AnalysisObj construction below -- AnalysisObj's constructor
    // calls SetFillColorAlpha(), which lazily creates+registers a new TColor
    // into ROOT's *global* color table the first time a given (color, alpha)
    // pair is used. That registration isn't thread-safe: two worker threads
    // racing to build their first AnalysisObj map concurrently can corrupt
    // the heap (confirmed via a crash report: SIGABRT inside
    // TColor::GetColorTransparent, malloc "pointer being freed was not
    // allocated"). Serializing the whole construction, not just the registry
    // push_back, means only one thread ever runs it at a time -- by the next
    // thread's turn, the colors it needs already exist and no race occurs.
    std::mutex analysisBuildMutex;

    auto t0 = std::chrono::steady_clock::now();
    proc.Process([&](TTreeReader& reader) {
      // Fresh BranchPointerWrapper per invocation: TTreeProcessorMT may hand
      // the same worker thread a different TTreeReader instance across task
      // ranges, so a stale thread_local BranchPointerWrapper bound to a
      // previous reader would be invalid. Binding is just pointer setup, so
      // this is cheap.
      BranchPointerWrapper branch(reader);

      // Histogram accumulation state IS safe to keep thread_local: it isn't
      // tied to any particular TTreeReader, and rebuilding ~7 AnalysisObj
      // (~60 histograms) on every task range would be wasteful.
      thread_local std::unique_ptr<std::map<Score, AnalysisObj>> tlAnalyses;
      if (!tlAnalyses) {
        std::lock_guard<std::mutex> lock(analysisBuildMutex);
        tlAnalyses = std::make_unique<std::map<Score, AnalysisObj>>(buildAnalysisMap());
      }

      auto resVec = tsResults.Get();  // resolve this thread's slot once, not per-event
      while (reader.Next()) {
        auto res = processEventData(&branch, false, true, *tlAnalyses);
        resVec->emplace_back(makeKey(reader), res);
      }
    });
    auto t1 = std::chrono::steady_clock::now();
    tPar = std::chrono::duration<double>(t1 - t0).count();

    Long64_t n = 0;
    for (unsigned i = 0; i < tsResults.GetNSlots(); ++i) {
      auto v = tsResults.GetAtSlot(i);
      if (!v) continue;
      for (auto& [k, r] : *v) {
        parallel.emplace(k, r);
        ++n;
      }
    }
    std::cout << "Parallel pass:   " << n << " events in " << tPar << " s ("
              << nThreads << " threads)\n";
  }

  // --- Compare ---------------------------------------------------------------
  std::cout << "\n=== Comparison ===\n";
  bool ok = true;

  if (baseline.size() != parallel.size()) {
    ok = false;
    std::cout << "EVENT COUNT MISMATCH: sequential=" << baseline.size()
              << " parallel=" << parallel.size() << "\n";
  }

  int nMismatch = 0;
  const int MAX_PRINT = 10;
  for (auto& [key, seqRes] : baseline) {
    auto it = parallel.find(key);
    if (it == parallel.end()) {
      ok = false;
      ++nMismatch;
      if (nMismatch <= MAX_PRINT)
        std::cout << "MISSING in parallel pass: " << key.file << ":" << key.entry << "\n";
      continue;
    }
    if (!sameResult(seqRes, it->second)) {
      ok = false;
      ++nMismatch;
      if (nMismatch <= MAX_PRINT) {
        std::cout << "MISMATCH " << key.file << ":" << key.entry << "\n";
        printResult("seq", seqRes);
        printResult("par", it->second);
      }
    }
  }

  if (ok)
    std::cout << baseline.size() << "/" << baseline.size() << " events match exactly.\n";
  else
    std::cout << nMismatch << " mismatching/missing events out of " << baseline.size() << ".\n";

  std::cout << "\nSequential: " << tSeq << " s\n";
  std::cout << "Parallel:   " << tPar << " s (" << nThreads << " threads)\n";
  std::cout << "Speedup:    " << (tPar > 0.0 ? tSeq / tPar : 0.0) << "x\n";

  return ok ? 0 : 1;
}
