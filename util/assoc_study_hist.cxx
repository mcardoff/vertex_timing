// ---------------------------------------------------------------------------
// assoc_study_hist.cxx
//   Histogramming stage of the track-to-vertex z-association study: runs the
//   clustering event loop ONCE and, for every event, re-selects the clustering
//   track list under each rule in assocRules() and processes the event
//   separately for each. Writes every histogram to a single ROOT file for
//   assoc_study_plot to turn into the report.
//
//   See util/assoc_study_common.h for what the study asks and why it is
//   structured as one loop over many rules rather than many runs.
//
//   Usage:
//     ./assoc_study_hist [--sample=<name>] [--threads=N] [--max-events=N]
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "histogram_io.h"
#include "plotting_utilities.h"
#include "sample_config.h"

#include "assoc_study_common.h"

#include <ROOT/TTreeProcessorMT.hxx>
#include <ROOT/TThreadedObject.hxx>

#include <atomic>
#include <mutex>

using namespace MyUtl;

// ---------------------------------------------------------------------------
// ThreadState — everything one worker thread accumulates.
//   `maps` is parallel to assocRules(): maps[i] is the analysis map for rule i.
//   The track-accounting vectors are parallel to it too, and give the
//   resolution numbers their context: a rule that improves σ(Δt) by throwing
//   away most of its tracks is a different result from one that improves it at
//   constant statistics, and only these counters distinguish the two.
// ---------------------------------------------------------------------------
struct ThreadState {
  std::vector<std::map<Score, AnalysisObj>> maps;
  std::vector<Long64_t> nKept, nKeptHS, nKeptPU;   // tracks the rule associates
  Long64_t nAccepted   = 0;   // events surviving the event selection
  Long64_t nFwdAvail   = 0;   // forward tracks passing kinematics, any z
  Long64_t nFwdAvailHS = 0;   // ... of which truth-HS; recall denominator

  ThreadState() {
    const auto& rules = assocRules();
    maps.reserve(rules.size());
    for (const auto& r : rules) maps.push_back(buildAssocMap(r));
    nKept  .assign(rules.size(), 0);
    nKeptHS.assign(rules.size(), 0);
    nKeptPU.assign(rules.size(), 0);
  }

  void mergeFrom(ThreadState& o) {
    for (size_t i = 0; i < maps.size(); ++i) {
      for (auto& [score, analysis] : maps[i])
        analysis.mergeFrom(o.maps[i].at(score));
      nKept  [i] += o.nKept  [i];
      nKeptHS[i] += o.nKeptHS[i];
      nKeptPU[i] += o.nKeptPU[i];
    }
    nAccepted   += o.nAccepted;
    nFwdAvail   += o.nFwdAvail;
    nFwdAvailHS += o.nFwdAvailHS;
  }
};

int main(int argc, char* argv[]) {
  // Detach histograms from gDirectory before ANY histogram is constructed --
  // every worker thread builds an identically-named set, and without this they
  // race on gDirectory's global registry. Same requirement as clustering_hist.
  TH1::AddDirectory(kFALSE);

  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);
  MyUtl::ENERGY_LABEL    = sample.energyLabel;
  MyUtl::OUTPUT_DIR      = sample.outputDir;
  MyUtl::SAMPLE_NAME     = sample.sampleName;
  MyUtl::OVERLAP_REMOVAL = sample.overlapRemoval;
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  if (MyUtl::SAMPLE_NAME.empty())
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/hists");

  const unsigned nThreads  = MyUtl::resolveThreads(argc, argv);
  const Long64_t maxEvents = MyUtl::resolveMaxEvents(argc, argv);
  const auto&    rules     = assocRules();

  std::cout << "[assoc-study] comparing " << rules.size()
            << " track-to-vertex association rules in one pass:\n";
  for (const auto& r : rules)
    std::cout << "    " << std::left << std::setw(10) << r.tag
              << ruleAscii(r) << '\n';
  std::cout << "[assoc-study] counting scan pinned at z0 significance < "
            << COUNTING_NSIGMA
            << " for every rule, so plot x-axes are common\n";
  std::cout << "[assoc-study] scores: ";
  for (const auto& s : assocScores()) std::cout << s.toStringShort() << ' ';
  std::cout << '\n';

  TChain chain("ntuple");
  setupChain(chain, sample.ntupleDir.c_str());
  ROOT::EnableImplicitMT(nThreads);
  gErrorIgnoreLevel = kFatal;

  // Per-thread state registry. The mutex serializes the whole lazy
  // construction, not just the push_back: AnalysisObj's constructor calls
  // SetFillColorAlpha(), which lazily registers a TColor into ROOT's global
  // color table, and that registration is not thread-safe (see the same
  // pattern and its crash history in clustering_hist.cxx).
  std::mutex stateMutex;
  std::vector<std::unique_ptr<ThreadState>> stateRegistry;

  std::atomic<Long64_t> progressCounter{0};

  const Long64_t N_EVENT       = chain.GetEntries();
  const Long64_t progressDenom = (maxEvents > 0) ? std::min(maxEvents, N_EVENT) : N_EVENT;
  if (maxEvents > 0)
    std::cout << "Restricting to first " << progressDenom << " events (--max-events)\n";
  std::cout << "Starting Event Loop over " << N_EVENT << " events\n";

  ROOT::TTreeProcessorMT proc(chain, nThreads);
  proc.Process([&](TTreeReader& reader) {
    BranchPointerWrapper branch(reader);

    thread_local ThreadState* st = nullptr;
    if (!st) {
      std::lock_guard<std::mutex> lock(stateMutex);
      stateRegistry.push_back(std::make_unique<ThreadState>());
      st = stateRegistry.back().get();
    }

    while (reader.Next()) {
      Long64_t n = ++progressCounter;
      if (maxEvents > 0 && n > maxEvents) break;
      if (n % 5000 == 0)
        std::cout << "Progress: " << n << "/" << progressDenom << "\r" << std::flush;

      // Process the SAME event once per rule. Event selection (step A of
      // processEventData) happens before any track is looked at and does not
      // depend on the track list, so every rule accepts or rejects the same
      // events -- which is what makes this a paired comparison and lets the
      // accounting below be done once for all rules.
      int code = 0;
      for (size_t i = 0; i < rules.size(); ++i)
        code = processEventData(&branch, false, true, st->maps[i], &rules[i]).code;

      if (code < 0) continue;
      ++st->nAccepted;

      // Track accounting, on accepted events only. Deliberately re-derived
      // here rather than returned from processEventData: it is study-specific
      // bookkeeping, and keeping it local means EventResult -- shared by every
      // other executable -- does not grow fields nothing else reads.
      for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
        if (!passTrackKinematics(trk, &branch, MIN_TRACK_PT, MAX_TRACK_PT)) continue;
        const bool isHS = (branch.trackToTruthvtx[trk] == 0);
        ++st->nFwdAvail;
        if (isHS) ++st->nFwdAvailHS;
        for (size_t i = 0; i < rules.size(); ++i) {
          if (!passTrackVertexAssociation((int)trk, 0, &branch, rules[i])) continue;
          ++st->nKept[i];
          (isHS ? st->nKeptHS[i] : st->nKeptPU[i])++;
        }
      }
    }
  });
  std::cout << "\nFINISHED PROCESSING\n";

  // --- Merge per-thread state ---
  ThreadState merged;
  for (auto& s : stateRegistry) merged.mergeFrom(*s);

  const Long64_t nEventsProcessed = progressCounter.load();
  std::cout << "\nEvents processed: " << nEventsProcessed
            << "   accepted by selection: " << merged.nAccepted << '\n';

  // --- Save ---
  const std::string histPath = MyUtl::histFilePath(ASSOC_HIST_BASENAME);
  MyUtl::HistWriter writer(histPath);
  for (auto& m : merged.maps)
    for (auto& [score, analysis] : m) analysis.saveTo(writer);

  for (size_t i = 0; i < rules.size(); ++i) {
    writer.WriteScalar(assocScalarKey(rules[i], "n_kept"),    merged.nKept  [i]);
    writer.WriteScalar(assocScalarKey(rules[i], "n_kept_hs"), merged.nKeptHS[i]);
    writer.WriteScalar(assocScalarKey(rules[i], "n_kept_pu"), merged.nKeptPU[i]);
  }
  writer.WriteScalar("meta_n_accepted",     merged.nAccepted);
  writer.WriteScalar("meta_n_fwd_avail",    merged.nFwdAvail);
  writer.WriteScalar("meta_n_fwd_avail_hs", merged.nFwdAvailHS);
  writer.WriteRunMeta(MyUtl::ENERGY_LABEL, nEventsProcessed);
  writer.Close();
  std::cout << "Wrote histograms to " << histPath << '\n';

  // --- Track accounting, printed here as well as saved, so a run's log alone
  //     says what each rule did to the track list ---
  std::cout << "\n=== TRACK ACCOUNTING (forward, pT/quality-selected, accepted events) ===\n";
  std::cout << "  available: " << merged.nFwdAvail << " tracks, of which "
            << merged.nFwdAvailHS << " truth-HS\n\n";
  std::cout << std::left << std::setw(16) << "rule"
            << std::right << std::setw(12) << "kept"
            << std::setw(12) << "HS kept"
            << std::setw(12) << "PU kept"
            << std::setw(10) << "purity"
            << std::setw(10) << "recall" << '\n';
  std::cout << std::string(72, '-') << '\n';
  for (size_t i = 0; i < rules.size(); ++i) {
    double purity = merged.nKept[i]  > 0 ? 100.0 * merged.nKeptHS[i] / merged.nKept[i]   : 0.0;
    double recall = merged.nFwdAvailHS > 0 ? 100.0 * merged.nKeptHS[i] / merged.nFwdAvailHS : 0.0;
    std::cout << std::left << std::setw(16) << ruleAscii(rules[i])
              << std::right << std::setw(12) << merged.nKept  [i]
              << std::setw(12) << merged.nKeptHS[i]
              << std::setw(12) << merged.nKeptPU[i]
              << std::fixed << std::setprecision(2)
              << std::setw(9) << purity << "%"
              << std::setw(9) << recall << "%" << '\n';
  }
  std::cout << std::string(72, '-') << '\n';

  return 0;
}
