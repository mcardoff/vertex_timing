#ifndef SAMPLE_CONFIG_H
#define SAMPLE_CONFIG_H

// ---------------------------------------------------------------------------
// sample_config.h
//   Sample selection via a --sample=<name> CLI flag: which ntuple directory
//   to read, which energy/process label to draw on plots, and which output
//   directory to save PDFs to. Keeps the main VBF H->Invisible ntuple as the
//   default when no flag is given.
// ---------------------------------------------------------------------------

#include <Rtypes.h>  // Long64_t

#include "clustering_constants.h"  // VBS_JET_D_ETA, VBS_JET_MJJ (runtime-settable)

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <cstdio>
#include <unistd.h>
#include <chrono>
#include <thread>

namespace MyUtl {

  struct SampleConfig {
    std::string ntupleDir;
    std::string energyLabel;
    std::string outputDir;
    std::string sampleName;  // "" for the local/default run; e.g. "vbf" when --sample= is given
    bool        overlapRemoval = false;  // lepton–jet overlap removal (Z+jets only; see OVERLAP_REMOVAL)
  };

  // Mutable globals read by plotting code (event-label text, PDF output root).
  // Set once in main() right after resolveSample(), before any plots are made.
  //
  // Default (no --sample, local dev) keeps the "../figs" convention: you
  // `cd build && ./clustering_hist`, and output lands in the repo-root
  // figs/, one level up. This one is never run on condor, so the "assumes
  // we're inside build/" convention is harmless here.
  //
  // The grid-sample entries (vbf/zjets/dijet) below are cwd-relative instead
  // ("vbf" not "../vbf") -- deliberately, since those are the ones actually
  // submitted to condor. The old "../<name>" convention only resolved
  // correctly there because run_analysis.sh faked being inside a build/
  // subdirectory (mkdir it, mv the binary in, cd into it) purely so "../"
  // would land back at the job's scratch root. If that step were ever
  // skipped -- e.g. running the transferred binary directly by hand on an
  // execute node -- "../<name>" would climb OUT of the scratch directory
  // into its parent, which may not be writable (or may belong to a
  // different job entirely). A cwd-relative path can't escape the directory
  // the job is actually confined to. See run_analysis.sh, which no longer
  // needs the build/ indirection as a result.
  inline std::string ENERGY_LABEL = "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.";
  inline std::string OUTPUT_DIR   = "../figs";

  // Sample name for the current run ("" for the local/default run, otherwise
  // e.g. "vbf"). Mirrors ENERGY_LABEL/OUTPUT_DIR: set once in main() right
  // after resolveSample(). Read by histFilePath()/plotFilePath() below to
  // decide whether output file names get a sample prefix.
  inline std::string SAMPLE_NAME  = "";

  // Non-default selection marker, e.g. "deta0p0". Empty for the standard
  // selection, so every existing output path is byte-identical to before.
  // Set by resolveSelection() and woven into histFilePath()/plotFilePath() so a
  // loosened run cannot silently overwrite the standard one's outputs.
  inline std::string SELECTION_TAG = "";

  // Default (no --sample flag): local VBF ntuple, VBF label, ../figs output.
  //
  // PAIRING the mu=0 samples with a mu=200 counterpart:
  //   vbf_mu0     600026 r15697  pairs cleanly with "vbf" (600026 r15700).
  //   zeejets_mu0 601189 Zee     pairs with "zeejets", NOT with "zjets".
  //                              The zjets directory holds TWO DSIDs side by
  //                              side -- 601189 (Zee) and 601190 (Zmumu) -- and
  //                              setupChain descends one level, so --sample=zjets
  //                              reads BOTH channels. Every earlier zjets number
  //                              in this study is therefore Zee+Zmumu combined.
  //                              Since only Zee exists at mu=0, comparing it
  //                              against the combined sample would confound
  //                              pileup with channel mix, so "zeejets" below
  //                              points at the 601189 subdirectory alone.
  //   ttbar_mu0   601229 r15697  has no mu=200 counterpart yet; the ttbar
  //                              directory exists but is empty.
  inline SampleConfig resolveSample(int argc, char** argv) {
    static const std::map<std::string, SampleConfig> registry = {
      {"vbf",   {"/data/mcardiff/exotic_superntuples/highstats_vbf/",
                 "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "vbf",   "vbf",   false}},
      {"zjets", {"/data/mcardiff/exotic_superntuples/zjets/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Z+jets",               "zjets", "zjets", true}},
      {"dijet", {"/data/mcardiff/exotic_superntuples/dijet/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Dijet",                "dijet", "dijet", false}},

      // Zee ONLY at mu=200, i.e. the channel-matched partner for zeejets_mu0.
      // "zjets" above deliberately keeps reading both DSIDs, so every previously
      // published zjets number stays reproducible; this entry exists purely so
      // the mu0-vs-mu200 comparison is like-for-like.
      {"zeejets", {"/data/mcardiff/exotic_superntuples/zjets/"
                   "user.mcardiff.mc21_14TeV.601189.PhPy8EG_AZNLO_Zee.SuperNtuple."
                   "e8481_s4290_r15700.20260712_Output/",
                   "#sqrt{s} = 14 TeV, HL-LHC, Z(ee)+jets", "zeejets", "zeejets", true}},

      // ---- mu = 0 (PhaseIINoPileUp, r15697) --------------------------------
      // These exist to measure the INTRINSIC timing floor: with no pileup, the
      // core fraction a sample can reach is set by per-track time resolution
      // alone. That is the number which turns "the residual needs new detector
      // information" from an elimination argument into a measurement.
      //
      // Expect the SELECTION EFFICIENCY to collapse on zeejets_mu0 and
      // ttbar_mu0, and that collapse is itself the result: passJetPtCut needs a
      // forward jet (2.38 < |eta| < 4.0, pT > 30 GeV), and neither Z+jets nor
      // ttbar produces one on its own -- at mu=200 that requirement is largely
      // satisfied by PILEUP jets. vbf_mu0 should pass normally, since forward
      // tagging jets are the VBF signature, which is why it is the one that can
      // support a timing measurement.
      {"vbf_mu0",     {"/data/mcardiff/exotic_superntuples/vbf_mu0/",
                       "#sqrt{s} = 14 TeV, #mu = 0, VBF H#rightarrowinv.",
                       "vbf_mu0", "vbf_mu0", false}},
      // Z->ee, so the OS-SF pair in passLeptonSelection is found among electrons.
      // Compare against "zeejets" (same DSID at mu=200), NOT "zjets", which is
      // Zee+Zmumu combined -- see the pairing note in the header comment.
      {"zeejets_mu0", {"/data/mcardiff/exotic_superntuples/zeejets_mu0/",
                       "#sqrt{s} = 14 TeV, #mu = 0, Z(ee)+jets",
                       "zeejets_mu0", "zeejets_mu0", true}},
      // ttbar SingleLep has exactly ONE lepton, so overlapRemoval is left false:
      // it currently also switches on passLeptonSelection, whose OS-SF PAIR
      // requirement would reject every event. Separating those two behaviours is
      // a prerequisite for ever running ttbar with lepton-jet overlap removal.
      {"ttbar_mu0",   {"/data/mcardiff/exotic_superntuples/ttbar_mu0/",
                       "#sqrt{s} = 14 TeV, #mu = 0, t#bar{t} (1 lep)",
                       "ttbar_mu0", "ttbar_mu0", false}},
    };

    const std::string prefix = "--sample=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;
      std::string name = arg.substr(prefix.size());
      auto it = registry.find(name);
      if (it == registry.end()) {
        std::cerr << "Unknown --sample value '" << name
                  << "'. Valid options: vbf, zjets, zeejets, dijet, vbf_mu0, zeejets_mu0, ttbar_mu0.\n";
        std::exit(1);
      }
      return it->second;
    }

    return {"/Users/mcard/project/ntuple-hgtd/", "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "../figs", "", false};
  }

  // ---------------------------------------------------------------------------
  // resolveSelection
  //   Optional selection overrides for the VBS candidate-pair cuts:
  //     --vbs-deta=<x>  VBS_JET_D_ETA (default 3.0); pass 0 to disable while
  //                     keeping the >=1 forward jet requirement, i.e. "still a
  //                     forward-jet event, but not forced into a VBS topology".
  //     --vbs-mjj=<x>   VBS_JET_MJJ (default 0, a no-op).
  //   Independent knobs, not alternatives -- passJetPtCut requires both. See
  //   the comments on VBS_JET_D_ETA / VBS_JET_MJJ in clustering_constants.h
  //   for why each needs to be measurable/settable rather than assumed.
  //
  //   Each sets a SELECTION_TAG component whenever its value differs from its
  //   default, joined with '_' if both are non-default, so a loosened run
  //   writes e.g. <sample>_deta0p0_mjj200p0_<name>.root instead of colliding
  //   with the standard <sample>_<name>.root. Neither flag given => no tag =>
  //   every existing path unchanged.
  //
  //   Call once in main() alongside resolveSample(). Executables that never
  //   see either flag are unaffected.
  // ---------------------------------------------------------------------------
  inline void resolveSelection(int argc, char** argv) {
    // Parses --<name>=<x> if present; returns false (default, val untouched)
    // if the flag never appears. Shared by both knobs below so a bad/missing
    // value is reported the same way for each.
    auto parseFlag = [&](const char* name, double& val) -> bool {
      const std::string prefix = std::string("--") + name + "=";
      for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind(prefix, 0) != 0) continue;
        try {
          val = std::stod(arg.substr(prefix.size()));
        } catch (const std::exception&) {
          std::cerr << "Bad --" << name << " value '" << arg.substr(prefix.size())
                    << "'; expected a number (e.g. --" << name << "=0).\n";
          std::exit(1);
        }
        if (val < 0.0) {
          std::cerr << "--" << name << " must be >= 0 (got " << val << ").\n";
          std::exit(1);
        }
        return true;
      }
      return false;
    };

    auto addTag = [&](const char* label, double v) {
      std::ostringstream os;
      os << label << std::fixed << std::setprecision(1) << v;
      std::string tag = os.str();
      std::replace(tag.begin(), tag.end(), '.', 'p');
      SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
    };

    double deta = VBS_JET_D_ETA, mjj = VBS_JET_MJJ;
    if (parseFlag("vbs-deta", deta)) {
      VBS_JET_D_ETA = deta;
      if (std::abs(deta - VBS_JET_D_ETA_DEFAULT) > 1e-9) addTag("deta", deta);
    }
    if (parseFlag("vbs-mjj", mjj)) {
      VBS_JET_MJJ = mjj;
      if (std::abs(mjj - VBS_JET_MJJ_DEFAULT) > 1e-9) addTag("mjj", mjj);
    }

    std::cout << "[selection] VBS |Deta| >= " << VBS_JET_D_ETA
              << ", m_jj >= " << VBS_JET_MJJ << " GeV"
              << (SELECTION_TAG.empty() ? "  (default)"
                                        : "  (tag: " + SELECTION_TAG + ")")
              << '\n';
  }

  // ---------------------------------------------------------------------------
  // resolveThreads
  //   Thread count for TTreeProcessorMT via a --threads=<N> CLI flag. Default
  //   is capped (not raw hardware_concurrency()) so it doesn't silently drift
  //   away from whatever request_cpus a condor job was submitted with --
  //   threads beyond what condor's cgroup granted just get throttled, which
  //   would eat the whole parallelization benefit without any visible error.
  // ---------------------------------------------------------------------------
  inline unsigned resolveThreads(int argc, char** argv) {
    unsigned nThreads = std::min(std::thread::hardware_concurrency(), 8u);

    const std::string prefix = "--threads=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;
      nThreads = static_cast<unsigned>(std::stoul(arg.substr(prefix.size())));
      break;
    }

    return nThreads > 0 ? nThreads : 1u;
  }

  // ---------------------------------------------------------------------------
  // Shard: which slice of a sample's files this process handles.
  //   index in [0, count), count >= 1. The default {0, 1} means "all files",
  //   so any executable or invocation that never mentions sharding is
  //   unchanged.
  //
  //   Sharding exists because a single job pays the whole per-file startup
  //   cost serially (0.3-1.2 s per file on the AF's /data) AND runs the whole
  //   event loop. Splitting a sample N ways divides both. See CLAUDE.md's
  //   "Job cost" section for the measurements that motivated it.
  // ---------------------------------------------------------------------------
  struct Shard {
    int index = 0;
    int count = 1;
    // Whether --file-shard was actually passed. Deliberately separate from
    // count > 1: a submit file that templates the output name on the shard
    // must get a shard-tagged name for EVERY shard count it might be run at,
    // including N = 1. Tying the tag to count > 1 instead would make
    // --file-shard=0/1 write an untagged file that the .sub then fails to
    // transfer back -- a foot-gun that only fires at the least-tested value.
    bool given = false;

    // True when the rule actually partitions the file list. N = 1 selects
    // everything, so it does not.
    bool active() const { return count > 1; }

    // Filename-safe marker woven into the output path so shards of one sample
    // cannot overwrite each other.
    std::string tag() const {
      return given ? ".shard" + std::to_string(index) + "of" + std::to_string(count)
                   : std::string();
    }
  };

  // ---------------------------------------------------------------------------
  // resolveShard
  //   --file-shard=<i>/<N>: this process handles files i, i+N, i+2N, ... of the
  //   sample's sorted file list. Absent => {0, 1} => every file.
  // ---------------------------------------------------------------------------
  inline Shard resolveShard(int argc, char** argv) {
    const std::string prefix = "--file-shard=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;
      std::string v = arg.substr(prefix.size());
      auto slash = v.find('/');
      if (slash == std::string::npos) {
        std::cerr << "Bad --file-shard value '" << v << "'; expected <i>/<N>, e.g. 0/10.\n";
        std::exit(1);
      }
      Shard s;
      s.given = true;
      try {
        s.index = std::stoi(v.substr(0, slash));
        s.count = std::stoi(v.substr(slash + 1));
      } catch (const std::exception&) {
        std::cerr << "Bad --file-shard value '" << v << "'; expected <i>/<N>, e.g. 0/10.\n";
        std::exit(1);
      }
      if (s.count < 1 || s.index < 0 || s.index >= s.count) {
        std::cerr << "--file-shard out of range: got " << s.index << "/" << s.count
                  << "; need 0 <= i < N and N >= 1.\n";
        std::exit(1);
      }
      return s;
    }
    return {};
  }

  // Set once in main() alongside SAMPLE_NAME, read by histFilePath() so every
  // shard writes a distinct file. Mirrors SELECTION_TAG's role.
  inline Shard FILE_SHARD = {};

  // ---------------------------------------------------------------------------
  // resolveMaxEvents
  //   Optional event cap via a --max-events=<N> CLI flag, for quick local
  //   sanity checks (e.g. a diagnostic breakdown) without waiting on a full
  //   multi-million-event condor-scale run. Returns -1 (no cap; default) when
  //   the flag is absent, so every existing invocation is unaffected.
  // ---------------------------------------------------------------------------
  inline Long64_t resolveMaxEvents(int argc, char** argv) {
    const std::string prefix = "--max-events=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;
      return static_cast<Long64_t>(std::stoll(arg.substr(prefix.size())));
    }
    return -1;
  }

  // ---------------------------------------------------------------------------
  // resolveHistFile
  //   Path to the histogram ROOT file a *_plot executable reads, via an
  //   optional --hist-file=<path> CLI override (e.g. when a file transferred
  //   back from condor lands somewhere other than the default --sample=
  //   convention). Defaults to `defaultPath` (typically
  //   OUTPUT_DIR + "/hists/<name>.root", matching what the corresponding
  //   *_hist executable wrote).
  // ---------------------------------------------------------------------------
  inline std::string resolveHistFile(int argc, char** argv, const std::string& defaultPath) {
    const std::string prefix = "--hist-file=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;
      return arg.substr(prefix.size());
    }
    return defaultPath;
  }

  // ---------------------------------------------------------------------------
  // histFilePath / plotFilePath
  //   Sample-aware output path builders, read by both stages of the
  //   *_hist/*_plot split. When a --sample was given (SAMPLE_NAME non-empty),
  //   every output file name is prefixed with "<SAMPLE_NAME>_" so files stay
  //   distinguishable if copied out of their per-sample directory into a
  //   shared one. For histograms this also flattens away the old "hists/"
  //   subdirectory (a single sample-prefixed file directly under OUTPUT_DIR
  //   no longer needs it to stay unambiguous). Plots keep their existing
  //   subdirectory layout (comparisons/, inclusive/, rpt_plots/, ...) -- only
  //   the file name gets prefixed. For the local/default run (no --sample,
  //   SAMPLE_NAME == ""), both fall back to the original unprefixed
  //   conventions.
  // ---------------------------------------------------------------------------
  //   SELECTION_TAG (non-empty only for a non-default selection) is inserted
  //   ahead of the base name so a loosened run lands beside the standard one
  //   rather than overwriting it.
  inline std::string selectionTagged(const std::string& baseName) {
    return SELECTION_TAG.empty() ? baseName : SELECTION_TAG + "_" + baseName;
  }

  // True when stdout is a terminal. Progress output that overwrites its own
  // line with '\r' is right for a terminal and wrong for a condor log, where
  // it collapses the whole run into one unreadable line; callers branch on
  // this rather than emitting terminal control characters into a file.
  inline const bool STDOUT_IS_TTY = isatty(fileno(stdout)) != 0;

  // ---------------------------------------------------------------------------
  // PhaseTimer
  //   Prints "[t+SSS.S s] <label>" for each milestone of a long-running job,
  //   flushed immediately.
  //
  //   Both halves matter on condor. The timestamps make the pre-event-loop
  //   cost attributable instead of guessable -- the startup dead time on the
  //   AF has been measured at 5-30 minutes and varies ~4x run-to-run on the
  //   same sample, so it has to be measured per run, not assumed. And the
  //   flush matters because a redirected stdout is block-buffered: without it
  //   a `condor_tail` of a job that is genuinely working shows nothing at
  //   all, which is indistinguishable from a hang.
  // ---------------------------------------------------------------------------
  class PhaseTimer {
  public:
    PhaseTimer() : start_(std::chrono::steady_clock::now()) {}
    // Composed into one string and written with a single insertion: mark() is
    // called from TTreeProcessorMT worker threads as well as main, and
    // streaming several fragments would let two threads interleave mid-line.
    void mark(const std::string& label) const {
      double t = std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - start_).count();
      std::ostringstream os;
      os << "[t+" << std::fixed << std::setprecision(1) << std::setw(7)
         << t << " s] " << label << '\n';
      std::cout << os.str() << std::flush;
    }
    double elapsed() const {
      return std::chrono::duration<double>(
               std::chrono::steady_clock::now() - start_).count();
    }
  private:
    std::chrono::steady_clock::time_point start_;
  };

  //   FILE_SHARD (non-default only when --file-shard was given) inserts a
  //   ".shard<i>of<N>" marker before the extension, so the N shards of one
  //   sample land side by side instead of overwriting each other. Empty for
  //   an unsharded run, leaving every existing path byte-identical.
  inline std::string shardTagged(const std::string& baseName) {
    if (FILE_SHARD.tag().empty()) return baseName;
    auto dot = baseName.rfind('.');
    return dot == std::string::npos
             ? baseName + FILE_SHARD.tag()
             : baseName.substr(0, dot) + FILE_SHARD.tag() + baseName.substr(dot);
  }

  inline std::string histFilePath(const std::string& baseName) {
    const std::string base = shardTagged(selectionTagged(baseName));
    if (!SAMPLE_NAME.empty())
      return OUTPUT_DIR + "/" + SAMPLE_NAME + "_" + base;
    return OUTPUT_DIR + "/hists/" + base;
  }

  inline std::string plotFilePath(const std::string& subdir, const std::string& baseName) {
    const std::string base = selectionTagged(baseName);
    const std::string name = SAMPLE_NAME.empty() ? base : SAMPLE_NAME + "_" + base;
    return OUTPUT_DIR + "/" + subdir + "/" + name;
  }

}

#endif
