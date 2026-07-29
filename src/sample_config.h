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

#include "clustering_constants.h"  // VBS_JET_D_ETA (runtime-settable)

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
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
  inline SampleConfig resolveSample(int argc, char** argv) {
    static const std::map<std::string, SampleConfig> registry = {
      {"vbf",   {"/data/mcardiff/exotic_superntuples/highstats_vbf/",
                 "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "vbf",   "vbf",   false}},
      {"zjets", {"/data/mcardiff/exotic_superntuples/zjets/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Z+jets",               "zjets", "zjets", true}},
      {"dijet", {"/data/mcardiff/exotic_superntuples/dijet/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Dijet",                "dijet", "dijet", false}},
    };

    const std::string prefix = "--sample=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;
      std::string name = arg.substr(prefix.size());
      auto it = registry.find(name);
      if (it == registry.end()) {
        std::cerr << "Unknown --sample value '" << name
                  << "'. Valid options: vbf, zjets, dijet.\n";
        std::exit(1);
      }
      return it->second;
    }

    return {"/Users/mcard/project/ntuple-hgtd/", "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "../figs", "", false};
  }

  // ---------------------------------------------------------------------------
  // resolveSelection
  //   Optional selection override via --vbs-deta=<x>, which sets the VBS
  //   candidate-pair |Deta| requirement (default VBS_JET_D_ETA_DEFAULT = 3.0;
  //   pass 0 to disable the cut entirely while keeping the >=1 forward jet
  //   requirement, i.e. "still a forward-jet event, but not forced into a VBS
  //   topology"). See the comment on VBS_JET_D_ETA in clustering_constants.h for
  //   why this needs to be measurable rather than assumed.
  //
  //   Also sets SELECTION_TAG whenever the value differs from the default, so
  //   the loosened run writes <sample>_deta0p0_<name>.root instead of colliding
  //   with the standard <sample>_<name>.root. Absent flag => no tag => every
  //   existing path unchanged.
  //
  //   Call once in main() alongside resolveSample(). Executables that never see
  //   the flag are unaffected.
  // ---------------------------------------------------------------------------
  inline void resolveSelection(int argc, char** argv) {
    const std::string prefix = "--vbs-deta=";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.rfind(prefix, 0) != 0) continue;

      double v = 0.0;
      try {
        v = std::stod(arg.substr(prefix.size()));
      } catch (const std::exception&) {
        std::cerr << "Bad --vbs-deta value '" << arg.substr(prefix.size())
                  << "'; expected a number (e.g. --vbs-deta=0).\n";
        std::exit(1);
      }
      if (v < 0.0) {
        std::cerr << "--vbs-deta must be >= 0 (got " << v << ").\n";
        std::exit(1);
      }

      VBS_JET_D_ETA = v;
      if (std::abs(v - VBS_JET_D_ETA_DEFAULT) > 1e-9) {
        std::ostringstream os;
        os << "deta" << std::fixed << std::setprecision(1) << v;
        std::string tag = os.str();
        std::replace(tag.begin(), tag.end(), '.', 'p');
        SELECTION_TAG = tag;
      }
      break;
    }
    std::cout << "[selection] VBS |Deta| >= " << VBS_JET_D_ETA
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

  inline std::string histFilePath(const std::string& baseName) {
    const std::string base = selectionTagged(baseName);
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
