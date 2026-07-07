#ifndef SAMPLE_CONFIG_H
#define SAMPLE_CONFIG_H

// ---------------------------------------------------------------------------
// sample_config.h
//   Sample selection via a --sample=<name> CLI flag: which ntuple directory
//   to read, which energy/process label to draw on plots, and which output
//   directory to save PDFs to. Keeps the main VBF H->Invisible ntuple as the
//   default when no flag is given.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <thread>

namespace MyUtl {

  struct SampleConfig {
    std::string ntupleDir;
    std::string energyLabel;
    std::string outputDir;
    std::string sampleName;  // "" for the local/default run; e.g. "vbf" when --sample= is given
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

  // Default (no --sample flag): local VBF ntuple, VBF label, ../figs output.
  inline SampleConfig resolveSample(int argc, char** argv) {
    static const std::map<std::string, SampleConfig> registry = {
      {"vbf",   {"/data/mcardiff/exotic_superntuples/highstats_vbf/",
                 "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "vbf",   "vbf"}},
      {"zjets", {"/data/mcardiff/exotic_superntuples/zjets/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Z+jets",               "zjets", "zjets"}},
      {"dijet", {"/data/mcardiff/exotic_superntuples/dijet/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Dijet",                "dijet", "dijet"}},
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

    return {"/Users/mcard/project/ntuple-hgtd/", "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "../figs", ""};
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
  inline std::string histFilePath(const std::string& baseName) {
    if (!SAMPLE_NAME.empty())
      return OUTPUT_DIR + "/" + SAMPLE_NAME + "_" + baseName;
    return OUTPUT_DIR + "/hists/" + baseName;
  }

  inline std::string plotFilePath(const std::string& subdir, const std::string& baseName) {
    const std::string name = SAMPLE_NAME.empty() ? baseName : SAMPLE_NAME + "_" + baseName;
    return OUTPUT_DIR + "/" + subdir + "/" + name;
  }

}

#endif
