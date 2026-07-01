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
  };

  // Mutable globals read by plotting code (event-label text, PDF output root).
  // Set once in main() right after resolveSample(), before any plots are made.
  inline std::string ENERGY_LABEL = "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.";
  inline std::string OUTPUT_DIR   = "../figs";

  // Default (no --sample flag): local VBF ntuple, VBF label, ../figs output.
  inline SampleConfig resolveSample(int argc, char** argv) {
    static const std::map<std::string, SampleConfig> registry = {
      {"vbf",   {"/data/mcardiff/exotic_superntuples/highstats_vbf/",
                 "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "../vbf"}},
      {"zjets", {"/data/mcardiff/exotic_superntuples/zjets/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Z+jets",               "../zjets"}},
      {"dijet", {"/data/mcardiff/exotic_superntuples/dijet/",
                 "#sqrt{s} = 14 TeV, HL-LHC, Dijet",                "../dijet"}},
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

    return {"../../ntuple-hgtd/", "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.", "../figs"};
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

}

#endif
