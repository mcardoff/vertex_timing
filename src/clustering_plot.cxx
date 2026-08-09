#include "clustering_constants.h"
#include "clustering_dt_common.h"
#include "sample_config.h"
#include "plotting_utilities.h"
#include "histogram_io.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// clustering_plot.cxx
//   Plotting stage of the clustering_dt split: reads the histogram file
//   written by clustering_hist.cxx (no ntuple access, no event loop) and
//   reproduces the same comparison/inclusive/shape/fullplots PDF output the
//   former monolithic clustering_dt.cxx produced -- see CLAUDE.md's "Main
//   Executables" section.
// ---------------------------------------------------------------------------

auto main(int argc, char** argv) -> int {
  SetAtlasStyle();  // ATLAS publication style (fonts, margins, ticks, …)

  // --- Sample selection (--sample=vbf|zjets|dijet; default: local VBF ntuple) ---
  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);  // --vbs-deta=<x>; sets SELECTION_TAG
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  for (const char* sub : {"comparisons", "inclusive", "fullplots", "diagnostics"})
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/" + sub);

  gErrorIgnoreLevel = kFatal;

  // --- Canvas ---
  TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);

  // --- Load histograms written by clustering_hist.cxx ---
  const std::string histFile = MyUtl::resolveHistFile(
    argc, argv, MyUtl::histFilePath("clustering_hist.root"));
  std::cout << "Reading histograms from " << histFile << '\n';

  std::map<Score, AnalysisObj> mapHGTD = buildAnalysisMap(Scenario::HGTD);
  {
    MyUtl::HistReader reader(histFile);
    for (auto& [score, analysis] : mapHGTD) analysis.loadFrom(reader);
  }

  auto allMaps = { &mapHGTD };

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

  inclusivePlot(MyUtl::plotFilePath("inclusive", "inclusivereso_logscale.pdf").c_str(),
		true,  false, -400, 400, canvas, RESO_SET, resoGetter, &FIT_TRPGAUS, diffLabel, false, INCL_RESO_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(MyUtl::plotFilePath("inclusive", "inclusivereso_linscale.pdf").c_str(),
		false, false, -200, 200, canvas, RESO_SET, resoGetter, &FIT_TRPGAUS, diffLabel, false, INCL_RESO_YMAX_LIN, 0.0);


  auto pullGetter = [](AnalysisObj* a) -> ResoTriple {
    return { a->inclusivePullSig.get(),
             a->inclusivePullMix.get(),
             a->inclusivePullBkg.get() }; };
  const char* pullLabel = "(t_{0}-t_{truth})/#sigma_{t}";
  inclusivePlot(MyUtl::plotFilePath("inclusive", "inclusivepull_logscale.pdf").c_str(),
		true,  false, -10, 10, canvas, RESO_SET,
		pullGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LOG, INCL_YMIN_LOG);
  inclusivePlot(MyUtl::plotFilePath("inclusive", "inclusivepull_linscale.pdf").c_str(),
		false, false, -10, 10, canvas, RESO_SET,
		pullGetter, &FIT_PULLGAUS, pullLabel, false, INCL_PULL_YMAX_LIN, 0.0);

  // --- Normalized shape comparison: TRKPTZ vs oracle scores ---
  const std::vector<AnalysisObj*> SHAPE_SET = {
    &mapHGTD.at(Score::TRKPTZ),
    &mapHGTD.at(Score::WAVES),
    &mapHGTD.at(Score::WAVES_MISAS)
  };
  shapeComparisonPlot(
    MyUtl::plotFilePath("inclusive", "shape_comparison_linscale.pdf").c_str(),
    false, -200, 200, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlot(
    MyUtl::plotFilePath("inclusive", "shape_comparison_logscale.pdf").c_str(),
    true, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlotPair(
    MyUtl::plotFilePath("inclusive", "shape_comparison.pdf").c_str(),
    -200, 200, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);

  std::cout << "FINISHED PLOT PRINTING\n";

  return 0;
}
