#include <fstream>
#include <sstream>
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
  MyUtl::resolveTimingInputs(argc, argv);  // --clean-times / --assoc-nsigma: the file's tag
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  for (const char* sub : {"comparisons", "inclusive", "fullplots", "diagnostics", "fixes", "talk"})
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
    reader.CheckSelection(MyUtl::VBS_JET_D_ETA, MyUtl::VBS_JET_MJJ);
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
  // The talk's t0s overlaid (slides 8/9): WAVeS against Athena (user, 2026-11-11).
  const std::vector<AnalysisObj*> TALK_SET = {&mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::WAVES)};
  shapeComparisonPlot(MyUtl::plotFilePath("talk", "talk_t0_shape_logscale.pdf").c_str(),
    true, -400, 400, canvas, TALK_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlot(MyUtl::plotFilePath("talk", "talk_t0_shape_linscale.pdf").c_str(),
    false, -200, 200, canvas, TALK_SET, resoGetter, nullptr, diffLabel);
  shapeComparisonPlotPair(
    MyUtl::plotFilePath("inclusive", "shape_comparison.pdf").c_str(),
    -200, 200, -400, 400, canvas, SHAPE_SET, resoGetter, nullptr, diffLabel);

  // --- Failure-mode fixes: inclusive core fraction and resolution per row ---
  // The numbers the talk's summary figure reads (python/gp_failure_fixes.py):
  // core fraction = |t0 - t_HS| < 60 ps over all selected events (no t0 = a
  // failure), resolution = FIT_PEAKTAIL on the inclusive t0 - t_HS.
  if (mapHGTD.count(Score::WFIX_RETIME)) {
    const std::string jpath = MyUtl::plotFilePath("fixes", "failure_fixes.json");
    std::ofstream js(jpath);
    js << "{\n  \"sample\": \"" << (MyUtl::SAMPLE_NAME.empty() ? "local" : MyUtl::SAMPLE_NAME)
       << "\",\n  \"tag\": \"" << MyUtl::SELECTION_TAG << "\",\n  \"rows\": [\n";
    std::printf("\n%-14s %-40s %9s %9s %11s %11s\n", "row", "", "core [%]", "+-", "sigPeak[ps]", "sigCore[ps]");
    bool first = true;
    for (const Score* sc : {&Score::HGTD, &Score::TRKPTZ, &Score::WAVES, &Score::WFIX_RETIME, &Score::WFIX_MISTIMED,
                            &Score::WFIX_SELECT, &Score::WFIX_MISCL, &Score::WFIX_TEFF, &Score::WFIX_C2,
                            &Score::WFIX_C3, &Score::WFIX_C4, &Score::WFIX_C5, &Score::WFIX_C6}) {
      AnalysisObj& a = mapHGTD.at(*sc);
      PlotObj* po = a.dataObjects.at("hs_track").get();
      const int nb = po->effTotal->GetNbinsX();
      const double tot = po->effTotal->Integral(0, nb + 1), pas = po->effPass->Integral(0, nb + 1);
      const double eff = tot > 0 ? pas / tot : 0.0, err = tot > 0 ? std::sqrt(eff * (1 - eff) / tot) : 0.0;
      std::unique_ptr<TH1D> h((TH1D*)a.inclusiveResoSig->Clone(TString::Format("fixsum_%d", sc->id)));
      h->Add(a.inclusiveResoMix.get()); h->Add(a.inclusiveResoBkg.get());
      double sp = -1, spe = -1, scr = -1, fsh = -1, ntail = 0;
      if (h->GetEntries() > 50) {
        TF1* f = makePeakTailFit(h.get(), false);
        sp = f->GetParameter("SigmaPeak"); spe = f->GetParError(f->GetParNumber("SigmaPeak"));
        fsh = f->GetParameter("FracSh");
        scr = sp / std::sqrt(1.0 - fsh + fsh * std::pow(f->GetParameter("RatioSh"), 2));
        ntail = f->GetParameter("NormTail") * f->GetParameter("SigmaTail") * std::sqrt(2 * M_PI) / h->GetBinWidth(1);
        delete f;
      }
      std::printf("%-14s %-40s %9.2f %9.2f %11.2f %11.2f\n", sc->toStringShort(), sc->toString(), 100 * eff, 100 * err, sp, scr);
      js << (first ? "" : ",\n") << "    {\"short\": \"" << sc->toStringShort() << "\", \"name\": \"" << sc->toString()
         << "\", \"fixes\": " << sc->fixes << ", \"n\": " << tot << ", \"pass\": " << pas
         << ", \"core\": " << eff << ", \"core_err\": " << err << ", \"sigma_peak\": " << sp
         << ", \"sigma_peak_err\": " << spe << ", \"sigma_core\": " << scr << ", \"frac_shoulder\": " << fsh
         << ", \"n_tail_fit\": " << ntail << ", \"n_reso\": " << h->GetEntries();
      // Per bin of forward HS tracks with a time, the underflow (n <= 2, below
      // HS_TRACK_MIN, where most failures are) included as the first bin; the
      // last bin is the fold (>= FOLD_HS_TRACK). sigma_peak: FIT_PEAKTAIL on the
      // bin's t0 - t_HS, as the main program's slices (-1 below 100 entries).
      js << ", \"nhs\": {\"x\": [";
      const int nx = po->hist->GetNbinsX();
      std::ostringstream xs, ps, ts, ss, es;
      for (int ib = 0; ib <= nx; ++ib) {
        const char* sep = ib ? ", " : "";
        xs << sep << (ib == 0 ? -1.0 : po->hist->GetXaxis()->GetBinCenter(ib));
        ps << sep << po->effPass->GetBinContent(ib);
        ts << sep << po->effTotal->GetBinContent(ib);
        std::unique_ptr<TH1D> sl(po->hist->ProjectionY(TString::Format("fixsl_%d_%d", sc->id, ib), ib, ib));
        double v = -1, e = -1;
        if (sl->GetEntries() >= 100) {
          TF1* f = makePeakTailFit(sl.get(), false);
          v = f->GetParameter("SigmaPeak"); e = f->GetParError(f->GetParNumber("SigmaPeak"));
          delete f;
        }
        ss << sep << v; es << sep << e;
      }
      js << xs.str() << "], \"pass\": [" << ps.str() << "], \"total\": [" << ts.str()
         << "], \"sigma_peak\": [" << ss.str() << "], \"sigma_peak_err\": [" << es.str() << "]}}";
      first = false;
    }
    js << "\n  ]\n}\n";
    std::cout << "wrote " << jpath << '\n';
  }

  std::cout << "FINISHED PLOT PRINTING\n";

  return 0;
}
