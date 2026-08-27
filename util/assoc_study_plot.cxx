// ---------------------------------------------------------------------------
// assoc_study_plot.cxx
//   Plotting/report stage of the track-to-vertex z-association study. Reads
//   the histogram file written by assoc_study_hist.cxx (no ntuple access, no
//   event loop) and produces:
//     * console tables of every metric, for both selection scores
//     * a machine-readable block (prefix "CSV|") for pasting into the report
//     * overlay PDFs under <OUTPUT_DIR>/assoc_study/
//
//   Metric definitions, all reusing the project's existing machinery so these
//   numbers are directly comparable with the rest of the analysis:
//
//     core sigma         ACTIVE_FIT_MODEL's Sigma1 -- the same fit and the
//                        same parameter PlotObj::plotPostProcessing writes
//                        into params->sigmaDist, so the inclusive number and
//                        the per-bin numbers are on one footing.
//     core fraction      passed / total from the efficiency histograms, with
//                        under- and overflow included (Integral(0, n+1)).
//                        Including them is not cosmetic: PlotObj's n-HS-track
//                        axis starts at 2.5, so every event with fewer than
//                        three forward HS tracks lands in the underflow bin --
//                        and per the project's failure decomposition that is
//                        where most failures live. Excluding it would quietly
//                        report a core fraction for the easy events only.
//     truncated RMS      RMS of Delta t inside |Delta t| < PASS_SIGMA. Carried
//                        as a fit-model-independent cross-check, since a
//                        single Gaussian fitted over the full +/-1000 ps range
//                        is pulled by the tail and can move for reasons that
//                        have nothing to do with the core.
//
//   Usage:
//     ./assoc_study_plot [--sample=<name>] [--hist-file=<path>]
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "plotting_utilities.h"
#include "histogram_io.h"
#include "sample_config.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

#include "assoc_study_common.h"

#include <iomanip>
#include <iostream>
#include <sstream>

using namespace MyUtl;

namespace {

// One rule x one score, fully reduced to the numbers the report quotes.
struct Summary {
  double sigma = 0.0, sigmaErr = 0.0;      // inclusive core sigma [ps]
  double coreFrac = 0.0, coreFracErr = 0.0; // inclusive core fraction [%]
  double truncRMS = 0.0;                    // RMS inside |dt| < PASS_SIGMA [ps]
  double sigmaLow = 0.0;                    // core sigma, nHS <= 5 events
  double nEvents = 0.0;
};

// Sum the three purity layers into the one inclusive Delta t distribution.
// Caller owns the result.
TH1D* stackTotal(TH1D* sig, TH1D* mix, TH1D* bkg, const char* name) {
  TH1D* t = (TH1D*)bkg->Clone(name);
  t->Add(mix);
  t->Add(sig);
  return t;
}

Summary summarize(AnalysisObj& a) {
  Summary s;

  std::unique_ptr<TH1D> total(stackTotal(a.inclusiveResoSig.get(),
                                         a.inclusiveResoMix.get(),
                                         a.inclusiveResoBkg.get(),
                                         TString::Format("tot_%p", (void*)&a)));
  if (total->GetEntries() > 0) {
    std::unique_ptr<TF1> fit(ACTIVE_FIT_MODEL.makeFit(total.get()));
    s.sigma    = fit->GetParameter("Sigma1");
    s.sigmaErr = fit->GetParError("Sigma1");

    std::unique_ptr<TH1D> trunc((TH1D*)total->Clone(TString::Format("trc_%p", (void*)&a)));
    trunc->GetXaxis()->SetRangeUser(-PASS_SIGMA, PASS_SIGMA);
    s.truncRMS = trunc->GetRMS();
  }

  std::unique_ptr<TH1D> low(stackTotal(a.inclusiveResoLowTrackSig.get(),
                                       a.inclusiveResoLowTrackMix.get(),
                                       a.inclusiveResoLowTrackBkg.get(),
                                       TString::Format("low_%p", (void*)&a)));
  if (low->GetEntries() > 0) {
    std::unique_ptr<TF1> fitLow(ACTIVE_FIT_MODEL.makeFit(low.get()));
    s.sigmaLow = fitLow->GetParameter("Sigma1");
  }

  // Core fraction: same definition inclusivePlot() uses, under/overflow in.
  PlotObj* p = a.ptrHSTrack;
  if (p && p->effTotal && p->effPass) {
    int n = p->effTotal->GetNbinsX();
    double tot  = p->effTotal->Integral(0, n + 1);
    double pass = p->effPass ->Integral(0, n + 1);
    if (tot > 0.0) {
      double f = pass / tot;
      s.coreFrac    = 100.0 * f;
      s.coreFracErr = 100.0 * std::sqrt(f * (1.0 - f) / tot);
      s.nEvents     = tot;
    }
  }
  return s;
}

// Per-bin core fraction as a TH1D with binomial errors. Cheaper to overlay
// than TEfficiency (which needs the graph-conversion dance to suppress
// ROOT's repaint-time x-error restoration) and carries the same numbers.
TH1D* coreFractionHist(PlotObj* p, const char* name) {
  TH1D* h = (TH1D*)p->effPass->Clone(name);
  h->Divide(p->effPass.get(), p->effTotal.get(), 1.0, 1.0, "B");
  h->SetTitle("");
  // Bins with an empty denominator are undefined, not zero. The event
  // selection requires at least one forward jet, so the n-forward-jets = 0 bin
  // is always empty, and Divide leaves it at 0 -- which would draw a marker at
  // 0% and read as "this rule fails every such event". Park those bins far
  // below the drawn y-range instead so they are simply absent from the plot;
  // the tables print them as "--" by testing effTotal directly.
  for (int b = 1; b <= h->GetNbinsX(); ++b)
    if (p->effTotal->GetBinContent(b) <= 0.0) {
      h->SetBinContent(b, -999.0);
      h->SetBinError  (b, 0.0);
    }
  return h;
}

// Last bin of a PlotObj's x-axis that plotPostProcessing actually filled: it
// stops once it reaches the folded overflow bin, so bins past the fold are
// empty by construction rather than by lack of statistics.
int lastFilledBin(PlotObj* p) {
  TH1D* h = p->params->sigmaDist;
  for (int i = 1; i <= h->GetNbinsX(); ++i)
    if (h->GetBinLowEdge(i) >= p->foldValue) return i;
  return h->GetNbinsX();
}

// ---------------------------------------------------------------------------
// overlay
//   One page: every rule's curve for a single (score, variable, metric),
//   labelled by rule. moneyPlot cannot be reused here -- it labels its legend
//   entries with the SCORE name, which is identical for all seven curves in
//   this study and would produce seven indistinguishable legend rows.
// ---------------------------------------------------------------------------
void overlayPage(const char* fname, TCanvas* canvas,
                 const std::vector<TH1D*>& curves,
                 const std::vector<std::string>& labels,
                 const char* title, const char* xtitle, const char* ytitle,
                 double xMin, double xMax, double yMin, double yMax,
                 const char* scoreLabel, const char* pageLabel,
                 bool legendBottom = false) {
  canvas->Clear();
  // Core-fraction curves rise towards the right and fill the top of the pad;
  // sigma curves fall and fill the left. One legend corner does not suit both.
  TLegend* leg = legendBottom ? new TLegend(0.58, 0.18, 0.93, 0.48)
                              : new TLegend(0.58, 0.62, 0.93, 0.92);
  StyleLegend(leg, 0.024);

  bool first = true;
  for (size_t i = 0; i < curves.size(); ++i) {
    TH1D* h = curves[i];
    if (!h) continue;
    Color_t c = COLORS[i % COLORS.size()];
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(20 + (int)(i % 8));
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->SetTitle(title);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetRangeUser(xMin, xMax);
    h->GetYaxis()->SetRangeUser(yMin, yMax);
    h->Draw(first ? "E1" : "E1 SAME");
    first = false;
    leg->AddEntry(h, labels[i].c_str(), "lep");
  }
  if (first) return;  // nothing drawable

  ATLASLabel(0.18, 0.88, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
  TLatex latex;
  latex.SetTextSize(0.035);
  latex.SetTextAlign(13);
  latex.DrawLatexNDC(0.20, 0.76, scoreLabel);
  latex.DrawLatexNDC(0.20, 0.70, pageLabel);
  leg->Draw("SAME");
  canvas->Print(fname);
}

// Highest bin content across a set of curves, ignoring empty bins, so the
// four differential pages can each pick a sensible ceiling instead of
// inheriting one from whichever curve happened to be drawn first.
double curvesMax(const std::vector<TH1D*>& curves) {
  double m = 0.0;
  for (TH1D* h : curves) {
    if (!h) continue;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
      m = std::max(m, h->GetBinContent(i) + h->GetBinError(i));
  }
  return m;
}

// Lowest drawn value, skipping the -999 sentinels that mark undefined bins.
double curvesMin(const std::vector<TH1D*>& curves) {
  double m = 1e30;
  for (TH1D* h : curves) {
    if (!h) continue;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      double v = h->GetBinContent(i);
      if (v < -900.0) continue;
      m = std::min(m, v - h->GetBinError(i));
    }
  }
  return (m > 1e29) ? 0.0 : m;
}

}  // namespace

int main(int argc, char** argv) {
  // Detach histograms from gDirectory. Not cosmetic here: this stage books 14
  // AnalysisObj (7 rules x 2 scores), and with the default AddDirectory(kTRUE)
  // any of them constructed while a TFile is open become owned by that file
  // and are deleted when it closes -- leaving every unique_ptr in the map
  // dangling, which surfaces as a segfault inside plotPostProcessing rather
  // than anywhere near the actual mistake. The booking below is also kept
  // outside the reader's scope (as clustering_plot.cxx does), so this is a
  // belt-and-braces guard against re-introducing the same bug.
  TH1::AddDirectory(kFALSE);

  SetAtlasStyle();

  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/assoc_study");

  gErrorIgnoreLevel = kFatal;

  const auto& rules  = assocRules();
  const auto& scores = assocScores();

  // --- Load ---
  const std::string histFile = MyUtl::resolveHistFile(
    argc, argv, MyUtl::histFilePath(ASSOC_HIST_BASENAME));
  std::cout << "Reading histograms from " << histFile << '\n';

  // Construct-empty-then-load: book every histogram with the SAME code the
  // hist stage used (so names and binning cannot drift, and HistReader::LoadInto
  // catches it if they ever do), all of it before the file is opened.
  std::vector<std::map<Score, AnalysisObj>> maps;
  maps.reserve(rules.size());
  for (const auto& r : rules) maps.push_back(buildAssocMap(r));

  std::vector<Long64_t> nKept(rules.size()), nKeptHS(rules.size()), nKeptPU(rules.size());
  Long64_t nAccepted = 0, nFwdAvail = 0, nFwdAvailHS = 0;
  {
    MyUtl::HistReader reader(histFile);
    reader.CheckSelection(MyUtl::VBS_JET_D_ETA, MyUtl::VBS_JET_MJJ);
    MyUtl::ENERGY_LABEL = reader.ReadEnergyLabel();
    for (auto& m : maps)
      for (auto& [score, analysis] : m) analysis.loadFrom(reader);
    for (size_t i = 0; i < rules.size(); ++i) {
      nKept  [i] = reader.ReadScalarLong(assocScalarKey(rules[i], "n_kept"));
      nKeptHS[i] = reader.ReadScalarLong(assocScalarKey(rules[i], "n_kept_hs"));
      nKeptPU[i] = reader.ReadScalarLong(assocScalarKey(rules[i], "n_kept_pu"));
    }
    nAccepted   = reader.ReadScalarLong("meta_n_accepted");
    nFwdAvail   = reader.ReadScalarLong("meta_n_fwd_avail");
    nFwdAvailHS = reader.ReadScalarLong("meta_n_fwd_avail_hs");
  }

  for (auto& m : maps)
    for (auto& [score, analysis] : m) analysis.postProcessing();

  // --- Inclusive summaries: [rule][score] ---
  std::vector<std::vector<Summary>> sums(rules.size());
  for (size_t i = 0; i < rules.size(); ++i)
    for (const auto& s : scores)
      sums[i].push_back(summarize(maps[i].at(s)));

  // ==========================================================================
  // Console report
  // ==========================================================================
  const std::string BAR(96, '=');
  std::cout << '\n' << BAR << '\n'
            << "TRACK-TO-VERTEX z-ASSOCIATION STUDY";
  if (!MyUtl::SAMPLE_NAME.empty()) std::cout << "   sample: " << MyUtl::SAMPLE_NAME;
  std::cout << "\n" << BAR << '\n'
            << "Events passing selection: " << nAccepted << '\n'
            << "Counting scan pinned at z0 significance < " << COUNTING_NSIGMA
            << " for every rule (common plot x-axes)\n";

  // --- Track accounting ---
  std::cout << "\n--- TRACK ACCOUNTING (forward, pT/quality-selected) ---\n"
            << "available: " << nFwdAvail << " tracks, of which "
            << nFwdAvailHS << " truth-HS\n\n"
            << std::left << std::setw(24) << "rule"
            << std::right << std::setw(13) << "kept"
            << std::setw(13) << "HS kept"
            << std::setw(13) << "PU kept"
            << std::setw(10) << "purity"
            << std::setw(10) << "recall" << '\n'
            << std::string(75, '-') << '\n';
  for (size_t i = 0; i < rules.size(); ++i) {
    double pur = nKept[i]    > 0 ? 100.0 * nKeptHS[i] / nKept[i]    : 0.0;
    double rec = nFwdAvailHS > 0 ? 100.0 * nKeptHS[i] / nFwdAvailHS : 0.0;
    std::cout << std::left << std::setw(24) << ruleAscii(rules[i])
              << std::right << std::setw(13) << nKept[i]
              << std::setw(13) << nKeptHS[i]
              << std::setw(13) << nKeptPU[i]
              << std::fixed << std::setprecision(2)
              << std::setw(9) << pur << "%"
              << std::setw(9) << rec << "%" << '\n';
  }

  // --- Inclusive metrics, per score ---
  for (size_t si = 0; si < scores.size(); ++si) {
    std::cout << "\n--- INCLUSIVE:  " << scores[si].toStringShort() << " ---\n"
              << std::left << std::setw(24) << "rule"
              << std::right << std::setw(16) << "core sigma[ps]"
              << std::setw(16) << "core frac[%]"
              << std::setw(14) << "truncRMS[ps]"
              << std::setw(16) << "sigma(nHS<=5)"
              << std::setw(12) << "N events" << '\n'
              << std::string(90, '-') << '\n';
    for (size_t i = 0; i < rules.size(); ++i) {
      const Summary& s = sums[i][si];
      std::ostringstream sig, cf;
      sig << std::fixed << std::setprecision(2) << s.sigma << " +/- " << s.sigmaErr;
      cf  << std::fixed << std::setprecision(2) << s.coreFrac << " +/- " << s.coreFracErr;
      std::cout << std::left << std::setw(24) << ruleAscii(rules[i])
                << std::right << std::setw(16) << sig.str()
                << std::setw(16) << cf.str()
                << std::fixed << std::setprecision(2)
                << std::setw(14) << s.truncRMS
                << std::setw(16) << s.sigmaLow
                << std::setprecision(0)
                << std::setw(12) << s.nEvents << '\n';
    }
    // Deltas against the incumbent (rule 0 = significance < 3.0)
    std::cout << "\n  vs. incumbent (" << ruleAscii(rules[0]) << "):\n";
    for (size_t i = 1; i < rules.size(); ++i) {
      double dS = sums[i][si].sigma    - sums[0][si].sigma;
      double dF = sums[i][si].coreFrac - sums[0][si].coreFrac;
      std::cout << "    " << std::left << std::setw(24) << ruleAscii(rules[i])
                << std::right << std::showpos << std::fixed << std::setprecision(2)
                << std::setw(9) << dS << " ps"
                << std::setw(9) << dF << " %pt" << std::noshowpos << '\n';
    }
  }

  // --- Differential tables + plots ---
  struct VarSpec { const char* key; const char* xtitle; };
  const std::vector<VarSpec> VARS = {
    { "hs_track", "n Forward HS Tracks" },
    { "fjet",     "n Forward Jets"      },
  };

  TCanvas* canvas = new TCanvas("canvas", "assoc study", 800, 600);

  for (size_t si = 0; si < scores.size(); ++si) {
    const Score& sc = scores[si];
    for (const auto& v : VARS) {
      PlotObj* ref = maps[0].at(sc).get(v.key).get();
      const int nb = lastFilledBin(ref);

      // ---- core sigma vs variable ----
      std::cout << "\n--- CORE SIGMA [ps] vs " << v.xtitle
                << "   (" << sc.toStringShort() << ") ---\n"
                << std::left << std::setw(24) << "rule";
      for (int b = 1; b <= nb; ++b) {
        std::ostringstream lab;
        int lo = (int)std::lround(ref->params->sigmaDist->GetBinLowEdge(b) + 0.5);
        lab << (b == nb ? ">=" : "") << lo;
        std::cout << std::right << std::setw(9) << lab.str();
      }
      std::cout << '\n' << std::string(24 + 9 * nb, '-') << '\n';

      // Denominator per bin, printed once per table. The counting scan is
      // pinned, so this row is identical for every rule -- which is the point:
      // it shows where a per-bin difference between rules can be believed and
      // where it is a handful of events talking.
      std::cout << std::left << std::setw(24) << "  [N events]" << std::right;
      for (int b = 1; b <= nb; ++b)
        std::cout << std::setw(9) << (Long64_t)ref->effTotal->GetBinContent(b);
      std::cout << '\n' << std::string(24 + 9 * nb, '-') << '\n';

      std::vector<TH1D*> sigCurves, cfCurves;
      std::vector<std::string> labels;
      for (size_t i = 0; i < rules.size(); ++i) {
        PlotObj* p = maps[i].at(sc).get(v.key).get();
        TH1D* sd = p->params->sigmaDist;
        std::cout << std::left << std::setw(24) << ruleAscii(rules[i]) << std::right;
        for (int b = 1; b <= nb; ++b) {
          double val = sd->GetBinContent(b);
          if (val == 0.0) std::cout << std::setw(9) << "--";
          else            std::cout << std::fixed << std::setprecision(1)
                                    << std::setw(9) << val;
        }
        std::cout << '\n';
        // Clone before blanking: sigmaDist is owned by the PlotObj's FitParams,
        // and a bin left empty by plotPostProcessing's fit-quality gate means
        // "no usable fit here", not "sigma = 0". Drawn as-is it would put a
        // marker on the axis and pull the y-range to zero.
        TH1D* sdDraw = (TH1D*)sd->Clone(
          TString::Format("sig_%s_%s_%zu", sc.toStringShort(), v.key, i));
        for (int b = 1; b <= sdDraw->GetNbinsX(); ++b)
          if (sdDraw->GetBinContent(b) == 0.0) sdDraw->SetBinContent(b, -999.0);
        sigCurves.push_back(sdDraw);
        cfCurves.push_back(coreFractionHist(
          p, TString::Format("cf_%s_%s_%zu", sc.toStringShort(), v.key, i)));
        labels.push_back(rules[i].label);
      }

      // ---- core fraction vs variable ----
      std::cout << "\n--- CORE FRACTION [%] vs " << v.xtitle
                << "   (" << sc.toStringShort() << ") ---\n"
                << std::left << std::setw(24) << "rule";
      for (int b = 1; b <= nb; ++b) {
        std::ostringstream lab;
        int lo = (int)std::lround(ref->effTotal->GetBinLowEdge(b) + 0.5);
        lab << (b == nb ? ">=" : "") << lo;
        std::cout << std::right << std::setw(9) << lab.str();
      }
      std::cout << '\n' << std::string(24 + 9 * nb, '-') << '\n';
      std::cout << std::left << std::setw(24) << "  [N events]" << std::right;
      for (int b = 1; b <= nb; ++b)
        std::cout << std::setw(9) << (Long64_t)ref->effTotal->GetBinContent(b);
      std::cout << '\n' << std::string(24 + 9 * nb, '-') << '\n';
      for (size_t i = 0; i < rules.size(); ++i) {
        PlotObj* p = maps[i].at(sc).get(v.key).get();
        std::cout << std::left << std::setw(24) << ruleAscii(rules[i]) << std::right;
        for (int b = 1; b <= nb; ++b) {
          if (p->effTotal->GetBinContent(b) <= 0.0) std::cout << std::setw(9) << "--";
          else std::cout << std::fixed << std::setprecision(1)
                         << std::setw(9) << 100.0 * cfCurves[i]->GetBinContent(b);
        }
        std::cout << '\n';
      }
      // Core fraction histograms are 0-1; scale to percent for drawing. The
      // -999 sentinels set by coreFractionHist stay far out of the drawn range.
      for (TH1D* h : cfCurves) h->Scale(100.0);

      // Eleven curves on one pad are unreadable, and the two rule families are
      // the natural split: within a family the curves form an ordered
      // tightening sequence, so the reader is looking for a trend, not for one
      // curve among eleven. Third page pairs the incumbent against the best
      // point of each family (best = highest inclusive core fraction for THIS
      // score, computed once from the inclusive table above).
      // Grouping keys on (kind, orGhost), not kind alone: a union rule shares
      // its base kind with the pure z-rules but is answering a different
      // question -- widen the list rather than tighten it -- so it belongs on
      // its own page next to the base rules it is built from, not buried in
      // an ordered tightening sequence it does not belong to.
      auto pageIdx = [&](AssocRule::Kind k, bool ghost) {
        std::vector<size_t> v2;
        for (size_t i = 0; i < rules.size(); ++i)
          if (rules[i].kind == k && rules[i].orGhost == ghost) v2.push_back(i);
        return v2;
      };
      auto unionPage = [&]() {
        std::vector<size_t> v2;
        // Each union rule preceded by the base rule it widens, so the cost of
        // the union is readable off the page without cross-referencing.
        for (size_t i = 0; i < rules.size(); ++i)
          if (rules[i].orGhost) {
            for (size_t b = 0; b < rules.size(); ++b)
              if (!rules[b].orGhost && rules[b].kind == rules[i].kind &&
                  rules[b].cut == rules[i].cut) v2.push_back(b);
            v2.push_back(i);
          }
        return v2;
      };
      auto bestOf = [&](AssocRule::Kind k) {
        size_t best = 0; double bestCF = -1.0;
        for (size_t i = 0; i < rules.size(); ++i)
          if (rules[i].kind == k && !rules[i].orGhost && sums[i][si].coreFrac > bestCF) {
            bestCF = sums[i][si].coreFrac; best = i;
          }
        return best;
      };
      struct PageSpec { std::string title; std::vector<size_t> idx; };
      std::vector<PageSpec> pages = {
        { "z_{0}-significance family",  pageIdx(AssocRule::Kind::SIGNIFICANCE, false) },
        { "parameterisation family",    pageIdx(AssocRule::Kind::DZ_PARA,      false) },
        { "union with ghost assoc.",    unionPage()                                   },
        { "incumbent vs. best of each",
          { 0, bestOf(AssocRule::Kind::SIGNIFICANCE), bestOf(AssocRule::Kind::DZ_PARA) } },
      };

      const double xLo = ref->xMin, xHi = ref->foldMax;
      // y-ceilings computed over ALL curves, not per page, so the three pages
      // of one PDF share a scale and can be flipped between.
      const double sigMax = 1.35 * curvesMax(sigCurves);

      auto emit = [&](const char* stem, const std::vector<TH1D*>& all,
                      const char* ytitle, double yLo, double yHi, bool legBottom) {
        TString fname = MyUtl::plotFilePath(
          "assoc_study",
          TString::Format("%s_%s_%s.pdf", stem, sc.toStringShort(), v.key).Data()).c_str();
        canvas->Print(fname + "[");
        for (const auto& pg : pages) {
          std::vector<TH1D*> sub;
          std::vector<std::string> subLabels;
          for (size_t i : pg.idx) { sub.push_back(all[i]); subLabels.push_back(labels[i]); }
          overlayPage(fname, canvas, sub, subLabels, "", v.xtitle, ytitle,
                      xLo, xHi, yLo, yHi, sc.toStringShort(), pg.title.c_str(), legBottom);
        }
        canvas->Print(fname + "]");
      };

      // Core fraction lives in a narrow band near the top (roughly 72-100% vs
      // n HS tracks, 88-95% vs n forward jets). Drawn on a 0-100 axis every
      // rule looks identical, which is exactly the difference being measured,
      // so the axis is scaled to the data with a little headroom instead.
      const double cfLo = std::max(0.0, curvesMin(cfCurves) - 3.0);
      const double cfHi = std::min(103.0, curvesMax(cfCurves) + 2.0);

      emit("assoc_sigma",    sigCurves, "Core #sigma(#Delta t) [ps]", 0.0, sigMax, false);
      emit("assoc_corefrac", cfCurves,  "Core Fraction [%]",          cfLo, cfHi,  true);
    }
  }

  // --- Inclusive Delta t shape overlay, one page per score ---
  for (const auto& sc : scores) {
    std::vector<TH1D*> curves;
    std::vector<std::string> labels;
    for (size_t i = 0; i < rules.size(); ++i) {
      AnalysisObj& a = maps[i].at(sc);
      TH1D* t = stackTotal(a.inclusiveResoSig.get(), a.inclusiveResoMix.get(),
                           a.inclusiveResoBkg.get(),
                           TString::Format("shape_%s_%zu", sc.toStringShort(), i));
      if (t->Integral() > 0) t->Scale(1.0 / t->Integral());
      curves.push_back(t);
      labels.push_back(rules[i].label);
    }
    TString fname = MyUtl::plotFilePath(
      "assoc_study",
      TString::Format("assoc_dt_shape_%s.pdf", sc.toStringShort()).Data()).c_str();
    const double yHi = 1.35 * curvesMax(curves);
    canvas->Print(fname + "[");
    struct ShapePage { AssocRule::Kind kind; bool ghost; const char* title; };
    for (auto pg : { ShapePage{AssocRule::Kind::SIGNIFICANCE, false, "z_{0}-significance family"},
                     ShapePage{AssocRule::Kind::DZ_PARA,      false, "parameterisation family"},
                     ShapePage{AssocRule::Kind::SIGNIFICANCE, true,  "union with ghost assoc."},
                     ShapePage{AssocRule::Kind::DZ_PARA,      true,  "union with ghost assoc."} }) {
      std::vector<TH1D*> sub;
      std::vector<std::string> subLabels;
      for (size_t i = 0; i < rules.size(); ++i)
        if (rules[i].kind == pg.kind && rules[i].orGhost == pg.ghost) {
          sub.push_back(curves[i]); subLabels.push_back(labels[i]);
        }
      if (sub.empty()) continue;
      overlayPage(fname, canvas, sub, subLabels, "", "#Delta t [ps]", "Normalised Entries",
                  -200.0, 200.0, 0.0, yHi, sc.toStringShort(), pg.title);
    }
    canvas->Print(fname + "]");
  }

  // ==========================================================================
  // Machine-readable dump — one row per (rule, score), for the write-up.
  // ==========================================================================
  std::cout << "\nCSV|rule,score,kept,kept_hs,kept_pu,purity_pct,recall_pct,"
               "sigma_ps,sigma_err_ps,corefrac_pct,corefrac_err_pct,trunc_rms_ps,sigma_lowmult_ps,n_events\n";
  for (size_t i = 0; i < rules.size(); ++i) {
    double pur = nKept[i]    > 0 ? 100.0 * nKeptHS[i] / nKept[i]    : 0.0;
    double rec = nFwdAvailHS > 0 ? 100.0 * nKeptHS[i] / nFwdAvailHS : 0.0;
    for (size_t si = 0; si < scores.size(); ++si) {
      const Summary& s = sums[i][si];
      std::cout << "CSV|" << rules[i].tag << ',' << scores[si].toStringShort() << ','
                << nKept[i] << ',' << nKeptHS[i] << ',' << nKeptPU[i] << ','
                << std::fixed << std::setprecision(4)
                << pur << ',' << rec << ','
                << s.sigma << ',' << s.sigmaErr << ','
                << s.coreFrac << ',' << s.coreFracErr << ','
                << s.truncRMS << ',' << s.sigmaLow << ','
                << std::setprecision(0) << s.nEvents << '\n';
    }
  }

  std::cout << "\nPlots written to " << MyUtl::OUTPUT_DIR << "/assoc_study/\n";
  return 0;
}
