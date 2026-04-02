#ifndef SWEEP_UTILITIES_H
#define SWEEP_UTILITIES_H

// sweep_utilities.h
//
// Shared includes, constants, and plotting helpers for the 1D parameter sweep
// executables (distcut_sweep, mintrackpt_sweep, maxtrackpt_sweep, …).
//
// Contents:
//   Constants        — N_REGIONS, RES_MIN/MAX/BWID, OUT_DIR
//   fitCoreSigma     — double-Gaussian core-sigma fit for residual histograms
//   LegendBox        — NDC bounding box for a legend
//   bestLegendCorner — choose the least-occupied plot corner for the legend
//   saveGraphPair    — style, draw, and print two TGraphErrors to a PDF page

#include <TROOT.h>
#include <TStyle.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

// ── Shared sweep constants ────────────────────────────────────────────────────
static constexpr int    N_REGIONS = 2;        // 0 = low HS, 1 = high HS
static constexpr double RES_MIN   = -1000.0;  // residual histogram range [ps]
static constexpr double RES_MAX   =  1000.0;
static constexpr double RES_BWID  =     2.0;  // bin width [ps]
static const char*      OUT_DIR   = "../parameter_sweep";
// ─────────────────────────────────────────────────────────────────────────────

// ── fitCoreSigma ──────────────────────────────────────────────────────────────
// Fit a residual histogram with a double Gaussian and return
// {core sigma, core sigma error}.  Returns {0,0} for fewer than 20 entries.
// ─────────────────────────────────────────────────────────────────────────────
inline std::pair<double,double> fitCoreSigma(TH1D* h, int idx) {
  if (h->GetEntries() < 20) return {0.0, 0.0};
  TF1* fit = new TF1(
    Form("dgaus_fit_%d", idx),
    "[1]*TMath::Exp(-0.5*((x-[0])/[3])^2) + [2]*TMath::Exp(-0.5*((x-[0])/[4])^2)",
    RES_MIN, RES_MAX);
  fit->SetParameters(
    0,
    0.8  * h->GetMaximum(),
    0.01 * h->GetMaximum(),
    26.0,
    MyUtl::PILEUP_SMEAR);
  fit->FixParameter(0, 0.0);
  fit->SetParLimits(1, 0, 1e6);
  fit->SetParLimits(2, 0, 1e6);
  fit->SetParLimits(3, 0.1, MyUtl::PILEUP_SMEAR);
  fit->FixParameter(4, MyUtl::PILEUP_SMEAR);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(2);
  fit->SetNpx(1000);
  h->Fit(fit, "RQ");
  double sigma    = fit->GetParameter(3);
  double sigmaErr = fit->GetParError(3);
  delete fit;
  return {sigma, sigmaErr};
}

// ── bestLegendCorner ──────────────────────────────────────────────────────────
// Returns the NDC bounding box (x1, y1, x2, y2) for the plot corner with the
// fewest data points (including ±1σ error bars) overlapping the legend area.
// Examines four corners: top-right (conventional default), top-left,
// bottom-right, bottom-left.  Ties are broken in that order of preference.
//
//   gr1, gr2     — the two data series to inspect
//   xAxisLo/Hi   — x-axis limits passed to GetXaxis()->SetLimits()
//   yMin, yMax   — y-axis range passed to GetYaxis()->SetRangeUser()
//   hasRefLine   — true when a reference-line entry will be in the legend
//                  (adds a small amount of extra height)
// ─────────────────────────────────────────────────────────────────────────────
struct LegendBox { double x1, y1, x2, y2; };

inline LegendBox bestLegendCorner(
  TGraphErrors* gr1, TGraphErrors* gr2,
  double xAxisLo, double xAxisHi,
  double yMin, double yMax,
  bool hasRefLine = false)
{
  // Canvas margins assumed by all sweep plots (canvas->SetLeftMargin(0.15))
  const double leftM = 0.15, rightM = 0.05;
  const double botM  = 0.10, topM   = 0.05;
  const double padW  = 1.0 - leftM - rightM;
  const double padH  = 1.0 - botM  - topM;

  // Legend dimensions matching the user's preferred size from distcut_sweep:
  //   width = 0.88 - 0.55 = 0.33,  height = 0.88 - 0.68 = 0.20 (2 entries)
  //                                          0.88 - 0.66 = 0.22 (3 entries)
  const double legW = 0.33;
  const double legH = hasRefLine ? 0.22 : 0.20;

  // Four candidate placements ordered by preference
  const LegendBox CORNERS[4] = {
    {0.55,        1.0-topM-legH, 0.88,        1.0-topM       },  // top-right
    {leftM+0.01,  1.0-topM-legH, leftM+0.34,  1.0-topM       },  // top-left
    {0.55,        botM+0.02,     0.88,        botM+legH+0.02  },  // bottom-right
    {leftM+0.01,  botM+0.02,     leftM+0.34,  botM+legH+0.02 },  // bottom-left
  };

  // Convert from data coordinates to NDC, clamped to the pad extent
  auto toNDCx = [&](double x) -> double {
    return leftM + (x - xAxisLo) / (xAxisHi - xAxisLo) * padW;
  };
  auto toNDCy = [&](double y) -> double {
    double v = botM + (y - yMin) / (yMax - yMin) * padH;
    return std::max(botM, std::min(1.0 - topM, v));
  };

  int bestIdx   = 0;
  int bestScore = INT_MAX;

  for (int c = 0; c < 4; ++c) {
    const auto& box = CORNERS[c];
    int score = 0;
    for (auto* gr : {gr1, gr2}) {
      for (int i = 0; i < gr->GetN(); ++i) {
        double gx, gy;
        gr->GetPoint(i, gx, gy);
        double nx   = toNDCx(gx);
        double nyLo = toNDCy(gy - gr->GetErrorY(i));
        double nyHi = toNDCy(gy + gr->GetErrorY(i));
        bool xHit = (nx >= box.x1 && nx <= box.x2);
        bool yHit = (nyHi >= box.y1 && nyLo <= box.y2);
        if (xHit && yHit) ++score;
      }
    }
    if (score < bestScore) { bestScore = score; bestIdx = c; }
  }
  return CORNERS[bestIdx];
}

// ── saveGraphPair ─────────────────────────────────────────────────────────────
// Style, draw, and overlay two TGraphErrors (low/high HS regions) on canvas,
// place the legend in the least-occupied corner, and print the page to fname.
// A horizontal reference line is added when refVal >= 0 and refLabel != nullptr.
//
//   xtitle         — x-axis label (e.g. "Distance Cut [#sigma]")
//   xAxisLo/Hi     — axis limits (used for SetLimits and the reference line)
//   hsSplit        — HS_TRACK_SPLIT value for the legend labels
// ─────────────────────────────────────────────────────────────────────────────
inline void saveGraphPair(
  const std::string& fname,
  TCanvas*      canvas,
  TGraphErrors* grLow,    // region 0 (low HS)
  TGraphErrors* grHigh,   // region 1 (high HS)
  const char*   title,
  const char*   xtitle,
  const char*   ytitle,
  double        xAxisLo,
  double        xAxisHi,
  double        yMin,
  double        yMax,
  int           hsSplit,
  double        refVal    = -1,      // horizontal reference line (y value); <0 = off
  const char*   refLabel  = nullptr,
  double        vRefVal   = -1,      // vertical reference line (x value); <0 = off
  const char*   vRefLabel = nullptr)
{
  canvas->cd();

  grLow->SetTitle(Form("%s;%s;%s", title, xtitle, ytitle));
  grLow->SetLineColor(MyUtl::C01);
  grLow->SetMarkerColor(MyUtl::C01);
  grLow->SetMarkerStyle(20);
  grLow->SetMarkerSize(0.9);
  grLow->SetLineWidth(2);

  grHigh->SetLineColor(MyUtl::C02);
  grHigh->SetMarkerColor(MyUtl::C02);
  grHigh->SetMarkerStyle(21);
  grHigh->SetMarkerSize(0.9);
  grHigh->SetLineWidth(2);

  grLow->Draw("APE1");
  grLow->GetXaxis()->SetLimits(xAxisLo, xAxisHi);
  grLow->GetYaxis()->SetRangeUser(yMin, yMax);
  grLow->GetXaxis()->SetNdivisions(510);
  grHigh->Draw("PE1 SAME");

  bool hasRef = (refVal >= 0 && refLabel != nullptr);
  LegendBox lb = bestLegendCorner(grLow, grHigh, xAxisLo, xAxisHi, yMin, yMax, hasRef);
  TLegend* leg = new TLegend(lb.x1, lb.y1, lb.x2, lb.y2);
  leg->AddEntry(grLow,  Form("nFTrackHS #leq %d", hsSplit), "lp");
  leg->AddEntry(grHigh, Form("nFTrackHS > %d",    hsSplit), "lp");

  if (hasRef) {
    TLine* refLine = new TLine(xAxisLo, refVal, xAxisHi, refVal);
    refLine->SetLineColor(kRed);
    refLine->SetLineWidth(2);
    refLine->SetLineStyle(4);
    refLine->Draw("SAME");
    leg->AddEntry(refLine, refLabel, "l");
  }

  bool hasVRef = (vRefVal >= 0 && vRefLabel != nullptr);
  if (hasVRef) {
    TLine* vLine = new TLine(vRefVal, yMin, vRefVal, yMax);
    vLine->SetLineColor(kGray + 2);
    vLine->SetLineWidth(2);
    vLine->SetLineStyle(7);
    vLine->Draw("SAME");
    leg->AddEntry(vLine, vRefLabel, "l");
  }

  leg->Draw("SAME");

  canvas->Print(fname.c_str());
}

// ── bestLegendCorner (TGraph overload) ────────────────────────────────────────
// Same four-corner selection as above but for plain TGraph (no error bars).
// An optional list of extra data-coordinate positions (e.g. text-label anchor
// points) are also scored against each candidate box so the chosen corner
// avoids both markers and labels.
// ─────────────────────────────────────────────────────────────────────────────
inline LegendBox bestLegendCorner(
  TGraph* gr1, TGraph* gr2,
  double xAxisLo, double xAxisHi,
  double yMin, double yMax,
  const std::vector<std::pair<double,double>>& extraPts = {})
{
  const double leftM = 0.15, rightM = 0.05;
  const double botM  = 0.10, topM   = 0.05;
  const double padW  = 1.0 - leftM - rightM;
  const double padH  = 1.0 - botM  - topM;
  const double legH  = 0.20;

  const LegendBox CORNERS[4] = {
    {0.55,        1.0-topM-legH, 0.88,        1.0-topM       },  // top-right
    {leftM+0.01,  1.0-topM-legH, leftM+0.34,  1.0-topM       },  // top-left
    {0.55,        botM+0.02,     0.88,        botM+legH+0.02  },  // bottom-right
    {leftM+0.01,  botM+0.02,     leftM+0.34,  botM+legH+0.02 },  // bottom-left
  };

  auto toNDCx = [&](double x) -> double {
    return leftM + (x - xAxisLo) / (xAxisHi - xAxisLo) * padW;
  };
  auto toNDCy = [&](double y) -> double {
    double v = botM + (y - yMin) / (yMax - yMin) * padH;
    return std::max(botM, std::min(1.0 - topM, v));
  };
  auto hitBox = [&](const LegendBox& box, double gx, double gy) -> bool {
    double nx = toNDCx(gx), ny = toNDCy(gy);
    return (nx >= box.x1 && nx <= box.x2 && ny >= box.y1 && ny <= box.y2);
  };

  int bestIdx = 0, bestScore = INT_MAX;
  for (int c = 0; c < 4; ++c) {
    const auto& box = CORNERS[c];
    int score = 0;
    for (auto* gr : {gr1, gr2}) {
      for (int i = 0; i < gr->GetN(); ++i) {
        double gx, gy;
        gr->GetPoint(i, gx, gy);
        if (hitBox(box, gx, gy)) ++score;
      }
    }
    for (const auto& [px, py] : extraPts)
      if (hitBox(box, px, py)) ++score;
    if (score < bestScore) { bestScore = score; bestIdx = c; }
  }
  return CORNERS[bestIdx];
}

// ── bestLegendCornerHist ──────────────────────────────────────────────────────
// Select the least-occupied canvas corner for a legend on a 1D histogram
// overlay plot.  Each candidate corner is scored by summing the normalised bin
// content of every supplied histogram whose bin-centre NDC coordinate falls
// inside the candidate box.  The corner with the lowest cumulative score wins.
//
//   hists        — all histograms that will be drawn (raw counts; normalised
//                  internally per-histogram)
//   nLegEntries  — number of text entries in the legend (controls box height)
//   logY         — true when the plot uses a log Y axis (changes NDC mapping)
// ─────────────────────────────────────────────────────────────────────────────
inline LegendBox bestLegendCornerHist(
  const std::vector<TH1D*>& hists,
  double xAxisLo, double xAxisHi,
  double yMin,    double yMax,
  int    nLegEntries = 2,
  bool   logY        = false)
{
  const double leftM = 0.15, rightM = 0.05;
  const double botM  = 0.10, topM   = 0.05;
  const double padW  = 1.0 - leftM - rightM;
  const double padH  = 1.0 - botM  - topM;
  const double legH  = 0.055 + nLegEntries * 0.038;  // ~0.20 for 4 entries

  const LegendBox CORNERS[4] = {
    {0.55,        1.0-topM-legH, 0.88,        1.0-topM       },  // top-right
    {leftM+0.01,  1.0-topM-legH, leftM+0.34,  1.0-topM       },  // top-left
    {0.55,        botM+0.02,     0.88,        botM+legH+0.02  },  // bottom-right
    {leftM+0.01,  botM+0.02,     leftM+0.34,  botM+legH+0.02 },  // bottom-left
  };

  auto toNDCx = [&](double x) -> double {
    return leftM + (x - xAxisLo) / (xAxisHi - xAxisLo) * padW;
  };
  auto toNDCy = [&](double y) -> double {
    double v;
    if (logY && yMin > 0 && yMax > yMin)
      v = botM + std::log(y / yMin) / std::log(yMax / yMin) * padH;
    else
      v = botM + (y - yMin) / (yMax - yMin) * padH;
    return std::max(botM, std::min(1.0 - topM, v));
  };

  int    bestIdx   = 0;
  double bestScore = 1e30;
  for (int c = 0; c < 4; ++c) {
    const auto& box = CORNERS[c];
    double score = 0.0;
    for (auto* h : hists) {
      double norm = h->Integral();
      if (norm <= 0.0) continue;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double bx = h->GetBinCenter(b);
        double by = h->GetBinContent(b) / norm;
        double nx  = toNDCx(bx);
        double ny  = toNDCy(by);
        if (nx >= box.x1 && nx <= box.x2 && ny >= box.y1 && ny <= box.y2)
          score += by;
      }
    }
    if (score < bestScore) { bestScore = score; bestIdx = c; }
  }
  return CORNERS[bestIdx];
}

// ── graphYVals ────────────────────────────────────────────────────────────────
// Extract the Y values from a TGraphErrors into a std::vector<double>.
// ─────────────────────────────────────────────────────────────────────────────
inline std::vector<double> graphYVals(TGraphErrors* gr) {
  std::vector<double> v(gr->GetN());
  for (int i = 0; i < gr->GetN(); ++i) {
    double x, y;
    gr->GetPoint(i, x, y);
    v[i] = y;
  }
  return v;
}

// ── saveRocPlot ───────────────────────────────────────────────────────────────
// Draw a pseudo-ROC (metric-vs-metric) curve for the two HS regions.
// Each point corresponds to one sweep step and is labelled with its
// parameter value string.  Both regions are overlaid on one canvas page.
//
//   xLow / yLow   — (x,y) metric pairs for region 0 (low HS)
//   xHigh / yHigh — (x,y) metric pairs for region 1 (high HS)
//   labels        — one ROOT-LaTeX label string per point (shared by both)
//   hsSplit       — HS track boundary for the legend entries
// ─────────────────────────────────────────────────────────────────────────────
inline void saveRocPlot(
  const std::string&              fname,
  TCanvas*                        canvas,
  const std::vector<double>&      xLow,
  const std::vector<double>&      yLow,
  const std::vector<double>&      xHigh,
  const std::vector<double>&      yHigh,
  const std::vector<std::string>& labels,
  const char*                     title,
  const char*                     xtitle,
  const char*                     ytitle,
  int                             hsSplit)
{
  canvas->cd();
  const int nPts = (int)xLow.size();

  // Compute axis ranges across both regions
  auto bothMinMax = [](const std::vector<double>& a, const std::vector<double>& b)
    -> std::pair<double,double> {
    double lo = std::min(*std::min_element(a.begin(), a.end()),
                         *std::min_element(b.begin(), b.end()));
    double hi = std::max(*std::max_element(a.begin(), a.end()),
                         *std::max_element(b.begin(), b.end()));
    return {lo, hi};
  };

  auto [xLo, xHi] = bothMinMax(xLow, xHigh);
  auto [yLo, yHi] = bothMinMax(yLow, yHigh);
  const double dxPad = 0.10 * (xHi - xLo + 1e-9);
  const double dyPad = 0.10 * (yHi - yLo + 1e-9);
  xLo -= dxPad;
  xHi += 2.0 * dxPad;   // extra right margin for end-point label
  yLo -= dyPad;
  yHi += 1.5 * dyPad;   // extra top  margin for first-point label

  TGraph* grLow  = new TGraph(nPts, xLow.data(),  yLow.data());
  TGraph* grHigh = new TGraph(nPts, xHigh.data(), yHigh.data());

  grLow->SetTitle(Form("%s;%s;%s", title, xtitle, ytitle));
  grLow->SetLineColor(MyUtl::C01);
  grLow->SetMarkerColor(MyUtl::C01);
  grLow->SetMarkerStyle(20);
  grLow->SetMarkerSize(0.9);
  grLow->SetLineWidth(2);

  grHigh->SetLineColor(MyUtl::C02);
  grHigh->SetMarkerColor(MyUtl::C02);
  grHigh->SetMarkerStyle(21);
  grHigh->SetMarkerSize(0.9);
  grHigh->SetLineWidth(2);

  grLow->Draw("APL");
  grLow->GetXaxis()->SetLimits(xLo, xHi);
  grLow->GetYaxis()->SetRangeUser(yLo, yHi);
  grHigh->Draw("PL SAME");

  // Label only the first and last points to avoid clutter on dense curves
  TLatex tex;
  tex.SetTextSize(0.022);
  tex.SetTextFont(42);
  const double offY = 0.022 * (yHi - yLo);
  const double offX = 0.008 * (xHi - xLo);

  // first point: label above (+1.0);  last point: label below (-1.8)
  const int    labelIdx[2] = { 0, nPts - 1 };
  const double labelSign[2] = { 1.0, -1.8 };

  // Collect label anchor positions so the legend avoids both markers and labels
  std::vector<std::pair<double,double>> labelPts;
  for (int ii = 0; ii < 2; ++ii) {
    const int i = labelIdx[ii];
    if (i < 0 || i >= nPts) continue;
    const double s = labelSign[ii];
    labelPts.push_back({xLow[i]  + offX, yLow[i]  + s * offY});
    labelPts.push_back({xHigh[i] + offX, yHigh[i] + s * offY});
  }

  for (int ii = 0; ii < 2; ++ii) {
    const int i = labelIdx[ii];
    if (i < 0 || i >= nPts) continue;
    const double s = labelSign[ii];
    tex.SetTextColor(MyUtl::C01);
    tex.DrawLatex(xLow[i]  + offX, yLow[i]  + s * offY, labels[i].c_str());
    tex.SetTextColor(MyUtl::C02);
    tex.DrawLatex(xHigh[i] + offX, yHigh[i] + s * offY, labels[i].c_str());
  }

  LegendBox lb = bestLegendCorner(grLow, grHigh, xLo, xHi, yLo, yHi, labelPts);
  TLegend* leg = new TLegend(lb.x1, lb.y1, lb.x2, lb.y2);
  leg->AddEntry(grLow,  Form("nFTrackHS #leq %d", hsSplit), "lp");
  leg->AddEntry(grHigh, Form("nFTrackHS > %d",    hsSplit), "lp");
  leg->Draw("SAME");

  canvas->Print(fname.c_str());
  delete grLow;
  delete grHigh;
}

#endif // SWEEP_UTILITIES_H
