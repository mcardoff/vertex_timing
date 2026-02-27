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
  double        refVal   = -1,
  const char*   refLabel = nullptr)
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
  leg->Draw("SAME");

  canvas->Print(fname.c_str());
}

#endif // SWEEP_UTILITIES_H
