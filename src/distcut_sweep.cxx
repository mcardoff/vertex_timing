// distcut_sweep.cxx
//
// 1D parameter sweep over the cone clustering distance cut.
// For each value of distCut in [DIST_MIN, DIST_MAX] (step DIST_STEP) we:
//   - Run the full track selection and TRKPTZ cone clustering
//   - Evaluate overall efficiency (passEfficiency) and collect timing residuals
//   - Count the mean number of final clusters per event
//
// After the event loop three summary plots are saved to ../parameter_sweep/:
//   1. Overall efficiency vs. distCut
//   2. Core sigma (double-Gaussian fit, bkg fixed at 175 ps) vs. distCut
//   3. Mean number of clusters per event vs. distCut
//
// ── Configuration ────────────────────────────────────────────────────────────
//   MAX_EVENTS   – number of events to process (set to -1 for the full sample)
//   DIST_MIN     – first distance cut value (exclusive: algorithm needs > 0)
//   DIST_MAX     – last  distance cut value (inclusive)
//   DIST_STEP    – step size between cut values
// ─────────────────────────────────────────────────────────────────────────────

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

#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "suppress_stdout.h"

using namespace MyUtl;

// ── Sweep configuration ───────────────────────────────────────────────────────
static constexpr Long64_t MAX_EVENTS = 10000; // -1 = full sample
static constexpr double   DIST_MIN   = 0.2;
static constexpr double   DIST_MAX   = 4.0;
static constexpr double   DIST_STEP  = 0.2;
// ─────────────────────────────────────────────────────────────────────────────

// Residual histogram range and bin width (ps)
static constexpr double RES_MIN   = -1000.0;
static constexpr double RES_MAX   =  1000.0;
static constexpr double RES_BWID  =    2.0;

// Output directory (must exist — created by CMake / Makefile if needed)
static const char* OUT_DIR = "../parameter_sweep";

// ── Helper: build the list of cut values ─────────────────────────────────────
std::vector<double> makeCutValues() {
  std::vector<double> cuts;
  // Round to avoid floating-point drift
  int n = (int)std::round((DIST_MAX - DIST_MIN) / DIST_STEP) + 1;
  for (int i = 0; i < n; ++i)
    cuts.push_back(DIST_MIN + i * DIST_STEP);
  return cuts;
}

// ── Helper: save a TGraphErrors with cosmetics ────────────────────────────────
void saveGraph(
  const char*      fname,
  TCanvas*         canvas,
  TGraphErrors*    gr,
  const char*      title,
  const char*      ytitle,
  double           yMin,
  double           yMax,
  Color_t          color,
  double           refVal   = -1,
  const char*      refLabel = nullptr
) {
  canvas->cd();

  gr->SetTitle(Form("%s;Distance Cut [#sigma];%s", title, ytitle));
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.9);
  gr->SetLineWidth(2);

  gr->Draw("APE1");
  gr->GetXaxis()->SetLimits(DIST_MIN - DIST_STEP * 0.5, DIST_MAX + DIST_STEP * 0.5);
  gr->GetYaxis()->SetRangeUser(yMin, yMax);
  gr->GetXaxis()->SetNdivisions(510);

  if (refVal >= 0 && refLabel) {
    TLine* refLine = new TLine(
      DIST_MIN - DIST_STEP * 0.5, refVal,
      DIST_MAX + DIST_STEP * 0.5, refVal);
    refLine->SetLineColor(kRed);
    refLine->SetLineWidth(2);
    refLine->SetLineStyle(4);
    refLine->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.15, 0.88, 0.30);
    leg->AddEntry(gr,       "TRKPTZ",  "lp");
    leg->AddEntry(refLine,  refLabel,  "l" );
    leg->Draw("SAME");
  }

  canvas->Print(fname);
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main() {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kFatal;

  // --- Data ---
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  ROOT::EnableImplicitMT();

  const Long64_t nTotal = chain.GetEntries();
  const Long64_t nProcess = (MAX_EVENTS < 0) ? nTotal
                                              : std::min(MAX_EVENTS, nTotal);

  std::cout << "Events to process : " << nProcess << " / " << nTotal << '\n';

  // --- Build cut list ---
  const std::vector<double> cuts = makeCutValues();
  const int nCuts = (int)cuts.size();

  // --- Per-cut accumulators ---
  const int nResBins = (int)((RES_MAX - RES_MIN) / RES_BWID);

  // Residual histograms (one per cut value)
  std::vector<TH1D*> resHists(nCuts);
  // Pass / total counts for efficiency
  std::vector<Long64_t> nPass(nCuts, 0);
  std::vector<Long64_t> nTot (nCuts, 0);
  // Cluster count accumulators
  std::vector<double>   sumNClusters(nCuts, 0.0);
  std::vector<Long64_t> nEventsWithClusters(nCuts, 0);

  for (int ic = 0; ic < nCuts; ++ic) {
    resHists[ic] = new TH1D(
      Form("res_cut%.2f", cuts[ic]),
      Form("TRKPTZ #Deltat (distCut=%.2f);#Delta t [ps];Entries", cuts[ic]),
      nResBins, RES_MIN, RES_MAX);
  }

  // --- Event loop ---
  std::cout << "Starting sweep event loop\n";
  Long64_t evProcessed = 0;

  while (reader.Next()) {
    const Long64_t readNum = chain.GetReadEntry() + 1;

    if (readNum % 200 == 0)
      std::cout << "Progress: " << readNum << "/" << nProcess << "\r" << std::flush;

    if (readNum > nProcess) break;

    // ── Event selection (same as main program) ──────────────────────────────
    if (!branch.passBasicCuts()) continue;
    if (!branch.passJetPtCut())  continue;

    // Track collection (looser sigma for counting, then tightened)
    std::vector<int> tracks =
      getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);

    if (MAX_NSIGMA != 3.0)
      tracks.erase(
        std::remove_if(tracks.begin(), tracks.end(), [&](int trk) {
          return !passTrackVertexAssociation(trk, 0, &branch, MAX_NSIGMA);
        }),
        tracks.end());

    // ── Per-cut clustering ───────────────────────────────────────────────────
    for (int ic = 0; ic < nCuts; ++ic) {
      // HGTD times, cone clustering, no z0 (1D time only, matches main program)
      std::vector<Cluster> clusters =
        clusterTracksInTime(
          tracks, &branch, cuts[ic],
          /*useSmearedTimes=*/false,
          /*checkTimeValid=*/true,
          /*smearRes=*/10.0,
          /*useCone=*/true,
          /*usez0=*/false,
          /*calcPurityFlag=*/false);

      if (clusters.empty()) continue;

      // Select best cluster by TRKPTZ
      Cluster best = chooseCluster(clusters, Score::TRKPTZ);

      // Fill residual histogram
      double diff = best.values.at(0) - branch.truthVtxTime[0];
      resHists[ic]->Fill(diff);

      // Efficiency
      nTot[ic]++;
      if (best.passEfficiency(&branch))
        nPass[ic]++;

      // Cluster count
      sumNClusters[ic]  += (double)clusters.size();
      nEventsWithClusters[ic]++;
    }

    evProcessed++;
  }

  std::cout << "\nFinished event loop (" << evProcessed << " events processed)\n";
  ROOT::DisableImplicitMT();

  // --- Build summary TGraphErrors ---
  TGraphErrors* grEff    = new TGraphErrors(nCuts);
  TGraphErrors* grSigma  = new TGraphErrors(nCuts);
  TGraphErrors* grNClusters = new TGraphErrors(nCuts);

  TLatex latex;
  latex.SetTextSize(0.035);
  latex.SetTextAlign(13);

  // Fit each residual histogram and fill graphs
  // Suppress ROOT's OBJ printf during fitting
  {
    SuppressStdout suppress;

    for (int ic = 0; ic < nCuts; ++ic) {
      double dc = cuts[ic];

      // ── Efficiency ─────────────────────────────────────────────────────────
      double eff    = (nTot[ic] > 0) ? (double)nPass[ic] / (double)nTot[ic] : 0.0;
      // Wilson score interval (standard binomial CI for efficiency)
      double effErr = (nTot[ic] > 0)
        ? std::sqrt(eff * (1.0 - eff) / (double)nTot[ic])
        : 0.0;
      grEff->SetPoint(ic, dc, eff);
      grEff->SetPointError(ic, 0.0, effErr);

      // ── Core sigma from double-Gaussian fit ────────────────────────────────
      double sigma    = 0.0;
      double sigmaErr = 0.0;
      if (resHists[ic]->GetEntries() >= 20) {
        // Reuse createDblFit with fixbkg=true (bkg Gaussian fixed at 175 ps)
        // Mean fixed at 0, core sigma is free, bkg sigma fixed at PILEUP_SMEAR
        TF1* fit = new TF1(
          Form("dgaus_fit_%d", ic),
          "[1]*TMath::Exp(-0.5*((x-[0])/[3])^2) + [2]*TMath::Exp(-0.5*((x-[0])/[4])^2)",
          RES_MIN, RES_MAX);
        fit->SetParameters(
          0,                                    // [0] mean — fixed at 0
          0.8  * resHists[ic]->GetMaximum(),   // [1] core norm
          0.01 * resHists[ic]->GetMaximum(),   // [2] bkg norm
          26.0,                                 // [3] core sigma (free)
          PILEUP_SMEAR                          // [4] bkg sigma — fixed
        );
        fit->FixParameter(0, 0.0);
        fit->SetParLimits(1, 0, 1e6);
        fit->SetParLimits(2, 0, 1e6);
        fit->SetParLimits(3, 0.1, PILEUP_SMEAR);
        fit->FixParameter(4, PILEUP_SMEAR);
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->SetNpx(1000);
        resHists[ic]->Fit(fit, "RQ");

        // Parameter 3 is always the core sigma (free), parameter 4 is bkg (fixed)
        sigma    = fit->GetParameter(3);
        sigmaErr = fit->GetParError(3);
        delete fit;
      }
      grSigma->SetPoint(ic, dc, sigma);
      grSigma->SetPointError(ic, 0.0, sigmaErr);

      // ── Mean cluster count ─────────────────────────────────────────────────
      double meanN    = (nEventsWithClusters[ic] > 0)
        ? sumNClusters[ic] / (double)nEventsWithClusters[ic]
        : 0.0;
      double meanNErr = 0.0; // statistical uncertainty negligible for mean
      grNClusters->SetPoint(ic, dc, meanN);
      grNClusters->SetPointError(ic, 0.0, meanNErr);
    }
  } // SuppressStdout scope

  // --- Plotting ---
  TCanvas* canvas = new TCanvas("canvas", "Sweep Results", 800, 600);
  canvas->SetLeftMargin(0.15);

  const char* outFile = Form("%s/distcut_sweep.pdf", OUT_DIR);
  canvas->Print(Form("%s[", outFile));

  // Plot 1: Efficiency vs distCut
  saveGraph(
    outFile, canvas, grEff,
    "TRKPTZ Efficiency vs. Distance Cut",
    "Overall Efficiency",
    0.0, 1.15,
    C01,
    0.99, "99% Efficiency");

  // Plot 2: Core sigma vs distCut
  saveGraph(
    outFile, canvas, grSigma,
    "TRKPTZ Core #sigma vs. Distance Cut",
    "Core #sigma [ps]",
    0.0, 60.0,
    C02,
    15.0, "15 ps");

  // Plot 3: Mean cluster count vs distCut
  saveGraph(
    outFile, canvas, grNClusters,
    "Mean Cluster Count vs. Distance Cut",
    "Mean N_{clusters} / event",
    0.0, grNClusters->GetYaxis()->GetXmax() * 1.2,
    C03);

  // Bonus: draw all residual histograms (normalised) on one page to show shape evolution
  {
    canvas->SetLogy(true);
    TLegend* leg = new TLegend(0.65, 0.60, 0.90, 0.90);
    leg->SetTextSize(0.025);
    bool first = true;
    // Only draw every other cut to avoid clutter
    for (int ic = 0; ic < nCuts; ++ic) {
      if (ic % 2 != 0 && ic != nCuts - 1) continue;
      TH1D* h = (TH1D*)resHists[ic]->Clone(Form("res_clone_%d", ic));
      if (h->Integral() == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[ic % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-400, 400);
      h->GetYaxis()->SetTitle("Normalised Entries");
      h->SetTitle("TRKPTZ #Deltat shape vs. Distance Cut;#Delta t [ps];Normalised Entries");
      if (first) { h->Draw("HIST"); first = false; }
      else        h->Draw("HIST SAME");
      leg->AddEntry(h, Form("distCut=%.1f", cuts[ic]), "l");
    }
    leg->Draw("SAME");
    canvas->Print(outFile);
    canvas->SetLogy(false);
  }

  canvas->Print(Form("%s]", outFile));

  std::cout << "Saved: " << outFile << '\n';

  // Print a quick summary table to stdout
  std::cout << "\n--- Summary Table ---\n";
  std::cout << Form("%-10s  %-10s  %-10s  %-12s\n",
                    "distCut", "Efficiency", "CoreSigma", "MeanNClusters");
  for (int ic = 0; ic < nCuts; ++ic) {
    double dc, eff, sig, ncl, dummy;
    grEff->GetPoint(ic, dc, eff);
    grSigma->GetPoint(ic, dc, sig);
    grNClusters->GetPoint(ic, dc, ncl);
    (void)dummy;
    std::cout << Form("%-10.2f  %-10.4f  %-10.2f  %-12.2f\n",
                      dc, eff, sig, ncl);
  }

  return 0;
}
