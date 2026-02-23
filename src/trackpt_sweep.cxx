// trackpt_sweep.cxx
//
// 1D parameter sweep over the maximum track pT cut applied during track
// selection.  For each value of maxTrackPt in [PT_MIN, PT_MAX] (step PT_STEP)
// we:
//   - Run the full track selection (using that pT ceiling) and TRKPTZ cone
//     clustering with the nominal distance cut
//   - Evaluate overall efficiency (passEfficiency) and collect timing residuals
//   - Count the mean number of final clusters per event
//
// After the event loop three summary plots are saved to ../parameter_sweep/:
//   1. Overall efficiency vs. maxTrackPt
//   2. Core sigma (double-Gaussian fit, bkg fixed at 175 ps) vs. maxTrackPt
//   3. Mean number of clusters per event vs. maxTrackPt
//   4. Normalised residual shape overlay (every other cut value)
//
// ── Configuration ────────────────────────────────────────────────────────────
//   MAX_EVENTS  – number of events to process (set to -1 for the full sample)
//   PT_MIN      – first max-pT cut value (GeV)
//   PT_MAX      – last  max-pT cut value (GeV, inclusive)
//   PT_STEP     – step size between cut values (GeV)
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

using namespace MyUtl;

// ── Sweep configuration ───────────────────────────────────────────────────────
static constexpr Long64_t MAX_EVENTS = -1;   // -1 = full sample
static constexpr double   PT_MIN     =  2.0; // GeV
static constexpr double   PT_MAX     = 50.0; // GeV
static constexpr double   PT_STEP    =  1.0; // GeV
// ─────────────────────────────────────────────────────────────────────────────

// Residual histogram range and bin width (ps)
static constexpr double RES_MIN  = -1000.0;
static constexpr double RES_MAX  =  1000.0;
static constexpr double RES_BWID =     2.0;

// Output directory (must exist — created by CMake / Makefile if needed)
static const char* OUT_DIR = "../parameter_sweep";

// ── Helper: build the list of cut values ─────────────────────────────────────
std::vector<double> makeCutValues() {
  std::vector<double> cuts;
  // Round to avoid floating-point drift
  int n = (int)std::round((PT_MAX - PT_MIN) / PT_STEP) + 1;
  for (int i = 0; i < n; ++i)
    cuts.push_back(PT_MIN + i * PT_STEP);
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
  
  gr->SetTitle(Form("%s;Max Track p_{T} [GeV];%s", title, ytitle));
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.9);
  gr->SetLineWidth(2);

  gr->Draw("APE1");
  gr->GetXaxis()->SetLimits(PT_MIN - PT_STEP * 0.5, PT_MAX + PT_STEP * 0.5);
  gr->GetYaxis()->SetRangeUser(yMin, yMax);
  gr->GetXaxis()->SetNdivisions(510);

  if (refVal >= 0 && refLabel) {
    TLine* refLine = new TLine(
      PT_MIN - PT_STEP * 0.5, refVal,
      PT_MAX + PT_STEP * 0.5, refVal);
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

  const Long64_t nTotal   = chain.GetEntries();
  const Long64_t nProcess = (MAX_EVENTS < 0) ? nTotal
                                              : std::min(MAX_EVENTS, nTotal);

  std::cout << "Events to process : " << nProcess << " / " << nTotal << '\n';

  // --- Build cut list ---
  const std::vector<double> cuts = makeCutValues();
  const int nCuts = (int)cuts.size();

  // --- Per-cut accumulators ---
  const int nResBins = (int)((RES_MAX - RES_MIN) / RES_BWID);

  std::vector<TH1D*>     resHists          (nCuts, nullptr);
  std::vector<Long64_t>  nPass             (nCuts, 0);
  std::vector<Long64_t>  nTot              (nCuts, 0);
  std::vector<double>    sumNClusters      (nCuts, 0.0);
  std::vector<Long64_t>  nEventsWithClusters(nCuts, 0);

  for (int ic = 0; ic < nCuts; ++ic) {
    resHists[ic] = new TH1D(
      Form("res_ptmax%.1f", cuts[ic]),
      Form("TRKPTZ #Deltat (maxTrkPt=%.1f GeV);#Delta t [ps];Entries", cuts[ic]),
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

    // ── Per-cut clustering ───────────────────────────────────────────────────
    for (int ic = 0; ic < nCuts; ++ic) {
      const double maxPt = cuts[ic];

      // Collect tracks at looser 3σ with this pT ceiling
      std::vector<int> tracks =
        getAssociatedTracks(&branch, MIN_TRACK_PT, maxPt, 3.0);

      // Tighten to MAX_NSIGMA if needed
      if (MAX_NSIGMA != 3.0)
        tracks.erase(
          std::remove_if(tracks.begin(), tracks.end(), [&](int trk) {
            return !passTrackVertexAssociation(trk, 0, &branch, MAX_NSIGMA);
          }),
          tracks.end());

      // HGTD times, cone clustering with the nominal distance cut (3σ)
      std::vector<Cluster> clusters =
        clusterTracksInTime(
          tracks, &branch, 3.0,
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
      sumNClusters[ic]   += (double)clusters.size();
      nEventsWithClusters[ic]++;
    }

    evProcessed++;
  }

  std::cout << "\nFinished event loop (" << evProcessed << " events processed)\n";
  ROOT::DisableImplicitMT();

  // --- Build summary TGraphErrors ---
  TGraphErrors* grEff       = new TGraphErrors(nCuts);
  TGraphErrors* grSigma     = new TGraphErrors(nCuts);
  TGraphErrors* grNClusters = new TGraphErrors(nCuts);

  TLatex latex;
  latex.SetTextSize(0.035);
  latex.SetTextAlign(13);

  // Fit each residual histogram and fill graphs
  for (int ic = 0; ic < nCuts; ++ic) {
    double pt = cuts[ic];

    // ── Efficiency ─────────────────────────────────────────────────────────
    double eff    = (nTot[ic] > 0) ? (double)nPass[ic] / (double)nTot[ic] : 0.0;
    double effErr = (nTot[ic] > 0)
      ? std::sqrt(eff * (1.0 - eff) / (double)nTot[ic])
      : 0.0;
    grEff->SetPoint(ic, pt, eff);
    grEff->SetPointError(ic, 0.0, effErr);

    // ── Core sigma from double-Gaussian fit ────────────────────────────────
    double sigma    = 0.0;
    double sigmaErr = 0.0;
    if (resHists[ic]->GetEntries() >= 20) {
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

      sigma    = fit->GetParameter(3);
      sigmaErr = fit->GetParError(3);
      delete fit;
    }
    grSigma->SetPoint(ic, pt, sigma);
    grSigma->SetPointError(ic, 0.0, sigmaErr);

    // ── Mean cluster count ─────────────────────────────────────────────────
    double meanN = (nEventsWithClusters[ic] > 0)
      ? sumNClusters[ic] / (double)nEventsWithClusters[ic]
      : 0.0;
    grNClusters->SetPoint(ic, pt, meanN);
    grNClusters->SetPointError(ic, 0.0, 0.0);
  }

  // --- Plotting ---
  TCanvas *canvas = new TCanvas("canvas", "pT Sweep Results", 800, 600);
  canvas->SetLogy(true);  
  canvas->SetLeftMargin(0.15);

  const char* outFile = Form("%s/trackpt_sweep.pdf", OUT_DIR);
  canvas->Print(Form("%s[", outFile));

  // Plot 1: Efficiency vs maxTrackPt
  saveGraph(
    outFile, canvas, grEff,
    "TRKPTZ Efficiency vs. Max Track p_{T}",
    "Overall Efficiency",
    0.0, 1.15,
    C01,
    0.99, "99% Efficiency");

  // Plot 2: Core sigma vs maxTrackPt
  saveGraph(
    outFile, canvas, grSigma,
    "TRKPTZ Core #sigma vs. Max Track p_{T}",
    "Core #sigma [ps]",
    0.0, 60.0,
    C02,
    15.0, "15 ps");

  // Plot 3: Mean cluster count vs maxTrackPt
  saveGraph(
    outFile, canvas, grNClusters,
    "Mean Cluster Count vs. Max Track p_{T}",
    "Mean N_{clusters} / event",
    0.0, grNClusters->GetYaxis()->GetXmax() * 1.2,
    C03);

  // Plot 4: Normalised residual shape overlay (every other cut to avoid clutter)
  {
    canvas->SetLogy(true);
    TLegend* leg = new TLegend(0.65, 0.60, 0.90, 0.90);
    leg->SetTextSize(0.025);
    bool first = true;
    for (int ic = 0; ic < nCuts; ++ic) {
      if (ic % 2 != 0 && ic != nCuts - 1) continue;
      TH1D* h = (TH1D*)resHists[ic]->Clone(Form("res_clone_%d", ic));
      if (h->Integral() == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[ic % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-400, 400);
      h->GetYaxis()->SetTitle("Normalised Entries");
      h->SetTitle("TRKPTZ #Deltat shape vs. Max Track p_{T};#Delta t [ps];Normalised Entries");
      if (first) { h->Draw("HIST"); first = false; }
      else        h->Draw("HIST SAME");
      leg->AddEntry(h, Form("maxPt=%.0f GeV", cuts[ic]), "l");
    }
    leg->Draw("SAME");
    canvas->Print(outFile);
    canvas->SetLogy(false);
  }

  canvas->Print(Form("%s]", outFile));

  std::cout << "Saved: " << outFile << '\n';

  // Print a quick summary table to stdout
  std::cout << "\n--- Summary Table ---\n";
  std::cout << Form("%-12s  %-10s  %-10s  %-14s\n",
                    "maxTrackPt", "Efficiency", "CoreSigma", "MeanNClusters");
  for (int ic = 0; ic < nCuts; ++ic) {
    double pt, eff, sig, ncl, dummy;
    grEff->GetPoint(ic, pt, eff);
    grSigma->GetPoint(ic, dummy, sig);
    grNClusters->GetPoint(ic, dummy, ncl);
    (void)dummy;
    std::cout << Form("%-12.1f  %-10.4f  %-10.2f  %-14.2f\n",
                      pt, eff, sig, ncl);
  }

  return 0;
}
