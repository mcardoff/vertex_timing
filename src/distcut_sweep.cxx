// distcut_sweep.cxx
//
// 1D parameter sweep over the cone clustering distance cut.
// For each value of distCut in [DIST_MIN, DIST_MAX] (step DIST_STEP) we:
//   - Run the full track selection and TRKPTZ cone clustering
//   - Evaluate overall efficiency (passEfficiency) and collect timing residuals
//   - Count the mean number of final clusters per event
//
// Events are split into two regions by the number of forward hard-scatter
// tracks (nFTrackHS from countForwardTracks):
//   Region 0 — Low HS:  nFTrackHS <= HS_TRACK_SPLIT
//   Region 1 — High HS: nFTrackHS >  HS_TRACK_SPLIT
//
// For each metric (efficiency, core sigma, mean cluster count, mean purity,
// mean track count) the two regions are overlaid on the same canvas page.
//
// Output PDF: ../parameter_sweep/distcut_sweep.pdf
//
// ── Configuration ────────────────────────────────────────────────────────────
//   MAX_EVENTS    – number of events to process (-1 = full sample)
//   DIST_MIN      – first distance cut value (> 0)
//   DIST_MAX      – last  distance cut value (inclusive)
//   DIST_STEP     – step size between cut values
//   HS_TRACK_SPLIT – boundary between the two HS-track regions (inclusive low)
// ─────────────────────────────────────────────────────────────────────────────

#include "sweep_utilities.h"

using namespace MyUtl;

// ── Sweep configuration ───────────────────────────────────────────────────────
static constexpr Long64_t MAX_EVENTS    = -1;  // -1 = full sample
static constexpr double   DIST_MIN      = 0.2;
static constexpr double   DIST_MAX      = 4.0;
static constexpr double   DIST_STEP     = 0.1;
static constexpr int      HS_TRACK_SPLIT = 5;  // ≤ this → "low HS" region
// ─────────────────────────────────────────────────────────────────────────────

// ── Helper: build the list of cut values ─────────────────────────────────────
std::vector<double> makeCutValues() {
  std::vector<double> cuts;
  int n = (int)std::round((DIST_MAX - DIST_MIN) / DIST_STEP) + 1;
  for (int i = 0; i < n; ++i)
    cuts.push_back(DIST_MIN + i * DIST_STEP);
  return cuts;
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main() {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  // gErrorIgnoreLevel = kFatal;

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
  std::cout << "HS track split    : nFTrackHS <= " << HS_TRACK_SPLIT
            << "  (low) / > " << HS_TRACK_SPLIT << " (high)\n";

  // --- Build cut list ---
  const std::vector<double> cuts = makeCutValues();
  const int nCuts = (int)cuts.size();

  // --- Per-cut, per-region accumulators ---
  // [region][cut_index]
  const int nResBins = (int)((RES_MAX - RES_MIN) / RES_BWID);

  std::vector<std::vector<TH1D*>>    resHists           (N_REGIONS, std::vector<TH1D*>(nCuts, nullptr));
  std::vector<std::vector<Long64_t>> nPass              (N_REGIONS, std::vector<Long64_t>(nCuts, 0));
  std::vector<std::vector<Long64_t>> nTot               (N_REGIONS, std::vector<Long64_t>(nCuts, 0));
  std::vector<std::vector<double>>   sumNClusters       (N_REGIONS, std::vector<double>(nCuts, 0.0));
  std::vector<std::vector<Long64_t>> nEventsWithClusters(N_REGIONS, std::vector<Long64_t>(nCuts, 0));
  std::vector<std::vector<double>>   sumPurity          (N_REGIONS, std::vector<double>(nCuts, 0.0));
  std::vector<std::vector<double>>   sumTrackCount      (N_REGIONS, std::vector<double>(nCuts, 0.0));

  const std::string regionLabel[N_REGIONS] = {
    Form("nFTrackHS <= %d", HS_TRACK_SPLIT),
    Form("nFTrackHS > %d",  HS_TRACK_SPLIT)
  };

  for (int r = 0; r < N_REGIONS; ++r) {
    for (int ic = 0; ic < nCuts; ++ic) {
      resHists[r][ic] = new TH1D(
        Form("res_r%d_cut%.2f", r, cuts[ic]),
        Form("TRKPTZ #Deltat (%s, distCut=%.2f);#Delta t [ps];Entries",
             regionLabel[r].c_str(), cuts[ic]),
        nResBins, RES_MIN, RES_MAX);
    }
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

    // ── Track collection ─────────────────────────────────────────────────────
    std::vector<int> tracks =
      getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);

    if (MAX_NSIGMA != 3.0)
      tracks.erase(
        std::remove_if(tracks.begin(), tracks.end(), [&](int trk) {
          return !passTrackVertexAssociation(trk, 0, &branch, MAX_NSIGMA);
        }),
        tracks.end());

    // ── Determine HS region ──────────────────────────────────────────────────
    int nFTrack = 0, nFTrackHS = 0, nFTrackPU = 0;
    branch.countForwardTracks(nFTrack, nFTrackHS, nFTrackPU, tracks, /*checkTimeValid=*/true);
    const int region = (nFTrackHS <= HS_TRACK_SPLIT) ? 0 : 1;

    // ── Per-cut clustering ───────────────────────────────────────────────────
    for (int ic = 0; ic < nCuts; ++ic) {
      std::vector<Cluster> clusters =
        clusterTracksInTime(
          tracks, &branch, cuts[ic],
          /*useSmearedTimes=*/false,
          /*checkTimeValid=*/true,
          /*smearRes=*/10.0,
          ClusteringMethod::CONE,
          /*usez0=*/false,
          /*sortTracks=*/false,
          /*calcPurityFlag=*/true);

      if (clusters.empty()) continue;

      Cluster best = chooseCluster(clusters, Score::TRKPTZ);

      double diff = best.values.at(0) - branch.truthVtxTime[0];
      resHists[region][ic]->Fill(diff);

      nTot[region][ic]++;
      if (best.passEfficiency(&branch))
        nPass[region][ic]++;

      sumNClusters[region][ic]        += (double)clusters.size();
      nEventsWithClusters[region][ic]++;
      sumPurity[region][ic]           += best.purity;
      sumTrackCount[region][ic]       += (double)best.trackIndices.size();
    }

    evProcessed++;
  }

  std::cout << "\nFinished event loop (" << evProcessed << " events processed)\n";
  ROOT::DisableImplicitMT();

  // --- Build per-region TGraphErrors and fit ---
  TGraphErrors* grEff        [N_REGIONS];
  TGraphErrors* grSigma      [N_REGIONS];
  TGraphErrors* grNClusters  [N_REGIONS];
  TGraphErrors* grPurity     [N_REGIONS];
  TGraphErrors* grTrackCount [N_REGIONS];

  for (int r = 0; r < N_REGIONS; ++r) {
    grEff       [r] = new TGraphErrors(nCuts);
    grSigma     [r] = new TGraphErrors(nCuts);
    grNClusters [r] = new TGraphErrors(nCuts);
    grPurity    [r] = new TGraphErrors(nCuts);
    grTrackCount[r] = new TGraphErrors(nCuts);

    for (int ic = 0; ic < nCuts; ++ic) {
      double dc = cuts[ic];

      // Efficiency
      double eff    = (nTot[r][ic] > 0) ? (double)nPass[r][ic] / (double)nTot[r][ic] : 0.0;
      double effErr = (nTot[r][ic] > 0)
        ? std::sqrt(eff * (1.0 - eff) / (double)nTot[r][ic]) : 0.0;
      grEff[r]->SetPoint(ic, dc, eff);
      grEff[r]->SetPointError(ic, 0.0, effErr);

      // Core sigma
      auto [sigma, sigmaErr] = fitCoreSigma(resHists[r][ic], r * nCuts + ic);
      grSigma[r]->SetPoint(ic, dc, sigma);
      grSigma[r]->SetPointError(ic, 0.0, sigmaErr);

      // Mean cluster count
      double meanN = (nEventsWithClusters[r][ic] > 0)
        ? sumNClusters[r][ic] / (double)nEventsWithClusters[r][ic] : 0.0;
      grNClusters[r]->SetPoint(ic, dc, meanN);
      grNClusters[r]->SetPointError(ic, 0.0, 0.0);

      // Mean purity of selected cluster
      double meanPurity = (nEventsWithClusters[r][ic] > 0)
        ? sumPurity[r][ic] / (double)nEventsWithClusters[r][ic] : 0.0;
      grPurity[r]->SetPoint(ic, dc, meanPurity);
      grPurity[r]->SetPointError(ic, 0.0, 0.0);

      // Mean track count of selected cluster
      double meanTrackCount = (nEventsWithClusters[r][ic] > 0)
        ? sumTrackCount[r][ic] / (double)nEventsWithClusters[r][ic] : 0.0;
      grTrackCount[r]->SetPoint(ic, dc, meanTrackCount);
      grTrackCount[r]->SetPointError(ic, 0.0, 0.0);
    }
  }

  const int nPlot[N_REGIONS] = {nCuts, nCuts};

  // --- Plotting ---
  TCanvas* canvas = new TCanvas("canvas", "Sweep Results", 800, 600);
  canvas->SetLeftMargin(0.15);

  const std::string outFile = std::string(OUT_DIR) + "/distcut_sweep.pdf";
  canvas->Print((outFile + "[").c_str());

  const double xLo = DIST_MIN - DIST_STEP * 0.5;
  const double xHi = DIST_MAX + DIST_STEP * 0.5;
  const char*  xtitle = "Distance Cut [#sigma]";

  // Plot 1: Efficiency vs distCut (both regions overlaid)
  saveGraphPair(outFile, canvas, grEff[0], grEff[1],
    "TRKPTZ Efficiency vs. Distance Cut",
    xtitle, "Overall Efficiency",
    xLo, xHi, 0.6, 1.3, HS_TRACK_SPLIT,
    0.99, "99% Efficiency");

  // Plot 2: Core sigma vs distCut (both regions overlaid)
  saveGraphPair(outFile, canvas, grSigma[0], grSigma[1],
    "TRKPTZ Core #sigma vs. Distance Cut",
    xtitle, "Core #sigma [ps]",
    xLo, xHi, 0.0, 60.0, HS_TRACK_SPLIT,
    15.0, "15 ps");

  // Plot 3: Mean cluster count vs distCut (both regions overlaid)
  double nclYMax = 1.0;
  for (int r = 0; r < N_REGIONS; ++r) {
    for (int ic = 0; ic < nPlot[r]; ++ic) {
      double xv, yv;
      grNClusters[r]->GetPoint(ic, xv, yv);
      if (yv > nclYMax) nclYMax = yv;
    }
  }
  saveGraphPair(outFile, canvas, grNClusters[0], grNClusters[1],
    "Mean Cluster Count vs. Distance Cut",
    xtitle, "Mean N_{clusters} / event",
    xLo, xHi, 0.0, nclYMax * 1.5, HS_TRACK_SPLIT);

  // Plot 4: Mean purity of selected cluster vs distCut (both regions overlaid)
  saveGraphPair(outFile, canvas, grPurity[0], grPurity[1],
    "Mean Selected-Cluster Purity vs. Distance Cut",
    xtitle, "Mean Purity (HS p_{T} fraction)",
    xLo, xHi, 0.5, 1.1, HS_TRACK_SPLIT,
    0.75, "75% purity");

  // Plot 5: Mean track count of selected cluster vs distCut (both regions overlaid)
  double tcYMax = 1.0;
  for (int r = 0; r < N_REGIONS; ++r) {
    for (int ic = 0; ic < nPlot[r]; ++ic) {
      double xv, yv;
      grTrackCount[r]->GetPoint(ic, xv, yv);
      if (yv > tcYMax) tcYMax = yv;
    }
  }
  saveGraphPair(outFile, canvas, grTrackCount[0], grTrackCount[1],
    "Mean Selected-Cluster Track Count vs. Distance Cut",
    xtitle, "Mean N_{tracks} in selected cluster",
    xLo, xHi, 0.0, tcYMax * 1.2, HS_TRACK_SPLIT);

  // Plot 6: Normalised residual shape overlay per region (two pages)
  for (int r = 0; r < N_REGIONS; ++r) {
    canvas->SetLogy(true);
    TLegend* leg = new TLegend(0.62, 0.58, 0.90, 0.88);
    leg->SetTextSize(0.025);
    leg->SetHeader(regionLabel[r].c_str(), "C");
    bool first = true;
    for (int ic = 0; ic < nPlot[r]; ++ic) {
      if (ic % 2 != 0 && ic != nPlot[r] - 1) continue;
      TH1D* h = (TH1D*)resHists[r][ic]->Clone(Form("res_clone_r%d_%d", r, ic));
      if (h->Integral() == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[ic % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-400, 400);
      h->GetYaxis()->SetTitle("Normalised Entries");
      h->SetTitle(Form("TRKPTZ #Deltat shape (%s);#Delta t [ps];Normalised Entries",
                       regionLabel[r].c_str()));
      if (first) { h->Draw("HIST"); first = false; }
      else        h->Draw("HIST SAME");
      leg->AddEntry(h, Form("distCut=%.1f", cuts[ic]), "l");
    }
    if (!first) {
      leg->Draw("SAME");
      canvas->Print(outFile.c_str());
    }
    canvas->SetLogy(false);
  }

  canvas->Print((outFile + "]").c_str());
  std::cout << "Saved: " << outFile << '\n';

  // --- Summary tables ---
  for (int r = 0; r < N_REGIONS; ++r) {
    std::cout << "\n--- Summary Table [" << regionLabel[r] << "] ---\n";
    std::cout << Form("%-10s  %-10s  %-10s  %-12s  %-10s  %-12s\n",
                      "distCut", "Efficiency", "CoreSigma", "MeanNClusters", "MeanPurity", "MeanNTracks");
    for (int ic = 0; ic < nPlot[r]; ++ic) {
      double dc, eff, sig, ncl, pur, ntk, dummy;
      grEff       [r]->GetPoint(ic, dc,    eff);
      grSigma     [r]->GetPoint(ic, dummy, sig);
      grNClusters [r]->GetPoint(ic, dummy, ncl);
      grPurity    [r]->GetPoint(ic, dummy, pur);
      grTrackCount[r]->GetPoint(ic, dummy, ntk);
      (void)dummy;
      std::cout << Form("%-10.2f  %-10.4f  %-10.2f  %-12.2f  %-10.4f  %-12.2f\n",
                        dc, eff, sig, ncl, pur, ntk);
    }
  }

  return 0;
}
