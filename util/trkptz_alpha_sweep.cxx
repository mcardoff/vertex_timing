// trkptz_alpha_sweep.cxx
//
// 1D parameter sweep over the TRKPTZ exponential weighting coefficient α.
// The score is Σpₜ · exp(−α|Δz|); the production code uses α = 1.5 but the
// "natural" choice α = 1 (and the full range 1–5) has never been evaluated.
//
// Clustering is performed once per event at DIST_CUT_CONE (identical to the
// main analysis).  For each α the TRKPTZ score is recomputed from the
// already-stored TRKPT score and a fresh precision-weighted Δz (only two
// branch arrays, trackZ0 and trackVarZ0, are re-read), then the best cluster
// is selected and the standard sweep metrics are accumulated.
//
// Events are split into two regions by the number of forward hard-scatter
// tracks (nFTrackHS from countForwardTracks):
//   Region 0 — Low HS:  nFTrackHS <= HS_TRACK_SPLIT
//   Region 1 — High HS: nFTrackHS >  HS_TRACK_SPLIT
//
// Output PDF: ../parameter_sweep/trkptz_alpha_sweep.pdf
//
// ── Configuration ────────────────────────────────────────────────────────────
//   MAX_EVENTS     – number of events to process (-1 = full sample)
//   ALPHA_MIN      – first α value
//   ALPHA_MAX      – last  α value (inclusive)
//   ALPHA_STEP     – step size between α values
//   ALPHA_NOM      – current production value (drawn as a reference line)
//   HS_TRACK_SPLIT – boundary between the two HS-track regions (inclusive low)
// ─────────────────────────────────────────────────────────────────────────────

#include "sweep_utilities.h"

using namespace MyUtl;

// ── Sweep configuration ───────────────────────────────────────────────────────
static constexpr Long64_t MAX_EVENTS    = -1;   // -1 = full sample
static constexpr double   ALPHA_MIN     =  0.0;
static constexpr double   ALPHA_MAX     =  10.0;
static constexpr double   ALPHA_STEP    =  0.1;
static constexpr double   ALPHA_NOM     =  1.5; // current production value
static constexpr int      HS_TRACK_SPLIT =  5;  // ≤ this → "low HS" region
// ─────────────────────────────────────────────────────────────────────────────

// ── clusterDeltaZ ─────────────────────────────────────────────────────────────
// Precision-weighted z₀ average for the cluster's constituent tracks, minus
// the reco primary vertex z.  Only reads trackZ0 and trackVarZ0.
// ─────────────────────────────────────────────────────────────────────────────
static double clusterDeltaZ(const Cluster& c, const BranchPointerWrapper& b) {
  float znum = 0.f, zden = 0.f;
  for (int trk : c.trackIndices) {
    float varZ = b.trackVarZ0[trk];
    if (varZ <= 0.f) continue;
    znum += b.trackZ0[trk] / varZ;
    zden += 1.f / varZ;
  }
  return (zden > 0.f) ? (double)(znum / zden) - (double)b.recoVtxZ[0] : 0.0;
}

// ── chooseTRKPTZAlpha ─────────────────────────────────────────────────────────
// Select the cluster with the highest Σpₜ · exp(−α|Δz|) for the given α.
// ─────────────────────────────────────────────────────────────────────────────
static Cluster chooseTRKPTZAlpha(const std::vector<Cluster>& clusters,
                                  const BranchPointerWrapper& b, double alpha) {
  const Cluster* best = &clusters[0];
  double bestScore = clusters[0].scores.at(Score::TRKPT.id)
                   * std::exp(-alpha * std::abs(clusterDeltaZ(clusters[0], b)));
  for (const auto& c : clusters) {
    double s = c.scores.at(Score::TRKPT.id)
             * std::exp(-alpha * std::abs(clusterDeltaZ(c, b)));
    if (s > bestScore) { bestScore = s; best = &c; }
  }
  return *best;
}

// ── makeAlphaValues ───────────────────────────────────────────────────────────
static std::vector<double> makeAlphaValues() {
  std::vector<double> alphas;
  int n = (int)std::round((ALPHA_MAX - ALPHA_MIN) / ALPHA_STEP) + 1;
  for (int i = 0; i < n; ++i)
    alphas.push_back(ALPHA_MIN + i * ALPHA_STEP);
  return alphas;
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main() {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

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
  std::cout << "Alpha range       : [" << ALPHA_MIN << ", " << ALPHA_MAX
            << "] step " << ALPHA_STEP
            << "  (nominal = " << ALPHA_NOM << ")\n";

  // --- Build alpha list ---
  const std::vector<double> alphas = makeAlphaValues();
  const int nAlpha = (int)alphas.size();

  // --- Per-alpha, per-region accumulators ---
  const int nResBins = (int)((RES_MAX - RES_MIN) / RES_BWID);

  std::vector<std::vector<TH1D*>>    resHists           (N_REGIONS, std::vector<TH1D*>(nAlpha, nullptr));
  std::vector<std::vector<Long64_t>> nPass              (N_REGIONS, std::vector<Long64_t>(nAlpha, 0));
  std::vector<std::vector<Long64_t>> nTot               (N_REGIONS, std::vector<Long64_t>(nAlpha, 0));
  std::vector<std::vector<double>>   sumNClusters       (N_REGIONS, std::vector<double>(nAlpha, 0.0));
  std::vector<std::vector<Long64_t>> nEventsWithClusters(N_REGIONS, std::vector<Long64_t>(nAlpha, 0));
  std::vector<std::vector<double>>   sumPurity          (N_REGIONS, std::vector<double>(nAlpha, 0.0));
  std::vector<std::vector<double>>   sumTrackCount      (N_REGIONS, std::vector<double>(nAlpha, 0.0));

  const std::string regionLabel[N_REGIONS] = {
    Form("nFTrackHS <= %d", HS_TRACK_SPLIT),
    Form("nFTrackHS > %d",  HS_TRACK_SPLIT)
  };

  for (int r = 0; r < N_REGIONS; ++r) {
    for (int ia = 0; ia < nAlpha; ++ia) {
      resHists[r][ia] = new TH1D(
        Form("res_r%d_a%.2f", r, alphas[ia]),
        Form("TRKPTZ #Deltat (%s, #alpha=%.2f);#Delta t [ps];Entries",
             regionLabel[r].c_str(), alphas[ia]),
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

    if (tracks.empty()) continue;

    // ── Determine HS region ──────────────────────────────────────────────────
    int nFTrack = 0, nFTrackHS = 0, nFTrackPU = 0;
    branch.countForwardTracks(nFTrack, nFTrackHS, nFTrackPU, tracks, /*checkTimeValid=*/true);
    const int region = (nFTrackHS <= HS_TRACK_SPLIT) ? 0 : 1;

    // ── Cluster once per event ────────────────────────────────────────────────
    // Identical to the main analysis: iterative at DIST_CUT_CONE.
    std::vector<Cluster> clusters =
      clusterTracksInTime(
        tracks, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false,
        /*checkTimeValid=*/true,
        /*smearRes=*/IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE,
        /*usez0=*/false,
        /*sortTracks=*/false,
        /*calcPurityFlag=*/true);

    if (clusters.empty()) continue;

    const int nClusters = (int)clusters.size();

    // ── Per-alpha selection and fill ─────────────────────────────────────────
    for (int ia = 0; ia < nAlpha; ++ia) {
      Cluster best = chooseTRKPTZAlpha(clusters, branch, alphas[ia]);

      double diff = best.values.at(0) - branch.truthVtxTime[0];
      resHists[region][ia]->Fill(diff);

      nTot[region][ia]++;
      if (best.passEfficiency(&branch))
        nPass[region][ia]++;

      sumNClusters       [region][ia] += (double)nClusters;
      nEventsWithClusters[region][ia]++;
      sumPurity          [region][ia] += best.purity;
      sumTrackCount      [region][ia] += (double)best.trackIndices.size();
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
    grEff       [r] = new TGraphErrors(nAlpha);
    grSigma     [r] = new TGraphErrors(nAlpha);
    grNClusters [r] = new TGraphErrors(nAlpha);
    grPurity    [r] = new TGraphErrors(nAlpha);
    grTrackCount[r] = new TGraphErrors(nAlpha);

    for (int ia = 0; ia < nAlpha; ++ia) {
      double a = alphas[ia];

      double eff    = (nTot[r][ia] > 0) ? (double)nPass[r][ia] / (double)nTot[r][ia] : 0.0;
      double effErr = (nTot[r][ia] > 0)
        ? std::sqrt(eff * (1.0 - eff) / (double)nTot[r][ia]) : 0.0;
      grEff[r]->SetPoint(ia, a, eff);
      grEff[r]->SetPointError(ia, 0.0, effErr);

      auto [sigma, sigmaErr] = fitCoreSigma(resHists[r][ia], r * nAlpha + ia);
      grSigma[r]->SetPoint(ia, a, sigma);
      grSigma[r]->SetPointError(ia, 0.0, sigmaErr);

      double meanN = (nEventsWithClusters[r][ia] > 0)
        ? sumNClusters[r][ia] / (double)nEventsWithClusters[r][ia] : 0.0;
      grNClusters[r]->SetPoint(ia, a, meanN);
      grNClusters[r]->SetPointError(ia, 0.0, 0.0);

      double meanPurity = (nEventsWithClusters[r][ia] > 0)
        ? sumPurity[r][ia] / (double)nEventsWithClusters[r][ia] : 0.0;
      grPurity[r]->SetPoint(ia, a, meanPurity);
      grPurity[r]->SetPointError(ia, 0.0, 0.0);

      double meanTrackCount = (nEventsWithClusters[r][ia] > 0)
        ? sumTrackCount[r][ia] / (double)nEventsWithClusters[r][ia] : 0.0;
      grTrackCount[r]->SetPoint(ia, a, meanTrackCount);
      grTrackCount[r]->SetPointError(ia, 0.0, 0.0);
    }
  }

  // --- Plotting ---
  TCanvas* canvas = new TCanvas("canvas", "Sweep Results", 800, 600);
  canvas->SetLeftMargin(0.15);

  const std::string outFile = std::string(OUT_DIR) + "/trkptz_alpha_sweep.pdf";
  canvas->Print((outFile + "[").c_str());

  const double xLo    = ALPHA_MIN - ALPHA_STEP * 0.5;
  const double xHi    = ALPHA_MAX + ALPHA_STEP * 0.5;
  const char*  xtitle = "TRKPTZ Weighting Coefficient #alpha";

  // Plot 1: Efficiency vs alpha
  saveGraphPair(outFile, canvas, grEff[0], grEff[1],
    "TRKPTZ Efficiency vs. #alpha",
    xtitle, "Overall Efficiency",
    xLo, xHi, 0.6, 1.05, HS_TRACK_SPLIT,
    0.99, "99% Efficiency",
    ALPHA_NOM, "Current (#alpha = 1.5)");

  // Plot 2: Core sigma vs alpha
  saveGraphPair(outFile, canvas, grSigma[0], grSigma[1],
    "TRKPTZ Core #sigma vs. #alpha",
    xtitle, "Core #sigma [ps]",
    xLo, xHi, 0.0, 60.0, HS_TRACK_SPLIT,
    15.0, "15 ps",
    ALPHA_NOM, "Current (#alpha = 1.5)");

  // Plot 3: Mean cluster count vs alpha
  double nclYMax = 1.0;
  for (int r = 0; r < N_REGIONS; ++r)
    for (int ia = 0; ia < nAlpha; ++ia) {
      double xv, yv; grNClusters[r]->GetPoint(ia, xv, yv);
      if (yv > nclYMax) nclYMax = yv;
    }
  saveGraphPair(outFile, canvas, grNClusters[0], grNClusters[1],
    "Mean Cluster Count vs. #alpha",
    xtitle, "Mean N_{clusters} / event",
    xLo, xHi, 0.0, nclYMax * 1.5, HS_TRACK_SPLIT,
    -1.0, "",
    ALPHA_NOM, "Current (#alpha = 1.5)");

  // Plot 4: Mean purity vs alpha
  saveGraphPair(outFile, canvas, grPurity[0], grPurity[1],
    "Mean Selected-Cluster Purity vs. #alpha",
    xtitle, "Mean Purity (HS p_{T} fraction)",
    xLo, xHi, 0.5, 1.1, HS_TRACK_SPLIT,
    0.75, "75% purity",
    ALPHA_NOM, "Current (#alpha = 1.5)");

  // Plot 5: Mean track count vs alpha
  double tcYMax = 1.0;
  for (int r = 0; r < N_REGIONS; ++r)
    for (int ia = 0; ia < nAlpha; ++ia) {
      double xv, yv; grTrackCount[r]->GetPoint(ia, xv, yv);
      if (yv > tcYMax) tcYMax = yv;
    }
  saveGraphPair(outFile, canvas, grTrackCount[0], grTrackCount[1],
    "Mean Selected-Cluster Track Count vs. #alpha",
    xtitle, "Mean N_{tracks} in selected cluster",
    xLo, xHi, 0.0, tcYMax * 1.2, HS_TRACK_SPLIT,
    -1.0, "",
    ALPHA_NOM, "Current (#alpha = 1.5)");

  // Plot 6: Normalised residual shape overlay per region (sample every 5 steps)
  for (int r = 0; r < N_REGIONS; ++r) {
    std::vector<TH1D*> resClones;
    std::vector<std::string> resLabels;
    for (int ia = 0; ia < nAlpha; ++ia) {
      if (ia % 5 != 0 && ia != nAlpha - 1) continue;
      TH1D* h = (TH1D*)resHists[r][ia]->Clone(Form("res_clone_r%d_%d", r, ia));
      if (h->Integral() == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[resClones.size() % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-400, 400);
      h->SetTitle(Form("TRKPTZ #Deltat shape (%s);#Delta t [ps];Normalised Entries",
                       regionLabel[r].c_str()));
      resLabels.push_back(Form("#alpha = %.1f", alphas[ia]));
      resClones.push_back(h);
    }
    if (resClones.empty()) continue;
    double yMax = 0.0, yMin = 1e30;
    for (auto* h : resClones)
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double v = h->GetBinContent(b);
        if (v > 0.0) { yMax = std::max(yMax, v); yMin = std::min(yMin, v); }
      }
    yMax *= 3.0;  yMin = std::max(yMin * 0.3, 1e-6);
    resClones[0]->GetYaxis()->SetRangeUser(yMin, yMax);
    canvas->SetLogy(true);
    LegendBox lb = bestLegendCornerHist(
      resClones, -400.0, 400.0, yMin, yMax, (int)resClones.size() + 1, true);
    TLegend* leg = new TLegend(lb.x1, lb.y1, lb.x2, lb.y2);
    leg->SetTextSize(0.025);
    leg->SetHeader(regionLabel[r].c_str(), "C");
    bool first = true;
    for (size_t i = 0; i < resClones.size(); ++i) {
      if (first) { resClones[i]->Draw("HIST"); first = false; }
      else          resClones[i]->Draw("HIST SAME");
      leg->AddEntry(resClones[i], resLabels[i].c_str(), "l");
    }
    leg->Draw("SAME");
    canvas->Print(outFile.c_str());
    canvas->SetLogy(false);
    for (auto* h : resClones) delete h;
  }

  // --- Pseudo-ROC plots ---
  std::vector<std::string> alphaLabels(nAlpha);
  for (int ia = 0; ia < nAlpha; ++ia)
    alphaLabels[ia] = Form("%.1f", alphas[ia]);

  const auto effLow  = graphYVals(grEff[0]);
  const auto effHigh = graphYVals(grEff[1]);
  const auto sigLow  = graphYVals(grSigma[0]);
  const auto sigHigh = graphYVals(grSigma[1]);
  const auto purLow  = graphYVals(grPurity[0]);
  const auto purHigh = graphYVals(grPurity[1]);

  saveRocPlot(outFile, canvas,
    effLow, purLow, effHigh, purHigh, alphaLabels,
    "Purity vs. Efficiency (#alpha sweep)",
    "Overall Efficiency", "Mean Cluster Purity",
    HS_TRACK_SPLIT);

  saveRocPlot(outFile, canvas,
    sigLow, purLow, sigHigh, purHigh, alphaLabels,
    "Purity vs. Core #sigma (#alpha sweep)",
    "Core #sigma [ps]", "Mean Cluster Purity",
    HS_TRACK_SPLIT);

  saveRocPlot(outFile, canvas,
    sigLow, effLow, sigHigh, effHigh, alphaLabels,
    "Efficiency vs. Core #sigma (#alpha sweep)",
    "Core #sigma [ps]", "Overall Efficiency",
    HS_TRACK_SPLIT);

  canvas->Print((outFile + "]").c_str());
  std::cout << "Saved: " << outFile << '\n';

  // --- Summary tables ---
  for (int r = 0; r < N_REGIONS; ++r) {
    std::cout << "\n--- Summary Table [" << regionLabel[r] << "] ---\n";
    std::cout << Form("%-8s  %-10s  %-10s  %-12s  %-10s  %-12s\n",
                      "alpha", "Efficiency", "CoreSigma", "MeanNClusters", "MeanPurity", "MeanNTracks");
    for (int ia = 0; ia < nAlpha; ++ia) {
      double a, eff, sig, ncl, pur, ntk, dummy;
      grEff       [r]->GetPoint(ia, a,     eff);
      grSigma     [r]->GetPoint(ia, dummy, sig);
      grNClusters [r]->GetPoint(ia, dummy, ncl);
      grPurity    [r]->GetPoint(ia, dummy, pur);
      grTrackCount[r]->GetPoint(ia, dummy, ntk);
      (void)dummy;
      // mark the nominal production value
      const char* tag = (std::abs(a - ALPHA_NOM) < ALPHA_STEP * 0.5) ? " <-- current" : "";
      std::cout << Form("%-8.2f  %-10.4f  %-10.2f  %-12.2f  %-10.4f  %-12.2f%s\n",
                        a, eff, sig, ncl, pur, ntk, tag);
    }
  }

  return 0;
}
