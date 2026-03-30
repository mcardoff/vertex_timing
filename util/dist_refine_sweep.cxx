// dist_refine_sweep.cxx
//
// 1D parameter sweep over DIST_CUT_REFINE: the iterative clustering distance
// cut used in the REFINED algorithm.  For each value in [REFINE_MIN, REFINE_MAX]
// (step REFINE_STEP) we:
//   - Run cone clustering once per event at DIST_CUT_CONE (k-independent)
//   - Pool the trackIndices of the top CONE_ITER_K qualifying cone clusters
//   - Re-cluster the pooled tracks with ClusteringMethod::ITERATIVE at the
//     sweep distance cut value
//   - Apply MIN_CLUSTER_TRACKS filter; evaluate efficiency and collect residuals
//
// The sweep range should bracket DIST_CUT_CONE (3σ) from below: values above
// the cone cut are redundant since the tracks already came from a 3σ cone.
// Typical range: 0.5σ – 3.0σ in steps of 0.25σ.
//
// Events are split into two regions by the number of forward hard-scatter
// tracks (nFTrackHS):
//   Region 0 — Low HS:  nFTrackHS <= HS_TRACK_SPLIT
//   Region 1 — High HS: nFTrackHS >  HS_TRACK_SPLIT
//
// Output PDF: ../parameter_sweep/dist_refine_sweep.pdf
//
// ── Configuration ─────────────────────────────────────────────────────────────
//   MAX_EVENTS     – number of events to process (-1 = full sample)
//   REFINE_MIN     – first distance cut value (> 0)
//   REFINE_MAX     – last  distance cut value (inclusive)
//   REFINE_STEP    – step size between cut values
//   HS_TRACK_SPLIT – boundary between the two HS-track regions (inclusive low)
// ─────────────────────────────────────────────────────────────────────────────

#include "sweep_utilities.h"

using namespace MyUtl;

// ── Sweep configuration ───────────────────────────────────────────────────────
static constexpr Long64_t MAX_EVENTS     = -1;   // -1 = full sample
static constexpr double   REFINE_MIN     = 0.5;
static constexpr double   REFINE_MAX     = 3.0;
static constexpr double   REFINE_STEP    = 0.25;
static constexpr int      HS_TRACK_SPLIT = 5;    // ≤ this → "low HS" region
// ─────────────────────────────────────────────────────────────────────────────

// ── Helper: build the list of cut values ─────────────────────────────────────
std::vector<double> makeCutValues() {
  std::vector<double> cuts;
  int n = (int)std::round((REFINE_MAX - REFINE_MIN) / REFINE_STEP) + 1;
  for (int i = 0; i < n; ++i)
    cuts.push_back(REFINE_MIN + i * REFINE_STEP);
  return cuts;
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

  const Long64_t nTotal   = chain.GetEntries();
  const Long64_t nProcess = (MAX_EVENTS < 0) ? nTotal
                                              : std::min(MAX_EVENTS, nTotal);

  const std::vector<double> cuts = makeCutValues();
  const int nCuts = (int)cuts.size();
  const int nResBins = (int)((RES_MAX - RES_MIN) / RES_BWID);

  std::cout << "Events to process : " << nProcess << " / " << nTotal << '\n';
  std::cout << "Refinement cut sweep: " << REFINE_MIN << " .. " << REFINE_MAX
            << " (step " << REFINE_STEP << ", " << nCuts << " values)\n";
  std::cout << "CONE_ITER_K       : " << CONE_ITER_K << '\n';
  std::cout << "HS track split    : nFTrackHS <= " << HS_TRACK_SPLIT
            << "  (low) / > " << HS_TRACK_SPLIT << " (high)\n";

  const std::string regionLabel[N_REGIONS] = {
    Form("nFTrackHS <= %d", HS_TRACK_SPLIT),
    Form("nFTrackHS > %d",  HS_TRACK_SPLIT)
  };

  // --- Per-cut, per-region accumulators ---
  // [region][cut_index]
  std::vector<std::vector<TH1D*>>    resHists           (N_REGIONS, std::vector<TH1D*>(nCuts, nullptr));
  std::vector<std::vector<Long64_t>> nPass              (N_REGIONS, std::vector<Long64_t>(nCuts, 0));
  std::vector<std::vector<Long64_t>> nTot               (N_REGIONS, std::vector<Long64_t>(nCuts, 0));
  std::vector<std::vector<Long64_t>> nNullSel           (N_REGIONS, std::vector<Long64_t>(nCuts, 0));
  std::vector<std::vector<double>>   sumTrackCount      (N_REGIONS, std::vector<double>(nCuts, 0.0));
  std::vector<std::vector<double>>   sumPurity          (N_REGIONS, std::vector<double>(nCuts, 0.0));
  std::vector<std::vector<Long64_t>> nEventsWithClusters(N_REGIONS, std::vector<Long64_t>(nCuts, 0));

  for (int r = 0; r < N_REGIONS; ++r)
    for (int ic = 0; ic < nCuts; ++ic)
      resHists[r][ic] = new TH1D(
        Form("res_r%d_cut%.2f", r, cuts[ic]),
        Form("REFINED #Deltat (%s, refCut=%.2f#sigma);#Delta t [ps];Entries",
             regionLabel[r].c_str(), cuts[ic]),
        nResBins, RES_MIN, RES_MAX);

  std::vector<std::vector<TH1D*>> pullHist(N_REGIONS, std::vector<TH1D*>(nCuts, nullptr));
  for (int r = 0; r < N_REGIONS; ++r)
    for (int ic = 0; ic < nCuts; ++ic)
      pullHist[r][ic] = new TH1D(
        Form("pull_r%d_cut%.2f", r, cuts[ic]),
        Form("Track Time Pull (%s, refCut=%.2f#sigma);(t_{trk}-t_{clust})/#sigma_{t,trk};Normalised Entries",
             regionLabel[r].c_str(), cuts[ic]),
        80, -8.0, 8.0);

  // --- Event loop ---
  std::cout << "Starting sweep event loop\n";
  Long64_t evProcessed = 0;

  while (reader.Next()) {
    const Long64_t readNum = chain.GetReadEntry() + 1;

    if (readNum % 200 == 0)
      std::cout << "Progress: " << readNum << "/" << nProcess << "\r" << std::flush;

    if (readNum > nProcess) break;

    // ── Event selection ──────────────────────────────────────────────────────
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

    // ── Cone clustering — once per event (cut-independent) ───────────────────
    std::vector<Cluster> coneClusters =
      clusterTracksInTime(
        tracks, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false,
        /*checkTimeValid=*/true,
        /*smearRes=*/10.0,
        ClusteringMethod::CONE,
        /*usez0=*/false);

    // Apply MIN_CLUSTER_TRACKS filter and sort by TRKPTZ descending
    std::vector<Cluster> qualClusters;
    for (const auto& c : coneClusters)
      if (c.nConstituents >= MIN_CLUSTER_TRACKS) qualClusters.push_back(c);

    if (qualClusters.empty()) {
      for (int ic = 0; ic < nCuts; ++ic) {
        nTot   [region][ic]++;
        nNullSel[region][ic]++;
      }
      continue;
    }

    std::sort(qualClusters.begin(), qualClusters.end(), [](const Cluster& a, const Cluster& b) {
      return a.scores.at(Score::TRKPTZ.id) > b.scores.at(Score::TRKPTZ.id);
    });

    // Pool trackIndices from top CONE_ITER_K qualifying clusters (cut-independent)
    std::vector<int> pooledTracks;
    int k = std::min((int)qualClusters.size(), CONE_ITER_K);
    for (int i = 0; i < k; ++i) {
      const auto& idx = qualClusters[i].trackIndices;
      pooledTracks.insert(pooledTracks.end(), idx.begin(), idx.end());
    }

    // ── Per-cut iterative refinement ─────────────────────────────────────────
    for (int ic = 0; ic < nCuts; ++ic) {
      nTot[region][ic]++;

      std::vector<Cluster> refined =
        clusterTracksInTime(
          pooledTracks, &branch, cuts[ic],
          /*useSmearedTimes=*/false,
          /*checkTimeValid=*/true,
          /*smearRes=*/10.0,
          ClusteringMethod::ITERATIVE,
          /*usez0=*/false,
          /*sortTracks=*/false,
          /*calcPurityFlag=*/true);

      // Apply MIN_CLUSTER_TRACKS filter
      std::vector<Cluster> qualRefined;
      for (const auto& c : refined)
        if (c.nConstituents >= MIN_CLUSTER_TRACKS) qualRefined.push_back(c);

      if (qualRefined.empty()) {
        nNullSel[region][ic]++;
        continue;
      }

      Cluster best = chooseCluster(qualRefined, Score::TRKPTZ);

      double diff = best.values.at(0) - branch.truthVtxTime[0];
      resHists[region][ic]->Fill(diff);

      nEventsWithClusters[region][ic]++;
      sumTrackCount      [region][ic] += (double)best.trackIndices.size();
      sumPurity          [region][ic] += best.purity;

      // Track time pull within selected cluster
      for (int idx : best.trackIndices) {
        if (branch.trackTimeValid[idx] != 1) continue;
        double pull = (branch.trackTime[idx] - best.values.at(0)) / branch.trackTimeRes[idx];
        pullHist[region][ic]->Fill(pull);
      }

      if (best.passEfficiency(&branch))
        nPass[region][ic]++;
    }

    evProcessed++;
  }

  std::cout << "\nFinished event loop (" << evProcessed << " events processed)\n";

  // --- Build per-region TGraphErrors ---
  TGraphErrors* grEff       [N_REGIONS];
  TGraphErrors* grSigma     [N_REGIONS];
  TGraphErrors* grNullSel   [N_REGIONS];
  TGraphErrors* grTrackCount[N_REGIONS];
  TGraphErrors* grPurity    [N_REGIONS];

  for (int r = 0; r < N_REGIONS; ++r) {
    grEff       [r] = new TGraphErrors(nCuts);
    grSigma     [r] = new TGraphErrors(nCuts);
    grNullSel   [r] = new TGraphErrors(nCuts);
    grTrackCount[r] = new TGraphErrors(nCuts);
    grPurity    [r] = new TGraphErrors(nCuts);

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

      // Null selection rate
      double nullFrac = (nTot[r][ic] > 0)
        ? (double)nNullSel[r][ic] / (double)nTot[r][ic] : 0.0;
      grNullSel[r]->SetPoint(ic, dc, nullFrac);
      grNullSel[r]->SetPointError(ic, 0.0, 0.0);

      // Mean track count in selected cluster
      double meanTk = (nEventsWithClusters[r][ic] > 0)
        ? sumTrackCount[r][ic] / (double)nEventsWithClusters[r][ic] : 0.0;
      grTrackCount[r]->SetPoint(ic, dc, meanTk);
      grTrackCount[r]->SetPointError(ic, 0.0, 0.0);

      // Mean purity of selected cluster
      double meanPurity = (nEventsWithClusters[r][ic] > 0)
        ? sumPurity[r][ic] / (double)nEventsWithClusters[r][ic] : 0.0;
      grPurity[r]->SetPoint(ic, dc, meanPurity);
      grPurity[r]->SetPointError(ic, 0.0, 0.0);
    }
  }

  // --- Plotting ---
  TCanvas* canvas = new TCanvas("canvas", "Sweep Results", 800, 600);
  canvas->SetLeftMargin(0.15);

  const std::string outFile = std::string(OUT_DIR) + "/dist_refine_sweep.pdf";
  canvas->Print((outFile + "[").c_str());

  const double xLo = REFINE_MIN - REFINE_STEP * 0.5;
  const double xHi = REFINE_MAX + REFINE_STEP * 0.5;
  const char*  xtitle = "Refinement Distance Cut [#sigma]";

  // Plot 1: Efficiency vs refinement cut
  saveGraphPair(outFile, canvas, grEff[0], grEff[1],
    Form("REFINED Efficiency vs. Refinement Cut (k=%d)", CONE_ITER_K),
    xtitle, "Overall Efficiency",
    xLo, xHi, 0.0, 1.3, HS_TRACK_SPLIT,
    0.99, "99% Efficiency");

  // Plot 2: Core sigma vs refinement cut
  saveGraphPair(outFile, canvas, grSigma[0], grSigma[1],
    Form("REFINED Core #sigma vs. Refinement Cut (k=%d)", CONE_ITER_K),
    xtitle, "Core #sigma [ps]",
    xLo, xHi, 0.0, 60.0, HS_TRACK_SPLIT,
    15.0, "15 ps");

  // Plot 3: Null selection rate vs refinement cut
  saveGraphPair(outFile, canvas, grNullSel[0], grNullSel[1],
    Form("Null Selection Rate vs. Refinement Cut (k=%d)", CONE_ITER_K),
    xtitle, "Fraction of events with no qualifying refined cluster",
    xLo, xHi, 0.0, 1.0, HS_TRACK_SPLIT);

  // Plot 4: Mean selected cluster track count vs refinement cut
  double tcYMax = 1.0;
  for (int r = 0; r < N_REGIONS; ++r) {
    for (int ic = 0; ic < nCuts; ++ic) {
      double xv, yv;
      grTrackCount[r]->GetPoint(ic, xv, yv);
      if (yv > tcYMax) tcYMax = yv;
    }
  }
  saveGraphPair(outFile, canvas, grTrackCount[0], grTrackCount[1],
    Form("Mean Refined Cluster Track Count vs. Refinement Cut (k=%d)", CONE_ITER_K),
    xtitle, "Mean N_{tracks} in selected refined cluster",
    xLo, xHi, 0.0, tcYMax * 1.2, HS_TRACK_SPLIT);

  // Plot 5: Mean purity of selected cluster vs refinement cut
  saveGraphPair(outFile, canvas, grPurity[0], grPurity[1],
    Form("Mean Selected-Cluster Purity vs. Refinement Cut (k=%d)", CONE_ITER_K),
    xtitle, "Mean Purity (HS p_{T} fraction)",
    xLo, xHi, 0.5, 1.1, HS_TRACK_SPLIT,
    0.75, "75% purity");

  // Plot 6: Normalised residual shape overlay per region (two pages)
  for (int r = 0; r < N_REGIONS; ++r) {
    std::vector<TH1D*> resClones;
    std::vector<std::string> resLabels6;
    for (int ic = 0; ic < nCuts; ++ic) {
      if (ic % 2 != 0 && ic != nCuts - 1) continue;
      TH1D* h = (TH1D*)resHists[r][ic]->Clone(Form("res_clone_r%d_%d", r, ic));
      if (h->Integral() == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[resClones.size() % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-400, 400);
      h->SetTitle(Form("REFINED #Deltat shape (%s);#Delta t [ps];Normalised Entries",
                       regionLabel[r].c_str()));
      resLabels6.push_back(Form("%.2f#sigma", cuts[ic]));
      resClones.push_back(h);
    }
    if (resClones.empty()) continue;
    double yMax6 = 0.0, yMin6 = 1e30;
    for (auto* h : resClones)
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double v = h->GetBinContent(b);
        if (v > 0.0) { yMax6 = std::max(yMax6, v); yMin6 = std::min(yMin6, v); }
      }
    yMax6 *= 3.0;  yMin6 = std::max(yMin6 * 0.3, 1e-6);
    resClones[0]->GetYaxis()->SetRangeUser(yMin6, yMax6);
    canvas->SetLogy(true);
    LegendBox lb6 = bestLegendCornerHist(
      resClones, -400.0, 400.0, yMin6, yMax6, (int)resClones.size() + 1, true);
    TLegend* leg6 = new TLegend(lb6.x1, lb6.y1, lb6.x2, lb6.y2);
    leg6->SetTextSize(0.025);
    leg6->SetHeader(regionLabel[r].c_str(), "C");
    bool first6 = true;
    for (size_t i = 0; i < resClones.size(); ++i) {
      if (first6) { resClones[i]->Draw("HIST"); first6 = false; }
      else          resClones[i]->Draw("HIST SAME");
      leg6->AddEntry(resClones[i], resLabels6[i].c_str(), "l");
    }
    leg6->Draw("SAME");
    canvas->Print(outFile.c_str());
    canvas->SetLogy(false);
    for (auto* h : resClones) delete h;
  }

  // Plot 7: Track time pull within selected cluster, overlay per region (two pages)
  for (int r = 0; r < N_REGIONS; ++r) {
    std::vector<TH1D*> pullClones;
    std::vector<std::string> pullLabels;
    for (int ic = 0; ic < nCuts; ++ic) {
      if (ic % 2 != 0 && ic != nCuts - 1) continue;
      if (pullHist[r][ic]->Integral() == 0) continue;
      TH1D* h = (TH1D*)pullHist[r][ic]->Clone(Form("pull_clone_r%d_%d", r, ic));
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[pullClones.size() % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-8.0, 8.0);
      h->SetTitle(Form("Track Time Pull (%s);(t_{trk}-t_{clust})/#sigma_{t,trk};Normalised Entries",
                       regionLabel[r].c_str()));
      pullLabels.push_back(Form("%.2f#sigma", cuts[ic]));
      pullClones.push_back(h);
    }
    if (pullClones.empty()) continue;
    double yMaxP = 0.0, yMinP = 1e30;
    for (auto* h : pullClones)
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double v = h->GetBinContent(b);
        if (v > 0.0) { yMaxP = std::max(yMaxP, v); yMinP = std::min(yMinP, v); }
      }
    yMaxP *= 3.0;  yMinP = std::max(yMinP * 0.3, 1e-6);
    pullClones[0]->GetYaxis()->SetRangeUser(yMinP, yMaxP);
    canvas->SetLogy(true);
    LegendBox lbP = bestLegendCornerHist(
      pullClones, -8.0, 8.0, yMinP, yMaxP, (int)pullClones.size() + 1, true);
    TLegend* legP = new TLegend(lbP.x1, lbP.y1, lbP.x2, lbP.y2);
    legP->SetTextSize(0.025);
    legP->SetHeader(regionLabel[r].c_str(), "C");
    bool firstP = true;
    for (size_t i = 0; i < pullClones.size(); ++i) {
      if (firstP) { pullClones[i]->Draw("HIST"); firstP = false; }
      else          pullClones[i]->Draw("HIST SAME");
      legP->AddEntry(pullClones[i], pullLabels[i].c_str(), "l");
    }
    legP->Draw("SAME");
    canvas->Print(outFile.c_str());
    canvas->SetLogy(false);
    for (auto* h : pullClones) delete h;
  }

  // --- Pseudo-ROC plots (pages 7-9) ---
  // Build label strings: "0.50σ", "0.75σ", ...
  std::vector<std::string> cutLabels(nCuts);
  for (int ic = 0; ic < nCuts; ++ic)
    cutLabels[ic] = Form("%.2f#sigma", cuts[ic]);

  const auto effLow   = graphYVals(grEff[0]);
  const auto effHigh  = graphYVals(grEff[1]);
  const auto sigLow   = graphYVals(grSigma[0]);
  const auto sigHigh  = graphYVals(grSigma[1]);
  const auto purLow   = graphYVals(grPurity[0]);
  const auto purHigh  = graphYVals(grPurity[1]);

  // ROC page 7: Purity vs Efficiency
  saveRocPlot(outFile, canvas,
    effLow, purLow, effHigh, purHigh, cutLabels,
    Form("Purity vs. Efficiency (k=%d)", CONE_ITER_K),
    "Overall Efficiency", "Mean Cluster Purity",
    HS_TRACK_SPLIT);

  // ROC page 8: Purity vs Core sigma
  saveRocPlot(outFile, canvas,
    sigLow, purLow, sigHigh, purHigh, cutLabels,
    Form("Purity vs. Core #sigma (k=%d)", CONE_ITER_K),
    "Core #sigma [ps]", "Mean Cluster Purity",
    HS_TRACK_SPLIT);

  // ROC page 9: Efficiency vs Core sigma
  saveRocPlot(outFile, canvas,
    sigLow, effLow, sigHigh, effHigh, cutLabels,
    Form("Efficiency vs. Core #sigma (k=%d)", CONE_ITER_K),
    "Core #sigma [ps]", "Overall Efficiency",
    HS_TRACK_SPLIT);

  canvas->Print((outFile + "]").c_str());
  std::cout << "Saved: " << outFile << '\n';

  // --- Summary tables ---
  for (int r = 0; r < N_REGIONS; ++r) {
    std::cout << "\n--- Summary Table [" << regionLabel[r] << "] ---\n";
    std::cout << Form("%-10s  %-10s  %-10s  %-12s  %-12s  %-12s\n",
                      "refCut", "Efficiency", "CoreSigma", "NullSelRate", "MeanNTracks", "MeanPurity");
    for (int ic = 0; ic < nCuts; ++ic) {
      double dc, eff, sig, nul, ntk, pur, dummy;
      grEff       [r]->GetPoint(ic, dc,    eff);
      grSigma     [r]->GetPoint(ic, dummy, sig);
      grNullSel   [r]->GetPoint(ic, dummy, nul);
      grTrackCount[r]->GetPoint(ic, dummy, ntk);
      grPurity    [r]->GetPoint(ic, dummy, pur);
      (void)dummy;
      std::cout << Form("%-10.2f  %-10.4f  %-10.2f  %-12.4f  %-12.2f  %-12.4f\n",
                        dc, eff, sig, nul, ntk, pur);
    }
  }

  return 0;
}
