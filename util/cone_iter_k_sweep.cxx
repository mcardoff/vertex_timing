// cone_iter_k_sweep.cxx
//
// 1D parameter sweep over CONE_ITER_K: the number of top cone clusters
// whose tracks are pooled before running a tighter iterative pass (the
// REFINED algorithm).  For each k in [K_MIN, K_MAX] we:
//   - Run cone clustering once per event (k-independent, expensive step)
//   - Pool the trackIndices of the top-k qualifying cone clusters by TRKPTZ
//   - Re-cluster the pooled tracks with ClusteringMethod::ITERATIVE at DIST_CUT_REFINE
//   - Apply MIN_CLUSTER_TRACKS filter; evaluate efficiency and collect residuals
//
// A key diagnostic is the HS cone-cluster rank CDF:
//   "For what fraction of events is the HS cluster in the top-k by TRKPTZ?"
// This is the theoretical maximum efficiency REFINED can achieve at a given k.
// Comparing it to actual efficiency shows how much loss comes from pool size
// (k too small) vs. from the iterative step itself (fragmentation, PU winning).
//
// Events are split into two regions by the number of forward hard-scatter
// tracks (nFTrackHS):
//   Region 0 — Low HS:  nFTrackHS <= HS_TRACK_SPLIT
//   Region 1 — High HS: nFTrackHS >  HS_TRACK_SPLIT
//
// Output PDF: ../parameter_sweep/cone_iter_k_sweep.pdf
//
// ── Configuration ─────────────────────────────────────────────────────────────
//   MAX_EVENTS    – number of events to process (-1 = full sample)
//   K_MIN         – smallest k to test (≥ 1)
//   K_MAX         – largest  k to test (inclusive)
//   HS_TRACK_SPLIT – boundary between the two HS-track regions (inclusive low)
// ─────────────────────────────────────────────────────────────────────────────

#include "sweep_utilities.h"

using namespace MyUtl;

// ── Sweep configuration ───────────────────────────────────────────────────────
static constexpr Long64_t MAX_EVENTS     = -1;  // -1 = full sample
static constexpr int      K_MIN          = 1;
static constexpr int      K_MAX          = 6;
static constexpr int      HS_TRACK_SPLIT = 5;  // ≤ this → "low HS" region
// ─────────────────────────────────────────────────────────────────────────────

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

  std::cout << "Events to process : " << nProcess << " / " << nTotal << '\n';
  std::cout << "k sweep           : " << K_MIN << " .. " << K_MAX << '\n';
  std::cout << "HS track split    : nFTrackHS <= " << HS_TRACK_SPLIT
            << "  (low) / > " << HS_TRACK_SPLIT << " (high)\n";

  const int nK = K_MAX - K_MIN + 1;
  const int nResBins = (int)((RES_MAX - RES_MIN) / RES_BWID);

  const std::string regionLabel[N_REGIONS] = {
    Form("nFTrackHS <= %d", HS_TRACK_SPLIT),
    Form("nFTrackHS > %d",  HS_TRACK_SPLIT)
  };

  // --- Per-k, per-region accumulators ---
  // [region][k_index]
  std::vector<std::vector<TH1D*>>    resHists           (N_REGIONS, std::vector<TH1D*>(nK, nullptr));
  std::vector<std::vector<Long64_t>> nPass              (N_REGIONS, std::vector<Long64_t>(nK, 0));
  std::vector<std::vector<Long64_t>> nTot               (N_REGIONS, std::vector<Long64_t>(nK, 0));
  std::vector<std::vector<Long64_t>> nNullSel           (N_REGIONS, std::vector<Long64_t>(nK, 0));
  std::vector<std::vector<double>>   sumTrackCount      (N_REGIONS, std::vector<double>(nK, 0.0));
  std::vector<std::vector<double>>   sumPurity          (N_REGIONS, std::vector<double>(nK, 0.0));
  std::vector<std::vector<Long64_t>> nEventsWithClusters(N_REGIONS, std::vector<Long64_t>(nK, 0));

  // HS cone-cluster rank histogram — filled once per event, not per k.
  // Bin i corresponds to rank i (1-based); overflow catches rank > K_MAX.
  TH1D* rankHist[N_REGIONS];
  for (int r = 0; r < N_REGIONS; ++r) {
    rankHist[r] = new TH1D(
      Form("hs_rank_r%d", r),
      Form("HS Cone Cluster Rank by TRKPTZ (%s);Rank;Fraction of Events",
           regionLabel[r].c_str()),
      K_MAX + 1, 0.5, K_MAX + 1.5);  // bins 1..K_MAX, overflow = rank > K_MAX
  }

  for (int r = 0; r < N_REGIONS; ++r)
    for (int ik = 0; ik < nK; ++ik)
      resHists[r][ik] = new TH1D(
        Form("res_r%d_k%d", r, K_MIN + ik),
        Form("REFINED #Deltat (%s, k=%d);#Delta t [ps];Entries",
             regionLabel[r].c_str(), K_MIN + ik),
        nResBins, RES_MIN, RES_MAX);

  // Refined sub-cluster count distribution — filled once per (event, k)
  std::vector<std::vector<TH1D*>> nSubclustHist(N_REGIONS, std::vector<TH1D*>(nK, nullptr));
  for (int r = 0; r < N_REGIONS; ++r)
    for (int ik = 0; ik < nK; ++ik)
      nSubclustHist[r][ik] = new TH1D(
        Form("nsubclust_r%d_k%d", r, K_MIN + ik),
        Form("Refined Sub-cluster Count (%s, k=%d);N qualified refined clusters;Events",
             regionLabel[r].c_str(), K_MIN + ik),
        K_MAX + 2, -0.5, K_MAX + 1.5);

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

    // ── Cone clustering — once per event (k-independent) ────────────────────
    std::vector<Cluster> coneClusters =
      clusterTracksInTime(
        tracks, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false,
        /*checkTimeValid=*/true,
        /*smearRes=*/10.0,
        ClusteringMethod::CONE,
        /*usez0=*/false);

    // Apply MIN_CLUSTER_TRACKS filter
    std::vector<Cluster> qualClusters;
    for (const auto& c : coneClusters)
      if (c.nConstituents >= MIN_CLUSTER_TRACKS) qualClusters.push_back(c);

    if (qualClusters.empty()) {
      // No qualifying cone clusters: fill denominator for all k then skip
      for (int ik = 0; ik < nK; ++ik) {
        nTot   [region][ik]++;
        nNullSel[region][ik]++;
      }
      continue;
    }

    // Sort qualifying clusters by TRKPTZ descending — this order defines ranks
    std::sort(qualClusters.begin(), qualClusters.end(), [](const Cluster& a, const Cluster& b) {
      return a.scores.at(Score::TRKPTZ) > b.scores.at(Score::TRKPTZ);
    });

    // ── HS cone-cluster rank (filled once per event) ─────────────────────────
    // Find the qualifying cluster with the most HS-tagged tracks.
    // HS tracks: branch.trackToTruthvtx[idx] == 0  (truth vertex index 0 = HS)
    {
      int hsRank    = -1;   // 1-based rank of the HS-richest cluster; -1 = not found
      int maxHSTracks = 0;
      for (int ic = 0; ic < (int)qualClusters.size(); ++ic) {
        int nHS = 0;
        for (int idx : qualClusters[ic].trackIndices)
          if (branch.trackToTruthvtx[idx] == 0) ++nHS;
        if (nHS > maxHSTracks) {
          maxHSTracks = nHS;
          hsRank = ic + 1;  // 1-based
        }
      }
      if (hsRank > 0)
        rankHist[region]->Fill(hsRank);
    }

    // ── Per-k inner loop ─────────────────────────────────────────────────────
    for (int ik = 0; ik < nK; ++ik) {
      const int k = K_MIN + ik;

      nTot[region][ik]++;

      // Pool trackIndices from top-k qualifying cone clusters
      std::vector<int> pooledTracks;
      int limit = std::min(k, (int)qualClusters.size());
      for (int i = 0; i < limit; ++i) {
        const auto& idx = qualClusters[i].trackIndices;
        pooledTracks.insert(pooledTracks.end(), idx.begin(), idx.end());
      }

      // Re-cluster with tighter iterative pass
      std::vector<Cluster> refined =
        clusterTracksInTime(
          pooledTracks, &branch, DIST_CUT_REFINE,
          /*useSmearedTimes=*/false,
          /*checkTimeValid=*/true,
          /*smearRes=*/10.0,
          ClusteringMethod::ITERATIVE,
          /*usez0=*/false,
          /*sortTracks=*/false,
          /*calcPurityFlag=*/true);

      // Apply MIN_CLUSTER_TRACKS filter to refined sub-clusters
      std::vector<Cluster> qualRefined;
      for (const auto& c : refined)
        if (c.nConstituents >= MIN_CLUSTER_TRACKS) qualRefined.push_back(c);

      nSubclustHist[region][ik]->Fill((int)qualRefined.size());

      if (qualRefined.empty()) {
        nNullSel[region][ik]++;
        continue;
      }

      Cluster best = chooseCluster(qualRefined, Score::TRKPTZ);

      double diff = best.values.at(0) - branch.truthVtxTime[0];
      resHists[region][ik]->Fill(diff);

      nEventsWithClusters[region][ik]++;
      sumTrackCount      [region][ik] += (double)best.trackIndices.size();
      sumPurity          [region][ik] += best.purity;

      if (best.passEfficiency(&branch))
        nPass[region][ik]++;
    }

    evProcessed++;
  }

  std::cout << "\nFinished event loop (" << evProcessed << " events processed)\n";

  // --- Build per-region TGraphErrors ---
  TGraphErrors* grEff       [N_REGIONS];
  TGraphErrors* grSigma     [N_REGIONS];
  TGraphErrors* grCoverage  [N_REGIONS];
  TGraphErrors* grNullSel   [N_REGIONS];
  TGraphErrors* grTrackCount[N_REGIONS];
  TGraphErrors* grPurity    [N_REGIONS];

  for (int r = 0; r < N_REGIONS; ++r) {
    grEff       [r] = new TGraphErrors(nK);
    grSigma     [r] = new TGraphErrors(nK);
    grCoverage  [r] = new TGraphErrors(nK);
    grNullSel   [r] = new TGraphErrors(nK);
    grTrackCount[r] = new TGraphErrors(nK);
    grPurity    [r] = new TGraphErrors(nK);

    // HS coverage CDF: cumulative fraction of events with HS rank ≤ k
    double totalRankEvents = rankHist[r]->Integral(0, K_MAX + 2);  // include overflow

    for (int ik = 0; ik < nK; ++ik) {
      const double xk = (double)(K_MIN + ik);

      // Efficiency
      double eff    = (nTot[r][ik] > 0) ? (double)nPass[r][ik] / (double)nTot[r][ik] : 0.0;
      double effErr = (nTot[r][ik] > 0)
        ? std::sqrt(eff * (1.0 - eff) / (double)nTot[r][ik]) : 0.0;
      grEff[r]->SetPoint(ik, xk, eff);
      grEff[r]->SetPointError(ik, 0.0, effErr);

      // Core sigma
      auto [sigma, sigmaErr] = fitCoreSigma(resHists[r][ik], r * nK + ik);
      grSigma[r]->SetPoint(ik, xk, sigma);
      grSigma[r]->SetPointError(ik, 0.0, sigmaErr);

      // HS coverage CDF: integral of rankHist up to and including bin k
      double cumulative = (totalRankEvents > 0)
        ? rankHist[r]->Integral(1, K_MIN + ik) / totalRankEvents : 0.0;
      grCoverage[r]->SetPoint(ik, xk, cumulative);
      grCoverage[r]->SetPointError(ik, 0.0, 0.0);

      // Null selection rate
      double nullFrac = (nTot[r][ik] > 0)
        ? (double)nNullSel[r][ik] / (double)nTot[r][ik] : 0.0;
      grNullSel[r]->SetPoint(ik, xk, nullFrac);
      grNullSel[r]->SetPointError(ik, 0.0, 0.0);

      // Mean refined track count
      double meanTk = (nEventsWithClusters[r][ik] > 0)
        ? sumTrackCount[r][ik] / (double)nEventsWithClusters[r][ik] : 0.0;
      grTrackCount[r]->SetPoint(ik, xk, meanTk);
      grTrackCount[r]->SetPointError(ik, 0.0, 0.0);

      // Mean purity of selected refined cluster
      double meanPurity = (nEventsWithClusters[r][ik] > 0)
        ? sumPurity[r][ik] / (double)nEventsWithClusters[r][ik] : 0.0;
      grPurity[r]->SetPoint(ik, xk, meanPurity);
      grPurity[r]->SetPointError(ik, 0.0, 0.0);
    }
  }

  // --- Plotting ---
  TCanvas* canvas = new TCanvas("canvas", "Sweep Results", 800, 600);
  canvas->SetLeftMargin(0.15);

  const std::string outFile = std::string(OUT_DIR) + "/cone_iter_k_sweep.pdf";
  canvas->Print((outFile + "[").c_str());

  const double xLo = K_MIN - 0.5;
  const double xHi = K_MAX + 0.5;
  const char*  xtitle = "k (top cone clusters refined)";

  // Plot 1: Efficiency vs k
  saveGraphPair(outFile, canvas, grEff[0], grEff[1],
    "REFINED Efficiency vs. k",
    xtitle, "Overall Efficiency",
    xLo, xHi, 0.0, 1.3, HS_TRACK_SPLIT,
    0.99, "99% Efficiency");

  // Plot 2: Core sigma vs k
  saveGraphPair(outFile, canvas, grSigma[0], grSigma[1],
    "REFINED Core #sigma vs. k",
    xtitle, "Core #sigma [ps]",
    xLo, xHi, 0.0, 60.0, HS_TRACK_SPLIT,
    15.0, "15 ps");

  // Plot 3: HS cone-cluster rank CDF (theoretical max efficiency)
  saveGraphPair(outFile, canvas, grCoverage[0], grCoverage[1],
    "HS Coverage (Cone Rank CDF) vs. k",
    xtitle, "Fraction of events with HS rank #leq k",
    xLo, xHi, 0.0, 1.3, HS_TRACK_SPLIT);

  // Plot 4: Null selection rate vs k
  saveGraphPair(outFile, canvas, grNullSel[0], grNullSel[1],
    "Null Selection Rate vs. k",
    xtitle, "Fraction of events with no qualifying refined cluster",
    xLo, xHi, 0.0, 1.0, HS_TRACK_SPLIT);

  // Plot 5: Mean refined cluster track count vs k
  double tcYMax = 1.0;
  for (int r = 0; r < N_REGIONS; ++r) {
    for (int ik = 0; ik < nK; ++ik) {
      double xv, yv;
      grTrackCount[r]->GetPoint(ik, xv, yv);
      if (yv > tcYMax) tcYMax = yv;
    }
  }
  saveGraphPair(outFile, canvas, grTrackCount[0], grTrackCount[1],
    "Mean Refined Cluster Track Count vs. k",
    xtitle, "Mean N_{tracks} in selected refined cluster",
    xLo, xHi, 0.0, tcYMax * 1.2, HS_TRACK_SPLIT);

  // Plot 6: Normalised residual shape overlay per region (two pages)
  for (int r = 0; r < N_REGIONS; ++r) {
    std::vector<TH1D*> resClones;
    std::vector<std::string> resLabels6;
    for (int ik = 0; ik < nK; ++ik) {
      TH1D* h = (TH1D*)resHists[r][ik]->Clone(Form("res_clone_r%d_%d", r, ik));
      if (h->Integral() == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[resClones.size() % COLORS.size()]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetRangeUser(-400, 400);
      h->SetTitle(Form("REFINED #Deltat shape (%s);#Delta t [ps];Normalised Entries",
                       regionLabel[r].c_str()));
      resLabels6.push_back(Form("k=%d", K_MIN + ik));
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

  // Plot 7: Refined sub-cluster count distribution, overlay per k per region (two pages)
  for (int r = 0; r < N_REGIONS; ++r) {
    std::vector<TH1D*> sClones;
    std::vector<std::string> sLabels;
    for (int ik = 0; ik < nK; ++ik) {
      if (nSubclustHist[r][ik]->Integral() == 0) continue;
      TH1D* h = (TH1D*)nSubclustHist[r][ik]->Clone(Form("nsubclust_clone_r%d_%d", r, ik));
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(COLORS[sClones.size() % COLORS.size()]);
      h->SetLineWidth(2);
      h->SetTitle(Form("Refined Sub-cluster Count (%s);N qualified refined clusters;Normalised Events",
                       regionLabel[r].c_str()));
      sLabels.push_back(Form("k=%d", K_MIN + ik));
      sClones.push_back(h);
    }
    if (sClones.empty()) continue;
    double yMaxS = 0.0;
    for (auto* h : sClones)
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        yMaxS = std::max(yMaxS, h->GetBinContent(b));
    yMaxS *= 1.3;
    if (yMaxS <= 0.0) { for (auto* h : sClones) delete h; continue; }
    sClones[0]->GetYaxis()->SetRangeUser(0.0, yMaxS);
    LegendBox lbS = bestLegendCornerHist(
      sClones, -0.5, K_MAX + 1.5, 0.0, yMaxS, (int)sClones.size() + 1, false);
    TLegend* legS = new TLegend(lbS.x1, lbS.y1, lbS.x2, lbS.y2);
    legS->SetTextSize(0.025);
    legS->SetHeader(regionLabel[r].c_str(), "C");
    bool firstS = true;
    for (size_t i = 0; i < sClones.size(); ++i) {
      if (firstS) { sClones[i]->Draw("HIST"); firstS = false; }
      else          sClones[i]->Draw("HIST SAME");
      legS->AddEntry(sClones[i], sLabels[i].c_str(), "l");
    }
    legS->Draw("SAME");
    canvas->Print(outFile.c_str());
    for (auto* h : sClones) delete h;
  }

  // --- Pseudo-ROC plots (pages 7-9) ---
  // Build label strings: "k=1", "k=2", ...
  std::vector<std::string> kLabels(nK);
  for (int ik = 0; ik < nK; ++ik)
    kLabels[ik] = Form("k=%d", K_MIN + ik);

  const auto effLow   = graphYVals(grEff[0]);
  const auto effHigh  = graphYVals(grEff[1]);
  const auto sigLow   = graphYVals(grSigma[0]);
  const auto sigHigh  = graphYVals(grSigma[1]);
  const auto purLow   = graphYVals(grPurity[0]);
  const auto purHigh  = graphYVals(grPurity[1]);

  // ROC page 7: Purity vs Efficiency
  saveRocPlot(outFile, canvas,
    effLow, purLow, effHigh, purHigh, kLabels,
    "Purity vs. Efficiency (k sweep)",
    "Overall Efficiency", "Mean Cluster Purity",
    HS_TRACK_SPLIT);

  // ROC page 8: Purity vs Core sigma
  saveRocPlot(outFile, canvas,
    sigLow, purLow, sigHigh, purHigh, kLabels,
    "Purity vs. Core #sigma (k sweep)",
    "Core #sigma [ps]", "Mean Cluster Purity",
    HS_TRACK_SPLIT);

  // ROC page 9: Efficiency vs Core sigma
  saveRocPlot(outFile, canvas,
    sigLow, effLow, sigHigh, effHigh, kLabels,
    "Efficiency vs. Core #sigma (k sweep)",
    "Core #sigma [ps]", "Overall Efficiency",
    HS_TRACK_SPLIT);

  canvas->Print((outFile + "]").c_str());
  std::cout << "Saved: " << outFile << '\n';

  // --- Summary tables ---
  for (int r = 0; r < N_REGIONS; ++r) {
    std::cout << "\n--- Summary Table [" << regionLabel[r] << "] ---\n";
    std::cout << Form("%-6s  %-10s  %-10s  %-12s  %-12s  %-12s\n",
                      "k", "Efficiency", "CoreSigma", "HSCoverage", "NullSelRate", "MeanNTracks");
    for (int ik = 0; ik < nK; ++ik) {
      const int k = K_MIN + ik;
      double xv, eff, sig, cov, nul, ntk, dummy;
      grEff       [r]->GetPoint(ik, xv,    eff);
      grSigma     [r]->GetPoint(ik, dummy, sig);
      grCoverage  [r]->GetPoint(ik, dummy, cov);
      grNullSel   [r]->GetPoint(ik, dummy, nul);
      grTrackCount[r]->GetPoint(ik, dummy, ntk);
      (void)dummy;
      std::cout << Form("%-6d  %-10.4f  %-10.2f  %-12.4f  %-12.4f  %-12.2f\n",
                        k, eff, sig, cov, nul, ntk);
    }
  }

  return 0;
}
