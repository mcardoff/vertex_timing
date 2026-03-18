// ---------------------------------------------------------------------------
// rate_diagnostics.cxx
//
//   Computes two selection-agnostic diagnostic rates:
//
//   1. Misclustering rate — for HGTD and Ideal Res. scenarios separately:
//        Among all qualifying clusters in an event, does ANY cluster achieve
//        purity > 75%?  Events where none do are "misclustered" — no selection
//        algorithm can recover them.  The HGTD rate includes timing-induced
//        contamination; the Ideal Res. rate is the spatial/kinematic floor.
//
//   2. Misassignment rate — PASS oracle comparison:
//        How often does an event fail the PASS oracle with real HGTD times
//        but pass with Ideal Res. smeared times?  This is fully selection-
//        agnostic: it asks whether fixing timing converts a timing failure
//        into a timing success, regardless of which cluster is chosen.
//
//   Output: three tables binned by n Forward HS Tracks (with valid HGTD time),
//   matching the style of clustering_dt.cxx.
// ---------------------------------------------------------------------------

#include <RtypesCore.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include "clustering_constants.h"
#include "event_processing.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Per-bin accumulator
// ---------------------------------------------------------------------------

struct BinData {
  double total           = 0; // all events passing selection in this bin
  double hgtd_pure       = 0; // events where ≥1 HGTD cluster has purity > 0.75
  double ideal_pure      = 0; // events where ≥1 Ideal Res. cluster has purity > 0.75
  double misassign       = 0; // events where HGTD |dt| − ideal |dt| > 60 ps (best-purity cluster)
  double total_pure      = 0; // subset: events with HGTD pure cluster (no misclustering)
  double misassign_pure  = 0; // misassignment events within the no-misclustering subset
};

static const int NBINS = static_cast<int>(FOLD_HS_TRACK) + 1;

// ---------------------------------------------------------------------------
// Misclustering table printer
//   Shows n_hasPure | n_noPure | Total | Misclustering rate (%)
//   Misclustering rate = n_noPure / Total = fraction without any pure cluster
// ---------------------------------------------------------------------------

void printMisclusteringTable(
  const std::string&          header,
  const std::vector<BinData>& bins,
  bool                        useIdeal  // false=HGTD, true=Ideal Res.
) {
  std::cout << "\n----------------------------------------------------------------\n";
  std::cout << header << "\n";
  std::cout << "  vs. n Forward HS Tracks (with valid HGTD time)\n";
  std::cout << "----------------------------------------------------------------\n";
  printf("%-24s %10s %10s %10s %12s\n",
         "Bin", "n_hasPure", "n_noPure", "Total", "Misclust (%)");
  std::cout << "----------------------------------------------------------------\n";

  double tot_pure = 0, tot_all = 0;
  for (int b = 0; b < NBINS; ++b) {
    double den    = bins[b].total;
    double n_pure = useIdeal ? bins[b].ideal_pure : bins[b].hgtd_pure;
    double n_no   = den - n_pure;
    tot_pure += n_pure;
    tot_all  += den;
    if (den == 0) continue;
    double rate = 100.0 * n_no / den;
    if (b < static_cast<int>(FOLD_HS_TRACK))
      printf("[%5.2f, %5.2f)  %10.1f %10.1f %10.1f %12.1f\n",
             static_cast<double>(b), static_cast<double>(b + 1),
             n_pure, n_no, den, rate);
    else
      printf("[%5.2f, %6.2f) %10.1f %10.1f %10.1f %12.1f\n",
             static_cast<double>(b), static_cast<double>(b + 1),
             n_pure, n_no, den, rate);
  }
  std::cout << "----------------------------------------------------------------\n";
  double tot_no = tot_all - tot_pure;
  printf("Total                  %10.1f %10.1f %10.1f %12.1f\n",
         tot_pure, tot_no, tot_all,
         tot_all > 0 ? 100.0 * tot_no / tot_all : 0.0);
  std::cout << "----------------------------------------------------------------\n";
}

// ---------------------------------------------------------------------------
// Misassignment table printer
//   Shows n_improved | Total | Rate (%)
//   Rate = events where Ideal Res. PASS passes but HGTD PASS fails / Total
// ---------------------------------------------------------------------------

void printMisassignmentTable(const std::vector<BinData>& bins) {
  std::cout << "\n----------------------------------------------------------------\n";
  std::cout << "Misassignment Rate  (HGTD |dt| - Ideal Res. |dt| > 60 ps, best-purity cluster)\n";
  std::cout << "  vs. n Forward HS Tracks (with valid HGTD time)\n";
  std::cout << "----------------------------------------------------------------\n";
  printf("%-24s %12s %10s %12s\n", "Bin", "n_improved", "Total", "Rate (%)");
  std::cout << "----------------------------------------------------------------\n";

  double tot_imp = 0, tot_all = 0;
  for (int b = 0; b < NBINS; ++b) {
    double den = bins[b].total;
    double imp = bins[b].misassign;
    tot_imp += imp;
    tot_all += den;
    if (den == 0) continue;
    double rate = 100.0 * imp / den;
    if (b < static_cast<int>(FOLD_HS_TRACK))
      printf("[%5.2f, %5.2f)  %12.1f %10.1f %12.1f\n",
             static_cast<double>(b), static_cast<double>(b + 1),
             imp, den, rate);
    else
      printf("[%5.2f, %6.2f) %12.1f %10.1f %12.1f\n",
             static_cast<double>(b), static_cast<double>(b + 1),
             imp, den, rate);
  }
  std::cout << "----------------------------------------------------------------\n";
  printf("Total                  %12.1f %10.1f %12.1f\n",
         tot_imp, tot_all,
         tot_all > 0 ? 100.0 * tot_imp / tot_all : 0.0);
  std::cout << "----------------------------------------------------------------\n";
}

// ---------------------------------------------------------------------------
// Misassignment table — no-misclustering subset (purple curve)
//   Denominator: events where HGTD has ≥1 cluster with purity > 0.75
//   Numerator:   those events where HGTD |dt| − ideal |dt| > 60 ps
// ---------------------------------------------------------------------------

void printMisassignmentPureTable(const std::vector<BinData>& bins) {
  std::cout << "\n----------------------------------------------------------------\n";
  std::cout << "Misassignment Rate  (no misclustering subset, purity > 75%)\n";
  std::cout << "  HGTD |dt| - Ideal Res. |dt| > 60 ps  |  denom = events with pure HGTD cluster\n";
  std::cout << "  vs. n Forward HS Tracks (with valid HGTD time)\n";
  std::cout << "----------------------------------------------------------------\n";
  printf("%-24s %12s %10s %12s\n", "Bin", "n_improved", "Total", "Rate (%)");
  std::cout << "----------------------------------------------------------------\n";

  double tot_imp = 0, tot_all = 0;
  for (int b = 0; b < NBINS; ++b) {
    double den = bins[b].total_pure;
    double imp = bins[b].misassign_pure;
    tot_imp += imp;
    tot_all += den;
    if (den == 0) continue;
    double rate = 100.0 * imp / den;
    if (b < static_cast<int>(FOLD_HS_TRACK))
      printf("[%5.2f, %5.2f)  %12.1f %10.1f %12.1f\n",
             static_cast<double>(b), static_cast<double>(b + 1), imp, den, rate);
    else
      printf("[%5.2f, %6.2f) %12.1f %10.1f %12.1f\n",
             static_cast<double>(b), static_cast<double>(b + 1), imp, den, rate);
  }
  std::cout << "----------------------------------------------------------------\n";
  printf("Total                  %12.1f %10.1f %12.1f\n",
         tot_imp, tot_all,
         tot_all > 0 ? 100.0 * tot_imp / tot_all : 0.0);
  std::cout << "----------------------------------------------------------------\n";
}

// ---------------------------------------------------------------------------
// plotRates — produce ../figs/rate_diagnostics.pdf (2 pages)
//   Page 1: Misclustering rate — HGTD (blue) and Ideal Res. (red) overlaid
//   Page 2: Misassignment rate  — Δ|dt| > 60 ps (orange)
//
//   Style matches the efficiency-vs-nHSTrack plots in plotting_utilities.h:
//   TEfficiency with kFNormal statistic, COLORS palette, linewidth 2,
//   legend at (0.55, 0.65, 0.9, 0.88), x-axis folded at FOLD_HS_TRACK.
// ---------------------------------------------------------------------------

void plotRates(TH1D* h_total,
               TH1D* h_hgtd_nopu, TH1D* h_ideal_nopu,
               TH1D* h_misassign,
               TH1D* h_total_pure, TH1D* h_misassign_pure)
{
  using namespace MyUtl;

  gStyle->SetOptStat(0);

  const double xMin = -0.5;
  const double xMax =  static_cast<double>(FOLD_HS_TRACK) + 0.5;  // 10.5
  const double yMin =  0.0;
  const double yMax =  1.1;  // headroom above 1.0 for legend

  // ── Build TEfficiency objects ──────────────────────────────────────────────
  auto* eff_hgtd     = new TEfficiency(*h_hgtd_nopu,      *h_total);
  auto* eff_ideal    = new TEfficiency(*h_ideal_nopu,     *h_total);
  auto* eff_mis      = new TEfficiency(*h_misassign,      *h_total);
  auto* eff_mis_pure = new TEfficiency(*h_misassign_pure, *h_total_pure);

  for (auto* e : {eff_hgtd, eff_ideal, eff_mis, eff_mis_pure})
    e->SetStatisticOption(TEfficiency::kFNormal);

  eff_hgtd    ->SetLineColor(COLORS[0]);  // blue   — HGTD
  eff_ideal   ->SetLineColor(COLORS[1]);  // red    — Ideal Res.
  eff_mis     ->SetLineColor(COLORS[6]);  // orange — all events
  eff_mis_pure->SetLineColor(COLORS[4]);  // purple — no-misclustering subset
  for (auto* e : {eff_hgtd, eff_ideal, eff_mis, eff_mis_pure}) e->SetLineWidth(2);

  // Titles follow ROOT's "title;xlabel;ylabel" convention.
  // These are also picked up by GetPaintedGraph() after the first Draw+Update.
  const char* xLabel = "n Forward HS Tracks (with valid HGTD time)";
  eff_hgtd    ->SetTitle(Form("Misclustering Rate vs n Forward HS Tracks;%s;Misclustering Rate", xLabel));
  eff_ideal   ->SetTitle(Form("Misclustering Rate vs n Forward HS Tracks;%s;Misclustering Rate", xLabel));
  eff_mis     ->SetTitle(Form("Misassignment Rate vs n Forward HS Tracks;%s;Misassignment Rate", xLabel));
  eff_mis_pure->SetTitle(Form("Misassignment Rate vs n Forward HS Tracks;%s;Misassignment Rate", xLabel));

  // ── Helper: draw scaled total-events distribution in the bottom ~15% ────
  //   Mirrors the pattern in plotting_utilities.h::moneyPlot().
  //   A unique clone name is required each call to avoid ROOT name collisions.
  int _distCallCount = 0;
  auto drawTotalDist = [&](double yMaxVal) {
    auto* htot = static_cast<TH1D*>(
      h_total->Clone(Form("h_total_disp_%d", _distCallCount++)));
    htot->SetDirectory(nullptr);
    htot->SetLineColorAlpha(COLORS[8], 0.6);
    htot->SetFillColorAlpha(COLORS[8], 0.3);
    htot->SetLineWidth(2);
    if (htot->Integral() > 0) {
      htot->Scale(1.0 / htot->Integral());
      htot->Scale(0.15 * yMaxVal / htot->GetMaximum());
    }
    htot->Draw("HIST SAME");
  };

  // ── Helper: fix axes after gPad::Update() paints the underlying TGraph ───
  auto styleAxes = [&](TEfficiency* e, const char* ytitle, double yMaxVal) {
    auto* g = e->GetPaintedGraph();
    g->SetTitle(e->GetTitle());  // propagate title to the painted graph
    g->GetXaxis()->SetTitle(xLabel);
    g->GetXaxis()->SetLimits(xMin, xMax);
    g->GetXaxis()->SetRangeUser(xMin, xMax);
    g->GetXaxis()->SetNdivisions(510);
    g->GetYaxis()->SetTitle(ytitle);
    g->GetYaxis()->SetRangeUser(yMin, yMaxVal);
  };

  TCanvas* c = new TCanvas("c_rates", "Rate diagnostics", 800, 600);
  const char* fname = "../figs/rate_diagnostics.pdf";

  // ── Page 1: Misclustering rate ────────────────────────────────────────────
  eff_hgtd->Draw("AP");
  gPad->Update();
  styleAxes(eff_hgtd, "Misclustering Rate", yMax);
  gPad->Update();

  eff_ideal->Draw("P SAME");
  gPad->Update();

  drawTotalDist(yMax);
  gPad->Update();

  auto* leg1 = new TLegend(0.55, 0.65, 0.9, 0.88);
  leg1->AddEntry(eff_hgtd,  "HGTD",      "l");
  leg1->AddEntry(eff_ideal, "Ideal Res.", "l");
  leg1->Draw();

  c->Print(Form("%s[", fname));  // open multi-page PDF
  c->Print(fname);
  c->Clear();

  // ── Page 2: Misassignment rate (all events + no-misclustering subset) ───────
  eff_mis->Draw("AP");
  gPad->Update();
  styleAxes(eff_mis, "Misassignment Rate", 0.5);
  gPad->Update();

  eff_mis_pure->Draw("P SAME");
  gPad->Update();

  drawTotalDist(0.5);
  gPad->Update();

  auto* leg2 = new TLegend(0.5, 0.65, 0.9, 0.88);
  leg2->AddEntry(eff_mis,      "All events",                    "l");
  leg2->AddEntry(eff_mis_pure, "No misclustering (purity > 75%)", "l");
  leg2->Draw();

  c->Print(fname);
  c->Print(Form("%s]", fname));  // close PDF

  std::cout << "Saved plots to " << fname << "\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

auto main() -> int {
  // ── Data source ──────────────────────────────────────────────────────────
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = kFatal;

  // ── Accumulators ─────────────────────────────────────────────────────────
  std::vector<BinData> bins(NBINS);

  // ── Histograms for plotting (11 bins, one per folded nHSTrack value 0-10) ─
  auto makeH = [](const char* name) {
    return new TH1D(name, "", 11, -0.5, 10.5);
  };
  TH1D* h_total           = makeH("h_total");
  TH1D* h_hgtd_nopu       = makeH("h_hgtd_nopu");       // no-pure-cluster events (HGTD)
  TH1D* h_ideal_nopu      = makeH("h_ideal_nopu");      // no-pure-cluster events (Ideal Res.)
  TH1D* h_misassign       = makeH("h_misassign");        // misassignment events (all)
  TH1D* h_total_pure      = makeH("h_total_pure");      // denom: events with HGTD pure cluster (no misclustering)
  TH1D* h_misassign_pure  = makeH("h_misassign_pure");  // misassignment events | no misclustering

  // ── filterClusters: keep clusters with >= MIN_CLUSTER_TRACKS tracks ──────
  auto filterClusters = [](const std::vector<Cluster>& col) {
    std::vector<Cluster> out;
    for (const auto& c : col)
      if (c.nConstituents >= MIN_CLUSTER_TRACKS) out.push_back(c);
    return out;
  };

  // ── Event loop ───────────────────────────────────────────────────────────
  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  while (reader.Next()) {
    const Long64_t READ_NUM = chain.GetReadEntry() + 1;
    if (READ_NUM % 100 == 0)
      std::cout << "Progress: " << READ_NUM << "/" << N_EVENT << "\r" << std::flush;

    // ── A. Event selection ────────────────────────────────────────────────
    if (!branch.passBasicCuts()) continue;
    if (!branch.passJetPtCut())  continue;

    // ── B. Track selection ────────────────────────────────────────────────
    std::vector<int> tracks =
      getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);

    // ── C. Bin index: n HS tracks with valid HGTD time ────────────────────
    int nHSTrack = 0;
    for (int idx : tracks)
      if (branch.trackToTruthvtx[idx] == 0 && branch.trackTimeValid[idx] == 1)
        ++nHSTrack;
    const int bin = folded(nHSTrack, static_cast<int>(FOLD_HS_TRACK));

    // ── D. HGTD clustering (real times, purity enabled) ───────────────────
    auto hgtdRaw = clusterTracksInTime(
      tracks, &branch, DIST_CUT_CONE,
      /*useSmearedTimes=*/false, /*checkTimeValid=*/true, /*smearRes=*/-1,
      ClusteringMethod::CONE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);

    // ── E. Ideal Res. clustering (truth-centred 10 ps smear, purity) ──────
    auto idealRaw = clusterTracksInTime(
      tracks, &branch, DIST_CUT_CONE,
      /*useSmearedTimes=*/true, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
      ClusteringMethod::CONE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);

    const auto qualHGTD  = filterClusters(hgtdRaw);
    const auto qualIdeal = filterClusters(idealRaw);

    // ── F. Compute per-event flags ────────────────────────────────────────
    bool hgtd_hasPure  = false, ideal_hasPure = false;

    // For misassignment: track best (highest-purity) qualifying cluster per scenario
    // and its |cluster_time − truth_time|.  Negative sentinel means no cluster found.
    const float truthTime      = branch.truthVtxTime[0];
    float best_hgtd_purity     = -1.f,  best_ideal_purity = -1.f;
    float best_hgtd_dt         = -1.f,  best_ideal_dt     = -1.f;

    for (const auto& c : qualHGTD) {
      if (c.purity > 0.75f) hgtd_hasPure = true;
      if (c.purity > best_hgtd_purity) {
        best_hgtd_purity = c.purity;
        best_hgtd_dt     = std::abs(c.values[0] - truthTime);
      }
    }
    for (const auto& c : qualIdeal) {
      if (c.purity > 0.75f) ideal_hasPure = true;
      if (c.purity > best_ideal_purity) {
        best_ideal_purity = c.purity;
        best_ideal_dt     = std::abs(c.values[0] - truthTime);
      }
    }

    // ── G. Accumulate ─────────────────────────────────────────────────────
    bins[bin].total++;
    if (hgtd_hasPure)  bins[bin].hgtd_pure++;
    if (ideal_hasPure) bins[bin].ideal_pure++;
    // Misassignment: ideal timing reduces best-purity cluster |dt| by > 60 ps.
    // Require both scenarios have ≥1 qualifying cluster so the comparison is fair.
    const bool isMisassign =
      best_hgtd_dt >= 0.f && best_ideal_dt >= 0.f &&
      (best_hgtd_dt - best_ideal_dt) > 60.f;
    if (isMisassign) bins[bin].misassign++;

    // ── H. Fill histograms ────────────────────────────────────────────────
    h_total->Fill(bin);
    if (!hgtd_hasPure)  h_hgtd_nopu->Fill(bin);
    if (!ideal_hasPure) h_ideal_nopu->Fill(bin);
    if (isMisassign)    h_misassign->Fill(bin);
    // Conditioned on no misclustering: HGTD has ≥1 cluster with purity > 0.75
    if (hgtd_hasPure) {
      h_total_pure->Fill(bin);
      if (isMisassign) h_misassign_pure->Fill(bin);
      bins[bin].total_pure++;
      if (isMisassign) bins[bin].misassign_pure++;
    }
  }

  std::cout << "\nFINISHED PROCESSING\n";

  // ── Print tables ──────────────────────────────────────────────────────────
  printMisclusteringTable(
    "Misclustering Rate — HGTD  (any cluster purity > 75%)", bins, false);
  printMisclusteringTable(
    "Misclustering Rate — Ideal Res.  (any cluster purity > 75%)", bins, true);
  printMisassignmentTable(bins);
  printMisassignmentPureTable(bins);

  plotRates(h_total, h_hgtd_nopu, h_ideal_nopu, h_misassign,
            h_total_pure, h_misassign_pure);

  return 0;
}
