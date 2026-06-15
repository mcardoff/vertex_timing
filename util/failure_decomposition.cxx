// ---------------------------------------------------------------------------
// failure_decomposition.cxx
//
// For every event where HGTD TRKPTZ fails the efficiency test, classifies
// the failure into four (potentially overlapping) diagnostic categories:
//
//   1. Selection failure    — another qualifying cluster in the same HGTD
//                             cone-clustering collection passes
//                             passEfficiency(); TRKPTZ chose suboptimally.
//
//   2. Timing misassignment — the Ideal Resolution scenario (flat
//                             IDEAL_TRACK_RES per-track smear, real track
//                             efficiency: checkTimeValid=true) has TRKPTZ
//                             pass.  Fixing per-track timing resolves it.
//
//   3. Misclustering        — no qualifying HGTD cluster achieves purity
//                             > 75%; no selection algorithm could recover
//                             this event regardless of which cluster is
//                             chosen.
//
//   4. Track efficiency loss — the Ideal Res. + Ideal Eff. scenario
//                             (checkTimeValid=false) has TRKPTZ pass, but
//                             Ideal Res. alone does NOT.  Giving currently
//                             time-less HS tracks an ideal smeared time is
//                             what fixes the event.
//
// Note on construction:
//   Category 4 is defined as "Ideal+Eff. fixes it AND Ideal Res. does not",
//   so categories 2 and 4 are mutually exclusive by construction.  All other
//   pairwise overlaps are measured and printed.
//
// Output:
//   • Summary table over all events (counts, % of fail, % of total)
//   • Pairwise and triple overlap counts
//   • Exclusive counts (only one category applies)
//   • Per-bin breakdown vs n Forward HS Tracks (with valid HGTD time)
// ---------------------------------------------------------------------------

#include <TROOT.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPie.h>
#include <TStyle.h>
#include "clustering_constants.h"
#include "event_processing.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Per-bin accumulator
// ---------------------------------------------------------------------------
struct BinData {
  // Event totals
  double total      = 0;   // events passing event selection
  double n_pass     = 0;   // TRKPTZ passes efficiency
  double n_fail     = 0;   // TRKPTZ fails efficiency (our population)

  // The four diagnostic flags (among failing events)
  double cat_sel    = 0;   // 1. Another cluster passes → selection failure
  double cat_timing = 0;   // 2. Ideal Res. TRKPTZ passes → timing misassignment
  double cat_miscl  = 0;   // 3. No pure cluster (>75%) OR REFINED fixes it → misclustering
  double cat_teff   = 0;   // 4. Ideal+Eff. fixes it, Ideal Res. alone doesn't

  double cat_none   = 0;   // None of the above

  // Pairwise overlaps (2∩4 and any triple/quad involving 2∩4 are always 0)
  double ov_s_t     = 0;   // 1 ∩ 2
  double ov_s_m     = 0;   // 1 ∩ 3
  double ov_s_e     = 0;   // 1 ∩ 4
  double ov_t_m     = 0;   // 2 ∩ 3
  double ov_m_e     = 0;   // 3 ∩ 4

  // Triple overlaps
  double ov_stm     = 0;   // 1 ∩ 2 ∩ 3
  double ov_sme     = 0;   // 1 ∩ 3 ∩ 4
};

static const int NBINS = static_cast<int>(FOLD_HS_TRACK) + 1;

// ---------------------------------------------------------------------------
// chooseTRKPTZ / chooseWAVeS — return pointer to the highest-score qualifying
//                               cluster, or nullptr if the collection is empty.
// ---------------------------------------------------------------------------
static const Cluster* chooseTRKPTZ(const std::vector<Cluster>& qual) {
  if (qual.empty()) return nullptr;
  const Cluster* best = &qual[0];
  double maxScore = best->scores.at(Score::TRKPTZ.id);
  for (const auto& c : qual) {
    double s = c.scores.at(Score::TRKPTZ.id);
    if (s > maxScore) { maxScore = s; best = &c; }
  }
  return best;
}

static const Cluster* chooseWAVeS(const std::vector<Cluster>& qual) {
  if (qual.empty()) return nullptr;
  const Cluster* best = &qual[0];
  double maxScore = best->scores.at(Score::WAVES.id);
  for (const auto& c : qual) {
    double s = c.scores.at(Score::WAVES.id);
    if (s > maxScore) { maxScore = s; best = &c; }
  }
  return best;
}

// ---------------------------------------------------------------------------
// Aggregate bins into totals.  Returns a single BinData summed over all bins.
// ---------------------------------------------------------------------------
static BinData aggregate(const std::vector<BinData>& bins) {
  BinData t;
  for (const auto& b : bins) {
    t.total     += b.total;
    t.n_pass    += b.n_pass;
    t.n_fail    += b.n_fail;
    t.cat_sel   += b.cat_sel;
    t.cat_timing += b.cat_timing;
    t.cat_miscl += b.cat_miscl;
    t.cat_teff  += b.cat_teff;
    t.cat_none  += b.cat_none;
    t.ov_s_t    += b.ov_s_t;
    t.ov_s_m    += b.ov_s_m;
    t.ov_s_e    += b.ov_s_e;
    t.ov_t_m    += b.ov_t_m;
    t.ov_m_e    += b.ov_m_e;
    t.ov_stm    += b.ov_stm;
    t.ov_sme    += b.ov_sme;
  }
  return t;
}

// ---------------------------------------------------------------------------
// printSummary — full summary table for a BinData aggregate
// ---------------------------------------------------------------------------
static void printSummary(const BinData& g, const char* label = "TRKPTZ") {
  const double nf = g.n_fail;
  const double nt = g.total;
  auto pf = [&](double n) { return nf > 0 ? 100.*n/nf : 0.; };  // % of fail
  auto pt = [&](double n) { return nt > 0 ? 100.*n/nt : 0.; };  // % of total

  // Exclusive counts via inclusion-exclusion (categories 2 and 4 are disjoint)
  //   |only A| = |A| - |A∩B| - |A∩C| - |A∩D| + |A∩B∩C| + |A∩C∩D|
  double excl_sel    = g.cat_sel  - g.ov_s_t - g.ov_s_m - g.ov_s_e + g.ov_stm + g.ov_sme;
  double excl_timing = g.cat_timing - g.ov_s_t - g.ov_t_m + g.ov_stm;
  double excl_miscl  = g.cat_miscl  - g.ov_s_m - g.ov_t_m - g.ov_m_e + g.ov_stm + g.ov_sme;
  double excl_teff   = g.cat_teff   - g.ov_s_e - g.ov_m_e + g.ov_sme;

  std::cout << "\n";
  std::cout << "================================================================\n";
  printf("  %s Failure Decomposition\n", label);
  std::cout << "================================================================\n";
  printf("  Total events (passing selection): %9.0f\n", g.total);
  printf("  %s pass:                      %9.0f  (%5.2f%%)\n", label, g.n_pass, pt(g.n_pass));
  printf("  %s fail:                      %9.0f  (%5.2f%%)\n", label, g.n_fail, pt(g.n_fail));

  std::cout << "\n  --- Among failing events (N = " << (int)g.n_fail << ") ---\n\n";
  printf("  %-45s %8s %9s %9s\n",
         "Category", "Count", "%fail", "%total");
  std::cout << "  " << std::string(73, '-') << "\n";
  printf("  %-45s %8.0f %8.1f%% %8.1f%%\n",
         "1. Selection failure (another cluster passes)",
         g.cat_sel, pf(g.cat_sel), pt(g.cat_sel));
  printf("  %-45s %8.0f %8.1f%% %8.1f%%\n",
         "2. Timing misassignment (Ideal Res. fixes)",
         g.cat_timing, pf(g.cat_timing), pt(g.cat_timing));
  printf("  %-45s %8.0f %8.1f%% %8.1f%%\n",
         "3. Misclustering (no pure >75% OR REFINED fixes)",
         g.cat_miscl, pf(g.cat_miscl), pt(g.cat_miscl));
  printf("  %-45s %8.0f %8.1f%% %8.1f%%\n",
         "4. Track efficiency loss (Ideal+Eff. fixes)",
         g.cat_teff, pf(g.cat_teff), pt(g.cat_teff));
  std::cout << "  " << std::string(73, '-') << "\n";
  printf("  %-45s %8.0f %8.1f%% %8.1f%%\n",
         "None of the above (fundamental failures)",
         g.cat_none, pf(g.cat_none), pt(g.cat_none));

  std::cout << "\n  --- Pairwise overlaps ---\n\n";
  printf("  %-45s %8.0f %8.1f%%\n", "1 ∩ 2  Selection ∩ Timing",       g.ov_s_t, pf(g.ov_s_t));
  printf("  %-45s %8.0f %8.1f%%\n", "1 ∩ 3  Selection ∩ Miscl",        g.ov_s_m, pf(g.ov_s_m));
  printf("  %-45s %8.0f %8.1f%%\n", "1 ∩ 4  Selection ∩ TrkEff",       g.ov_s_e, pf(g.ov_s_e));
  printf("  %-45s %8.0f %8.1f%%\n", "2 ∩ 3  Timing ∩ Miscl",           g.ov_t_m, pf(g.ov_t_m));
  printf("  %-45s %8.0f %8.1f%%\n", "2 ∩ 4  Timing ∩ TrkEff",          0.0,      0.0);  // always 0
  printf("  %-45s %8.0f %8.1f%%\n", "3 ∩ 4  Miscl ∩ TrkEff",           g.ov_m_e, pf(g.ov_m_e));

  std::cout << "\n  --- Triple overlaps ---\n\n";
  printf("  %-45s %8.0f %8.1f%%\n", "1 ∩ 2 ∩ 3  Selection ∩ Timing ∩ Miscl", g.ov_stm, pf(g.ov_stm));
  printf("  %-45s %8.0f %8.1f%%\n", "1 ∩ 3 ∩ 4  Selection ∩ Miscl ∩ TrkEff", g.ov_sme, pf(g.ov_sme));
  printf("  %-45s %8.0f %8.1f%%\n", "1 ∩ 2 ∩ 4  (impossible: 2∩4=∅)",         0.0,      0.0);
  printf("  %-45s %8.0f %8.1f%%\n", "2 ∩ 3 ∩ 4  (impossible: 2∩4=∅)",         0.0,      0.0);

  std::cout << "\n  --- Exclusive contributions (only this category) ---\n\n";
  printf("  %-45s %8.0f %8.1f%%\n", "Only Selection failure",      excl_sel,    pf(excl_sel));
  printf("  %-45s %8.0f %8.1f%%\n", "Only Timing misassignment",   excl_timing, pf(excl_timing));
  printf("  %-45s %8.0f %8.1f%%\n", "Only Misclustering",          excl_miscl,  pf(excl_miscl));
  printf("  %-45s %8.0f %8.1f%%\n", "Only Track efficiency loss",  excl_teff,   pf(excl_teff));

  // Sanity: |A∪B∪C∪D| + |None| should equal n_fail.
  // Inclusion-exclusion with |B∩D|=0 (categories 2 and 4 are mutually exclusive):
  //   |A∪B∪C∪D| = Σ|single| - Σ|pairs| + Σ|triples|
  double union_abcd = g.cat_sel + g.cat_timing + g.cat_miscl + g.cat_teff
                    - g.ov_s_t - g.ov_s_m - g.ov_s_e - g.ov_t_m - g.ov_m_e
                    + g.ov_stm + g.ov_sme;
  double check = union_abcd + g.cat_none;
  printf("\n  Sanity check (should equal n_fail = %.0f):  %.0f\n", g.n_fail, check);
  std::cout << "================================================================\n";
}

// ---------------------------------------------------------------------------
// printPerBinTable — show per-bin failure rates and category breakdowns
// ---------------------------------------------------------------------------
static void printPerBinTable(const std::vector<BinData>& bins) {
  const BinData g = aggregate(bins);

  std::cout << "\n";
  std::cout << "---";
  std::cout << std::string(105, '-') << "\n";
  std::cout << "  Per-bin breakdown vs n Forward HS Tracks (with valid HGTD time)\n";
  std::cout << "  Category percentages are expressed as % of failing events in that bin.\n";
  std::cout << "---";
  std::cout << std::string(105, '-') << "\n";
  printf("  %-14s %7s %7s %7s | %9s %9s %9s %9s | %9s\n",
         "Bin", "Total", "Fail", "Fail%",
         "Sel%", "Timing%", "Miscl%", "TrkEff%", "None%");
  std::cout << "---";
  std::cout << std::string(105, '-') << "\n";

  auto rowFmt = [](const char* label, const BinData& b) {
    double fp   = b.total > 0  ? 100.*b.n_fail/b.total    : 0.;
    double selp = b.n_fail > 0 ? 100.*b.cat_sel/b.n_fail  : 0.;
    double timp = b.n_fail > 0 ? 100.*b.cat_timing/b.n_fail: 0.;
    double mclp = b.n_fail > 0 ? 100.*b.cat_miscl/b.n_fail: 0.;
    double tefp = b.n_fail > 0 ? 100.*b.cat_teff/b.n_fail : 0.;
    double nonp = b.n_fail > 0 ? 100.*b.cat_none/b.n_fail : 0.;
    printf("  %-14s %7.0f %7.0f %6.1f%% | %8.1f%% %8.1f%% %8.1f%% %8.1f%% | %8.1f%%\n",
           label, b.total, b.n_fail, fp,
           selp, timp, mclp, tefp, nonp);
  };

  for (int b = 0; b < NBINS; ++b) {
    if (bins[b].total == 0) continue;
    char label[32];
    if (b < (int)FOLD_HS_TRACK)
      std::snprintf(label, sizeof(label), "[%d, %d)", b, b+1);
    else
      std::snprintf(label, sizeof(label), "[%d, inf)", b);
    rowFmt(label, bins[b]);
  }

  std::cout << "---";
  std::cout << std::string(105, '-') << "\n";
  rowFmt("Total", g);
  std::cout << "---";
  std::cout << std::string(105, '-') << "\n";
}

// ---------------------------------------------------------------------------
// plotFailures — produce ../figs/failure_decomposition.pdf (2 pages)
//
//   Page 1: Overlaid TEfficiency rate curves (% of all events) vs nHSTrack.
//     Five curves: total failure (gray dashed), misclustering (blue),
//     timing misassignment (red), selection failure (green), track
//     efficiency loss (violet).  Same style as the money plots.
//
//   Page 2: Stacked bar chart showing the exclusive decomposition of the
//     failing events per bin.  Each bin sums to 1.0 (100% of failures).
//     Six exclusive groups, bottom to top:
//       1. Two-way overlaps (orange) — exactly two categories, ~57%
//       2. Three-way overlaps (cyan) — all three main categories, ~26%
//       3. Only misclustering (blue)
//       4. Only timing misassignment (red)
//       5. Only selection failure (green)
//       6. Only track efficiency loss (violet, tiny)
// ---------------------------------------------------------------------------
static void plotFailures(const std::vector<BinData>& bins) {
  gStyle->SetOptStat(0);

  const double xMin  = -0.5;
  const double xMax  = static_cast<double>(FOLD_HS_TRACK) + 0.5;  // 10.5
  const char*  xLbl  = "n Forward HS Tracks (with valid HGTD time)";
  const char*  fname = "../figs/failure_decomposition.pdf";

  // ── Build histograms from BinData ─────────────────────────────────────────
  auto makeH = [](const char* name) -> TH1D* {
    auto* h = new TH1D(name, "", 11, -0.5, 10.5);
    h->Sumw2(false);
    return h;
  };

  // Page 1: TEfficiency denominators / numerators (bin content = event count)
  TH1D* h_total  = makeH("hf_total");
  TH1D* h_fail   = makeH("hf_fail");
  TH1D* h_sel    = makeH("hf_sel");
  TH1D* h_timing = makeH("hf_timing");
  TH1D* h_miscl  = makeH("hf_miscl");
  TH1D* h_teff   = makeH("hf_teff");

  for (int b = 0; b < NBINS; ++b) {
    const auto& B  = bins[b];
    const int   rb = b + 1;  // ROOT bin index (1-based)

    h_total ->SetBinContent(rb, B.total);
    h_fail  ->SetBinContent(rb, B.n_fail);
    h_sel   ->SetBinContent(rb, B.cat_sel);
    h_timing->SetBinContent(rb, B.cat_timing);
    h_miscl ->SetBinContent(rb, B.cat_miscl);
    h_teff  ->SetBinContent(rb, B.cat_teff);

  }

  // ── Axis helper (applied after TEfficiency is drawn + gPad updated) ───────
  auto styleAxes = [&](TEfficiency* e, const char* ytitle, double ymax) {
    auto* g = e->GetPaintedGraph();
    g->GetXaxis()->SetTitle(xLbl);
    g->GetXaxis()->SetLimits(xMin, xMax);
    g->GetXaxis()->SetRangeUser(xMin, xMax);
    g->GetXaxis()->SetNdivisions(510);
    g->GetYaxis()->SetTitle(ytitle);
    g->GetYaxis()->SetRangeUser(0.5e-3, ymax);
  };

  // ── Scaled total distribution overlay (bottom 15%), same as money plots ───
  int _nc = 0;
  auto drawTotDist = [&](double ymax) {
    auto* hd = static_cast<TH1D*>(h_total->Clone(Form("hf_tot_disp_%d", _nc++)));
    hd->SetDirectory(nullptr);
    hd->SetLineColorAlpha(COLORS[8], 0.6);
    hd->SetFillColorAlpha(COLORS[8], 0.3);
    hd->SetLineWidth(2);
    if (hd->Integral() > 0) {
      hd->Scale(1.0 / hd->Integral());
      hd->Scale(0.15 * ymax / hd->GetMaximum());
    }
    hd->Draw("HIST SAME");
  };

  // ── Build TEfficiency objects for page 1 ──────────────────────────────────
  auto* eff_fail   = new TEfficiency(*h_fail,   *h_total);  // overall rate (/ all events)
  auto* eff_miscl  = new TEfficiency(*h_miscl,  *h_fail);   // fraction of failing events
  auto* eff_timing = new TEfficiency(*h_timing, *h_fail);
  auto* eff_sel    = new TEfficiency(*h_sel,    *h_fail);
  auto* eff_teff   = new TEfficiency(*h_teff,   *h_fail);

  // for (auto* e : {eff_fail, eff_miscl, eff_timing, eff_sel, eff_teff})
    // e->SetStatisticOption(TEfficiency::kFNormal);

  // Total failure: gray dashed
  eff_fail  ->SetLineColor(COLORS[3]);  eff_fail  ->SetLineStyle(7); eff_fail  ->SetLineWidth(3);
  // Categories: consistent colours with page 2 stacked bars
  eff_miscl ->SetLineColor(COLORS[0]);  eff_miscl ->SetLineWidth(2);  // blue
  eff_timing->SetLineColor(COLORS[1]);  eff_timing->SetLineWidth(2);  // red
  eff_sel   ->SetLineColor(COLORS[7]);  eff_sel   ->SetLineWidth(2);  // green
  eff_teff  ->SetLineColor(COLORS[4]);  eff_teff  ->SetLineWidth(2);  // violet

  eff_fail->SetTitle(Form(
    "TRKPTZ Failure Rate Decomposition;%s;Rate (fraction of failing events)", xLbl));

  // ── Canvas: 3:1 split (top 75% = log-scale rates, bottom 25% = nHS dist.) ────
  TCanvas* c = new TCanvas("c_fail_decomp", "Failure decomposition", 800, 700);

  // Explicit pads — identical left/right margins so y-axes align
  const double split = 0.25;  // fraction given to bottom pad
  const double lm = 0.12, rm = 0.05;

  auto* padTop = new TPad("padTop", "", 0, split, 1, 1);
  padTop->SetTopMargin(0.12);
  padTop->SetBottomMargin(0.00);
  padTop->SetLeftMargin(lm);
  padTop->SetRightMargin(rm);
  padTop->Draw();

  auto* padBot = new TPad("padBot", "", 0, 0, 1, split);
  padBot->SetTopMargin(0.00);
  padBot->SetBottomMargin(0.38);  // proportionally large since pad is short
  padBot->SetLeftMargin(lm);
  padBot->SetRightMargin(rm);
  padBot->Draw();

  // ── Top pad: TEfficiency curves with log Y-scale ──────────────────────────────
  padTop->cd();
  // padTop->SetLogy(true);

  const double yPage1 = 1.5;

  eff_fail->Draw("AP");      gPad->Update();

  for (auto* e : {eff_miscl, eff_timing, eff_sel, eff_teff}) {
    e->Draw("P SAME"); gPad->Update();
  }

  // drawTotDist(yPage1);

  auto* leg1 = new TLegend(0.60, 0.65, 0.95, 0.88);
  leg1->AddEntry(eff_fail,   "Overall failure rate (/ all events)",      "l");
  leg1->AddEntry(eff_sel,    "1. Selection failure (other cluster)",     "l");
  leg1->AddEntry(eff_timing, "2. Timing misassignment (Ideal Res.)",     "l");
  leg1->AddEntry(eff_miscl,  "3. Misclustering (no pure >75% OR REFINED)", "l");
  leg1->AddEntry(eff_teff,   "4. Track efficiency loss",                 "l");
  leg1->Draw();
  gPad->Update();  // final update before axis customisation

  // Apply axis styling LAST — per ROOT docs, GetPaintedGraph()->GetXaxis()->
  // SetRangeUser() must be called after Draw+Update and NOT followed by
  // another Update() or ROOT auto-expands the TGraph frame again.
  
  gPad->Update();
  styleAxes(eff_fail,   "Rate (fraction of failing events)", yPage1);
  styleAxes(eff_miscl,  "Rate (fraction of failing events)", yPage1);
  styleAxes(eff_timing, "Rate (fraction of failing events)", yPage1);
  styleAxes(eff_sel,    "Rate (fraction of failing events)", yPage1);
  styleAxes(eff_teff,   "Rate (fraction of failing events)", yPage1);
  // Suppress x-axis on top pad — shared axis shown on bottom pad only
  eff_fail->GetPaintedGraph()->GetXaxis()->SetLabelSize(0);
  eff_fail->GetPaintedGraph()->GetXaxis()->SetTitleSize(0);
  gPad->Update();
  padTop->SetLogy(false);

  // ── Bottom pad: nHS distribution with linear Y-scale ──────────────────────────
  padBot->cd();

  // Scale font sizes up to compensate for the shorter pad (factor ~3)
  const double sf = 1.0 / split;  // = 4.0 (NDC font sizes scale with pad height)
  h_total->SetLineColor(COLORS[8]);
  h_total->SetFillColor(COLORS[8]);
  h_total->SetFillStyle(3004);  // diagonal hatch
  h_total->SetLineWidth(2);

  // Use DrawFrame to force EXACT x-range matching the top pad's TEfficiency frame.
  // Axis styling goes on the frame, then h_total is drawn on top with "HIST SAME".
  double yBotMax = h_total->GetMaximum() * 1.15;
  auto* hframe = static_cast<TH1*>(gPad->DrawFrame(xMin, 0, xMax, yBotMax));
  hframe->GetXaxis()->SetTitle(xLbl);
  hframe->GetXaxis()->SetNdivisions(510);
  hframe->GetXaxis()->SetLabelSize(0.035 * sf);
  hframe->GetXaxis()->SetTitleSize(0.035 * sf);
  hframe->GetXaxis()->SetTitleOffset(0.85);
  hframe->GetYaxis()->SetTitle("Events");
  hframe->GetYaxis()->SetLabelSize(0.030 * sf);
  hframe->GetYaxis()->SetTitleSize(0.030 * sf);
  hframe->GetYaxis()->SetTitleOffset(0.4);
  hframe->GetYaxis()->SetNdivisions(504);  // fewer y-axis ticks
  gPad->Update();
  h_total->Draw("HIST SAME");
  gPad->Update();

  // ── Print page 1 (both pads together)
  c->cd();
  c->Print(Form("%s[", fname));
  c->Print(fname);
  c->Clear();

  // ── Page 2: pie chart — named exclusive decomposition ─────────────────────
  const BinData g = aggregate(bins);

  // All mutually-exclusive intersection pieces (sum = g.n_fail)
  double excl_sel    = g.cat_sel    - g.ov_s_t - g.ov_s_m - g.ov_s_e + g.ov_stm + g.ov_sme;
  double excl_timing = g.cat_timing - g.ov_s_t - g.ov_t_m             + g.ov_stm;
  double excl_miscl  = g.cat_miscl  - g.ov_s_m - g.ov_t_m - g.ov_m_e + g.ov_stm + g.ov_sme;
  double excl_teff   = g.cat_teff   - g.ov_s_e             - g.ov_m_e             + g.ov_sme;
  double two_s_t     = g.ov_s_t - g.ov_stm;
  double two_s_m     = g.ov_s_m - g.ov_stm - g.ov_sme;
  double two_t_m     = g.ov_t_m - g.ov_stm;
  double two_s_e     = g.ov_s_e - g.ov_sme;
  double two_m_e     = g.ov_m_e - g.ov_sme;
  double three_stm   = g.ov_stm;
  double three_sme   = g.ov_sme;

  // Original ordering is preserved except excl_teff (Pure Track Eff., ~0%) is moved
  // from index 3 (where it was directly adjacent to cat_none/None) to index 11 (end).
  // This puts ~274° of arc between the None and Pure Track Eff. labels so they
  // no longer overlap. The SetEntryRadiusOffset pops the tiny slice outward.
  const int NS = 12;
  double pieVals[NS] = {
    excl_sel,    // 0   Pure Selection              3.9%
    excl_timing, // 1   Pure Misassignment          5.4%
    excl_miscl,  // 2   Pure Misclustering          7.4%
    two_s_t,     // 4   Sel. ∩ Misassign            0.2%
    two_s_m,     // 5   Sel. ∩ Miscluster          19.2%
    two_t_m,     // 6   Misassign ∩ Miscluster      9.7%
    two_s_e,     // 7   Sel. ∩ Track Eff.          19.1%
    two_m_e,     // 8   Miscluster ∩ Track Eff.     3.2%
    three_stm,   // 9   Sel. ∩ Miassign ∩ Miscl    5.8%
    g.cat_none,  // 3   None                        6.9%  (excl_teff removed from here)
    three_sme,   // 10  Sel. ∩ Miscluster ∩ TrkEff 18.9%
    excl_teff,   // 11  Pure Track Eff.            ~0%   ← moved to end, now 274° from None
  };

  // ── Verify pie slices sum to n_fail; clamp any negative exclusive pieces ──
  {
    double pieSum = 0;
    for (int i = 0; i < NS; ++i) pieSum += pieVals[i];
    printf("Pie chart check: sum of slices = %.0f  (n_fail = %.0f)%s\n",
           pieSum, g.n_fail,
           std::abs(pieSum - g.n_fail) < 0.5 ? " OK" : "  *** MISMATCH ***");
    for (int i = 0; i < NS; ++i) {
      if (pieVals[i] < 0) {
        printf("  WARNING: slice %d clamped from %.2f to 0\n", i, pieVals[i]);
        pieVals[i] = 0;
      }
    }
  }

  // Each slice gets its own distinct color from the standard P10/P6 palette.
  Int_t pieCols[NS] = {
    COLORS[7],   //  0  excl_sel    — green
    COLORS[1],   //  1  excl_timing — red
    COLORS[0],   //  2  excl_miscl  — blue
    COLORS[5],   //  4  two_s_t     — brown
    COLORS[2],   //  5  two_s_m     — yellow
    COLORS[9],   //  6  two_t_m     — cyan
    COLORS[3],   //  7  two_s_e     — gray
    COLORS[10],  //  8  two_m_e     — P6Blue
    COLORS[8],   //  9  three_stm   — ash
    COLORS[6],   //  3  cat_none    — orange
    (Int_t)kGray,// 10  three_sme   — medium gray
    COLORS[4],   // 11  excl_teff   — violet
  };
  auto pct = [&](double v) { return Form("%.1f%%", 100.*v/g.n_fail); };
  const char* pieLabels[NS] = {
    "Pure Selection",
    "Pure Misassignment",
    "Pure Misclustering",
    "Sel. #cap Misassign",
    "Sel. #cap Miscluster",
    "Misassign #cap Miscluster",
    "Sel. #cap Track Eff.",
    "Miscluster #cap Track Eff.",
    "Sel. #cap Miassign #cap Miscluster",
    "None",
    "Sel. #cap Miscluster #cap Track Eff.",
    "Pure Track Eff.",
  };

  auto* pie = new TPie("pie_decomp",
    "TRKPTZ Failure - Exclusive Category Decomposition",
    NS, pieVals, pieCols, pieLabels);
  pie->SetValueFormat("%4.0f/2076");
  pie->SetLabelFormat("#splitline{%txt}{%val %perc}");
  pie->SetCircle(0.45, 0.35, 0.28);
  pie->SetTextSize(0.018);
  pie->SetLabelsOffset(.01);
  pie->SetAngularOffset(20.);
  pie->SetEntryRadiusOffset(9, 0.05);
  pie->SetEntryRadiusOffset(11, 0.1);
  c->SetTopMargin(0.08);
  pie->Draw();

  auto* leg2 = pie->MakeLegend(0.02, 0.7, 0.98, 0.9);
  leg2->SetNColumns(4);
  leg2->SetTextSize(0.018);
  gPad->Update();

  c->Print(fname);
  c->Print(Form("%s]", fname));

  std::cout << "Saved plots to " << fname << "\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
auto main() -> int {
  // ── Data source ───────────────────────────────────────────────────────────
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = kFatal;

  std::vector<BinData> bins(NBINS);    // TRKPTZ path
  std::vector<BinData> bins_w(NBINS);  // WAVeS path (same oracle categories)

  // ── filterClusters: keep clusters with ≥ MIN_CLUSTER_TRACKS tracks ────────
  auto filterClusters = [](const std::vector<Cluster>& col) {
    std::vector<Cluster> out;
    for (const auto& c : col)
      if (c.nConstituents >= MIN_CLUSTER_TRACKS) out.push_back(c);
    return out;
  };

  // ── Event loop ────────────────────────────────────────────────────────────
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
    // Use ITERATIVE (production method) so WAVeS scores are populated.
    auto hgtdRaw = clusterTracksInTime(
      tracks, &branch, DIST_CUT_CONE,
      /*useSmearedTimes=*/false, /*checkTimeValid=*/true, /*smearRes=*/-1.0,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);
    const auto qualHGTD = filterClusters(hgtdRaw);

    // ── E. Per-event pass/fail for TRKPTZ and WAVeS ───────────────────────
    const Cluster* hgtdBest   = chooseTRKPTZ(qualHGTD);
    const Cluster* hgtdBestW  = chooseWAVeS(qualHGTD);
    const bool hgtd_trkptz_pass = hgtdBest  && hgtdBest->passEfficiency(&branch);
    const bool hgtd_waves_pass  = hgtdBestW && hgtdBestW->passEfficiency(&branch);

    bins[bin].total++;
    bins_w[bin].total++;
    if (hgtd_trkptz_pass) bins[bin].n_pass++;   else bins[bin].n_fail++;
    if (hgtd_waves_pass)  bins_w[bin].n_pass++;  else bins_w[bin].n_fail++;

    // Skip shared oracle work entirely if both passed
    if (hgtd_trkptz_pass && hgtd_waves_pass) continue;

    // ── F. Shared cluster-level properties ───────────────────────────────
    bool any_hgtd_pass = false;
    bool hgtd_hasPure  = false;
    for (const auto& c : qualHGTD) {
      if (c.passEfficiency(&branch)) any_hgtd_pass = true;
      if (c.purity > 0.75f)         hgtd_hasPure  = true;
    }

    // Category 1 (selection failure): depends on which cluster was chosen.
    // If the chosen cluster failed but another cluster in the collection passes,
    // the algorithm made a wrong selection choice.
    const bool cat_sel_T = !hgtd_trkptz_pass && any_hgtd_pass;
    const bool cat_sel_W = !hgtd_waves_pass  && any_hgtd_pass;

    // ── F2. Category 3: misclustering (score-independent) ─────────────────
    bool refined_fixes = false;
    {
      if (!qualHGTD.empty()) {
        auto bestIt = std::max_element(qualHGTD.begin(), qualHGTD.end(),
            [](const Cluster& a, const Cluster& b) {
              return a.scores.at(Score::TRKPTZ.id) < b.scores.at(Score::TRKPTZ.id);
            });
        Cluster refined = refineClusterTiming(*bestIt, &branch, DIST_CUT_REFINE);
        refined_fixes = refined.passEfficiency(&branch);
      }
    }
    const bool cat_miscl = !hgtd_hasPure || refined_fixes;

    // ── G. Category 2: timing misassignment (oracle — score-independent) ──
    auto idealRaw = clusterTracksInTime(
      tracks, &branch, DIST_CUT_CONE,
      /*useSmearedTimes=*/true, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/false);
    const auto qualIdeal = filterClusters(idealRaw);
    const Cluster* idealBest = chooseTRKPTZ(qualIdeal);
    const bool ideal_pass = idealBest && idealBest->passEfficiency(&branch);
    const bool cat_timing = ideal_pass;

    // ── H. Category 4: track efficiency loss (oracle — score-independent) ─
    bool cat_teff = false;
    if (!ideal_pass) {
      auto idealEffRaw = clusterTracksInTime(
        tracks, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/true, /*checkTimeValid=*/false, IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
        /*sortTracks=*/false, /*calcPurityFlag=*/false);
      const auto qualIdealEff = filterClusters(idealEffRaw);
      const Cluster* idealEffBest = chooseTRKPTZ(qualIdealEff);
      cat_teff = idealEffBest && idealEffBest->passEfficiency(&branch);
    }

    // ── I. Accumulate TRKPTZ path ─────────────────────────────────────────
    if (!hgtd_trkptz_pass) {
      const bool cat_none_T = !cat_sel_T && !cat_timing && !cat_miscl && !cat_teff;
      if (cat_none_T) {
        static constexpr const char* EVTDISPLAY_FMT =
          "python3 event_display.py --file_num %s --event_num %lld --extra_time %.2f";
        std::string fname   = chain.GetCurrentFile()->GetName();
        std::string fileNum = "?";
        auto p1 = fname.rfind("._");
        auto p2 = fname.rfind(".SuperNtuple");
        if (p1 != std::string::npos && p2 != std::string::npos)
          fileNum = fname.substr(p1 + 2, p2 - p1 - 2);
        const Long64_t localEntry = chain.GetTree()->GetReadEntry();
        std::cout << "\n===== cat_none event (TRKPTZ) =====\n";
        printf((std::string(EVTDISPLAY_FMT) + "\n").c_str(),
               fileNum.c_str(), localEntry, 0.0);
        for (const auto& c : hgtdRaw) {
          std::cout << "---------\n";
          std::cout << "t: " << c.values.at(0) << "\n";
          std::cout << "score: " << c.scores.at(Score::TRKPTZ.id) << "\n";
          std::cout << "purity: " << c.purity << " (" << c.nConstituents << " tracks)\n";
          std::cout << "passes? " << c.passEfficiency(&branch) << "\n";
          std::cout << "---------\n";
        }
        std::cout << "===== end cat_none event =====\n";
      }
      auto& B = bins[bin];
      if (cat_sel_T)  B.cat_sel++;
      if (cat_timing) B.cat_timing++;
      if (cat_miscl)  B.cat_miscl++;
      if (cat_teff)   B.cat_teff++;
      if (cat_none_T) B.cat_none++;
      if (cat_sel_T  && cat_timing) B.ov_s_t++;
      if (cat_sel_T  && cat_miscl)  B.ov_s_m++;
      if (cat_sel_T  && cat_teff)   B.ov_s_e++;
      if (cat_timing && cat_miscl)  B.ov_t_m++;
      if (cat_miscl  && cat_teff)   B.ov_m_e++;
      if (cat_sel_T && cat_timing && cat_miscl) B.ov_stm++;
      if (cat_sel_T && cat_miscl  && cat_teff)  B.ov_sme++;
    }

    // ── J. Accumulate WAVeS path ──────────────────────────────────────────
    if (!hgtd_waves_pass) {
      const bool cat_none_W = !cat_sel_W && !cat_timing && !cat_miscl && !cat_teff;
      auto& W = bins_w[bin];
      if (cat_sel_W)  W.cat_sel++;
      if (cat_timing) W.cat_timing++;
      if (cat_miscl)  W.cat_miscl++;
      if (cat_teff)   W.cat_teff++;
      if (cat_none_W) W.cat_none++;
      if (cat_sel_W  && cat_timing) W.ov_s_t++;
      if (cat_sel_W  && cat_miscl)  W.ov_s_m++;
      if (cat_sel_W  && cat_teff)   W.ov_s_e++;
      if (cat_timing && cat_miscl)  W.ov_t_m++;
      if (cat_miscl  && cat_teff)   W.ov_m_e++;
      if (cat_sel_W && cat_timing && cat_miscl) W.ov_stm++;
      if (cat_sel_W && cat_miscl  && cat_teff)  W.ov_sme++;
    }
  }

  std::cout << "\nFINISHED PROCESSING\n";

  // ── Print results ─────────────────────────────────────────────────────────
  const BinData g  = aggregate(bins);
  const BinData gw = aggregate(bins_w);
  printSummary(g,  "TRKPTZ");
  printSummary(gw, "WAVeS");
  printPerBinTable(bins);

  // ── Generate plots ────────────────────────────────────────────────────────
  plotFailures(bins);

  return 0;
}
