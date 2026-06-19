// ---------------------------------------------------------------------------
// failure_decomposition.cxx
//
// For every event where the selection score fails the efficiency test, assigns
// the failure to a SINGLE principal cause via a fixed precedence ladder.  An
// event can satisfy several diagnostic conditions at once (e.g. timing
// misassignment AND the hard-scatter pT split across clusters); rather than
// reporting every overlap, we attribute the event to the most prominent /
// most actionable cause — the first rung of this ladder that "fixes" it:
//
//   1. Timing misassignment — the HS timing is the culprit, by either route:
//                             (a) the Ideal Resolution scenario passes
//                             (resolution-limited), or (b) a majority of the HS
//                             pT is mis-timed vs the truth vertex time
//                             (|pull| >= 3sigma), which ideal resolution can't
//                             fix because the times are offset, not just smeared.
//
//   2. Selection failure    — another qualifying cluster in the same real-time
//                             collection passes passEfficiency(); the right
//                             answer was present but not chosen.
//
//   3. Misclustering        — grouping ALL truth-HS tracks into one cluster
//                             (real times, no PU contamination) passes the
//                             window, but the algorithm's cluster does not: the
//                             HS energy was split / contaminated.
//
//   4. Track efficiency loss — Ideal Res. + Ideal Eff. passes but Ideal Res.
//                             alone does not: giving currently time-less HS
//                             tracks an ideal time is what recovers the event.
//
//   5. Irreducible          — none of the above; even a perfectly-grouped HS
//                             cluster with ideal timing+efficiency cannot reach
//                             the window.  Not an algorithm failure.
//
// Precedence (first match wins):  Timing > Selection > Misclustering > TrkEff.
//
// The PDF decomposition is produced for WAVeS (the primary analysis score);
// the TRKPTZ baseline is still printed as a console table.
//
// Output:
//   • Console summary table (WAVeS + TRKPTZ): principal-cause counts, the five
//     buckets are mutually exclusive and sum to n_fail.
//   • Per-bin breakdown vs n Forward HS Tracks (composition sums to 100%).
//   • PDF page 1: per-bin stacked bar of principal-cause composition (WAVeS).
//   • PDF page 2: overall pie of the five principal causes (WAVeS).
// ---------------------------------------------------------------------------

#include <TROOT.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include "clustering_constants.h"
#include "event_processing.h"

using namespace MyUtl;

// Track pT floor (GeV) for the truth-agnostic PU-suppression cross-check oracle.
static constexpr double PU_SUPPRESS_PT = 2.0;

// Fraction of HS pT that must be mis-timed (vs the truth vertex time) for an
// event to count as a timing-misassignment failure on that basis alone.
static constexpr double BAD_HS_TIME_FRAC = 0.5;

// Precedence cross-check: the 6 orderings of the top-3 failure modes
// (0 = timing, 1 = selection, 2 = misclustering).  TrkEff / Irreducible are
// rungs 4-5 and invariant to this reshuffle.
static const int   TOP3_ORD[6][3]   = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
static const char* TOP3_ORD_LBL[6]  = {"Timing > Sel > Miscl  (current)",
                                       "Timing > Miscl > Sel",
                                       "Sel > Timing > Miscl",
                                       "Sel > Miscl > Timing",
                                       "Miscl > Timing > Sel",
                                       "Miscl > Sel > Timing"};

// ---------------------------------------------------------------------------
// Per-bin accumulator — five mutually-exclusive principal-cause buckets
// ---------------------------------------------------------------------------
struct BinData {
  double total   = 0;   // events passing event selection
  double n_pass  = 0;   // selection score passes efficiency
  double n_fail  = 0;   // selection score fails efficiency (our population)

  // Principal cause (exactly one per failing event; sum == n_fail)
  double pc_timing = 0; // 1. timing misassignment
  double pc_sel    = 0; // 2. selection failure
  double pc_miscl  = 0; // 3. misclustering
  double pc_teff   = 0; // 4. track efficiency loss
  double pc_none   = 0; // 5. fundamental

  // Overlap-aware composition: comb[m] = #events whose set of applicable top-3
  // modes is the bitmask m (bit0=timing, bit1=selection, bit2=misclustering).
  // comb[0] (no top-3 mode) splits into pc_teff + pc_none.  Used by the plots.
  double comb[8] = {};
};

static const int NBINS = static_cast<int>(FOLD_HS_TRACK) + 1;

// ---------------------------------------------------------------------------
// chooseBest — highest-scoring qualifying cluster for a given score, or nullptr
// ---------------------------------------------------------------------------
static const Cluster* chooseBest(const std::vector<Cluster>& qual, int scoreId) {
  if (qual.empty()) return nullptr;
  const Cluster* best = &qual[0];
  double maxScore = best->scores.at(scoreId);
  for (const auto& c : qual) {
    double s = c.scores.at(scoreId);
    if (s > maxScore) { maxScore = s; best = &c; }
  }
  return best;
}

// ---------------------------------------------------------------------------
// hsCorePass — keeping a cluster's grouping, drop its truth-PU tracks and
//   recompute the time from the truth-HS core alone (grouped as one cluster,
//   real times).  Returns whether that decontaminated time passes the window;
//   the HS-core time is returned via tOut.  Empty HS core → false.
// ---------------------------------------------------------------------------
static bool hsCorePass(const Cluster& c, BranchPointerWrapper* branch, double& tOut) {
  std::vector<int> core;
  for (int idx : c.trackIndices)
    if (branch->trackToTruthvtx[idx] == 0 && branch->trackTimeValid[idx] == 1)
      core.push_back(idx);
  tOut = 0.;
  if (core.empty()) return false;
  auto cl = clusterTracksInTime(core, branch, /*distCut=*/1e9,
              /*useSmearedTimes=*/false, /*checkTimeValid=*/true, /*smearRes=*/-1.0,
              ClusteringMethod::ITERATIVE, /*useZ0=*/false,
              /*sortTracks=*/false, /*calcPurityFlag=*/false);
  bool pass = false;
  for (const auto& cc : cl) { tOut = cc.values.at(0); if (cc.passEfficiency(branch)) pass = true; }
  return pass;
}

// ---------------------------------------------------------------------------
// assignPrincipal — attribute a failing event to one bucket via the ladder
//   Timing > Selection > Misclustering > TrkEff > Fundamental
// ---------------------------------------------------------------------------
static void assignPrincipal(BinData& B,
                            bool timing, bool sel, bool miscl, bool teff) {
  if      (timing) B.pc_timing++;
  else if (sel)    B.pc_sel++;
  else if (miscl)  B.pc_miscl++;
  else if (teff)   B.pc_teff++;
  else             B.pc_none++;
}

// ---------------------------------------------------------------------------
// aggregate — sum bins into a single BinData total
// ---------------------------------------------------------------------------
static BinData aggregate(const std::vector<BinData>& bins) {
  BinData t;
  for (const auto& b : bins) {
    t.total     += b.total;
    t.n_pass    += b.n_pass;
    t.n_fail    += b.n_fail;
    t.pc_timing += b.pc_timing;
    t.pc_sel    += b.pc_sel;
    t.pc_miscl  += b.pc_miscl;
    t.pc_teff   += b.pc_teff;
    t.pc_none   += b.pc_none;
    for (int m = 0; m < 8; ++m) t.comb[m] += b.comb[m];
  }
  return t;
}

// ---------------------------------------------------------------------------
// printSummary — principal-cause table for a BinData aggregate
// ---------------------------------------------------------------------------
static void printSummary(const BinData& g, const char* label) {
  const double nf = g.n_fail;
  auto pf = [&](double n) { return nf > 0 ? 100.*n/nf : 0.; };       // % of fail
  auto pt = [&](double n) { return g.total > 0 ? 100.*n/g.total : 0.; }; // % of total

  std::cout << "\n";
  std::cout << "================================================================\n";
  printf("  %s Failure Decomposition (principal cause per event)\n", label);
  std::cout << "================================================================\n";
  printf("  Total events (passing selection): %9.0f\n", g.total);
  printf("  %s pass:  %9.0f  (%5.2f%%)\n", label, g.n_pass, pt(g.n_pass));
  printf("  %s fail:  %9.0f  (%5.2f%%)\n", label, g.n_fail, pt(g.n_fail));

  std::cout << "\n  --- Principal cause among failing events (N = "
            << (int)g.n_fail << ") ---\n\n";
  printf("  %-42s %8s %9s %9s\n", "Principal cause", "Count", "%fail", "%total");
  std::cout << "  " << std::string(70, '-') << "\n";
  printf("  %-42s %8.0f %8.1f%% %8.1f%%\n",
         "1. Timing misassignment",  g.pc_timing, pf(g.pc_timing), pt(g.pc_timing));
  printf("  %-42s %8.0f %8.1f%% %8.1f%%\n",
         "2. Selection failure",     g.pc_sel,    pf(g.pc_sel),    pt(g.pc_sel));
  printf("  %-42s %8.0f %8.1f%% %8.1f%%\n",
         "3. Misclustering",         g.pc_miscl,  pf(g.pc_miscl),  pt(g.pc_miscl));
  printf("  %-42s %8.0f %8.1f%% %8.1f%%\n",
         "4. Track efficiency loss", g.pc_teff,   pf(g.pc_teff),   pt(g.pc_teff));
  printf("  %-42s %8.0f %8.1f%% %8.1f%%\n",
         "5. Irreducible (no oracle recovers)", g.pc_none, pf(g.pc_none), pt(g.pc_none));
  std::cout << "  " << std::string(70, '-') << "\n";

  double sum = g.pc_timing + g.pc_sel + g.pc_miscl + g.pc_teff + g.pc_none;
  printf("  Sanity check (should equal n_fail = %.0f):  %.0f%s\n",
         g.n_fail, sum, std::abs(sum - g.n_fail) < 0.5 ? "  OK" : "  *** MISMATCH ***");
  std::cout << "================================================================\n";
}

// ---------------------------------------------------------------------------
// printPerBinTable — per-bin principal-cause composition (sums to 100%)
// ---------------------------------------------------------------------------
static void printPerBinTable(const std::vector<BinData>& bins, const char* label) {
  const BinData g = aggregate(bins);

  std::cout << "\n---" << std::string(98, '-') << "\n";
  printf("  %s per-bin principal-cause composition vs n Forward HS Tracks\n", label);
  std::cout << "  Cause percentages are % of failing events in that bin (sum to 100%).\n";
  std::cout << "---" << std::string(98, '-') << "\n";
  printf("  %-14s %7s %7s %7s | %8s %8s %8s %8s %8s\n",
         "Bin", "Total", "Fail", "Fail%",
         "Timing%", "Sel%", "Miscl%", "TrkEff%", "Irred%");
  std::cout << "---" << std::string(98, '-') << "\n";

  auto rowFmt = [](const char* lbl, const BinData& b) {
    double fp  = b.total  > 0 ? 100.*b.n_fail/b.total     : 0.;
    double tip = b.n_fail > 0 ? 100.*b.pc_timing/b.n_fail : 0.;
    double sep = b.n_fail > 0 ? 100.*b.pc_sel/b.n_fail    : 0.;
    double mcp = b.n_fail > 0 ? 100.*b.pc_miscl/b.n_fail  : 0.;
    double tep = b.n_fail > 0 ? 100.*b.pc_teff/b.n_fail   : 0.;
    double nop = b.n_fail > 0 ? 100.*b.pc_none/b.n_fail   : 0.;
    printf("  %-14s %7.0f %7.0f %6.1f%% | %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n",
           lbl, b.total, b.n_fail, fp, tip, sep, mcp, tep, nop);
  };

  for (int b = 0; b < NBINS; ++b) {
    if (bins[b].total == 0) continue;
    char lbl[32];
    if (b < (int)FOLD_HS_TRACK) std::snprintf(lbl, sizeof(lbl), "[%d, %d)", b, b+1);
    else                        std::snprintf(lbl, sizeof(lbl), "[%d, inf)", b);
    rowFmt(lbl, bins[b]);
  }
  std::cout << "---" << std::string(98, '-') << "\n";
  rowFmt("Total", g);
  std::cout << "---" << std::string(98, '-') << "\n";
}

// ---------------------------------------------------------------------------
// drawStriped — fill [x0,x1]x[y0,y1] in the current pad's user coords: one
//   solid colour for a single-mode event class, or a woven multi-colour hatch
//   for a class fixable equally by several modes (one colour per stripe).
// ---------------------------------------------------------------------------
static void drawStriped(double x0, double y0, double x1, double y1,
                        const std::vector<Int_t>& cols) {
  auto* base = new TBox(x0, y0, x1, y1);
  base->SetFillColor(cols[0]); base->SetFillStyle(1001); base->Draw();
  const int pat[2] = {3144, 3145};               // opposite diagonals → weave
  for (size_t i = 1; i < cols.size() && i <= 2; ++i) {
    auto* h = new TBox(x0, y0, x1, y1);
    h->SetFillColor(cols[i]); h->SetFillStyle(pat[i-1]); h->Draw();
  }
  auto* bd = new TBox(x0, y0, x1, y1);           // black outline
  bd->SetFillStyle(0); bd->SetLineColor(kBlack); bd->SetLineWidth(2); bd->Draw();
}

// ---------------------------------------------------------------------------
// plotFailures — ../figs/failure_decomposition.pdf (2 pages)
//   Overlap-aware composition: events are grouped by the SET of top-3 modes
//   that can fix them.  Single-mode events are solid (timing = red,
//   selection = tan, misclustering = blue); events equally fixable by several
//   modes are woven stripes of the applicable colours.  TrkEff (violet) and
//   Irreducible (gray) are solid.
//   Page 1: per-bin composition (sums to 1 per bin) over an n_fail dist pad.
//   Page 2: overall composition column + full legend with counts.
// ---------------------------------------------------------------------------
static void plotFailures(const std::vector<BinData>& bins, const char* label) {
  gStyle->SetOptStat(0);
  gStyle->SetHatchesLineWidth(2);   // thicker, more legible hatch lines
  gStyle->SetHatchesSpacing(3.0);   // slightly wider spacing so they don't crowd

  const double xData = static_cast<double>(FOLD_HS_TRACK) + 0.5;  // 20.5
  const double xMin  = -0.5;
  const double xMax  = xData + 13.5;            // right margin reserved for legend
  const char*  xLbl  = "n Forward HS Tracks (with valid HGTD time)";
  const char*  fname = "../figs/failure_decomposition.pdf";

  const Int_t COL_TIMING = COLORS[1]; // red
  const Int_t COL_MISCL  = COLORS[0]; // blue
  const Int_t COL_SEL    = COLORS[2]; // tan
  const Int_t COL_TEFF   = COLORS[4]; // violet
  const Int_t COL_NONE   = COLORS[3]; // gray

  // Stack segments (bottom→top). kind 0..7 = top-3 bitmask (bit0=timing,
  // bit1=selection, bit2=misclustering); kind 8 = TrkEff; kind 9 = Irreducible.
  struct Seg { int kind; std::vector<Int_t> cols; const char* name; };
  const std::vector<Seg> segs = {
    {1, {COL_TIMING},                     "Timing"},
    {5, {COL_TIMING, COL_MISCL},          "Timing #cap Miscl."},
    {3, {COL_TIMING, COL_SEL},            "Timing #cap Selection"},
    {7, {COL_TIMING, COL_SEL, COL_MISCL}, "All three"},
    {2, {COL_SEL},                        "Selection"},
    {6, {COL_SEL, COL_MISCL},             "Selection #cap Miscl."},
    {4, {COL_MISCL},                      "Misclustering"},
    {8, {COL_TEFF},                       "Track eff. loss"},
    {9, {COL_NONE},                       "Irreducible"},
  };
  auto segCount = [](const BinData& B, const Seg& s) -> double {
    if (s.kind < 8) return B.comb[s.kind];
    return (s.kind == 8) ? B.pc_teff : B.pc_none;
  };

  auto* h_fail = new TH1D("pc_faildist", "", NBINS, xMin, xData);
  for (int b = 0; b < NBINS; ++b) h_fail->SetBinContent(b + 1, bins[b].n_fail);

  TCanvas* c = new TCanvas("c_fail_decomp", "Failure decomposition", 950, 700);
  const double split = 0.25;
  const double lm = 0.10, rm = 0.02;

  auto* padTop = new TPad("padTop", "", 0, split, 1, 1);
  padTop->SetTopMargin(0.10); padTop->SetBottomMargin(0.00);
  padTop->SetLeftMargin(lm);  padTop->SetRightMargin(rm);
  padTop->Draw();
  auto* padBot = new TPad("padBot", "", 0, 0, 1, split);
  padBot->SetTopMargin(0.00); padBot->SetBottomMargin(0.38);
  padBot->SetLeftMargin(lm);  padBot->SetRightMargin(rm);
  padBot->Draw();

  // ── Top pad: per-bin striped composition ─────────────────────────────────
  padTop->cd();
  auto* fr = static_cast<TH1*>(gPad->DrawFrame(xMin, 0, xMax, 1.0));
  fr->SetTitle(Form("%s Failure - Composition by Recoverability;;Fraction of failing events", label));
  fr->GetXaxis()->SetLabelSize(0);
  fr->GetXaxis()->SetTickLength(0);
  fr->GetYaxis()->SetTitleOffset(0.7);
  fr->GetYaxis()->SetNdivisions(505);

  for (int b = 0; b < NBINS; ++b) {
    const double nf = bins[b].n_fail;
    if (nf <= 0) continue;
    double yc = 0.;
    const double x0 = b - 0.45, x1 = b + 0.45;
    for (const auto& s : segs) {
      double h = segCount(bins[b], s) / nf;
      if (h <= 0) continue;
      drawStriped(x0, yc, x1, yc + h, s.cols);
      yc += h;
    }
  }

  // legend keys in the reserved right margin (data ends at xData)
  {
    const double kx0 = xData + 1.0, kx1 = xData + 2.2, tx = xData + 2.6;
    auto key = [&](double y, const std::vector<Int_t>& cols, const char* txt) {
      drawStriped(kx0, y - 0.03, kx1, y + 0.03, cols);
      auto* t = new TLatex(tx, y, txt);
      t->SetTextSize(0.045); t->SetTextAlign(12); t->Draw();
    };
    auto* hdr = new TLatex(kx0, 0.97, "#bf{Recoverable by}");
    hdr->SetTextSize(0.045); hdr->SetTextAlign(12); hdr->Draw();
    key(0.88, {COL_TIMING}, "Timing");
    key(0.78, {COL_SEL},    "Selection");
    key(0.68, {COL_MISCL},  "Misclustering");
    key(0.58, {COL_TEFF},   "Track eff.");
    key(0.48, {COL_NONE},   "Irreducible");
    key(0.32, {COL_TIMING, COL_MISCL}, "#geq2 modes");
    auto* note = new TLatex(kx0, 0.22, "#splitline{(stripes = event fixable}{ by multiple modes)}");
    note->SetTextSize(0.038); note->SetTextAlign(12); note->Draw();
  }
  gPad->Update();

  // ── Bottom pad: n_fail distribution ──────────────────────────────────────
  padBot->cd();
  const double sf = 1.0 / split;
  double yBotMax = h_fail->GetMaximum() * 1.15;
  auto* hframe = static_cast<TH1*>(gPad->DrawFrame(xMin, 0, xMax, yBotMax));
  hframe->GetXaxis()->SetTitle(xLbl);
  hframe->GetXaxis()->SetNdivisions(510);
  hframe->GetXaxis()->SetLabelSize(0.035 * sf);
  hframe->GetXaxis()->SetTitleSize(0.035 * sf);
  hframe->GetXaxis()->SetTitleOffset(0.85);
  hframe->GetYaxis()->SetTitle("Failing events");
  hframe->GetYaxis()->SetLabelSize(0.030 * sf);
  hframe->GetYaxis()->SetTitleSize(0.030 * sf);
  hframe->GetYaxis()->SetTitleOffset(0.35);
  hframe->GetYaxis()->SetNdivisions(503);
  h_fail->SetFillColor(COLORS[8]);
  h_fail->SetFillStyle(3004);
  h_fail->SetLineWidth(2);
  h_fail->Draw("HIST SAME");
  gPad->Update();

  c->cd();
  c->Print(Form("%s[", fname));
  c->Print(fname);
  c->Clear();

  // ── Page 2: overall composition column + detailed legend ─────────────────
  const BinData g = aggregate(bins);
  c->cd();
  auto* p2 = new TPad("p2", "", 0, 0, 1, 1);
  p2->SetTopMargin(0.10); p2->SetLeftMargin(0.05); p2->SetRightMargin(0.03);
  p2->Draw(); p2->cd();
  auto* fr2 = static_cast<TH1*>(gPad->DrawFrame(0, 0, 1, 1));
  fr2->SetTitle(Form("%s Failure - Overall Composition (N = %.0f failing events)", label, g.n_fail));
  fr2->GetXaxis()->SetTickLength(0); fr2->GetXaxis()->SetLabelSize(0);
  fr2->GetYaxis()->SetTickLength(0); fr2->GetYaxis()->SetLabelSize(0);

  // one tall stacked column on the left
  double yc = 0.;
  const double bx0 = 0.12, bx1 = 0.34;
  for (const auto& s : segs) {
    double h = g.n_fail > 0 ? segCount(g, s) / g.n_fail : 0.;
    if (h <= 0) continue;
    drawStriped(bx0, yc, bx1, yc + h, s.cols);
    yc += h;
  }
  { auto* t0 = new TLatex(bx0 - 0.01, 0.0, "0%");
    t0->SetTextSize(0.028); t0->SetTextAlign(32); t0->Draw();
    auto* t1 = new TLatex(bx0 - 0.01, 1.0, "100%");
    t1->SetTextSize(0.028); t1->SetTextAlign(32); t1->Draw(); }

  // detailed legend on the right: every present class with count + %
  {
    double y = 0.90; const double dy = 0.085;
    auto* hdr = new TLatex(0.46, 0.985, "#bf{Recoverable by  (count, % of failures):}");
    hdr->SetTextSize(0.030); hdr->SetTextAlign(13); hdr->Draw();
    for (const auto& s : segs) {
      double cnt = segCount(g, s);
      if (cnt <= 0) continue;
      drawStriped(0.46, y - 0.025, 0.53, y + 0.025, s.cols);
      auto* t = new TLatex(0.55, y, Form("%s:  %.0f  (%.1f%%)",
                  s.name, cnt, g.n_fail > 0 ? 100.*cnt/g.n_fail : 0.));
      t->SetTextSize(0.028); t->SetTextAlign(12); t->Draw();
      y -= dy;
    }
  }
  gPad->Update();

  c->Print(fname);
  c->Print(Form("%s]", fname));
  std::cout << "Saved plots to " << fname << "\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
auto main() -> int {
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = kFatal;

  std::vector<BinData> bins(NBINS);    // TRKPTZ path (console only)
  std::vector<BinData> bins_w(NBINS);  // WAVeS path (console + plots)

  // Misclustering-oracle cross-check tallies (among WAVeS failures)
  double xc_truth = 0, xc_ptcut = 0, xc_both = 0;
  double xc_badtime = 0;  // failures with > BAD_HS_TIME_FRAC of HS pT mis-timed
  double n_irred = 0, n_irred_zerohs = 0;  // irreducible WAVeS fails; subset w/ 0 HS tracks
  double ord_cnt[6][3] = {};  // WAVeS top-3 split per ordering: [ordering][timing/sel/miscl]

  auto filterClusters = [](const std::vector<Cluster>& col) {
    std::vector<Cluster> out;
    for (const auto& c : col)
      if (c.nConstituents >= MIN_CLUSTER_TRACKS) out.push_back(c);
    return out;
  };

  std::cout << "Starting Event Loop\n";
  const Long64_t N_EVENT = chain.GetEntries();

  while (reader.Next()) {
    const Long64_t READ_NUM = chain.GetReadEntry() + 1;
    if (READ_NUM % 100 == 0)
      std::cout << "Progress: " << READ_NUM << "/" << N_EVENT << "\r" << std::flush;

    // ── A. Event + track selection ────────────────────────────────────────
    if (!branch.passBasicCuts()) continue;
    if (!branch.passJetPtCut())  continue;
    std::vector<int> tracks =
      getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);

    // ── B. Truth-HS tracks (valid HGTD time) + bin index ──────────────────
    std::vector<int> hsTracks;
    for (int idx : tracks)
      if (branch.trackToTruthvtx[idx] == 0 && branch.trackTimeValid[idx] == 1)
        hsTracks.push_back(idx);
    const int nHSTrack = static_cast<int>(hsTracks.size());
    const int bin = folded(nHSTrack, static_cast<int>(FOLD_HS_TRACK));

    // ── C. HGTD clustering (real times, purity enabled) ───────────────────
    auto hgtdRaw = clusterTracksInTime(
      tracks, &branch, DIST_CUT_CONE,
      /*useSmearedTimes=*/false, /*checkTimeValid=*/true, /*smearRes=*/-1.0,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/true);
    const auto qualHGTD = filterClusters(hgtdRaw);

    // ── D. Per-event pass/fail for TRKPTZ and WAVeS ───────────────────────
    const Cluster* hgtdBestT = chooseBest(qualHGTD, Score::TRKPTZ.id);
    const Cluster* hgtdBestW = chooseBest(qualHGTD, Score::WAVES.id);
    const bool pass_T = hgtdBestT && hgtdBestT->passEfficiency(&branch);
    const bool pass_W = hgtdBestW && hgtdBestW->passEfficiency(&branch);

    bins[bin].total++;   bins_w[bin].total++;
    if (pass_T) bins[bin].n_pass++;   else bins[bin].n_fail++;
    if (pass_W) bins_w[bin].n_pass++; else bins_w[bin].n_fail++;

    if (pass_T && pass_W) continue;  // nothing to decompose

    // ── E. Selection failure (does another cluster already pass?) ─────────
    bool any_hgtd_pass = false;
    for (const auto& cc : qualHGTD)
      if (cc.passEfficiency(&branch)) { any_hgtd_pass = true; break; }
    const bool cat_sel_T = !pass_T && any_hgtd_pass;
    const bool cat_sel_W = !pass_W && any_hgtd_pass;

    // ── F. Misclustering oracles (score-independent) ──────────────────────
    // Headline (truth-PU decontamination, keeps the algorithm's grouping): for
    // each real cluster, recompute the time from its truth-HS tracks only.  If
    // any decontaminated cluster passes, a good HS core existed and PU
    // contamination is what dragged the measured time out of the window —
    // i.e. better PU rejection (not better clustering granularity) recovers it.
    bool cat_miscl = false;
    for (const auto& c : qualHGTD) {
      double tcore;
      if (hsCorePass(c, &branch, tcore)) { cat_miscl = true; break; }
    }
    // Cross-check (truth-agnostic): re-cluster only tracks above PU_SUPPRESS_PT
    // to suppress (predominantly low-pT) pile-up.  If a passing cluster emerges,
    // a realistic label-free pT cut would also have recovered the event.  Used
    // for comparison only — it does not drive the classification.
    bool miscl_ptcut = false;
    {
      std::vector<int> hiPt;
      for (int idx : tracks)
        if (branch.trackPt[idx] >= PU_SUPPRESS_PT) hiPt.push_back(idx);
      if (!hiPt.empty()) {
        auto reCl = clusterTracksInTime(
          hiPt, &branch, DIST_CUT_CONE,
          /*useSmearedTimes=*/false, /*checkTimeValid=*/true, /*smearRes=*/-1.0,
          ClusteringMethod::ITERATIVE, /*useZ0=*/false,
          /*sortTracks=*/false, /*calcPurityFlag=*/false);
        for (const auto& cc : reCl)
          if (cc.passEfficiency(&branch)) { miscl_ptcut = true; break; }
      }
    }

    // ── G. Timing misassignment ──────────────────────────────────────────
    // Two routes, either of which makes the HS *timing* the culprit:
    //  (a) Ideal-resolution oracle — re-smearing every track to its truth time
    //      with 1 ps resolution recovers the event (resolution-limited).
    //  (b) Bad-HS-timing gate — a majority of the HS pT has a measured time
    //      inconsistent with the truth VERTEX time (|pull| ≥ 3σ vs truthVtxTime,
    //      the reference passEfficiency grades against).  This catches events
    //      whose HS tracks are systematically offset from the vertex time, which
    //      ideal *resolution* cannot fix (it smears around the track truth time,
    //      not the vertex time).
    double hs_pt_tot = 0., hs_pt_bad = 0.;
    for (int idx : hsTracks) {
      double pull = std::abs(branch.trackTime[idx] - branch.truthVtxTime[0])
                    / branch.trackTimeRes[idx];
      double pT = branch.trackPt[idx];
      hs_pt_tot += pT;
      if (pull >= TRUTH_PULL_CUT) hs_pt_bad += pT;
    }
    const bool cat_badtime =
      hs_pt_tot > 0. && (hs_pt_bad / hs_pt_tot) > BAD_HS_TIME_FRAC;

    auto idealRaw = clusterTracksInTime(
      tracks, &branch, DIST_CUT_CONE,
      /*useSmearedTimes=*/true, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
      ClusteringMethod::ITERATIVE, /*useZ0=*/false,
      /*sortTracks=*/false, /*calcPurityFlag=*/false);
    const auto qualIdeal = filterClusters(idealRaw);
    const Cluster* idealBestT = chooseBest(qualIdeal, Score::TRKPTZ.id);
    const Cluster* idealBestW = chooseBest(qualIdeal, Score::WAVES.id);
    const bool cat_timing_T = cat_badtime || (idealBestT && idealBestT->passEfficiency(&branch));
    const bool cat_timing_W = cat_badtime || (idealBestW && idealBestW->passEfficiency(&branch));

    // ── H. Track efficiency loss oracle (Ideal Res. + Ideal Eff.) ─────────
    bool cat_teff_T = false, cat_teff_W = false;
    if (!cat_timing_T || !cat_timing_W) {
      auto idealEffRaw = clusterTracksInTime(
        tracks, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/true, /*checkTimeValid=*/false, IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
        /*sortTracks=*/false, /*calcPurityFlag=*/false);
      const auto qualIdealEff = filterClusters(idealEffRaw);
      const Cluster* effBestT = chooseBest(qualIdealEff, Score::TRKPTZ.id);
      const Cluster* effBestW = chooseBest(qualIdealEff, Score::WAVES.id);
      cat_teff_T = !cat_timing_T && effBestT && effBestT->passEfficiency(&branch);
      cat_teff_W = !cat_timing_W && effBestW && effBestW->passEfficiency(&branch);
    }

    // ── I. Assign principal cause (Timing > Selection > Miscl > TrkEff) ───
    if (!pass_T) {
      assignPrincipal(bins[bin],   cat_timing_T, cat_sel_T, cat_miscl, cat_teff_T);
      bins[bin].comb[(cat_timing_T?1:0)|(cat_sel_T?2:0)|(cat_miscl?4:0)]++;
    }

    if (!pass_W) {
      assignPrincipal(bins_w[bin], cat_timing_W, cat_sel_W, cat_miscl, cat_teff_W);
      bins_w[bin].comb[(cat_timing_W?1:0)|(cat_sel_W?2:0)|(cat_miscl?4:0)]++;

      // Cross-check: how often each misclustering oracle fires among WAVeS fails
      xc_truth += cat_miscl;
      xc_ptcut += miscl_ptcut;
      xc_both  += (cat_miscl && miscl_ptcut);
      xc_badtime += cat_badtime;

      // Precedence cross-check: attribute this event under all 6 top-3 orderings.
      // (If none of the top-3 fire, the event is TrkEff/Irreducible — invariant.)
      const bool top3[3] = { cat_timing_W, cat_sel_W, cat_miscl };
      for (int o = 0; o < 6; ++o)
        for (int r = 0; r < 3; ++r) {
          int cid = TOP3_ORD[o][r];
          if (top3[cid]) { ord_cnt[o][cid]++; break; }
        }

      // Irreducible WAVeS failures: not recovered by selection, ideal timing,
      // ideal efficiency, or PU decontamination.  Emit a ready-to-run event
      // display command plus per-cluster HS-core diagnostics.
      const bool waves_irreducible =
        !cat_timing_W && !cat_sel_W && !cat_miscl && !cat_teff_W;
      if (waves_irreducible) {
        n_irred++;
        if (nHSTrack == 0) n_irred_zerohs++;
        std::string fpath   = chain.GetCurrentFile()->GetName();
        std::string fileNum = "?";
        auto p1 = fpath.rfind("._");
        auto p2 = fpath.rfind(".SuperNtuple");
        if (p1 != std::string::npos && p2 != std::string::npos)
          fileNum = fpath.substr(p1 + 2, p2 - p1 - 2);
        const Long64_t localEntry = chain.GetTree()->GetReadEntry();
        std::cout << "\n===== irreducible WAVeS failure =====\n";
        printf("python3 event_display.py --file_num %s --event_num %lld --extra_time 0.00\n",
               fileNum.c_str(), localEntry);
        printf("truth HS time: %.2f ps | HS pT mis-timed frac: %.2f | pT-cut recovers? %d\n",
               branch.truthVtxTime[0],
               hs_pt_tot > 0. ? hs_pt_bad / hs_pt_tot : -1., (int)miscl_ptcut);
        for (const auto& cc : hgtdRaw) {
          double tcore; bool corePass = hsCorePass(cc, &branch, tcore);
          std::cout << "--------- t: " << cc.values.at(0)
                    << "  waves: " << cc.scores.at(Score::WAVES.id)
                    << "  purity: " << cc.purity << " (" << cc.nConstituents << " trk)"
                    << "  HScore t: " << tcore << " (pass? " << corePass << ")"
                    << "  passes? " << cc.passEfficiency(&branch) << "\n";
        }
        std::cout << "===== end irreducible WAVeS failure =====\n";
      }
    }
  }

  std::cout << "\nFINISHED PROCESSING\n";

  printSummary(aggregate(bins_w), "WAVeS");
  printSummary(aggregate(bins),   "TRKPTZ");

  // Misclustering-oracle cross-check: how much does truth dependence inflate it?
  {
    const double nf = aggregate(bins_w).n_fail;
    auto p = [&](double n) { return nf > 0 ? 100.*n/nf : 0.; };
    std::cout << "\n  --- WAVeS misclustering oracle cross-check (among "
              << (int)nf << " failures) ---\n";
    printf("    truth-PU decontamination recovers : %6.0f  (%.1f%%)   [headline]\n",
           xc_truth, p(xc_truth));
    printf("    pT > %.1f GeV re-clustering recovers: %6.0f  (%.1f%%)   [label-free]\n",
           PU_SUPPRESS_PT, xc_ptcut, p(xc_ptcut));
    printf("    both                              : %6.0f  (%.1f%%)\n",
           xc_both, p(xc_both));
    printf("  --- of which timing: > %.0f%% HS pT mis-timed vs truth vtx : %6.0f  (%.1f%%) ---\n",
           100.*BAD_HS_TIME_FRAC, xc_badtime, p(xc_badtime));
    printf("  --- irreducible: %.0f total, %.0f with zero valid-time HS tracks (no HS timing) ---\n",
           n_irred, n_irred_zerohs);
  }

  // Precedence cross-check: how the top-3 split shifts under all 6 orderings.
  {
    const BinData g = aggregate(bins_w);
    const double top3tot = g.pc_timing + g.pc_sel + g.pc_miscl;  // invariant
    std::cout << "\n  --- WAVeS precedence cross-check: top-3 split under all 6 orderings ---\n";
    printf("    (TrkEff = %.0f and Irreducible = %.0f are rungs 4-5, invariant; "
           "top-3 always sums to %.0f)\n\n", g.pc_teff, g.pc_none, top3tot);
    printf("    %-34s %9s %10s %14s\n", "Ordering", "Timing", "Selection", "Misclustering");
    std::cout << "    " << std::string(70, '-') << "\n";
    auto pc = [&](double n) { return top3tot > 0 ? 100.*n/top3tot : 0.; };
    for (int o = 0; o < 6; ++o)
      printf("    %-34s %5.0f(%4.1f%%) %5.0f(%4.1f%%) %6.0f(%4.1f%%)\n",
             TOP3_ORD_LBL[o],
             ord_cnt[o][0], pc(ord_cnt[o][0]),
             ord_cnt[o][1], pc(ord_cnt[o][1]),
             ord_cnt[o][2], pc(ord_cnt[o][2]));
    std::cout << "    " << std::string(70, '-') << "\n";
  }

  printPerBinTable(bins_w, "WAVeS");

  plotFailures(bins_w, "WAVeS");
  return 0;
}
