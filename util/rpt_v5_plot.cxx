// rpt_v5_plot.cxx — Plotting stage of the rpt_v5 split.
//
// Reads the Hard-Scatter/Pile-Up RpT histograms and scalar accumulators
// written by rpt_v5_hist.cxx (no ntuple access, no event loop) and
// reproduces the same ROC/console-summary/PDF output the former monolithic
// rpt_v5.cxx produced -- see CLAUDE.md's "Main Executables" section.
//
// Structured in four lexical sections, in order:
//   1. Load histograms (raw, unstyled) from the saved ROOT file.
//   2. Derive every Integral()/GetMean()-based number (ROC graphs,
//      the harmonized [0.6, 1.0] ROC window, the console summary, the
//      rejection table).
//      styleScen (cosmetic axis restriction) is NOT YET DEFINED here --
//      see the sentinel comment below marking the end of this section.
//   3. Define + call styleScen. Because styleScen is a local lambda whose
//      *definition* now lives lexically after section 2, calling it early
//      would be a compile error (undeclared identifier), not just a
//      discipline violation -- see the bug this guards against, noted at
//      generate_roc()/rejAtEff() below.
//   4. PDF drawing (ROC+ratio, per-scenario pages, overlays).
//
// Output: <OUTPUT_DIR>/rpt_plots/rpt_v5.pdf

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#include <boost/filesystem.hpp>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "clustering_constants.h"
#include "sample_config.h"
#include "clustering_includes.h"
#include "rpt_v5_common.h"
#include "histogram_io.h"

using namespace MyUtl;

// -----------------------------------------------------------------------------
// ROC: HS efficiency vs PU rejection (1 / mistag).
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// generate_roc_grid — ROC sampled at fixed HS-efficiency working points.
//
// Same construction rpt_v6 uses to match the TDR figure. generate_roc below
// emits one point per histogram bin, which samples wherever the R_pT
// distribution happens to put thresholds: points bunch unevenly, and wherever
// the distribution has an atom (a large spike of jets sharing one R_pT value,
// e.g. jets with no associated tracks at all) no threshold exists inside it, so
// the curve shows a hard break.
//
// Sampling on efficiency instead walks the cumulative distributions and
// linearly interpolates 1/mistag at each target efficiency. Interpolating
// across an atom is the physically right thing -- it is what breaking ties at
// random inside that bin would deliver -- and because every scenario then
// shares identical x-points, ratio panels become exact rather than relying on
// TGraph::Eval.
//
// Same explicit-bin-range discipline as generate_roc: Integral(i, bin+1)
// throughout, never the no-arg form, so a prior SetRangeUser cannot silently
// truncate the denominators.
// -----------------------------------------------------------------------------
TGraph* generate_roc_grid(TH1D* PU_hist, TH1D* HS_hist,
                          const std::vector<double>& targets) {
  const int bin = PU_hist->GetNbinsX();
  const double HS_tot = HS_hist->Integral(1, bin + 1);
  const double PU_tot = PU_hist->Integral(1, bin + 1);
  std::vector<double> vx, vy;
  if (HS_tot <= 0 || PU_tot <= 0) return new TGraph();

  std::vector<double> S(bin + 2), B(bin + 2);
  for (int i = 1; i <= bin + 1; ++i) {
    S[i] = HS_hist->Integral(i, bin + 1) / HS_tot;
    B[i] = PU_hist->Integral(i, bin + 1) / PU_tot;
  }
  for (double e : targets) {
    if (e <= 0 || e > S[1]) continue;              // efficiency unreachable
    int i = 1;
    while (i <= bin + 1 && S[i] > e) ++i;
    if (i <= 1 || i > bin + 1) continue;
    const double s_hi = S[i - 1], s_lo = S[i];
    const double b_hi = B[i - 1], b_lo = B[i];
    const double f = (s_hi > s_lo) ? (s_hi - e) / (s_hi - s_lo) : 0.0;
    const double b = b_hi + f * (b_lo - b_hi);
    if (b <= 1e-9) continue;                       // background exhausted
    vx.push_back(e);
    vy.push_back(1.0 / b);
  }
  return new TGraph((int)vx.size(), vx.data(), vy.data());
}

// Working points matching the TDR figure's marker spacing.
inline std::vector<double> rocEffGrid() {
  std::vector<double> g;
  for (double e = 0.800; e <= 0.9751; e += 0.025) g.push_back(e);
  return g;
}

TGraph* generate_roc(TH1D* PU_hist, TH1D* HS_hist) {
  int bin = PU_hist->GetNbinsX();
  // Explicit bin range: Integral() with no arguments silently respects any
  // prior GetXaxis()->SetRangeUser() restriction (e.g. the display zoom
  // applied elsewhere in this file), while the per-bin Integral(i, bin+1)
  // calls below always use explicit bins and ignore it — so a no-arg call
  // here would silently truncate the denominator and inflate every
  // efficiency/rejection number. Match the per-bin calls' bin range instead.
  double HS_tot = HS_hist->Integral(1, bin + 1);
  double PU_tot = PU_hist->Integral(1, bin + 1);
  std::vector<double> vx, vy;
  for (int i = 1; i <= bin; ++i) {
    double HS_eff    = HS_hist->Integral(i, bin + 1) / HS_tot;
    double PU_mistag = PU_hist->Integral(i, bin + 1) / PU_tot;
    if (std::abs(PU_mistag) > 1e-6 && HS_eff < 0.99) {
      vx.push_back(HS_eff);
      vy.push_back(1.0 / PU_mistag);
    }
  }
  return new TGraph((int)vx.size(), vx.data(), vy.data());
}

int main(int argc, char** argv) {
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  // --- Sample selection (--sample=vbf|zjets|dijet; default: local VBF ntuple) ---
  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);  // --vbs-deta=<x>; sets SELECTION_TAG
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;

  // ===========================================================================
  // SECTION 1: Load histograms (raw, unstyled) from the saved ROOT file.
  // ===========================================================================
  const std::string histFile = MyUtl::resolveHistFile(
    argc, argv, MyUtl::histFilePath("rpt_v5_hist.root"));
  std::cout << "Reading histograms from " << histFile << '\n';

  std::vector<Scenario> scen_lo = makeScenarios("_lo");
  std::vector<Scenario> scen_hi = makeScenarios("_hi");
  // Central (|eta| < 2.4) counterparts of the same two pT slices -- the
  // ITk-only baseline the forward numbers are quoted against.
  std::vector<Scenario> scen_lo_cen = makeScenarios("_lo_cen");
  std::vector<Scenario> scen_hi_cen = makeScenarios("_hi_cen");
  // VBS-topology regions, see rpt_v5_hist.cxx's region block for definitions.
  // scen_r2's h_hs is EMPTY by construction (region 2 has no forward HS jet);
  // its ROC borrows scen_r1's HS side -- see SECTION 2b.
  std::vector<Scenario> scen_r1 = makeScenarios("_r1");
  std::vector<Scenario> scen_r2 = makeScenarios("_r2");
  long n_total, n_pass_basic, n_hgtd_valid;
  double pu_tot_pt, pu_floor_pt, hs_tot_pt, hs_floor_pt;
  double pu_tot_lo, pu_floor_lo, hs_tot_lo, hs_floor_lo;
  {
    MyUtl::HistReader reader(histFile);
    reader.CheckSelection(MyUtl::VBS_JET_D_ETA, MyUtl::VBS_JET_MJJ);
    loadScenarios(reader, scen_lo);
    loadScenarios(reader, scen_hi);
    loadScenarios(reader, scen_r1);
    loadScenarios(reader, scen_r2);
    loadScenarios(reader, scen_lo_cen);
    loadScenarios(reader, scen_hi_cen);
    n_total      = reader.ReadScalarLong("meta_n_total");
    n_pass_basic = reader.ReadScalarLong("meta_n_pass_basic");
    n_hgtd_valid = reader.ReadScalarLong("meta_n_hgtd_valid");
    pu_tot_pt    = reader.ReadScalarDouble("meta_pu_tot_pt");
    pu_floor_pt  = reader.ReadScalarDouble("meta_pu_floor_pt");
    hs_tot_pt    = reader.ReadScalarDouble("meta_hs_tot_pt");
    hs_floor_pt  = reader.ReadScalarDouble("meta_hs_floor_pt");
    pu_tot_lo    = reader.ReadScalarDouble("meta_pu_tot_lo");
    pu_floor_lo  = reader.ReadScalarDouble("meta_pu_floor_lo");
    hs_tot_lo    = reader.ReadScalarDouble("meta_hs_tot_lo");
    hs_floor_lo  = reader.ReadScalarDouble("meta_hs_floor_lo");
  }

  // ── Output paths ─────────────────────────────────────────────────────────────
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/rpt_plots");
  const TString out_pdf = MyUtl::plotFilePath("rpt_plots", "rpt_v5.pdf").c_str();

  TCanvas* canvas = new TCanvas("canvas", "RpT v5", 800, 700);
  canvas->Print(out_pdf + "[");

  auto drawLabels = [](const char* extra = nullptr) {
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
    if (extra) {
      TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
      t.DrawLatex(0.18, 0.76, extra);
    }
  };
  const char* lbl_lo = "30 < p_{T}^{jet} < 40 GeV@@2.4 < |#eta| < 3.8";
  const char* lbl_hi = "p_{T}^{jet} > 40 GeV@@2.4 < |#eta| < 3.8";
  // VBS-topology region labels (see rpt_v5_hist.cxx for the definitions).
  const char* lbl_lo_cen = "30 < p_{T}^{jet} < 40 GeV@@|#eta| < 2.4@@[central baseline]";
  const char* lbl_hi_cen = "p_{T}^{jet} > 40 GeV@@|#eta| < 2.4@@[central baseline]";
  const char* lbl_r1 = "VBS R1: both tags fwd@@HS leg vs PU leg@@p_{T}^{jet} > 30 GeV";
  const char* lbl_r2 = "VBS R2: fwd PU + cen. HS@@signal side from R1";

  // ===========================================================================
  // SECTION 2: Derive every Integral()/GetMean()-based number. styleScen is
  // NOT defined yet -- see the sentinel comment at the end of this section.
  // ===========================================================================

  // Harmonized "high-efficiency working point" window: fixed for every sample
  // and pT-slice (rather than a per-sample adaptive range) so plots are
  // directly comparable across samples. [0.80, 1.0] matches the x-range the
  // reference study's ROCs use, so our panels can be read side by side with
  // theirs; it also starts exactly where rocEffGrid()'s first sampled working
  // point sits, so no computed point falls outside the drawn window.
  // (An earlier [0.85, 1.0] zoom cropped the first grid point; the still
  // wider [0.6, 1.0] view was only needed while Z+jets curves fell short of
  // the window, which the lepton overlap removal fixed.)
  const double roc_xmin_default = 0.80, roc_xmax_default = 1.0;
  std::vector<TGraph*> rocs_lo, rocs_hi;
  // Fixed-efficiency working points (see generate_roc_grid): evenly spaced,
  // no breaks across R_pT atoms, and identical x-points across scenarios so
  // the ratio panels are exact.
  const std::vector<double> effGrid = rocEffGrid();
  for (auto& s : scen_lo) rocs_lo.push_back(generate_roc_grid(s.h_pu, s.h_hs, effGrid));
  for (auto& s : scen_hi) rocs_hi.push_back(generate_roc_grid(s.h_pu, s.h_hs, effGrid));
  std::vector<TGraph*> rocs_lo_cen, rocs_hi_cen;
  for (auto& s : scen_lo_cen) rocs_lo_cen.push_back(generate_roc_grid(s.h_pu, s.h_hs, effGrid));
  for (auto& s : scen_hi_cen) rocs_hi_cen.push_back(generate_roc_grid(s.h_pu, s.h_hs, effGrid));

  const double roc_xmin_lo = roc_xmin_default, roc_xmax_lo = roc_xmax_default;
  const double roc_xmin_hi = roc_xmin_default, roc_xmax_hi = roc_xmax_default;

  // ── SECTION 2b: VBS-topology region ROCs ──────────────────────────────────
  // R1 is self-contained: its own forward HS leg vs its own forward PU leg.
  //
  // R2 is NOT. By construction its only forward jet is the PU leg (the HS leg
  // is central, outside HGTD acceptance), so scen_r2[i].h_hs is empty and
  // generate_roc against it would divide by zero. The physics question R2 asks
  // is "how well can timing reject a forward PU jet", which needs a forward-HS
  // reference to trade off against -- and R1's HS leg is exactly that: a
  // genuine forward HS jet under the same track selection and timing gate. So
  // R2's ROC pairs r2.h_pu against r1.h_hs. This is a cross-region ROC and is
  // labelled as such on the plot; do not "fix" it to use r2.h_hs.
  std::vector<TGraph*> rocs_r1, rocs_r2;
  for (auto& s : scen_r1) rocs_r1.push_back(generate_roc_grid(s.h_pu, s.h_hs, effGrid));
  for (size_t i = 0; i < scen_r2.size(); ++i)
    rocs_r2.push_back(generate_roc_grid(scen_r2[i].h_pu, scen_r1[i].h_hs, effGrid));
  // (styled below, once styleRoc is in scope)

  auto styleRoc = [&](TGraph* g, Color_t col, double xmin, double xmax) {
    g->SetTitle("R_{pT} Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->SetLineWidth(2);
    g->GetXaxis()->SetLimits(xmin, xmax);
    g->GetXaxis()->SetNdivisions(810);
    g->SetMinimum(1.0);
  };
  for (size_t i = 0; i < rocs_lo.size(); ++i) {
    styleRoc(rocs_lo[i], scen_lo[i].color, roc_xmin_lo, roc_xmax_lo);
    styleRoc(rocs_hi[i], scen_hi[i].color, roc_xmin_hi, roc_xmax_hi);
  }
  for (size_t i = 0; i < rocs_lo_cen.size(); ++i) {
    styleRoc(rocs_lo_cen[i], scen_lo_cen[i].color, roc_xmin_default, roc_xmax_default);
    styleRoc(rocs_hi_cen[i], scen_hi_cen[i].color, roc_xmin_default, roc_xmax_default);
  }
  for (size_t i = 0; i < rocs_r1.size(); ++i) {
    styleRoc(rocs_r1[i], scen_r1[i].color, roc_xmin_default, roc_xmax_default);
    styleRoc(rocs_r2[i], scen_r2[i].color, roc_xmin_default, roc_xmax_default);
  }

  // ── Console summary ──────────────────────────────────────────────────────────
  // Computed here, BEFORE styleScen() below restricts the R_pT histograms'
  // axis range for display — see the note at styleScen's definition (Section 3).
  std::cout << "\n=== EVENT COUNTS ===\n";
  std::cout << "Total events        : " << n_total      << '\n';
  std::cout << "Pass basic cuts     : " << n_pass_basic << '\n';
  std::cout << "HGTD vtx valid      : " << n_hgtd_valid << " / " << n_pass_basic
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * n_hgtd_valid / n_pass_basic) << "%)\n";

  std::cout << "\n=== ENTRIES / MEAN RpT PER SCENARIO (30-40 GeV) ===\n";
  for (auto& s : scen_lo)
    std::cout << "  " << std::setw(12) << std::left << s.name
              << "  HS: " << (long)s.h_hs->Integral()
              << " <RpT>=" << std::setprecision(3) << s.h_hs->GetMean()
              << "   PU: " << (long)s.h_pu->Integral()
              << " <RpT>=" << s.h_pu->GetMean() << '\n';

  std::cout << "\n=== ENTRIES PER SCENARIO (>40 GeV) ===\n";
  for (auto& s : scen_hi)
    std::cout << "  " << std::setw(12) << std::left << s.name
              << "  HS: " << (long)s.h_hs->Integral()
              << "  PU: " << (long)s.h_pu->Integral() << '\n';

  // ── PU rejection at fixed HS-efficiency working points (>40 GeV slice) ──────
  auto rejAtEff = [](TH1D* PU, TH1D* HS, double targetEff) -> double {
    int bin = HS->GetNbinsX();
    // Explicit bin range — see the note in generate_roc() above.
    double hsTot = HS->Integral(1, bin + 1), puTot = PU->Integral(1, bin + 1);
    for (int i = 1; i <= bin; ++i) {
      double eff = HS->Integral(i, bin + 1) / hsTot;     // decreases with i
      if (eff <= targetEff) {
        double mis = PU->Integral(i, bin + 1) / puTot;
        return mis > 0 ? 1.0 / mis : 0.0;
      }
    }
    return 0.0;
  };
  auto printRejTable = [&](const char* hdr, std::vector<Scenario>& sv) {
    std::cout << "\n=== PU REJECTION (1/mistag) AT FIXED HS EFF " << hdr << " ===\n";
    std::cout << "  scenario       eff=0.85  0.90  0.93  0.95  0.97\n";
    for (auto& s : sv)
      std::printf("  %-12s   %7.1f %6.1f %6.1f %6.1f %6.1f\n", s.name.c_str(),
                  rejAtEff(s.h_pu, s.h_hs, 0.85), rejAtEff(s.h_pu, s.h_hs, 0.90),
                  rejAtEff(s.h_pu, s.h_hs, 0.93), rejAtEff(s.h_pu, s.h_hs, 0.95),
                  rejAtEff(s.h_pu, s.h_hs, 0.97));
  };
  printRejTable("(30-40 GeV, 2.4<|eta|<3.8)", scen_lo);
  printRejTable("(>40 GeV, 2.4<|eta|<3.8)",   scen_hi);
  printRejTable("(30-40 GeV, |eta|<2.4 CENTRAL BASELINE)", scen_lo_cen);
  printRejTable("(>40 GeV, |eta|<2.4 CENTRAL BASELINE)",   scen_hi_cen);
  std::cout << "\n  Central is outside HGTD acceptance, so its four scenarios are\n"
               "  expected to agree to within the few jets whose tracks reach\n"
               "  |eta| > 2.4; divergence there indicates a leak, not a gain.\n";

  std::cout << "\n=== UNTIMED-TRACK FLOOR (ghost & z-selected) ===\n";
  std::printf("  30-40 GeV  PU untimed: %.1f%%   HS untimed: %.1f%%\n",
              pu_tot_lo > 0 ? 100.0 * pu_floor_lo / pu_tot_lo : 0.0,
              hs_tot_lo > 0 ? 100.0 * hs_floor_lo / hs_tot_lo : 0.0);
  std::printf("  >40 GeV    PU untimed: %.1f%%   HS untimed: %.1f%%\n",
              pu_tot_pt > 0 ? 100.0 * pu_floor_pt / pu_tot_pt : 0.0,
              hs_tot_pt > 0 ? 100.0 * hs_floor_pt / hs_tot_pt : 0.0);

  // === END OF DERIVED-NUMBER COMPUTATION. Do not add Integral()/GetMean()
  //     calls below this line without moving them above styleScen(). ===

  // ===========================================================================
  // SECTION 3: Cosmetic styling. styleScen is defined (and immediately
  // called) only here, after every Integral()/GetMean()-based number above
  // has already been computed and printed -- calling it earlier is now a
  // compile error, not just a discipline violation. See the note at
  // generate_roc() above for why the ordering matters: SetRangeUser()
  // silently truncates any subsequent no-arg Integral()/GetMean() call.
  // ===========================================================================
  auto styleScen = [](std::vector<Scenario>& sv) {
    for (auto& s : sv) {
      if (!s.h_hs || !s.h_pu) continue;
      for (auto* h : {s.h_hs, s.h_pu}) {
        h->GetXaxis()->SetRangeUser(0.0, 1.5);
        // 5 primary divisions, not 15: over the [0, 1.5] display window the
        // latter puts a label every 0.1, which run together into an unreadable
        // "0.10.20.30.4..." strip at this canvas width.
        h->GetXaxis()->SetNdivisions(505);
        h->SetLineWidth(2);
      }
      s.h_hs->SetLineColor(s.color);
      s.h_pu->SetLineColor(s.color);
      s.h_pu->SetLineStyle(2);
    }
  };
  styleScen(scen_lo);
  styleScen(scen_hi);
  // Regions too, so their R_pT distributions get the same [0, 1.5] display
  // window. Without this they render against the full booked axis, whose last
  // bin runs to 375 (rpt_bins' single wide overflow bin) -- squashing every
  // real entry into the leftmost 0.4% of the frame. Safe here for the same
  // reason as the two calls above: the region ROCs were built back in
  // SECTION 2b via generate_roc, which uses explicit bin ranges and is immune
  // to SetRangeUser, and the region yield printouts use GetEntries().
  styleScen(scen_lo_cen);
  styleScen(scen_hi_cen);
  styleScen(scen_r1);
  styleScen(scen_r2);

  // ===========================================================================
  // SECTION 4: PDF drawing.
  // ===========================================================================

  // ROC page with ratio-to-zonly panel.
  auto drawRocWithRatio = [&](std::vector<TGraph*>& gs,
                               std::vector<Scenario>& sc,
                               double ymax, double ratio_ymax,
                               double roc_xmin, double roc_xmax,
                               const char* extra_label) {
    canvas->Clear();
    canvas->SetLogy(false);

    TPad* pad_main = new TPad("pad_main", "", 0.0, 0.30, 1.0, 1.0);
    pad_main->SetBottomMargin(0.02);
    pad_main->SetTopMargin(0.07);
    pad_main->SetLeftMargin(0.16);
    pad_main->SetRightMargin(0.05);
    pad_main->Draw();
    pad_main->cd();

    gs[0]->SetMaximum(ymax);
    gs[0]->GetXaxis()->SetLabelSize(0);
    gs[0]->GetXaxis()->SetTitleSize(0);
    gs[0]->Draw("ALP");
    for (size_t i = 1; i < gs.size(); ++i) gs[i]->Draw("LP SAME");

    TLegend* L = new TLegend(0.58, 0.60, 0.93, 0.88);
    StyleLegend(L);
    for (size_t i = 0; i < gs.size(); ++i) L->AddEntry(gs[i], sc[i].legend.c_str());
    L->Draw();

    {
      TLatex t; t.SetNDC(); t.SetTextColor(1);
      t.SetTextFont(72); t.SetTextSize(0.05); t.DrawLatex(0.18, 0.88, "ATLAS");
      t.SetTextFont(42); t.SetTextSize(0.05); t.DrawLatex(0.28, 0.88, "Simulation Internal");
      t.SetTextSize(0.044); t.DrawLatex(0.18, 0.80, MyUtl::ENERGY_LABEL.c_str());
      // Selection annotation, wrapped. The legend is an opaque box whose left
      // edge is at NDC 0.58 and which spans the same y range as this text, so a
      // single long line is painted over by it. Lines are split on "@@" by the
      // caller -- not on a width estimate, since TLatex markup means the glyph
      // count is not the string length. "|" cannot be the delimiter: it already
      // appears in every |#eta| label.
      if (extra_label) {
        t.SetTextSize(0.038);
        const std::string el(extra_label);
        double ly = 0.72;
        for (size_t pos = 0;;) {
          const size_t d = el.find("@@", pos);
          t.DrawLatex(0.18, ly, el.substr(pos, d == std::string::npos ? d : d - pos).c_str());
          ly -= 0.055;
          if (d == std::string::npos) break;
          pos = d + 2;
        }
      }
    }

    canvas->cd();

    TPad* pad_ratio = new TPad("pad_ratio", "", 0.0, 0.0, 1.0, 0.30);
    pad_ratio->SetTopMargin(0.02);
    pad_ratio->SetBottomMargin(0.38);
    pad_ratio->SetLeftMargin(0.16);
    pad_ratio->SetRightMargin(0.05);
    pad_ratio->Draw();
    pad_ratio->cd();

    TGraph* ref = gs[0];
    std::vector<TGraph*> ratios;
    for (size_t i = 0; i < gs.size(); ++i) {
      TGraph* g = gs[i];
      int n = g->GetN();
      std::vector<double> rx, ry;
      for (int j = 0; j < n; ++j) {
        double x = g->GetX()[j];
        if (x < roc_xmin || x > roc_xmax) continue;
        double ref_y = ref->Eval(x);
        if (ref_y <= 0) continue;
        rx.push_back(x);
        ry.push_back(g->GetY()[j] / ref_y);
      }
      TGraph* r = new TGraph(rx.size(), rx.data(), ry.data());
      r->SetLineColor(g->GetLineColor());
      r->SetLineWidth(2);
      ratios.push_back(r);
    }

    TH1F* frame = pad_ratio->DrawFrame(roc_xmin, 0.5, roc_xmax, ratio_ymax,
        ";Hard Scatter Efficiency;Ratio to ITk-only");
    frame->GetXaxis()->SetNdivisions(810);
    frame->GetXaxis()->SetLabelSize(0.13);
    frame->GetXaxis()->SetTitleSize(0.13);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetLabelSize(0.11);
    frame->GetYaxis()->SetTitleSize(0.11);
    frame->GetYaxis()->SetTitleOffset(0.65);
    frame->GetYaxis()->SetNdivisions(504);

    TLine* line1 = new TLine(roc_xmin, 1.0, roc_xmax, 1.0);
    line1->SetLineColor(kGray + 1);
    line1->SetLineStyle(2);
    line1->Draw();
    for (auto* r : ratios) r->Draw("L SAME");

    canvas->cd();
    canvas->Print(out_pdf);
  };

  // ROC y-maximum, auto-scaled to the curves actually visible in the x-window.
  // This was previously hard-coded to 300, which only suited the old [0.6, 1.0]
  // view: pile-up rejection falls steeply with efficiency, so a zoomed
  // high-efficiency window leaves every curve squashed into the bottom of a
  // fixed 0-300 frame. Computed once across BOTH pT slices so they keep a
  // shared, directly-comparable scale (as the fixed 300 did). The 1.5 headroom
  // keeps the topmost curve near ~2/3 pad height, clear of the ATLAS/sample
  // labels. Reads the graphs via GetPoint, so it is unaffected by any
  // SetRangeUser applied for display.
  auto rocYMaxIn = [](const std::vector<TGraph*>& gs, double xmin, double xmax) {
    double m = 0.0;
    for (auto* g : gs)
      for (int i = 0; i < g->GetN(); ++i) {
        double x, y;
        g->GetPoint(i, x, y);
        if (x >= xmin && x <= xmax) m = std::max(m, y);
      }
    return m;
  };
  double roc_ymax = 1.5 * std::max(rocYMaxIn(rocs_lo, roc_xmin_lo, roc_xmax_lo),
                                   rocYMaxIn(rocs_hi, roc_xmin_hi, roc_xmax_hi));
  if (!(roc_ymax > 0.0)) roc_ymax = 300.0;  // fallback: no points in window

  // (1) ROC — 30–40 GeV.  Linear, shared y maxima across slices; ratio ymax = 4.
  drawRocWithRatio(rocs_lo, scen_lo, roc_ymax, 4.0, roc_xmin_lo, roc_xmax_lo, lbl_lo);

  // (2) ROC — >40 GeV.
  drawRocWithRatio(rocs_hi, scen_hi, roc_ymax, 4.0, roc_xmin_hi, roc_xmax_hi, lbl_hi);

  // (3) ROC — VBS region R1: both candidate legs forward, HS leg vs PU leg.
  // (4) ROC — VBS region R2: forward PU leg, signal side borrowed from R1
  //     (R2 has no forward HS jet by construction -- see SECTION 2b).
  // Placed directly after the inclusive ROCs so the four discrimination plots
  // sit together at the front of the PDF. drawRocWithRatio Print()s the canvas
  // itself -- do not add another Print here.
  // (2b) ROC — central baseline, both pT slices. Own y-scale: central rejection
  //      sits well below forward, so sharing the forward maximum would squash
  //      these into the axis.
  double roc_ymax_cen = 1.5 * std::max(rocYMaxIn(rocs_lo_cen, roc_xmin_default, roc_xmax_default),
                                       rocYMaxIn(rocs_hi_cen, roc_xmin_default, roc_xmax_default));
  drawRocWithRatio(rocs_lo_cen, scen_lo_cen, roc_ymax_cen, 4.0,
                   roc_xmin_default, roc_xmax_default, lbl_lo_cen);
  drawRocWithRatio(rocs_hi_cen, scen_hi_cen, roc_ymax_cen, 4.0,
                   roc_xmin_default, roc_xmax_default, lbl_hi_cen);

  // R1/R2 get their own y-scale rather than the inclusive slices'. With the
  // truth-t0 and clean-timing rows added, both regions overshoot the inclusive
  // maximum badly -- R1's clean-timing peaks around 880 and R2's truth around
  // 1180 against an inclusive ceiling near 640 -- so sharing it clipped the
  // top off both panels.
  double roc_ymax_reg = 1.5 * std::max(rocYMaxIn(rocs_r1, roc_xmin_default, roc_xmax_default),
                                       rocYMaxIn(rocs_r2, roc_xmin_default, roc_xmax_default));
  drawRocWithRatio(rocs_r1, scen_r1, roc_ymax_reg, 4.0,
                   roc_xmin_default, roc_xmax_default, lbl_r1);
  drawRocWithRatio(rocs_r2, scen_r2, roc_ymax_reg, 4.0,
                   roc_xmin_default, roc_xmax_default, lbl_r2);

  // (3+) Per-scenario HS vs PU, 30–40 GeV slice, log-y.
  canvas->Clear();
  canvas->SetLogy(true);
  for (auto& s : scen_lo) {
    TLegend* L = new TLegend(0.60, 0.75, 0.92, 0.88);
    StyleLegend(L);
    L->AddEntry(s.h_hs, "Hard Scatter");
    L->AddEntry(s.h_pu, "Pile-Up");
    double ymax = 50.0 * std::max(s.h_hs->GetMaximum(), s.h_pu->GetMaximum());
    s.h_pu->SetMaximum(ymax);
    s.h_pu->SetTitle((std::string("R_{pT}: ") + s.legend).c_str());
    s.h_pu->Draw("HIST");
    s.h_hs->Draw("HIST SAME");
    drawLabels(lbl_lo);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // All-scenario HS overlay — 30–40 GeV.
  {
    TLegend* L = new TLegend(0.55, 0.68, 0.92, 0.88);
    StyleLegend(L);
    double ymax = 0;
    for (auto& s : scen_lo) ymax = std::max(ymax, s.h_hs->GetMaximum());
    scen_lo[0].h_hs->SetMaximum(50.0 * ymax);
    scen_lo[0].h_hs->SetTitle("Hard Scatter R_{pT} — all scenarios");
    for (size_t i = 0; i < scen_lo.size(); ++i) {
      scen_lo[i].h_hs->SetLineStyle(1);
      scen_lo[i].h_hs->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(scen_lo[i].h_hs, scen_lo[i].legend.c_str());
    }
    drawLabels(lbl_lo);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // All-scenario PU overlay — 30–40 GeV.
  {
    TLegend* L = new TLegend(0.55, 0.68, 0.92, 0.88);
    StyleLegend(L);
    double ymax = 0;
    for (auto& s : scen_lo) ymax = std::max(ymax, s.h_pu->GetMaximum());
    scen_lo[0].h_pu->SetMaximum(50.0 * ymax);
    scen_lo[0].h_pu->SetTitle("Pile-Up R_{pT} — all scenarios");
    for (size_t i = 0; i < scen_lo.size(); ++i) {
      scen_lo[i].h_pu->SetLineStyle(1);
      scen_lo[i].h_pu->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(scen_lo[i].h_pu, scen_lo[i].legend.c_str());
    }
    drawLabels(lbl_lo);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // R_pT distributions for the CENTRAL baseline, same treatment as the forward
  // slices above: per-scenario HS-vs-PU, then the all-scenario HS and PU
  // overlays. The four central scenarios are expected to sit on top of one
  // another (|eta| < 2.4 is outside HGTD acceptance, so the time gate is a
  // no-op there) -- the overlays are where that is easiest to see, which is
  // why they are worth drawing even though they carry no timing gain.
  auto drawShapes = [&](std::vector<Scenario>& sv, const char* lbl) {
    for (auto& s : sv) {
      TLegend* L = new TLegend(0.60, 0.75, 0.92, 0.88);
      StyleLegend(L);
      L->AddEntry(s.h_hs, "Hard Scatter");
      L->AddEntry(s.h_pu, "Pile-Up");
      s.h_pu->SetMaximum(50.0 * std::max(s.h_hs->GetMaximum(), s.h_pu->GetMaximum()));
      s.h_pu->SetTitle((std::string("R_{pT}: ") + s.legend).c_str());
      s.h_pu->Draw("HIST");
      s.h_hs->Draw("HIST SAME");
      drawLabels(lbl);
      L->Draw();
      canvas->Print(out_pdf);
    }
    for (int pass = 0; pass < 2; ++pass) {          // 0 = HS overlay, 1 = PU overlay
      TLegend* L = new TLegend(0.55, 0.68, 0.92, 0.88);
      StyleLegend(L);
      auto pick = [&](Scenario& s) { return pass == 0 ? s.h_hs : s.h_pu; };
      double ymax = 0;
      for (auto& s : sv) ymax = std::max(ymax, pick(s)->GetMaximum());
      pick(sv[0])->SetMaximum(50.0 * ymax);
      pick(sv[0])->SetTitle(pass == 0 ? "Hard Scatter R_{pT} — all scenarios"
                                      : "Pile-Up R_{pT} — all scenarios");
      for (size_t i = 0; i < sv.size(); ++i) {
        pick(sv[i])->SetLineStyle(1);
        pick(sv[i])->Draw(i == 0 ? "HIST" : "HIST SAME");
        L->AddEntry(pick(sv[i]), sv[i].legend.c_str());
      }
      drawLabels(lbl);
      L->Draw();
      canvas->Print(out_pdf);
    }
  };
  drawShapes(scen_lo_cen, lbl_lo_cen);
  drawShapes(scen_hi_cen, lbl_hi_cen);

  // ===========================================================================
  // SECTION 5: VBS-topology regions.
  //   R1  both VBS candidate legs forward, one truth-HS + one truth-PU --
  //       timing has to pick which forward jet is the hard-scatter one.
  //   R2  forward truth-PU leg + central truth-HS leg -- can timing reject the
  //       forward PU jet? Its ROC's signal side is R1's forward HS leg (R2 has
  //       no forward HS jet by construction); see SECTION 2b.
  // ===========================================================================
  // (The two region ROCs are drawn earlier, as pages 3-4 alongside the
  // inclusive ROCs.) What remains here is the per-region R_pT shapes.
  auto drawRegionShapes = [&](std::vector<Scenario>& sc, bool useHS,
                              const char* title, const char* lbl) {
    double ymax = 0.0;
    for (auto& s : sc) ymax = std::max(ymax, (useHS ? s.h_hs : s.h_pu)->GetMaximum());
    if (ymax <= 0.0) return;  // nothing filled (e.g. R2's empty HS side)
    canvas->Clear();
    canvas->SetLogy(true);
    // Bottom-right: the region labels are long enough to run under a top-right
    // legend, and these distributions fall off well before the right edge, so
    // that corner is empty.
    TLegend* L = new TLegend(0.58, 0.16, 0.93, 0.40);
    StyleLegend(L);
    for (size_t i = 0; i < sc.size(); ++i) {
      TH1D* h = useHS ? sc[i].h_hs : sc[i].h_pu;
      h->SetLineStyle(1);
      if (i == 0) { h->SetMaximum(50.0 * ymax); h->SetTitle(title); }
      h->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(h, sc[i].legend.c_str());
    }
    drawLabels(lbl);
    L->Draw();
    canvas->Print(out_pdf);
  };
  drawRegionShapes(scen_r1, true,  "R1 forward HS leg R_{pT};R_{pT};Entries", lbl_r1);
  drawRegionShapes(scen_r1, false, "R1 forward PU leg R_{pT};R_{pT};Entries", lbl_r1);
  drawRegionShapes(scen_r2, false, "R2 forward PU leg R_{pT};R_{pT};Entries", lbl_r2);

  // Region yields -- these topologies are rare, so print the counts that the
  // ROCs above are built from. A region with very few entries makes its ROC
  // unreliable no matter how clean the curve looks.
  std::printf("\n=== VBS-topology region yields (zonly scenario) ===\n");
  std::printf("  R1 forward HS legs : %8.0f\n", scen_r1[0].h_hs->GetEntries());
  std::printf("  R1 forward PU legs : %8.0f\n", scen_r1[0].h_pu->GetEntries());
  std::printf("  R2 forward PU legs : %8.0f\n", scen_r2[0].h_pu->GetEntries());
  std::printf("  (R2 HS side is empty by construction; its ROC uses R1's HS legs)\n");

  canvas->Print(out_pdf + "]");

  std::cout << "\nWrote " << out_pdf << "\n";
  return 0;
}
