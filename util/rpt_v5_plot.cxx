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
//      pickXRange's window, the console summary, the rejection table).
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
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;

  // ===========================================================================
  // SECTION 1: Load histograms (raw, unstyled) from the saved ROOT file.
  // ===========================================================================
  const std::string histFile = MyUtl::resolveHistFile(
    argc, argv, MyUtl::OUTPUT_DIR + "/hists/rpt_v5_hist.root");
  std::cout << "Reading histograms from " << histFile << '\n';

  std::vector<Scenario> scen_lo = makeScenarios("_lo");
  std::vector<Scenario> scen_hi = makeScenarios("_hi");
  long n_total, n_pass_basic, n_hgtd_valid;
  double pu_tot_pt, pu_floor_pt, hs_tot_pt, hs_floor_pt;
  double pu_tot_lo, pu_floor_lo, hs_tot_lo, hs_floor_lo;
  {
    MyUtl::HistReader reader(histFile);
    loadScenarios(reader, scen_lo);
    loadScenarios(reader, scen_hi);
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
  const TString out_pdf = TString::Format("%s/rpt_plots/rpt_v5.pdf", MyUtl::OUTPUT_DIR.c_str());

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
  const char* lbl_lo = "30 < p_{T}^{jet} < 40 GeV, 2.4 < |#eta| < 3.8";
  const char* lbl_hi = "p_{T}^{jet} > 40 GeV, 2.4 < |#eta| < 3.8";

  // ===========================================================================
  // SECTION 2: Derive every Integral()/GetMean()-based number. styleScen is
  // NOT defined yet -- see the sentinel comment at the end of this section.
  // ===========================================================================

  // Preferred "high-efficiency working point" window, as used in the paper.
  const double roc_xmin_default = 0.8, roc_xmax_default = 1.0;
  std::vector<TGraph*> rocs_lo, rocs_hi;
  for (auto& s : scen_lo) rocs_lo.push_back(generate_roc(s.h_pu, s.h_hs));
  for (auto& s : scen_hi) rocs_hi.push_back(generate_roc(s.h_pu, s.h_hs));

  // A "high-efficiency working point" plot is meaningless below 50% HS
  // efficiency (e.g. Z+jets forward "HS" jets with ~no tracks can otherwise
  // drag the fallback window's lower edge down arbitrarily far).
  static constexpr double MIN_HS_EFF = 0.5;

  // Some samples (e.g. Z+jets, where a large fraction of forward "HS" jets
  // carry ~no tracks) collapse the ROC curve to efficiencies well below the
  // default window, leaving zero points in [roc_xmin_default, roc_xmax_default]
  // and an empty page.  Detect that per pT-slice and fall back to the actual
  // observed x-range instead of hardcoding a sample-specific override.
  auto pickXRange = [](const std::vector<TGraph*>& gs, double lo_default, double hi_default) {
    double data_lo = 1.0, data_hi = 0.0, window_max = 0.0;
    bool any_in_window = false;
    for (auto* g : gs) {
      int n = g->GetN();
      for (int j = 0; j < n; ++j) {
        double x = g->GetX()[j];
        data_lo = std::min(data_lo, x);
        data_hi = std::max(data_hi, x);
        if (x >= lo_default && x <= hi_default) {
          any_in_window = true;
          window_max = std::max(window_max, x);
        }
      }
    }
    // Requiring merely "one point somewhere in the window" isn't enough: a
    // spike near R_pT=0 can make efficiency crash from ~1.0 to e.g. 0.85 in
    // a single bin step, leaving only a stub hugging the window's low edge
    // (e.g. zjets 30-40 GeV, cut off at eff~0.89) — technically non-empty,
    // but the "high-efficiency working point" framing is meaningless there.
    // Require the curve to actually reach near the window's top edge too.
    bool reaches_near_top = window_max >= hi_default - 0.05;
    if ((any_in_window && reaches_near_top) || data_hi <= data_lo)
      return std::make_pair(lo_default, hi_default);
    double pad = 0.02 * (data_hi - data_lo);
    double lo  = std::max(MIN_HS_EFF, data_lo - pad);
    return std::make_pair(lo, std::min(1.0, data_hi + pad));
  };
  auto [roc_xmin_lo, roc_xmax_lo] = pickXRange(rocs_lo, roc_xmin_default, roc_xmax_default);
  auto [roc_xmin_hi, roc_xmax_hi] = pickXRange(rocs_hi, roc_xmin_default, roc_xmax_default);

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
  printRejTable("(30-40 GeV)", scen_lo);
  printRejTable("(>40 GeV)",   scen_hi);

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
        h->GetXaxis()->SetNdivisions(515);
        h->SetLineWidth(2);
      }
      s.h_hs->SetLineColor(s.color);
      s.h_pu->SetLineColor(s.color);
      s.h_pu->SetLineStyle(2);
    }
  };
  styleScen(scen_lo);
  styleScen(scen_hi);

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
      if (extra_label) { t.SetTextSize(0.044); t.DrawLatex(0.18, 0.72, extra_label); }
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

  // (1) ROC — 30–40 GeV.  Linear, shared y maxima: ROC ymax = 300, ratio ymax = 4.
  drawRocWithRatio(rocs_lo, scen_lo, 300.0, 4.0, roc_xmin_lo, roc_xmax_lo, lbl_lo);

  // (2) ROC — >40 GeV.
  drawRocWithRatio(rocs_hi, scen_hi, 300.0, 4.0, roc_xmin_hi, roc_xmax_hi, lbl_hi);

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

  canvas->Print(out_pdf + "]");

  std::cout << "\nWrote " << out_pdf << "\n";
  return 0;
}
