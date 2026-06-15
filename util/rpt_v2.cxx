// rpt_v2.cxx — RpT discrimination study using the main-analysis event selection.
//
// Five scenarios (HS vs PU RpT, ROC):
//   1. zonly       — z-significance tracks only, no time gate
//   2. hgtd        — ntuple RecoVtx_time, RecoVtx_timeRes  (per-track gate; no event-level recoVtxValid gate)
//   3. mine        — TRKPTZ-selected cluster time
//   4. mine_pure   — same time as #3, only events with TRKPTZ cluster purity > 0.75   (TEST_MISCL gate)
//   5. mine_misas  — same time as #3, only events with HS timing purity >= 0.95       (TEST_MISAS gate)
//
// Event/track selection mirrors clustering_dt (passBasicCuts + passJetPtCut +
// getAssociatedTracks).  Jets for RpT are forward + above MIN_JET_PT.
//
// Output: figs/rpt_plots/rpt_v2.pdf

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <limits>
#include <string>
#include <vector>

#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "clustering_constants.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

#define debug false

using namespace MyUtl;

// -----------------------------------------------------------------------------
// ROC: HS efficiency vs PU rejection (1 / mistag).  Copied verbatim from
// util/generate_rpt.cxx so the new script reproduces the historical curve shape.
// -----------------------------------------------------------------------------
TGraph* generate_roc(TH1D* PU_hist, TH1D* HS_hist) {
  int bin = PU_hist->GetNbinsX();
  std::vector<float> vx, vy;
  for (int i = 1; i <= bin; ++i) {
    double HS_eff   = HS_hist->Integral(i, bin + 1) / HS_hist->Integral();
    double PU_mistag = PU_hist->Integral(i, bin + 1) / PU_hist->Integral();
    if (std::abs(PU_mistag) > 1e-6 && HS_eff < 0.99) {
      vx.push_back(HS_eff);
      vy.push_back(1.0 / PU_mistag);
    }
  }
  return new TGraph((int)vx.size(), vx.data(), vy.data());
}

// -----------------------------------------------------------------------------
// Scenario container — five histograms in parallel.
// -----------------------------------------------------------------------------
struct Scenario {
  std::string name;
  std::string legend;
  Color_t color;
  TH1D* h_hs;
  TH1D* h_pu;
};

// Cone-style ΔR association used by the original RpT computation.
static inline double dR(double j_eta, double j_phi, double t_eta, double t_phi) {
  double deta = j_eta - t_eta;
  double dphi = TVector2::Phi_mpi_pi(j_phi - t_phi);
  return std::sqrt(deta * deta + dphi * dphi);
}

// Compute RpT for one jet given a track index set (numerator includes any
// track with dR < 0.3, regardless of pT — same convention as generate_rpt).
static double computeRpT(BranchPointerWrapper* b,
                         const std::vector<int>& ghost_indices,
                         double j_pt,
                         const std::unordered_set<int>& assoc_set) {
  double sumpt = 0.0;
  for (int idx : ghost_indices)
    if (assoc_set.count(idx)) sumpt += b->trackPt[idx];
  return sumpt / j_pt;
}

int main() {
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found.  Aborting.\n";
    return 1;
  }

  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  // Paper jet-label helper (ATL-HGTD-PUB-2022-001 Sec. 3).
  // HS  : dR(reco, truthHS) < 0.3  AND  truthHS_pT > 10 GeV
  // PU  : dR > 0.6 from ALL truth jets (HS + ITPU + OOTPU) with pT > 4 GeV
  auto paperIsHS = [&](double j_eta, double j_phi) {
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      if (branch.truthHSJetPt[t] < 10.0) continue;
      if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.3)
        return true;
    }
    return false;
  };
  // PU: dR > 0.6 from any truth HS jet with pT > 4 GeV (Sec. 3).
  // Only HS truth jets enter — ITPU/OOTPU are irrelevant for this cut.
  auto paperIsPU = [&](double j_eta, double j_phi) {
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      if (branch.truthHSJetPt[t] < 4.0) continue;
      if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.6)
        return false;
    }
    return true;
  };

  // Non-uniform binning: 0.01-wide bins in [0, 2.5] for ROC granularity,
  // then one wide bin [2.5, 375] to capture the tail without bloating memory.
  static std::vector<double> rpt_bins = []() {
    std::vector<double> b;
    b.reserve(202);
    for (int i = 0; i <= 200; ++i) b.push_back(0.0125 * i);
    b.push_back(375.0);
    return b;
  }();
  const int nbin = (int)rpt_bins.size() - 1;  // 201

  auto makeHist = [&](const char* name, const char* title) {
    return new TH1D(name, title, nbin, rpt_bins.data());
  };

  auto makeScenarios = [&](const std::string& suffix) {
    std::vector<Scenario> s = {
      {"zonly",       "z-only (no timing)",              C05, nullptr, nullptr},
      {"hgtd",        "HGTD ntuple t_{0}",               C01, nullptr, nullptr},
      {"mine",        "My algo (TRKPTZ)",                   C02, nullptr, nullptr},
      {"mine_pure",   "My algo (TRKPTZ) + pure cluster",  C03, nullptr, nullptr},
      {"mine_misas",  "My algo (TRKPTZ) + clean timing",  C04, nullptr, nullptr},
    };
    for (auto& sc : s) {
      sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                         ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
      sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                         ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
    }
    return s;
  };

  std::vector<Scenario> scen      = makeScenarios("");        // pT > 30 GeV
  std::vector<Scenario> scen_hipt = makeScenarios("_hipt");   // pT > 50 GeV

  // Flagging vs Correcting scenarios (advisor request).
  auto makeFlaggedScenarios = [&](const std::string& suffix) {
    std::vector<Scenario> s = {
      {"zonly_f",     "z-only (no timing)",                              C05, nullptr, nullptr},
      {"correcting",  "My algo + clean timing (correcting showers)",     C04, nullptr, nullptr},
      {"flagging",    "My algo + clean timing (flagging showers)",       C02, nullptr, nullptr},
    };
    for (auto& sc : s) {
      sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                         ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
      sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                         ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
    }
    return s;
  };
  std::vector<Scenario> flag      = makeFlaggedScenarios("_f");      // pT > 30 GeV
  std::vector<Scenario> flag_hipt = makeFlaggedScenarios("_fhipt");  // pT > 50 GeV

  // Diagnostic: distribution of (t_mine − t_hgtd) per selected event,
  // filled only when both vertex times are available.
  TH1D* h_dt_diag = new TH1D("h_dt_diag",
      "TRKPTZ #minus HGTD vertex time;t_{mine} #minus t_{HGTD} [ps];Events",
      500, -500.0, 500.0);
  h_dt_diag->SetLineColor(C01);
  h_dt_diag->SetLineWidth(2);

  // Diagnostic: vertex time resolution — HGTD reco vs TRKPTZ cluster sigma.
  TH1D* h_res_hgtd = new TH1D("h_res_hgtd",
      "Vertex time resolution;#sigma_{vtx} [ps];Events",
      200, 0.0, 40.0);
  h_res_hgtd->SetLineColor(C01);
  h_res_hgtd->SetLineWidth(2);
  TH1D* h_res_mine = new TH1D("h_res_mine", "", 200, 0.0, 40.0);
  h_res_mine->SetLineColor(C02);
  h_res_mine->SetLineWidth(2);

  long n_total = 0, n_pass_sel = 0, n_pure_evt = 0, n_misas_evt = 0;

  while (reader.Next()) {
    ++n_total;

    // ── Event selection: same as clustering_dt. ─────────────────────────────
    if (!branch.passBasicCuts())  continue;
    if (!branch.passJetPtCut())   continue;
    ++n_pass_sel;

    // ── Track selection: same as clustering_dt (3σ z-sig at PV). ────────────
    std::vector<int> trk_z = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 2.5);

    // ── TRKPTZ clustering for scenarios 3/4/5 (purity flag on for #4). ──────
    auto clusters = clusterTracksInTime(
        trk_z, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
        /*sortTracks=*/false, /*calcPurityFlag=*/true);

    // TRKPTZ selection: highest-score cluster (no oracle gate).
    double t_trkptz = 0.0, var_trkptz = 0.0;
    float  cluster_purity = 0.0f;
    bool   trkptz_ok = false;
    if (!clusters.empty()) {
      auto best      = chooseCluster(clusters, Score::TRKPTZ);
      t_trkptz       = best.values[0];
      var_trkptz     = best.sigmas[0] * best.sigmas[0];
      cluster_purity = best.purity;
      trkptz_ok      = true;
    }

    // ── HGTD ntuple vertex time. ────────────────────────────────────────────
    double t_hgtd       = branch.recoVtxTime[0];
    double var_hgtd     = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
    bool   hgtd_vtx_valid = (branch.recoVtxValid[0] == 1);

    // ── Diagnostics. ─────────────────────────────────────────────────────────
    if (trkptz_ok && hgtd_vtx_valid)
      h_dt_diag->Fill(t_trkptz - t_hgtd);
    if (hgtd_vtx_valid)
      h_res_hgtd->Fill(branch.recoVtxTimeRes[0]);
    if (trkptz_ok)
      h_res_mine->Fill(std::sqrt(var_trkptz));

    // ── Event-level oracle gates for scenarios 4 / 5 (TRKPTZ-based). ────────
    bool gate_pure  = trkptz_ok && (cluster_purity > 0.50f);          // TEST_MISCL
    float hs_tp     = calcHSTimingPurity(trk_z, &branch);
    bool gate_misas = (hs_tp >= 0.9f);                                // TEST_MISAS
    if (gate_pure)  ++n_pure_evt;
    if (gate_misas) ++n_misas_evt;

    // ── Per-scenario track lists.  Z-sig is the baseline; timing scenarios
    //    apply a 3σ time gate only when both the vertex time AND the track
    //    time are valid.  If either is invalid, the track falls through to
    //    the z-only inclusion (matches generate_rpt.cxx:275 behavior). ─────
    auto applyTimeGate = [&](double t_vtx, double var_vtx, bool vtx_valid) {
      std::vector<int> out;
      out.reserve(trk_z.size());
      for (int idx : trk_z) {
        bool apply = vtx_valid && branch.trackTimeValid[idx] == 1;
        if (!apply) { out.push_back(idx); continue; }
        double dt    = branch.trackTime[idx] - t_vtx;
        double var_t = branch.trackTimeRes[idx] * branch.trackTimeRes[idx];
        // pull width ~1.5 from fit → inflate var_vtx by 1.5² = 2.25
        double pull  = std::abs(dt) / std::sqrt(2.25 * var_vtx + var_t);
        if (pull < 2.0) out.push_back(idx);
      }
      return out;
    };

    std::vector<int> trk_hgtd   = applyTimeGate(t_hgtd,   var_hgtd,   hgtd_vtx_valid);
    std::vector<int> trk_trkptz = applyTimeGate(t_trkptz, var_trkptz, trkptz_ok);

    // Build per-scenario sets once per event for O(1) ghost-index lookup.
    std::unordered_set<int> set_z     (trk_z.begin(),      trk_z.end());
    std::unordered_set<int> set_hgtd  (trk_hgtd.begin(),   trk_hgtd.end());
    std::unordered_set<int> set_trkptz(trk_trkptz.begin(), trk_trkptz.end());

    // Fill RpT histograms for one scenario set at a given min jet pT.
    auto fillScenario = [&](Scenario& s, const std::unordered_set<int>& s_set,
                            bool eventOk, double min_pt) {
      if (!eventOk) return;
      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
        double j_pt  = branch.topoJetPt[j];
        double j_eta = branch.topoJetEta[j];
        double j_phi = branch.topoJetPhi[j];
        if (j_pt < min_pt) continue;
        if (std::abs(j_eta) < MIN_ABS_ETA_JET || std::abs(j_eta) > MAX_ABS_ETA_JET) continue;
        bool isHS = paperIsHS(j_eta, j_phi);
        bool isPU = paperIsPU(j_eta, j_phi);
        if (!isHS && !isPU) continue;
        double r = computeRpT(&branch, branch.topoJetGhostTrackIdx[j], j_pt, s_set);
        if (isHS) s.h_hs->Fill(r);
        else      s.h_pu->Fill(r);
      }
    };

    // pT > 30 GeV (inclusive)
    fillScenario(scen[0], set_z,       true,       MIN_JET_PT);
    fillScenario(scen[1], set_hgtd,    true,       MIN_JET_PT);
    fillScenario(scen[2], set_trkptz,  true,       MIN_JET_PT);
    fillScenario(scen[3], set_trkptz,  gate_pure,  MIN_JET_PT);
    fillScenario(scen[4], set_trkptz,  gate_misas, MIN_JET_PT);

    // pT > 50 GeV (high-pT check)
    fillScenario(scen_hipt[0], set_z,      true,       50.0);
    fillScenario(scen_hipt[1], set_hgtd,   true,       50.0);
    fillScenario(scen_hipt[2], set_trkptz, true,       50.0);
    fillScenario(scen_hipt[3], set_trkptz, gate_pure,  50.0);
    fillScenario(scen_hipt[4], set_trkptz, gate_misas, 50.0);

    // Flagging vs Correcting scenarios.
    // sv[0]=z-only  sv[1]=correcting(mine_misas)  sv[2]=flagging(same gate, z-only tracks)
    fillScenario(flag[0],      set_z,      true,       MIN_JET_PT);
    fillScenario(flag[1],      set_trkptz, gate_misas, MIN_JET_PT);  // correcting
    fillScenario(flag[2],      set_z,      gate_misas, MIN_JET_PT);  // flagging
    fillScenario(flag_hipt[0], set_z,      true,       50.0);
    fillScenario(flag_hipt[1], set_trkptz, gate_misas, 50.0);        // correcting
    fillScenario(flag_hipt[2], set_z,      gate_misas, 50.0);        // flagging
  }

  // ── Plotting ───────────────────────────────────────────────────────────────
  boost::filesystem::create_directories("../figs/rpt_plots");
  const TString out_pdf = "../figs/rpt_plots/rpt_v2.pdf";

  TCanvas* canvas = new TCanvas("canvas", "RpT", 800, 700);
  canvas->Print(out_pdf + "[");

  // ATLAS label helper: "ATLAS Simulation Internal" + energy/channel line +
  // optional third line for the jet selection.
  auto drawLabels = [](const char* extra = nullptr) {
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82);  // "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv."
    if (extra) {
      TLatex t;
      t.SetNDC();
      t.SetTextFont(42);
      t.SetTextSize(0.032);
      t.DrawLatex(0.18, 0.76, extra);
    }
  };
  const char* jetcut    = "p_{T}^{jet} > 30 GeV, 2.38 < |#eta| < 4.0";
  const char* jetcut_hi = "p_{T}^{jet} > 50 GeV, 2.38 < |#eta| < 4.0";

  // Common cosmetics helper.
  auto styleScen = [](std::vector<Scenario>& sv) {
    for (auto& s : sv) {
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
  styleScen(scen);
  styleScen(scen_hipt);

  // Build ROC graphs for both pT slices now (histos already filled).
  const double roc_xmin = 0.8, roc_xmax = 1.0;
  std::vector<TGraph*> rocs, rocs_hipt;
  for (auto& s : scen)      rocs.push_back(generate_roc(s.h_pu, s.h_hs));
  for (auto& s : scen_hipt) rocs_hipt.push_back(generate_roc(s.h_pu, s.h_hs));

  std::vector<TGraph*> rocs_f, rocs_fhipt;
  for (auto& s : flag)      rocs_f.push_back(generate_roc(s.h_pu, s.h_hs));
  for (auto& s : flag_hipt) rocs_fhipt.push_back(generate_roc(s.h_pu, s.h_hs));

  auto styleRoc = [&](TGraph* g, Color_t col, const char* title) {
    g->SetTitle(title);
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->SetLineWidth(2);
    g->GetXaxis()->SetLimits(roc_xmin, roc_xmax);
    g->GetXaxis()->SetNdivisions(810);
    g->SetMinimum(1.0);
    g->SetMaximum(180);
  };
  for (size_t i = 0; i < rocs.size(); ++i) {
    styleRoc(rocs[i],      scen[i].color,
             "R_{pT} Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    styleRoc(rocs_hipt[i], scen_hipt[i].color,
             "R_{pT} Discriminant (p_{T} > 50 GeV);Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    rocs[i]->SetMaximum(300);
    rocs_hipt[i]->SetMaximum(500);
  }

  // Helper: draw one ROC page with a ratio-to-zonly panel below.
  auto drawRocWithRatio = [&](std::vector<TGraph*>& gs,
                               std::vector<Scenario>& sc,
                               double ymax,
                               double ratio_ymax,
                               const char* extra_label) {
    canvas->Clear();
    canvas->SetLogy(false);

    // Main pad: 70% of height.
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

    TLegend* L = new TLegend(0.58, 0.52, 0.93, 0.88);
    StyleLegend(L);
    for (size_t i = 0; i < gs.size(); ++i) L->AddEntry(gs[i], sc[i].legend.c_str());
    L->Draw();

    // Draw ATLAS labels manually to avoid the gap formula misbehaving in sub-pads.
    {
      TLatex t; t.SetNDC(); t.SetTextColor(1);
      t.SetTextFont(72); t.SetTextSize(0.05); t.DrawLatex(0.18, 0.88, "ATLAS");
      t.SetTextFont(42); t.SetTextSize(0.05); t.DrawLatex(0.28, 0.88, "Simulation Internal");
      t.SetTextSize(0.044); t.DrawLatex(0.18, 0.80, "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.");
      if (extra_label) { t.SetTextSize(0.044); t.DrawLatex(0.18, 0.72, extra_label); }
    }

    canvas->cd();

    // Ratio pad: bottom 30%.
    TPad* pad_ratio = new TPad("pad_ratio", "", 0.0, 0.0, 1.0, 0.30);
    pad_ratio->SetTopMargin(0.02);
    pad_ratio->SetBottomMargin(0.38);
    pad_ratio->SetLeftMargin(0.16);
    pad_ratio->SetRightMargin(0.05);
    pad_ratio->Draw();
    pad_ratio->cd();

    // Build ratio graphs: each curve divided by z-only at same x.
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

    // Draw ratio frame.
    TH1F* frame = pad_ratio->DrawFrame(roc_xmin, 0.5, roc_xmax, ratio_ymax,
        ";Hard Scatter Efficiency;Ratio to z-only");
    frame->GetXaxis()->SetNdivisions(810);
    frame->GetXaxis()->SetLabelSize(0.13);
    frame->GetXaxis()->SetTitleSize(0.13);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetLabelSize(0.11);
    frame->GetYaxis()->SetTitleSize(0.11);
    frame->GetYaxis()->SetTitleOffset(0.65);
    frame->GetYaxis()->SetNdivisions(504);

    // Reference line at 1.
    TLine* line1 = new TLine(roc_xmin, 1.0, roc_xmax, 1.0);
    line1->SetLineColor(kGray + 1);
    line1->SetLineStyle(2);
    line1->Draw();

    for (auto* r : ratios) r->Draw("L SAME");

    canvas->cd();
    canvas->Print(out_pdf);
  };

  // (1) ROC — pT > 30 GeV.
  drawRocWithRatio(rocs,      scen,      300.0, 4.0, jetcut);

  // (2) ROC — pT > 50 GeV.
  drawRocWithRatio(rocs_hipt, scen_hipt, 500.0, 5.0, jetcut_hi);

  // (3) Diagnostic: t_mine − t_hgtd.
  h_dt_diag->Draw("HIST");
  drawLabels();
  canvas->Print(out_pdf);

  // (3b) Diagnostic: vertex time resolution — HGTD reco vs TRKPTZ cluster sigma.
  {
    canvas->SetLogy(false);
    double ymax = 1.3 * std::max(h_res_hgtd->GetMaximum(), h_res_mine->GetMaximum());
    h_res_hgtd->SetMaximum(ymax);
    h_res_hgtd->Draw("HIST");
    h_res_mine->Draw("HIST SAME");
    TLegend* L = new TLegend(0.55, 0.72, 0.92, 0.84);
    StyleLegend(L);
    L->AddEntry(h_res_hgtd, "HGTD reco #sigma_{vtx}");
    L->AddEntry(h_res_mine, "TRKPTZ cluster #sigma_{vtx}");
    drawLabels();
    L->Draw();
    canvas->Print(out_pdf);
  }

  // (4) Per-scenario HS vs PU overlays, pT > 30 GeV only.
  canvas->SetLogy(true);
  for (auto& s : scen) {
    TLegend* L = new TLegend(0.60, 0.75, 0.92, 0.88);
    StyleLegend(L);
    L->AddEntry(s.h_hs, "Hard Scatter");
    L->AddEntry(s.h_pu, "Pile-Up");
    double ymax = 50.0 * std::max(s.h_hs->GetMaximum(), s.h_pu->GetMaximum());
    s.h_pu->SetMaximum(ymax);
    s.h_pu->SetTitle((std::string("R_{pT}: ") + s.legend).c_str());
    s.h_pu->Draw("HIST");
    s.h_hs->Draw("HIST SAME");
    drawLabels(jetcut);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // (5) All-scenario HS overlay.
  {
    TLegend* L = new TLegend(0.55, 0.58, 0.92, 0.88);
    StyleLegend(L);
    double ymax = 0;
    for (auto& s : scen) ymax = std::max(ymax, s.h_hs->GetMaximum());
    scen[0].h_hs->SetMaximum(50.0 * ymax);
    scen[0].h_hs->SetTitle("Hard Scatter R_{pT} — all scenarios");
    for (size_t i = 0; i < scen.size(); ++i) {
      scen[i].h_hs->SetLineStyle(1);
      scen[i].h_hs->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(scen[i].h_hs, scen[i].legend.c_str());
    }
    drawLabels(jetcut);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // (6) All-scenario PU overlay.
  {
    TLegend* L = new TLegend(0.55, 0.58, 0.92, 0.88);
    StyleLegend(L);
    double ymax = 0;
    for (auto& s : scen) ymax = std::max(ymax, s.h_pu->GetMaximum());
    scen[0].h_pu->SetMaximum(50.0 * ymax);
    scen[0].h_pu->SetTitle("Pile-Up R_{pT} — all scenarios");
    for (size_t i = 0; i < scen.size(); ++i) {
      scen[i].h_pu->SetLineStyle(1);
      scen[i].h_pu->Draw(i == 0 ? "HIST" : "HIST SAME");
      L->AddEntry(scen[i].h_pu, scen[i].legend.c_str());
    }
    drawLabels(jetcut);
    L->Draw();
    canvas->Print(out_pdf);
  }

  // ── Flagging vs Correcting pages ─────────────────────────────────────────────
  for (size_t i = 0; i < rocs_f.size(); ++i) {
    styleRoc(rocs_f[i],      flag[i].color,
             "R_{pT} Flagging vs Correcting;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    styleRoc(rocs_fhipt[i],  flag_hipt[i].color,
             "R_{pT} Flagging vs Correcting (p_{T} > 50 GeV);Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    rocs_f[i]->SetMaximum(300);
    rocs_fhipt[i]->SetMaximum(500);
  }
  // correcting = solid, flagging = dashed
  rocs_f[1]->SetLineStyle(1);     rocs_fhipt[1]->SetLineStyle(1);
  rocs_f[2]->SetLineStyle(2);     rocs_fhipt[2]->SetLineStyle(2);

  drawRocWithRatio(rocs_f,     flag,      300.0, 4.0, jetcut);
  drawRocWithRatio(rocs_fhipt, flag_hipt, 500.0, 5.0, jetcut_hi);

  canvas->Print(out_pdf + "]");

  // ── Console summary ────────────────────────────────────────────────────────
  std::cout << "\n=== EVENT SELECTION ===\n";
  std::cout << "Total events                   : " << n_total    << '\n';
  std::cout << "Pass basic+jet selection       : " << n_pass_sel << '\n';
  std::cout << "  with pure TRKPTZ cluster     : " << n_pure_evt  << '\n';
  std::cout << "  with clean HS timing (>=95%) : " << n_misas_evt << '\n';

  std::cout << "\n=== ENTRIES PER SCENARIO ===\n";
  for (auto& s : scen) {
    std::cout << "  " << s.name
              << "  HS: " << (long)s.h_hs->Integral()
              << "  PU: " << (long)s.h_pu->Integral() << '\n';
  }

  std::cout << "\nWrote " << out_pdf << "\n";
  return 0;
}
