// rpt_v5.cxx — Jet-level RpT analysis using the WAVeS algorithm.
//
// Clean comparison of five scenarios (no per-case diagnostics — see rpt_v4 for those).
// The "mine" time now comes from the WAVeS-selected cluster (Score::WAVES) with the
// in-jet timing refinement, replacing the TRKPTZ selection used in rpt_v4.
//
// Scenarios:
//   1. zonly        — z-significance tracks only, no time gate   ("ITk-only")
//   2. hgtd         — ntuple RecoVtx_time / RecoVtx_timeRes gate  ("HGTD t_{0}")
//   3. waves        — WAVeS-selected cluster time gate            ("WAVeS t_{0}")
//   4. waves_misas  — WAVeS time gate, events gated on HS timing  ("WAVeS t_{0} + clean timing")
//                     purity ≥ MISAS_PURITY_CUT (event filter)
//   5. truth        — reco track times gated vs truth vertex t_{0} ("Truth t_{0}")
//
// No event-level selection besides the MISAS filter on scenario 4 — every forward
// jet in the acceptance contributes an independent RpT measurement.
//
// Two jet pT windows:
//   Slice A: 30 < pT < 40 GeV
//   Slice B: pT > 40 GeV
// Jet eta acceptance: 2.4 < |eta| < 3.8.
//
// Output: figs/rpt_plots/rpt_v5.pdf

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <boost/filesystem.hpp>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
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

// Jet acceptance for this script (paper values, distinct from rpt_v2).
static constexpr double JET_ETA_MIN = 2.4;
static constexpr double JET_ETA_MAX = 3.8;

// Event-level HS-timing-purity threshold for the "clean timing" WAVeS scenario:
// a jet's RpT is filled only in events where ≥ this fraction of HS pT is timed
// within |pull| < 3σ (calcHSTimingPurity).  Mirrors the Score::WAVES_MISAS
// oracle.
static constexpr float MISAS_PURITY_CUT = 0.75f;

// Per-track time-gate half-width in σ.  A 2σ cut over-trims genuine HS tracks
// when the vertex time is slightly mis-estimated, dragging the high-efficiency
// end of the ROC below ITk-only (worst in the >40 GeV slice).  Loosening it
// recovers that region at a small low-efficiency cost.
static constexpr double GATE_SIGMA = 2.5;

// -----------------------------------------------------------------------------
// ROC: HS efficiency vs PU rejection (1 / mistag).
// -----------------------------------------------------------------------------
TGraph* generate_roc(TH1D* PU_hist, TH1D* HS_hist) {
  int bin = PU_hist->GetNbinsX();
  double HS_tot = HS_hist->Integral();
  double PU_tot = PU_hist->Integral();
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

// -----------------------------------------------------------------------------
// Scenario container.
// -----------------------------------------------------------------------------
struct Scenario {
  std::string name;
  std::string legend;
  Color_t color;
  TH1D* h_hs;
  TH1D* h_pu;
};

static inline double dR(double j_eta, double j_phi, double t_eta, double t_phi) {
  double deta = j_eta - t_eta;
  double dphi = TVector2::Phi_mpi_pi(j_phi - t_phi);
  return std::sqrt(deta * deta + dphi * dphi);
}

// RpT using ghost association (paper definition): sum pT of tracks that are
// both ghost-associated to the jet AND in the caller's z₀-selected set.
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
  setupChain(chain, "../../highstats-ntuple/");
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found.  Aborting.\n";
    return 1;
  }

  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  // Paper jet-label helper (ATL-HGTD-PUB-2022-001 Sec. 3).
  // HS  : dR(reco, truthHS) < 0.3  AND  truthHS_pT > 10 GeV
  auto paperIsHS = [&](double j_eta, double j_phi) {
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      if (branch.truthHSJetPt[t] < 10.0) continue;
      if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.3)
        return true;
    }
    return false;
  };
  // PU: dR > 0.6 from any truth HS jet with pT > 4 GeV (Sec. 3).
  auto paperIsPU = [&](double j_eta, double j_phi) {
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      if (branch.truthHSJetPt[t] < 4.0) continue;
      if (dR(j_eta, j_phi, branch.truthHSJetEta[t], branch.truthHSJetPhi[t]) < 0.6)
        return false;
    }
    return true;
  };

  // Non-uniform binning: fine bins in [0, 2.5] for ROC granularity, then one
  // wide bin to capture the tail without bloating memory.  250 bins of width
  // 0.01 — coarse enough that the ROC + error bars aren't overcrowded.
  static std::vector<double> rpt_bins = []() {
    std::vector<double> b;
    b.reserve(252);
    for (int i = 0; i <= 250; ++i) b.push_back(0.01 * i);
    b.push_back(375.0);
    return b;
  }();
  const int nbin = (int)rpt_bins.size() - 1;  // 251

  auto makeHist = [&](const char* name, const char* title) {
    return new TH1D(name, title, nbin, rpt_bins.data());
  };

  auto makeScenarios = [&](const std::string& suffix) {
    std::vector<Scenario> s = {
      {"zonly",       "ITk-only",                    C05, nullptr, nullptr},
      {"hgtd",        "HGTD t_{0}",                  C01, nullptr, nullptr},
      {"waves",       "WAVeS t_{0}",                 C03, nullptr, nullptr},
      {"waves_misas", "WAVeS t_{0} + clean timing",  C04, nullptr, nullptr},
      {"truth",       "Truth t_{0}",                 C06, nullptr, nullptr},
    };
    for (auto& sc : s) {
      sc.h_hs = makeHist(("HS_" + sc.name + suffix).c_str(),
                         ("Hard Scatter R_{pT}: " + sc.legend + ";R_{pT};Entries").c_str());
      sc.h_pu = makeHist(("PU_" + sc.name + suffix).c_str(),
                         ("Pile-Up R_{pT}: "      + sc.legend + ";R_{pT};Entries").c_str());
    }
    return s;
  };

  std::vector<Scenario> scen_lo = makeScenarios("_lo");  // 30–40 GeV
  std::vector<Scenario> scen_hi = makeScenarios("_hi");  // >40 GeV

  long n_total = 0, n_pass_basic = 0, n_hgtd_valid = 0;

  // Untimed-track-floor diagnostic (>40 GeV slice): of the ghost-associated,
  // z-selected track pT in each jet, how much is carried by tracks with NO valid
  // HGTD time?  Those tracks pass every timing gate unconditionally, so they set
  // the irreducible RpT floor that no timing method (even Truth) can remove.
  double pu_tot_pt = 0, pu_floor_pt = 0, hs_tot_pt = 0, hs_floor_pt = 0;        // >40
  double pu_tot_lo = 0, pu_floor_lo = 0, hs_tot_lo = 0, hs_floor_lo = 0;        // 30-40

  while (reader.Next()) {
    ++n_total;

    // ── Require only vertex quality (paper Sec. 3: |z_reco − z_truth| < 2 mm).
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    if (std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]) > MAX_VTX_DZ) continue;
    ++n_pass_basic;

    // ── Track selection: all tracks by z-significance (no eta cut) for the
    //    z-only baseline, matching the paper's ITk-only scenario. ────────────
    std::vector<int> trk_all;
    for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
      double trkPt = branch.trackPt[trk];
      if (trkPt < MIN_TRACK_PT || trkPt > MAX_TRACK_PT) continue;
      if (!branch.trackQuality[trk]) continue;
      if (passTrackVertexAssociation((int)trk, 0, &branch, 2.5))
        trk_all.push_back((int)trk);
    }

    // ── HGTD-acceptance tracks only (used for WAVeS clustering). ────────────
    std::vector<int> trk_z = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 2.5);

    // ── WAVeS clustering + selection. ────────────────────────────────────────
    auto clusters = clusterTracksInTime(
        trk_z, &branch, DIST_CUT_CONE,
        /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
        ClusteringMethod::ITERATIVE, /*useZ0=*/false,
        /*sortTracks=*/false, /*calcPurityFlag=*/true);

    // WAVeS selection: highest WAVeS-score cluster; time via in-jet refinement.
    double t_waves = 0.0, var_waves = 0.0;
    bool   waves_ok = false;
    if (!clusters.empty()) {
      auto best   = chooseCluster(clusters, Score::WAVES);
      t_waves     = best.calculateTime(Score::WAVES, &branch);  // in-jet refined
      var_waves   = best.sigmas[0] * best.sigmas[0];
      waves_ok    = true;
    }

    // ── HGTD ntuple vertex time. ─────────────────────────────────────────────
    double t_hgtd         = branch.recoVtxTime[0];
    double var_hgtd       = branch.recoVtxTimeRes[0] * branch.recoVtxTimeRes[0];
    bool   hgtd_vtx_valid = (branch.recoVtxValid[0] == 1);
    if (hgtd_vtx_valid) ++n_hgtd_valid;

    // ── Per-track time gate.  pull width ~1.5 → var_vtx ×2.25. ──────────────
    auto applyTimeGate = [&](const std::vector<int>& base,
                              double t_vtx, double var_vtx, bool vtx_valid,
                              double sigma = 2.0) {
      std::vector<int> out;
      out.reserve(base.size());
      for (int idx : base) {
        bool apply = vtx_valid && branch.trackTimeValid[idx] == 1;
        if (!apply) { out.push_back(idx); continue; }
        double dt    = branch.trackTime[idx] - t_vtx;
        double var_t = branch.trackTimeRes[idx] * branch.trackTimeRes[idx];
        double pull  = std::abs(dt) / std::sqrt(2.25 * var_vtx + var_t);
        if (pull < sigma) out.push_back(idx);
      }
      return out;
    };

    std::vector<int> trk_hgtd  = applyTimeGate(trk_all, t_hgtd,  var_hgtd,  hgtd_vtx_valid, GATE_SIGMA);
    std::vector<int> trk_waves = applyTimeGate(trk_all, t_waves, var_waves, waves_ok,       GATE_SIGMA);

    // ── Truth-vertex-t₀: gate the reco track times against the perfect HS
    //    vertex time (var_vtx = 0). ──────────────────────────────────────────
    double t_truth = branch.truthVtxTime[0];
    std::vector<int> trk_truth = applyTimeGate(trk_all, t_truth, 0.0, true, GATE_SIGMA);

    // ── "Clean timing" event filter: only fill the waves_misas scenario in
    //    events whose HS timing purity clears MISAS_PURITY_CUT.  Same WAVeS time
    //    gate otherwise. ───────────────────────────────────────────────────────
    bool gate_misas = (calcHSTimingPurity(trk_z, &branch) >= MISAS_PURITY_CUT);

    // Build per-scenario sets once per event for O(1) ghost-index lookup.
    std::unordered_set<int> set_all  (trk_all.begin(),   trk_all.end());
    std::unordered_set<int> set_hgtd (trk_hgtd.begin(),  trk_hgtd.end());
    std::unordered_set<int> set_waves(trk_waves.begin(), trk_waves.end());
    std::unordered_set<int> set_truth(trk_truth.begin(), trk_truth.end());

    // ── Fill jets into pT slices. ─────────────────────────────────────────────
    auto fillJets = [&](std::vector<Scenario>& sv, double pt_lo, double pt_hi) {
      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
        double j_pt  = branch.topoJetPt[j];
        double j_eta = branch.topoJetEta[j];
        double j_phi = branch.topoJetPhi[j];
        if (j_pt <= pt_lo || j_pt >= pt_hi) continue;
        if (std::abs(j_eta) < JET_ETA_MIN || std::abs(j_eta) > JET_ETA_MAX) continue;
        bool isHS = paperIsHS(j_eta, j_phi);
        bool isPU = paperIsPU(j_eta, j_phi);
        if (!isHS && !isPU) continue;
        const auto& ghost = branch.topoJetGhostTrackIdx[j];
        auto fill = [&](Scenario& s, const std::unordered_set<int>& s_set, bool ok = true) {
          if (!ok) return;
          double r = computeRpT(&branch, ghost, j_pt, s_set);
          if (isHS) s.h_hs->Fill(r);
          else      s.h_pu->Fill(r);
        };
        fill(sv[0], set_all);                       // ITk-only
        fill(sv[1], set_hgtd);                      // HGTD t0
        fill(sv[2], set_waves);                     // WAVeS t0
        fill(sv[3], set_waves, gate_misas);         // WAVeS t0 + clean timing (event filter)
        fill(sv[4], set_truth);                     // Truth t0

        // Untimed floor accounting, per slice.
        for (int idx : ghost) {
          if (!set_all.count(idx)) continue;
          double pt = branch.trackPt[idx];
          bool untimed = (branch.trackTimeValid[idx] != 1);
          if (pt_lo >= 40.0) {
            if (isHS) { hs_tot_pt += pt; if (untimed) hs_floor_pt += pt; }
            else      { pu_tot_pt += pt; if (untimed) pu_floor_pt += pt; }
          } else {
            if (isHS) { hs_tot_lo += pt; if (untimed) hs_floor_lo += pt; }
            else      { pu_tot_lo += pt; if (untimed) pu_floor_lo += pt; }
          }
        }
      }
    };

    fillJets(scen_lo, 30.0, 40.0);
    fillJets(scen_hi, 40.0, 1e9);
  }

  // ── Output paths ─────────────────────────────────────────────────────────────
  boost::filesystem::create_directories("../figs/rpt_plots");
  const TString out_pdf = "../figs/rpt_plots/rpt_v5.pdf";

  TCanvas* canvas = new TCanvas("canvas", "RpT v5", 800, 700);
  canvas->Print(out_pdf + "[");

  auto drawLabels = [](const char* extra = nullptr) {
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82);
    if (extra) {
      TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
      t.DrawLatex(0.18, 0.76, extra);
    }
  };
  const char* lbl_lo = "30 < p_{T}^{jet} < 40 GeV, 2.4 < |#eta| < 3.8";
  const char* lbl_hi = "p_{T}^{jet} > 40 GeV, 2.4 < |#eta| < 3.8";

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

  const double roc_xmin = 0.8, roc_xmax = 1.0;
  std::vector<TGraph*> rocs_lo, rocs_hi;
  for (auto& s : scen_lo) rocs_lo.push_back(generate_roc(s.h_pu, s.h_hs));
  for (auto& s : scen_hi) rocs_hi.push_back(generate_roc(s.h_pu, s.h_hs));

  auto styleRoc = [&](TGraph* g, Color_t col) {
    g->SetTitle("R_{pT} Discriminant;Hard Scatter Efficiency;Pile-Up Rejection (1 / Mistag Rate)");
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->SetLineWidth(2);
    g->GetXaxis()->SetLimits(roc_xmin, roc_xmax);
    g->GetXaxis()->SetNdivisions(810);
    g->SetMinimum(1.0);
  };
  for (size_t i = 0; i < rocs_lo.size(); ++i) {
    styleRoc(rocs_lo[i], scen_lo[i].color);
    styleRoc(rocs_hi[i], scen_hi[i].color);
  }

  // ROC page with ratio-to-zonly panel.
  auto drawRocWithRatio = [&](std::vector<TGraph*>& gs,
                               std::vector<Scenario>& sc,
                               double ymax, double ratio_ymax,
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
      t.SetTextSize(0.044); t.DrawLatex(0.18, 0.80, "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.");
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
  drawRocWithRatio(rocs_lo, scen_lo, 300.0, 4.0, lbl_lo);

  // (2) ROC — >40 GeV.
  drawRocWithRatio(rocs_hi, scen_hi, 300.0, 4.0, lbl_hi);

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

  // ── Console summary ───────────────────────────────────────────────────────────
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
    double hsTot = HS->Integral(), puTot = PU->Integral();
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

  std::cout << "\nWrote " << out_pdf << "\n";
  return 0;
}
