// track_dt_nhits.cxx — per-track timing residual split by number of HGTD hits.
//
// track_dt.cxx's quantity, (HGTD track time - truth particle prodVtx time), for
// every track with a valid HGTD time and a linked truth particle, binned by
// Track_nHGTDHits (1, 2, 3, >= 4). No event selection, no track cuts.
//
// Also split by Track_nHGTDPrimaryHits == 0 against >= 1: a track none of whose
// hits come from its own particle carries somebody else's time.
//
// Output: <OUTPUT_DIR>/track_dt_nhits.pdf (3 pages: residual, pull, the
// no-primary-hit split) and a console table.
//   ./track_dt_nhits [--sample=<name>] [--max-events=N]

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>

#include <boost/filesystem.hpp>
#include <iostream>

#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "sample_config.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

namespace {
struct Acc {
  TH1D *dt = nullptr, *pull = nullptr, *dtNoPrim = nullptr, *dtPrim = nullptr;
  long n = 0, n30 = 0, n60 = 0, n90 = 0, n3s = 0, nNoPrim = 0, nNoPrimBad = 0, nPrimBad = 0;
  double sumRes = 0, sumDt = 0;
};

// shared-mean double Gaussian, as track_dt.cxx; returns (core sigma, tail sigma, core fraction of the fit area)
std::array<double,3> fitDG(TH1D* h, const char* name) {
  TF1* f = new TF1(name, "[1]*TMath::Exp(-0.5*((x-[0])/[3])^2)+[2]*TMath::Exp(-0.5*((x-[0])/[4])^2)",
                   h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  f->SetParameters(h->GetMean(), 0.8 * h->GetMaximum(), 1e-2 * h->GetMaximum(), 25.0, h->GetStdDev());
  f->SetParLimits(0, -50.0, 50.0); f->SetParLimits(1, 0, 1e12); f->SetParLimits(2, 0, 1e12);
  f->SetParLimits(3, 5.0, 1e6);    f->SetParLimits(4, 5.0, 1e6);
  f->SetNpx(1000);
  h->Fit(f, "RQ0");
  double s1 = f->GetParameter(3), s2 = f->GetParameter(4), n1 = f->GetParameter(1), n2 = f->GetParameter(2);
  if (s1 > s2) { std::swap(s1, s2); std::swap(n1, n2); }
  const double a1 = n1 * s1, a2 = n2 * s2;
  return {s1, s2, (a1 + a2 > 0) ? a1 / (a1 + a2) : 0.0};
}
}  // namespace

int main(int argc, char** argv) {
  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  const auto cfg = MyUtl::resolveSample(argc, argv);
  MyUtl::SAMPLE_NAME  = cfg.sampleName;
  MyUtl::OUTPUT_DIR   = cfg.outputDir;
  MyUtl::ENERGY_LABEL = cfg.energyLabel;
  const Long64_t maxEvents = MyUtl::resolveMaxEvents(argc, argv);

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  const int NB = 4;
  const char* lab[NB] = {"1 hit", "2 hits", "3 hits", "#geq 4 hits"};
  const Color_t col[NB] = {C02, C01, C08, C04};
  Acc a[NB];
  for (int i = 0; i < NB; ++i) {
    a[i].dt       = new TH1D(Form("dt_%d", i), ";t_{track}^{HGTD} #minus t_{truth particle} [ps];Fraction of tracks / 4 ps", 200, -400, 400);
    a[i].pull     = new TH1D(Form("pull_%d", i), ";(t_{track}^{HGTD} #minus t_{truth particle}) / #sigma_{t};Fraction of tracks", 200, -10, 10);
    a[i].dtNoPrim = new TH1D(Form("dtnp_%d", i), ";t_{track}^{HGTD} #minus t_{truth particle} [ps];Tracks / 4 ps", 200, -400, 400);
    a[i].dtPrim   = new TH1D(Form("dtp_%d", i), ";t_{track}^{HGTD} #minus t_{truth particle} [ps];Tracks / 4 ps", 200, -400, 400);
  }

  Long64_t nEv = 0;
  while (reader.Next()) {
    if (maxEvents > 0 && nEv >= maxEvents) break;
    ++nEv;
    for (size_t idx = 0; idx < branch.trackTime.GetSize(); ++idx) {
      if (branch.trackTimeValid[idx] != 1) continue;
      const int p = branch.trackToParticle[idx];
      if (p < 0 || p >= (int)branch.particleT.GetSize()) continue;
      const int nh = branch.trackHgtdHits[idx];
      if (nh < 1) continue;
      Acc& x = a[std::min(nh, NB) - 1];
      const double dt = branch.trackTime[idx] - branch.particleT[p];
      const double res = branch.trackTimeRes[idx];
      const bool noPrim = (branch.trackPrimHits[idx] == 0);
      const bool bad = (res > 0) && (std::abs(dt) >= 3.0 * res);
      x.dt->Fill(dt);
      if (res > 0) x.pull->Fill(dt / res);
      (noPrim ? x.dtNoPrim : x.dtPrim)->Fill(dt);
      ++x.n; x.sumRes += res; x.sumDt += dt;
      if (std::abs(dt) < 30) ++x.n30;
      if (std::abs(dt) < 60) ++x.n60;
      if (std::abs(dt) < 90) ++x.n90;
      if (!bad) ++x.n3s;
      if (noPrim) { ++x.nNoPrim; if (bad) ++x.nNoPrimBad; }
      else if (bad) ++x.nPrimBad;
    }
  }

  long nAll = 0;
  for (auto& x : a) nAll += x.n;
  std::printf("\n=== HGTD track time - truth particle time, by number of HGTD hits (%lld events, %ld tracks) ===\n", (long long)nEv, nAll);
  std::printf("%-9s %10s %7s %8s %9s %9s %9s %7s %7s %7s %8s | %12s %14s %14s\n", "n hits", "tracks", "share", "<sig_t>", "core sig", "tail sig", "core frac",
              "<30ps", "<60ps", "<90ps", "<3 sig_t", "no prim. hit", "bad | no prim", "bad | >=1 prim");
  for (int i = 0; i < NB; ++i) {
    Acc& x = a[i];
    if (x.n == 0) continue;
    auto f = fitDG(x.dt, Form("fit_%d", i));
    std::printf("%-9s %10ld %6.1f%% %7.1f %8.1f %9.1f %8.1f%% %6.1f%% %6.1f%% %6.1f%% %7.1f%% | %11.1f%% %13.1f%% %13.1f%%\n",
                i == 3 ? ">=4" : Form("%d", i + 1), x.n, 100.0 * x.n / nAll, x.sumRes / x.n, f[0], f[1], 100 * f[2],
                100.0 * x.n30 / x.n, 100.0 * x.n60 / x.n, 100.0 * x.n90 / x.n, 100.0 * x.n3s / x.n,
                100.0 * x.nNoPrim / x.n, x.nNoPrim ? 100.0 * x.nNoPrimBad / x.nNoPrim : 0.0,
                (x.n - x.nNoPrim) ? 100.0 * x.nPrimBad / (x.n - x.nNoPrim) : 0.0);
  }
  std::printf("  core/tail sig: shared-mean double Gaussian over +-400 ps.  'bad' = |dt| >= 3 sigma_t.\n"
              "  no prim. hit: Track_nHGTDPrimaryHits == 0, i.e. none of the track's hits come from its own particle.\n");

  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  const std::string out = MyUtl::plotFilePath("", "track_dt_nhits.pdf");
  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Print((out + "[").c_str());
  auto drawSet = [&](bool pull, bool logy) {
    c->SetLogy(logy);
    double mx = 0;
    std::vector<TH1D*> hs;
    for (int i = 0; i < NB; ++i) {
      TH1D* h = (TH1D*)(pull ? a[i].pull : a[i].dt)->Clone(Form("n_%d_%d_%d", i, pull, logy));
      if (h->Integral(0, -1) > 0) h->Scale(1.0 / h->Integral(0, -1));
      h->SetLineColor(col[i]); h->SetLineWidth(2); h->SetMarkerSize(0);
      mx = std::max(mx, h->GetMaximum()); hs.push_back(h);
    }
    hs[0]->SetMaximum(logy ? 30 * mx : 1.35 * mx);
    if (logy) hs[0]->SetMinimum(2e-5);
    TLegend* leg = new TLegend(0.66, 0.70, 0.92, 0.90); StyleLegend(leg);
    for (int i = 0; i < NB; ++i) {
      hs[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
      leg->AddEntry(hs[i], Form("%s (%.0f%%)", lab[i], 100.0 * a[i].n / nAll), "l");
    }
    leg->Draw();
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.032);
    tl.DrawLatex(0.18, 0.76, "Tracks with an HGTD time and a truth particle");
    tl.DrawLatex(0.18, 0.71, "Unit area; legend: share of tracks");
    c->Print(out.c_str());
  };
  drawSet(false, true);
  drawSet(false, false);
  drawSet(true, true);
  // page 4: 1-hit and 2-hit tracks, split by whether any hit is from the track's own particle
  c->SetLogy(true);
  {
    TH1D* h[4] = {a[0].dtPrim, a[0].dtNoPrim, a[1].dtPrim, a[1].dtNoPrim};
    const Color_t cc[4] = {C02, C02, C01, C01};
    const char* ll[4] = {"1 hit, from own particle", "1 hit, not from own particle", "2 hits, #geq 1 from own particle", "2 hits, none from own particle"};
    double mx = 0; for (auto* x : h) mx = std::max(mx, x->GetMaximum());
    TLegend* leg = new TLegend(0.56, 0.70, 0.92, 0.90); StyleLegend(leg);
    for (int i = 0; i < 4; ++i) {
      h[i]->SetLineColor(cc[i]); h[i]->SetLineWidth(2); h[i]->SetLineStyle(i % 2 ? 2 : 1);
      if (i == 0) { h[i]->SetMaximum(30 * mx); h[i]->SetMinimum(0.5); }
      h[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
      leg->AddEntry(h[i], ll[i], "l");
    }
    leg->Draw();
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
    c->Print(out.c_str());
  }
  c->Print((out + "]").c_str());
  std::cout << "Wrote " << out << "\n";
  return 0;
}
