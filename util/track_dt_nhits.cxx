// track_dt_nhits.cxx — per-track timing residual split by number of HGTD hits.
//
// track_dt.cxx's quantity, (HGTD track time - truth particle prodVtx time), for
// every track with a valid HGTD time and a linked truth particle, binned by
// Track_nHGTDHits (1, 2, 3, >= 4). No event selection, no track cuts.
//
// Also split by Track_nHGTDPrimaryHits == 0 against >= 1: a track none of whose
// hits come from its own particle carries somebody else's time.
//
// --by=pt bins the same quantities in track pT instead (< 1, 1-2, 2-5, 5-10,
// 10-30, > 30 GeV) and writes track_dt_pt.pdf.
//
// Output: <OUTPUT_DIR>/track_dt_nhits.pdf or track_dt_pt.pdf (4 pages: residual
// log and linear, pull, the no-primary-hit split) and a console table.
//   ./track_dt_nhits [--by=nhits|pt] [--sample=<name>] [--max-events=N]

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
  long n = 0, n30 = 0, n60 = 0, n90 = 0, n3s = 0, nNoPrim = 0, nNoPrimBad = 0, nPrimBad = 0, nOneHit = 0;
  double sumRes = 0, sumDt = 0;
};

// Core width: Gaussian + flat floor over +-120 ps. The tails are flat out to +-400 ps (mis-matched hits),
// not Gaussian, so a double-Gaussian tail width is ill-defined here and its fit unstable from bin to bin;
// the tail is reported as counted fractions instead.
std::array<double,3> fitDG(TH1D* h, const char* name) {
  TF1* f = new TF1(name, "[1]*TMath::Exp(-0.5*((x-[0])/[2])^2)+[3]", -120.0, 120.0);
  f->SetParameters(h->GetBinCenter(h->GetMaximumBin()), 0.9 * h->GetMaximum(), 28.0, 0.02 * h->GetMaximum());
  // Limits scaled to the histogram: Minuit's bounded-parameter transform loses all precision when a
  // value of ~1e3 sits in [0, 1e12], which is what broke the fit in the smallest bin.
  const double mx = h->GetMaximum();
  f->SetParLimits(0, -50.0, 50.0); f->SetParLimits(1, 0, 10 * mx);
  f->SetParLimits(2, 5.0, 100.0);  f->SetParLimits(3, 0, mx);
  h->Fit(f, "RQ0L");   // likelihood: the > 30 GeV bin has ~28k tracks and a chi2 fit lands in a false minimum
  return {f->GetParameter(2), f->GetParameter(0), f->GetParError(2)};
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

  bool byPt = false;
  for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--by=pt") byPt = true;
  const int NB = byPt ? 6 : 4;
  const std::vector<const char*> lab = byPt
    ? std::vector<const char*>{"p_{T} < 1 GeV", "1 #minus 2 GeV", "2 #minus 5 GeV", "5 #minus 10 GeV", "10 #minus 30 GeV", "> 30 GeV"}
    : std::vector<const char*>{"1 hit", "2 hits", "3 hits", "#geq 4 hits"};
  const std::vector<const char*> row = byPt
    ? std::vector<const char*>{"<1", "1-2", "2-5", "5-10", "10-30", ">30"}
    : std::vector<const char*>{"1", "2", "3", ">=4"};
  const std::vector<Color_t> col = byPt ? std::vector<Color_t>{C02, C07, C03, C08, C01, C04}
                                        : std::vector<Color_t>{C02, C01, C08, C04};
  const double ptEdge[5] = {1.0, 2.0, 5.0, 10.0, 30.0};
  // the two bins drawn on the own-particle page
  const int splitA = 0, splitB = byPt ? 3 : 1;
  std::vector<Acc> a(NB);
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
      int bin = std::min(nh, 4) - 1;
      if (byPt) { const double pt = branch.trackPt[idx]; bin = 0; while (bin < 5 && pt >= ptEdge[bin]) ++bin; }
      Acc& x = a[bin];
      if (nh == 1) ++x.nOneHit;
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
  std::printf("\n=== HGTD track time - truth particle time, by %s (%lld events, %ld tracks) ===\n",
              byPt ? "track pT [GeV]" : "number of HGTD hits", (long long)nEv, nAll);
  std::printf("%-9s %10s %7s %8s %9s %9s %9s %7s %7s %7s %8s | %12s %14s %14s %9s\n", byPt ? "pT" : "n hits", "tracks", "share", "<sig_t>", "core sig", "peak", "sig/quoted",
              "<30ps", "<60ps", "<90ps", "<3 sig_t", "no prim. hit", "bad | no prim", "bad | >=1 prim", "1-hit");
  for (int i = 0; i < NB; ++i) {
    Acc& x = a[i];
    if (x.n == 0) continue;
    auto f = fitDG(x.dt, Form("fit_%d", i));
    std::printf("%-9s %10ld %6.2f%% %7.1f %8.1f %9.1f %9.2f %6.1f%% %6.1f%% %6.1f%% %7.1f%% | %11.1f%% %13.1f%% %13.1f%% %8.1f%%\n",
                row[i], x.n, 100.0 * x.n / nAll, x.sumRes / x.n, f[0], f[1], f[0] / (x.sumRes / x.n),
                100.0 * x.n30 / x.n, 100.0 * x.n60 / x.n, 100.0 * x.n90 / x.n, 100.0 * x.n3s / x.n,
                100.0 * x.nNoPrim / x.n, x.nNoPrim ? 100.0 * x.nNoPrimBad / x.nNoPrim : 0.0,
                (x.n - x.nNoPrim) ? 100.0 * x.nPrimBad / (x.n - x.nNoPrim) : 0.0, 100.0 * x.nOneHit / x.n);
  }
  std::printf("  core sig / peak: Gaussian + flat floor fitted over +-120 ps.  'bad' = |dt| >= 3 sigma_t.\n"
              "  no prim. hit: Track_nHGTDPrimaryHits == 0, i.e. none of the track's hits come from its own particle.\n");

  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  const std::string out = MyUtl::plotFilePath("", byPt ? "track_dt_pt.pdf" : "track_dt_nhits.pdf");
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
    TLegend* leg = new TLegend(0.64, byPt ? 0.62 : 0.70, 0.92, 0.90); StyleLegend(leg);
    for (int i = 0; i < NB; ++i) {
      hs[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
      leg->AddEntry(hs[i], Form(100.0 * a[i].n / nAll < 1 ? "%s (%.2f%%)" : "%s (%.0f%%)", lab[i], 100.0 * a[i].n / nAll), "l");
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
  // page 4: two of the bins, split by whether any hit is from the track's own particle. Unit area per
  // BIN (both of its curves share one normalisation), so bins of very different size can be compared.
  c->SetLogy(true);
  {
    const int bb[2] = {splitA, splitB};
    std::vector<TH1D*> h; double mx = 0;
    TLegend* leg = new TLegend(0.50, 0.70, 0.92, 0.90); StyleLegend(leg);
    for (int k = 0; k < 2; ++k) {
      const Acc& x = a[bb[k]];
      for (int np = 0; np < 2; ++np) {
        TH1D* g = (TH1D*)(np ? x.dtNoPrim : x.dtPrim)->Clone(Form("sp_%d_%d", k, np));
        if (x.n > 0) g->Scale(1.0 / x.n);
        g->SetLineColor(col[bb[k]]); g->SetLineWidth(2); g->SetLineStyle(np ? 2 : 1);
        g->GetYaxis()->SetTitle("Fraction of the bin's tracks / 4 ps");
        mx = std::max(mx, g->GetMaximum()); h.push_back(g);
        leg->AddEntry(g, Form("%s, %s", lab[bb[k]], np ? "no hit from own particle" : "#geq 1 hit from own particle"), "l");
      }
    }
    h[0]->SetMaximum(30 * mx); h[0]->SetMinimum(2e-6);
    for (size_t i = 0; i < h.size(); ++i) h[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
    leg->Draw();
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
    c->Print(out.c_str());
  }
  c->Print((out + "]").c_str());
  std::cout << "Wrote " << out << "\n";
  return 0;
}
