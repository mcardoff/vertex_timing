// track_dt_nhits.cxx — per-track timing residual, binned three ways in one pass.
//
// track_dt.cxx's quantity, (HGTD track time - truth particle prodVtx time), for
// every track with a valid HGTD time and a linked truth particle. No event
// selection, no track cuts. Binned by
//   number of HGTD hits   1, 2, 3, >= 4
//   track pT              < 1, 1-2, 2-5, 5-10, 10-30, > 30 GeV
//   track |eta|           2.4-2.6, 2.6-2.8, 2.8-3.0, 3.0-3.2, 3.2-3.5, 3.5-4.0
// and for each binning: residual (log and linear), pull, and the split by
// whether any of the track's hits comes from its own particle
// (Track_nHGTDPrimaryHits == 0 against >= 1) for two of the bins.
//
// Core width: Gaussian + flat floor fitted over +-120 ps. The tails are flat
// out to +-400 ps (mis-matched hits), not Gaussian, so a double-Gaussian tail
// is ill-defined here; the tail is reported as counted fractions.
//
// Mass correction: the reconstruction's time-of-flight correction takes the
// particle at the speed of light. A particle of mass m and momentum p arrives
// later by (L/c)(1/beta - 1), L = z_HGTD / |cos theta|. The "corrected"
// residual subtracts that, using the TRUTH mass (TruthPart_m) and the track's
// momentum pT cosh(eta): a test of the hypothesis, not a reco-level correction.
// Each binning gets a corrected page, the tables get before/after columns, and
// a species split (pion / kaon / proton / e / mu / other) is added at the end.
//
// Output: <OUTPUT_DIR>/track_dt_binned.pdf (20 pages: per binning residual log, linear, pull, mass-corrected, hit provenance) and the console tables.
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
  TH1D *dt = nullptr, *pull = nullptr, *dtNoPrim = nullptr, *dtPrim = nullptr, *dtc = nullptr;
  long n = 0, n30 = 0, n60 = 0, n90 = 0, n3s = 0, nNoPrim = 0, nNoPrimBad = 0, nPrimBad = 0, nOneHit = 0;
  long n60c = 0, n3sc = 0;
  double sumRes = 0, sumDt = 0;
};
constexpr double Z_HGTD = 3500.0;   // mm, nominal front face
constexpr double C_MM_PS = 0.299792458;  // mm / ps

// pure Gaussian over +-90 ps -> (sigma, mean, chi2/ndf); the fit-quality comparison before/after
std::array<double,3> fitGaus(TH1D* h, const char* name) {
  TF1* f = new TF1(name, "gaus", -90.0, 90.0);
  f->SetParameters(h->GetMaximum(), h->GetBinCenter(h->GetMaximumBin()), 28.0);
  h->Fit(f, "RQ0");
  return {f->GetParameter(2), f->GetParameter(1), f->GetNDF() > 0 ? f->GetChisquare() / f->GetNDF() : 0.0};
}

struct Binning {
  const char* tag;                 // histogram-name stem
  const char* title;               // console / page heading
  std::vector<const char*> lab;    // legend
  std::vector<const char*> row;    // console
  std::vector<Color_t> col;
  int splitA, splitB;              // the two bins drawn on the own-particle page
  std::vector<Acc> a;
};

// Gaussian + flat floor over +-120 ps -> (core sigma, peak, error on sigma)
std::array<double,3> fitCore(TH1D* h, const char* name) {
  TF1* f = new TF1(name, "[1]*TMath::Exp(-0.5*((x-[0])/[2])^2)+[3]", -120.0, 120.0);
  const double mx = h->GetMaximum();
  f->SetParameters(h->GetBinCenter(h->GetMaximumBin()), 0.9 * mx, 28.0, 0.02 * mx);
  // Limits scaled to the histogram: Minuit's bounded-parameter transform loses
  // its precision when a value of ~1e3 sits in [0, 1e12].
  f->SetParLimits(0, -50.0, 50.0); f->SetParLimits(1, 0, 10 * mx);
  f->SetParLimits(2, 5.0, 100.0);  f->SetParLimits(3, 0, mx);
  h->Fit(f, "RQ0L");
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
  // Bound here rather than in BranchPointerWrapper: the grid skims do not carry
  // these, and a wrapper branch missing from a file makes TTreeReader iterate
  // zero entries for every executable.
  TTreeReaderArray<float> particleM(reader, "TruthPart_m");
  TTreeReaderArray<float> particlePdgId(reader, "TruthPart_pdgId");   // stored as float

  std::vector<Binning> B = {
    {"nhits", "number of HGTD hits",
     {"1 hit", "2 hits", "3 hits", "#geq 4 hits"}, {"1", "2", "3", ">=4"}, {C02, C01, C08, C04}, 0, 1, {}},
    {"pt", "track p_{T}",
     {"p_{T} < 1 GeV", "1 #minus 2 GeV", "2 #minus 5 GeV", "5 #minus 10 GeV", "10 #minus 30 GeV", "> 30 GeV"},
     {"<1", "1-2", "2-5", "5-10", "10-30", ">30"}, {C02, C07, C03, C08, C01, C04}, 0, 3, {}},
    {"eta", "track |#eta|",
     {"|#eta| 2.4 #minus 2.6", "2.6 #minus 2.8", "2.8 #minus 3.0", "3.0 #minus 3.2", "3.2 #minus 3.5", "3.5 #minus 4.0"},
     {"2.4-2.6", "2.6-2.8", "2.8-3.0", "3.0-3.2", "3.2-3.5", "3.5-4.0"}, {C02, C07, C03, C08, C01, C04}, 0, 5, {}},
  };
  Binning SP = {"species", "truth particle species",
     {"#pi^{#pm}", "K^{#pm}", "p / #bar{p}", "e^{#pm}", "#mu^{#pm}", "other"},
     {"pion", "kaon", "proton", "e", "mu", "other"}, {C02, C07, C03, C08, C01, C04}, 0, 2, {}};
  B.push_back(SP);
  const double ptEdge[5]  = {1.0, 2.0, 5.0, 10.0, 30.0};
  const double etaEdge[5] = {2.6, 2.8, 3.0, 3.2, 3.5};
  for (auto& b : B) {
    b.a.resize(b.lab.size());
    for (size_t i = 0; i < b.a.size(); ++i) {
      b.a[i].dt       = new TH1D(Form("dt_%s_%zu", b.tag, i), ";t_{track}^{HGTD} #minus t_{truth particle} [ps];Fraction of tracks / 4 ps", 200, -400, 400);
      b.a[i].pull     = new TH1D(Form("pull_%s_%zu", b.tag, i), ";(t_{track}^{HGTD} #minus t_{truth particle}) / #sigma_{t};Fraction of tracks", 200, -10, 10);
      b.a[i].dtNoPrim = new TH1D(Form("dtnp_%s_%zu", b.tag, i), ";t_{track}^{HGTD} #minus t_{truth particle} [ps];Fraction of the bin's tracks / 4 ps", 200, -400, 400);
      b.a[i].dtPrim   = new TH1D(Form("dtp_%s_%zu", b.tag, i), ";t_{track}^{HGTD} #minus t_{truth particle} [ps];Fraction of the bin's tracks / 4 ps", 200, -400, 400);
      b.a[i].dtc      = new TH1D(Form("dtc_%s_%zu", b.tag, i), ";t_{track}^{HGTD} #minus t_{truth particle} #minus #Deltat_{TOF}(m, p) [ps];Fraction of tracks / 4 ps", 200, -400, 400);
    }
  }

  Long64_t nEv = 0; long nAll = 0;
  while (reader.Next()) {
    if (maxEvents > 0 && nEv >= maxEvents) break;
    ++nEv;
    for (size_t idx = 0; idx < branch.trackTime.GetSize(); ++idx) {
      if (branch.trackTimeValid[idx] != 1) continue;
      const int p = branch.trackToParticle[idx];
      if (p < 0 || p >= (int)branch.particleT.GetSize()) continue;
      const int nh = branch.trackHgtdHits[idx];
      if (nh < 1) continue;
      const double pt = branch.trackPt[idx], aeta = std::abs(branch.trackEta[idx]);
      int bPt = 0;  while (bPt < 5 && pt >= ptEdge[bPt]) ++bPt;
      int bEta = 0; while (bEta < 5 && aeta >= etaEdge[bEta]) ++bEta;
      if (aeta < MIN_HGTD_ETA || aeta > MAX_HGTD_ETA) bEta = -1;   // a few timed tracks sit just outside
      const int pdg = std::abs((int)std::lround(particlePdgId[p]));
      const int bSp = pdg == 211 ? 0 : pdg == 321 ? 1 : pdg == 2212 ? 2 : pdg == 11 ? 3 : pdg == 13 ? 4 : 5;
      const int bins[4] = {std::min(nh, 4) - 1, bPt, bEta, bSp};
      const double dt = branch.trackTime[idx] - branch.particleT[p];
      // expected extra flight time of a massive particle against the beta = 1 hypothesis
      const double pmom = pt * std::cosh(aeta), m = particleM[p];
      const double L = Z_HGTD * std::cosh(aeta) / std::sinh(aeta);           // z / |cos theta|
      const double beta = pmom / std::sqrt(pmom * pmom + m * m);
      const double dtTof = (L / C_MM_PS) * (1.0 / beta - 1.0);
      const double dtc = dt - dtTof;
      const double res = branch.trackTimeRes[idx];
      const bool noPrim = (branch.trackPrimHits[idx] == 0);
      const bool bad = (res > 0) && (std::abs(dt) >= 3.0 * res);
      ++nAll;
      for (int k = 0; k < 4; ++k) {
        if (bins[k] < 0) continue;
        Acc& x = B[k].a[bins[k]];
        x.dt->Fill(dt);
        x.dtc->Fill(dtc);
        if (std::abs(dtc) < 60) ++x.n60c;
        if (!((res > 0) && (std::abs(dtc) >= 3.0 * res))) ++x.n3sc;
        if (res > 0) x.pull->Fill(dt / res);
        (noPrim ? x.dtNoPrim : x.dtPrim)->Fill(dt);
        ++x.n; x.sumRes += res; x.sumDt += dt;
        if (nh == 1) ++x.nOneHit;
        if (std::abs(dt) < 30) ++x.n30;
        if (std::abs(dt) < 60) ++x.n60;
        if (std::abs(dt) < 90) ++x.n90;
        if (!bad) ++x.n3s;
        if (noPrim) { ++x.nNoPrim; if (bad) ++x.nNoPrimBad; }
        else if (bad) ++x.nPrimBad;
      }
    }
  }

  // ── Console tables ──────────────────────────────────────────────────────────
  std::printf("\n%lld events, %ld tracks with an HGTD time and a truth particle\n", (long long)nEv, nAll);
  for (auto& b : B) {
    std::printf("\n=== HGTD track time - truth particle time, by %s ===\n", b.title);
    std::printf("%-9s %10s %7s %8s %9s %7s %10s %7s %7s %7s %8s | %12s %14s %14s %7s\n", "bin", "tracks", "share", "<sig_t>", "core sig", "peak", "sig/quoted",
                "<30ps", "<60ps", "<90ps", "<3 sig_t", "no prim. hit", "bad | no prim", "bad | >=1 prim", "1-hit");
    for (size_t i = 0; i < b.a.size(); ++i) {
      Acc& x = b.a[i];
      if (x.n == 0) continue;
      auto f = fitCore(x.dt, Form("fit_%s_%zu", b.tag, i));
      std::printf("%-9s %10ld %6.2f%% %7.1f %8.1f %7.1f %10.2f %6.1f%% %6.1f%% %6.1f%% %7.1f%% | %11.1f%% %13.1f%% %13.1f%% %6.1f%%\n",
                  b.row[i], x.n, 100.0 * x.n / nAll, x.sumRes / x.n, f[0], f[1], f[0] / (x.sumRes / x.n),
                  100.0 * x.n30 / x.n, 100.0 * x.n60 / x.n, 100.0 * x.n90 / x.n, 100.0 * x.n3s / x.n,
                  100.0 * x.nNoPrim / x.n, x.nNoPrim ? 100.0 * x.nNoPrimBad / x.nNoPrim : 0.0,
                  (x.n - x.nNoPrim) ? 100.0 * x.nPrimBad / (x.n - x.nNoPrim) : 0.0, 100.0 * x.nOneHit / x.n);
    }
  }
  std::printf("  core sig / peak: Gaussian + flat floor fitted over +-120 ps.  'bad' = |dt| >= 3 sigma_t.\n"
              "  no prim. hit: Track_nHGTDPrimaryHits == 0, i.e. none of the track's hits comes from its own particle.\n");
  std::printf("\n=== before vs after subtracting the expected mass delay, (L/c)(1/beta - 1) with the truth mass ===\n");
  std::printf("  pure Gaussian over +-90 ps: sigma, mean, chi2/ndf;  then the counted fractions\n");
  for (auto& b : B) {
    std::printf("  -- by %s --\n", b.title);
    std::printf("  %-9s %10s | %7s %7s %9s | %7s %7s %9s | %7s %7s | %8s %8s\n", "bin", "tracks", "sig", "mean", "chi2/ndf", "sig'", "mean'", "chi2/ndf'", "<60ps", "<60ps'", "<3sig", "<3sig'");
    for (size_t i = 0; i < b.a.size(); ++i) {
      Acc& x = b.a[i];
      if (x.n == 0) continue;
      auto f0 = fitGaus(x.dt, Form("g0_%s_%zu", b.tag, i));
      auto f1 = fitGaus(x.dtc, Form("g1_%s_%zu", b.tag, i));
      std::printf("  %-9s %10ld | %7.1f %7.1f %9.0f | %7.1f %7.1f %9.0f | %6.1f%% %6.1f%% | %7.1f%% %7.1f%%\n",
                  b.row[i], x.n, f0[0], f0[1], f0[2], f1[0], f1[1], f1[2],
                  100.0 * x.n60 / x.n, 100.0 * x.n60c / x.n, 100.0 * x.n3s / x.n, 100.0 * x.n3sc / x.n);
    }
  }

  // ── Pages ───────────────────────────────────────────────────────────────────
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  const std::string out = MyUtl::plotFilePath("", "track_dt_binned.pdf");
  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Print((out + "[").c_str());
  auto labels = [&](const Binning& b, const char* what) {
    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.032);
    tl.DrawLatex(0.18, 0.76, Form("%s, by %s", what, b.title));
    tl.DrawLatex(0.18, 0.71, "Tracks with an HGTD time and a truth particle");
  };
  for (auto& b : B) {
    const int NB = (int)b.a.size();
    auto drawSet = [&](bool pull, bool logy) {
      c->SetLogy(logy);
      double mx = 0; std::vector<TH1D*> hs;
      for (int i = 0; i < NB; ++i) {
        TH1D* h = (TH1D*)(pull ? b.a[i].pull : b.a[i].dt)->Clone(Form("n_%s_%d_%d_%d", b.tag, i, pull, logy));
        if (h->Integral(0, -1) > 0) h->Scale(1.0 / h->Integral(0, -1));
        h->SetLineColor(b.col[i]); h->SetLineWidth(2); h->SetMarkerSize(0);
        mx = std::max(mx, h->GetMaximum()); hs.push_back(h);
      }
      hs[0]->SetMaximum(logy ? 30 * mx : 1.35 * mx);
      if (logy) hs[0]->SetMinimum(2e-5);
      TLegend* leg = new TLegend(0.62, NB > 4 ? 0.62 : 0.70, 0.92, 0.90); StyleLegend(leg);
      for (int i = 0; i < NB; ++i) {
        hs[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        const double sh = 100.0 * b.a[i].n / nAll;
        leg->AddEntry(hs[i], Form(sh < 1 ? "%s (%.2f%%)" : "%s (%.0f%%)", b.lab[i], sh), "l");
      }
      leg->Draw();
      labels(b, pull ? "Pull, unit area" : "Residual, unit area");
      c->Print(out.c_str());
    };
    drawSet(false, true);
    drawSet(false, false);
    drawSet(true, true);
    // corrected residual, log
    {
      c->SetLogy(true);
      double mx = 0; std::vector<TH1D*> hs;
      for (int i = 0; i < NB; ++i) {
        TH1D* h = (TH1D*)b.a[i].dtc->Clone(Form("nc_%s_%d", b.tag, i));
        if (h->Integral(0, -1) > 0) h->Scale(1.0 / h->Integral(0, -1));
        h->SetLineColor(b.col[i]); h->SetLineWidth(2); h->SetMarkerSize(0);
        mx = std::max(mx, h->GetMaximum()); hs.push_back(h);
      }
      hs[0]->SetMaximum(30 * mx); hs[0]->SetMinimum(2e-5);
      TLegend* leg = new TLegend(0.62, NB > 4 ? 0.62 : 0.70, 0.92, 0.90); StyleLegend(leg);
      for (int i = 0; i < NB; ++i) {
        hs[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        const double sh = 100.0 * b.a[i].n / nAll;
        leg->AddEntry(hs[i], Form(sh < 1 ? "%s (%.2f%%)" : "%s (%.0f%%)", b.lab[i], sh), "l");
      }
      leg->Draw();
      labels(b, "Mass-corrected residual, unit area");
      c->Print(out.c_str());
    }
    // own-particle split for two bins, unit area per BIN (both curves of a bin share one normalisation)
    c->SetLogy(true);
    {
      const int bb[2] = {b.splitA, b.splitB};
      std::vector<TH1D*> h; double mx = 0;
      TLegend* leg = new TLegend(0.50, 0.70, 0.92, 0.90); StyleLegend(leg);
      for (int k = 0; k < 2; ++k) {
        const Acc& x = b.a[bb[k]];
        for (int np = 0; np < 2; ++np) {
          TH1D* g = (TH1D*)(np ? x.dtNoPrim : x.dtPrim)->Clone(Form("sp_%s_%d_%d", b.tag, k, np));
          if (x.n > 0) g->Scale(1.0 / x.n);
          g->SetLineColor(b.col[bb[k]]); g->SetLineWidth(2); g->SetLineStyle(np ? 2 : 1);
          mx = std::max(mx, g->GetMaximum()); h.push_back(g);
          leg->AddEntry(g, Form("%s, %s", b.lab[bb[k]], np ? "no hit from own particle" : "#geq 1 hit from own particle"), "l");
        }
      }
      h[0]->SetMaximum(30 * mx); h[0]->SetMinimum(2e-6);
      for (size_t i = 0; i < h.size(); ++i) h[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
      leg->Draw();
      labels(b, "Residual by hit provenance");
      c->Print(out.c_str());
    }
  }
  c->Print((out + "]").c_str());
  std::cout << "Wrote " << out << "\n";
  return 0;
}
