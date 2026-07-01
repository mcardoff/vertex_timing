// track_dt.cxx — Inclusive per-track timing residual.
//
// Track-level analogue of the inclusive vertex-t0 Delta t distributions:
// for every track that has a valid HGTD time AND a linked truth particle,
// fill (HGTD track time - truth particle prodVtx time).  No event selection
// and no track pT / quality cuts are applied — every qualifying track in the
// sample contributes one entry.
//
// The distribution is fit with the shared-mean double Gaussian (FIT_DBLGAUS)
// reused from plotting_utilities.h, reproducing the sigma_1 / sigma_2 dgaus
// annotations of the original plot.
//
// Output: figs/track_dt.pdf

#include <TCanvas.h>
#include <TChain.h>
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
#include "AtlasStyle.h"
#include "AtlasLabels.h"

#define debug false

using namespace MyUtl;

int main() {
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found in ../../ntuple-hgtd/. Aborting.\n";
    return 1;
  }

  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TH1D* h = new TH1D("h_track_dt",
                     "Inclusive Track t - Truth t;#Delta t[ps];Entries",
                     400, -400.0, 400.0);

  long n_tracks = 0, n_filled = 0;
  long n_acc = 0, n_acc_timed = 0;  // HGTD-acceptance tracks; those assigned a time
  long n_within90 = 0;              // of timed+truth tracks, |dt| < 90 ps

  while (reader.Next()) {
    for (size_t idx = 0; idx < branch.trackTime.GetSize(); ++idx) {
      ++n_tracks;

      // Timing efficiency in HGTD acceptance (|eta| in [MIN,MAX]_HGTD_ETA):
      // fraction of acceptance tracks that get a valid HGTD time.
      double aeta = std::abs(branch.trackEta[idx]);
      bool inAcc = (aeta >= MIN_HGTD_ETA && aeta <= MAX_HGTD_ETA);
      if (inAcc) {
        ++n_acc;
        if (branch.trackTimeValid[idx] == 1) ++n_acc_timed;
      }

      if (branch.trackTimeValid[idx] != 1) continue;  // require an HGTD time
      int part_idx = branch.trackToParticle[idx];
      if (part_idx < 0 || part_idx >= (int)branch.particleT.GetSize())
        continue;                                      // require a truth time
      double dt = branch.trackTime[idx] - branch.particleT[part_idx];
      h->Fill(dt);
      ++n_filled;
      if (std::abs(dt) < 90.0) ++n_within90;           // within 90 ps of truth
    }
  }

  double eff_timed  = n_acc    > 0 ? (double)n_acc_timed / n_acc : 0.0;
  double frac_w90   = n_filled > 0 ? (double)n_within90  / n_filled : 0.0;

  // ── Double-Gaussian fit (shared mean free, not displayed). ───────────────────
  double xlo = h->GetXaxis()->GetXmin(), xhi = h->GetXaxis()->GetXmax();
  TF1* fit = new TF1("track_dt_fit",
                     "[1]*TMath::Exp(-0.5*((x-[0])/[3])^2)+"
                     "[2]*TMath::Exp(-0.5*((x-[0])/[4])^2)", xlo, xhi);
  fit->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
  fit->SetParameters(h->GetMean(), 0.8 * h->GetMaximum(), 1e-2 * h->GetMaximum(),
                     13.0, h->GetStdDev());
  fit->SetParLimits(0, -50.0, 50.0);  // mean free, constrained +/-50 ps
  fit->SetParLimits(1, 0, 1.E6);
  fit->SetParLimits(2, 0, 1.E6);
  fit->SetParLimits(3, 5.0, 1.E6);    // core sigma >= 5 ps
  fit->SetParLimits(4, 5.0, 1.E6);    // tail sigma >= 5 ps
  fit->SetNpx(1000);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(2);
  h->Fit(fit, "RQI");

  // Order so Sigma1 is the narrow core (swap only if the fit inverted them).
  bool   swap    = fit->GetParameter("Sigma1") > fit->GetParameter("Sigma2");
  double sigCore = fit->GetParameter(swap ? "Sigma2" : "Sigma1");
  double sigTail = fit->GetParameter(swap ? "Sigma1" : "Sigma2");
  double errCore = fit->GetParError(swap ? "Sigma2" : "Sigma1");
  double errTail = fit->GetParError(swap ? "Sigma1" : "Sigma2");

  std::vector<TString> sigmaLines = {
    TString::Format("#sigma_{1}^{dgaus}=%.2f #pm (%.2f%%)", sigCore, 100 * errCore / sigCore),
    TString::Format("#sigma_{2}^{dgaus}=%.2f #pm (%.2f%%)", sigTail, 100 * errTail / sigTail),
    TString::Format("Timing eff. = %.1f%%", 100 * eff_timed),
    TString::Format("Core frac. = %.1f%%",  100 * frac_w90),
  };

  // ── Draw ─────────────────────────────────────────────────────────────────────
  boost::filesystem::create_directories("../figs");
  const TString out_pdf = "../figs/track_dt.pdf";

  TCanvas* canvas = new TCanvas("canvas", "Inclusive Track Delta t", 800, 600);
  canvas->SetLogy(true);

  h->SetLineColor(kBlue + 2);
  h->SetLineWidth(2);
  h->SetMaximum(5.0 * h->GetMaximum());
  h->Draw("HIST");
  fit->Draw("SAME");

  TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.90);
  StyleLegend(leg);
  leg->AddEntry(h,   "Histogram",            "l");
  leg->AddEntry(fit, "Double Gaussian Fit",  "l");
  leg->Draw();

  ATLASLabel(0.18, 0.88, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.82);

  // sigma_1 / sigma_2 dgaus annotations, top-left below the labels.
  TLatex tl;
  tl.SetNDC();
  tl.SetTextFont(42);
  tl.SetTextSize(0.035);
  double y = 0.74;
  for (const auto& line : sigmaLines) {
    tl.DrawLatex(0.18, y, line);
    y -= 0.055;
  }

  canvas->Print(out_pdf);

  // ── Console summary ──────────────────────────────────────────────────────────
  std::cout << "\n=== Inclusive Track t - Truth t ===\n";
  std::cout << "Tracks seen              : " << n_tracks << '\n';
  std::cout << "Tracks in HGTD acceptance: " << n_acc << '\n';
  std::cout << "  ... assigned a time    : " << n_acc_timed
            << " (timing eff. = " << 100 * eff_timed << "%)\n";
  std::cout << "Tracks filled            : " << n_filled
            << " (valid HGTD time + truth particle link)\n";
  std::cout << "  ... |dt| < 90 ps       : " << n_within90
            << " (" << 100 * frac_w90 << "%)\n";
  std::cout << "sigma_1 (core)           : " << sigCore << " ps\n";
  std::cout << "sigma_2 (tail)           : " << sigTail << " ps\n";
  std::cout << "Wrote " << out_pdf << "\n";
  return 0;
}
