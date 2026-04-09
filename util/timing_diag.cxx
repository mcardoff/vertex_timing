// timing_diag.cxx
// Diagnostic: print trackTimeRes distribution and |dt| for HS vs PU tracks
// to understand why Ideal t0 ROC doesn't match ATL-HGTD-PUB-2022-001.

#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp>

#include "clustering_includes.h"
#include "clustering_structs.h"

#define debug false

using namespace MyUtl;

int main() {
  gStyle->SetOptStat(0);

  TChain chain("ntuple");
  for (const auto& entry : boost::filesystem::directory_iterator("../../ntuple-hgtd")) {
    chain.Add(entry.path().c_str());
  }

  TTreeReader reader(&chain);
  MyUtl::BranchPointerWrapper branch(reader);

  // Histograms for timing resolution
  TH1D* h_timeRes_all  = new TH1D("h_timeRes_all",  "Track time resolution (all HGTD tracks);#sigma_{t} [ps];Entries",  200, 0, 200);
  TH1D* h_timeRes_hs   = new TH1D("h_timeRes_hs",   "Track time resolution (HS tracks);#sigma_{t} [ps];Entries",         200, 0, 200);
  TH1D* h_timeRes_pu   = new TH1D("h_timeRes_pu",   "Track time resolution (PU tracks);#sigma_{t} [ps];Entries",         200, 0, 200);

  // Histograms for |dt| = |trackTime - truthVtxTime|
  TH1D* h_dt_hs = new TH1D("h_dt_hs", "|t_{track} - t_{0}^{true}| (HS tracks);|#Delta t| [ps];Entries", 300, 0, 600);
  TH1D* h_dt_pu = new TH1D("h_dt_pu", "|t_{track} - t_{0}^{true}| (PU tracks);|#Delta t| [ps];Entries", 300, 0, 600);

  // Histograms for nsigma = |dt| / sigma_t
  TH1D* h_nsig_hs = new TH1D("h_nsig_hs", "|#Delta t|/#sigma_{t} (HS tracks);n#sigma;Entries", 100, 0, 10);
  TH1D* h_nsig_pu = new TH1D("h_nsig_pu", "|#Delta t|/#sigma_{t} (PU tracks);n#sigma;Entries", 100, 0, 10);

  // Check units: print a few raw values
  int n_printed = 0;
  int n_events = 0;

  // Counters
  long n_hs_tracks = 0, n_pu_tracks = 0;
  double sum_timeRes_hs = 0, sum_timeRes_pu = 0;
  double sum_dt_hs = 0, sum_dt_pu = 0;

  while (reader.Next()) {
    double tru_vtx_t = branch.truthVtxTime[0];
    double pri_vtx_z = branch.recoVtxZ[0];

    // Skip if reco vertex is far from truth
    if (std::abs(branch.truthVtxZ[0] - branch.recoVtxZ[0]) > 2.0) continue;

    if (n_printed < 5) {
      std::cout << "Event " << n_events
                << "  truthVtxTime[0]=" << tru_vtx_t
                << "  recoVtxTime[0]=" << branch.recoVtxTime[0]
                << "  recoVtxTimeRes[0]=" << branch.recoVtxTimeRes[0]
                << std::endl;
      // Print first few tracks
      int nprint_trk = 0;
      for (int i = 0; i < (int)branch.trackEta.GetSize() && nprint_trk < 3; i++) {
        double eta = std::abs(branch.trackEta[i]);
        if (eta < MIN_ABS_ETA_TRACK || eta > MAX_ABS_ETA_TRACK) continue;
        if (branch.trackTimeValid[i] != 1) continue;
        std::cout << "  trk[" << i << "]  time=" << branch.trackTime[i]
                  << "  timeRes=" << branch.trackTimeRes[i]
                  << "  truthVtx_idx=" << branch.trackToTruthvtx[i]
                  << "  |dt|=" << std::abs(branch.trackTime[i] - tru_vtx_t)
                  << "  |dt|/sig=" << std::abs(branch.trackTime[i] - tru_vtx_t) / branch.trackTimeRes[i]
                  << std::endl;
        nprint_trk++;
      }
      n_printed++;
    }

    for (int i = 0; i < (int)branch.trackEta.GetSize(); i++) {
      double eta = std::abs(branch.trackEta[i]);
      double pt  = branch.trackPt[i];

      if (eta < MIN_ABS_ETA_TRACK || eta > MAX_ABS_ETA_TRACK) continue;
      if (pt < MIN_TRACK_PT || pt > MAX_TRACK_PT) continue;
      if (!branch.trackQuality[i]) continue;
      if (branch.trackTimeValid[i] != 1) continue;

      double t   = branch.trackTime[i];
      double sig = branch.trackTimeRes[i];
      double dt  = std::abs(t - tru_vtx_t);
      double nsig = (sig > 0) ? dt / sig : -1;
      bool isHS = (branch.trackToTruthvtx[i] == 0);

      h_timeRes_all->Fill(sig);

      if (isHS) {
        h_timeRes_hs->Fill(sig);
        h_dt_hs->Fill(dt);
        if (nsig >= 0) h_nsig_hs->Fill(nsig);
        sum_timeRes_hs += sig;
        sum_dt_hs += dt;
        n_hs_tracks++;
      } else {
        h_timeRes_pu->Fill(sig);
        h_dt_pu->Fill(dt);
        if (nsig >= 0) h_nsig_pu->Fill(nsig);
        sum_timeRes_pu += sig;
        sum_dt_pu += dt;
        n_pu_tracks++;
      }
    }
    n_events++;
  }

  // Print summary
  std::cout << "\n========================================" << std::endl;
  std::cout << "TIMING DIAGNOSTIC SUMMARY" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Events processed: " << n_events << std::endl;
  std::cout << "\nHS HGTD tracks:  " << n_hs_tracks << std::endl;
  if (n_hs_tracks > 0) {
    std::cout << "  mean timeRes:  " << sum_timeRes_hs / n_hs_tracks << " ps" << std::endl;
    std::cout << "  mean |dt|:     " << sum_dt_hs / n_hs_tracks << " ps" << std::endl;
    std::cout << "  timeRes peak:  " << h_timeRes_hs->GetBinCenter(h_timeRes_hs->GetMaximumBin()) << " ps" << std::endl;
    std::cout << "  |dt| peak:     " << h_dt_hs->GetBinCenter(h_dt_hs->GetMaximumBin()) << " ps" << std::endl;
  }
  std::cout << "\nPU HGTD tracks:  " << n_pu_tracks << std::endl;
  if (n_pu_tracks > 0) {
    std::cout << "  mean timeRes:  " << sum_timeRes_pu / n_pu_tracks << " ps" << std::endl;
    std::cout << "  mean |dt|:     " << sum_dt_pu / n_pu_tracks << " ps" << std::endl;
    std::cout << "  timeRes peak:  " << h_timeRes_pu->GetBinCenter(h_timeRes_pu->GetMaximumBin()) << " ps" << std::endl;
    std::cout << "  |dt| peak:     " << h_dt_pu->GetBinCenter(h_dt_pu->GetMaximumBin()) << " ps" << std::endl;
  }

  // Fraction of HS tracks passing 2.5 sigma cut
  double hs_pass = h_nsig_hs->Integral(1, h_nsig_hs->FindBin(2.5)) / h_nsig_hs->Integral();
  double pu_pass = h_nsig_pu->Integral(1, h_nsig_pu->FindBin(2.5)) / h_nsig_pu->Integral();
  std::cout << "\nFraction passing |dt|/sig < 2.5:" << std::endl;
  std::cout << "  HS: " << hs_pass * 100.0 << "%" << std::endl;
  std::cout << "  PU: " << pu_pass * 100.0 << "%" << std::endl;

  // Save diagnostic plots
  TCanvas* c = new TCanvas("c", "Timing Diagnostics", 1200, 800);
  c->Divide(3, 2);

  c->cd(1);
  h_timeRes_all->SetLineColor(kBlack);
  h_timeRes_hs->SetLineColor(kRed);
  h_timeRes_pu->SetLineColor(kBlue);
  h_timeRes_all->Draw("hist");
  h_timeRes_hs->Draw("hist same");
  h_timeRes_pu->Draw("hist same");
  auto leg1 = new TLegend(0.5, 0.6, 0.9, 0.9);
  leg1->AddEntry(h_timeRes_all, "All HGTD tracks", "l");
  leg1->AddEntry(h_timeRes_hs,  "HS tracks",       "l");
  leg1->AddEntry(h_timeRes_pu,  "PU tracks",       "l");
  leg1->Draw();

  c->cd(2);
  h_dt_hs->SetLineColor(kRed);
  h_dt_pu->SetLineColor(kBlue);
  h_dt_pu->Draw("hist");
  h_dt_hs->Draw("hist same");
  auto leg2 = new TLegend(0.5, 0.6, 0.9, 0.9);
  leg2->AddEntry(h_dt_hs, "HS tracks", "l");
  leg2->AddEntry(h_dt_pu, "PU tracks", "l");
  leg2->Draw();

  c->cd(3);
  h_nsig_hs->SetLineColor(kRed);
  h_nsig_pu->SetLineColor(kBlue);
  // normalize
  if (h_nsig_hs->Integral() > 0) h_nsig_hs->Scale(1.0 / h_nsig_hs->Integral());
  if (h_nsig_pu->Integral() > 0) h_nsig_pu->Scale(1.0 / h_nsig_pu->Integral());
  h_nsig_pu->Draw("hist");
  h_nsig_hs->Draw("hist same");
  auto leg3 = new TLegend(0.5, 0.6, 0.9, 0.9);
  leg3->AddEntry(h_nsig_hs, "HS tracks", "l");
  leg3->AddEntry(h_nsig_pu, "PU tracks", "l");
  leg3->Draw();

  c->SaveAs("../figs/timing_diag.pdf");
  std::cout << "\nSaved: ../figs/timing_diag.pdf" << std::endl;

  return 0;
}
