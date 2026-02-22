// evaluate_ml_lowtrack.cxx
// Evaluates neural network scores on clusters in events with <= 4 hard scatter tracks
// This helps understand ML model performance in low-track regimes

#include <RtypesCore.h>
#include <TROOT.h>
#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

using namespace MyUtl;

auto main() -> int {
  gStyle->SetOptStat(0);

  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  ROOT::EnableImplicitMT();

  // Histograms for ML score distributions
  TH1D* h_ml_score_all = new TH1D("h_ml_score_all",
    "ML Score (All Clusters, N_{HS} #leq 4);ML Score;Clusters", 50, 0, 1);
  TH1D* h_ml_score_hs = new TH1D("h_ml_score_hs",
    "ML Score (HS Clusters, N_{HS} #leq 4);ML Score;Clusters", 50, 0, 1);
  TH1D* h_ml_score_pu = new TH1D("h_ml_score_pu",
    "ML Score (PU Clusters, N_{HS} #leq 4);ML Score;Clusters", 50, 0, 1);
  TH1D* h_ml_score_chosen = new TH1D("h_ml_score_chosen",
    "ML Score (Chosen Cluster, N_{HS} #leq 4);ML Score;Events", 50, 0, 1);

  // Histograms by number of HS tracks
  TH1D* h_ml_score_1hs = new TH1D("h_ml_score_1hs", "ML Score (N_{HS} = 1);ML Score;Clusters", 50, 0, 1);
  TH1D* h_ml_score_2hs = new TH1D("h_ml_score_2hs", "ML Score (N_{HS} = 2);ML Score;Clusters", 50, 0, 1);
  TH1D* h_ml_score_3hs = new TH1D("h_ml_score_3hs", "ML Score (N_{HS} = 3);ML Score;Clusters", 50, 0, 1);
  TH1D* h_ml_score_4hs = new TH1D("h_ml_score_4hs", "ML Score (N_{HS} = 4);ML Score;Clusters", 50, 0, 1);

  // 2D: ML score vs cluster purity
  TH2D* h_ml_vs_purity = new TH2D("h_ml_vs_purity",
    "ML Score vs Cluster Purity (N_{HS} #leq 4);Cluster Purity;ML Score",
    50, 0, 1, 50, 0, 1);

  // 2D: ML score vs N_HS tracks
  TH2D* h_ml_vs_nhs = new TH2D("h_ml_vs_nhs",
    "ML Score vs N_{HS} Tracks;N_{HS} Tracks;ML Score",
    5, -0.5, 5.5, 50, 0, 1);

  // Profile histograms: average ML score vs purity and N_HS
  TProfile* p_ml_vs_purity = new TProfile("p_ml_vs_purity",
    "Average ML Score vs Cluster Purity;Cluster Purity;#LT ML Score #GT",
    20, 0, 1, 0, 1);
  TProfile* p_ml_vs_nhs = new TProfile("p_ml_vs_nhs",
    "Average ML Score vs N_{HS} Tracks;N_{HS} Tracks;#LT ML Score #GT",
    4, -0.5, 4.5, 0, 1);

  // Efficiency histograms: correct selection when ML score > 0.5
  TH1D* h_eff_total = new TH1D("h_eff_total", "Total Events;N_{HS} Tracks;Events", 5, 0.5, 5.5);
  TH1D* h_eff_pass = new TH1D("h_eff_pass", "Correct Selection (ML > 0.5);N_{HS} Tracks;Events", 5, 0.5, 5.5);
  TH1D* h_eff_pass_timing = new TH1D("h_eff_pass_timing",
    "Correct Timing (ML > 0.5 && |#Deltat| < 60 ps);N_{HS} Tracks;Events", 5, 0.5, 5.5);

  // Counters
  Long64_t nProcessed = 0;
  Long64_t nLowTrack = 0;
  Long64_t nCorrectSelection = 0;
  Long64_t nMLAboveThreshold = 0;
  Long64_t nCorrectWithMLThreshold = 0;

  std::cout << "Starting Event Loop for ML Evaluation (N_HS <= 4)\n";
  Long64_t nEvent = chain.GetEntries();

  while (reader.Next()) {
    Long64_t readNum = chain.GetReadEntry() + 1;
    if (readNum % 5000 == 0)
      std::cout << "Progress: " << readNum << "/" << nEvent << "\n";

    // Apply basic cuts
    if (not branch.passBasicCuts()) continue;
    if (not branch.passJetPtCut()) continue;

    nProcessed++;

    // Get associated tracks
    std::vector<int> tracks = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT,3.0);

    // Count forward tracks
    int nForwardTrack = 0, nForwardTrackHS = 0, nForwardTrackPU = 0;
    branch.countForwardTracks(nForwardTrack, nForwardTrackHS, nForwardTrackPU, tracks, true);

    // Only process events with <= 4 HS tracks
    if (nForwardTrackHS > 4) continue;

    nLowTrack++;

    // Perform clustering
    std::vector<Cluster> clusters = clusterTracksInTime(
      tracks, &branch, 3.0,
      false, true, 10.0,  // useSmearedTimes=false, checkValidTimes=true
      true, false);        // useCone=true, useZ0=false

    if (clusters.empty()) continue;

    // Fill total efficiency histogram
    h_eff_total->Fill(nForwardTrackHS);

    // Find the chosen cluster (highest ML score)
    Cluster chosenCluster = clusters[0];
    double maxMLScore = chosenCluster.scores.at(Score::TESTML);
    for (auto& cluster : clusters) {
      double mlScore = cluster.scores.at(Score::TESTML);
      if (mlScore > maxMLScore) {
        maxMLScore = mlScore;
        chosenCluster = cluster;
      }
    }

    // Fill ML score for chosen cluster
    h_ml_score_chosen->Fill(maxMLScore);

    // Check if chosen cluster is correct (passes efficiency)
    bool isCorrect = chosenCluster.passEfficiency(&branch);
    if (isCorrect) nCorrectSelection++;

    // Check ML threshold
    if (maxMLScore > 0.5) {
      nMLAboveThreshold++;
      if (isCorrect) {
        nCorrectWithMLThreshold++;
        h_eff_pass->Fill(nForwardTrackHS);

        // Also check timing
        double timeDiff = std::abs(chosenCluster.values.at(0) - branch.truthVtxTime[0]);
        if (timeDiff < 60.0) {
          h_eff_pass_timing->Fill(nForwardTrackHS);
        }
      }
    }

    // Fill histograms for all clusters
    for (auto& cluster : clusters) {
      double mlScore = cluster.scores.at(Score::TESTML);
      double purity = cluster.purity;

      h_ml_score_all->Fill(mlScore);
      h_ml_vs_purity->Fill(purity, mlScore);
      h_ml_vs_nhs->Fill(nForwardTrackHS, mlScore);
      p_ml_vs_purity->Fill(purity, mlScore);
      p_ml_vs_nhs->Fill(nForwardTrackHS, mlScore);

      // Classify cluster as HS or PU based on purity
      if (purity > 0.5) {
        h_ml_score_hs->Fill(mlScore);
      } else {
        h_ml_score_pu->Fill(mlScore);
      }

      // Fill by N_HS
      switch (nForwardTrackHS) {
        case 1: h_ml_score_1hs->Fill(mlScore); break;
        case 2: h_ml_score_2hs->Fill(mlScore); break;
        case 3: h_ml_score_3hs->Fill(mlScore); break;
        case 4: h_ml_score_4hs->Fill(mlScore); break;
      }
    }
  }

  // Print summary statistics
  std::cout << "\n========================================\n";
  std::cout << "ML Evaluation Summary (N_HS <= 4)\n";
  std::cout << "========================================\n";
  std::cout << "Total events processed: " << nProcessed << "\n";
  std::cout << "Events with 1-4 HS tracks: " << nLowTrack << "\n";
  std::cout << "Correct cluster selection (all): " << nCorrectSelection
            << " (" << 100.0 * nCorrectSelection / nLowTrack << "%)\n";
  std::cout << "Events with ML score > 0.5: " << nMLAboveThreshold
            << " (" << 100.0 * nMLAboveThreshold / nLowTrack << "%)\n";
  std::cout << "Correct selection with ML > 0.5: " << nCorrectWithMLThreshold
            << " (" << 100.0 * nCorrectWithMLThreshold / nMLAboveThreshold << "% of ML>0.5)\n";
  std::cout << "========================================\n";

  // Create output plots
  TCanvas* canvas = new TCanvas("canvas", "ML Evaluation", 1200, 800);
  gErrorIgnoreLevel = kWarning;

  // Plot 1: ML score distributions (HS vs PU)
  canvas->Clear();
  canvas->SetLogy(true);
  h_ml_score_hs->SetLineColor(kBlue);
  h_ml_score_hs->SetLineWidth(2);
  h_ml_score_pu->SetLineColor(kRed);
  h_ml_score_pu->SetLineWidth(2);

  h_ml_score_hs->SetMaximum(1.5 * std::max(h_ml_score_hs->GetMaximum(), h_ml_score_pu->GetMaximum()));
  h_ml_score_hs->Draw("HIST");
  h_ml_score_pu->Draw("HIST SAME");

  TLegend* leg1 = new TLegend(0.15, 0.75, 0.45, 0.88);
  leg1->AddEntry(h_ml_score_hs, "HS Clusters (purity > 0.5)", "l");
  leg1->AddEntry(h_ml_score_pu, "PU Clusters (purity < 0.5)", "l");
  leg1->Draw();

  canvas->SaveAs("../figs/ml_score_hs_vs_pu_lowtrack.pdf");

  // Plot 2: ML score by N_HS tracks
  canvas->Clear();
  canvas->SetLogy(false);
  h_ml_score_1hs->SetLineColor(kBlue);
  h_ml_score_2hs->SetLineColor(kGreen+2);
  h_ml_score_3hs->SetLineColor(kOrange+1);
  h_ml_score_4hs->SetLineColor(kRed);

  h_ml_score_1hs->SetLineWidth(2);
  h_ml_score_2hs->SetLineWidth(2);
  h_ml_score_3hs->SetLineWidth(2);
  h_ml_score_4hs->SetLineWidth(2);

  // Normalize to unity
  if (h_ml_score_1hs->Integral() > 0) h_ml_score_1hs->Scale(1.0/h_ml_score_1hs->Integral());
  if (h_ml_score_2hs->Integral() > 0) h_ml_score_2hs->Scale(1.0/h_ml_score_2hs->Integral());
  if (h_ml_score_3hs->Integral() > 0) h_ml_score_3hs->Scale(1.0/h_ml_score_3hs->Integral());
  if (h_ml_score_4hs->Integral() > 0) h_ml_score_4hs->Scale(1.0/h_ml_score_4hs->Integral());

  double maxY = 0;
  maxY = std::max(maxY, h_ml_score_1hs->GetMaximum());
  maxY = std::max(maxY, h_ml_score_2hs->GetMaximum());
  maxY = std::max(maxY, h_ml_score_3hs->GetMaximum());
  maxY = std::max(maxY, h_ml_score_4hs->GetMaximum());

  h_ml_score_1hs->SetMaximum(1.3 * maxY);
  h_ml_score_1hs->GetYaxis()->SetTitle("Normalized");
  h_ml_score_1hs->Draw("HIST");
  h_ml_score_2hs->Draw("HIST SAME");
  h_ml_score_3hs->Draw("HIST SAME");
  h_ml_score_4hs->Draw("HIST SAME");

  TLegend* leg2 = new TLegend(0.15, 0.65, 0.40, 0.88);
  leg2->AddEntry(h_ml_score_1hs, "N_{HS} = 1", "l");
  leg2->AddEntry(h_ml_score_2hs, "N_{HS} = 2", "l");
  leg2->AddEntry(h_ml_score_3hs, "N_{HS} = 3", "l");
  leg2->AddEntry(h_ml_score_4hs, "N_{HS} = 4", "l");
  leg2->Draw();

  canvas->SaveAs("../figs/ml_score_by_nhs.pdf");

  // Plot 3: 2D ML score vs purity
  canvas->Clear();
  canvas->SetLogy(false);
  canvas->SetRightMargin(0.15);
  h_ml_vs_purity->Draw("COLZ");
  canvas->SaveAs("../figs/ml_score_vs_purity_lowtrack.pdf");

  // Plot 4: 2D ML score vs N_HS
  canvas->Clear();
  h_ml_vs_nhs->Draw("COLZ");
  canvas->SaveAs("../figs/ml_score_vs_nhs.pdf");

  // Plot 4a: Profile - average ML score vs purity
  canvas->Clear();
  canvas->SetRightMargin(0.05);
  p_ml_vs_purity->SetLineColor(kBlue);
  p_ml_vs_purity->SetLineWidth(2);
  p_ml_vs_purity->SetMarkerStyle(20);
  p_ml_vs_purity->SetMarkerColor(kBlue);
  p_ml_vs_purity->GetYaxis()->SetRangeUser(0, 1);
  p_ml_vs_purity->Draw("E1");

  // Add horizontal line at 0.5 threshold
  TLine* hline1 = new TLine(0, 0.5, 1, 0.5);
  hline1->SetLineColor(kRed);
  hline1->SetLineStyle(2);
  hline1->SetLineWidth(2);
  hline1->Draw();

  canvas->SaveAs("../figs/ml_score_profile_vs_purity.pdf");

  // Plot 4b: Profile - average ML score vs N_HS
  canvas->Clear();
  p_ml_vs_nhs->SetLineColor(kBlue);
  p_ml_vs_nhs->SetLineWidth(2);
  p_ml_vs_nhs->SetMarkerStyle(20);
  p_ml_vs_nhs->SetMarkerColor(kBlue);
  p_ml_vs_nhs->GetYaxis()->SetRangeUser(0, 1);
  p_ml_vs_nhs->Draw("E1");

  // Add horizontal line at 0.5 threshold
  TLine* hline2 = new TLine(-0.5, 0.5, 4.5, 0.5);
  hline2->SetLineColor(kRed);
  hline2->SetLineStyle(2);
  hline2->SetLineWidth(2);
  hline2->Draw();

  canvas->SaveAs("../figs/ml_score_profile_vs_nhs.pdf");

  // Plot 5: ML score for chosen clusters
  canvas->Clear();
  canvas->SetRightMargin(0.05);
  h_ml_score_chosen->SetLineColor(kBlack);
  h_ml_score_chosen->SetLineWidth(2);
  h_ml_score_chosen->Draw("HIST");

  // Add vertical line at 0.5 threshold
  TLine* line = new TLine(0.5, 0, 0.5, h_ml_score_chosen->GetMaximum());
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  canvas->SaveAs("../figs/ml_score_chosen_cluster_lowtrack.pdf");

  // Plot 6: Efficiency vs N_HS
  canvas->Clear();
  TH1D* h_eff = (TH1D*)h_eff_pass->Clone("h_eff");
  h_eff->Divide(h_eff_total);
  h_eff->SetTitle("Selection Efficiency (ML > 0.5);N_{HS} Tracks;Efficiency");
  h_eff->SetLineColor(kBlue);
  h_eff->SetLineWidth(2);
  h_eff->SetMarkerStyle(20);
  h_eff->SetMarkerColor(kBlue);
  h_eff->GetYaxis()->SetRangeUser(0, 1.2);
  h_eff->Draw("E1");

  TH1D* h_eff_time = (TH1D*)h_eff_pass_timing->Clone("h_eff_time");
  h_eff_time->Divide(h_eff_total);
  h_eff_time->SetLineColor(kRed);
  h_eff_time->SetLineWidth(2);
  h_eff_time->SetMarkerStyle(21);
  h_eff_time->SetMarkerColor(kRed);
  h_eff_time->Draw("E1 SAME");

  TLegend* leg3 = new TLegend(0.15, 0.75, 0.55, 0.88);
  leg3->AddEntry(h_eff, "Correct selection (ML > 0.5)", "lep");
  leg3->AddEntry(h_eff_time, "Correct timing (|#Deltat| < 60 ps)", "lep");
  leg3->Draw();

  canvas->SaveAs("../figs/ml_efficiency_vs_nhs.pdf");

  // Save histograms to ROOT file
  TFile* outFile = new TFile("ml_evaluation_lowtrack.root", "RECREATE");
  h_ml_score_all->Write();
  h_ml_score_hs->Write();
  h_ml_score_pu->Write();
  h_ml_score_chosen->Write();
  h_ml_score_1hs->Write();
  h_ml_score_2hs->Write();
  h_ml_score_3hs->Write();
  h_ml_score_4hs->Write();
  h_ml_vs_purity->Write();
  h_ml_vs_nhs->Write();
  p_ml_vs_purity->Write();
  p_ml_vs_nhs->Write();
  h_eff_total->Write();
  h_eff_pass->Write();
  h_eff_pass_timing->Write();
  outFile->Close();

  std::cout << "\nOutput saved to:\n";
  std::cout << "  - ../figs/ml_score_hs_vs_pu_lowtrack.pdf\n";
  std::cout << "  - ../figs/ml_score_by_nhs.pdf\n";
  std::cout << "  - ../figs/ml_score_vs_purity_lowtrack.pdf\n";
  std::cout << "  - ../figs/ml_score_vs_nhs.pdf\n";
  std::cout << "  - ../figs/ml_score_profile_vs_purity.pdf\n";
  std::cout << "  - ../figs/ml_score_profile_vs_nhs.pdf\n";
  std::cout << "  - ../figs/ml_score_chosen_cluster_lowtrack.pdf\n";
  std::cout << "  - ../figs/ml_efficiency_vs_nhs.pdf\n";
  std::cout << "  - ml_evaluation_lowtrack.root\n";

  return 0;
}
