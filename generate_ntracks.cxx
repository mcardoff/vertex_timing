#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TStyle.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <boost/filesystem.hpp>

#define debug false

// colors
auto c1 = TColor::GetColor("#3f90da");
auto c2 = TColor::GetColor("#ffa90e");
auto c3 = TColor::GetColor("#bd1f01");
auto c4 = TColor::GetColor("#94a4a2");
auto c5 = TColor::GetColor("#832db6");

// cut variables
auto min_jets          =  2;    // min number of jets
auto min_jetpt         =  30.0; // self explanatory
auto min_abs_eta       =  2.0;  // min eta for a vtx to be considered "forward"
auto min_abs_track_eta =  2.4;  // min eta for a track to be considered "forward"
auto min_track_pt      =  1.0;  // track_pt > 1.0 GeV
auto max_vtx_dz        =  2.0;  // max error for reco HS vertex z
auto max_nsigma        =  3.0;  // how close a track can be to PV

void setup_chain(TChain &chain) {
  // VBF H->Invisible sample
  for (const auto& entry : boost::filesystem::directory_iterator("./ntuple")) {
    if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
    chain.Add(entry.path().c_str());
    // break; // add 1
  }
  
  if (chain.GetEntries() == 0) {
    std::cerr << "No ROOT files found in directory: " << std::endl;
    return;
  }
}

void draw_class(TH1I* h1, TH1I* h2, TH1I* h3, TLegend* legend) {
  h1->SetMaximum(1.1*std::max({h1->GetMaximum(),h2->GetMaximum(),h3->GetMaximum()}));

  h1->SetLineColor(c1);
  h1->SetLineWidth(2);
  h1->Draw("HIST");

  h2->SetLineColor(c2);
  h2->SetLineWidth(2);
  h2->Draw("HIST SAME");

  h3->SetLineColor(c3);
  h3->SetLineWidth(2);
  h3->Draw("HIST SAME");
  legend->Draw("SAME");
}

void draw_set(TH1I* h1, TH1I* h2, TH1I* h3, TH1I* h4, TLegend* legend) {
  h1->SetMaximum(1.1*std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum(), h4->GetMaximum()}));
  
  h1->SetLineColor(c1);
  h1->SetLineWidth(2);
  h1->Draw("HIST");

  h2->SetLineColor(c2);
  h2->SetLineWidth(2);
  h2->Draw("HIST SAME");

  h3->SetLineColor(c3);
  h3->SetLineWidth(2);
  h3->Draw("HIST SAME");
  
  h4->SetLineColor(c5);
  h4->SetLineWidth(2);
  h4->Draw("HIST SAME");
  legend->Draw("SAME");
}

void generate_ntracks() {
  gStyle->SetOptStat(0);

  TChain chain("ntuple");
  setup_chain(chain);
  TTreeReader reader(&chain);
  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);

  // jet variables
  TTreeReaderArray<float> jet_pt
    // (reader, "TruthHSJet_pt");
    (reader, "AntiKt4EMTopoJets_pt");
  TTreeReaderArray<float> jet_eta
    // (reader, "TruthHSJet_eta");
    (reader, "AntiKt4EMTopoJets_eta");
  TTreeReaderArray<std::vector<int>> jet_track_indices
    (reader, "AntiKt4EMTopoJets_track_idx");

  // vertex variables
  /// truth vertex
  TTreeReaderArray<float> truth_vtx_z
    (reader, "TruthVtx_z");
  TTreeReaderArray<float> truth_vtx_time
    (reader, "TruthVtx_time");

  /// reco vertex
  TTreeReaderArray<float> reco_vtx_z
    (reader, "RecoVtx_z");
  TTreeReaderArray<float> reco_vtx_time
    (reader, "RecoVtx_time");
  TTreeReaderArray<float> reco_vtx_timeRes
    (reader, "RecoVtx_timeRes");
  TTreeReaderArray<int>   reco_vtx_valid
    (reader, "RecoVtx_hasValidTime");

  // track variables
  TTreeReaderArray<float> track_z0
    (reader, "Track_z0");
  TTreeReaderArray<float> track_pt
    (reader, "Track_pt");
  TTreeReaderArray<float> track_eta
    (reader, "Track_eta");
  TTreeReaderArray<float> track_time
    (reader, "Track_time");
  TTreeReaderArray<float> track_time_res
    (reader, "Track_timeRes");
  TTreeReaderArray<float> track_var_z0
    (reader, "Track_var_z0");
  TTreeReaderArray<int>   track_to_truthvtx
    (reader, "Track_truthVtx_idx");
  TTreeReaderArray<int>   track_to_particle
    (reader, "Track_truthPart_idx");
  TTreeReaderArray<int>   track_time_valid
    (reader, "Track_hasValidTime");

  // particle vars
  TTreeReaderArray<float> prod_vtx_z
    (reader, "TruthPart_prodVtx_z");

  int bins = 50;

  /// Various flavors of ntracks
  //// Number of tracks with valid time associated to HS vertex
  TH1I *nTracks1 = new TH1I("nTracks1", "nTracks valid time in Event,Forward Jet=0;nTracks;Entries", 80,-0.5,79.5);
  TH1I *nTracks2 = new TH1I("nTracks2", "nTracks valid time in Event,Forward Jet=1;nTracks;Entries", 80,-0.5,79.5);
  TH1I *nTracks3 = new TH1I("nTracks3", "nTracks valid time in Event,Forward Jet=2;nTracks;Entries", 80,-0.5,79.5);
  TH1I *nTracks4 = new TH1I("nTracks4", "nTracks valid time in Event,Forward Jet>2;nTracks;Entries", 80,-0.5,79.5);

  //// number of forward tracks with valid time associated to HS Vertex 
  TH1I *nForTracks1 = new TH1I("nForTracks1", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);
  TH1I *nForTracks2 = new TH1I("nForTracks2", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);
  TH1I *nForTracks3 = new TH1I("nForTracks3", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);
  TH1I *nForTracks4 = new TH1I("nForTracks4", "n Forward Tracks with valid time in HS Vertex;nTracks;Entries", 40,-0.5,39.5);

  //// ntracks with valid time, including their HS classifications
  ///// 0 forward Jet
  TH1I *nHSTracks1 = new TH1I("nHSTracks1", "Track Classification Forward Jet=0;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks1 = new TH1I("nNonHSTracks1", "", 50,-0.5,49.5);
  TH1I *nUnkHSTracks1 = new TH1I("nUnkHSTracks1", "", 50,-0.5,49.5);

  ///// 1 Forward Jet
  TH1I *nHSTracks2 = new TH1I("nHSTracks2", "Track Classification Forward Jet=1;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks2 = new TH1I("nNonHSTracks2", "", 50,-0.5,49.5);
  TH1I *nUnkHSTracks2 = new TH1I("nUnkHSTracks2", "", 50,-0.5,49.5);

  ///// 2 Forward Jet
  TH1I *nHSTracks3 = new TH1I("nHSTracks3", "Track Classification Forward Jet=2;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks3 = new TH1I("nNonHSTracks3", "", 50,-0.5,49.5);
  TH1I *nUnkHSTracks3 = new TH1I("nUnkHSTracks3", "", 50,-0.5,49.5);
  
  TH1I *nHSTracks4 = new TH1I("nHSTracks4", "Track Classification Forward Jet>2;nTracks;Entries", 50,-0.5,49.5);
  TH1I *nNonHSTracks4 = new TH1I("nNonHSTracks4", "", 50,-0.5,49.5);
  TH1I *nUnkHSTracks4 = new TH1I("nUnkHSTracks4", "", 50,-0.5,49.5);

  while (reader.Next()) {
    if (jet_pt.GetSize() < min_jets ) {
      if(debug)
	std::cout << "Skipping low jet event" << std::endl;
      continue;  // Skip if no jets are present
    }

    if(reco_vtx_valid[0] == 0) {
      if(debug)
	std::cout << "Skipping event where reco vertex has no valid time" << std::endl;
      continue;  // Skip if reco vertex has no valid time
    }
    
    // check reco HS vertex is with 2mm of truth HS vertex
    if(std::abs(truth_vtx_z[0] - reco_vtx_z[0]) > max_vtx_dz ) {
      if(debug)
	std::cout << "Skipping event due to incorrect HS vertex" << std::endl;
      continue;
    }
    
    // check if jet pt > 30 GeV
    float min_pt = *std::min_element(jet_pt.begin(), jet_pt.end());
    if(min_pt < min_jetpt) {
      if(debug)
	std::cout << "Skipping event due to low jet pt" << std::endl;
      continue;
    }

    int nForwardJet = 0;
    for(auto eta : jet_eta) {
      if (std::abs(eta) > min_abs_eta)
	nForwardJet++;
    }

    if (debug)
      std::cout << "nForwardJet = " << nForwardJet << std::endl;

    if (debug)
      std::cout << "nForwardJet = " << nForwardJet << std::endl;

    // ntracks calculation
    // counter for number of tracks with valid time within max_nsigma of HS Vertex
    int nvalidtime    = 0;
    // counter for number of forward tracks with valid time within max_nsigma of HS Vertex
    int nvalidftracks = 0;
    // counter for number of tracks with valid time (tagged as HS)
    int nTruthHS      = 0;
    // counter for number of tracks with valid time (tagged as assoc with a non-HS vertex)
    int nTruthNonHS   = 0;
    // counter for number of tracks with valid time (Not linked to a vertex)
    int nTruthUnk     = 0;
    // track loops
    for (int idx = 0; idx < track_z0.GetSize(); ++idx) {
      if(debug)
	std::cout << "track: " << idx << " hasvalidtime==" << track_time_valid[idx] << std::endl;

      if (!(track_time_valid[idx] == 1 && track_pt[idx] > 1.0)) // skip ALL low pt tracks with invalid time
	continue;
      
      float nsigma = (track_z0[idx] - reco_vtx_z[0])/std::sqrt(track_var_z0[idx]);

      if(std::abs(nsigma) > max_nsigma) {
	if(debug)
	  std::cout << "Skipping Track due to High n sigma" << std::endl;
	continue;
      }
      
      nvalidtime++;

      if(track_eta[idx] > min_abs_track_eta)
	nvalidftracks++;

      if(track_to_truthvtx[idx] == 0) {
	// associated to truth hs
	nTruthHS++;
      } else if (track_to_truthvtx[idx] == -1) {
	// unknown associated vertex (not used for vertex fitting)
	// In this case we should determine the closest truth vertex to the associated partciles production vertex
	if (track_to_particle[idx] != -1) {
	  double abs_dz=1000000;
	  int nearest_idx = -1000;
	  for(int vtxIndex = 0; vtxIndex < truth_vtx_z.GetSize(); vtxIndex++) {
	    auto abs_diff = std::abs(truth_vtx_z[vtxIndex] - prod_vtx_z[track_to_particle[idx]]);
	    if(abs_diff < abs_dz) {
	      abs_dz=abs_diff;
	      nearest_idx = vtxIndex;
	    }
	  }
	  if (nearest_idx == 0) {
	    nTruthHS++;
	  } else if (abs_dz < max_vtx_dz) { // we can reasonably say it is associated to another non HS vertex
	    nTruthNonHS++;
	  } else {
	    nTruthUnk++;
	  }
	} else {
	  nTruthUnk++;
	}
      } else {
	// associated to other primary vtx
	nTruthNonHS++;
      }
    }

    // exactly 0 forward jets
    if(nForwardJet == 0){
      nTracks1->Fill(nvalidtime);
      nForTracks1->Fill(nvalidftracks);
      nHSTracks1->Fill(nTruthHS);
      nNonHSTracks1->Fill(nTruthNonHS);
      nUnkHSTracks1->Fill(nTruthUnk);
    }
    // exactly 1 forward jet
    if(nForwardJet == 1){
      nTracks2->Fill(nvalidtime);
      nForTracks2->Fill(nvalidftracks);
      nHSTracks2->Fill(nTruthHS);
      nNonHSTracks2->Fill(nTruthNonHS);
      nUnkHSTracks2->Fill(nTruthUnk);
    }
    // exactly 2 forward jet
    if(nForwardJet == 2) {
      nTracks3->Fill(nvalidtime);
      nForTracks3->Fill(nvalidftracks);
      nHSTracks3->Fill(nTruthHS);
      nNonHSTracks3->Fill(nTruthNonHS);
      nUnkHSTracks3->Fill(nTruthUnk);
    }
    // > 2 forward jet
    if(nForwardJet > 2) {
      nTracks4->Fill(nvalidtime);
      nForTracks4->Fill(nvalidftracks);
      nHSTracks4->Fill(nTruthHS);
      nNonHSTracks4->Fill(nTruthNonHS);
      nUnkHSTracks4->Fill(nTruthUnk);
    }
  }

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);  // Align left-top
  
  // ntracks plot
  TLegend *nTracksLegend = new TLegend(0.65, 0.65, 0.9, 0.9);
  nTracksLegend->AddEntry(nTracks1, Form("=0 Forward Jet"), "l");
  nTracksLegend->AddEntry(nTracks2, Form("=1 Forward Jet"), "l");
  nTracksLegend->AddEntry(nTracks3, Form("=2 Forward Jet"), "l");
  nTracksLegend->AddEntry(nTracks4, Form(">2 Forward Jet"), "l");
  
  draw_set(nTracks1,nTracks2,nTracks3,nTracks4, nTracksLegend);
  canvas->Print("figs/ntracks.pdf(", "pdf");

  draw_set(nForTracks1,nForTracks2,nForTracks3,nForTracks4, nTracksLegend);
  canvas->Print("figs/ntracks.pdf", "pdf");

  TLegend *classLegend = new TLegend(0.65, 0.65, 0.9, 0.9);
  classLegend->AddEntry(nHSTracks1, Form("Truth HS"), "l");
  classLegend->AddEntry(nNonHSTracks1, Form("Truth Non-HS"), "l");
  classLegend->AddEntry(nUnkHSTracks1, Form("Truth Unknown"), "l");
  classLegend->Draw();

  draw_class(nHSTracks1, nNonHSTracks1, nUnkHSTracks1, classLegend);
  canvas->Print("figs/ntracks.pdf", "pdf");

  draw_class(nHSTracks2, nNonHSTracks2, nUnkHSTracks2, classLegend);
  canvas->Print("figs/ntracks.pdf", "pdf");

  draw_class(nHSTracks3, nNonHSTracks3, nUnkHSTracks3, classLegend);
  canvas->Print("figs/ntracks.pdf", "pdf");

  draw_class(nHSTracks4, nNonHSTracks4, nUnkHSTracks4, classLegend);
  canvas->Print("figs/ntracks.pdf)", "pdf");

}
