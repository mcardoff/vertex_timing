#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLine.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <boost/filesystem.hpp>

// Explicitly load the dictionary for ROOT to recognize nested vectors
#ifdef __CLING__
#pragma link C++ class std::vector<std::vector<int>>+;
#endif

#define debug false

TTree* setup_tree(std::string Number) {
  TFile *file = TFile::Open(Form("./ntuple/user.scheong.42871997.Output._%s.SuperNtuple.root", Number.c_str()));
  TTree *tree = (TTree*)file->Get("ntuple");
  return tree;
}

void generate_onefile_timeplot() {
  gStyle->SetOptStat(0);
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  // colors
  auto c1 = TColor::GetColor("#3f90da");
  auto c2 = TColor::GetColor("#ffa90e");
  auto c3 = TColor::GetColor("#bd1f01");
  auto c4 = TColor::GetColor("#94a4a2");
  auto c5 = TColor::GetColor("#832db6");

  auto Numbers = {"000031", "000032", "000033", "000034", "000036", "000037",
		  "000045", "000046", "001445", "001446", "001447", "001449",
		  "001450", "001451", "001452", "001528", "001537", "001538",
		  "001550", "001975", "001982", "002016", "002017", "002029"};
  std::vector<std::string> eventDisplayCommands;
  for(auto Number: Numbers){
    TCanvas *canvas = new TCanvas("canvas", "Track Time Histogram", 800, 600);
    TTree * tree = setup_tree(Number);
  
    TTreeReader reader(tree);
    // jet variables
    TTreeReaderArray<float> jet_pt
      (reader, "TruthHSJet_pt");
      // (reader, "AntiKt4EMTopoJets_pt");
    TTreeReaderArray<float> jet_eta
      (reader, "TruthHSJet_eta");
      // (reader, "AntiKt4EMTopoJets_eta");
    TTreeReaderArray<std::vector<int>> jet_track_indices
      (reader, "AntiKt4EMTopoJets_track_idx");

    // vertex variable
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
    TTreeReaderArray<float> reco_vtx_time_res
      (reader, "RecoVtx_timeRes");
    TTreeReaderArray<int>   reco_vtx_valid
      (reader, "RecoVtx_hasValidTime");

    // track Variables
    TTreeReaderArray<float> track_z0
      (reader, "Track_z0");
    TTreeReaderArray<float> track_pt
      (reader, "Track_pt");
    TTreeReaderArray<float> track_eta
      (reader, "Track_eta");
    TTreeReaderArray<float> track_time
      (reader, "Track_time");
    TTreeReaderArray<float> track_var_z0
      (reader, "Track_var_z0");
    TTreeReaderArray<int>   track_time_valid
      (reader, "Track_hasValidTime");

    // int num_events = (int)tree->GetEntries();
    int numpages = 1;
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_0fjet.pdf("  ,Number), "pdf");
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_1fjet.pdf("  ,Number), "pdf");
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_2fjet.pdf("  ,Number), "pdf");
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_2pfjet.pdf(" ,Number), "pdf");
    // std::vector<int> nHSTracksValidTime;

    // cut variables
    int    min_jets     = 2;    // min number of jets
    double min_jetpt    = 30.0; // self explanatory
    double min_abs_eta  = 2.0;  // min eta for a vtx to be considered "forward"
    double min_abs_eta_track = 2.4;  // min eta for a track to be considered "forward"
    double min_track_pt = 1.0;  // track_pt > 1.0 GeV
    double max_vtx_dz   = 2.0;  // max error for reco HS vertex z
    double max_nsigma   = 2.5;  // how close a track can be to PV

    while (reader.Next()) {
      if (jet_pt.GetSize() < min_jets) {
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
      std::cout << "Min_pt: " << min_pt << std::endl;
      if(min_pt < min_jetpt) {
	if(debug)
	  std::cout << "Skipping event due to low jet pt" << std::endl;
	continue;
      }

      float truth_vertex_time = truth_vtx_time[0];
      float reco_vertex_time = reco_vtx_time[0];
      float reco_vertex_time_res = reco_vtx_time_res[0];

      int nForwardJet = 0;
      for(auto eta : jet_eta) {
	if (std::abs(eta) - min_abs_eta > 1e-6)
	  nForwardJet++;
      }
    
      // select valid times
      std::vector<float> selected_times;
      int nvalidtime = 0;
      for (int idx = 0; idx < track_z0.GetSize(); ++idx) {
	if(debug)
	  std::cout << "track: " << idx << " hasvalidtime==" << track_time_valid[idx] << std::endl;

	if (!(track_time_valid[idx] == 1 && track_pt[idx] > min_track_pt))
	  continue; // skip low pt tracks

	float nsigma = (track_z0[idx] - reco_vtx_z[0])/std::sqrt(track_var_z0[idx]+reco_vertex_time_res);
	if(std::abs(nsigma) > max_nsigma) {
	  if(debug)
	    std::cout << "Skipping track " << idx << " due to high nsigma" << std::endl;
	  continue; // skip tracks too far from vertex
	}

	selected_times.push_back(track_time[idx]);
	if(track_eta[idx] > min_abs_eta_track)
	  nvalidtime++; // nvalidtime is now number of forward tracks
     
      }

      // nHSTracksValidTime.push_back(nvalidtime);

      if (!selected_times.empty()) {
	float min_time = std::min({*std::min_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
	float max_time = std::max({*std::max_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
      
	float range = max_time - min_time;
	float extended_min_time = min_time - 0.05 * range;
	float extended_max_time = max_time + 0.05 * range;
      
	TH1F *hist = new TH1F("track_time_hist_event",
			      Form("Track Time Distribution Event %lld", reader.GetCurrentEntry()), 50, extended_min_time, extended_max_time);
	hist->GetXaxis()->SetTitle("Time (ps)");
	hist->GetYaxis()->SetTitle("A.U.");
	for (float time : selected_times)
	  hist->Fill(time);
      
	hist->SetFillColorAlpha(kBlue, 0.35);
	// hist->Scale(1/hist->Integral());
	hist->Draw("HIST F");

	TLine *truthHSLine = new TLine(truth_vertex_time, 0, truth_vertex_time, hist->GetMaximum());
	truthHSLine->SetLineColor(kRed);
	truthHSLine->SetLineWidth(2);
	truthHSLine->Draw("SAME");

	TLine *recoHSLine = new TLine(reco_vertex_time, 0, reco_vertex_time, hist->GetMaximum());
	recoHSLine->SetLineColor(kGreen);
	recoHSLine->SetLineWidth(2);
	recoHSLine->Draw("SAME");

	TLatex *latex = new TLatex();
	latex->SetNDC();  // Use Normalized Device Coordinates (0 to 1)
	latex->SetTextSize(0.03);
	latex->SetTextAlign(13);  // Align text left
	latex->DrawLatex(0.15, 0.88,
			 Form("Truth Vertex Time: %.2f ps, %d fjet(s), %d Forward Tracks",
			      truth_vertex_time, nForwardJet, nvalidtime));
	latex->DrawLatex(0.15, 0.85,
			 Form("#Deltat: %.2f ps, #Deltat/res: %.2f", truth_vertex_time-reco_vertex_time, (truth_vertex_time-reco_vertex_time)/reco_vertex_time_res));
	latex->DrawLatex(0.15, 0.82, Form("Jet 0 #eta, p^{T}: %.2f, %.2f GeV", jet_eta[0],jet_pt[0]));
	if(jet_eta.GetSize() > 1) {
	  latex->DrawLatex(0.15, 0.79, Form("Jet 1 #eta, p^{T}: %.2f, %.2f GeV", jet_eta[1],jet_pt[1]));
	  if(jet_eta.GetSize() > 2)
	    latex->DrawLatex(0.15, 0.76, Form("Jet 2 #eta, p^{T}: %.2f, %.2f GeV", jet_eta[2],jet_pt[2]));
	}
      
	hist->SetMaximum(1.25*hist->GetMaximum());
	// separate into nforjet regions
	if (nForwardJet == 0)
	  canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_0fjet.pdf",Number), "pdf");
	if (nForwardJet == 1) 
	  canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_1fjet.pdf",Number), "pdf");
	if (nForwardJet == 2) 
	  canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_2fjet.pdf",Number), "pdf");
	if (nForwardJet >= 2) 
	  canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_2pfjet.pdf",Number), "pdf");

	if (nForwardJet == 0)
	  eventDisplayCommands.push_back(Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0",Number,reader.GetCurrentEntry()));
	
	delete hist;
	delete latex;
	delete truthHSLine;
	delete recoHSLine;
      }
    }

    // may result in dupes
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_0fjet.pdf)",Number), "pdf");
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_1fjet.pdf)",Number), "pdf");
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_2fjet.pdf)",Number), "pdf");
    canvas->Print(Form("figs/trackhists/%s/truthjettrackhists_2pfjet.pdf)",Number), "pdf");

    // int minN = *std::min_element(nHSTracksValidTime.begin(), nHSTracksValidTime.end());
    // int maxN = *std::max_element(nHSTracksValidTime.begin(), nHSTracksValidTime.end());
    // int rangeN = maxN - minN;
  
    // TH1F *newHist = new TH1F("newHist", "Number of tracks in HS Vertex with valid time", rangeN, minN-0.5, maxN-0.5);
    // for(auto n : nHSTracksValidTime)
    //   newHist->Fill(n);
    // newHist->SetLineColor(kBlue);
    // newHist->SetFillColorAlpha(kBlue, 0.35);
    // newHist->SetLineWidth(2);
    // newHist->GetXaxis()->SetTitle("nTracks");
    // newHist->GetYaxis()->SetTitle("Entries");
  
    // newHist->Draw("HIST F");

    // double mean = newHist->GetMean();

    // TLine *meanLine = new TLine(mean, 0, mean, newHist->GetMaximum());
    // meanLine->SetLineColor(kRed);
    // meanLine->SetLineWidth(2);
    // meanLine->Draw("SAME");
    // TLatex *latex = new TLatex();
    // latex->SetNDC();  // Use Normalized Device Coordinates (0 to 1)
    // latex->SetTextSize(0.03);
    // latex->SetTextAlign(13);  // Align text left
    // latex->DrawLatex(0.68, 0.88, Form("Avg: %.2f", newHist->GetMean()));
    // canvas->SaveAs("nTracks.pdf");
  }
  for(auto e: eventDisplayCommands) {
    std::cout << e << std::endl;
  }
}
