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

void test_track_particle_association() {
  gStyle->SetOptStat(0);
  // colors
  auto c1 = TColor::GetColor("#3f90da");
  auto c2 = TColor::GetColor("#ffa90e");
  auto c3 = TColor::GetColor("#bd1f01");
  auto c4 = TColor::GetColor("#94a4a2");
  auto c5 = TColor::GetColor("#832db6");
  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  TChain chain("ntuple");

  // ttbar sample
  for (const auto& entry : boost::filesystem::directory_iterator("./ntuple")) {
    // if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
    chain.Add(entry.path().c_str());
    break; // add 1
  }

  TTreeReader reader(&chain);
  // gives truth particle index from track index
  TTreeReaderArray<int> track_to_particle(reader, "Track_truthPart_idx");
  TTreeReaderArray<int> track_to_truthvtx(reader, "Track_truthVtx_idx");
  TTreeReaderArray<float> truth_vertex_z (reader, "TruthVtx_z");
  TTreeReaderArray<float> prod_vertex_z  (reader, "TruthPart_prodVtx_z");

  TH1F *test1 = new TH1F("test1", "Particle ProdVtx_{z}-TruthVtx_{z};#Deltaz(mm);Entries", 20, -0.5, 19.5);

  // first find truth particle associated 
  while(reader.Next()) {
    // loop over partiles
    for(int i = 0; i<track_to_particle.GetSize(); ++i) {
      // i is track index
      int particle_index = track_to_particle[i];
      int truthVtx_index = track_to_truthvtx[i];
      // before anything, see what dz is

      //find closest truth vertex to prodvtx
      if (particle_index != -1) {
	double abs_dz=1000000;
	int nearest_idx = -1000;
	for(int vtxIndex = 0; vtxIndex < truth_vertex_z.GetSize(); vtxIndex++) {
	  auto abs_diff = std::abs(truth_vertex_z[vtxIndex] - prod_vertex_z[particle_index]);
	  if(abs_diff < abs_dz) {
	    abs_dz=abs_diff;
	    nearest_idx = vtxIndex;
	  }
	}
	if(truthVtx_index == -1){
	  std::cout << "----------------------------------------" << std::endl;
	  std::cout << "TruthVtx Index == " << truthVtx_index << std::endl;
	  std::cout << "TruthParticle Index == " << particle_index << std::endl;
	  std::cout << "Closest TruthVtx to ProdVtx for Particle " <<  nearest_idx << std::endl;
	  std::cout << "----------------------------------------" << std::endl;
	  test1->Fill(nearest_idx);
	}
      }
    }
  }

  test1->SetLineColor(c1);
  test1->SetLineWidth(2);
  test1->Draw("HIST");
  canvas->Print("test.pdf");
  // test1->Show();
}
