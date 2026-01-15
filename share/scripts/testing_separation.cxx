#include "common_includes.h"
#include <TLegend.h>

using namespace myutl;

void testing_separation() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple-hgtd");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  double s_min = -5, s_max = 5, s_width = 0.1;
  int s_nbins = (int)((s_max-s_min)/(s_width));

  // Pileup histogram, only tracks where isHS[track_to_truthvtx] == false
  TH1D * s_pu = new TH1D("s_pu", "Significance of PU Tracks Associated to PV;s(z);Entries",
			 s_nbins, s_min, s_max);
  
  // Hard Scatter histogram, only tracks where isHS[track_to_truthvtx] == false
  TH1D * s_hs = new TH1D("s_hs", "Significance of HS Tracks Associated to PV;s(z);Entries",
			 s_nbins, s_min, s_max);

  bool check_valid_times = true;
  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  while (reader.Next()) {
    if (not branch.pass_basic_cuts()) continue;

    if (not branch.pass_jet_pt_cut()) continue;
    
    std::vector<int> tracks = getAssociatedTracks(&branch, min_track_pt);

    // int nForwardTrack=0, nForwardTrack_HS=0, nForwardTrack_PU=0;
    // branch.count_forward_tracks(nForwardTrack,nForwardTrack_HS,nForwardTrack_PU,tracks, check_valid_times);

    // if (not branch.pass_forward_hs_tracks(nForwardTrack_HS)) continue;
    for (const auto i_trk: tracks) {
      // find the nearest vertex in z, if the track is PU, store that in s_pu,
      double z0 = branch.track_z0[i_trk], var_z0 = branch.track_var_z0[i_trk];
      int vtx = branch.track_to_truthvtx[i_trk];
      bool isHS = vtx != -1 ? branch.truth_vtx_ishs[vtx] : false;

      // idx of nearest RECO vtx
      int nearest_vtx = 0;
      double dzNearest = z0-branch.reco_vtx_z[nearest_vtx];
      double absdzNearest = std::abs(dzNearest);
      for (int i_vtx = 0; i_vtx < branch.reco_vtx_z.GetSize(); ++i_vtx) {
	double dz = z0-branch.reco_vtx_z[i_vtx];
	double absdz = std::abs(dz);
	if (absdz < absdzNearest) {
	  nearest_vtx = i_vtx;
	  dzNearest = dz;
	  absdzNearest = absdz;
	}
      }
      double s_vtx_trk = dzNearest / std::sqrt(var_z0);
      // double s_vtx_trk = branch.track_near_sig[i_trk];
      if (isHS)
	s_hs->Fill(s_vtx_trk);
      else
	s_pu->Fill(s_vtx_trk);
    }
  }

  const char* fname = "figs/sigcomp.pdf";
  canvas->Print(Form("%s[",fname));

  // set colors/aesthetics
  s_hs->SetLineColor(c01);
  s_pu->SetLineColor(c02);
  s_hs->SetLineWidth(2);
  s_pu->SetLineWidth(2);
  s_pu->Scale(1/s_pu->Integral());
  s_hs->Scale(1/s_hs->Integral());
  
  // draw only HS
  s_hs->Draw("HIST");
  canvas->Print(fname);
  // draw only PU
  s_pu->Draw("HIST");
  canvas->Print(fname);
  // draw both overlayed with new title
  s_hs->SetMaximum(1.1*s_pu->GetMaximum());
  s_hs->Draw("HIST");
  s_pu->Draw("HIST SAME");
  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(s_hs, "Nearest is HS");
  legend->AddEntry(s_pu, "Nearest is PU");
  legend->Draw("SAME");
  canvas->Print(fname);
  canvas->Print(Form("%s]",fname));
}
