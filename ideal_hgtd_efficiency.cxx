#include <Rtypes.h>
#include <TRandom1.h>
#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLatex.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGaxis.h>

#include <boost/filesystem.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>

#include "./my_utilities.h"

#define debug false

using namespace myutl;

double getSmearedTrackTime(int idx, double res, BranchPointerWrapper *branch) {
  // know that it has a valid time in HGTD and has a valid truth link
  
  int particle_idx = branch->track_to_particle[idx]; // this is anything
  int truthvtx_idx = branch->track_to_truthvtx[idx]; // this is anything
  if (debug) std::cout << "Particle_idx: " << particle_idx << std::endl;
  if (debug) std::cout << "TruthVtx_idx: " << truthvtx_idx << std::endl;
  double t_part, reso = res;
  if (particle_idx == -1) { // no particle link
    if (truthvtx_idx != -1) { // some valid truth link, smear that vertex time
      if (debug) std::cout << "PATH A" << std::endl;
      t_part = branch->truth_vtx_time[truthvtx_idx];
    } else { // some random pileup, assume
      if (debug) std::cout << "PATH B" << std::endl;
      t_part = gRandom->Gaus(branch->truth_vtx_time[0],175);
    }
  } else {
    if (debug) std::cout << "PATH C" << std::endl;
    t_part = branch->particle_t[particle_idx];
  }
  if (debug) std::cout << "DONE" << std::endl;
  return gRandom->Gaus(t_part,reso);
}

std::vector<Cluster> doConeClustering(std::vector<int> indices, BranchPointerWrapper *branch) {
  std::vector<bool> checked(branch->track_pt.GetSize(), false);
  std::vector<Cluster> clusters;
  double fixed_res = 10.0;
  
  if (debug)
    std::cout << "Got: " << indices.size() << std::endl;

  std::vector<int> timefiltered_indices;
  std::map<int,double> smeared_times, smeared_resos;
  for (auto idx: indices) {
    timefiltered_indices.push_back(idx);
    smeared_times.emplace(idx, getSmearedTrackTime(idx, fixed_res, branch));
  }

  while (true) {
    // Step 1: Find the highest-pt unchecked track
    int max_pt_idx = -1;
    double max_pt = -1.0;

    for (int idx : timefiltered_indices) {
      if (checked[idx]) continue;
      double this_pt = branch->track_pt[idx];
      if (this_pt > max_pt) {
        max_pt = this_pt;
        max_pt_idx = idx;
      }
    }

    if (max_pt_idx == -1) break; // all tracks have been checked

    // Step 2: Build cluster around max_pt_idx
    checked[max_pt_idx] = true;

    double t1 = smeared_times[max_pt_idx];
    double s1 = fixed_res;

    std::vector<int> to_merge = {max_pt_idx};

    for (int idx : timefiltered_indices) {
      if (checked[idx]) continue;

      double t2 = smeared_times[idx];
      double s2 = fixed_res;
      double nsigma = std::abs(t1 - t2) / std::sqrt(s1 * s1 + s2 * s2);

      if (nsigma < 3.0) {
        to_merge.push_back(idx);
        checked[idx] = true;
      }
    }

    // Step 3: Form cluster from selected indices
    std::vector<Cluster> simpleclusters;
    for (int idx : to_merge) {
      
      double trk_time = smeared_times[idx];
      double trk_pt = branch->track_pt[idx];
      double res = fixed_res;
      double significance = std::abs(branch->track_z0[idx] - branch->reco_vtx_z[0]) / std::sqrt(branch->track_var_z0[idx]);
      bool ismaxpt = idx == max_pt_idx;
      
      std::map<ScoreType,double> scores;
      scores[ScoreType::HGTD] = 0; // do not use this
      scores[ScoreType::SIG] = exp(-significance);
      scores[ScoreType::MAXPT] = 0;
      scores[ScoreType::TRKPT] = trk_pt;
      scores[ScoreType::SUMPT2] = trk_pt*trk_pt;
      scores[ScoreType::SIGTRKPT] = trk_pt*exp(-significance);
      scores[ScoreType::IDEAL] = trk_time;

      Cluster cluster = {{trk_time}, {res}, {trk_time}, {idx}, scores};
      cluster.max_pt_cluster = ismaxpt;
      if (not simpleclusters.empty()) {
	Cluster newCluster = mergeClusters(simpleclusters[0],cluster);
	simpleclusters.at(0) = newCluster;
      } else {
	simpleclusters.push_back(cluster);
      }
    }
    clusters.push_back(simpleclusters.at(0));
  }

  for (Cluster& cluster: clusters)
    cluster.purity = calc_pt_purity(cluster, branch);

  return clusters;
}

std::vector<Cluster> clusterTracksInTime(std::vector<int> track_indices, double cluster_distance, BranchPointerWrapper *branch) {
  if (track_indices.empty() && debug)
    std::cout << "EMPTY!!!!" << std::endl;
  std::vector<Cluster> collection;
  collection = doConeClustering(track_indices, branch);

  for (Cluster& cluster: collection)
    cluster.purity = calc_pt_purity(cluster, branch);

  return collection;
}

void ideal_hgtd_efficiency() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);
  // gPad->SetLeftMargin(0.1);

  double diff_min = -1000.0, diff_max = 1000.0;
  double diff_width = 2.0;

  double purity_min = 0, purity_max = 1;
  double purity_width = 0.05;

  double fjet_min = 1, fjet_max = 31, fold_fjet = 5;
  double fjet_width = 1.0;

  double track_min = 0, track_max = 72, fold_track = 20;
  double track_width = 2.0;

  double pu_track_min = track_min, pu_track_max = track_max, fold_hs_track = fold_track;
  double pu_track_width = track_width;

  double hs_track_min = track_min, hs_track_max = track_max, fold_pu_track = fold_track;
  double hs_track_width = track_width;

  double pu_frac_min = 0, pu_frac_max = 1 , fold_pu_frac = pu_frac_max;
  double pu_frac_width = 0.05;

  double z_min = -200, z_max = 200, fold_z = 100;
  double z_width = 10.0;

  // All Needed Histogramming information
  std::map<ScoreType,TH1D*> inclusive_resos;
  std::map<ScoreType,TH1D*> inclusive_purity;
  for (ScoreType score: enum_vec) {
    inclusive_resos[score] = new TH1D(Form("inclusivereso_%s",toString(score)),
				      Form("Inclusive %s t_{0} - TruthVtx t_{0} (Ideal Eff. HGTD);#Delta t[ps];Entries",toString(score)),
				      (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    inclusive_purity[score] = new TH1D(Form("inclusivepurity_%s",toString(score)),
				       Form("Inclusive %s Purity (Ideal Eff. HGTD);Purity;Entries",toString(score)),
				       (int)((purity_max-purity_min)/purity_width), purity_min, purity_max);
  }

  PlotObj fjet     ("n Forward Jets", "fjet", "figs/idealeff_nfjet.pdf",
		    fjet_min     , fjet_max    , fjet_width    ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_fjet    , fold_fjet+fjet_width        );
  
  PlotObj ftrack   ("n Forward Tracks", "track", "figs/idealeff_ntrack.pdf",
		    track_min    , track_max   , track_width   ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_track   , fold_track+track_width      );
  
  PlotObj pu_frac  ("Pile Up Fraction", "pu_frac", "figs/idealeff_pufrac.pdf",
		    pu_frac_min  , pu_frac_max , pu_frac_width ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_pu_frac , pu_frac_max                 );
  
  PlotObj hs_track ("n Forward HS Tracks", "hs_track", "figs/idealeff_nhstrack.pdf",
		    hs_track_min , hs_track_max, hs_track_width,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_hs_track, fold_hs_track+hs_track_width);
  
  PlotObj pu_track ("n Forward PU Tracks", "pu_track", "figs/idealeff_nputrack.pdf",
		    pu_track_min , pu_track_max, pu_track_width,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    fold_pu_track, fold_pu_track+pu_track_width);
  
  PlotObj recovtx_z("(Reco Vtx z) (mm)", "z", "figs/idealeff_vtx_z.pdf",
		    z_min        , z_max       , z_width       ,
		    diff_min     , diff_max    , diff_width    ,
		    purity_min   , purity_max  , purity_width  ,
		    z_max        , z_max                       );

  TH2D* hs_pu_inclusive = new TH2D("hs_vs_pu","N HardScatter vs. N Pile-Up;n HS Tracks;n PU Tracks",
				   (int)((hs_track_max-hs_track_min)/hs_track_width), hs_track_min, hs_track_max,
				   (int)((pu_track_max-pu_track_min)/pu_track_width), pu_track_min, pu_track_max);

  std::cout << "Starting Event Loop" << std::endl;
  bool progress = true;
  int n_incorrect = 0, n_chosen = 0;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed

    auto readnum = chain.GetReadEntry()+1;
    if (progress and readnum % 1000 == 0)
      std::cout << "Progress: " << readnum << "/" << chain.GetEntries() << "\n";

    if (not branch.pass_basic_cuts()) // check n_jet cut and vertex dz cut
      continue;

    // check if there is one forward jet with pt > 30 GeV
    int nForwardJet;
    branch.count_forward_jets(nForwardJet);

    if (not branch.pass_jet_pt_cut()) {
      if (debug) std::cout << "Skipping pt cut event" << std::endl;
      continue;
    }

    std::vector<int> tracks = getAssociatedTracks(&branch, min_track_pt);

    int nForwardTrack, nForwardTrack_HS, nForwardTrack_PU;
    branch.count_forward_tracks(nForwardJet,nForwardTrack_HS,nForwardTrack_PU,tracks);

    hs_pu_inclusive->Fill(nForwardTrack_HS,nForwardTrack_PU);
    
    double pu_ratio = (double)nForwardTrack_PU / (double)nForwardTrack;
    double reco_z = branch.reco_vtx_z[0];

    auto eff_fill_val_fjet     = folded(nForwardJet     , (int)fold_fjet) ;
    auto eff_fill_val_track    = folded(nForwardTrack   , (int)fold_track);
    auto eff_fill_val_hs_track = folded(nForwardTrack_HS, (int)fold_hs_track);
    auto eff_fill_val_pu_track = folded(nForwardTrack_PU, (int)fold_pu_track);
    auto eff_fill_val_z        = reco_z;
    auto eff_fill_val_pu_ratio = pu_ratio;

    if (debug) {
      std::cout << "nForwardJet = " << nForwardJet << std::endl;
      std::cout << "nForwardTrack = " << nForwardTrack << std::endl;
      std::cout << "Vertex_z = " << reco_z << std::endl;
    }

    std::vector<Cluster> cluster = clusterTracksInTime(tracks, 3.0, &branch);

    fjet.eff_total->     Fill(eff_fill_val_fjet    );
    ftrack.eff_total->   Fill(eff_fill_val_track   );
    pu_frac.eff_total->  Fill(eff_fill_val_pu_ratio);
    hs_track.eff_total-> Fill(eff_fill_val_hs_track);
    pu_track.eff_total-> Fill(eff_fill_val_pu_track);
    recovtx_z.eff_total->Fill(eff_fill_val_z       );
    
    // std::cout << "Acessing cluster size" << std::endl;
    if (cluster.size() == 0)
      continue;

    std::map<ScoreType,Cluster> chosen = chooseCluster(cluster, &branch);

    // run HGTD Clustering
    // std::vector<Cluster> hgtd_clusters = clusterTracksInTime(tracks, 3.0, &branch, false);

    // if (branch.reco_vtx_valid[0] == 1)
    //   chosen[ScoreType::HGTD] = chooseHGTDCluster(hgtd_clusters, &branch);

    for (ScoreType score : enum_vec) {
      
      if (branch.reco_vtx_valid[0] == 0)
	continue;
      
      if (fjet.eff_pass.count(score) != 1)
	continue;
      
      if (score == ScoreType::INVALID)
	continue;
      
      Cluster scored = chosen[score];
      double score_based_time = scored.values.at(0);
      double cluster_purity = scored.purity;
      double diff = score_based_time - branch.truth_vtx_time[0];
      inclusive_resos.at(score)->Fill(diff);
      inclusive_purity.at(score)->Fill(cluster_purity);

      // if (score == TRKPT and branch.reco_vtx_valid[0] == 1 and
      //     diff < std::abs(branch.reco_vtx_time[0] - branch.truth_vtx_time[0]) and
      //     60 < diff)
        // std::cout << Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0 --extra_time %.2f\n", filename.substr(40,6).c_str(),this_evnt,score_based_time);
      
      if (passEfficiency(diff, scored, &branch)) {
	fjet.eff_pass.at(score)->     Fill(eff_fill_val_fjet    );
	ftrack.eff_pass.at(score)->   Fill(eff_fill_val_track   );
	pu_frac.eff_pass.at(score)->  Fill(eff_fill_val_pu_ratio);
	hs_track.eff_pass.at(score)-> Fill(eff_fill_val_hs_track);
	pu_track.eff_pass.at(score)-> Fill(eff_fill_val_pu_track);
	recovtx_z.eff_pass.at(score)->Fill(eff_fill_val_z       );
      } else if (score == TRKPT) {
	n_chosen++;
      }

      // fill diff hists
      fjet.hist.at(score)->     Fill(nForwardJet     , diff);
      ftrack.hist.at(score)->   Fill(nForwardTrack   , diff);
      pu_frac.hist.at(score)->  Fill(pu_ratio        , diff);
      hs_track.hist.at(score)-> Fill(nForwardTrack_HS, diff);
      pu_track.hist.at(score)-> Fill(nForwardTrack_PU, diff);
      recovtx_z.hist.at(score)->Fill(reco_z          , diff);

      // fill purities
      fjet.purity.at(score)->     Fill(nForwardJet     , cluster_purity);
      ftrack.purity.at(score)->   Fill(nForwardTrack   , cluster_purity);
      pu_frac.purity.at(score)->  Fill(pu_ratio        , cluster_purity);
      hs_track.purity.at(score)-> Fill(nForwardTrack_HS, cluster_purity);
      pu_track.purity.at(score)-> Fill(nForwardTrack_PU, cluster_purity);
      recovtx_z.purity.at(score)->Fill(reco_z          , cluster_purity);
    }
  }
  
  // add all objects to plotobjs before drawing everything
  std::vector<PlotObj> plots = {fjet, ftrack, hs_track, pu_track, pu_frac, recovtx_z};
  for (int i = 0; i < plots.size(); i++) {
    
  }

  std::cout << "FINISHED CREATING " << std::endl;
  
  for (auto plot: plots) {
    plot.plot_logic(canvas);
  }  

  plot_inclusive("figs/idealeff_resos_logscale.pdf", true, true, -400, 400,
		 canvas, inclusive_resos);

  plot_inclusive("figs/idealeff_resos_linscale.pdf", false, true, -150, 150,
		 canvas, inclusive_resos);

  auto purity_fname = "figs/idealeff_purity.pdf";
  canvas->Print(Form("%s[",purity_fname));
  inclusive_purity.at(ScoreType::TRKPT)->Draw("HIST");
  double maxval = 0.0;
  for (auto pair: inclusive_purity) {
    if (pair.first == ScoreType::HGTD) continue;
    TH1D *hist = pair.second;
    hist->SetLineColor(colors[pair.first % colors.size()]);
    hist->SetLineWidth(2);
    hist->SetBinContent(hist->GetNbinsX(), 
			hist->GetBinContent(hist->GetNbinsX()) + hist->GetBinContent(hist->GetNbinsX()+1));
    hist->SetBinContent(hist->GetNbinsX()+1, 0); // clear overflow
    hist->Scale(1/hist->Integral());
    double thismax = hist->GetMaximum();
    if (thismax > maxval)
      maxval = thismax;
  }
  for (auto pair: inclusive_purity) {
    if (pair.first == ScoreType::HGTD) continue;
    TH1D *hist = pair.second;
    hist->SetMaximum(maxval/(1-(fjet.algolegend->GetY2()-fjet.algolegend->GetY1())));
    hist->SetMinimum(0.00);
    hist->Draw("HIST SAME");
  }
  fjet.algolegend->Draw("SAME");
  canvas->Print(Form("%s",purity_fname));
  canvas->Print(Form("%s]",purity_fname));

  auto hs_pu_fname = "figs/hgtdtimes_hs_v_pu.pdf";
  canvas->Print(Form("%s[",hs_pu_fname));
  hs_pu_inclusive->GetXaxis()->SetRangeUser(0, 35);
  hs_pu_inclusive->GetYaxis()->SetRangeUser(0, 50);
  hs_pu_inclusive->Draw("COLZ");
  canvas->Print(Form("%s",hs_pu_fname));
  canvas->Print(Form("%s]",hs_pu_fname));

  // std::cout << "N incorrect = " << n_incorrect << std::endl;
  // std::cout << "N BAD = " << n_chosen << std::endl;
}
