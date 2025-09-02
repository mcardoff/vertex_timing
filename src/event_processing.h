#ifndef EVENT_PROCESSING_H
#define EVENT_PROCESSING_H

#include "clustering_includes.h"
#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "plotting_utilities.h"

using boost::filesystem::directory_iterator;

namespace myutl {

  static void setup_chain(
    TChain &chain, const char* ntuple_dir
  ) {
    // VBF H->Invisible sample
    for (const auto& entry : directory_iterator(ntuple_dir)) {
      if (entry.is_directory()) continue;
      if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
      chain.Add(entry.path().c_str());
      // break;
    }
  
    if (chain.GetEntries() == 0) {
      std::cerr << "No ROOT files found in directory: " << std::endl;
      return;
    }
  }

  static void setup_chain(
    TChain &chain, std::string Number
  ) {
    chain.Add(Form("../ntuple-hgtd/user.mcardiff.45809429.Output._%s.SuperNtuple.root", Number.c_str()));
  }

  static bool passTrackVertexAssociation(
    int track_idx, int vertex_idx, 
    BranchPointerWrapper *branch,
    double significance_cut
  ) {
    double
      trk_z0    = branch->track_z0[track_idx],
      trk_z_var = branch->track_var_z0[track_idx],
      vx_z      = branch->reco_vtx_z[vertex_idx],
      vx_z_var  = 0.0;
    
    double nsigma_prim = std::abs(trk_z0 - vx_z) / std::sqrt(trk_z_var + vx_z_var);
    return nsigma_prim < significance_cut;
  }

  static std::vector<int> pileupRemoval(
    std::vector<int> tracks, 
    BranchPointerWrapper *branch,
    double significance_cut
  ) {
    std::vector<int> output;
    for (const auto& trk: tracks) {
       int near = (int)branch->track_near_idx[trk];
       double near_z0 = branch->track_near_z0sin[trk] / std::sin(branch->track_theta[trk]);
       double near_var_z0 = (pow(branch->track_near_z0sin_unc[trk],2)-pow(near_z0*std::cos(branch->track_theta[trk])*branch->track_var_theta[trk],2))/std::sin(branch->track_theta[trk])*std::sin(branch->track_theta[trk]);

       double sig_near = std::abs(near_z0 - branch->reco_vtx_z[near]) / std::sqrt(near_var_z0);
       double sig_prim = std::abs(branch->track_z0[trk] - branch->reco_vtx_z[0]) / std::sqrt(branch->track_var_z0[trk]);

       // if (near != 0 && sig_near < sig_prim && sig_near < significance_cut)
       if (near != 0 && sig_near < significance_cut)
	 continue;

       output.push_back(trk);
    }
    return output;
  }

  static std::vector<int> getAssociatedTracks(
    BranchPointerWrapper *branch,
    double min_trk_pt, double max_trk_pt
  ) {
    std::vector<int> good_tracks;

    for (int trk = 0; trk < branch->track_z0.GetSize(); ++trk) {
      double
	trk_eta = branch->track_eta[trk],
	trk_pt  = branch->track_pt[trk],
	trk_quality = branch->track_quality[trk];

      if (std::abs(trk_eta) < min_hgtd_eta or
	  std::abs(trk_eta) > max_hgtd_eta)
        continue;

      if (trk_pt < min_trk_pt or trk_pt > max_trk_pt)
	continue;

      if (not trk_quality)
	continue;

      if (passTrackVertexAssociation(trk, 0, branch, 3.0))
	good_tracks.push_back(trk);
    }
    
    return good_tracks;
  }

  static std::pair<int,double> process_event_data(
    BranchPointerWrapper *branch,
    bool use_smeared_times,
    bool check_valid_times,
    bool use_z0, // whether or not to use 2D Clustering
    std::map<ScoreType,AnalysisObj>& analyses
  ) {
    int return_code = 0;
    double return_val = -1.;
    // check if vertex selection is correct & number of jets
    if (not branch->pass_basic_cuts()) return std::make_pair(-1,1.);

    // check if there is one forward jet with pt > 30 GeV
    if (not branch->pass_jet_pt_cut()) return std::make_pair(-1,1.);;

    int nForwardJet=0;
    branch->count_forward_jets(nForwardJet);
    
    std::vector<int> tracks = getAssociatedTracks(branch, min_track_pt, max_track_pt);

    int nForwardTrack=0, nForwardTrack_HS=0, nForwardTrack_PU=0;
    branch->count_forward_tracks(nForwardTrack,nForwardTrack_HS,nForwardTrack_PU,tracks, check_valid_times);

    // if (not branch->pass_forward_hs_tracks(nForwardTrack_HS)) return std::make_pair(-1,1.);
    
    double pu_ratio = (double)nForwardTrack_PU / (double)nForwardTrack;
    double reco_z = branch->reco_vtx_z[0];

    auto eff_fill_val_fjet     = folded(nForwardJet     , (int)fold_fjet) ;
    auto eff_fill_val_track    = folded(nForwardTrack   , (int)fold_track);
    auto eff_fill_val_hs_track = folded(nForwardTrack_HS, (int)fold_hs_track);
    auto eff_fill_val_pu_track = folded(nForwardTrack_PU, (int)fold_pu_track);
    auto eff_fill_val_pu_ratio = pu_ratio;

    if (debug) {
      std::cout << "nForwardJet = " << nForwardJet << std::endl;
      std::cout << "nForwardTrack = " << nForwardTrack << std::endl;
      std::cout << "nForwardTrack_HS = " << nForwardTrack_HS << std::endl;
      std::cout << "nForwardTrack_PU = " << nForwardTrack_PU << std::endl;
      std::cout << "Vertex_z = " << reco_z << std::endl;
    }
    
    std::vector<Cluster> clusters =
      clusterTracksInTime(tracks, branch, 3.0, 30.0,
			  use_smeared_times, check_valid_times, true, use_z0);

    if (debug) std::cout << "LEFT CLUSTERING" << std::endl;
    for (auto& [score,analysis]: analyses) {
      if (debug) std::cout << "SCORE: " << toString(score) << std::endl;
      if (debug) std::cout << "ACCESSING FJET" << std::endl;
      analysis["fjet"]->     FillTotal(eff_fill_val_fjet    );
      if (debug) std::cout << "ACCESSING FTRK" << std::endl;
      analysis["ftrack"]->   FillTotal(eff_fill_val_track   );
      if (debug) std::cout << "ACCESSING PUFRAC" << std::endl;
      analysis["pu_frac"]->  FillTotal(eff_fill_val_pu_ratio);
      if (debug) std::cout << "ACCESSING HS TRK" << std::endl;
      analysis["hs_track"]-> FillTotal(eff_fill_val_hs_track);
      if (debug) std::cout << "ACCESSING PU TRK" << std::endl;
      analysis["pu_track"]-> FillTotal(eff_fill_val_pu_track);
    }
    
    std::map<ScoreType,Cluster> chosen;
    if (clusters.size() != 0)
      chosen = chooseCluster(clusters, branch);
    if (debug) std::cout << "Chose My Clusters" << std::endl;

    // run HGTD Clustering (simultaneous)
    std::vector<Cluster> hgtd_clusters =
      clusterTracksInTime(tracks, branch, 3.0, -1, false, true, false, false);

    if (branch->reco_vtx_valid[0] == 1 and analyses.count(ScoreType::HGTD))
      chosen[ScoreType::HGTD] = chooseHGTDCluster(hgtd_clusters, branch);

    if (debug) std::cout << "Chose clusters" << std::endl;

    bool passes_hgtd = false, passes_mine = false;
    
    for (auto& [score, analysis] : analyses) {
      if (debug) std::cout << "Filling Scores: " << toString(score) << std::endl;
      if (branch->reco_vtx_valid[0] == 0 and score == HGTD)
	continue;

      if (clusters.size() == 0 and score != HGTD)
	continue;

      if (!chosen.count(score) and score != HGTD)
	continue;
      
      if (debug) std::cout << "Attempting to access chosen[" << toString(score) << "]" << std::endl;
      Cluster scored = chosen[score];
      if (debug) std::cout << "Attempting to access scores values, size: " << scored.values.size() << std::endl;
      double score_based_time = score == HGTD ? branch->reco_vtx_time[0] : scored.values.at(0);
      if (debug) std::cout << "Attempting to access purity" << std::endl;
      double cluster_purity = clusters.size() != 0 ? scored.purity : 0;
      if (debug) std::cout << "Purity: " << cluster_purity << std::endl;
      double diff = score_based_time - branch->truth_vtx_time[0];
      if (debug) std::cout << "Diff: " << diff << std::endl;

      if (score == ScoreType::TRKPTZ)
	return_val = score_based_time;
      
      analysis.inclusive_reso->Fill(diff);
      analysis.inclusive_purity->Fill(cluster_purity);
      if (scored.passEfficiency(branch)) {
	analysis["fjet"]->     FillPass(eff_fill_val_fjet    );
	analysis["ftrack"]->   FillPass(eff_fill_val_track   );
	analysis["pu_frac"]->  FillPass(eff_fill_val_pu_ratio);
	analysis["hs_track"]-> FillPass(eff_fill_val_hs_track);
	analysis["pu_track"]-> FillPass(eff_fill_val_pu_track);
	if (score == HGTD)
	  passes_hgtd = true;
	if (score == TRKPTZ)
	  passes_mine = true;
      }
      
      // fill diff hists
      analysis["fjet"]->     FillDiff(nForwardJet     , diff);
      analysis["ftrack"]->   FillDiff(nForwardTrack   , diff);
      analysis["pu_frac"]->  FillDiff(pu_ratio        , diff);
      analysis["hs_track"]-> FillDiff(nForwardTrack_HS, diff);
      analysis["pu_track"]-> FillDiff(nForwardTrack_PU, diff);
      
      // fill purities
      analysis["fjet"]->     FillPurity(nForwardJet     , cluster_purity);
      analysis["ftrack"]->   FillPurity(nForwardTrack   , cluster_purity);
      analysis["pu_frac"]->  FillPurity(pu_ratio        , cluster_purity);
      analysis["hs_track"]-> FillPurity(nForwardTrack_HS, cluster_purity);
      analysis["pu_track"]-> FillPurity(nForwardTrack_PU, cluster_purity);
    }
    if ((not passes_hgtd) and passes_mine)
      return_code = 2;
    else if (passes_mine)
      return_code = 1;

    if (not passes_mine) return_code = 0;
    
    return std::make_pair(return_code, return_val);
  }
}
#endif // EVENT_PROCESSING_H
