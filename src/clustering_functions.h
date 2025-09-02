#ifndef CLUSTERING_FUNCTIONS_H
#define CLUSTERING_FUNCTIONS_H

#include "clustering_includes.h"
#include "clustering_constants.h"
#include "clustering_structs.h"

namespace myutl {

  static double getSmearedTrackTime(
    int idx, double res, BranchPointerWrapper *branch
  ) {
    int particle_idx = branch->track_to_particle[idx]; // this is anything
    int truthvtx_idx = branch->track_to_truthvtx[idx]; // this is anything
    if (debug) std::cout << "Particle_idx: " << particle_idx << std::endl;
    if (debug) std::cout << "TruthVtx_idx: " << truthvtx_idx << std::endl;
    double t_part, smear_res = res;
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
    return gRandom->Gaus(t_part,smear_res);
  }
  
  static double getDistanceBetweenClusters(
    Cluster a, Cluster b
  ) {
    std::vector<double> a_v = a.values;
    std::vector<double> a_s = a.sigmas;
    std::vector<double> b_v = b.values;
    std::vector<double> b_s = b.sigmas;

    std::vector<double> metric = {1.,-1.};

    if (a_v.size() != b_v.size())
      std::cout << "Uh ohhh" << std::endl;

    std::vector<double> distances(a_v.size());
    for (int i=0; i < a_v.size(); i++) {
    
      double diff_i = (a_v.at(i) - b_v.at(i));
      double denom_i = sqrt(a_s.at(i)*a_s.at(i) + b_s.at(i)*b_s.at(i));
      distances.at(i) = diff_i/denom_i;
    }

    double dsqr = 0.0;
    for (int i=0; i < distances.size(); i++){
      auto d = distances[i];
      dsqr += metric.at(i)*d*d;
    }
    return sqrt(dsqr);
  }

  static Cluster mergeClusters(
    Cluster a, Cluster b
  ) {
    Cluster merged_cluster;
    merged_cluster.was_merged = true;
  
    for (size_t i = 0; i < a.values.size(); i++) {
      double v1 = a.values.at(i);
      double v2 = b.values.at(i);
      double s1 = a.sigmas.at(i);
      double s2 = b.sigmas.at(i);

      double new_cluster_value =
        (v1 / (s1 * s1) + v2 / (s2 * s2)) / (1. / (s1 * s1) + 1. / (s2 * s2));
      double new_cluster_sigma =
	sqrt((s2 * s1) * (s2 * s1) / (s2 * s2 + s1 * s1));
      merged_cluster.values.push_back(new_cluster_value);
      merged_cluster.sigmas.push_back(new_cluster_sigma);
    }

    for (auto t: a.all_times)
      merged_cluster.all_times.push_back(t);
    for (auto t: b.all_times)
      merged_cluster.all_times.push_back(t);

    for (auto i: a.track_indices)
      merged_cluster.track_indices.push_back(i);
    for (auto i: b.track_indices)
      merged_cluster.track_indices.push_back(i);

    merged_cluster.n_constituents = a.n_constituents+b.n_constituents;
    merged_cluster.max_pt_cluster = a.max_pt_cluster || b.max_pt_cluster;

    for (ScoreType score: enum_vec)
      merged_cluster.scores[score] = a.scores[score]+b.scores[score];

    merged_cluster.scores[ScoreType::IDEAL] = merged_cluster.values.at(0);

    // if (merged_cluster.max_pt_cluster)
    //   merged_cluster.scores[ScoreType::MAXPT] = 1e6;

    return merged_cluster;
  }

  static std::vector<Cluster> makeSimpleClusters(
    std::vector<int> track_indices,
    BranchPointerWrapper *branch,
    bool use_smeared_times,                         // if true, use smeared_times map
    const std::map<int, double>& smeared_times_map, // map of smeared times
    const std::map<int, double>& smeared_timeRes_map, // map of resolution for smeared times if used
    double fixed_res,                               // fixed resolution for smeared times
    bool check_time_valid,                          // if true, apply branch->track_time_valid check
    bool usez0

  ) {
    std::vector<Cluster> simpleClusters;
    for (auto idx: track_indices) {
      if (!check_time_valid || branch->track_time_valid[idx] == 1) {
	double trk_time;
	double trk_res; 
	if (use_smeared_times) {
	  trk_time = smeared_times_map.at(idx);
	  trk_res  = smeared_timeRes_map.at(idx);
	} else {
	  trk_time = branch->track_time[idx];   
	  trk_res  = branch->track_time_res[idx];
	}

	double trk_pt = branch->track_pt[idx];
	double trk_z0 = branch->track_z0[idx];
	double trk_res_z0 = branch->track_var_z0[idx];
	double z_significance = std::abs(trk_z0 - branch->reco_vtx_z[0]) / trk_res_z0;

	std::map<ScoreType,double> scores;
	scores[ScoreType::HGTD] = 0; // do not use this
	// scores[ScoreType::SIG] = exp(-z_significance);
	// scores[ScoreType::MAXPT] = 0;
	scores[ScoreType::TRKPT] = trk_pt;
	scores[ScoreType::SUMPT2] = trk_pt*trk_pt;
	// scores[ScoreType::SIGTRKPT] = trk_pt*exp(-z_significance);
	// scores[ScoreType::MAXHS] = branch->track_to_truthvtx[idx] == 0 ? 1 : 0;
	// scores[ScoreType::JETPTDR] = branch->calc_jetpt_dr_score(idx);
	// scores[ScoreType::TRKPTDR] = branch->calc_trkpt_dr_score(idx);
	// scores[ScoreType::IDEAL] = trk_time;

	Cluster cluster;
	if (usez0)
	  cluster = {{trk_time,trk_z0}, {trk_res,trk_res_z0}, {trk_time}, {idx}, scores};
	else
	  cluster = {{trk_time}, {trk_res}, {trk_time}, {idx}, scores};
	simpleClusters.push_back(cluster);
      }
    }
    return simpleClusters;
  }

  static void doSimultaneousClustering(
    std::vector<Cluster> *collection, double dist_cut
  ) {
    double distance = 1.e30;
    while (collection->size() > 1) {
      // std::cout << "entering while loop" << std::endl;

      int i0 = 0;
      int j0 = 0;

      distance = getDistanceBetweenClusters(collection->at(0), collection->at(1));

      for (size_t i = 0; i < collection->size(); i++) {
	for (size_t j = i + 1; j < collection->size(); j++) {
	  Cluster a = collection->at(i);
	  Cluster b = collection->at(j);

	  if (a.was_merged or b.was_merged)
	    continue;
	
	  double current_distance =
	    getDistanceBetweenClusters(a, b);
	  if (current_distance <= distance) {
	    distance = current_distance;
	    i0 = i;
	    j0 = j;
	  }
	} // j loop
      } // i loop

      // fuse closest two vertices
      if (distance < dist_cut and i0 != j0) {
	Cluster new_cluster = mergeClusters(collection->at(i0),collection->at(j0));
	collection->erase(collection->begin()+j0);
	if (i0 < j0)
	  collection->erase(collection->begin()+i0);
	else
	  collection->erase(collection->begin()+(i0-1));

	collection->push_back(new_cluster);
      } else {
	if (std::find_if(collection->begin(), collection->end(),
			 [](Cluster a) {return a.was_merged;}) != collection->end()) {
	  for (int idx=0; idx < collection->size(); ++idx) {
	    collection->at(idx).was_merged = false;
	  }
	
	} else {
	  break;
	}
      }
      // break;
    } // While
  }

  static void doConeClustering(
    std::vector<Cluster> *collection, double dist_cut
  ) {
    if (debug) std::cout << "entering while loop" << std::endl;
    while (collection->size() > 1) {
      
      // find seed
      int i0 = -1;
      double seed_value = -1.;
      for (int i=0; i < collection->size(); ++i) {
        if (collection->at(i).was_merged)
	  continue;
	Cluster check = collection->at(i);
	double check_value = check.scores.at(TRKPT);
	
	if (seed_value < check_value) {
	  seed_value = check_value;
	  i0 = i;
	}
      }
      if (i0 == -1) break;
      
      if (debug) std::cout << "Found Seed" << std::endl;

      // find merge candidates
      Cluster seed = collection->at(i0);
      std::vector<int> to_merge;
      for (int j = 0; j < collection->size(); ++j) {
	if (j == i0) continue;
	double d = getDistanceBetweenClusters(seed, collection->at(j));
	if (d < dist_cut) {
	  to_merge.push_back(j);
	}
      }
      
      if (debug) std::cout << "Found " << to_merge.size() << " merge candidates" << std::endl;

      if (!to_merge.empty()) {
	// 1) Build the merged cone without erasing yet
	Cluster new_cluster = seed;
	for (int idx : to_merge) {
	  new_cluster = mergeClusters(new_cluster, collection->at(idx));
	}

	new_cluster.was_merged = true;

	// 2) Erase seed + all candidates in descending order
	to_merge.push_back(i0);
	std::sort(to_merge.begin(), to_merge.end(), std::greater<int>());
	for (int idx : to_merge) {
	  collection->erase(collection->begin() + idx);
	}

	// 3) Push back the final cone
	collection->push_back(new_cluster);
      } else {
	// no merges possible, seed is isolated → just move it to the back
	Cluster seed_only = collection->at(i0);
	seed_only.was_merged = true;
	collection->erase(collection->begin() + i0);
	collection->push_back(seed_only);
      }
    } // While
    if (debug) std::cout << "Finished Clustering" << std::endl;
  }

  static std::vector<Cluster> clusterTracksInTime(
     std::vector<int> track_indices,
     BranchPointerWrapper *branch,
     double distance_cut,    // Cone size/distance cut
     double smear_res,       // fixed resolution for smeared times
     bool use_smeared_times, // for ideal cases
     bool check_time_valid,  // for ideal efficiency
     bool useCone,
     bool usez0
  ) {
    if (track_indices.empty() && debug)
      std::cout << "EMPTY!!!!" << std::endl;

    std::map<int,double> smeared_times_map, smeared_timeRes_map;

    // if use smeared times, generated
    if (use_smeared_times) {
      for (auto idx: track_indices) {
	double res = smear_res;
	if (branch->track_time_valid[idx] == 1 and branch->track_hgtd_hits[idx] > 0 )
	  res /= std::sqrt(branch->track_hgtd_hits[idx]);
	if (!check_time_valid || branch->track_time_valid[idx] == 1) {
	  double smeared_time = getSmearedTrackTime(idx, res, branch);
	  smeared_times_map.emplace(idx, smeared_time);
	  smeared_timeRes_map.emplace(idx, res);
	}
      }
    }

    // std::sort(track_indices.begin(), track_indices.end(),
    // 	      [&branch](int a, int b){return branch->track_pt[a] > branch->track_pt[b];});

    std::vector<Cluster> collection = makeSimpleClusters(track_indices, branch, use_smeared_times,
							 smeared_times_map, smeared_timeRes_map, smear_res, check_time_valid, usez0);
    
    if (not useCone) {
      doSimultaneousClustering(&collection, distance_cut);
    } else {
      doConeClustering(&collection, distance_cut);
      // collection = doConeClustering(track_indices, branch, distance_cut, use_smeared_times, smeared_times_map, smear_res, check_time_valid, usez0);
    }

    // if using z0, have to update the trkpt score:
    if (usez0) {
      for (Cluster& cluster: collection) {
	double dz = std::abs(cluster.values.at(1)-branch->reco_vtx_z[0]);
	auto oldscore = std::pow(cluster.scores.at(TRKPT),0.9);
	cluster.scores[TRKPTZ] = oldscore*exp(-1.5*dz);
	cluster.scores[DZ] = exp(-dz);
      }
    } else {
      for (Cluster& cluster: collection) {
	double znum=0., zden=0.;
	double dnum=0., dden=0.;
	// auto max_pt_idx = -1;
	auto max_pt = -1.;
	for (auto trk: cluster.track_indices) {
	  auto
	    trk_z = branch->track_z0[trk],
	    trk_var_z = branch->track_var_z0[trk],
	    trk_d0 = branch->track_d0[trk],
	    trk_var_d0 = branch->track_var_d0[trk];
	   // std::cout << "trk_z[" << trk << "] = " << trk_z << std::endl;
	  // std::cout << "trk_var_z[" << trk << "] = " << trk_var_z << std::endl;
	  znum += trk_z/(trk_var_z);
	  zden += 1/(trk_var_z);
	  dnum += trk_d0/(trk_var_d0);
	  dden += 1/(trk_var_d0);
	  if (branch->track_pt[trk] > max_pt)
	    max_pt = branch->track_pt[trk];
	}

	double z = znum/zden, zsigma = 1/std::sqrt(zden);
	double dz = std::abs(z-branch->reco_vtx_z[0]);
	// std::cout << "dz: " << dz << std::endl;
	// std::cout << "d0: " << dnum/dden << std::endl;
	// std::cout << "sumpt: " << cluster.scores.at(TRKPT) << std::endl;
	auto oldscore = std::pow(cluster.scores.at(TRKPT),0.9);
	cluster.scores[TRKPTZ] = oldscore*exp(-1.5*dz);
      }
    }
    
    for (Cluster& cluster: collection)
      cluster.calcPurity(branch);
    if (debug) std::cout << "Finished Clustering" << std::endl;
    return collection;
  }

  static Cluster chooseHGTDCluster(
    std::vector<Cluster> collection,
    BranchPointerWrapper *branch
  ) {
    double min_diff = 1e50, reco_time = branch->reco_vtx_time[0];
    Cluster min_cluster;
    for (Cluster cluster: collection) {
      double this_diff = std::abs(cluster.values.at(0) - reco_time);
      if (this_diff < min_diff) {
	min_diff = this_diff;
	min_cluster = cluster;
      }
    }
    return min_cluster;
  }
  
  static std::map<ScoreType, Cluster> chooseCluster(
    std::vector<Cluster> collection,
    BranchPointerWrapper *branch
  ) {
    std::map<ScoreType,Cluster> output;
    if (debug) std::cout << "Choosing score" << std::endl;
    for (ScoreType score: enum_vec) {
      if (debug) std::cout << "SCORE: " << toString(score)<< std::endl;
      if (score == ScoreType::HGTD) {
	continue;
      }

      if (score == ScoreType::PASS) {
	if (debug) std::cout << "Choosing pass score" << std::endl;
	for (Cluster& cluster: collection) {
	  if (cluster.passEfficiency(branch))
	    output[score] = cluster;
	}
	continue;
      }
      
      output[score] = collection[0]; // final time we are giving to the user
      double max_score = score != ScoreType::IDEAL
	? output[score].scores.at(score)
	: std::abs(output[score].scores.at(score) - branch->truth_vtx_time[0]);
    
      for (Cluster& cluster: collection) {
	double comp_score = score != ScoreType::IDEAL
	  ? cluster.scores[score]
	  : std::abs(cluster.scores[score] - branch->truth_vtx_time[0]);

	bool query = score != ScoreType::IDEAL
	  ? comp_score > max_score
	  : comp_score < max_score;

	if (query) {
	  max_score = comp_score;
	  output[score] = cluster;
	}
      }
    }
  
    return output;
  }
}
#endif // CLUSTERING_FUNCTIONS_H
