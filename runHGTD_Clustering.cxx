#include <RtypesCore.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TH1.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>

struct Cluster {
  std::vector<double> values;
  std::vector<double> sigmas;
  std::vector<double> all_times;
  std::vector<int>    indices;
  double score;
  bool wasMerged = false;
  int nConstituents=1;
};

struct BranchPointerWrapper {
  // track variables
  std::vector<float> *track_z0 = nullptr;
  std::vector<float> *track_pt = nullptr;
  std::vector<float> *track_eta = nullptr;
  std::vector<float> *track_time = nullptr;
  std::vector<float> *track_time_res = nullptr;
  std::vector<float> *track_var_z0 = nullptr;
  std::vector<int>   *track_to_truthvtx = nullptr;
  std::vector<int>   *track_time_valid = nullptr;

  // vertex variables
  std::vector<float> *truth_vtx_z = nullptr;
  std::vector<float> *truth_vtx_time = nullptr;
  
  /// reco vertex
  std::vector<float> *reco_vtx_z = nullptr;
  std::vector<float> *reco_vtx_time = nullptr;
  std::vector<float> *reco_vtx_timeRes = nullptr;
  std::vector<int>   *reco_vtx_time_valid = nullptr;
};

bool passTrackVertexAssociation(
  int track_idx, int vertex_idx, BranchPointerWrapper *branch,
  double min_trk_pt, double significance_cut) {
  double
    trk_z0    = branch->track_z0->at(track_idx),
    trk_z_var = branch->track_var_z0->at(track_idx),
    trk_pt    = branch->track_pt->at(track_idx),
    trk_eta   = branch->track_eta->at(track_idx),
    vx_z      = branch->reco_vtx_z->at(vertex_idx),
    vx_z_var  = 0.0;
  
  if (trk_pt < min_trk_pt || std::abs(trk_eta) > 4.0)
    return false;

  double nsigma = std::abs(trk_z0 - vx_z) / std::sqrt(trk_z_var + vx_z_var);
  return nsigma < significance_cut;
}

double getDistanceBetweenClusters(Cluster a, Cluster b) {
  std::vector<double> a_v = a.values;
  std::vector<double> a_s = a.sigmas;
  std::vector<double> b_v = b.values;
  std::vector<double> b_s = b.sigmas;

  if (a_v.size() != b_v.size())
    std::cout << "Uh ohhh" << std::endl;

  std::vector<double> distances(a_v.size());
  for (int i=0; i < a_v.size(); i++) {
    
    double diff_i = (a_v.at(i) - b_v.at(i));
    double denom_i = sqrt(a_s.at(i)*a_s.at(i) + b_s.at(i)*b_s.at(i));
    distances.at(i) = diff_i/denom_i;
  }

  double dsqr = 0.0;
  for (double d : distances){
    dsqr += d*d;
  }
  return sqrt(dsqr);
}

Cluster mergeClusters(Cluster a, Cluster b) {
  Cluster merged_cluster;
  
  for (size_t i = 0; i < a.values.size(); i++) {
    double v1 = a.values.at(i);
    double v2 = b.values.at(i);
    double s1 = a.sigmas.at(i);
    double s2 = b.sigmas.at(i);

    double new_cluster_value =
        (v1 / (s1 * s1) + v2 / (s2 * s2)) / (1. / (s1 * s1) + 1. / (s2 * s2));
    double new_cluster_sigma =
      sqrt((s2 * s1) * (s2 * s1) / (s2 * s2 + s1 * s1));
    // std::cout << "in mergeClusters" << std::endl;
    merged_cluster.values.push_back(new_cluster_value);
    merged_cluster.sigmas.push_back(new_cluster_sigma);
  }

  for (auto t: a.all_times)
    merged_cluster.all_times.push_back(t);
  for (auto t: b.all_times)
    merged_cluster.all_times.push_back(t);

  for (auto t: a.indices)
    merged_cluster.indices.push_back(t);
  for (auto t: b.indices)
    merged_cluster.indices.push_back(t);

  merged_cluster.nConstituents = a.nConstituents+b.nConstituents;
  merged_cluster.score = a.score+b.score;
  merged_cluster.wasMerged = true;

  return merged_cluster;
}

void doSimpleConeClustering(std::vector<Cluster> *collection, BranchPointerWrapper *branch) {
  while (collection->size() > 1) {
    std::cout << "entering while loop" << std::endl;

    int i0 = 0; // this should be the cluster with highest pT (maybe)
    for (size_t idx = 0; idx < collection->size(); ++idx) {
      Cluster this_cluster = collection->at(idx);
      Cluster that_cluster = collection->at(i0);
      if (this_cluster.wasMerged) {
	continue;
      }
      double this_pt = branch->track_pt->at(this_cluster.indices.at(0)), max_pt = branch->track_pt->at(that_cluster.indices.at(0));
      if (this_pt > max_pt) {
	i0 = idx;
      }
    }

    std::cout << "i0 = " << i0 << std::endl;

    int j0 = 0;

    Cluster a = collection->at(i0);
    double sa = a.sigmas.at(0);
    double distance = 1000;

    for (size_t j = 0; j < collection->size(); j++) {
      if (j == i0)
	continue;
      Cluster b = collection->at(j);
      double sb = b.sigmas.at(0);

      if (a.wasMerged or b.wasMerged)
	continue;

      double current_distance = std::abs(a.values.at(0) - b.values.at(0)) / std::sqrt(sa*sa + sb*sb);
      std::cout << "Current Distance in t: " << current_distance << std::endl;
	
      if (current_distance <= distance) {
	distance = current_distance;
	j0 = j;
      }
      std::cout << "Saved Distance: " << distance << std::endl;
    } // j loop

    // fuse closest two vertices
    if (distance < 3.0 and i0 != j0) {
      Cluster new_cluster = mergeClusters(collection->at(i0),collection->at(j0));
      if (i0 > j0){
	collection->erase(collection->begin()+i0);
	collection->erase(collection->begin()+j0);
      } else {
	collection->erase(collection->begin()+j0);
	collection->erase(collection->begin()+i0);
      }
      collection->push_back(new_cluster);
    } else {
      if (std::find_if(collection->begin(), collection->end(),
		       [](Cluster a) {return a.wasMerged;}) != collection->end()) {
	for (int idx=0; idx < collection->size(); ++idx) {
	  collection->at(idx).wasMerged = false;
	}
	  
      } else {
	break;
      }
    }
    // break;
  } // While
}

std::vector<Cluster> doConeClustering(std::vector<int> indices, BranchPointerWrapper *branch) {
  std::vector<bool> checked(branch->track_pt->size(), false);
  std::vector<Cluster> clusters;

  std::vector<int> timefiltered_indices;
  for (auto idx: indices) {
    if (branch->track_time_valid->at(idx) == 1)
      timefiltered_indices.push_back(idx);
  }

  while (true) {
    // Step 1: Find the highest-pt unchecked track
    int max_pt_idx = -1;
    double max_pt = -1.0;

    for (int idx : timefiltered_indices) {
      if (checked[idx]) continue;
      double this_pt = branch->track_pt->at(idx);
      if (this_pt > max_pt) {
        max_pt = this_pt;
        max_pt_idx = idx;
      }
    }

    if (max_pt_idx == -1) break; // all tracks have been checked

    // Step 2: Build cluster around max_pt_idx
    checked[max_pt_idx] = true;
    double t1 = branch->track_time->at(max_pt_idx);
    double s1 = branch->track_time_res->at(max_pt_idx);

    std::vector<int> to_merge = {max_pt_idx};

    for (int idx : timefiltered_indices) {
      if (checked[idx]) continue;

      double t2 = branch->track_time->at(idx);
      double s2 = branch->track_time_res->at(idx);
      double nsigma = std::abs(t1 - t2) / std::sqrt(s1 * s1 + s2 * s2);

      if (nsigma < 3.0) {
        to_merge.push_back(idx);
        checked[idx] = true;
      }
    }

    // Step 3: Form cluster from selected indices
    std::vector<Cluster> simpleclusters;
    for (int idx : to_merge) {
      double trk_time = branch->track_time->at(idx);
      double res = branch->track_time_res->at(idx);
      double significance = std::abs(branch->track_z0->at(idx) - branch->reco_vtx_z->at(0)) / std::sqrt(branch->track_var_z0->at(idx));
      Cluster cluster = {{trk_time}, {res}, {trk_time}, {idx}, exp(-significance)};
      if (not simpleclusters.empty()) {
	Cluster newCluster = mergeClusters(simpleclusters[0],cluster);
	simpleclusters.at(0) = newCluster;
      } else {
	simpleclusters.push_back(cluster);
      }
    }
    clusters.push_back(simpleclusters.at(0));
  }

  return clusters;
}

void doEagerClustering(
  std::vector<Cluster> *collection, double dist_cut) {

  double distance = 1.e30;
  while (collection->size() > 1) {

    int i0 = 0;
    int j0 = 0;

    distance = getDistanceBetweenClusters(collection->at(0), collection->at(1));

    for (size_t i = 0; i < collection->size(); i++) {
      for (size_t j = i + 1; j < collection->size(); j++) {
	Cluster a = collection->at(i);
	Cluster b = collection->at(j);
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
    if (distance < dist_cut) {
      Cluster new_cluster = mergeClusters(collection->at(i0),collection->at(j0));
      
      collection->erase(collection->begin()+j0);
      if (i0 < j0)
	collection->erase(collection->begin()+i0);
      else
	collection->erase(collection->begin()+(i0-1));

      collection->push_back(new_cluster);
    } else {
      break;
    }
  } // While
}

void doSimultaneousClustering(
  std::vector<Cluster> *collection, double dist_cut) {

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

	if (a.wasMerged or b.wasMerged)
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
		       [](Cluster a) {return a.wasMerged;}) != collection->end()) {
	for (int idx=0; idx < collection->size(); ++idx) {
	  collection->at(idx).wasMerged = false;
	}
	  
      } else {
	break;
      }
    }
    // break;
  } // While
}

std::vector<Cluster> clusterTracksInTime(
   std::vector<int> track_indices, double cluster_distance,
   BranchPointerWrapper *branch) {
  std::vector<Cluster> collection;
  for (auto idx: track_indices) {
    if (branch->track_time_valid->at(idx)) {
      double trk_time = branch->track_time->at(idx);
      double res = branch->track_time_res->at(idx);
      double significance = std::abs(branch->track_z0->at(idx) - branch->reco_vtx_z->at(0)) / std::sqrt(branch->track_var_z0->at(idx));
      Cluster cluster = {{trk_time}, {res}, {trk_time}, {idx}, exp(-significance)};
      collection.push_back(cluster);
    } else {
      // std::cout << "No valid time" << std::endl;
    }
  }

  // doSimultaneousClustering(&collection, 3.0);
  // doSimpleConeClustering(&collection, branch);
  collection = doConeClustering(track_indices, branch);
  // doEagerClustering(&collection, 3.0);
  
  return collection;
}

TTree *setup_tree(std::string Number, BranchPointerWrapper* branch) {
  TFile *file = TFile::Open(Form("../ntuple/user.scheong.42871997.Output._%s.SuperNtuple.root", Number.c_str()));
  TTree *tree = (TTree*)file->Get("ntuple");

  tree->SetBranchAddress("Track_z0"            , &branch->track_z0);
  tree->SetBranchAddress("Track_var_z0"        , &branch->track_var_z0);
  tree->SetBranchAddress("Track_pt"            , &branch->track_pt);
  tree->SetBranchAddress("Track_eta"           , &branch->track_eta);
  tree->SetBranchAddress("Track_time"          , &branch->track_time);
  tree->SetBranchAddress("Track_timeRes"       , &branch->track_time_res);
  tree->SetBranchAddress("Track_hasValidTime"  , &branch->track_time_valid);
  tree->SetBranchAddress("TruthVtx_z"          , &branch->truth_vtx_z);
  tree->SetBranchAddress("TruthVtx_time"       , &branch->truth_vtx_time);
  tree->SetBranchAddress("RecoVtx_z"           , &branch->reco_vtx_z);
  tree->SetBranchAddress("RecoVtx_time"        , &branch->reco_vtx_time);
  tree->SetBranchAddress("RecoVtx_timeRes"     , &branch->reco_vtx_timeRes);
  tree->SetBranchAddress("RecoVtx_hasValidTime", &branch->reco_vtx_time_valid);
    
  return tree;
}

std::vector<int> getAssociatedTracks(
  BranchPointerWrapper *branch, double min_trk_pt) {
  std::vector<int> good_tracks;

  for (int trk = 0; trk < branch->track_z0->size(); ++trk) {
    double
      trk_eta = branch->track_eta->at(trk),
      trk_pt  = branch->track_pt->at(trk);
    if (std::abs(trk_eta) < 2.4 || std::abs(trk_eta) > 4.0) {
      // if (trk == 58)
	// std::cout << "Out of eta acceptance " << trk_eta << std::endl;
      continue;
    }
    
    if (trk_pt < 1.0) {
      // if (trk == 58)
	// std::cout << "Out of pt acceptance" << std::endl;
      continue;
    }
    
    if (passTrackVertexAssociation(trk, 0, branch, min_trk_pt, 3.0))
      good_tracks.push_back(trk);
  }

  return good_tracks;
}

void runHGTD_Clustering(std::string Number, Long64_t event_num) {
  BranchPointerWrapper *branch = new BranchPointerWrapper;
  TTree *tree = setup_tree(Number, branch);
  tree->GetEntry(event_num);

  std::vector<int> tracks = getAssociatedTracks(branch, 1.0);
  std::vector<float> selected_times;
  for (auto idx: tracks) {
    selected_times.push_back(branch->track_time->at(idx));
    if (idx == 58)
      std::cout << "BLEH" << std::endl;
  }

  float truth_vertex_time = branch->truth_vtx_time->at(0);
  float reco_vertex_time = branch->reco_vtx_time->at(0);
  float min_time = std::min({*std::min_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
  float max_time = std::max({*std::max_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
  
  float range = max_time - min_time;
  float extended_min_time = min_time - 0.05 * range;
  float extended_max_time = max_time + 0.05 * range;

  // TH1F* ex_hist = new TH1F("track_time_hist_event",
			   // Form("Track Time Distribution Event %lld;Time (ps);A.U.", event_num),
			   // 50, extended_min_time, extended_max_time);
  std::vector<Cluster> clusters = clusterTracksInTime(tracks, 3.0, branch);
  std::vector<TH1F*> cluster_hist;
  // std::cout << clusters.size() << "\n"; 
  for (int i = 0; i < clusters.size(); i++) {
    // std::cout << "idx: " << i << std::endl;
    auto cluster = clusters.at(i);
    if (cluster.nConstituents < 2)
      continue;
    // std::cout << "Cluster has " << cluster.indices.size() << " tracks\n";
    // print out min and max time
    std::cout << "---------\n";
    std::cout << "t: " << cluster.values.at(0) << "\n";
    for (auto idx: cluster.indices) {
      std::cout << idx << "\n";
    }
    std::cout << "---------\n";
  }
}
