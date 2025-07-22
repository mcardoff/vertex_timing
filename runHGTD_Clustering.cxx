#include "common_includes.h" // Includes clustering_utilities.h

using namespace myutl;

void runHGTD_Clustering(std::string Number, Long64_t event_num) {
  // Initialize TChain to read the ntuple files
  TChain chain("ntuple");
  // Setup the chain using the utility function from clustering_utilities.h
  setup_chain(chain, Number);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  // Set the specific entry to read
  // chain.GetEntry(event_num);
  reader.SetEntry(event_num);

  // Get tracks associated with the primary vertex, using the utility function
  // min_track_pt is a static variable defined in clustering_utilities.h
  std::vector<int> tracks = getAssociatedTracks(&branch, min_track_pt);

  // Collect selected track times for finding min/max range for histogram
  std::vector<float> selected_times;
  for (auto idx : tracks) {
    selected_times.push_back(branch.track_time[idx]);
  }

  // Get truth and reconstructed vertex times
  double truth_vertex_time = branch.truth_vtx_time[0];
  double reco_vertex_time = branch.reco_vtx_time[0];

  // Calculate the min and max time for histogram range, including vertex times
  float min_time = truth_vertex_time;
  float max_time = truth_vertex_time;

  if (!selected_times.empty()) {
    // min_time = std::min({*std::min_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
    // max_time = std::max({*std::max_element(selected_times.begin(), selected_times.end()), truth_vertex_time, reco_vertex_time});
  } else {
      // If no tracks, set a default range or handle appropriately
      min_time = truth_vertex_time - 100.0; // Example: 100 ps window around truth vertex
      max_time = truth_vertex_time + 100.0;
  }
  
  // Extend the range slightly for better visualization
  float range = max_time - min_time;
  float extended_min_time = min_time - 0.05 * range;
  float extended_max_time = max_time + 0.05 * range;

  // Perform clustering using the utility function `clusterTracksInTime`
  // Parameters: track_indices, BranchPointerWrapper, distance_cut, smear_res, use_smeared_times, check_time_valid, isCone
  // For this example, we'll use a distance cut of 3.0, no smearing, check time validity, and use cone clustering.
  std::vector<Cluster> clusters = clusterTracksInTime(
      tracks,            // Indices of tracks to cluster
      &branch,           // Pointer to the BranchPointerWrapper
      3.0,               // distance_cut (e.g., 3.0 sigma for cone clustering)
      -1,                // smear_res (not used if use_smeared_times is false)
      false,             // use_smeared_times (set to false for real data)
      true,              // check_time_valid (only cluster tracks with valid times)
      true               // isCone (true for cone clustering, false for simultaneous)
  );

  // Iterate through the found clusters and print their properties
  for (int i = 0; i < clusters.size(); i++) {
    auto cluster = clusters.at(i);
    std::cout << "---------\n";
    std::cout << "t: " << cluster.values.at(0) << "\n";
    for (auto idx: cluster.track_indices) {
      std::cout << idx << "\n";
    }
    std::cout << "---------\n";
  }
}
