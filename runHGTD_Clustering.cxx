#include "common_includes.h" // Includes clustering_utilities.h

using namespace myutl;

void runHGTD_Clustering(std::string Number, Long64_t event_num, bool smear_times) {
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

  // Perform clustering using the utility function `clusterTracksInTime`
  // Parameters: track_indices, BranchPointerWrapper, distance_cut, smear_res, use_smeared_times, check_time_valid, isCone
  // For this example, we'll use a distance cut of 3.0, no smearing, check time validity, and use cone clustering.
  std::vector<Cluster> clusters = clusterTracksInTime(
      tracks,            // Indices of tracks to cluster
      &branch,           // Pointer to the BranchPointerWrapper
      3.0,               // distance_cut (e.g., 3.0 sigma for cone clustering)
      30.0,              // smear_res (not used if use_smeared_times is false)
      smear_times,       // use_smeared_times (set to false for real data)
      true,              // check_time_valid (only cluster tracks with valid times)
      true               // isCone (true for cone clustering, false for simultaneous)
  );

  // Iterate through the found clusters and print their properties
  for (int i = 0; i < clusters.size(); i++) {
    auto cluster = clusters.at(i);
    std::cout << "---------\n";
    std::cout << "t: " << cluster.values.at(0) << "\n";
    for (int i=0; i < cluster.track_indices.size(); i++) {
      std::cout << cluster.track_indices[i] << "," << cluster.all_times[i] << "\n";
    }
    std::cout << "---------\n";
  }
}
