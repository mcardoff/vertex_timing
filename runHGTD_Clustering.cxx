#include "common_includes.h" // Includes clustering_utilities.h

using namespace myutl;

void runHGTD_Clustering(std::string Number, Long64_t event_num, bool smear_times, bool use_near) {
  // Initialize TChain to read the ntuple files
  TChain chain("ntuple");
  // Setup the chain using the utility function from clustering_utilities.h
  setup_chain(chain, Number);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  reader.SetEntry(event_num);
  
  std::vector<int> tracks = getAssociatedTracks(&branch, min_track_pt);

  std::vector<Cluster> clusters =
    clusterTracksInTime(tracks, &branch, 3.0, 10.0, false, true, true, true);

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
