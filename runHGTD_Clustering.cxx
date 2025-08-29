#include "clustering_includes.h"
#include "event_processing.h"

using namespace myutl;

void runHGTD_Clustering(std::string Number, Long64_t event_num, bool smear_times, bool use_near) {
  // Initialize TChain to read the ntuple files
  TChain chain("ntuple");
  // Setup the chain using the utility function from clustering_utilities.h
  setup_chain(chain, Number);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  reader.SetEntry(event_num);
  
  std::vector<int> tracks = getAssociatedTracks(&branch, min_track_pt,max_track_pt);

  bool
    use_smear_times = false,
    use_valid_times_only = true,
    use_cone_clustering = true,
    use_z0 = false;
  

  std::vector<Cluster> clusters =
    clusterTracksInTime(tracks, &branch, 3.0, 30.0,
			use_smear_times,
			use_valid_times_only,
			use_cone_clustering,
			use_z0);

  // auto clustermerge = mergeClusters(clusters[0], clusters[2]);
  // clusters.push_back(clustermerge);
  // clusters.erase(clusters.begin()+0);
  // clusters.erase(clusters.begin()+1);

  for (int j = 0; j < clusters.size(); j++) {
    auto cluster = clusters.at(j);
    auto score = cluster.scores.at(TRKPT);
    std::cout << "---------\n";
    std::cout << "t: " << cluster.values.at(0) << "\n";
    if (cluster.values.size() > 1) std::cout << "z: " << cluster.values.at(1) << "\n";
    std::cout << "score: " << score << "\n";
    for (int i=0; i < cluster.track_indices.size(); i++) {
      std::cout << cluster.track_indices[i] << "," << cluster.all_times[i] << "\n";
    }
    std::cout << "passes? " << cluster.passEfficiency(&branch) << std::endl;
    std::cout << "---------\n";
  }
}
