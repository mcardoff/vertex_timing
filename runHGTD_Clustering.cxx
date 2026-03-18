#include "src/clustering_constants.h"
R__ADD_INCLUDE_PATH(/opt/homebrew/opt/boost/include)
R__ADD_LIBRARY_PATH(/opt/homebrew/opt/boost/lib)
R__LOAD_LIBRARY(libboost_filesystem)

#include "src/clustering_functions.h"
#include "src/event_processing.h"
#include <TRandom.h>

using namespace MyUtl;

void runHGTD_Clustering(std::string number, Long64_t eventNum) {
  // Initialize TChain to read the ntuple files
  TChain chain("ntuple");
  // Setup the chain using the utility function from clustering_utilities.h
  setupChain(chain, number);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  reader.SetEntry(eventNum);
  
  std::vector<int> tracks = getAssociatedTracks(&branch, MIN_TRACK_PT,MAX_TRACK_PT, 3.0);

  bool useSmearTimes = false, useValidTimesOnly = true,
       useZ0 = false, usePURemoval = false;

  if (usePURemoval) {
    std::vector<int> pu_filtered_tracks = filterTruthMatchedTracks(tracks, &branch, 3.0);
    
    tracks = pu_filtered_tracks;
  }

  // gRandom->SetSeed(21);

  ClusteringMethod method = ClusteringMethod::CONE;

  std::vector<Cluster> clusters =
    clusterTracksInTime(tracks, &branch, 3.0,
			useSmearTimes,
			useValidTimesOnly,
			30.0,
			method,
			useZ0);

  // auto clustermerge = mergeClusters(clusters[0], clusters[1]);
  // clustermerge = mergeClusters(clustermerge, clusters[2]);
  // clusters = {clustermerge};
  // clusters.push_back(clustermerge);
  // clusters.erase(clusters.begin()+0);
  // clusters.erase(clusters.begin()+1);

  // ===== REFINED STAGE — comment out this block to revert to plain cone output =====
  // Two-pass timing refinement: find the TRKPTZ winner among the 3σ cone clusters,
  // then recompute its time using only tracks within DIST_CUT_REFINE σ of the centroid.
  {
    if (!clusters.empty()) {
      auto bestIt = std::max_element(clusters.begin(), clusters.end(),
          [](const Cluster& a, const Cluster& b) {
            return a.scores.at(Score::TRKPTZ) < b.scores.at(Score::TRKPTZ);
          });
      *bestIt = refineClusterTiming(*bestIt, &branch, DIST_CUT_REFINE);
    }
  }
  // ===== END REFINED STAGE =====

  for (int j = 0; j < clusters.size(); j++) {
    auto cluster = clusters.at(j);
    auto score = cluster.scores.at(Score::TRKPTZ);
    std::cout << "---------\n";
    std::cout << "t: " << cluster.values.at(0) << "\n";
    if (cluster.values.size() > 1) std::cout << "z: " << cluster.values.at(1) << "\n";
    std::cout << "score: " << score << "\n";
    for (int i=0; i < cluster.trackIndices.size(); i++) {
      std::cout << cluster.trackIndices[i] << "," << cluster.allTimes[i] << "\n";
      
    }
    std::cout << "passes? " << cluster.passEfficiency(&branch) << std::endl;
    std::cout << "---------\n";
  }
}
