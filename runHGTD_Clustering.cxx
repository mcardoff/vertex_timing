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
       useConeClustering = true, useZ0 = false, usePURemoval = false;

  if (usePURemoval) {
    std::vector<int> pu_filtered_tracks;
    for (auto trk : tracks) {
      // auto z0_trk = branch.trackZ0[trk];
      // auto var_z0_trk = branch.trackVarZ0[trk];
      // double closest_nsigma = 1.0e6;
      int closest_reco_vtx = branch.trackToTruthvtx[trk];
      
      // for (int i_vtx = 0; i_vtx < branch.recoVtxZ.GetSize(); ++i_vtx) {
      // 	double z_vtx = branch.recoVtxZ[i_vtx];
      //   double nsigma = std::abs(z_vtx - z0_trk) / std::sqrt(var_z0_trk);
      //   if (nsigma < closest_nsigma) {
      //     closest_nsigma = nsigma;
      // 	  closest_reco_vtx = i_vtx;
      // 	}
      // }

      if (closest_reco_vtx == 0) {
	pu_filtered_tracks.push_back(trk);
      }
    }

    tracks = pu_filtered_tracks;
  }

  // gRandom->SetSeed(21);

  std::vector<Cluster> clusters =
    clusterTracksInTime(tracks, &branch, 2.0,
			useSmearTimes,
			useValidTimesOnly,
			30.0,
			useConeClustering,
			useZ0);

  // auto clustermerge = mergeClusters(clusters[0], clusters[1]);
  // clustermerge = mergeClusters(clustermerge, clusters[2]);
  // clusters = {clustermerge};
  // clusters.push_back(clustermerge);
  // clusters.erase(clusters.begin()+0);
  // clusters.erase(clusters.begin()+1);

  for (int j = 0; j < clusters.size(); j++) {
    auto cluster = clusters.at(j);
    auto score = cluster.scores.at(TRKPTZ);
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
