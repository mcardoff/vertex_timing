#include "../src/clustering_constants.h"
R__ADD_INCLUDE_PATH(/opt/homebrew/opt/boost/include)
R__ADD_LIBRARY_PATH(/opt/homebrew/opt/boost/lib)
R__LOAD_LIBRARY(libboost_filesystem)

#include "../src/clustering_functions.h"
#include "../src/event_processing.h"
#include "../src/idealised_timing.h"
#include <TRandom.h>
#include <TRandom3.h>
#include <iomanip>

using namespace MyUtl;

// assocNSigma: z0-significance cut for the clustering input. 3.0 is the
// display's historical default; rpt_v5_hist clusters from a 2.5 list, so pass
// 2.5 to see the clusters its t0s came from.
//
// idealTimes: replay rpt_v5's idealised world (src/idealised_timing.h) --
// print every timed track's smeared time ("smeared: idx,time"), the smeared
// vertex time ("ttruthvtx:"), and cluster the smeared world exactly as the
// waves_ideal row does ("tclosest:" is the cluster that row selects). Same
// seed, same draw order, so these are the times the histograms were filled
// with, not a fresh draw.
void runHGTD_Clustering(std::string filePath, Long64_t eventNum,
                        double assocNSigma = 3.0, bool idealTimes = false) {
  // Initialize TChain to read the ntuple file. A single Add() of the full
  // path (any sample -- vbf/zjets/dijet/local) replaces the old
  // setupChain(chain, std::string) single-file overload, which only worked
  // for the local default VBF ntuple's fixed naming convention.
  TChain chain("ntuple");
  chain.Add(filePath.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  reader.SetEntry(eventNum);

  // Full precision: the display re-applies the analysis' time gate to these
  // numbers, and a truncated time can flip a track sitting on the gate edge.
  std::cout << std::setprecision(17);

  std::vector<int> tracks = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, assocNSigma);

  bool useSmearTimes = false, useValidTimesOnly = true, useZ0 = false;

  // gRandom->SetSeed(21);

  ClusteringMethod method = ClusteringMethod::ITERATIVE;

  std::vector<Cluster> clusters;
  if (idealTimes) {
    TRandom3 rng(1);  // re-seeded per event inside makeIdealisedTiming
    const IdealisedTiming ideal = makeIdealisedTiming(&branch, rng);
    std::vector<int> timed;
    for (const auto& kv : ideal.times) timed.push_back(kv.first);
    std::sort(timed.begin(), timed.end());
    std::cout << "smearres: " << TRUTH_TRK_SMEAR << "\n";
    for (int idx : timed) std::cout << "smeared: " << idx << "," << ideal.times.at(idx) << "\n";
    std::cout << "ttruthvtx: " << ideal.tTruthVtx << "\n";
    clusters = idealisedClusters(tracks, &branch, ideal);
    const int bi = closestToTruthCluster(clusters, branch.truthVtxTime[0]);
    if (bi >= 0)
      std::cout << "tclosest: " << clusters[bi].values[0] << " " << clusters[bi].sigmas[0] << "\n";
  } else {
    clusters = clusterTracksInTime(tracks, &branch, 3.0,
                                   useSmearTimes, useValidTimesOnly, 30.0,
                                   method, useZ0);
  }

  // ===== REFINED STAGE — comment out this block to revert to plain cone output =====
  // Two-pass timing refinement: find the TRKPTZ winner among the 3σ cone clusters,
  // then recompute its time using only tracks within DIST_CUT_REFINE σ of the centroid.
  // {
  //   if (!clusters.empty()) {
  //     auto bestIt = std::max_element(clusters.begin(), clusters.end(),
  //         [](const Cluster& a, const Cluster& b) {
  //           return a.scores.at(Score::TRKPTZ.id) < b.scores.at(Score::TRKPTZ.id);
  //         });
  //     *bestIt = refineClusterTiming(*bestIt, &branch, DIST_CUT_REFINE);
  //   }
  // }
  // ===== END REFINED STAGE =====

  for (int j = 0; j < clusters.size(); j++) {
    auto cluster = clusters.at(j);
    std::cout << "---------\n";
    std::cout << "t: " << cluster.values.at(0) << "\n";
    if (cluster.values.size() > 1) std::cout << "z: " << cluster.values.at(1) << "\n";
    // Idealised-world clusters never go through updateScores, so they carry
    // no TRKPTZ/WAVeS score; event_display.py recomputes one when absent.
    if (cluster.scores.count(Score::TRKPTZ.id))
      std::cout << "score_trkptz: " << cluster.scores.at(Score::TRKPTZ.id) << "\n";
    if (cluster.scores.count(Score::WAVES.id))
      std::cout << "score_waves: " << cluster.scores.at(Score::WAVES.id) << "\n";
    for (int i=0; i < cluster.trackIndices.size(); i++) {
      std::cout << cluster.trackIndices[i] << "," << cluster.allTimes[i] << "\n";

    }
    std::cout << "passes? " << cluster.passEfficiency(&branch) << std::endl;
    std::cout << "---------\n";
  }
}
