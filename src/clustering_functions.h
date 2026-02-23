#ifndef CLUSTERING_FUNCTIONS_H
#define CLUSTERING_FUNCTIONS_H

// ---------------------------------------------------------------------------
// clustering_functions.h
//   Core clustering algorithms and cluster-selection logic.  All functions
//   live inside the MyUtl namespace.  Sections:
//     1. getSmearedTrackTime   — smear a track time with a given resolution
//     2. getDistanceBetweenClusters — Mahalanobis-like distance metric
//     3. mergeClusters         — precision-weighted merge of two clusters
//     4. makeSimpleClusters    — seed one cluster per track
//     5. doSimultaneousClustering — agglomerative (minimum-distance) merge
//     6. doConeClustering      — seed-and-cone merge
//     7. clusterTracksInTime   — top-level clustering entry point
//     8. chooseHGTDCluster     — select cluster closest to reco vertex time
//     9. chooseCluster (all scores) — select best cluster for every Score
//    10. chooseCluster (single score) — select best cluster for one Score
// ---------------------------------------------------------------------------

#include <cmath>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "ml_model.h"

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // 1. getSmearedTrackTime
  //   Returns a Gaussian-smeared time for track idx using resolution res.
  //   The true underlying time is taken from (in priority order):
  //     a) the associated truth particle's production time, if a particle link
  //        exists;
  //     b) the associated truth vertex time, if a vertex link exists;
  //     c) a random draw from Gaus(truthVtxTime[0], PILEUP_SMEAR) to model
  //        an unlinked pileup track.
  //   Used when useSmearedTimes == true in makeSimpleClusters.
  // ---------------------------------------------------------------------------
  auto getSmearedTrackTime(
    int idx, double res, BranchPointerWrapper *branch
  ) -> double {
    int particleIdx = branch->trackToParticle[idx]; // this is anything
    int truthvtxIdx = branch->trackToTruthvtx[idx]; // this is anything
    double tPart = NAN;
    double smearRes = res;
    if (particleIdx == -1) { // no particle link
      if (truthvtxIdx != -1) { // some valid truth link, smear that vertex time
        tPart = branch->truthVtxTime[truthvtxIdx];
      } else { // some random pileup, assume
        tPart = gRandom->Gaus(branch->truthVtxTime[0],PILEUP_SMEAR);
      }
    } else {
      tPart = branch->particleT[particleIdx];
    }
    return gRandom->Gaus(tPart,smearRes);
  }
  
  // ---------------------------------------------------------------------------
  // 2. getDistanceBetweenClusters
  //   Computes the Mahalanobis-like distance between two clusters:
  //     d = sqrt( Σ_i  ((a.values[i] - b.values[i]) /
  //                      sqrt(a.sigmas[i]² + b.sigmas[i]²))² )
  //   Each dimension contributes a normalised squared difference; the combined
  //   uncertainty in the denominator accounts for both clusters' uncertainties.
  //   Used as the merge criterion in both clustering algorithms.
  // ---------------------------------------------------------------------------
  auto getDistanceBetweenClusters(
    const Cluster& a, const Cluster& b
  ) -> double {
    // Work directly on references — no copies, no intermediate vector.
    // Accumulate d^2 in a single pass.
    const size_t N = a.values.size();
    double dsqr = 0.0;
    for (size_t i = 0; i < N; ++i) {
      double diff  = a.values[i] - b.values[i];
      double denom = std::sqrt(a.sigmas[i]*a.sigmas[i] + b.sigmas[i]*b.sigmas[i]);
      double d = diff / denom;
      dsqr += d * d;
    }
    return std::sqrt(dsqr);
  }

  // ---------------------------------------------------------------------------
  // 3. mergeClusters
  //   Returns a new Cluster formed by precision-weighted averaging of clusters
  //   a and b in each dimension:
  //     new_value = (v1/σ1² + v2/σ2²) / (1/σ1² + 1/σ2²)
  //     new_sigma = σ1·σ2 / sqrt(σ1² + σ2²)
  //   allTimes and trackIndices are concatenated.  Scores are summed (additive
  //   combination) so that cluster-level scores reflect all constituent tracks.
  //   wasMerged is set true on the result so the clustering loop can identify
  //   freshly merged clusters and reset the flag before the next pass.
  // ---------------------------------------------------------------------------
  auto mergeClusters(
    Cluster a, Cluster b
  ) -> Cluster {
    Cluster mergedCluster;
    mergedCluster.wasMerged = true;
  
    for (size_t i = 0; i < a.values.size(); i++) {
      double v1 = a.values.at(i);
      double v2 = b.values.at(i);
      double s1 = a.sigmas.at(i);
      double s2 = b.sigmas.at(i);

      double newClusterValue =
        (v1 / (s1 * s1) + v2 / (s2 * s2)) / (1. / (s1 * s1) + 1. / (s2 * s2));
      double newClusterSigma =
	sqrt((s2 * s1) * (s2 * s1) / (s2 * s2 + s1 * s1));
      mergedCluster.values.push_back(newClusterValue);
      mergedCluster.sigmas.push_back(newClusterSigma);
    }

    mergedCluster.allTimes.reserve(a.allTimes.size() + b.allTimes.size());
    mergedCluster.allTimes.insert(mergedCluster.allTimes.end(), a.allTimes.begin(), a.allTimes.end());
    mergedCluster.allTimes.insert(mergedCluster.allTimes.end(), b.allTimes.begin(), b.allTimes.end());

    mergedCluster.trackIndices.reserve(a.trackIndices.size() + b.trackIndices.size());
    mergedCluster.trackIndices.insert(mergedCluster.trackIndices.end(), a.trackIndices.begin(), a.trackIndices.end());
    mergedCluster.trackIndices.insert(mergedCluster.trackIndices.end(), b.trackIndices.begin(), b.trackIndices.end());
    
    mergedCluster.nConstituents = a.nConstituents+b.nConstituents;
    mergedCluster.maxPtCluster = a.maxPtCluster || b.maxPtCluster;

    for (auto [k,v] : a.scores) {
      mergedCluster.scores[k] = v;
    }

    for (auto [k,v] : b.scores) {
      if (mergedCluster.scores.count(k))
	mergedCluster.scores[k] += v;
      else
	mergedCluster.scores[k] = v;
    }

    return mergedCluster;
  }

  // ---------------------------------------------------------------------------
  // 4. makeSimpleClusters
  //   Seeds one single-track Cluster per qualifying track in trackIndices.
  //   A track is included only if checkTimeValid is false or if its
  //   Track_hasValidTime flag is set.  When useSmearedTimes is true the
  //   pre-computed smearedTimesMap and smearedTimeResMap are used in place
  //   of the raw branch times.  When usez0 is true both time and z₀ are
  //   stored as cluster dimensions; otherwise only time is used.
  //   Initial scores HGTD (0) and TRKPT (track_pT) are set inline.
  // ---------------------------------------------------------------------------
  auto makeSimpleClusters(
    const std::vector<int>& trackIndices,
    BranchPointerWrapper *branch,
    bool useSmearedTimes,                         // if true, use smeared_times map
    const std::map<int, double>& smearedTimesMap, // map of smeared times
    const std::map<int, double>& smearedTimeResMap, // map of resolution for smeared times if used
    bool checkTimeValid,                          // if true, apply branch->track_time_valid check
    bool usez0
  ) -> std::vector<Cluster> {
    std::vector<Cluster> simpleClusters;
    for (auto idx: trackIndices) {
      if (!checkTimeValid || branch->trackTimeValid[idx] == 1) {
	double trkTime = NAN;
	double trkRes = NAN; 
	if (useSmearedTimes) {
	  trkTime = smearedTimesMap.at(idx);
	  trkRes  = smearedTimeResMap.at(idx);
	} else {
	  trkTime = branch->trackTime[idx];   
	  trkRes  = branch->trackTimeRes[idx];
	}

	double trkPt = branch->trackPt[idx];
	double trkZ0 = branch->trackZ0[idx];
	double trkResZ0 = std::sqrt(branch->trackVarZ0[idx]);

        // Initialise scores inline — avoids constructing a temporary map
        // and two separate tree-node insertions per track.
        Cluster cluster;
        if (usez0)
          cluster = {{trkTime,trkZ0}, {trkRes,trkResZ0}, {trkTime}, {idx},
                     {{Score::HGTD, 0.0}, {Score::TRKPT, trkPt}}};
        else
          cluster = {{trkTime}, {trkRes}, {trkTime}, {idx},
                     {{Score::HGTD, 0.0}, {Score::TRKPT, trkPt}}};
        simpleClusters.push_back(std::move(cluster));
      }
    }
    return simpleClusters;
  }

  // ---------------------------------------------------------------------------
  // 5. doSimultaneousClustering
  //   Agglomerative (bottom-up) clustering: at each iteration the globally
  //   closest pair of unmerged clusters is found and merged if their distance
  //   is below distCut.  The loop terminates when no pair is closer than
  //   distCut.  The wasMerged flag is used to skip freshly merged clusters in
  //   the inner distance search and is reset between passes.
  // ---------------------------------------------------------------------------
  void doSimultaneousClustering(
    std::vector<Cluster> *collection, double distCut
  ) {
    double distance = 1.e30;
    while (collection->size() > 1) {
      int i0 = 0;
      int j0 = 0;

      distance = getDistanceBetweenClusters(collection->at(0), collection->at(1));

      for (int i = 0; i < collection->size(); i++) {
	for (int j = i + 1; j < collection->size(); j++) {
	  const Cluster& a = collection->at(i);
	  const Cluster& b = collection->at(j);

	  if (a.wasMerged or b.wasMerged)
	    continue;
	
	  double currentDistance =
	    getDistanceBetweenClusters(a, b);
	  if (currentDistance <= distance) {
	    distance = currentDistance;
	    i0 = i;
	    j0 = j;
	  }
	} // j loop
      } // i loop

      // fuse closest two vertices
      if (distance < distCut and i0 != j0) {
	Cluster newCluster = mergeClusters(collection->at(i0),collection->at(j0));
	collection->erase(collection->begin()+j0);
	if (i0 < j0)
	  collection->erase(collection->begin()+i0);
	else
	  collection->erase(collection->begin()+(i0-1));
	collection->push_back(newCluster);
      } else {
	if (std::find_if(collection->begin(), collection->end(),
			 [](const Cluster& a) {return a.wasMerged;}) != collection->end()) {
	  for (auto & idx : *collection)
	    idx.wasMerged = false;
	} else
	  break;
      }
    } // While
  }

  // ---------------------------------------------------------------------------
  // 6. doConeClustering
  //   Seed-and-cone clustering: the unmerged cluster with the highest TRKPT
  //   score is selected as the seed, all other clusters within distCut of
  //   the seed are absorbed into it, and the merged cluster is marked as
  //   processed.  Repeats until every cluster has been merged or assigned.
  //   Produces a result biased toward the highest-pT cluster being the core.
  // ---------------------------------------------------------------------------
  void doConeClustering(
    std::vector<Cluster> *collection, double distCut
  ) {
    if (DEBUG) { std::cout << "entering while loop\n"; }
    while (collection->size() > 1) {
      // find seed
      int i0 = -1;
      double seedValue = -1.;
      for (int i=0; i < collection->size(); ++i) {
        if (collection->at(i).wasMerged) { continue; }
	const Cluster& check = collection->at(i);
	double checkValue = check.scores.at(Score::TRKPT);
	
	if (seedValue < checkValue) {
	  seedValue = checkValue;
	  i0 = i;
	}
      }
      if (i0 == -1) { break; }
      
      if (DEBUG) std::cout << "Found Seed\n";

      // find merge candidates
      Cluster seed = collection->at(i0);
      std::vector<int> toMerge;
      for (int j = 0; j < collection->size(); ++j) {
	if (j == i0) { continue; }
	double clusterDistance = getDistanceBetweenClusters(seed, collection->at(j));
	if (clusterDistance < distCut) {
	  toMerge.push_back(j);
	}
      }
      
      if (DEBUG) std::cout << "Found " << toMerge.size() << " merge candidates\n";

      if (!toMerge.empty()) {
	Cluster newCluster = seed;
	for (int idx : toMerge)
	  newCluster = mergeClusters(newCluster, collection->at(idx));

	newCluster.wasMerged = true;
	toMerge.push_back(i0);
	std::sort(toMerge.begin(), toMerge.end(), std::greater<int>());
	for (int idx : toMerge) {
	  collection->erase(collection->begin() + idx);
	}
	collection->push_back(newCluster);
      } else {
	Cluster seedOnly = collection->at(i0);
	seedOnly.wasMerged = true;
	collection->erase(collection->begin() + i0);
	collection->push_back(seedOnly);
      }
    } // While
    if (DEBUG) std::cout << "Finished Clustering\n";
  }

  // ---------------------------------------------------------------------------
  // 7. clusterTracksInTime
  //   Top-level clustering entry point called from event_processing.h.
  //   Orchestrates the full pipeline:
  //     a) If useSmearedTimes, generate a smeared time for every track (using
  //        getSmearedTrackTime).  The resolution is scaled by 1/sqrt(nHits)
  //        for tracks with valid HGTD hits.
  //     b) Seed one cluster per track with makeSimpleClusters.
  //     c) Cluster with doConeClustering (useCone == true) or
  //        doSimultaneousClustering (useCone == false).
  //     d) Load the ML model once via static initialisation (first call only).
  //     e) Optionally compute cluster purity (calcPurityFlag) and call
  //        updateScores on every cluster to fill the derived score map.
  // ---------------------------------------------------------------------------
  auto clusterTracksInTime(
     const std::vector<int>& trackIndices,
     BranchPointerWrapper *branch,
     double distanceCut,          // Cone size/distance cut
     bool useSmearedTimes,        // for ideal cases
     bool checkTimeValid,         // for ideal efficiency
     double smearRes,             // fixed resolution for smeared times
     bool useCone,                // whether or not to use improved cone clustering
     bool usez0,                  // Use z info in clustering as well
     bool calcPurityFlag = false  // only compute purity when TEST_MISCL is active
  ) -> std::vector<Cluster> {
    if (trackIndices.empty() && DEBUG)
      std::cout << "EMPTY!!!!\n";

    std::map<int,double> smearedTimesMap;
    std::map<int,double> smearedTimeResMap;

    if (useSmearedTimes) {
      for (auto idx: trackIndices) {
	double res = smearRes;
	bool timeQuery = branch->trackTimeValid[idx] == 1
	  and branch->trackHgtdHits[idx] > 0;
	if (timeQuery) { res /= std::sqrt(branch->trackHgtdHits[idx]); }
	if (!checkTimeValid || branch->trackTimeValid[idx] == 1) {
	  double smearedTime = getSmearedTrackTime(idx, res, branch);
	  smearedTimesMap.emplace(idx, smearedTime);
	  smearedTimeResMap.emplace(idx, res);
	}
      }
    }

    std::vector<Cluster> collection =
      makeSimpleClusters(trackIndices, branch, useSmearedTimes,
			 smearedTimesMap, smearedTimeResMap, checkTimeValid, usez0);
    
    if (not useCone)
      doSimultaneousClustering(&collection, distanceCut);
    else
      doConeClustering(&collection, distanceCut);

    // Load ML model only once using static initialization (lazy evaluation)
    static MLModel mlModel = []() {
        MLModel m;
        m.load_weights("../share/models/model_weights.json");
        std::cout << "✓ ML model loaded (one-time initialization)" << std::endl;
        return m;
    }();

    for (Cluster& cluster: collection) {
      if (calcPurityFlag)
        cluster.calcPurity(branch);  // only needed when TEST_MISCL is active
      cluster.updateScores(branch, &mlModel);
    }
    
    if (DEBUG) std::cout << "Finished Clustering\n";
    return collection;
  }

  // ---------------------------------------------------------------------------
  // 8. chooseHGTDCluster
  //   Selects the cluster whose weighted-mean time is closest (in absolute
  //   difference) to the reconstructed primary vertex time (recoVtxTime[0]).
  //   Used exclusively for the HGTD score, which compares directly against
  //   the reco-vertex timing rather than using a score ranking.
  // ---------------------------------------------------------------------------
  auto chooseHGTDCluster(
    const std::vector<Cluster>& collection,
    BranchPointerWrapper *branch
  ) -> Cluster {
    double minDiff = 1e50;
    double recoTime = branch->recoVtxTime[0];
    Cluster minCluster;
    for (Cluster cluster: collection) {
      double thisDiff = std::abs(cluster.values.at(0) - recoTime);
      if (thisDiff < minDiff) {
	minDiff = thisDiff;
	minCluster = cluster;
      }
    }
    return minCluster;
  }
  
  // ---------------------------------------------------------------------------
  // 9. chooseCluster  [all-scores overload]
  //   Iterates over every Score in ENUM_VEC and selects the highest-scoring
  //   cluster for each, returning a map keyed on Score.  Special handling:
  //     HGTD     — skipped here; handled separately by chooseHGTDCluster.
  //     FILTJET / FILT60 / FILT90 — skipped; use their own pre-filtered
  //                collection and are chosen by the single-score overload.
  //     JUST60 / JUST90 — a synthetic cluster with only a calorimeter time
  //                is inserted directly (no track-based cluster chosen).
  //     PASS    — uses passEfficiency() as the selection criterion instead
  //                of a numeric score.
  //     CALO90 / CALO60 — additionally require the cluster time to be within
  //                2σ of a smeared calorimeter time.
  //   Intended for the main analysis path where all active scores are needed
  //   simultaneously, avoiding repeated iteration over the cluster collection.
  // ---------------------------------------------------------------------------
  auto chooseCluster(
    std::vector<Cluster>& collection,
    BranchPointerWrapper *branch
  ) -> std::map<Score, Cluster> {
    std::map<Score,Cluster> output;
    if (DEBUG) std::cout << "Choosing score\n";

    double caloTime90, caloTime60;
    if (std::find(ENUM_VEC.begin(), ENUM_VEC.end(), CALO90) != ENUM_VEC.end())      
      caloTime90 = gRandom->Gaus(branch->truthVtxTime[0],90);

    if (std::find(ENUM_VEC.begin(), ENUM_VEC.end(), CALO60) != ENUM_VEC.end())      
      caloTime60 = gRandom->Gaus(branch->truthVtxTime[0],60);
    
    for (Score score: ENUM_VEC) {
      if (DEBUG)
        std::cout << "SCORE: " << toString(score) << '\n';
      if (score == HGTD)
	continue;

      // skip filtered tracks, they use a different collection
      if (score == FILTJET or score == FILT60 or score == FILT90)
	continue;
      if (score == Score::JUST60) {
	output[score] = {{caloTime60}, {60.0}, {}, {-1}, {}};
	continue;
      }
      
      if (score == Score::JUST90) {
	output[score] = {{caloTime90}, {90.0}, {}, {-1}, {}};
	continue;
      }

      if (score == Score::PASS) {        
	if (DEBUG) std::cout << "Choosing pass score\n";
        for (Cluster &cluster : collection)
          if (cluster.passEfficiency(branch))            
	    output[score] = cluster;
	continue;
      }

      output[score] = collection[0]; // final time we are giving to the user
      double maxScore = output[score].scores.at(score);
      // std::cout << ""<< maxScore << '\n';
    
      for (Cluster& cluster: collection) {
	double compScore = cluster.scores[score];
	bool query = compScore > maxScore;

	if (score == Score::CALO90) {
	  double caloDiff = std::abs(cluster.values.at(0)-caloTime90);
	  double caloSigma = std::hypot(cluster.sigmas.at(0), 90);
	  double nsigma = caloDiff/caloSigma;
	  if (nsigma > 2.0)
	    continue;
	}

	if (score == Score::CALO60) {
	  double caloDiff = std::abs(cluster.values.at(0)-caloTime60);
	  double caloSigma = std::hypot(cluster.sigmas.at(0), 60);
	  double nsigma = caloDiff/caloSigma;
	  if (nsigma > 2.0)
	    continue;
	}
	
	if (query) {
	  maxScore = compScore;
	  output[score] = cluster;
	}
      }
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 10. chooseCluster  [single-score overload]
  //   Selects the highest-scoring cluster for a single specified Score by
  //   linear scan.  Simpler than the all-scores overload: no calorimeter-time
  //   gating, no PASS/JUST special cases.  Used to choose the best cluster
  //   from pre-filtered collections (FILT60, FILT90, FILTJET) where only one
  //   score (TRKPTZ) is needed.
  // ---------------------------------------------------------------------------
  auto chooseCluster(
    std::vector<Cluster>& collection,
    Score score
  ) -> Cluster {
    if (DEBUG) std::cout << "Choosing score\n";

    // if (score >= 3) {
      // std::cout << "DONT CALL LIKE THIS\n";
      // return {};
    // }

    Cluster output = collection[0]; // final time we are giving to the user
    double maxScore = output.scores.at(score);
    
    for (Cluster& cluster: collection) {
      double compScore = cluster.scores[score];
      bool query = compScore > maxScore;
	
      if (query) {
	maxScore = compScore;
	output = cluster;
      }
    }
    return output;
  }
}
#endif // CLUSTERING_FUNCTIONS_H
