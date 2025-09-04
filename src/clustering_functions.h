#ifndef CLUSTERING_FUNCTIONS_H
#define CLUSTERING_FUNCTIONS_H

#include <cmath>

#include "clustering_constants.h"
#include "clustering_structs.h"

namespace MyUtl {

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
  
  auto getDistanceBetweenClusters(
    const Cluster& a, const Cluster& b
  ) -> double {
    std::vector<double> aVal = a.values;
    std::vector<double> aSig = a.sigmas;
    std::vector<double> bVal = b.values;
    std::vector<double> bSig = b.sigmas;

    std::vector<double> metric = {1.,-1.};

    if (aVal.size() != bVal.size()) { std::cout << "Uh ohhh\n"; }

    std::vector<double> distances(aVal.size());
    for (int i=0; i < aVal.size(); i++) {
    
      double diffI = (aVal.at(i) - bVal.at(i));
      double denomI = sqrt(aSig.at(i)*aSig.at(i) + bSig.at(i)*bSig.at(i));
      distances.at(i) = diffI/denomI;
    }

    double dsqr = 0.0;
    for (int i=0; i < distances.size(); i++){
      auto d = distances[i];
      dsqr += metric.at(i)*d*d;
    }
    return sqrt(dsqr);
  }

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

    for (auto time: a.allTimes)
      mergedCluster.allTimes.push_back(time);
    for (auto time: b.allTimes)
      mergedCluster.allTimes.push_back(time);
      
    for (auto idx: a.trackIndices)
      mergedCluster.trackIndices.push_back(idx);
    for (auto idx: b.trackIndices)
      mergedCluster.trackIndices.push_back(idx);
    
    mergedCluster.nConstituents = a.nConstituents+b.nConstituents;
    mergedCluster.maxPtCluster = a.maxPtCluster || b.maxPtCluster;

    for (Score score : ENUM_VEC) {
      mergedCluster.scores[score] = a.scores[score]+b.scores[score];
    }

    return mergedCluster;
  }

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
	double trkResZ0 = branch->trackVarZ0[idx];
	double zSignificance = std::abs(trkZ0 - branch->recoVtxZ[0]) / trkResZ0;

	std::map<Score,double> scores;
	scores[Score::HGTD] = 0;
	scores[Score::TRKPT] = trkPt;

	Cluster cluster;
	if (usez0) cluster = {{trkTime,trkZ0}, {trkRes,trkResZ0}, {trkTime}, {idx}, scores};
	else cluster = {{trkTime}, {trkRes}, {trkTime}, {idx}, scores};
	simpleClusters.push_back(cluster);
      }
    }
    return simpleClusters;
  }

  void doSimultaneousClustering(
    std::vector<Cluster> *collection, double distCut
  ) {
    double distance = 1.e30;
    while (collection->size() > 1) {
      // std::cout << "entering while loop" << std::endl;

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
      // break;
    } // While
  }

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

  auto clusterTracksInTime(
     const std::vector<int>& trackIndices,
     BranchPointerWrapper *branch,
     double distanceCut,    // Cone size/distance cut
     bool useSmearedTimes, // for ideal cases
     bool checkTimeValid,  // for ideal efficiency
     double smearRes,       // fixed resolution for smeared times
     bool useCone,
     bool usez0
  ) -> std::vector<Cluster> {
    if (trackIndices.empty() && DEBUG)
      std::cout << "EMPTY!!!!\n";

    std::map<int,double> smearedTimesMap;
    std::map<int,double> smearedTimeResMap;

    // if use smeared times, generated
    if (useSmearedTimes) {
      for (auto idx: trackIndices) {
	double res = smearRes;
	bool timeQuery = branch->trackTimeValid[idx] == 1 and branch->trackHgtdHits[idx] > 0;
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
        
    for (Cluster& cluster: collection) {
      cluster.calcPurity(branch);
      cluster.updateScores(branch);
    }
    
    if (DEBUG) std::cout << "Finished Clustering\n";
    return collection;
  }

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
  
  auto chooseCluster(
    std::vector<Cluster> collection,
    BranchPointerWrapper *branch
  ) -> std::map<Score, Cluster> {
    std::map<Score,Cluster> output;
    if (DEBUG) std::cout << "Choosing score\n";

    double caloTime90, caloTime60;
    if(std::find(ENUM_VEC.begin(), ENUM_VEC.end(), CALO90)!=ENUM_VEC.end())
      caloTime90 = gRandom->Gaus(branch->truthVtxTime[0],90);

    if(std::find(ENUM_VEC.begin(), ENUM_VEC.end(), CALO60)!=ENUM_VEC.end())
      caloTime60 = gRandom->Gaus(branch->truthVtxTime[0],60);
      
    for (Score score: ENUM_VEC) {
      if (DEBUG) std::cout << "SCORE: " << toString(score) << '\n';
      if (score == Score::HGTD)
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
	for (Cluster& cluster: collection)
	  if (cluster.passEfficiency(branch))
	    output[score] = cluster;
	continue;
      }

      output[score] = collection[0]; // final time we are giving to the user
      double maxScore = output[score].scores.at(score);
    
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
}
#endif // CLUSTERING_FUNCTIONS_H
