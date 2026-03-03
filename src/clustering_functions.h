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
//     7. clusterTracksInTime      — top-level clustering entry point
//     8. chooseHGTDCluster        — select cluster closest to reco vertex time
//     9. chooseHGTDSortCluster    — select highest-BDT cluster (HGTD_SORT score)
//    10. chooseCluster (all scores) — select best cluster for every Score
//    11. chooseCluster (single score) — select best cluster for one Score
// ---------------------------------------------------------------------------

#include <cmath>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "ml_model.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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
    const Cluster& a, const Cluster& b
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
    const std::unordered_map<int, double>& smearedTimesMap, // map of smeared times
    const std::unordered_map<int, double>& smearedTimeResMap, // map of resolution for smeared times if used
    bool checkTimeValid,                          // if true, apply branch->track_time_valid check
    bool usez0
  ) -> std::vector<Cluster> {
    std::vector<Cluster> simpleClusters;
    simpleClusters.reserve(trackIndices.size()); // avoid repeated reallocation
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
	Cluster newCluster = mergeClusters(collection->at(i0), collection->at(j0));
	// j0 > i0 always (inner loop starts at j = i+1).
	// Remove the higher index first so i0 stays valid, both in O(1).
	auto sz = collection->size();
	if ((size_t)j0 != sz - 1) (*collection)[j0] = std::move((*collection)[sz - 1]);
	collection->pop_back();
	sz = collection->size();
	if ((size_t)i0 != sz - 1) (*collection)[i0] = std::move((*collection)[sz - 1]);
	collection->pop_back();
	collection->push_back(std::move(newCluster));
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

      // find merge candidates — use a const-ref to avoid copying the seed
      // cluster before we know whether any merge will happen.
      const Cluster& seedRef = (*collection)[i0];
      std::vector<int> toMerge;
      for (int j = 0; j < (int)collection->size(); ++j) {
	if (j == i0) { continue; }
	double clusterDistance = getDistanceBetweenClusters(seedRef, (*collection)[j]);
	if (clusterDistance < distCut) {
	  toMerge.push_back(j);
	}
      }

      if (DEBUG) std::cout << "Found " << toMerge.size() << " merge candidates\n";

      if (!toMerge.empty()) {
	Cluster newCluster = (*collection)[i0]; // copy only when we know we need it
	for (int idx : toMerge)
	  newCluster = mergeClusters(newCluster, collection->at(idx));

	newCluster.wasMerged = true;
	toMerge.push_back(i0);
	// Remove all consumed clusters in descending-index order via O(1) swap+pop_back.
	// Descending order ensures each swap targets a position that hasn't been
	// vacated yet by an earlier removal.
	std::sort(toMerge.begin(), toMerge.end(), std::greater<int>());
	for (int idx : toMerge) {
	  auto last = collection->size() - 1;
	  if ((size_t)idx != last) (*collection)[idx] = std::move((*collection)[last]);
	  collection->pop_back();
	}
	collection->push_back(std::move(newCluster));
      } else {
	// Seed has no merge candidates: mark it done in-place and rotate it to
	// the back so the seed-search loop skips it (wasMerged == true).
	(*collection)[i0].wasMerged = true;
	if ((size_t)i0 != collection->size() - 1)
	  std::swap((*collection)[i0], collection->back());
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
     bool useSmearedTimes,        // For Ideal Scenarios
     bool checkTimeValid,         // For Ideal Efficiency Scenario
     double smearRes,             // Fixed resolution for smeared times
     bool useCone,                // Whether or not to use improved cone clustering
     bool usez0,                  // Use z info in clustering as well
     bool sortTracks = false,      // Before clustering, sort tracks by pT
     bool calcPurityFlag = false  // only compute purity when TEST_MISCL is active
  ) -> std::vector<Cluster> {
    if (trackIndices.empty() && DEBUG)
      std::cout << "EMPTY!!!!\n";

    std::vector<int> sortedIndices;
    if (sortTracks) {
      sortedIndices = trackIndices;
      std::sort(sortedIndices.begin(), sortedIndices.end(), [&](int a, int b) {
        return branch->trackPt[a] > branch->trackPt[b];
      });
    }
    const std::vector<int>& indices = sortTracks ? sortedIndices : trackIndices;

    std::unordered_map<int,double> smearedTimesMap;
    std::unordered_map<int,double> smearedTimeResMap;

    if (useSmearedTimes) {
      for (auto idx: indices) {
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
      makeSimpleClusters(indices, branch, useSmearedTimes,
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
    for (const Cluster& cluster: collection) {
      double thisDiff = std::abs(cluster.values.at(0) - recoTime);
      if (thisDiff < minDiff) {
	minDiff = thisDiff;
	minCluster = cluster;
      }
    }
    return minCluster;
  }
  
  // ---------------------------------------------------------------------------
  // 9. chooseHGTDSortCluster
  //   Evaluates each cluster with the ATLAS HGTD TMVA BDT and returns the
  //   cluster with the highest BDT output.  The BDT score is stored in
  //   cluster.scores[Score::HGTD_SORT] so the caller can apply the 0.3
  //   confidence threshold.
  //
  //   TMVA variables (order must match weights XML exactly):
  //     m_delta_z          — cluster z minus primary vertex z (mm)
  //     m_z_sigma          — cluster z uncertainty (mm)
  //     m_q_over_p         — precision-weighted charge-over-momentum
  //     m_q_over_p_sigma   — q/p uncertainty
  //     m_delta_z_resunits — delta_z in cluster z-resolution units
  //     m_cluster_sumpt2   — sum of constituent track pT² (GeV²)
  //     m_d0               — precision-weighted transverse impact parameter (mm)
  //     m_d0_sigma         — d0 uncertainty (mm)
  //
  //   The TMVA::Reader is initialised once (static) — safe for the serial
  //   TTreeReader event loop in clustering_dt.cxx.
  // ---------------------------------------------------------------------------
  auto chooseHGTDSortCluster(
    const std::vector<Cluster>& collection,
    BranchPointerWrapper* branch
  ) -> Cluster {
    if (collection.empty()) return Cluster{};

    // Static floats: TMVA Reader holds float* pointers to these, so they must
    // have stable addresses.  Initialised to 0 by the C++ standard.
    static float m_delta_z, m_z_sigma, m_q_over_p, m_q_over_p_sigma;
    static float m_delta_z_resunits, m_cluster_sumpt2, m_d0, m_d0_sigma;

    // Initialise TMVA Reader once; variable order must match the weights XML.
    static TMVA::Reader* reader = []() {
      TMVA::Tools::Instance();
      TMVA::Reader* r = new TMVA::Reader("!Color:!Silent");
      r->AddVariable("m_delta_z",          &m_delta_z);
      r->AddVariable("m_z_sigma",          &m_z_sigma);
      r->AddVariable("m_q_over_p",         &m_q_over_p);
      r->AddVariable("m_q_over_p_sigma",   &m_q_over_p_sigma);
      r->AddVariable("m_delta_z_resunits", &m_delta_z_resunits);
      r->AddVariable("m_cluster_sumpt2",   &m_cluster_sumpt2);
      r->AddVariable("m_d0",               &m_d0);
      r->AddVariable("m_d0_sigma",         &m_d0_sigma);
      r->BookMVA("BDT", "../share/models/TMVA.VBFinv.mu200.Step3p1.8var.weights.xml");
      std::cout << "✓ HGTD BDT loaded (one-time initialization)" << std::endl;
      return r;
    }();

    Cluster best;
    double bestScore = -1e50;
    float refVtxZ = branch->recoVtxZ[0];

    for (const Cluster& cluster : collection) {
      float znum = 0., zden = 0., dnum = 0., dden = 0., qpnum = 0., qpden = 0.;
      float sumpt2 = 0.;

      for (auto trk : cluster.trackIndices) {
        float trkZ   = branch->trackZ0[trk],  trkVarZ = branch->trackVarZ0[trk];
        float trkD   = branch->trackD0[trk],  trkVarD = branch->trackVarD0[trk];
        float trkQ   = branch->trackQP[trk],  trkVarQ = branch->trackVarQp[trk];
        float trkPt  = branch->trackPt[trk];
        znum  += trkZ / trkVarZ;  zden  += 1.f / trkVarZ;
        dnum  += trkD / trkVarD;  dden  += 1.f / trkVarD;
        qpnum += trkQ / trkVarQ;  qpden += 1.f / trkVarQ;
        sumpt2 += trkPt * trkPt;
      }

      float clusterZ        = znum / zden;
      float clusterZSigma   = 1.f / std::sqrt(zden);
      float deltaZ          = clusterZ - refVtxZ;

      m_delta_z          = deltaZ;
      m_z_sigma          = clusterZSigma;
      m_q_over_p         = qpnum / qpden;
      m_q_over_p_sigma   = 1.f / std::sqrt(qpden);
      m_delta_z_resunits = deltaZ / clusterZSigma;
      m_cluster_sumpt2   = sumpt2;
      m_d0               = dnum / dden;
      m_d0_sigma         = 1.f / std::sqrt(dden);

      double bdtScore = reader->EvaluateMVA("BDT");
      if (bdtScore > bestScore) {
        bestScore = bdtScore;
        best = cluster;
      }
    }

    best.scores[Score::HGTD_SORT] = bestScore;
    return best;
  }

  // ---------------------------------------------------------------------------
  // 10. chooseCluster  [all-scores overload]
  //   Iterates over every Score in SCORE_REGISTRY and selects the highest-
  //   scoring cluster for each, returning a map keyed on Score.
  //   Special handling:
  //     HGTD / HGTD_SORT — skipped; each uses its own dedicated collection
  //                         (handled by chooseHGTDCluster / chooseHGTDSortCluster).
  //     FILTJET — skipped; uses a pre-filtered jet-cone track collection and
  //               is chosen by the single-score overload.
  //     PASS    — uses passEfficiency() as the selection criterion instead
  //               of a numeric score.
  //   Intended for the main analysis path where all active scores are needed
  //   simultaneously, avoiding repeated iteration over the cluster collection.
  // ---------------------------------------------------------------------------
  auto chooseCluster(
    const std::vector<Cluster>& collection,
    BranchPointerWrapper *branch
  ) -> std::map<Score, Cluster> {
    std::map<Score,Cluster> output;
    if (DEBUG) std::cout << "Choosing score\n";

    for (Score score: SCORE_REGISTRY) {
      if (DEBUG)
        std::cout << "SCORE: " << toString(score) << '\n';
      if (score.usesOwnCollection)
	continue;

      // FILTJET uses its own pre-filtered collection; chosen by single-score overload
      if (score == Score::FILTJET)
	continue;

      if (score == Score::PASS) {
	if (DEBUG) std::cout << "Choosing pass score\n";
        for (const Cluster& cluster : collection)
          if (cluster.passEfficiency(branch))
	    output[score] = cluster;
	continue;
      }

      // Track the best cluster via pointer — only copy the winner once at the end.
      const Cluster* best = &collection[0];
      double maxScore = best->scores.at(score);

      for (const Cluster& cluster: collection) {
	double compScore = cluster.scores.at(score);
	if (compScore > maxScore) {
	  maxScore = compScore;
	  best = &cluster;
	}
      }
      output[score] = *best;
    }
    return output;
  }

  // ---------------------------------------------------------------------------
  // 11. chooseCluster  [single-score overload]
  //   Selects the highest-scoring cluster for a single specified Score by
  //   linear scan.  Simpler than the all-scores overload: no PASS special
  //   case.  Used to choose the best cluster from the pre-filtered FILTJET
  //   collection where only TRKPTZ scoring is needed.
  // ---------------------------------------------------------------------------
  auto chooseCluster(
    const std::vector<Cluster>& collection,
    Score score
  ) -> Cluster {
    if (DEBUG) std::cout << "Choosing score\n";

    // if (score >= 3) {
      // std::cout << "DONT CALL LIKE THIS\n";
      // return {};
    // }

    Cluster output = collection[0]; // final time we are giving to the user
    double maxScore = output.scores.at(score);
    
    for (const Cluster& cluster: collection) {
      double compScore = cluster.scores.at(score);
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
