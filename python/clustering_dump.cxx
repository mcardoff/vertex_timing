// clustering_dump.cxx -- the C++ clustering and WAVeS choice of one event, for
// slide figures (python/gp_clustering_figure.py; the Global Performance talk's
// clustering slide, 2026-11-12). Same headers as runHGTD_Clustering.cxx, so the
// clusters, scores and t0 are exactly the analysis's, with the talk's inputs:
// TDR track-time cleaning (CLEAN_TIMES) and the 2.5 sigma z association.
//
//   root -l -b -q 'clustering_dump.cxx("<file>", -1)'              scan: list candidate events
//   root -l -b -q 'clustering_dump.cxx("<file>", <entry>, "<out.json>")'   dump one event
//
// A candidate (scan) is an event where WAVeS chooses the hard scatter's cluster
// (hsClusterOf, the selection fix's identification) and passes, while the plain
// pT sum would choose a pileup cluster, with >= 3 clusters of >= 2 tracks: the
// case that shows what the score is for.
#include "../src/clustering_constants.h"
R__ADD_INCLUDE_PATH(/opt/homebrew/opt/boost/include)
R__ADD_LIBRARY_PATH(/opt/homebrew/opt/boost/lib)
R__LOAD_LIBRARY(libboost_filesystem)

#include "../src/clustering_functions.h"
#include "../src/event_processing.h"
#include <fstream>
#include <iomanip>

using namespace MyUtl;

void clustering_dump(std::string filePath, Long64_t eventNum = -1, std::string out = "") {
  CLEAN_TIMES = true;                       // binds Track_hitInLastLayer
  TChain chain("ntuple");
  chain.Add(filePath.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  auto process = [&](Long64_t entry, bool dump) {
    reader.SetEntry(entry);
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection() || branch.vetoLeptonOverlap()) return;
    if (!branch.passBasicCuts() || !branch.passJetPtCut()) return;
    std::vector<int> tracks = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 2.5);
    tracks.erase(std::remove_if(tracks.begin(), tracks.end(), [&](int i) {
                   return branch.trackTimeValid[i] == 1 && !branch.timeSurvivesCleaning(i); }), tracks.end());
    if (tracks.empty()) return;
    auto clusters = clusterTracksInTime(tracks, &branch, DIST_CUT_CONE, false, true, -1.0,
                                        ClusteringMethod::ITERATIVE, false, false, /*purity=*/true);
    if (clusters.empty()) return;
    const double tHS = branch.truthVtxTime[0];
    int iW = 0, iP = 0, iH = -1, nBig = 0;
    std::vector<double> sumPt(clusters.size(), 0.0);
    for (size_t k = 0; k < clusters.size(); ++k) {
      for (int t : clusters[k].trackIndices) sumPt[k] += branch.trackPt[t];
      if (clusters[k].scores.at(Score::WAVES.id) > clusters[iW].scores.at(Score::WAVES.id)) iW = k;
      if (sumPt[k] > sumPt[iP]) iP = k;
      if (clusters[k].trackIndices.size() >= 2) ++nBig;
    }
    if (const Cluster* h = hsClusterOf(clusters, &branch))
      for (size_t k = 0; k < clusters.size(); ++k) if (&clusters[k] == h) iH = k;
    const double t0W = clusters[iW].calculateTime(Score::WAVES, &branch);
    if (!dump) {
      const bool cand = iH >= 0 && iW == iH && iP != iH && nBig >= 3 && std::abs(t0W - tHS) < PASS_SIGMA;
      if (cand)
        std::cout << "CAND " << entry << " nclus " << clusters.size() << " nbig " << nBig
                  << " hs_ntrk " << clusters[iH].trackIndices.size() << " sumpt_hs " << sumPt[iH]
                  << " sumpt_max " << sumPt[iP] << " waves_ratio "
                  << clusters[iW].scores.at(Score::WAVES.id) / std::max(1e-9, clusters[iP].scores.at(Score::WAVES.id))
                  << "\n";
      return;
    }
    std::ofstream js(out);
    js << std::setprecision(8) << "{\n  \"file\": \"" << filePath << "\", \"entry\": " << entry
       << ", \"t_hs\": " << tHS << ", \"z_pv\": " << branch.recoVtxZ[0]
       << ", \"athena_t0\": " << branch.recoVtxTime[0] << ", \"athena_valid\": " << (int)branch.recoVtxValid[0]
       << ", \"t0_waves\": " << t0W << ", \"i_waves\": " << iW << ", \"i_sumpt\": " << iP << ", \"i_hs\": " << iH
       << ",\n  \"clusters\": [\n";
    for (size_t k = 0; k < clusters.size(); ++k) {
      const auto& c = clusters[k];
      js << (k ? ",\n" : "") << "    {\"t\": " << c.values[0] << ", \"sigma\": " << c.sigmas[0] << ", \"sumpt\": " << sumPt[k]
         << ", \"purity\": " << c.purity << ", \"waves\": " << c.scores.at(Score::WAVES.id) << ", \"tracks\": [";
      for (size_t i = 0; i < c.trackIndices.size(); ++i) {
        const int t = c.trackIndices[i];
        js << (i ? ", " : "") << "{\"idx\": " << t << ", \"t\": " << c.allTimes[i] << ", \"res\": " << branch.trackTimeRes[t]
           << ", \"pt\": " << branch.trackPt[t] << ", \"eta\": " << branch.trackEta[t] << ", \"phi\": " << branch.trackPhi[t]
           << ", \"hs\": " << (branch.trackToTruthvtx[t] == 0 ? 1 : 0) << "}";
      }
      js << "]}";
    }
    js << "\n  ],\n  \"jets\": [";
    bool first = true;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j) || branch.topoJetPt[j] < MIN_JET_PT) continue;
      const double ae = std::abs((double)branch.topoJetEta[j]);
      if (ae < MIN_ABS_ETA_JET || ae > MAX_ABS_ETA_JET) continue;
      js << (first ? "" : ", ") << "{\"pt\": " << branch.topoJetPt[j] << ", \"eta\": " << branch.topoJetEta[j]
         << ", \"phi\": " << branch.topoJetPhi[j] << "}";
      first = false;
    }
    js << "]\n}\n";
    std::cout << "wrote " << out << "\n";
  };

  if (eventNum >= 0) { process(eventNum, true); return; }
  const Long64_t n = reader.GetEntries(true);
  for (Long64_t e = 0; e < n; ++e) process(e, false);
}
