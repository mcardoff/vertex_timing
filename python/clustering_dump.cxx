// clustering_dump.cxx -- the C++ clustering and WAVeS choice of one event, for
// slide figures (python/gp_clustering_figure.py; the Global Performance talk's
// clustering slide, 2026-11-12). Same headers as runHGTD_Clustering.cxx, so the
// clusters, scores and t0 are exactly the analysis's, with the talk's inputs:
// TDR track-time cleaning (CLEAN_TIMES) and the 2.5 sigma z association.
//
//   root -l -b -q 'clustering_dump.cxx("<file>", -1)'              scan: list candidate events
//   root -l -b -q 'clustering_dump.cxx("<file>", <entry>, "<out.json>")'   dump one event
//   root -l -b -q -e '.L clustering_dump.cxx' -e 'clustering_dump_list("<list>", "<outdir>")'
//       dump many: each line of <list> is "<file> <entry> <name>", written to <outdir>/<name>.json
//       (the failure-mode displays, python/gp_failure_display.py)
//
// A candidate (scan) is an event where WAVeS chooses the hard scatter's cluster
// (hsClusterOf, the selection fix's identification) and passes, while the plain
// pT sum would choose a pileup cluster, with >= 3 clusters of >= 2 tracks: the
// case that shows what the score is for.
//
// Per track the dump also carries "injet" (within dR < 0.4 of a qualifying forward
// jet: the tracks WAVeS's in-jet time uses, as Cluster::calculateTime) and
// "mistimed" (a hard-scatter track >= 3 sigma from the truth vertex time:
// util/failure_decomposition's definition); per cluster "t_hscore", its HS tracks'
// inverse-variance mean (failure_decomposition's misclustering oracle); per event
// "t0_waves_fix", the WAVeS t0 after the mistimed HS tracks are dropped, the event
// re-clustered and WAVeS re-chosen (its misassigned-times oracle).
#include "../src/clustering_constants.h"
R__ADD_INCLUDE_PATH(/opt/homebrew/opt/boost/include)
R__ADD_LIBRARY_PATH(/opt/homebrew/opt/boost/lib)
R__LOAD_LIBRARY(libboost_filesystem)

#include "../src/clustering_functions.h"
#include "../src/event_processing.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

using namespace MyUtl;

// the event at the reader's current entry; scan (out empty) prints a candidate line, else writes JSON
static void dumpCurrent(BranchPointerWrapper& branch, const std::string& filePath, Long64_t entry,
                        const std::string& out) {
  const bool dump = !out.empty();
  branch.computeOverlapRemoval();
  if (!branch.passLeptonSelection() || branch.vetoLeptonOverlap()) return;
  if (!branch.passBasicCuts() || !branch.passJetPtCut()) return;
  std::vector<int> tracks = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 2.5);
  int nHsAssoc = 0, nHsNoTime = 0, nHsCleaned = 0;   // the hard scatter's associated forward tracks, and why some lack a time
  std::vector<int> hsUntimed;
  for (int i : tracks) {
    if (branch.trackToTruthvtx[i] != 0) continue;
    ++nHsAssoc;
    if (branch.trackTimeValid[i] != 1) { ++nHsNoTime; hsUntimed.push_back(i); }
    else if (!branch.timeSurvivesCleaning(i)) { ++nHsCleaned; hsUntimed.push_back(i); }
  }
  tracks.erase(std::remove_if(tracks.begin(), tracks.end(), [&](int i) {
                 return branch.trackTimeValid[i] == 1 && !branch.timeSurvivesCleaning(i); }), tracks.end());
  if (tracks.empty()) { if (dump) std::cout << "no tracks: " << entry << "\n"; return; }
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
  // qualifying forward jets, as Cluster::calculateTime's in-jet subset
  std::vector<std::pair<double, double>> fj; std::vector<int> fjIdx;
  for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
    if (branch.isJetRemoved(j) || branch.topoJetPt[j] < MIN_JET_PT) continue;
    const double ae = std::abs((double)branch.topoJetEta[j]);
    if (ae < MIN_ABS_ETA_JET || ae > MAX_ABS_ETA_JET) continue;
    fj.push_back({branch.topoJetEta[j], branch.topoJetPhi[j]}); fjIdx.push_back(j);
  }
  auto jetOf = [&](int t) {   // the first qualifying forward jet within dR < 0.4 (ntuple index), else -1
    for (size_t q = 0; q < fj.size(); ++q)
      if (std::hypot(branch.trackEta[t] - fj[q].first, TVector2::Phi_mpi_pi(branch.trackPhi[t] - fj[q].second)) < 0.4)
        return fjIdx[q];
    return -1;
  };
  auto inJet = [&](int t) { return jetOf(t) >= 0; };
  auto mistimed = [&](int t) {
    return branch.trackToTruthvtx[t] == 0 && branch.trackTimeValid[t] == 1 &&
           std::abs(branch.trackTime[t] - tHS) / branch.trackTimeRes[t] >= 3.0;
  };
  // the misassigned-times oracle: drop the mistimed HS tracks, re-cluster, re-choose
  int nMis = 0; double t0Fix = -99999.0;
  { std::vector<int> keep;
    for (int t : tracks) { if (mistimed(t)) ++nMis; else keep.push_back(t); }
    if (nMis > 0 && !keep.empty()) {
      auto fixCl = clusterTracksInTime(keep, &branch, DIST_CUT_CONE, false, true, -1.0,
                                       ClusteringMethod::ITERATIVE, false, false, false);
      int jW = -1;
      for (size_t k = 0; k < fixCl.size(); ++k)
        if (jW < 0 || fixCl[k].scores.at(Score::WAVES.id) > fixCl[jW].scores.at(Score::WAVES.id)) jW = k;
      if (jW >= 0) t0Fix = fixCl[jW].calculateTime(Score::WAVES, &branch);
    } }
  std::ofstream js(out);
  js << std::setprecision(8) << "{\n  \"file\": \"" << filePath << "\", \"entry\": " << entry
     << ", \"t_hs\": " << tHS << ", \"z_pv\": " << branch.recoVtxZ[0]
     << ", \"athena_t0\": " << branch.recoVtxTime[0] << ", \"athena_valid\": " << (int)branch.recoVtxValid[0]
     << ", \"t0_waves\": " << t0W << ", \"i_waves\": " << iW << ", \"i_sumpt\": " << iP << ", \"i_hs\": " << iH
     << ", \"n_mistimed_hs\": " << nMis << ", \"t0_waves_fix\": " << t0Fix
     << ", \"n_hs_assoc\": " << nHsAssoc << ", \"n_hs_notime\": " << nHsNoTime << ", \"n_hs_cleaned\": " << nHsCleaned
     << ", \"hs_untimed\": [";
  for (size_t i = 0; i < hsUntimed.size(); ++i) js << (i ? ", " : "") << hsUntimed[i];
  js << "]"
     << ",\n  \"clusters\": [\n";
  for (size_t k = 0; k < clusters.size(); ++k) {
    const auto& c = clusters[k];
    double sw = 0.0, swt = 0.0; int nhs = 0;
    for (int t : c.trackIndices)
      if (branch.trackToTruthvtx[t] == 0 && branch.trackTimeValid[t] == 1) {
        const double w = 1.0 / (branch.trackTimeRes[t] * branch.trackTimeRes[t]);
        sw += w; swt += w * branch.trackTime[t]; ++nhs;
      }
    js << (k ? ",\n" : "") << "    {\"t\": " << c.values[0] << ", \"sigma\": " << c.sigmas[0] << ", \"sumpt\": " << sumPt[k]
       << ", \"purity\": " << c.purity << ", \"waves\": " << c.scores.at(Score::WAVES.id)
       << ", \"nhs\": " << nhs << ", \"t_hscore\": " << (sw > 0 ? swt / sw : -99999.0) << ", \"tracks\": [";
    for (size_t i = 0; i < c.trackIndices.size(); ++i) {
      const int t = c.trackIndices[i];
      js << (i ? ", " : "") << "{\"idx\": " << t << ", \"t\": " << c.allTimes[i] << ", \"res\": " << branch.trackTimeRes[t]
         << ", \"pt\": " << branch.trackPt[t] << ", \"eta\": " << branch.trackEta[t] << ", \"phi\": " << branch.trackPhi[t]
         << ", \"hs\": " << (branch.trackToTruthvtx[t] == 0 ? 1 : 0)
         << ", \"injet\": " << (inJet(t) ? 1 : 0) << ", \"jet\": " << jetOf(t) << ", \"mistimed\": " << (mistimed(t) ? 1 : 0) << "}";
    }
    js << "]}";
  }
  js << "\n  ],\n  \"jets\": [";
  bool first = true;
  for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
    if (branch.isJetRemoved(j) || branch.topoJetPt[j] < MIN_JET_PT) continue;
    const double ae = std::abs((double)branch.topoJetEta[j]);
    if (ae < MIN_ABS_ETA_JET || ae > MAX_ABS_ETA_JET) continue;
    js << (first ? "" : ", ") << "{\"idx\": " << j << ", \"pt\": " << branch.topoJetPt[j] << ", \"eta\": " << branch.topoJetEta[j]
       << ", \"phi\": " << branch.topoJetPhi[j] << "}";
    first = false;
  }
  js << "]\n}\n";
  std::cout << "wrote " << out << "\n";
}

void clustering_dump(std::string filePath, Long64_t eventNum = -1, std::string out = "") {
  CLEAN_TIMES = true;                       // binds Track_hitInLastLayer
  TChain chain("ntuple");
  chain.Add(filePath.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  if (eventNum >= 0) { reader.SetEntry(eventNum); dumpCurrent(branch, filePath, eventNum, out); return; }
  const Long64_t n = reader.GetEntries(true);
  for (Long64_t e = 0; e < n; ++e) { reader.SetEntry(e); dumpCurrent(branch, filePath, e, ""); }
}

void clustering_dump_list(std::string listPath, std::string outDir) {
  CLEAN_TIMES = true;
  std::map<std::string, std::vector<std::pair<Long64_t, std::string>>> byFile;
  std::ifstream in(listPath);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line); std::string f, name; Long64_t e;
    if (ls >> f >> e >> name) byFile[f].push_back({e, name});
  }
  for (auto& [f, evs] : byFile) {
    TChain chain("ntuple");
    chain.Add(f.c_str());
    TTreeReader reader(&chain);
    BranchPointerWrapper branch(reader);
    for (auto& [e, name] : evs) { reader.SetEntry(e); dumpCurrent(branch, f, e, outDir + "/" + name + ".json"); }
  }
}
