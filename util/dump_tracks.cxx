// dump_tracks: loose per-track dump for OFFLINE track-list studies.
//
// export_training_data stores only the tracks that pass the nominal clustering
// selection (quality, 1 < pT < 30 GeV, z0 significance < MAX_NSIGMA), so nothing
// about a looser or different selection can be learned from its output. This
// writes every HGTD-timed track in acceptance within a LOOSE z window (8 sigma;
// 4.5 below 1 GeV; pT > 0.5 GeV, no ceiling, no quality cut), with the
// quantities the selection cuts on stored as columns, so python/kde_study can
// apply any rule and rerun the whole clustering + scoring chain.
//
// Event selection is export_training_data's (reco only). `n_nominal` is the
// number of dumped tracks that pass the nominal list, so the nominal event
// population (n_nominal > 0) can be reproduced exactly.
//
//   trees:  events (file_idx, event_num, vz, ttruth, n_nominal, weight)
//           tracks (file_idx, event_num, pt, eta, phi, z0, sigma_z0, sigma_d0,
//                   time, timeRes, nhgtd_hits, quality, nsig, dr_nearest_fwdjet,
//                   truth_is_hs)
// Sharded like the exporter; merge with hadd (TTrees only).

#include <TFile.h>
#include <TTree.h>
#include "clustering_constants.h"
#include "event_processing.h"
#include "sample_config.h"

using namespace MyUtl;

namespace {
constexpr double LOOSE_NSIG    = 8.0;   // z0 significance window of the dump
constexpr double LOOSE_NSIG_LO = 4.5;   // ... for tracks below 1 GeV (there are ~6x more of them)
constexpr double LOOSE_PT_MIN  = 0.5;   // GeV
}

auto main(int argc, char** argv) -> int {
  gErrorIgnoreLevel = kWarning;
  const auto cfg = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);
  MyUtl::SAMPLE_NAME  = cfg.sampleName;
  MyUtl::OUTPUT_DIR   = cfg.outputDir;
  MyUtl::ENERGY_LABEL = cfg.energyLabel;
  MyUtl::OVERLAP_REMOVAL = cfg.overlapRemoval;
  const Long64_t maxEvents = MyUtl::resolveMaxEvents(argc, argv);
  MyUtl::FILE_SHARD = MyUtl::resolveShard(argc, argv);
  const MyUtl::Shard shard = MyUtl::FILE_SHARD;

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), shard);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  if (MyUtl::SAMPLE_NAME.empty())
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/hists");
  const std::string outPath = MyUtl::histFilePath("trackdump.root");
  std::unique_ptr<TFile> out(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) { std::cerr << "ERROR: cannot open " << outPath << "\n"; return 1; }

  auto* evTree = new TTree("events", "one row per selected event");
  auto* trTree = new TTree("tracks", "one row per loosely selected timed track");
  float e_file, e_evt, e_vz, e_tt, e_nnom, e_w;
  evTree->Branch("file_idx", &e_file); evTree->Branch("event_num", &e_evt);
  evTree->Branch("vz", &e_vz); evTree->Branch("ttruth", &e_tt);
  evTree->Branch("n_nominal", &e_nnom); evTree->Branch("weight", &e_w);
  float t_pt, t_eta, t_phi, t_z0, t_sz0, t_sd0, t_time, t_res, t_nh, t_q, t_nsig, t_dr, t_hs;
  trTree->Branch("file_idx", &e_file); trTree->Branch("event_num", &e_evt);
  trTree->Branch("pt", &t_pt); trTree->Branch("eta", &t_eta); trTree->Branch("phi", &t_phi);
  trTree->Branch("z0", &t_z0); trTree->Branch("sigma_z0", &t_sz0); trTree->Branch("sigma_d0", &t_sd0);
  trTree->Branch("time", &t_time); trTree->Branch("timeRes", &t_res);
  trTree->Branch("nhgtd_hits", &t_nh); trTree->Branch("quality", &t_q); trTree->Branch("nsig", &t_nsig);
  trTree->Branch("dr_nearest_fwdjet", &t_dr); trTree->Branch("truth_is_hs", &t_hs);

  Long64_t nSeen = 0, nEv = 0, nTrk = 0;
  while (reader.Next()) {
    const Long64_t evtNum = chain.GetReadEntry() - chain.GetChainOffset();
    const int fileIdx = shard.active() ? shard.index + chain.GetTreeNumber() * shard.count
                                       : chain.GetTreeNumber();
    if (maxEvents > 0 && nSeen >= maxEvents) break;
    ++nSeen;

    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    if (!branch.passBasicCuts())       continue;
    if (!branch.passJetPtCut())        continue;

    // Qualifying forward jets, as Cluster::calculateTime's in-jet test.
    std::vector<std::pair<double,double>> jets;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j)) continue;
      if (branch.topoJetPt[j] < MIN_JET_PT) continue;
      const double je = branch.topoJetEta[j];
      if (std::abs(je) < MIN_ABS_ETA_JET || std::abs(je) > MAX_ABS_ETA_JET) continue;
      jets.push_back({je, branch.topoJetPhi[j]});
    }

    const double vz = branch.recoVtxZ[0];
    e_file = (float)fileIdx; e_evt = (float)evtNum; e_vz = (float)vz;
    e_tt = (float)branch.truthVtxTime[0]; e_w = (float)*branch.weight;
    int nNom = 0;
    std::vector<int> keep;
    for (size_t trk = 0; trk < branch.trackZ0.GetSize(); ++trk) {
      const double eta = std::abs(branch.trackEta[trk]);
      if (eta < MIN_HGTD_ETA || eta > MAX_HGTD_ETA) continue;
      if (branch.trackTimeValid[trk] != 1) continue;
      const double vzz = branch.trackVarZ0[trk];
      if (!(vzz > 0.0)) continue;
      const double dz = std::abs(branch.trackZ0[trk] - vz);
      const double nsig = dz / std::sqrt(Z0_VAR_INFLATION * vzz);
      const double pt = branch.trackPt[trk];
      if (pt < LOOSE_PT_MIN) continue;
      if (nsig > (pt < 1.0 ? LOOSE_NSIG_LO : LOOSE_NSIG)) continue;
      keep.push_back((int)trk);
      if (passTrackKinematics(trk, &branch, MIN_TRACK_PT, MAX_TRACK_PT) && nsig < MAX_NSIGMA) ++nNom;
    }
    e_nnom = (float)nNom;
    evTree->Fill(); ++nEv;
    for (int trk : keep) {
      t_pt = branch.trackPt[trk]; t_eta = branch.trackEta[trk]; t_phi = branch.trackPhi[trk];
      t_z0 = branch.trackZ0[trk]; t_sz0 = std::sqrt(branch.trackVarZ0[trk]);
      t_sd0 = std::sqrt(branch.trackVarD0[trk]);
      t_time = branch.trackTime[trk]; t_res = branch.trackTimeRes[trk];
      t_nh = (float)branch.trackHgtdHits[trk]; t_q = branch.trackQuality[trk] ? 1.f : 0.f;
      t_nsig = (float)(std::abs(branch.trackZ0[trk] - vz) / std::sqrt(Z0_VAR_INFLATION * branch.trackVarZ0[trk]));
      double best = 99.0;
      for (auto& [je, jp] : jets)
        best = std::min(best, std::hypot(je - (double)t_eta, TVector2::Phi_mpi_pi(jp - (double)t_phi)));
      t_dr = (float)best;
      t_hs = (branch.trackToTruthvtx[trk] == 0) ? 1.f : 0.f;
      trTree->Fill(); ++nTrk;
    }
  }
  out->cd(); evTree->Write(); trTree->Write(); out->Close();
  std::cout << "dump_tracks: " << nSeen << " read, " << nEv << " events, " << nTrk << " tracks -> " << outPath << "\n";
  return 0;
}
