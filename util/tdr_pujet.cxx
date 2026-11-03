// -----------------------------------------------------------------------------
// tdr_pujet -- TDR Sec. 3.3.3, "Suppression of pileup jets", per jet.
//
// One row per reconstructed jet (EMTopo, pT > 30 GeV, |eta| < 4.0) in every event
// passing vertex quality (|z_reco - z_truth| < MAX_VTX_DZ) and, for Z+jets, the
// Z->ll selection with lepton-jet overlap removal -- no other event selection, as
// in the TDR's per-jet performance. Each row carries the jet's kinematics, its
// truth label (TDR: hard scatter within dR < 0.3 of a truth jet with pT > 10 GeV,
// pileup dR > 0.6 from every one with pT > 4 GeV), and R_pT for every tagger of
// util/tdr_common.h: ITk only, t0-gated (HGTD/Athena, TRKPTZ, WAVeS, truth),
// self-tagging and combined. Working points, ROCs (TDR Fig 3.20), the pileup-jet
// rate vs |eta| (Fig 3.21) and the 85% thresholds for tdr_vbf come from the tree
// offline, with python/tdr_pujet_plots.py.
//
// Separate from rpt_v5 on purpose (its forward window is 2.4-3.8 and its gate
// 3 sigma; this program follows the TDR); the ingredients are copied from it.
//
// Usage:  tdr_pujet [--sample=<s>] [--ntuple-dir=<d>] [--file-shard=i/N]
//                   [--max-events=N] [--gate-sigma=X (default 2)] [--clean-times=all|t0|off (default all)] [--out=<path>]
// Output: <prefix>tdr_pujet.root (TTrees "jets", "meta"; merge shards with hadd -f).
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <boost/filesystem.hpp>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "sample_config.h"
#include "tdr_common.h"

using namespace MyUtl;

int main(int argc, char** argv) {
  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL    = cfg.energyLabel;
  OUTPUT_DIR      = cfg.outputDir;
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  FILE_SHARD      = resolveShard(argc, argv);
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);

  double gate = tdr::GATE_SIGMA_TDR;
  std::string outOverride;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--gate-sigma=", 0) == 0) gate = std::stod(a.substr(13));
    if (a.rfind("--out=", 0) == 0) outOverride = a.substr(6);
  }
  if (!(gate > 0.0)) { std::cerr << "[tdr_pujet] --gate-sigma must be > 0\n"; return 1; }
  if (gate != tdr::GATE_SIGMA_TDR) {   // a different gate is a different tagger: tagged
    std::ostringstream os; os << "gate" << std::fixed << std::setprecision(1) << gate;
    std::string tag = os.str(); std::replace(tag.begin(), tag.end(), '.', 'p');
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  const tdr::CleanMode cleanMode = tdr::cleaningMode(argc, argv);   // before the output path: it tags the file
  const Inflation INFL = inflationFor(SAMPLE_NAME);
  const std::string outPath = outOverride.empty() ? histFilePath("tdr_pujet.root") : outOverride;
  {
    const auto dir = boost::filesystem::path(outPath).parent_path();
    if (!dir.empty()) boost::filesystem::create_directories(dir);
  }
  std::printf("[tdr_pujet] sample=%s gate=%.1f sigma, t0 inflation hgtd %.2f trkptz %.2f waves %.2f\n",
              SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME.c_str(), gate, INFL.hgtd, INFL.trkptz, INFL.waves);
  std::cout << "[tdr_pujet] track-time cleaning: " << tdr::cleanModeName(cleanMode) << "\n";
  std::cout << "[tdr_pujet] out=" << outPath << "\n";

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), FILE_SHARD);
  recordAvailableBranches(chain);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  const tdr::TimeCleaning CLEAN = tdr::setupCleaning(reader, cleanMode);

  TFile out(outPath.c_str(), "RECREATE");
  if (out.IsZombie()) { std::cerr << "[tdr_pujet] cannot open " << outPath << "\n"; return 1; }
  TTree jets("jets", "TDR 3.3.3: one row per jet, pT > 30 GeV, |eta| < 4.0");
  int file_idx = -1, hs = 0, pu = 0, ntrk = 0, ntimed = 0;
  Long64_t entry = -1;
  float weight = 0.f, pt = 0.f, eta = 0.f, phi = 0.f;
  bool okHgtd = false;
  double rItk = 0.0, rSelf = 0.0, rComb = 0.0;
  double rT0[tdr::N_T0] = {};
  jets.Branch("file_idx", &file_idx, "file_idx/I");
  jets.Branch("entry", &entry, "entry/L");
  jets.Branch("weight", &weight, "weight/F");
  jets.Branch("pt", &pt, "pt/F");
  jets.Branch("eta", &eta, "eta/F");
  jets.Branch("phi", &phi, "phi/F");
  jets.Branch("hs", &hs, "hs/I");
  jets.Branch("pu", &pu, "pu/I");
  jets.Branch("ntrk", &ntrk, "ntrk/I");
  jets.Branch("ntimed", &ntimed, "ntimed/I");
  jets.Branch("ok_hgtd", &okHgtd, "ok_hgtd/O");
  jets.Branch("rpt_itk", &rItk, "rpt_itk/D");
  for (int k = 0; k < tdr::N_T0; ++k) {
    const std::string n = std::string("rpt_t0_") + tdr::T0_NAMES[k];
    jets.Branch(n.c_str(), &rT0[k], (n + "/D").c_str());
  }
  jets.Branch("rpt_self", &rSelf, "rpt_self/D");
  jets.Branch("rpt_comb", &rComb, "rpt_comb/D");

  Long64_t nRead = 0, nVtx = 0, nLep = 0, nJetRows = 0;
  while (reader.Next()) {
    if (maxEvents > 0 && nRead >= maxEvents) break;
    ++nRead;
    if (nRead % 100000 == 0) std::cout << "[tdr_pujet] read " << nRead << ", jets " << nJetRows << std::endl;
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    if (std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]) > MAX_VTX_DZ) continue;
    ++nVtx;
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;   // no-op except Z+jets
    if (branch.vetoLeptonOverlap()) continue;
    ++nLep;

    // the jets first: most events have none in reach, and the t0s cost a clustering
    std::vector<int> js;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j)) continue;
      if (branch.topoJetPt[j] < tdr::JET_PT_MIN) continue;
      if (std::abs((double)branch.topoJetEta[j]) >= tdr::FWD_ETA_MAX) continue;
      js.push_back(j);
    }
    if (js.empty()) continue;
    const tdr::T0Set T = tdr::computeT0s(branch, INFL, CLEAN);
    const tdr::TrackLists L = tdr::rptTrackLists(branch);

    const std::string fpath = chain.GetCurrentFile()->GetName();
    { const size_t p = fpath.rfind('_'), d = (p == std::string::npos) ? p : fpath.find('.', p);
      file_idx = (p != std::string::npos && d != std::string::npos) ? std::atoi(fpath.substr(p + 1, d - p - 1).c_str()) : -1; }
    entry = chain.GetTree()->GetReadEntry();
    weight = (float)*branch.weight;
    okHgtd = T.ok[tdr::T0_HGTD];
    for (int j : js) {
      pt = branch.topoJetPt[j];  eta = branch.topoJetEta[j];  phi = branch.topoJetPhi[j];
      hs = branch.isJetPaperHS(eta, phi) ? 1 : 0;
      pu = branch.isJetPaperPU(eta, phi) ? 1 : 0;
      const tdr::JetRpt r = tdr::jetRpt(branch, j, L, T, gate, CLEAN);
      ntrk = r.ntrk;  ntimed = r.ntimed;
      rItk = r.itk;  rSelf = r.self;  rComb = r.comb;
      for (int k = 0; k < tdr::N_T0; ++k) rT0[k] = r.t0[k];
      jets.Fill();
      ++nJetRows;
    }
  }
  if (nRead == 0) {
    std::cerr << "[tdr_pujet] ERROR: read 0 events -- check stderr for TTreeReader branch errors\n";
    return 1;
  }

  TTree meta("meta", "run record (one row per job; sum over shards)");
  std::string sample = SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME, ndir = cfg.ntupleDir;
  double fwdMin = tdr::FWD_ETA_MIN, fwdMax = tdr::FWD_ETA_MAX, ptMin = tdr::JET_PT_MIN;
  double iH = INFL.hgtd, iT = INFL.trkptz, iW = INFL.waves;
  meta.Branch("n_read", &nRead, "n_read/L");
  meta.Branch("n_pass_vtx", &nVtx, "n_pass_vtx/L");
  meta.Branch("n_pass_lepton", &nLep, "n_pass_lepton/L");
  meta.Branch("n_jets", &nJetRows, "n_jets/L");
  meta.Branch("gate_sigma", &gate, "gate_sigma/D");
  std::string cleanName = tdr::cleanModeName(cleanMode);
  meta.Branch("time_cleaning", &cleanName);
  meta.Branch("jet_pt_min", &ptMin, "jet_pt_min/D");
  meta.Branch("fwd_eta_min", &fwdMin, "fwd_eta_min/D");
  meta.Branch("fwd_eta_max", &fwdMax, "fwd_eta_max/D");
  meta.Branch("infl_hgtd", &iH, "infl_hgtd/D");
  meta.Branch("infl_trkptz", &iT, "infl_trkptz/D");
  meta.Branch("infl_waves", &iW, "infl_waves/D");
  meta.Branch("sample", &sample);
  meta.Branch("ntuple_dir", &ndir);
  meta.Fill();
  out.Write();
  out.Close();
  std::printf("[tdr_pujet] read %lld, vertex %lld, lepton/OR %lld; %lld jets written\n",
              (long long)nRead, (long long)nVtx, (long long)nLep, (long long)nJetRows);
  std::cout << "[tdr_pujet] wrote " << outPath << "\n";
  return 0;
}
