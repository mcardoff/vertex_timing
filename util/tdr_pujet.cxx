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
//                   [--lojo]   (+ the WAVeS t0 with the tested jet left out, todo vt-9zl; tags the output "lojo")
//                   [--max-events=N] [--gate-sigma=X (default 2)] [--clean-times=all|t0|off (default all)] [--rpt-match=ghost|dr02 (default ghost)] [--rpt-trk-ptmax=X|none (default 30)] [--out=<path>]
// Output: <prefix>tdr_pujet.root (TTrees "jets", "meta"; merge shards with hadd -f).
// The jvt_rpt* columns need Track_recoVtx_idx: an input without it gets -1 there for every jet, with a warning.
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <boost/filesystem.hpp>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
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
  bool lojo = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--gate-sigma=", 0) == 0) gate = std::stod(a.substr(13));
    if (a.rfind("--out=", 0) == 0) outOverride = a.substr(6);
    if (a == "--lojo") lojo = true;
  }
  if (!(gate > 0.0)) { std::cerr << "[tdr_pujet] --gate-sigma must be > 0\n"; return 1; }
  if (gate != tdr::GATE_SIGMA_TDR) {   // a different gate is a different tagger: tagged
    std::ostringstream os; os << "gate" << std::fixed << std::setprecision(1) << gate;
    std::string tag = os.str(); std::replace(tag.begin(), tag.end(), '.', 'p');
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  if (lojo) SELECTION_TAG = SELECTION_TAG.empty() ? "lojo" : SELECTION_TAG + "_lojo";
  const tdr::CleanMode cleanMode = tdr::cleaningMode(argc, argv);   // before the output path: it tags the file
  const tdr::RptMatch rptMatch = tdr::rptMatchMode(argc, argv);     // likewise
  tdr::rptTrkPtMaxMode(argc, argv);                                  // likewise
  const Inflation INFL = inflationFor(SAMPLE_NAME);
  const std::string outPath = outOverride.empty() ? histFilePath("tdr_pujet.root") : outOverride;
  {
    const auto dir = boost::filesystem::path(outPath).parent_path();
    if (!dir.empty()) boost::filesystem::create_directories(dir);
  }
  std::printf("[tdr_pujet] sample=%s gate=%.1f sigma, t0 inflation hgtd %.2f trkptz %.2f waves %.2f\n",
              SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME.c_str(), gate, INFL.hgtd, INFL.trkptz, INFL.waves);
  std::cout << "[tdr_pujet] track-time cleaning: " << tdr::cleanModeName(cleanMode) << "; R_pT matching: " << tdr::rptMatchName(rptMatch) << "\n";
  std::cout << "[tdr_pujet] out=" << outPath << "\n";

  EXTENDED_BRANCHES = true;                      // JVT's R_pT needs the vertex fit's Track_recoVtx_idx
  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), FILE_SHARD);
  recordAvailableBranches(chain);
  // Without the vertex fit's assignment recoVtxOf is a sentinel for every track and JVT's R_pT would read 0 for every
  // jet, the value of a jet with no PV0-fitted track. The TDR R_pT columns do not need it, so the run goes on and the
  // JVT columns are written as -1 (an R_pT is never negative) instead.
  const bool haveVtxFit = hasBranch("Track_recoVtx_idx");
  if (!haveVtxFit)
    std::cerr << "[tdr_pujet] WARNING: this sample lacks Track_recoVtx_idx: jvt_rpt and jvt_rpt_t0_*"
              << (lojo ? " (and jvt_rpt_t0_waves_lojo)" : "") << " are written as -1 for every jet\n";
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
  double rJvt = 0.0, rJvtT0[tdr::N_T0] = {};   // JVT's R_pT (vertex-fit PV0 tracks), plain and t0-gated
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
  // the jet's R_pT tracks that are hard-scatter with a misassigned time (>= TRUTH_PULL_CUT sigma from t_HS, the
  // wfix_time fix's definition): the fix drops them from the t0 input, but the gate still sees them on the jet
  int nMisHs = 0;
  jets.Branch("nmistimed_hs", &nMisHs, "nmistimed_hs/I");
  jets.Branch("ok_hgtd", &okHgtd, "ok_hgtd/O");
  jets.Branch("rpt_itk", &rItk, "rpt_itk/D");
  for (int k = 0; k < tdr::N_T0; ++k) {
    const std::string n = std::string("rpt_t0_") + tdr::T0_NAMES[k];
    jets.Branch(n.c_str(), &rT0[k], (n + "/D").c_str());
  }
  jets.Branch("rpt_self", &rSelf, "rpt_self/D");
  jets.Branch("rpt_comb", &rComb, "rpt_comb/D");
  jets.Branch("jvt_rpt", &rJvt, "jvt_rpt/D");
  for (int k = 0; k < tdr::N_T0; ++k) {
    const std::string n = std::string("jvt_rpt_t0_") + tdr::T0_NAMES[k];
    jets.Branch(n.c_str(), &rJvtT0[k], (n + "/D").c_str());
  }

  // --lojo: the WAVeS t0 with the tested jet left out (util/tdr_common.h's wavesT0LeaveOut), forward jets only --
  // a central jet is outside WAVeS's jet list and its tracks outside the t0 input, so for it these copy the WAVeS
  // columns. dt = t0 - t_HS (NaN without a t0).
  double rLojo = 0.0, rLojet = 0.0, rJvtLojo = 0.0;
  double dtW = 0.0, dtLojo = 0.0, dtLojet = 0.0;
  int nInLojo = 0, nInLojet = 0;   // the leave-out cluster's tracks inside another forward jet (-1: central jet)
  if (lojo) {
    jets.Branch("rpt_t0_waves_lojo", &rLojo, "rpt_t0_waves_lojo/D");     // jet and its tracks left out
    jets.Branch("rpt_t0_waves_lojet", &rLojet, "rpt_t0_waves_lojet/D");  // jet left out of WAVeS's jet list only
    jets.Branch("jvt_rpt_t0_waves_lojo", &rJvtLojo, "jvt_rpt_t0_waves_lojo/D");
    jets.Branch("dt_waves", &dtW, "dt_waves/D");
    jets.Branch("dt_waves_lojo", &dtLojo, "dt_waves_lojo/D");
    jets.Branch("dt_waves_lojet", &dtLojet, "dt_waves_lojet/D");
    jets.Branch("ninjet_lojo", &nInLojo, "ninjet_lojo/I");
    jets.Branch("ninjet_lojet", &nInLojet, "ninjet_lojet/I");
  }

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
      const tdr::JetRpt r = tdr::jetRpt(branch, j, L, T, gate, CLEAN, rptMatch);
      ntrk = r.ntrk;  ntimed = r.ntimed;
      nMisHs = 0;
      for (int idx : tdr::jetRptTracks(branch, j, std::abs((double)eta) >= tdr::FWD_ETA_MIN ? L.fwd : L.cen, rptMatch))
        if (branch.trackToTruthvtx[idx] == 0 && tdr::timed(branch, idx, CLEAN, CLEAN.jets()) &&
            std::abs(branch.trackTime[idx] - branch.truthVtxTime[0]) / branch.trackTimeRes[idx] >= TRUTH_PULL_CUT) ++nMisHs;
      rItk = r.itk;  rSelf = r.self;  rComb = r.comb;
      for (int k = 0; k < tdr::N_T0; ++k) rT0[k] = r.t0[k];
      if (haveVtxFit) {
        const tdr::JvtRpt jr = tdr::jvtRpt(branch, j, T, gate, CLEAN);
        rJvt = jr.plain;
        for (int k = 0; k < tdr::N_T0; ++k) rJvtT0[k] = jr.t0[k];
      } else {
        rJvt = -1.0;
        for (int k = 0; k < tdr::N_T0; ++k) rJvtT0[k] = -1.0;
      }
      if (lojo) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double tHS = branch.truthVtxTime[0];
        dtW = T.ok[tdr::T0_WAVES] ? T.t[tdr::T0_WAVES] - tHS : nan;
        rLojo = rLojet = rT0[tdr::T0_WAVES];  rJvtLojo = rJvtT0[tdr::T0_WAVES];  dtLojo = dtLojet = dtW;
        nInLojo = nInLojet = -1;
        const bool fwd = std::abs((double)eta) >= tdr::FWD_ETA_MIN;
        if (fwd) {
          const std::vector<int> trk = tdr::jetRptTracks(branch, j, L.fwd, rptMatch);
          std::vector<int> pv;
          for (int idx : branch.topoJetGhostTrackIdx[j]) if (branch.recoVtxOf(idx) == 0) pv.push_back(idx);
          const double iw = T.infl[tdr::T0_WAVES];
          const tdr::WavesT0 a = tdr::wavesT0LeaveOut(branch, j, tdr::LeaveOut::JET_AND_TRACKS, CLEAN);
          rLojo = tdr::sumPt(branch, tdr::timeGate(branch, trk, a.t, a.sig, a.ok, iw, gate, CLEAN)) / pt;
          rJvtLojo = haveVtxFit ? tdr::sumPt(branch, tdr::timeGate(branch, pv, a.t, a.sig, a.ok, iw, gate, CLEAN)) / pt
                                : -1.0;
          dtLojo = a.ok ? a.t - tHS : nan;
          nInLojo = a.nInJet;
          const tdr::WavesT0 c = tdr::wavesT0LeaveOut(branch, j, tdr::LeaveOut::JET, CLEAN);
          rLojet = tdr::sumPt(branch, tdr::timeGate(branch, trk, c.t, c.sig, c.ok, iw, gate, CLEAN)) / pt;
          dtLojet = c.ok ? c.t - tHS : nan;
          nInLojet = c.nInJet;
        }
      }
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
  std::string matchName = tdr::rptMatchName(rptMatch);
  meta.Branch("rpt_match", &matchName);
  double rptPtMax = tdr::RPT_TRK_PT_MAX;
  meta.Branch("rpt_trk_ptmax", &rptPtMax, "rpt_trk_ptmax/D");
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
