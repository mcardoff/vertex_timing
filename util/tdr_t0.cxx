// -----------------------------------------------------------------------------
// tdr_t0 -- the vertex t0 part of the TDR-structured presentation (ATLAS-TDR-031
// Sec. 3.2.2, "Determination of the time of the primary vertex"), with the TDR's
// track-time cleaning (Sec. 3.2.1) against the ntuples' times as they are.
//
// One row per event: t0 - t_HS for
//   hgtd           Athena's RecoVtx_time (computed upstream: cannot be re-cleaned)
//   trkptz_raw / trkptz_clean, waves_raw / waves_clean
//                  the TRKPTZ and WAVeS t0s of util/tdr_common.h -- the ones the
//                  presentation's R_pT plots use -- from the ntuple's track times
//                  and from the cleaned ones, in the same event
// and dt_<name>_clean for every later t0 source of util/tdr_common.h (wfix_*, the waves_* selection levers)
// (NaN where a t0 is absent), with the number of forward hard-scatter tracks with a
// valid time (n_hs_trk: the main analysis's EventCounts on the 3-sigma list, from the
// ntuple's times, so raw and cleaned share one x axis), of forward jets, and the pT of the
// hardest forward truth-pileup jet (pu_fwd_ptmax, 0 = none).
// Events: the main clustering analysis's gate -- Z+jets lepton selection and overlap
// removal, passBasicCuts (MIN_JETS, vertex within MAX_VTX_DZ), passJetPtCut.
//
// Fits: FIT_PEAKTAIL (src/plotting_utilities.h) of each t0's inclusive residual, run
// by --fit=<tree file> (writes <stem>_fits.root) -- automatically after an unsharded
// event loop; for a sharded run, hadd -f the trees first, then --fit the merge.
// The core fraction (|t0 - t_HS| < 60 ps over ALL events, a missing t0 a failure) is
// computed from the tree by python/tdr_t0_plots.py.
//
// Usage:  tdr_t0 [--sample=<s>] [--ntuple-dir=<d>] [--file-shard=i/N] [--max-events=N]
//         tdr_t0 --fit=<tdr_t0 tree file>
// Needs Track_hitInLastLayer: the skims drop it, so read the raw ntuples (--ntuple-dir=).
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <boost/filesystem.hpp>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "sample_config.h"
#include "tdr_common.h"
#include "plotting_utilities.h"

using namespace MyUtl;

namespace {
const std::vector<std::string> NAMES = {"hgtd", "trkptz_raw", "trkptz_clean", "waves_raw", "waves_clean"};

int fitTree(const std::string& in) {
  TFile f(in.c_str(), "READ");
  TTree* t = f.IsZombie() ? nullptr : (TTree*)f.Get("events");
  if (!t) { std::cerr << "[tdr_t0] no tree 'events' in " << in << "\n"; return 1; }
  const std::string out = in.substr(0, in.rfind(".root")) + "_fits.root";
  TFile fo(out.c_str(), "RECREATE");
  const Long64_t n = t->GetEntries();
  std::printf("[tdr_t0] fits of %s (%lld events)\n  %-14s %16s %14s\n", in.c_str(), (long long)n, "t0", "sigma_core [ps]", "peak RMS [ps]");
  for (const auto& name : NAMES) {
    double dt = 0.0;
    t->SetBranchAddress(("dt_" + name).c_str(), &dt);
    TH1D h(("h_" + name).c_str(), ";t_{0} - t_{HS} [ps];Events", 200, -400, 400);
    for (Long64_t i = 0; i < n; ++i) { t->GetEntry(i); if (std::isfinite(dt)) h.Fill(dt); }
    TF1* fit = makePeakTailFit(&h, false);
    const double sp = fit->GetParameter("SigmaPeak"), fs = fit->GetParameter("FracSh"), r = fit->GetParameter("RatioSh");
    const double sc = sp / std::sqrt(1.0 - fs + fs * r * r), se = fit->GetParError(2) * sc / sp;
    std::printf("  %-14s %9.1f +- %3.1f %8.1f +- %3.1f\n", name.c_str(), sc, se, sp, fit->GetParError(2));
    fo.cd(); h.Write();
    TParameter<double>(("sigcore_" + name).c_str(), sc).Write();
    TParameter<double>(("sigcore_err_" + name).c_str(), se).Write();
    TParameter<double>(("sigpeak_" + name).c_str(), sp).Write();
    t->ResetBranchAddresses();
  }
  fo.Close();
  std::cout << "[tdr_t0] wrote " << out << "\n";
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
  TH1::AddDirectory(kFALSE);
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--fit=", 0) == 0) return fitTree(a.substr(6));
  }
  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL = cfg.energyLabel; OUTPUT_DIR = cfg.outputDir; SAMPLE_NAME = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  FILE_SHARD = resolveShard(argc, argv);
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);
  const Inflation INFL = inflationFor(SAMPLE_NAME);
  const std::string outPath = histFilePath("tdr_t0.root");
  {
    const auto dir = boost::filesystem::path(outPath).parent_path();
    if (!dir.empty()) boost::filesystem::create_directories(dir);
  }
  std::printf("[tdr_t0] sample=%s: t0s with and without the TDR track-time cleaning\n",
              SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME.c_str());
  printInflation(SAMPLE_NAME);
  std::cout << "[tdr_t0] out=" << outPath << "\n";

  EXTENDED_BRANCHES = true;   // Track_recoVtx_idx: the waves_pv0 / waves_npu t0s (tdr_common.h)
  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), FILE_SHARD);
  recordAvailableBranches(chain);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  const tdr::TimeCleaning RAW = tdr::setupCleaning(reader, tdr::CleanMode::OFF);
  const tdr::TimeCleaning CLEAN = tdr::setupCleaning(reader, tdr::CleanMode::T0);

  TFile out(outPath.c_str(), "RECREATE");
  if (out.IsZombie()) { std::cerr << "[tdr_t0] cannot open " << outPath << "\n"; return 1; }
  TTree tree("events", "t0 - t_HS per t0 version (NaN: no t0), main-analysis gate");
  const int N = (int)NAMES.size();
  std::vector<double> dt(N);
  int file_idx = -1, nHsTrk = 0, nFwdJet = 0;
  Long64_t entry = -1;
  float weight = 0.f;
  tree.Branch("file_idx", &file_idx, "file_idx/I");
  tree.Branch("entry", &entry, "entry/L");
  tree.Branch("weight", &weight, "weight/F");
  tree.Branch("n_hs_trk", &nHsTrk, "n_hs_trk/I");
  tree.Branch("n_fwd_jet", &nFwdJet, "n_fwd_jet/I");
  for (int k = 0; k < N; ++k) tree.Branch(("dt_" + NAMES[k]).c_str(), &dt[k], ("dt_" + NAMES[k] + "/D").c_str());
  // the failure breakdown (cleaned times): what the clustering offered (tdr::T0Extra)
  int nClus = 0; double dtBest = 0.0, dtWavesFull = 0.0;
  tree.Branch("n_clusters_clean", &nClus, "n_clusters_clean/I");
  tree.Branch("dt_best_clean", &dtBest, "dt_best_clean/D");
  tree.Branch("dt_waves_full_clean", &dtWavesFull, "dt_waves_full_clean/D");
  double dtHs = 0.0, hsPur = 0.0;
  tree.Branch("dt_hs_clean", &dtHs, "dt_hs_clean/D");
  tree.Branch("hs_purity_clean", &hsPur, "hs_purity_clean/D");
  int nHsT = 0, nHsMis = 0;
  double hsPtT = 0.0, hsPtMis = 0.0, hsPtMisHit = 0.0, hsPtJet = 0.0;
  tree.Branch("n_hs_timed_clean", &nHsT, "n_hs_timed_clean/I");
  tree.Branch("n_hs_mistimed_clean", &nHsMis, "n_hs_mistimed_clean/I");
  tree.Branch("hs_pt_timed_clean", &hsPtT, "hs_pt_timed_clean/D");
  tree.Branch("hs_pt_mistimed_clean", &hsPtMis, "hs_pt_mistimed_clean/D");
  tree.Branch("hs_pt_misassigned_clean", &hsPtMisHit, "hs_pt_misassigned_clean/D");
  tree.Branch("hs_pt_injet_clean", &hsPtJet, "hs_pt_injet_clean/D");
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  // every later t0 source of tdr_common.h (the truth-fixed WAVeS t0s, the selection levers of todo vt-5t0), from the
  // cleaned t0 input as the R_pT plots use them: dt_<name>_clean
  std::vector<double> dtX(tdr::N_T0, NaN);
  for (int k = tdr::T0_WFIX_SEL; k < tdr::N_T0; ++k) {
    const std::string n = std::string("dt_") + tdr::T0_NAMES[k] + "_clean";
    tree.Branch(n.c_str(), &dtX[k], (n + "/D").c_str());
  }
  // the event's hardest forward truth-pileup jet (2.4 < |eta| < 4.0, pT > 30 GeV, isJetPaperPU; 0 = none): the
  // > 50 GeV capture's population (util/scratch/kmeans_split_test.cxx's category)
  double puFwdPtMax = 0.0;
  tree.Branch("pu_fwd_ptmax", &puFwdPtMax, "pu_fwd_ptmax/D");

  Long64_t nRead = 0, nSel = 0;
  while (reader.Next()) {
    if (maxEvents > 0 && nRead >= maxEvents) break;
    ++nRead;
    if (nRead % 100000 == 0) std::cout << "[tdr_t0] read " << nRead << ", selected " << nSel << std::endl;
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection() || branch.vetoLeptonOverlap()) continue;
    if (!branch.passBasicCuts() || !branch.passJetPtCut()) continue;
    ++nSel;
    const EventCounts counts(&branch, getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0), true);
    nHsTrk = counts.nForwardTrackHS;
    nFwdJet = counts.nForwardJet;
    const tdr::T0Set raw = tdr::computeT0s(branch, INFL, RAW);
    tdr::T0Extra ex;
    const tdr::T0Set cln = tdr::computeT0s(branch, INFL, CLEAN, &ex);
    nClus = ex.nClusters;  dtBest = ex.bestDt;  dtWavesFull = ex.wavesFullDt;  dtHs = ex.hsDt;  hsPur = ex.hsPurity;
    nHsT = ex.nHsTimed;  nHsMis = ex.nHsMistimed;  hsPtT = ex.hsPtTimed;  hsPtMis = ex.hsPtMistimed;
    hsPtMisHit = ex.hsPtMisassigned;  hsPtJet = ex.hsPtInJet;
    const double tHS = branch.truthVtxTime[0];
    auto put = [&](int k, const tdr::T0Set& s, int idx) { dt[k] = s.ok[idx] ? s.t[idx] - tHS : NaN; };
    put(0, raw, tdr::T0_HGTD);
    put(1, raw, tdr::T0_TRKPTZ); put(2, cln, tdr::T0_TRKPTZ);
    put(3, raw, tdr::T0_WAVES);  put(4, cln, tdr::T0_WAVES);
    for (int k = tdr::T0_WFIX_SEL; k < tdr::N_T0; ++k) dtX[k] = cln.ok[k] ? cln.t[k] - tHS : NaN;
    puFwdPtMax = 0.0;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j) || branch.topoJetPt[j] < tdr::JET_PT_MIN) continue;
      const double ae = std::abs((double)branch.topoJetEta[j]);
      if (ae < tdr::FWD_ETA_MIN || ae >= tdr::FWD_ETA_MAX) continue;
      if (branch.isJetPaperPU(branch.topoJetEta[j], branch.topoJetPhi[j])) puFwdPtMax = std::max(puFwdPtMax, (double)branch.topoJetPt[j]);
    }
    const std::string fpath = chain.GetCurrentFile()->GetName();
    { const size_t p = fpath.rfind('_'), d = (p == std::string::npos) ? p : fpath.find('.', p);
      file_idx = (p != std::string::npos && d != std::string::npos) ? std::atoi(fpath.substr(p + 1, d - p - 1).c_str()) : -1; }
    entry = chain.GetTree()->GetReadEntry();
    weight = (float)*branch.weight;
    tree.Fill();
  }
  if (nRead == 0) { std::cerr << "[tdr_t0] ERROR: read 0 events\n"; return 1; }

  TTree meta("meta", "run record (one row per job; sum over shards)");
  std::string sample = SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME, ndir = cfg.ntupleDir;
  meta.Branch("n_read", &nRead, "n_read/L");
  meta.Branch("n_selected", &nSel, "n_selected/L");
  meta.Branch("sample", &sample);
  meta.Branch("ntuple_dir", &ndir);
  meta.Fill();
  out.Write();
  out.Close();
  std::printf("[tdr_t0] read %lld, selected %lld (main-analysis gate)\n", (long long)nRead, (long long)nSel);
  std::cout << "[tdr_t0] wrote " << outPath << "\n";
  return FILE_SHARD.given ? 0 : fitTree(outPath);   // a shard is fitted after the merge
}
