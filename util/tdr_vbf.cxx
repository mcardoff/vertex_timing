// -----------------------------------------------------------------------------
// tdr_vbf -- TDR Sec. 3.4.1, "Vector boson fusion Higgs production": VBF H->inv
// against Z+jets after the TDR's loose VBF preselection, one row per event.
//
// Selection (TDR p. 52):
//   vertex |z_reco - z_truth| < MAX_VTX_DZ; Z+jets: Z->ll with the leptons made
//   invisible (lepton-jet overlap removal), standing in for the TDR's Z->nunu;
//   >= 2 jets (EMTopo, pT > 30 GeV, |eta| < 4.5), the leading pT > 75 GeV and the
//   second pT > 50 GeV, |Delta eta(j1, j2)| > 3, after "all jets are required to
//   pass an ITk-only R_pT pileup jet tagger operating at 85% hard-scatter
//   efficiency" -- read as a jet cleaning: a failing jet is dropped and the
//   selection runs on the jets that remain.
// The tagger is NOT applied here. Its thresholds come from tdr_pujet's working
// points, and the timing taggers replace it offline (python/tdr_vbf_plots.py), so
// the tree keeps every event that could pass once jets are dropped: some pair of
// jets with pT > 75 / > 50 GeV and |Delta eta| > 3. presel_all flags the events
// that pass on all jets, before any tagging; mjj, deta, dphi and cat describe that
// all-jets leading pair, and the offline step recomputes them on the survivors.
// Categories (Fig 3.24): central-central / central-forward / forward-forward by
// the leading two jets' |eta| > 2.4 -- cat = 0 / 1 / 2.
// Recoil proxies, stored only (the TDR preselection has no E_T^miss cut):
//   met_truth  |sum pT(nu)| of the truth record (raw ntuples only; -1 otherwise)
//   z_pt       the dilepton pT (Z+jets; -1 otherwise)
//
// Usage:  tdr_vbf [--sample=<s>] [--ntuple-dir=<d>] [--file-shard=i/N]
//                 [--max-events=N] [--gate-sigma=X (default 2)] [--clean-times=all|t0|off (default all)] [--rpt-match=ghost|dr02 (default ghost)] [--out=<path>]
// Output: <prefix>tdr_vbf.root (TTrees "events", "meta"; merge shards with hadd -f).
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <boost/filesystem.hpp>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "sample_config.h"
#include "tdr_common.h"

using namespace MyUtl;

namespace {
constexpr double J1_PT_MIN = 75.0, J2_PT_MIN = 50.0, DETA_MIN = 3.0;   // TDR p. 52
constexpr double JET_ETA_MAX = 4.5;                                      // the jets entering the event
}

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
  if (!(gate > 0.0)) { std::cerr << "[tdr_vbf] --gate-sigma must be > 0\n"; return 1; }
  if (gate != tdr::GATE_SIGMA_TDR) {
    std::ostringstream os; os << "gate" << std::fixed << std::setprecision(1) << gate;
    std::string tag = os.str(); std::replace(tag.begin(), tag.end(), '.', 'p');
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  const tdr::CleanMode cleanMode = tdr::cleaningMode(argc, argv);   // before the output path: it tags the file
  const tdr::RptMatch rptMatch = tdr::rptMatchMode(argc, argv);     // likewise
  const Inflation INFL = inflationFor(SAMPLE_NAME);
  const std::string outPath = outOverride.empty() ? histFilePath("tdr_vbf.root") : outOverride;
  {
    const auto dir = boost::filesystem::path(outPath).parent_path();
    if (!dir.empty()) boost::filesystem::create_directories(dir);
  }
  std::printf("[tdr_vbf] sample=%s: jets > %.0f / %.0f GeV, |Deta(j1,j2)| > %.0f; gate %.1f sigma\n",
              SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME.c_str(), J1_PT_MIN, J2_PT_MIN, DETA_MIN, gate);
  std::cout << "[tdr_vbf] track-time cleaning: " << tdr::cleanModeName(cleanMode) << "; R_pT matching: " << tdr::rptMatchName(rptMatch) << "\n";
  std::cout << "[tdr_vbf] out=" << outPath << "\n";

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), FILE_SHARD);
  recordAvailableBranches(chain);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  const tdr::TimeCleaning CLEAN = tdr::setupCleaning(reader, cleanMode);
  const bool hasTruthMet = hasBranch("TruthPart_pdgId") && hasBranch("TruthPart_pt") && hasBranch("TruthPart_phi");
  std::unique_ptr<TTreeReaderArray<float>> tpPdg, tpPt, tpPhi;
  if (hasTruthMet) {
    tpPdg = std::make_unique<TTreeReaderArray<float>>(reader, "TruthPart_pdgId");
    tpPt  = std::make_unique<TTreeReaderArray<float>>(reader, "TruthPart_pt");
    tpPhi = std::make_unique<TTreeReaderArray<float>>(reader, "TruthPart_phi");
  }
  std::cout << "[tdr_vbf] truth MET: " << (hasTruthMet ? "from the truth record" : "UNAVAILABLE (met_truth = -1)") << "\n";

  TFile out(outPath.c_str(), "RECREATE");
  if (out.IsZombie()) { std::cerr << "[tdr_vbf] cannot open " << outPath << "\n"; return 1; }
  TTree tree("events", "TDR 3.4.1: every event that can pass the VBF preselection once R_pT-tagged jets are dropped");
  int file_idx = -1, cat = -1, nJets = 0;
  bool presel = false;
  Long64_t entry = -1;
  float weight = 0.f, metTruth = -1.f, zPt = -1.f;
  double mjj = 0.0, deta = 0.0, dphi = 0.0;
  double t0[tdr::N_T0] = {}, sig0[tdr::N_T0] = {}, infl[tdr::N_T0] = {};
  bool ok0[tdr::N_T0] = {};
  std::vector<float> jPt, jEta, jPhi;
  std::vector<int> jHs, jPu, jNtrk, jNtimed;
  std::vector<double> jItk, jSelf, jComb;
  std::vector<double> jT0[tdr::N_T0];
  tree.Branch("file_idx", &file_idx, "file_idx/I");
  tree.Branch("entry", &entry, "entry/L");
  tree.Branch("weight", &weight, "weight/F");
  tree.Branch("met_truth", &metTruth, "met_truth/F");
  tree.Branch("z_pt", &zPt, "z_pt/F");
  tree.Branch("n_jets", &nJets, "n_jets/I");
  tree.Branch("mjj", &mjj, "mjj/D");
  tree.Branch("deta", &deta, "deta/D");
  tree.Branch("dphi", &dphi, "dphi/D");
  tree.Branch("cat", &cat, "cat/I");
  tree.Branch("presel_all", &presel, "presel_all/O");
  for (int k = 0; k < tdr::N_T0; ++k) {
    const std::string s = tdr::T0_NAMES[k];
    tree.Branch(("t0_" + s).c_str(), &t0[k], ("t0_" + s + "/D").c_str());
    tree.Branch(("sig0_" + s).c_str(), &sig0[k], ("sig0_" + s + "/D").c_str());
    tree.Branch(("infl_" + s).c_str(), &infl[k], ("infl_" + s + "/D").c_str());
    tree.Branch(("ok_" + s).c_str(), &ok0[k], ("ok_" + s + "/O").c_str());
  }
  tree.Branch("jet_pt", &jPt);
  tree.Branch("jet_eta", &jEta);
  tree.Branch("jet_phi", &jPhi);
  tree.Branch("jet_hs", &jHs);
  tree.Branch("jet_pu", &jPu);
  tree.Branch("jet_ntrk", &jNtrk);
  tree.Branch("jet_ntimed", &jNtimed);
  tree.Branch("jet_rpt_itk", &jItk);
  for (int k = 0; k < tdr::N_T0; ++k)
    tree.Branch((std::string("jet_rpt_t0_") + tdr::T0_NAMES[k]).c_str(), &jT0[k]);
  tree.Branch("jet_rpt_self", &jSelf);
  tree.Branch("jet_rpt_comb", &jComb);

  Long64_t nRead = 0, nVtx = 0, nLep = 0, nSel = 0, nPresel = 0;
  double sumw = 0.0;
  while (reader.Next()) {
    if (maxEvents > 0 && nRead >= maxEvents) break;
    ++nRead;
    sumw += *branch.weight;
    if (nRead % 100000 == 0) std::cout << "[tdr_vbf] read " << nRead << ", selected " << nSel << std::endl;
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    if (std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]) > MAX_VTX_DZ) continue;
    ++nVtx;
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap()) continue;
    ++nLep;

    std::vector<int> js;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j)) continue;
      if (branch.topoJetPt[j] < tdr::JET_PT_MIN || std::abs((double)branch.topoJetEta[j]) >= JET_ETA_MAX) continue;
      js.push_back(j);
    }
    if (js.size() < 2) continue;
    std::sort(js.begin(), js.end(), [&](int x, int y) { return branch.topoJetPt[x] > branch.topoJetPt[y]; });
    // could some pair still pass once jets are dropped? (a superset of every tagger's selection)
    bool anyPair = false;
    for (size_t a = 0; a < js.size() && !anyPair; ++a) {
      if (branch.topoJetPt[js[a]] <= J1_PT_MIN) break;          // sorted: no leading candidate left
      for (size_t c = a + 1; c < js.size() && !anyPair; ++c) {
        if (branch.topoJetPt[js[c]] <= J2_PT_MIN) break;
        anyPair = std::abs((double)branch.topoJetEta[js[a]] - (double)branch.topoJetEta[js[c]]) > DETA_MIN;
      }
    }
    if (!anyPair) continue;
    ++nSel;
    const int j1 = js[0], j2 = js[1];
    deta = std::abs((double)branch.topoJetEta[j1] - (double)branch.topoJetEta[j2]);
    presel = branch.topoJetPt[j1] > J1_PT_MIN && branch.topoJetPt[j2] > J2_PT_MIN && deta > DETA_MIN;
    nPresel += presel;

    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(branch.topoJetPt[j1], branch.topoJetEta[j1], branch.topoJetPhi[j1], 0.0);
    v2.SetPtEtaPhiM(branch.topoJetPt[j2], branch.topoJetEta[j2], branch.topoJetPhi[j2], 0.0);
    mjj = (v1 + v2).M();
    dphi = std::abs(TVector2::Phi_mpi_pi((double)branch.topoJetPhi[j1] - (double)branch.topoJetPhi[j2]));
    const bool f1 = std::abs((double)branch.topoJetEta[j1]) > tdr::FWD_ETA_MIN;
    const bool f2 = std::abs((double)branch.topoJetEta[j2]) > tdr::FWD_ETA_MIN;
    cat = (int)f1 + (int)f2;
    nJets = (int)js.size();

    // recoil proxies (util/vbs_time_veto.cxx's definitions)
    metTruth = -1.f;  zPt = -1.f;
    if (hasTruthMet) {
      double px = 0.0, py = 0.0;
      for (size_t p = 0; p < tpPdg->GetSize(); ++p) {
        const long a = std::labs(std::lround((*tpPdg)[p]));
        if (a != 12 && a != 14 && a != 16) continue;
        px += (*tpPt)[p] * std::cos((*tpPhi)[p]);
        py += (*tpPt)[p] * std::sin((*tpPhi)[p]);
      }
      metTruth = (float)std::hypot(px, py);
    }
    if (branch.trackLeptonID) {
      std::vector<int> lep;
      for (int t = 0; t < (int)branch.trackLeptonID->GetSize(); ++t)
        if (branch.isGoodLepton(t)) lep.push_back(t);
      std::sort(lep.begin(), lep.end(), [&](int x, int y) { return branch.trackPt[x] > branch.trackPt[y]; });
      bool done = false;
      for (size_t x = 0; x < lep.size() && !done; ++x)
        for (size_t y = x + 1; y < lep.size() && !done; ++y) {
          const int px_ = branch.leptonPdg(lep[x]), py_ = branch.leptonPdg(lep[y]);
          if (std::abs(px_) != std::abs(py_) || px_ * py_ >= 0) continue;
          TLorentzVector lx, ly;
          lx.SetPtEtaPhiM(branch.trackPt[lep[x]], branch.trackEta[lep[x]], branch.trackPhi[lep[x]], 0.0);
          ly.SetPtEtaPhiM(branch.trackPt[lep[y]], branch.trackEta[lep[y]], branch.trackPhi[lep[y]], 0.0);
          zPt = (float)(lx + ly).Pt();
          done = true;
        }
    }

    const tdr::T0Set T = tdr::computeT0s(branch, INFL, CLEAN);
    const tdr::TrackLists L = tdr::rptTrackLists(branch);
    for (int k = 0; k < tdr::N_T0; ++k) { t0[k] = T.t[k]; sig0[k] = T.sig[k]; infl[k] = T.infl[k]; ok0[k] = T.ok[k]; }
    jPt.clear(); jEta.clear(); jPhi.clear(); jHs.clear(); jPu.clear(); jNtrk.clear(); jNtimed.clear();
    jItk.clear(); jSelf.clear(); jComb.clear();
    for (auto& v : jT0) v.clear();
    for (int j : js) {
      const float e = branch.topoJetEta[j], p = branch.topoJetPhi[j];
      jPt.push_back(branch.topoJetPt[j]);  jEta.push_back(e);  jPhi.push_back(p);
      jHs.push_back(branch.isJetPaperHS(e, p) ? 1 : 0);
      jPu.push_back(branch.isJetPaperPU(e, p) ? 1 : 0);
      // R_pT needs tracks: beyond |eta| = 4.0 it is undefined (stored as -1)
      if (std::abs((double)e) < tdr::FWD_ETA_MAX) {
        const tdr::JetRpt r = tdr::jetRpt(branch, j, L, T, gate, CLEAN, rptMatch);
        jNtrk.push_back(r.ntrk);  jNtimed.push_back(r.ntimed);
        jItk.push_back(r.itk);  jSelf.push_back(r.self);  jComb.push_back(r.comb);
        for (int k = 0; k < tdr::N_T0; ++k) jT0[k].push_back(r.t0[k]);
      } else {
        jNtrk.push_back(-1);  jNtimed.push_back(-1);
        jItk.push_back(-1.0);  jSelf.push_back(-1.0);  jComb.push_back(-1.0);
        for (int k = 0; k < tdr::N_T0; ++k) jT0[k].push_back(-1.0);
      }
    }
    const std::string fpath = chain.GetCurrentFile()->GetName();
    { const size_t pp = fpath.rfind('_'), d = (pp == std::string::npos) ? pp : fpath.find('.', pp);
      file_idx = (pp != std::string::npos && d != std::string::npos) ? std::atoi(fpath.substr(pp + 1, d - pp - 1).c_str()) : -1; }
    entry = chain.GetTree()->GetReadEntry();
    weight = (float)*branch.weight;
    tree.Fill();
  }
  if (nRead == 0) {
    std::cerr << "[tdr_vbf] ERROR: read 0 events -- check stderr for TTreeReader branch errors\n";
    return 1;
  }

  TTree meta("meta", "run record (one row per job; sum over shards)");
  std::string sample = SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME, ndir = cfg.ntupleDir;
  double j1 = J1_PT_MIN, j2 = J2_PT_MIN, dEtaMin = DETA_MIN;
  bool truthMet = hasTruthMet;
  meta.Branch("n_read", &nRead, "n_read/L");
  meta.Branch("sumw_read", &sumw, "sumw_read/D");
  meta.Branch("n_pass_vtx", &nVtx, "n_pass_vtx/L");
  meta.Branch("n_pass_lepton", &nLep, "n_pass_lepton/L");
  meta.Branch("n_selected", &nSel, "n_selected/L");
  meta.Branch("n_presel_all", &nPresel, "n_presel_all/L");
  meta.Branch("gate_sigma", &gate, "gate_sigma/D");
  std::string cleanName = tdr::cleanModeName(cleanMode);
  meta.Branch("time_cleaning", &cleanName);
  std::string matchName = tdr::rptMatchName(rptMatch);
  meta.Branch("rpt_match", &matchName);
  meta.Branch("j1_pt_min", &j1, "j1_pt_min/D");
  meta.Branch("j2_pt_min", &j2, "j2_pt_min/D");
  meta.Branch("deta_min", &dEtaMin, "deta_min/D");
  meta.Branch("has_truth_met", &truthMet, "has_truth_met/O");
  meta.Branch("sample", &sample);
  meta.Branch("ntuple_dir", &ndir);
  meta.Fill();
  out.Write();
  out.Close();
  std::printf("[tdr_vbf] read %lld, vertex %lld, lepton/OR %lld; written %lld (some pair passes), of which %lld pass on all jets\n",
              (long long)nRead, (long long)nVtx, (long long)nLep, (long long)nSel, (long long)nPresel);
  std::cout << "[tdr_vbf] wrote " << outPath << "\n";
  return 0;
}
