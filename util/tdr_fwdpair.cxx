// -----------------------------------------------------------------------------
// tdr_fwdpair -- ATL-HGTD-PUB-2022-001 Sec. 6, "Suppression of events with forward
// pile-up jets": events whose two leading forward jets are both inside the HGTD
// acceptance, tagged as a PAIR. One row per selected event; the taggers, the cuts
// and the BDTs are applied offline (python/tdr_fwdpair_plots.py).
//
// Selection (the note's Sec. 3):
//   reco HS vertex within MAX_VTX_DZ (2 mm) of the truth HS vertex; Z+jets: Z->ll with
//   the leptons made invisible (lepton-jet overlap removal), as tdr_vbf does;
//   jets EMTopo, pT > 30 GeV, 2.4 < |eta| < 3.8 (inside the tracker: the note's
//   "fully contained" window); the two leading of THOSE jets must have opposite-sign
//   eta, the leading pT > 50 GeV and the second > 30 GeV. ("At least two jets in
//   2.4 < |eta| < 3.8 with opposite-sign eta; the leading (sub-leading) jet pT > 50
//   (30) GeV" -- read as the two leading jets of that window.)
// --pair=one (tag onefwd_; the user's extension, 2026-11-04): AT LEAST ONE of the two jets in HGTD.
//   The pair is the two leading jets with |eta| < 3.8 (still inside the tracker, so every leg has tracks),
//   > 50 / 30 GeV, opposite eta, and at least one of them in 2.4 < |eta| < 3.8: central-forward pairs join
//   the forward-forward ones (cat = 1 / 2). A central jet's JVF uses the central association (z0
//   significance < 5, tdr_common's CENTRAL_Z_SIGNIF), as its R_pT does: the forward parameterisation is
//   several times too tight at |eta| < 2.4. Central jets carry no HGTD time.
// --jvt=loose|tight (tag jvtLoose_ / jvtTight_; the user's analysis baseline, 2026-11-04: section 3 of
//   slides_tdr_timing gets "central JVT + fJVT loose"): the in-code JVT / fJVT taggers (src/jet_tagging.h,
//   as vbs_time_veto / rpt_v5_hist --jvt) remove failing jets BEFORE the pair is formed. Needs the vertex
//   fit's track assignment (Track_recoVtx_idx, an EXTENDED branch). Each selected jet's JVT R_pT and fJVT
//   are stored (j*_jvt_rpt, j*_fjvt).
// Labels (offline): signal = VBF events with both jets hard-scatter; background =
// events with at least one pileup jet (the TDR / PUB-note truth labels, as tdr_vbf).
//
// Per jet (j1, j2):
//   JVF  (the note's eq. 2) sum pT(PV0 tracks) / [sum pT(PV0 tracks) + sum pT(PU-vertex
//        tracks)] over the jet's ghost-associated tracks (1 < pT < 30 GeV, quality).
//        PV0 association is tdr_common's forward one (|dz| / getNewDzpara < 1.4, the
//        note's 2.5 sigma); a PU-vertex track is one NOT associated to PV0 that lies in
//        the same window around some other reco vertex (counted once). -1 when the jet
//        has no associated track (the note: such events are lost to a cut on both jets).
//        With a t0 (per t0 source): PV0 association also needs the 2-sigma t0 gate
//        (untimed tracks pass), and a track failing it may then join a PU vertex.
//   R_pT the TDR taggers' values (tdr_common's jetRpt), for comparison.
//   core jet time  the timing veto's estimator (util/vbs_time_veto.cxx jetTimeCore, ported): timed ghost
//        tracks with the standard kinematics (passTrackKinematics: 1 < pT < 30 GeV, quality, HGTD eta)
//        within dR < 0.2 of the axis, clustered in time at 2 sigma; the highest-sum-pT cluster's
//        inverse-variance mean and its sigma (uncalibrated: the veto scales it offline), and its track
//        count. Track times cleaned as everywhere here (vbs_time_veto does not clean). -999 / -1 / 0
//        when the jet has no such track (always for a central jet).
//   jet time  over the jet's ghost-associated tracks with a (cleaned) HGTD time, any
//        vertex: t_jet = the uncertainty-weighted mean, w = 1 / sigma_t^2 (the note: 35 ps /
//        sqrt(N_hits); here the ntuple's Track_timeRes), var_mean = 1 / sum w its variance -- the
//        clustering's cluster time, and the reading of the note's sigma_jet^2 that is used
//        (the user's, 2026-11-04); var_spread = the weighted variance of the times, the other
//        reading; and the leading-pT timed track's time and resolution. -999 / -1 when the jet
//        has no timed track. At mu = 200 about half the ghost-associated tracks are pileup, so
//        either reading separates weakly (local VBF, AUC alone: 0.599 mean, 0.585 spread;
//        leading track 0.810). Also stored, beyond the note: the same over the
//        jet's PV0-associated (z only, no t0) timed tracks -- t_hs, var_hs_spread, var_hs_mean,
//        ntimed_hs.
//
// Usage:  tdr_fwdpair [--sample=<s>] [--ntuple-dir=<d>] [--file-shard=i/N] [--max-events=N]
//                     [--gate-sigma=X (default 2)] [--clean-times=all|t0|off (default all)] [--pair=ff|one (default ff)]
//                     [--jvt=none|loose|tight|jvt4pt60|jvt4pt120|jvt4pt60t0|jvt4pt60t0tru|jvt4pt60f03|jvt4pt60f05|jvt4pt60f10 (default none)]
//   jvt4pt60t0(tru): JVT's cut on the t0-gated R_pT (WAVeS / truth t0, --gate-sigma, cleaned times; untimed tracks
//   count); jvt4pt60fNN: forward threshold 0.NN. j*_jvt_rpt is always the plain R_pT, j*_jvt_rpt_t0 the WAVeS-gated.
//                     [--out=<path>]
// Output: <prefix>tdr_fwdpair.root (TTrees "events", "meta"; merge shards with hadd -f).
// Track-time cleaning needs Track_hitInLastLayer: read the RAW ntuples (--ntuple-dir=).
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <boost/filesystem.hpp>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

#include "jet_tagging.h"
#include "sample_config.h"
#include "tdr_common.h"

using namespace MyUtl;

namespace {
constexpr double ETA_MIN = 2.4, ETA_MAX = 3.8;        // the note's jet window
constexpr double J1_PT_MIN = 50.0, J2_PT_MIN = 30.0;  // leading / sub-leading

struct JetOut {
  float pt = 0, eta = 0, phi = 0;
  int hs = 0, pu = 0, ntrk = 0, ntimed = 0;
  double jvf_itk = -1, rpt_itk = -1;
  std::array<double, tdr::N_T0> jvf_t0{}, rpt_t0{};
  double t_jet = -999, var_spread = -1, var_mean = -1, t_lead = -999, sig_lead = -1, pt_lead = -1;
  int ntimed_hs = 0;
  double t_hs = -999, var_hs_spread = -1, var_hs_mean = -1;
  int core_n = 0;
  double core_t = -999, core_sig = -1;
  // the same core time without the misassigned hard-scatter tracks (tdr_common's misassignedHs): what the
  // perfect-timing t0 sources' vetoes compare against
  int core_n_fix = 0;
  double core_t_fix = -999, core_sig_fix = -1;
  float jvt_rpt = -1, fjvt = -999;           // the JVT / fJVT discriminants (--jvt), -1 / -999 without
  float jvt_rpt_t0 = -1;                      // JVT's R_pT with the WAVeS t0 gate (--jvt), -1 without
  // the core-time variants (2026-11-14, todo vt-r2k.21 / vt-r2k.20): the same estimator on other track sets, to see
  // what a better jet time buys the veto (CV_NAMES below; -999 / -1 / 0 without a track, as core_*)
  static constexpr int N_CV = 5;
  std::array<int, N_CV> cv_n{};
  std::array<double, N_CV> cv_t{-999, -999, -999, -999, -999}, cv_sig{-1, -1, -1, -1, -1};
};
constexpr double CORE_DR = 0.2, CORE_DIST_CUT = 2.0;   // util/vbs_time_veto.cxx's JET_TIME_DR, JET_TIME_DIST_CUT
// the core-time variants, j*_core_{t,sig,n}_<name>:
//   r04       tracks within dR < 0.4 of the axis instead of 0.2
//   nocap     no 30 GeV upper track pT (MAX_TRACK_PT, copied from the t0 clustering, never chosen for jets)
//   r04nocap  both
//   zwin      the dR < 0.2 core restricted to tracks inside PV0's forward z window (getNewDzpara x 1.4, the R_pT
//             association): a hard-scatter core built from pileup tracks far in z (vt-r2k.20). It also strips the
//             pileup jets of their own tracks, so the veto goes blind (results/jet_core_variants.md)
//   zpref     the standard tracks, but the core is the time cluster with the most pT INSIDE that z window; the
//             highest-sum-pT cluster when no cluster has a track there (a pileup jet keeps its collision's time)
constexpr const char* CV_NAMES[JetOut::N_CV] = {"r04", "nocap", "r04nocap", "zwin", "zpref"};

// weighted mean, weighted variance of the times, and 1 / sum w (w = 1 / sigma_t^2)
void timeStats(const std::vector<std::pair<double, double>>& tw, double& t, double& vSpread, double& vMean) {
  double sw = 0.0, swt = 0.0;
  for (const auto& [ti, w] : tw) { sw += w; swt += w * ti; }
  if (!(sw > 0.0)) return;
  t = swt / sw;
  double v = 0.0;
  for (const auto& [ti, w] : tw) v += w * (ti - t) * (ti - t);
  vSpread = v / sw;
  vMean = 1.0 / sw;
}

// JVF numerator / denominator for one PV0-associated set
double jvfOf(BranchPointerWrapper& b, const std::vector<int>& trk, const std::vector<char>& inHs,
             const std::vector<char>& nearPu) {
  double hs = 0.0, pu = 0.0;
  for (size_t k = 0; k < trk.size(); ++k) {
    if (inHs[k]) hs += b.trackPt[trk[k]];
    else if (nearPu[k]) pu += b.trackPt[trk[k]];
  }
  return hs + pu > 0.0 ? hs / (hs + pu) : -1.0;
}

// PV0 / PU-vertex association: the forward parameterisation for forward jets, z0 significance centrally
bool assocZ(BranchPointerWrapper& b, int i, double z, bool fwd) {
  const double dz = std::abs(b.trackZ0[i] - z);
  if (fwd) return dz <= MyUtl::getNewDzpara(b.trackEta[i], b.trackPt[i]) * tdr::DZ0_PARA_SCALE;
  const double v = b.trackVarZ0[i];
  return v > 0.0 && dz / std::sqrt(v) <= tdr::CENTRAL_Z_SIGNIF;
}

JetOut jetQuantities(BranchPointerWrapper& b, int j, const tdr::TrackLists& L, const tdr::T0Set& T,
                     double gate, const tdr::TimeCleaning& C) {
  JetOut o;
  o.pt = b.topoJetPt[j];  o.eta = b.topoJetEta[j];  o.phi = b.topoJetPhi[j];
  o.hs = b.isJetPaperHS(o.eta, o.phi) ? 1 : 0;
  o.pu = b.isJetPaperPU(o.eta, o.phi) ? 1 : 0;
  // the jet's ghost-associated tracks with the R_pT kinematics
  std::vector<int> trk;
  const int nT = (int)b.trackPt.GetSize();
  for (int idx : b.topoJetGhostTrackIdx[j]) {
    if (idx < 0 || idx >= nT) continue;
    const double pt = b.trackPt[idx];
    if (pt < MIN_TRACK_PT || pt > MAX_TRACK_PT || !b.trackQuality[idx]) continue;
    trk.push_back(idx);
  }
  o.ntrk = (int)trk.size();
  const bool fwd = std::abs((double)o.eta) >= ETA_MIN;
  // association: PV0 (z) and "within some PU vertex's window"
  const int nVtx = (int)b.recoVtxZ.GetSize();
  std::vector<char> hsZ(trk.size(), 0), nearPu(trk.size(), 0);
  for (size_t k = 0; k < trk.size(); ++k) {
    const int i = trk[k];
    hsZ[k] = assocZ(b, i, b.recoVtxZ[0], fwd);
    for (int v = 1; v < nVtx && !nearPu[k]; ++v) nearPu[k] = assocZ(b, i, b.recoVtxZ[v], fwd);
  }
  o.jvf_itk = jvfOf(b, trk, hsZ, nearPu);
  for (int s = 0; s < tdr::N_T0; ++s) {
    std::vector<char> hsT(trk.size(), 0);
    for (size_t k = 0; k < trk.size(); ++k) {
      if (!hsZ[k]) continue;
      hsT[k] = !tdr::timeGate(b, {trk[k]}, T.t[s], T.sig[s], T.ok[s], T.infl[s], gate, C, tdr::fixesTimes(s)).empty();
    }
    o.jvf_t0[s] = jvfOf(b, trk, hsT, nearPu);
  }
  // R_pT, as the TDR taggers (ghost matching)
  const tdr::JetRpt r = tdr::jetRpt(b, j, L, T, gate, C, tdr::RptMatch::GHOST);
  o.rpt_itk = r.itk;
  for (int s = 0; s < tdr::N_T0; ++s) o.rpt_t0[s] = r.t0[s];
  // jet time over the timed tracks: any vertex (the note's), and PV0-associated only
  std::vector<std::pair<double, double>> tw, twHs;   // (time, weight)
  for (size_t k = 0; k < trk.size(); ++k) {
    const int i = trk[k];
    if (!tdr::timed(b, i, C, C.jets()) || !(b.trackTimeRes[i] > 0.0)) continue;
    const double w = 1.0 / (b.trackTimeRes[i] * b.trackTimeRes[i]);
    tw.emplace_back(b.trackTime[i], w);
    if (hsZ[k]) twHs.emplace_back(b.trackTime[i], w);
    if (b.trackPt[i] > o.pt_lead) { o.pt_lead = b.trackPt[i]; o.t_lead = b.trackTime[i]; o.sig_lead = b.trackTimeRes[i]; }
  }
  o.ntimed = (int)tw.size();
  o.ntimed_hs = (int)twHs.size();
  timeStats(tw, o.t_jet, o.var_spread, o.var_mean);
  timeStats(twHs, o.t_hs, o.var_hs_spread, o.var_hs_mean);
  // the timing veto's core jet time
  std::vector<int> cand;
  for (int idx : b.topoJetGhostTrackIdx[j]) {
    if (idx < 0 || idx >= nT) continue;
    if (!tdr::timed(b, idx, C, C.jets()) || !(b.trackTimeRes[idx] > 0.0)) continue;
    if (!passTrackKinematics((size_t)idx, &b, MIN_TRACK_PT, MAX_TRACK_PT)) continue;
    if (tdr::dR(o.eta, o.phi, b.trackEta[idx], b.trackPhi[idx]) >= CORE_DR) continue;
    cand.push_back(idx);
  }
  auto core = [&](const std::vector<int>& tracks, double& t, double& sig, int& n) {
    if (tracks.empty()) return;
    static const std::unordered_map<int, double> kNone;
    std::vector<Cluster> cl = makeSimpleClusters(tracks, &b, /*useSmearedTimes=*/false, kNone, kNone,
                                                 /*checkTimeValid=*/true, /*usez0=*/false);
    doIterativeClustering(&cl, CORE_DIST_CUT);
    const Cluster* best = nullptr;   // TRKPT score = sum pT
    for (const Cluster& c : cl)
      if (!best || c.scores.at(Score::TRKPT.id) > best->scores.at(Score::TRKPT.id)) best = &c;
    t = best->values[0];  sig = best->sigmas[0];  n = (int)best->trackIndices.size();
  };
  core(cand, o.core_t, o.core_sig, o.core_n);
  std::vector<int> candFix;
  for (int idx : cand) if (!tdr::misassignedHs(b, idx)) candFix.push_back(idx);
  if (candFix.size() == cand.size()) { o.core_t_fix = o.core_t;  o.core_sig_fix = o.core_sig;  o.core_n_fix = o.core_n; }
  else core(candFix, o.core_t_fix, o.core_sig_fix, o.core_n_fix);
  // the core-time variants (CV_NAMES)
  for (int v = 0; v < 4; ++v) {
    const double drMax = (v == 0 || v == 2) ? 0.4 : CORE_DR;
    const double ptMax = (v == 1 || v == 2) ? std::numeric_limits<double>::infinity() : MAX_TRACK_PT;
    std::vector<int> cv;
    for (int idx : b.topoJetGhostTrackIdx[j]) {
      if (idx < 0 || idx >= nT) continue;
      if (!tdr::timed(b, idx, C, C.jets()) || !(b.trackTimeRes[idx] > 0.0)) continue;
      if (!passTrackKinematics((size_t)idx, &b, MIN_TRACK_PT, ptMax)) continue;
      if (tdr::dR(o.eta, o.phi, b.trackEta[idx], b.trackPhi[idx]) >= drMax) continue;
      if (v == 3 && !assocZ(b, idx, b.recoVtxZ[0], /*fwd=*/true)) continue;
      cv.push_back(idx);
    }
    core(cv, o.cv_t[v], o.cv_sig[v], o.cv_n[v]);
  }
  if (!cand.empty()) {   // zpref: the standard candidates, the cluster with the most z-compatible pT preferred
    static const std::unordered_map<int, double> kNone;
    std::vector<Cluster> cl = makeSimpleClusters(cand, &b, /*useSmearedTimes=*/false, kNone, kNone,
                                                 /*checkTimeValid=*/true, /*usez0=*/false);
    doIterativeClustering(&cl, CORE_DIST_CUT);
    const Cluster* best = nullptr;
    double bestZ = 0.0;
    for (const Cluster& c : cl) {
      double zpt = 0.0;
      for (int t : c.trackIndices) if (assocZ(b, t, b.recoVtxZ[0], /*fwd=*/true)) zpt += b.trackPt[t];
      if (zpt > bestZ) { bestZ = zpt; best = &c; }
    }
    if (!best)   // no z-compatible track anywhere: the standard choice
      for (const Cluster& c : cl)
        if (!best || c.scores.at(Score::TRKPT.id) > best->scores.at(Score::TRKPT.id)) best = &c;
    o.cv_t[4] = best->values[0];  o.cv_sig[4] = best->sigmas[0];  o.cv_n[4] = (int)best->trackIndices.size();
  }
  return o;
}
}  // namespace

int main(int argc, char** argv) {
  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL    = cfg.energyLabel;
  OUTPUT_DIR      = cfg.outputDir;
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  FILE_SHARD      = resolveShard(argc, argv);
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);

  double gate = tdr::GATE_SIGMA_TDR;
  std::string outOverride, pairMode = "ff", jvtName = "none";
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--gate-sigma=", 0) == 0) gate = std::stod(a.substr(13));
    if (a.rfind("--out=", 0) == 0) outOverride = a.substr(6);
    if (a.rfind("--pair=", 0) == 0) pairMode = a.substr(7);
    if (a.rfind("--jvt=", 0) == 0) jvtName = a.substr(6);
  }
  if (!(gate > 0.0)) { std::cerr << "[tdr_fwdpair] --gate-sigma must be > 0\n"; return 1; }
  if (pairMode != "ff" && pairMode != "one") { std::cerr << "[tdr_fwdpair] --pair must be ff or one\n"; return 1; }
  const bool oneFwd = pairMode == "one";
  if (oneFwd) SELECTION_TAG = SELECTION_TAG.empty() ? "onefwd" : SELECTION_TAG + "_onefwd";
  const JvtWP* jwp = jvtWPByName(jvtName);
  if (!jwp || jvtName == "looseext") { std::cerr << "[tdr_fwdpair] --jvt must be none, loose, tight, jvt4pt60, jvt4pt120, jvt4pt60t0, jvt4pt60t0tru, jvt4pt60f03, jvt4pt60f05, jvt4pt60f10, jvt4pt60tdr, jvt4pt60tdrw or jvt4pt60tdrtru\n"; return 1; }
  const bool jvtOn = jwp->rptMin >= 0.0;
  EXTENDED_BRANCHES = true;   // Track_recoVtx_idx: the JVT discriminants, and the waves_pv0 / waves_npu t0s (tdr_common.h)
  if (jvtOn) {
    const std::string tag = std::string(jwp->tag).substr(1);
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  if (gate != tdr::GATE_SIGMA_TDR) {
    std::ostringstream os; os << "gate" << std::fixed << std::setprecision(1) << gate;
    std::string tag = os.str(); std::replace(tag.begin(), tag.end(), '.', 'p');
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  const tdr::CleanMode cleanMode = tdr::cleaningMode(argc, argv);   // before the output path: it tags the file
  tdr::rptTrkPtMaxMode(argc, argv);                                  // likewise
  if (jvtOn && jwp->tdrRpt > -2 && std::isfinite(tdr::RPT_TRK_PT_MAX)) {
    std::cerr << "[tdr_fwdpair] --jvt=" << jvtName << " is calibrated on the uncapped TDR R_pT: pass --rpt-trk-ptmax=none\n";
    return 1;
  }
  const Inflation INFL = inflationFor(SAMPLE_NAME);
  const std::string outPath = outOverride.empty() ? histFilePath("tdr_fwdpair.root") : outOverride;
  {
    const auto dir = boost::filesystem::path(outPath).parent_path();
    if (!dir.empty()) boost::filesystem::create_directories(dir);
  }
  std::printf("[tdr_fwdpair] sample=%s: two leading jets in %s|eta| < %.1f, opposite eta, > %.0f / %.0f GeV%s; gate %.1f sigma\n",
              SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME.c_str(), oneFwd ? "" : "2.4 < ", ETA_MAX, J1_PT_MIN, J2_PT_MIN,
              oneFwd ? ", at least one in 2.4 < |eta| < 3.8" : "", gate);
  printInflation(SAMPLE_NAME);
  std::cout << "[tdr_fwdpair] track-time cleaning: " << tdr::cleanModeName(cleanMode) << "; JVT/fJVT before pairing: " << jvtName << "\n";
  std::cout << "[tdr_fwdpair] out=" << outPath << "\n";

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), FILE_SHARD);
  recordAvailableBranches(chain);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  const tdr::TimeCleaning CLEAN = tdr::setupCleaning(reader, cleanMode);

  TFile out(outPath.c_str(), "RECREATE");
  if (out.IsZombie()) { std::cerr << "[tdr_fwdpair] cannot open " << outPath << "\n"; return 1; }
  TTree tree("events", "PUB-2022-001 Sec. 6: events with their two leading forward jets in 2.4 < |eta| < 3.8, opposite eta");
  int file_idx = -1, jz = -1;
  Long64_t entry = -1;
  float weight = 0.f;
  double mjj = 0.0, deta = 0.0;
  int cat = 0;
  double t0[tdr::N_T0] = {}, sig0[tdr::N_T0] = {}, infl[tdr::N_T0] = {};
  bool ok0[tdr::N_T0] = {};
  JetOut J[2];
  tree.Branch("file_idx", &file_idx, "file_idx/I");
  tree.Branch("jz", &jz, "jz/I");
  tree.Branch("entry", &entry, "entry/L");
  tree.Branch("weight", &weight, "weight/F");
  tree.Branch("mjj", &mjj, "mjj/D");
  tree.Branch("deta", &deta, "deta/D");
  tree.Branch("cat", &cat, "cat/I");   // forward jets in the pair: 1 = central-forward, 2 = forward-forward
  for (int k = 0; k < tdr::N_T0; ++k) {
    const std::string s = tdr::T0_NAMES[k];
    tree.Branch(("t0_" + s).c_str(), &t0[k], ("t0_" + s + "/D").c_str());
    tree.Branch(("sig0_" + s).c_str(), &sig0[k], ("sig0_" + s + "/D").c_str());
    tree.Branch(("infl_" + s).c_str(), &infl[k], ("infl_" + s + "/D").c_str());
    tree.Branch(("ok_" + s).c_str(), &ok0[k], ("ok_" + s + "/O").c_str());
  }
  for (int n = 0; n < 2; ++n) {
    const std::string p = "j" + std::to_string(n + 1) + "_";
    JetOut& o = J[n];
    tree.Branch((p + "pt").c_str(), &o.pt, (p + "pt/F").c_str());
    tree.Branch((p + "eta").c_str(), &o.eta, (p + "eta/F").c_str());
    tree.Branch((p + "phi").c_str(), &o.phi, (p + "phi/F").c_str());
    tree.Branch((p + "hs").c_str(), &o.hs, (p + "hs/I").c_str());
    tree.Branch((p + "pu").c_str(), &o.pu, (p + "pu/I").c_str());
    tree.Branch((p + "ntrk").c_str(), &o.ntrk, (p + "ntrk/I").c_str());
    tree.Branch((p + "ntimed").c_str(), &o.ntimed, (p + "ntimed/I").c_str());
    tree.Branch((p + "jvf_itk").c_str(), &o.jvf_itk, (p + "jvf_itk/D").c_str());
    tree.Branch((p + "rpt_itk").c_str(), &o.rpt_itk, (p + "rpt_itk/D").c_str());
    for (int k = 0; k < tdr::N_T0; ++k) {
      const std::string s = tdr::T0_NAMES[k];
      tree.Branch((p + "jvf_t0_" + s).c_str(), &o.jvf_t0[k], (p + "jvf_t0_" + s + "/D").c_str());
      tree.Branch((p + "rpt_t0_" + s).c_str(), &o.rpt_t0[k], (p + "rpt_t0_" + s + "/D").c_str());
    }
    tree.Branch((p + "t_jet").c_str(), &o.t_jet, (p + "t_jet/D").c_str());
    tree.Branch((p + "var_spread").c_str(), &o.var_spread, (p + "var_spread/D").c_str());
    tree.Branch((p + "var_mean").c_str(), &o.var_mean, (p + "var_mean/D").c_str());
    tree.Branch((p + "t_lead").c_str(), &o.t_lead, (p + "t_lead/D").c_str());
    tree.Branch((p + "sig_lead").c_str(), &o.sig_lead, (p + "sig_lead/D").c_str());
    tree.Branch((p + "pt_lead").c_str(), &o.pt_lead, (p + "pt_lead/D").c_str());
    tree.Branch((p + "ntimed_hs").c_str(), &o.ntimed_hs, (p + "ntimed_hs/I").c_str());
    tree.Branch((p + "t_hs").c_str(), &o.t_hs, (p + "t_hs/D").c_str());
    tree.Branch((p + "var_hs_spread").c_str(), &o.var_hs_spread, (p + "var_hs_spread/D").c_str());
    tree.Branch((p + "var_hs_mean").c_str(), &o.var_hs_mean, (p + "var_hs_mean/D").c_str());
    tree.Branch((p + "core_t").c_str(), &o.core_t, (p + "core_t/D").c_str());
    tree.Branch((p + "core_sig").c_str(), &o.core_sig, (p + "core_sig/D").c_str());
    tree.Branch((p + "core_n").c_str(), &o.core_n, (p + "core_n/I").c_str());
    tree.Branch((p + "core_t_fix").c_str(), &o.core_t_fix, (p + "core_t_fix/D").c_str());
    tree.Branch((p + "core_sig_fix").c_str(), &o.core_sig_fix, (p + "core_sig_fix/D").c_str());
    tree.Branch((p + "core_n_fix").c_str(), &o.core_n_fix, (p + "core_n_fix/I").c_str());
    for (int v = 0; v < JetOut::N_CV; ++v) {
      const std::string sfx = std::string("_") + CV_NAMES[v];
      tree.Branch((p + "core_t" + sfx).c_str(), &o.cv_t[v], (p + "core_t" + sfx + "/D").c_str());
      tree.Branch((p + "core_sig" + sfx).c_str(), &o.cv_sig[v], (p + "core_sig" + sfx + "/D").c_str());
      tree.Branch((p + "core_n" + sfx).c_str(), &o.cv_n[v], (p + "core_n" + sfx + "/I").c_str());
    }
    tree.Branch((p + "jvt_rpt").c_str(), &o.jvt_rpt, (p + "jvt_rpt/F").c_str());
    tree.Branch((p + "fjvt").c_str(), &o.fjvt, (p + "fjvt/F").c_str());
    tree.Branch((p + "jvt_rpt_t0").c_str(), &o.jvt_rpt_t0, (p + "jvt_rpt_t0/F").c_str());
  }

  const std::regex jzRe("JZ([0-9])");
  Long64_t nRead = 0, nVtx = 0, nLep = 0, nTwo = 0, nSel = 0, nRmJvt = 0, nRmFjvt = 0;
  std::vector<JetDisc> D;
  double sumw = 0.0;
  while (reader.Next()) {
    if (maxEvents > 0 && nRead >= maxEvents) break;
    ++nRead;
    sumw += *branch.weight;
    if (nRead % 100000 == 0) std::cout << "[tdr_fwdpair] read " << nRead << ", selected " << nSel << std::endl;
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    if (std::abs(branch.recoVtxZ[0] - branch.truthVtxZ[0]) > MAX_VTX_DZ) continue;
    ++nVtx;
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap()) continue;
    ++nLep;

    if (jvtOn) D = computeJetDiscriminants(branch);   // for every jet: fJVT needs the central jets' vertices
    // a t0-gated JVT needs the t0 BEFORE pairing (the t0s do not depend on which jets JVT keeps)
    tdr::T0Set T;
    bool haveT = false;
    std::vector<float> rptPlain;
    if (jvtOn && jwp->t0Gate >= 0) {
      T = tdr::computeT0s(branch, INFL, CLEAN);  haveT = true;
      rptPlain.resize(D.size());
      for (int j = 0; j < (int)D.size(); ++j) {
        rptPlain[j] = D[j].rpt;
        D[j].rpt = (float)tdr::jvtRpt(branch, j, T, gate, CLEAN).t0[jwp->t0Gate];
      }
    }
    if (jvtOn && jwp->tdrRpt > -2) {   // the TDR R_pT (ITk or t0-gated) on the forward window's jets, JVT's elsewhere
      if (!haveT) { T = tdr::computeT0s(branch, INFL, CLEAN);  haveT = true; }
      const tdr::TrackLists Lp = tdr::rptTrackLists(branch);
      rptPlain.resize(D.size());
      for (int j = 0; j < (int)D.size(); ++j) {
        rptPlain[j] = D[j].rpt;
        if (!D[j].rptFwdWin) continue;
        const tdr::JetRpt r = tdr::jetRpt(branch, j, Lp, T, gate, CLEAN, tdr::RptMatch::GHOST);
        D[j].rpt = (float)(jwp->tdrRpt < 0 ? r.itk : r.t0[jwp->tdrRpt]);
      }
    }
    std::vector<int> js;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j) || branch.topoJetPt[j] < tdr::JET_PT_MIN) continue;
      if (jvtOn) {                                  // the tagger before pairing: a failing jet is dropped
        if (jvtRemoves(D[j], *jwp))  { ++nRmJvt;  continue; }
        if (fjvtRemoves(D[j], *jwp)) { ++nRmFjvt; continue; }
      }
      const double ae = std::abs((double)branch.topoJetEta[j]);
      if (ae < ETA_MAX && (oneFwd || ae > ETA_MIN)) js.push_back(j);
    }
    if (js.size() < 2) continue;
    ++nTwo;
    std::sort(js.begin(), js.end(), [&](int x, int y) { return branch.topoJetPt[x] > branch.topoJetPt[y]; });
    const int j1 = js[0], j2 = js[1];
    if (branch.topoJetPt[j1] <= J1_PT_MIN || branch.topoJetPt[j2] <= J2_PT_MIN) continue;
    if ((double)branch.topoJetEta[j1] * (double)branch.topoJetEta[j2] >= 0.0) continue;
    cat = (int)(std::abs((double)branch.topoJetEta[j1]) > ETA_MIN) + (int)(std::abs((double)branch.topoJetEta[j2]) > ETA_MIN);
    if (cat == 0) continue;                       // --pair=one: at least one jet in HGTD (always true for ff)
    ++nSel;

    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(branch.topoJetPt[j1], branch.topoJetEta[j1], branch.topoJetPhi[j1], 0.0);
    v2.SetPtEtaPhiM(branch.topoJetPt[j2], branch.topoJetEta[j2], branch.topoJetPhi[j2], 0.0);
    mjj = (v1 + v2).M();
    deta = std::abs((double)branch.topoJetEta[j1] - (double)branch.topoJetEta[j2]);

    if (!haveT) T = tdr::computeT0s(branch, INFL, CLEAN);
    const tdr::TrackLists L = tdr::rptTrackLists(branch);
    for (int k = 0; k < tdr::N_T0; ++k) { t0[k] = T.t[k]; sig0[k] = T.sig[k]; infl[k] = T.infl[k]; ok0[k] = T.ok[k]; }
    J[0] = jetQuantities(branch, j1, L, T, gate, CLEAN);
    J[1] = jetQuantities(branch, j2, L, T, gate, CLEAN);
    if (jvtOn) {
      J[0].jvt_rpt = rptPlain.empty() ? D[j1].rpt : rptPlain[j1];  J[0].fjvt = D[j1].fjvt;
      J[1].jvt_rpt = rptPlain.empty() ? D[j2].rpt : rptPlain[j2];  J[1].fjvt = D[j2].fjvt;
      J[0].jvt_rpt_t0 = (float)tdr::jvtRpt(branch, j1, T, gate, CLEAN).t0[tdr::T0_WAVES];
      J[1].jvt_rpt_t0 = (float)tdr::jvtRpt(branch, j2, T, gate, CLEAN).t0[tdr::T0_WAVES];
    }

    const std::string fpath = chain.GetCurrentFile()->GetName();
    { const size_t pp = fpath.rfind('_'), d = (pp == std::string::npos) ? pp : fpath.find('.', pp);
      file_idx = (pp != std::string::npos && d != std::string::npos) ? std::atoi(fpath.substr(pp + 1, d - pp - 1).c_str()) : -1; }
    std::smatch m;
    jz = std::regex_search(fpath, m, jzRe) ? std::stoi(m[1].str()) : -1;   // dijet slice; -1 elsewhere
    entry = chain.GetTree()->GetReadEntry();
    weight = (float)*branch.weight;
    tree.Fill();
  }
  if (nRead == 0) {
    std::cerr << "[tdr_fwdpair] ERROR: read 0 events -- check stderr for TTreeReader branch errors\n";
    return 1;
  }

  TTree meta("meta", "run record (one row per job; sum over shards)");
  std::string sample = SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME, ndir = cfg.ntupleDir;
  std::string cleanName = tdr::cleanModeName(cleanMode);
  meta.Branch("n_read", &nRead, "n_read/L");
  meta.Branch("sumw_read", &sumw, "sumw_read/D");
  meta.Branch("n_pass_vtx", &nVtx, "n_pass_vtx/L");
  meta.Branch("n_pass_lepton", &nLep, "n_pass_lepton/L");
  meta.Branch("n_two_fwd", &nTwo, "n_two_fwd/L");
  meta.Branch("n_selected", &nSel, "n_selected/L");
  meta.Branch("gate_sigma", &gate, "gate_sigma/D");
  double rptPtMax = tdr::RPT_TRK_PT_MAX;
  meta.Branch("rpt_trk_ptmax", &rptPtMax, "rpt_trk_ptmax/D");
  meta.Branch("pair_mode", &pairMode);
  meta.Branch("jvt_wp", &jvtName);
  meta.Branch("n_rm_jvt", &nRmJvt, "n_rm_jvt/L");
  meta.Branch("n_rm_fjvt", &nRmFjvt, "n_rm_fjvt/L");
  meta.Branch("time_cleaning", &cleanName);
  meta.Branch("sample", &sample);
  meta.Branch("ntuple_dir", &ndir);
  meta.Fill();
  out.Write();
  out.Close();
  std::printf("[tdr_fwdpair] read %lld, vertex %lld, lepton/OR %lld, two jets in the window %lld, selected %lld\n",
              (long long)nRead, (long long)nVtx, (long long)nLep, (long long)nTwo, (long long)nSel);
  std::cout << "[tdr_fwdpair] wrote " << outPath << "\n";
  return 0;
}
