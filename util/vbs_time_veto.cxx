// -----------------------------------------------------------------------------
// vbs_time_veto -- event-level timing veto on the VBS tagging-jet pair:
// jet-vs-jet time compatibility against jet-vs-t0 compatibility.
//
// Feedback of 2026-09-25: in R1 (both tagging jets forward, one hard-scatter
// and one pileup) no vertex t0 is needed to reject the event. The two jets can
// be tested against EACH OTHER: give each a time from its own timed tracks and
// reject the event when
//
//     |t_A - t_B| / sqrt(sigma_A^2 + sigma_B^2) >= 3.
//
// That needs BOTH legs timed, so it can only act on pairs whose legs are both
// in HGTD acceptance (the R1 geometry). The vertex-t0 alternative tests EACH
// timed leg against the event t0 instead, and rejects when either leg fails
//
//     |t_leg - t0| / sqrt(sigma_leg^2 + (f sigma_t0)^2) >= 3,
//
// so it can also reject R2 (forward PU + central HS) and any other pair with a
// timed pileup leg -- and can also lose a genuine pair to a wrong t0.
//
// The question is what each costs in SIGNAL (VBF H->inv) and buys in
// BACKGROUND (Z+jets), weighted, after the analysis-level recoil cut. This
// program writes the INGREDIENTS of both tests for every event passing the
// selection -- never the decisions -- so the threshold, the jet-time
// calibration and the recoil cut are all applied offline, with weights, by
// python/vbs_time_veto_plot.py.
//
// ── Selection ────────────────────────────────────────────────────────────────
// Exactly the composition plots' (vbs_region_diag + vbs_region_stack.py), so an
// efficiency here describes the same events as a region-plot column:
//   vertex |z_reco - z_truth| < MAX_VTX_DZ, Z->ll + lepton-jet OR (Z+jets only),
//   >= 2 jets > 30 GeV and >= 1 of them in 2.38 < |eta| < 4.0, re-checked after
//   the --jvt tagger removes jets, then the max-m_jj opposite-hemisphere pair
//   of the surviving jets, m_jj >= VBS_JET_MJJ (default 200 GeV; --vbs-mjj=)
//   and |Deta| >= VBS_JET_D_ETA (default 0; --vbs-deta=).
// The recoil cut is NOT applied here; both recoil variables are stored:
//   met_truth  |sum pT(nu)| over the truth record (pdgId 12/14/16). The skims
//              DROP the truth record, so VBF must be read from the ORIGINAL
//              ntuples (--ntuple-dir=.../highstats_vbf/) to get it; -1 when
//              the sample has no truth record.
//   z_pt       the dilepton pT, vbs_region_diag's definition; -1 unless Z+jets.
//
// ── Jet time ─────────────────────────────────────────────────────────────────
// Two estimators per leg, both from the jet's ghost-associated tracks that
// carry a valid HGTD time:
//
//   all_   the prescription read literally: EVERY timed ghost track, inverse-
//          variance mean. Kept for the record, but at mu = 200 it is a pileup
//          average: a forward jet's ghost area holds ~23 timed tracks, and on
//          local VBF only 21% of an HS jet's are from the hard scatter
//          (chi2/ndf ~ 65). Only 25% of HS+HS pairs pass the 3-sigma test with
//          it, barely more than R1 pairs (21%).
//   core_  the default: timed ghost tracks passing the standard track
//          kinematics (passTrackKinematics: 1 < pT < 30 GeV, quality flag,
//          HGTD eta) within dR < JET_TIME_DR of the jet axis (the R_pT cone),
//          clustered in time (doIterativeClustering at JET_TIME_DIST_CUT); the
//          jet time is the inverse-variance mean of the highest-sum-pT cluster.
//          Raises the HS-track share of the estimate to ~73% and the HS+HS
//          pass rate to 83% (raw sigma; both legs in 2.4-4.0) on local VBF.
// The quoted sigma of either is a pure statistical error: it knows nothing of
// the pileup tracks left in the cluster or of mis-assigned hits, so -- like the
// vertex-t0 sigma -- it is too small, by ~1.5 for core_. The calibration is
// measured from the paper-HS legs (printed below, and by the plotting script)
// and applied OFFLINE, like rpt_v5's t0 inflation.
//
// ── t0 ───────────────────────────────────────────────────────────────────────
// Identical code and inputs to rpt_v5_hist, so every t0 here equals its region
// tree's t0 for the same event: hgtd (Athena RecoVtx_time, where valid), trkptz,
// waves (in-jet refined), tzp (guarded in-jet), each with sigma before
// inflation and the per-sample inflation rpt_v5 uses (inflationFor, now in
// util/rpt_v5_common.h). Plus truth: t0 = TruthVtx_time[0], sigma 0 -- a
// perfect vertex time, i.e. the jet-time-limited ceiling of the t0 method.
//
// ── Outputs (histFilePath: sample-, selection- and shard-tagged) ─────────────
//   <prefix>vbs_time_veto.root
//     events  one row per selected event (see Row): the legs' columns, and the
//             jet_* arrays of every jet the pair was chosen from, so a timing
//             tagger can drop jets and re-form the pair offline (--repair in
//             python/vbs_time_veto_plot.py)
//     meta    one row per job: cut flow and sum of weights of EVERY event read
// TTrees only: merge shards with hadd (NOT hist_merge), and sum meta's rows.
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "clustering_constants.h"
#include "sample_config.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "jet_tagging.h"
#include "rpt_v5_common.h"   // inflationFor(): rpt_v5's per-sample t0 calibration

using namespace MyUtl;

namespace {

// core_ estimator constants. 0.2 is the R_pT cone (rpt_v5_hist's
// RPT_TRACK_JET_DR); 2 sigma is the in-jet re-clustering distance WAVES_RECLUST
// uses (DIST_CUT_T_REFINED). Chosen on local VBF truth, 2026-09-28, against
// dR 0.4 and 3 sigma: dR 0.2 is worth ~6 points of HS+HS pass rate, the cut
// ~1 point, and picking the cluster by sum pT rather than track count ~1 point.
constexpr double JET_TIME_DR       = 0.2;
constexpr double JET_TIME_DIST_CUT = 2.0;

// The clustering input for the t0 rows -- rpt_v5_hist's, verbatim.
constexpr double T0_ASSOC_NSIGMA = 2.5;

struct JetTime {
  double t = 0.0, sig = -1.0;   // sig < 0: no time
  int    n = 0;                 // tracks in the estimate
  int    nCand = 0;             // timed tracks it was built from
  double sumPt = 0.0;           // sum pT of the n tracks
  double hsFrac = -1.0;         // TRUTH diagnostic: their pT share from truth vertex 0
};

// all_: every timed ghost track, inverse-variance mean.
JetTime jetTimeAllGhost(BranchPointerWrapper& b, int j) {
  JetTime out;
  const int nT = (int)b.trackPt.GetSize();
  double sw = 0.0, swt = 0.0, hs = 0.0;
  for (int idx : b.topoJetGhostTrackIdx[j]) {
    if (idx < 0 || idx >= nT || b.trackTimeValid[idx] != 1) continue;
    const double s = b.trackTimeRes[idx];
    if (!(s > 0.0)) continue;
    sw  += 1.0 / (s * s);
    swt += b.trackTime[idx] / (s * s);
    ++out.n;
    out.sumPt += b.trackPt[idx];
    if (b.trackToTruthvtx[idx] == 0) hs += b.trackPt[idx];
  }
  out.nCand = out.n;
  if (out.n == 0) return out;
  out.t = swt / sw;
  out.sig = std::sqrt(1.0 / sw);
  out.hsFrac = out.sumPt > 0.0 ? hs / out.sumPt : -1.0;
  return out;
}

// core_: kinematic + dR-cone selection, time clustering, highest-sum-pT cluster.
JetTime jetTimeCore(BranchPointerWrapper& b, int j) {
  JetTime out;
  const int nT = (int)b.trackPt.GetSize();
  const double je = b.topoJetEta[j], jp = b.topoJetPhi[j];
  std::vector<int> cand;
  for (int idx : b.topoJetGhostTrackIdx[j]) {
    if (idx < 0 || idx >= nT || b.trackTimeValid[idx] != 1) continue;
    if (!(b.trackTimeRes[idx] > 0.0)) continue;
    if (!passTrackKinematics((size_t)idx, &b, MIN_TRACK_PT, MAX_TRACK_PT)) continue;
    const double dr = std::hypot(b.trackEta[idx] - je, TVector2::Phi_mpi_pi(b.trackPhi[idx] - jp));
    if (dr >= JET_TIME_DR) continue;
    cand.push_back(idx);
  }
  out.nCand = (int)cand.size();
  if (cand.empty()) return out;
  static const std::unordered_map<int, double> kNone;
  std::vector<Cluster> cl = makeSimpleClusters(cand, &b, /*useSmearedTimes=*/false, kNone, kNone,
                                               /*checkTimeValid=*/true, /*usez0=*/false);
  doIterativeClustering(&cl, JET_TIME_DIST_CUT);
  const Cluster* best = nullptr;   // TRKPT score = sum pT (summed on every merge)
  for (const Cluster& c : cl)
    if (!best || c.scores.at(Score::TRKPT.id) > best->scores.at(Score::TRKPT.id)) best = &c;
  out.t     = best->values[0];
  out.sig   = best->sigmas[0];
  out.n     = (int)best->trackIndices.size();
  out.sumPt = best->scores.at(Score::TRKPT.id);
  double hs = 0.0;
  for (int idx : best->trackIndices)
    if (b.trackToTruthvtx[idx] == 0) hs += b.trackPt[idx];
  out.hsFrac = out.sumPt > 0.0 ? hs / out.sumPt : -1.0;
  return out;
}

// Per-leg columns.
struct Leg {
  int    idx = -1;
  float  pt = 0.f, eta = 0.f, phi = 0.f;
  bool   hs = false, pu = false;          // paper labels (isJetPaperHS / isJetPaperPU)
  float  jvt_rpt = -1.f, fjvt = -1.f;     // the tagger's discriminants; -1 with --jvt=none
  double core_t = 0.0, core_sig = -1.0, core_sumpt = 0.0, core_hsfrac = -1.0;
  int    core_n = 0, core_ncand = 0;
  double all_t = 0.0, all_sig = -1.0, all_hsfrac = -1.0;
  int    all_n = 0;
};

// The t0 sources, in column order. truth has sigma 0 and inflation 1.
constexpr int N_T0 = 5;
const char* const T0_NAMES[N_T0] = {"hgtd", "trkptz", "waves", "tzp", "truth"};

struct Row {
  std::string file_path;
  Long64_t entry = -1;         // entry within file_path (event displays, weight lookups)
  Long64_t chain_entry = -1;   // entry within this job's chain: vbs_region_diag's `entry`
  int      file_idx = -1;      // vbs_region_diag's file_idx (trailing number of the name)
  float    weight = 0.f;
  float    met_truth = -1.f, z_pt = -1.f, z_mass = -1.f;
  int      n_jets = 0, n_jets_fwd_acc = 0, n_rm_jvt = 0, n_rm_fjvt = 0;
  double   mjj = -1.0, deta = -1.0;
  Leg      a, b;
  double   t0[N_T0] = {}, sig0[N_T0] = {}, infl[N_T0] = {};
  bool     ok[N_T0] = {};
  // EVERY jet the pair was chosen from -- the pT-passing, not-overlap-removed
  // jets that survived --jvt -- so a timing tagger can remove jets and the pair
  // be re-formed offline. Removing jets can only lower the maximum m_jj and the
  // jet counts, so these events are a superset of any re-paired selection. The
  // legs are the entries whose jet_idx equals a_idx / b_idx.
  std::vector<int>    jet_idx, jet_hs, jet_pu, jet_n;
  std::vector<float>  jet_pt, jet_eta, jet_phi;
  std::vector<double> jet_t, jet_sig;   // core_ jet time; jet_n = 0: no time
};

// Jet-time calibration accumulator over paper-HS legs (see the printout).
struct Calib {
  std::vector<double> dt, sig;
  void add(const JetTime& jt, double tHS) {
    if (jt.n == 0) return;
    dt.push_back(jt.t - tHS);
    sig.push_back(jt.sig);
  }
  void print(const char* name) const {
    if (dt.size() < 50) { std::printf("  %-6s %8zu   (too few)\n", name, dt.size()); return; }
    // (1) rpt_v5's t0 definition: RMS of dt about its mean inside |dt| < 150 ps,
    //     over the RMS quoted sigma of the same entries.
    double s1 = 0, s2 = 0, sv = 0; long n = 0;
    for (size_t i = 0; i < dt.size(); ++i) {
      if (std::abs(dt[i]) >= 150.0) continue;
      s1 += dt[i]; s2 += dt[i] * dt[i]; sv += sig[i] * sig[i]; ++n;
    }
    const double mean = s1 / n, rms = std::sqrt(std::max(0.0, s2 / n - mean * mean));
    const double quo = std::sqrt(sv / n);
    // (2) per-jet pulls: 1.4826 x median |pull| (the Gaussian-core width,
    //     tail-insensitive) and the fraction inside 3 at the raw sigma.
    std::vector<double> ap(dt.size());
    long in3 = 0;
    for (size_t i = 0; i < dt.size(); ++i) {
      ap[i] = std::abs(dt[i]) / sig[i];
      if (ap[i] < 3.0) ++in3;
    }
    std::nth_element(ap.begin(), ap.begin() + ap.size() / 2, ap.end());
    const double mad = 1.4826 * ap[ap.size() / 2];
    std::printf("  %-6s %8zu %8.1fps %8.1fps %8.1fps %7.2f %9.2f %8.1f%% %8.1f%%\n",
                name, dt.size(), mean, rms, quo, quo > 0 ? rms / quo : 0.0, mad,
                100.0 * in3 / dt.size(), 100.0 * (dt.size() - n) / dt.size());
  }
};

}  // namespace

int main(int argc, char** argv) {
  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL    = cfg.energyLabel;
  OUTPUT_DIR      = cfg.outputDir;
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  resolveSelection(argc, argv);   // --vbs-mjj / --vbs-deta; tags the output name
  if (VBS_JET_D_ETA < 0.0) {
    std::cerr << "[veto] --vbs-deta < 0 drops the VBS pair, and both tests act ON the pair\n";
    return 1;
  }
  // MUST be set explicitly: setupChain reads it, and nothing else populates it.
  FILE_SHARD = resolveShard(argc, argv);
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);

  // --jvt=none|loose|tight: pileup-jet tagging before pairing, exactly as
  // vbs_region_diag / rpt_v5_hist (src/jet_tagging.h). Joins SELECTION_TAG.
  const JvtWP* wpSel = &JVT_WPS[0];
  std::string outOverride;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--out=", 0) == 0) outOverride = a.substr(6);
    if (a.rfind("--jvt=", 0) != 0) continue;
    const JvtWP* hit = jvtWPByName(a.substr(6));
    if (!hit) { std::cerr << "[veto] unknown " << a << " (none|loose|tight)\n"; return 1; }
    wpSel = hit;
  }
  const JvtWP& WP = *wpSel;
  const bool jvtOn = WP.rptMin >= 0.0;
  if (jvtOn) {
    // The discriminants need the vertex fit's track assignment (EXTENDED).
    EXTENDED_BRANCHES = true;
    const std::string tag = std::string(WP.tag).substr(1);
    SELECTION_TAG = SELECTION_TAG.empty() ? tag : SELECTION_TAG + "_" + tag;
  }
  const Inflation INFL = inflationFor(SAMPLE_NAME);

  const std::string outPath = outOverride.empty() ? histFilePath("vbs_time_veto.root") : outOverride;
  boost::filesystem::create_directories(boost::filesystem::path(outPath).parent_path().empty()
                                            ? boost::filesystem::path(".")
                                            : boost::filesystem::path(outPath).parent_path());
  std::cout << "[veto] sample=" << (SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME)
            << " jvt=" << WP.name << " m_jj >= " << VBS_JET_MJJ << " |Deta| >= " << VBS_JET_D_ETA
            << "  jet time: core = timed ghost tracks, " << MIN_TRACK_PT << "-" << MAX_TRACK_PT
            << " GeV, dR < " << JET_TIME_DR << ", " << JET_TIME_DIST_CUT
            << "-sigma time clustering, max sum-pT cluster\n";
  std::printf("[veto] t0 inflation (rpt_v5's, '%s'): hgtd %.2f trkptz %.2f waves %.2f tzp %.2f\n",
              SAMPLE_NAME.empty() ? "local (vbf)" : SAMPLE_NAME.c_str(),
              INFL.hgtd, INFL.trkptz, INFL.waves, INFL.tzp);
  std::cout << "[veto] out=" << outPath << "\n";

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), FILE_SHARD);
  // Always: hasBranch() decides whether the truth record is there to read, and
  // with EXTENDED_BRANCHES on a missing branch would otherwise make the reader
  // iterate ZERO entries and exit 0 (Track_btagIp_* on every grid sample).
  recordAvailableBranches(chain);
  if (jvtOn && !hasBranch("Track_recoVtx_idx")) {
    std::cerr << "[veto] --jvt needs Track_recoVtx_idx, which this sample lacks; a silent "
                 "fallback would give every jet a JVT R_pT of 0 and remove them all\n";
    return 1;
  }
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  // Truth record for the VBF recoil proxy: raw ntuples only (the skims drop it).
  // Read lazily, so it costs I/O only for the events that get this far.
  const bool hasTruthMet = hasBranch("TruthPart_pdgId") && hasBranch("TruthPart_pt") &&
                           hasBranch("TruthPart_phi");
  std::unique_ptr<TTreeReaderArray<float>> tpPdg, tpPt, tpPhi;
  if (hasTruthMet) {
    tpPdg = std::make_unique<TTreeReaderArray<float>>(reader, "TruthPart_pdgId");
    tpPt  = std::make_unique<TTreeReaderArray<float>>(reader, "TruthPart_pt");
    tpPhi = std::make_unique<TTreeReaderArray<float>>(reader, "TruthPart_phi");
  }
  std::cout << "[veto] truth MET: " << (hasTruthMet ? "from the truth record"
                                                    : "UNAVAILABLE (no truth record; met_truth = -1)")
            << "\n";

  TFile out(outPath.c_str(), "RECREATE");
  if (out.IsZombie()) { std::cerr << "[veto] cannot open " << outPath << "\n"; return 1; }
  TTree tree("events", "VBS time veto ingredients, one row per selected event");
  Row R;
  tree.Branch("file_path", &R.file_path);
  tree.Branch("entry", &R.entry, "entry/L");
  tree.Branch("chain_entry", &R.chain_entry, "chain_entry/L");
  tree.Branch("file_idx", &R.file_idx, "file_idx/I");
  tree.Branch("weight", &R.weight, "weight/F");
  tree.Branch("met_truth", &R.met_truth, "met_truth/F");
  tree.Branch("z_pt", &R.z_pt, "z_pt/F");
  tree.Branch("z_mass", &R.z_mass, "z_mass/F");
  tree.Branch("n_jets", &R.n_jets, "n_jets/I");
  tree.Branch("n_jets_fwd_acc", &R.n_jets_fwd_acc, "n_jets_fwd_acc/I");
  tree.Branch("n_rm_jvt", &R.n_rm_jvt, "n_rm_jvt/I");
  tree.Branch("n_rm_fjvt", &R.n_rm_fjvt, "n_rm_fjvt/I");
  tree.Branch("mjj", &R.mjj, "mjj/D");
  tree.Branch("deta", &R.deta, "deta/D");
  auto legBranches = [&](const char* p, Leg& L) {
    auto n = [&](const char* f) { return std::string(p) + f; };
    auto br = [&](const char* f, void* addr, const char* type) {
      tree.Branch(n(f).c_str(), addr, (n(f) + "/" + type).c_str());
    };
    br("idx", &L.idx, "I");
    br("pt", &L.pt, "F");  br("eta", &L.eta, "F");  br("phi", &L.phi, "F");
    br("hs", &L.hs, "O");  br("pu", &L.pu, "O");
    br("jvt_rpt", &L.jvt_rpt, "F");  br("fjvt", &L.fjvt, "F");
    br("core_t", &L.core_t, "D");  br("core_sig", &L.core_sig, "D");
    br("core_n", &L.core_n, "I");  br("core_ncand", &L.core_ncand, "I");
    br("core_sumpt", &L.core_sumpt, "D");  br("core_hsfrac", &L.core_hsfrac, "D");
    br("all_t", &L.all_t, "D");  br("all_sig", &L.all_sig, "D");
    br("all_n", &L.all_n, "I");  br("all_hsfrac", &L.all_hsfrac, "D");
  };
  legBranches("a_", R.a);
  legBranches("b_", R.b);
  tree.Branch("jet_idx", &R.jet_idx);
  tree.Branch("jet_pt",  &R.jet_pt);
  tree.Branch("jet_eta", &R.jet_eta);
  tree.Branch("jet_phi", &R.jet_phi);
  tree.Branch("jet_hs",  &R.jet_hs);
  tree.Branch("jet_pu",  &R.jet_pu);
  tree.Branch("jet_t",   &R.jet_t);
  tree.Branch("jet_sig", &R.jet_sig);
  tree.Branch("jet_n",   &R.jet_n);
  for (int k = 0; k < N_T0; ++k) {
    const std::string s = T0_NAMES[k];
    tree.Branch(("t0_" + s).c_str(),   &R.t0[k],   ("t0_" + s + "/D").c_str());
    tree.Branch(("sig0_" + s).c_str(), &R.sig0[k], ("sig0_" + s + "/D").c_str());
    tree.Branch(("infl_" + s).c_str(), &R.infl[k], ("infl_" + s + "/D").c_str());
    tree.Branch(("ok_" + s).c_str(),   &R.ok[k],   ("ok_" + s + "/O").c_str());
  }

  Long64_t nRead = 0, nVtx = 0, nLep = 0, nJets = 0, nTagger = 0, nPair = 0, nSel = 0;
  double sumw = 0.0, sumw2 = 0.0, sumwSel = 0.0;
  Long64_t nBothCore = 0, nOneCore = 0;
  Calib calCore, calAll;

  while (reader.Next()) {
    if (maxEvents > 0 && nRead >= maxEvents) break;
    ++nRead;
    const double w = *branch.weight;
    sumw += w; sumw2 += w * w;
    if (nRead % 100000 == 0) std::cout << "[veto] read " << nRead << ", selected " << nSel << std::endl;

    // ── vbs_region_diag's preamble ───────────────────────────────────────────
    if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    ++nLep;
    if (!branch.passBasicCuts())       continue;
    ++nVtx;
    std::vector<int> passPtIdx;
    int nPt = 0, nPtEta = 0;
    branch.collectPtPassingJets(passPtIdx, nPt, nPtEta);
    if (nPt < MIN_PASSPT_JETS || nPtEta < MIN_PASSETA_JETS) continue;
    ++nJets;

    // ── Tagger before pairing, then the jet requirements again ──────────────
    std::vector<JetDisc> D;
    int nRmJvt = 0, nRmFjvt = 0;
    if (jvtOn) {
      D = computeJetDiscriminants(branch);
      std::vector<int> kept;
      int nPtWP = 0, nPtEtaWP = 0;
      for (int j : passPtIdx) {
        if (jvtRemoves(D[j], WP))  { ++nRmJvt;  continue; }
        if (fjvtRemoves(D[j], WP)) { ++nRmFjvt; continue; }
        kept.push_back(j);
        ++nPtWP;
        const double ae = std::abs((double)branch.topoJetEta[j]);
        if (ae > MIN_ABS_ETA_JET && ae < MAX_ABS_ETA_JET) ++nPtEtaWP;
      }
      if (nPtWP < MIN_PASSPT_JETS || nPtEtaWP < MIN_PASSETA_JETS) continue;
      passPtIdx.swap(kept);
    }
    ++nTagger;

    // ── The pair: max m_jj, opposite hemispheres (vbs_region_diag's wide_) ───
    const BranchPointerWrapper::VbsPair pair = branch.calcBestVbsPair(passPtIdx);
    if (!pair.valid()) continue;
    ++nPair;
    if (pair.mjj < VBS_JET_MJJ || pair.dEta < VBS_JET_D_ETA) continue;
    ++nSel;
    sumwSel += w;

    R = Row{};
    R.file_path   = chain.GetCurrentFile()->GetName();
    R.entry       = chain.GetTree()->GetReadEntry();
    R.chain_entry = chain.GetReadEntry();
    { const size_t p = R.file_path.rfind('_');
      const size_t d = (p == std::string::npos) ? std::string::npos : R.file_path.find('.', p);
      if (p != std::string::npos && d != std::string::npos)
        R.file_idx = std::atoi(R.file_path.substr(p + 1, d - p - 1).c_str()); }
    R.weight = (float)w;
    R.mjj = pair.mjj;
    R.deta = pair.dEta;
    R.n_jets = (int)passPtIdx.size();
    for (int j : passPtIdx) {
      const double ae = std::abs((double)branch.topoJetEta[j]);
      if (ae > MIN_ABS_ETA_JET && ae < MAX_ABS_ETA_JET) ++R.n_jets_fwd_acc;
    }
    R.n_rm_jvt = nRmJvt;
    R.n_rm_fjvt = nRmFjvt;

    // Recoil proxies.
    if (hasTruthMet) {
      double px = 0.0, py = 0.0;
      for (size_t p = 0; p < tpPdg->GetSize(); ++p) {
        const long a = std::labs(std::lround((*tpPdg)[p]));
        if (a != 12 && a != 14 && a != 16) continue;
        px += (*tpPt)[p] * std::cos((*tpPhi)[p]);
        py += (*tpPt)[p] * std::sin((*tpPhi)[p]);
      }
      R.met_truth = (float)std::hypot(px, py);
    }
    if (branch.trackLeptonID) {
      // vbs_region_diag's z_pt: good leptons by pT, first OS-SF pair, massless.
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
          R.z_pt = (float)(lx + ly).Pt();
          R.z_mass = (float)(lx + ly).M();
          done = true;
        }
    }

    // ── Every pairing candidate: kinematics, labels, core jet time ──────────
    std::unordered_map<int, JetTime> coreTime;
    for (int j : passPtIdx) {
      const JetTime c = jetTimeCore(branch, j);
      coreTime[j] = c;
      const float eta = branch.topoJetEta[j], phi = branch.topoJetPhi[j];
      R.jet_idx.push_back(j);
      R.jet_pt.push_back(branch.topoJetPt[j]);
      R.jet_eta.push_back(eta);
      R.jet_phi.push_back(phi);
      R.jet_hs.push_back(branch.isJetPaperHS(eta, phi) ? 1 : 0);
      R.jet_pu.push_back(branch.isJetPaperPU(eta, phi) ? 1 : 0);
      R.jet_t.push_back(c.t);
      R.jet_sig.push_back(c.sig);
      R.jet_n.push_back(c.n);
    }

    // ── Legs ────────────────────────────────────────────────────────────────
    const double tHS = branch.truthVtxTime[0];
    auto fillLeg = [&](int j, Leg& L) {
      L.idx = j;
      L.pt  = branch.topoJetPt[j];
      L.eta = branch.topoJetEta[j];
      L.phi = branch.topoJetPhi[j];
      L.hs  = branch.isJetPaperHS(L.eta, L.phi);
      L.pu  = branch.isJetPaperPU(L.eta, L.phi);
      if (jvtOn) { L.jvt_rpt = D[j].rpt; L.fjvt = D[j].fjvt; }
      const JetTime& c = coreTime.at(j);
      L.core_t = c.t;  L.core_sig = c.sig;  L.core_n = c.n;  L.core_ncand = c.nCand;
      L.core_sumpt = c.sumPt;  L.core_hsfrac = c.hsFrac;
      const JetTime g = jetTimeAllGhost(branch, j);
      L.all_t = g.t;  L.all_sig = g.sig;  L.all_n = g.n;  L.all_hsfrac = g.hsFrac;
      if (L.hs) { calCore.add(c, tHS); calAll.add(g, tHS); }
    };
    fillLeg(pair.idxI, R.a);
    fillLeg(pair.idxJ, R.b);
    const int nTimed = (R.a.core_n > 0) + (R.b.core_n > 0);
    if (nTimed == 2) ++nBothCore; else if (nTimed == 1) ++nOneCore;

    // ── t0: rpt_v5_hist's code and inputs, verbatim ──────────────────────────
    double t0v[N_T0] = {}, s0v[N_T0] = {};
    bool   okv[N_T0] = {};
    t0v[0] = branch.recoVtxTime[0];
    s0v[0] = branch.recoVtxTimeRes[0];
    okv[0] = (branch.recoVtxValid[0] == 1);
    {
      const std::vector<int> trk_z = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, T0_ASSOC_NSIGMA);
      const auto clusters = clusterTracksInTime(
          trk_z, &branch, DIST_CUT_CONE,
          /*useSmearedTimes=*/false, /*checkTimeValid=*/true, IDEAL_TRACK_RES,
          ClusteringMethod::ITERATIVE, /*useZ0=*/false,
          /*sortTracks=*/false, /*calcPurityFlag=*/true);
      if (!clusters.empty()) {
        auto bT = chooseCluster(clusters, Score::TRKPTZ);
        t0v[1] = bT.calculateTime(Score::TRKPTZ, &branch);
        s0v[1] = bT.sigmas[0];  okv[1] = true;
        auto bW = chooseCluster(clusters, Score::WAVES);
        t0v[2] = bW.calculateTime(Score::WAVES, &branch);        // in-jet refined
        s0v[2] = bW.sigmas[0];  okv[2] = true;
        auto bP = chooseCluster(clusters, Score::TRKPTZ_TZQ);
        t0v[3] = bP.calculateTime(Score::TRKPTZ_TZQ, &branch);   // guarded in-jet
        s0v[3] = bP.sigmas[0];  okv[3] = true;
      }
    }
    t0v[4] = tHS;  s0v[4] = 0.0;  okv[4] = true;
    const double inflv[N_T0] = {INFL.hgtd, INFL.trkptz, INFL.waves, INFL.tzp, 1.0};
    for (int k = 0; k < N_T0; ++k) {
      R.t0[k] = t0v[k];  R.sig0[k] = s0v[k];  R.ok[k] = okv[k];  R.infl[k] = inflv[k];
    }
    tree.Fill();
  }

  if (nRead == 0) {
    // A TTreeReader with a missing branch reports nothing: Next() just returns
    // false. Fail rather than write an empty tree with exit code 0.
    std::cerr << "[veto] ERROR: read 0 events -- check stderr for "
                 "TTreeReaderArrayBase::GetBranchAndLeaf errors\n";
    return 2;
  }

  // ── meta: one row per job, so hadd'd shards sum ────────────────────────────
  TTree meta("meta", "per-job cut flow and sum of weights over EVERY event read");
  int    jvtWp = (int)(wpSel - JVT_WPS);
  double cutMjj = VBS_JET_MJJ, cutDeta = VBS_JET_D_ETA, jdr = JET_TIME_DR, jcut = JET_TIME_DIST_CUT;
  bool   truthMet = hasTruthMet;
  std::string sample = SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME, ndir = cfg.ntupleDir;
  meta.Branch("n_read", &nRead, "n_read/L");
  meta.Branch("n_pass_lepton", &nLep, "n_pass_lepton/L");
  meta.Branch("n_pass_vtx", &nVtx, "n_pass_vtx/L");
  meta.Branch("n_pass_jets", &nJets, "n_pass_jets/L");
  meta.Branch("n_pass_tagger", &nTagger, "n_pass_tagger/L");
  meta.Branch("n_pass_pair", &nPair, "n_pass_pair/L");
  meta.Branch("n_selected", &nSel, "n_selected/L");
  meta.Branch("sumw_read", &sumw, "sumw_read/D");
  meta.Branch("sumw2_read", &sumw2, "sumw2_read/D");
  meta.Branch("sumw_selected", &sumwSel, "sumw_selected/D");
  meta.Branch("jvt_wp", &jvtWp, "jvt_wp/I");
  meta.Branch("vbs_mjj_cut", &cutMjj, "vbs_mjj_cut/D");
  meta.Branch("vbs_deta_cut", &cutDeta, "vbs_deta_cut/D");
  meta.Branch("jet_time_dr", &jdr, "jet_time_dr/D");
  meta.Branch("jet_time_dist_cut", &jcut, "jet_time_dist_cut/D");
  meta.Branch("has_truth_met", &truthMet, "has_truth_met/O");
  meta.Branch("sample", &sample);
  meta.Branch("ntuple_dir", &ndir);
  meta.Fill();
  out.cd();
  tree.Write();
  meta.Write();
  out.Close();

  std::printf("\n[veto] read %lld (sum w %.6g); Z->ll/OR %lld; vertex %lld; jets %lld; after --jvt=%s %lld; "
              "pair %lld; selected %lld (m_jj >= %g)\n",
              nRead, sumw, nLep, nVtx, nJets, WP.name, nTagger, nPair, nSel, VBS_JET_MJJ);
  std::printf("[veto] core jet time on both legs %lld (%.1f%%), one leg %lld (%.1f%%), neither %lld\n",
              nBothCore, nSel ? 100.0 * nBothCore / nSel : 0.0,
              nOneCore, nSel ? 100.0 * nOneCore / nSel : 0.0, nSel - nBothCore - nOneCore);
  std::printf("\n=== JET-TIME CALIBRATION (paper-HS legs with a time; dt = t_jet - t_HS^truth) ===\n");
  std::printf("  %-6s %8s %10s %10s %10s %7s %9s %9s %9s\n", "est.", "n", "mean dt", "core rms",
              "quoted", "ratio", "MAD pull", "|pull|<3", "|dt|>150");
  calCore.print("core");
  calAll.print("all");
  std::printf("  ratio   : rpt_v5's t0 definition (RMS about the mean in |dt| < 150 ps / RMS quoted sigma)\n");
  std::printf("  MAD pull: 1.4826 x median |dt/sigma| -- the Gaussian-core width of the per-jet pull\n");
  std::printf("  (the plotting script measures and applies the calibration; nothing here uses it)\n");
  std::cout << "[veto] wrote " << outPath << "\n";
  return 0;
}
