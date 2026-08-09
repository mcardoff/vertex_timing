// ---------------------------------------------------------------------------
// export_training_data.cxx
//
//   Exports per-cluster and per-track training data for the cluster-SELECTION
//   ML study (see share/docs / the WAVeS selection analysis).
//
//   Motivation
//   ----------
//   Hand-designed selection scores (TRKPTZ, WAVeS and its LOJO/JETCAP/KERNEL
//   variants) plateau far below what the clustering makes available: on Z+jets a
//   cluster within 60 ps of truth exists in ~87% of events, but the best score
//   selects it only ~62% of the time. The gap is a SELECTION problem, so this
//   exporter dumps the raw ingredients for a learned selector.
//
//   Output: one ROOT file per sample containing two TTrees
//     clusters : one row per cluster  (fixed-width features -> GBDT ranking)
//     tracks   : one row per track    (enables a per-track P(HS) model, whose
//                aggregate Sum pT*P(HS) is the learned form of what WAVeS
//                hard-codes as pT_jet/dR)
//
//   ROOT rather than CSV/parquet: the track tree is ~9M rows across samples
//   (~2.7 GB as CSV, ~350 MB compressed); ROOT is already linked here, whereas
//   parquet would need Apache Arrow C++. uproot reads either into pandas with
//   per-column selection, so nothing is lost on the analysis side.
//
//   IMPORTANT - truth handling
//   --------------------------
//   Every truth-derived column is prefixed `truth_` and must be excluded from
//   the model inputs. `delta_t` (cluster time - truth HS time) is the TARGET;
//   derive labels from it offline (argmin |delta_t| within event = best cluster;
//   |delta_t| < 60 ps = "good"). Event SELECTION here uses reco quantities only,
//   so the training population is reproducible at inference time.
//
//   Usage:  ./export_training_data [--sample=vbf|zjets|dijet] [--max-events=N]
// ---------------------------------------------------------------------------

#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TVector2.h>

#include <boost/filesystem.hpp>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "sample_config.h"

// ---------------------------------------------------------------------------
// TODO (for the no-clustering transformer, python/train_transformer.py)
//
// 1. PER-TRACK HS-JET LABEL.  Add a truth_ column giving the identity of the jet
//    a track is associated with:
//       truth_nearest_fwdjet_is_hs  : 1 if the nearest forward jet is matched to
//                                     a truth HS jet, else 0.
//    Derivable from branches that ALREADY survive the ntuple slim: a reco jet is
//    HS when AntiKt4EMTopoJets_truthHSJet_idx (bound as topoJetTruthHSIdx) is
//    non-empty, and PU otherwise -- the same definition rpt_v5's paperIsPU uses.
//    Note this needs the truth_ prefix so the leakage guard holds it out of the
//    feature set: it is supervision, never an input.
//
//    (TruthITPUJet_*/TruthOOTPUJet_* are NOT needed for this. Only the HS jets
//    are labelled positively; everything unmatched is PU by construction.)
//
// 2. JETS AS TOKENS.  Add a third TTree, one row per forward jet:
//       sample_id, event_num, jet_idx, pt, eta, phi, is_forward,
//       n_ghost_tracks, sumpt_ghost, truth_is_hs
//    Today jets reach the model only through per-track association
//    (dr_nearest_fwdjet, pt_nearest_fwdjet, is_ghost_of_nearest) and event-level
//    scalars (lead_jet_pt, n_forward_jets). A transformer can attend over jets
//    directly if they are their own token sequence, which the per-track summary
//    cannot express.
//
// Both are additive: existing consumers are unaffected, and train_transformer.py
// picks them up automatically when present (it treats them as optional).
// ---------------------------------------------------------------------------


using namespace MyUtl;

namespace {

constexpr float NaNf = std::numeric_limits<float>::quiet_NaN();

// Nearest qualifying forward jet to a track: returns {dR, jetPt, jetIndex}.
// "Qualifying" mirrors the WAVeS definition (pT > MIN_JET_PT, |eta| inside the
// HGTD forward band, not lepton-overlap removed) so the exported jet features
// are the same ingredients WAVeS combines by hand -- unbundled, so a model can
// learn their weighting instead of inheriting the 1/dR form.
struct NearJet { double dr = 1e9; double pt = 0.0; int idx = -1; };

NearJet nearestForwardJet(int trk, BranchPointerWrapper* b) {
  NearJet nj;
  const double te = b->trackEta[trk], tp = b->trackPhi[trk];
  for (int j = 0; j < (int)b->topoJetEta.GetSize(); ++j) {
    if (b->isJetRemoved(j)) continue;
    const double je = b->topoJetEta[j], jp = b->topoJetPt[j];
    if (jp < MIN_JET_PT) continue;
    if (std::abs(je) < MIN_ABS_ETA_JET || std::abs(je) > MAX_ABS_ETA_JET) continue;
    const double deta = je - te;
    const double dphi = TVector2::Phi_mpi_pi(b->topoJetPhi[j] - tp);
    const double dr   = std::hypot(deta, dphi);
    if (dr < nj.dr) { nj.dr = dr; nj.pt = jp; nj.idx = j; }
  }
  return nj;
}

float median(std::vector<float> v) {
  if (v.empty()) return NaNf;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n % 2) ? v[n / 2] : 0.5f * (v[n / 2 - 1] + v[n / 2]);
}

float rms(const std::vector<float>& v) {
  if (v.size() < 2) return 0.f;
  const float m = std::accumulate(v.begin(), v.end(), 0.f) / v.size();
  float s = 0.f;
  for (float x : v) s += (x - m) * (x - m);
  return std::sqrt(s / v.size());
}

// One row of the cluster tree. Plain floats so every branch registers the same
// way; integer-valued quantities are stored as floats deliberately (GBDT does
// not care, and it keeps the schema uniform for uproot -> pandas).
struct ClusterRow {
  // --- bookkeeping -------------------------------------------------------
  float event_num, cluster_idx, sample_id, weight;
  // --- Tier 0: original cluster features ---------------------------------
  float cluster_time, delta_z, delta_z_resunits, cluster_z_sigma;
  float cluster_d0, cluster_d0_sigma, cluster_qOverP, cluster_qOverP_sigma;
  float sumpt, cluster_time_sigma, n_tracks;
  // --- Tier A: pT spectrum / hardness ------------------------------------
  float sumpt2, maxpt, pt_2nd, pt_3rd, lead_pt_frac, meanpt, medianpt;
  float n_pt_gt2, n_pt_gt5;
  // --- Tier B: timing coherence ------------------------------------------
  float time_chi2_ndf, max_abs_tpull, n_tpull_gt2, time_spread;
  float min_timeRes, median_timeRes, max_timeRes;
  float n_valid_time, frac_valid_time, sumpt_valid_time;
  // --- Tier C: z compatibility with the reco PV --------------------------
  float z_chi2_ndf, max_abs_zpull, n_zpull_gt3, z_spread, median_z0_pull;
  // --- Tier D: HGTD detector quality -------------------------------------
  // mean_nhgtd_primary is written under the TRUTH prefix: Track_nHGTDPrimaryHits
  // is built in Athena from per-hit m_isprime ("keep which of the hits associated
  // in reco were primary hits (truth info!)"), so it is MC-only and unavailable
  // in data. Kept for diagnostics; the truth_ prefix keeps it out of the features.
  float mean_nhgtd, min_nhgtd, max_nhgtd, mean_nhgtd_primary;
  float frac_nhgtd_ge2, sumpt_w_nhgtd;
  // --- Tier E: track quality / kinematics --------------------------------
  float mean_quality, min_quality, mean_abs_eta, max_abs_eta, min_abs_eta;
  float eta_spread, mean_d0_signif, max_abs_d0_signif;
  // --- Tier F: jet & lepton association ----------------------------------
  float n_tracks_in_fwdjet, sumpt_in_fwdjet, frac_pt_in_fwdjet;
  float min_dr_to_fwdjet, pt_of_nearest_fwdjet, max_matched_jetpt;
  float n_ghost_tracks, sumpt_ghost;
  float trkptz_score, waves_score;
  float n_lepton_tracks, sumpt_lepton, dz_to_lepton_z0, dz_to_lepton_signif;
  // --- Tier G: event context (relative; filled in a second pass) ----------
  float n_clusters, sumpt_rank, sumpt_frac_of_event, sumpt_ratio_to_max;
  float trkptz_rank, trkptz_ratio_to_max, is_max_sumpt, is_max_trkptz;
  float dt_to_nearest_cluster, dz_to_nearest_cluster, n_clusters_within_60ps;
  // --- Tier H: event / topology globals ----------------------------------
  float n_forward_jets, lead_jet_pt, lead_jet_abseta;
  float sublead_jet_pt, sublead_jet_abseta;
  float n_fwd_tracks_reco, event_sumpt, n_reco_vertices, reco_vtx_z;
  // --- Tier I: HGTD external prior ---------------------------------------
  float hgtd_time, hgtd_valid, hgtd_time_res, dt_cluster_to_hgtd;
  // --- TARGET + truth-only diagnostics (never model inputs) ---------------
  float delta_t;                 // TARGET: cluster time - truth HS time
  float truth_purity;            // HS pT fraction of the cluster
  float truth_n_hs_tracks;       // HS tracks in this cluster
  float truth_hs_frac_tracks;    // HS track fraction
};

struct TrackRow {
  float event_num, cluster_idx, track_idx, sample_id;
  float pt, eta, phi, theta, z0, d0, qOverP;
  float sigma_z0, sigma_d0, sigma_qOverP;
  float time, timeRes, time_valid;
  float quality, nhgtd_hits, nhgtd_primary;  // nhgtd_primary -> truth_ branch (MC-only)
  float z0_pull_pv, t_pull_cluster;
  float dr_nearest_fwdjet, pt_nearest_fwdjet, is_ghost_of_nearest;
  float is_lepton;
  float cluster_time, cluster_delta_z;   // context, so track rows stand alone
  float truth_is_hs;                     // LABEL for the per-track P(HS) model
};

}  // namespace

// ---------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

  const auto cfg = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);  // --vbs-deta=<x>; sets SELECTION_TAG
  MyUtl::SAMPLE_NAME  = cfg.sampleName;
  MyUtl::OUTPUT_DIR   = cfg.outputDir;
  MyUtl::ENERGY_LABEL = cfg.energyLabel;
  MyUtl::OVERLAP_REMOVAL = cfg.overlapRemoval;
  const Long64_t maxEvents = MyUtl::resolveMaxEvents(argc, argv);

  // Numeric sample tag travels with every row so concatenated files still group
  // correctly: the ranking group key is (sample_id, event_num), since event_num
  // restarts at 0 for each sample's run.
  const float sampleId =
      cfg.sampleName == "vbf"   ? 0.f :
      cfg.sampleName == "zjets" ? 1.f :
      cfg.sampleName == "dijet" ? 2.f : 3.f;   // 3 = local/default

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  // Output directory must be created before TFile::Open -- ROOT will not make it,
  // and on condor the job's scratch dir starts empty (mirrors clustering_hist).
  boost::filesystem::create_directories(MyUtl::OUTPUT_DIR);
  if (MyUtl::SAMPLE_NAME.empty())
    boost::filesystem::create_directories(MyUtl::OUTPUT_DIR + "/hists");

  // Pass the UNPREFIXED base name: histFilePath() adds "<sample>_" itself when
  // a --sample was given, so "<sample>_training.root" here would double-prefix
  // into e.g. zjets/zjets_zjets_training.root and break the condor output remap.
  const std::string outPath = MyUtl::histFilePath("training.root");
  std::unique_ptr<TFile> out(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "ERROR: cannot open " << outPath << " for writing\n";
    return 1;
  }

  auto* clusterTree = new TTree("clusters", "one row per cluster");
  auto* trackTree   = new TTree("tracks",   "one row per track");
  ClusterRow C{};
  TrackRow   T{};
  auto BC = [&](const char* n, float* p) { clusterTree->Branch(n, p, Form("%s/F", n)); };
  auto BT = [&](const char* n, float* p) { trackTree  ->Branch(n, p, Form("%s/F", n)); };

  BC("event_num",&C.event_num); BC("cluster_idx",&C.cluster_idx);
  BC("sample_id",&C.sample_id); BC("weight",&C.weight);
  BC("cluster_time",&C.cluster_time); BC("delta_z",&C.delta_z);
  BC("delta_z_resunits",&C.delta_z_resunits); BC("cluster_z_sigma",&C.cluster_z_sigma);
  BC("cluster_d0",&C.cluster_d0); BC("cluster_d0_sigma",&C.cluster_d0_sigma);
  BC("cluster_qOverP",&C.cluster_qOverP); BC("cluster_qOverP_sigma",&C.cluster_qOverP_sigma);
  BC("sumpt",&C.sumpt); BC("cluster_time_sigma",&C.cluster_time_sigma);
  BC("n_tracks",&C.n_tracks);
  BC("sumpt2",&C.sumpt2); BC("maxpt",&C.maxpt); BC("pt_2nd",&C.pt_2nd);
  BC("pt_3rd",&C.pt_3rd); BC("lead_pt_frac",&C.lead_pt_frac);
  BC("meanpt",&C.meanpt); BC("medianpt",&C.medianpt);
  BC("n_pt_gt2",&C.n_pt_gt2); BC("n_pt_gt5",&C.n_pt_gt5);
  BC("time_chi2_ndf",&C.time_chi2_ndf); BC("max_abs_tpull",&C.max_abs_tpull);
  BC("n_tpull_gt2",&C.n_tpull_gt2); BC("time_spread",&C.time_spread);
  BC("min_timeRes",&C.min_timeRes); BC("median_timeRes",&C.median_timeRes);
  BC("max_timeRes",&C.max_timeRes); BC("n_valid_time",&C.n_valid_time);
  BC("frac_valid_time",&C.frac_valid_time); BC("sumpt_valid_time",&C.sumpt_valid_time);
  BC("z_chi2_ndf",&C.z_chi2_ndf); BC("max_abs_zpull",&C.max_abs_zpull);
  BC("n_zpull_gt3",&C.n_zpull_gt3); BC("z_spread",&C.z_spread);
  BC("median_z0_pull",&C.median_z0_pull);
  BC("mean_nhgtd",&C.mean_nhgtd); BC("min_nhgtd",&C.min_nhgtd);
  BC("max_nhgtd",&C.max_nhgtd); BC("truth_mean_nhgtd_primary",&C.mean_nhgtd_primary);
  BC("frac_nhgtd_ge2",&C.frac_nhgtd_ge2); BC("sumpt_w_nhgtd",&C.sumpt_w_nhgtd);
  BC("mean_quality",&C.mean_quality); BC("min_quality",&C.min_quality);
  BC("mean_abs_eta",&C.mean_abs_eta); BC("max_abs_eta",&C.max_abs_eta);
  BC("min_abs_eta",&C.min_abs_eta); BC("eta_spread",&C.eta_spread);
  BC("mean_d0_signif",&C.mean_d0_signif); BC("max_abs_d0_signif",&C.max_abs_d0_signif);
  BC("n_tracks_in_fwdjet",&C.n_tracks_in_fwdjet); BC("sumpt_in_fwdjet",&C.sumpt_in_fwdjet);
  BC("frac_pt_in_fwdjet",&C.frac_pt_in_fwdjet); BC("min_dr_to_fwdjet",&C.min_dr_to_fwdjet);
  BC("pt_of_nearest_fwdjet",&C.pt_of_nearest_fwdjet);
  BC("max_matched_jetpt",&C.max_matched_jetpt);
  BC("n_ghost_tracks",&C.n_ghost_tracks); BC("sumpt_ghost",&C.sumpt_ghost);
  BC("trkptz_score",&C.trkptz_score); BC("waves_score",&C.waves_score);
  BC("n_lepton_tracks",&C.n_lepton_tracks); BC("sumpt_lepton",&C.sumpt_lepton);
  BC("dz_to_lepton_z0",&C.dz_to_lepton_z0); BC("dz_to_lepton_signif",&C.dz_to_lepton_signif);
  BC("n_clusters",&C.n_clusters); BC("sumpt_rank",&C.sumpt_rank);
  BC("sumpt_frac_of_event",&C.sumpt_frac_of_event);
  BC("sumpt_ratio_to_max",&C.sumpt_ratio_to_max);
  BC("trkptz_rank",&C.trkptz_rank); BC("trkptz_ratio_to_max",&C.trkptz_ratio_to_max);
  BC("is_max_sumpt",&C.is_max_sumpt); BC("is_max_trkptz",&C.is_max_trkptz);
  BC("dt_to_nearest_cluster",&C.dt_to_nearest_cluster);
  BC("dz_to_nearest_cluster",&C.dz_to_nearest_cluster);
  BC("n_clusters_within_60ps",&C.n_clusters_within_60ps);
  BC("n_forward_jets",&C.n_forward_jets); BC("lead_jet_pt",&C.lead_jet_pt);
  BC("lead_jet_abseta",&C.lead_jet_abseta); BC("sublead_jet_pt",&C.sublead_jet_pt);
  BC("sublead_jet_abseta",&C.sublead_jet_abseta);
  BC("n_fwd_tracks_reco",&C.n_fwd_tracks_reco); BC("event_sumpt",&C.event_sumpt);
  BC("n_reco_vertices",&C.n_reco_vertices); BC("reco_vtx_z",&C.reco_vtx_z);
  BC("hgtd_time",&C.hgtd_time); BC("hgtd_valid",&C.hgtd_valid);
  BC("hgtd_time_res",&C.hgtd_time_res); BC("dt_cluster_to_hgtd",&C.dt_cluster_to_hgtd);
  BC("delta_t",&C.delta_t);
  BC("truth_purity",&C.truth_purity); BC("truth_n_hs_tracks",&C.truth_n_hs_tracks);
  BC("truth_hs_frac_tracks",&C.truth_hs_frac_tracks);

  BT("event_num",&T.event_num); BT("cluster_idx",&T.cluster_idx);
  BT("track_idx",&T.track_idx); BT("sample_id",&T.sample_id);
  BT("pt",&T.pt); BT("eta",&T.eta); BT("phi",&T.phi); BT("theta",&T.theta);
  BT("z0",&T.z0); BT("d0",&T.d0); BT("qOverP",&T.qOverP);
  BT("sigma_z0",&T.sigma_z0); BT("sigma_d0",&T.sigma_d0); BT("sigma_qOverP",&T.sigma_qOverP);
  BT("time",&T.time); BT("timeRes",&T.timeRes); BT("time_valid",&T.time_valid);
  BT("quality",&T.quality); BT("nhgtd_hits",&T.nhgtd_hits);
  BT("truth_nhgtd_primary",&T.nhgtd_primary);
  BT("z0_pull_pv",&T.z0_pull_pv); BT("t_pull_cluster",&T.t_pull_cluster);
  BT("dr_nearest_fwdjet",&T.dr_nearest_fwdjet);
  BT("pt_nearest_fwdjet",&T.pt_nearest_fwdjet);
  BT("is_ghost_of_nearest",&T.is_ghost_of_nearest);
  BT("is_lepton",&T.is_lepton);
  BT("cluster_time",&T.cluster_time); BT("cluster_delta_z",&T.cluster_delta_z);
  BT("truth_is_hs",&T.truth_is_hs);

  const Long64_t nTotal = chain.GetEntries();
  const Long64_t nLoop  = (maxEvents > 0 && maxEvents < nTotal) ? maxEvents : nTotal;
  Long64_t nPassEvent = 0, nClusterRows = 0, nTrackRows = 0;
  std::cout << "Sample '" << (cfg.sampleName.empty() ? "local" : cfg.sampleName)
            << "' -- exporting over " << nLoop << " / " << nTotal << " events\n";

  while (reader.Next()) {
    const Long64_t evtNum = chain.GetReadEntry();
    if (maxEvents > 0 && evtNum >= maxEvents) break;
    if ((evtNum + 1) % 20000 == 0)
      std::cout << "Progress: " << evtNum + 1 << " / " << nLoop << "\r" << std::flush;

    // ---- Event selection: RECO ONLY --------------------------------------
    // Mirrors processEventData's reco cuts so the training population is
    // reproducible at inference. The previous truth-based `nFTrackHS < 3` cut
    // was removed deliberately: it cannot be applied in production and biased
    // the sample toward events the algorithm already handles well.
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    if (!branch.passBasicCuts())       continue;
    if (!branch.passJetPtCut())        continue;

    std::vector<int> tracks =
        getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, MAX_NSIGMA);
    if (tracks.empty()) continue;

    std::vector<Cluster> clusters = clusterTracksInTime(
        tracks, &branch, DIST_CUT_CONE,
        /*useSmearedTimes*/ false, /*checkTimeValid*/ true,
        /*smearRes*/ IDEAL_TRACK_RES, ClusteringMethod::ITERATIVE,
        /*usez0*/ false, /*sortTracks*/ false, /*calcPurityFlag*/ true);
    if (clusters.empty()) continue;
    ++nPassEvent;

    // ---- Event-level quantities (Tier H/I), computed once -----------------
    const double vtxZ = branch.recoVtxZ[0];
    int nFTrack = 0, nFTrackHS = 0, nFTrackPU = 0;
    branch.countForwardTracks(nFTrack, nFTrackHS, nFTrackPU, tracks, true);

    double leadPt = 0, leadEta = 0, subPt = 0, subEta = 0;
    int nFwdJets = 0;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j)) continue;
      const double jp = branch.topoJetPt[j], je = std::abs(branch.topoJetEta[j]);
      if (jp > leadPt) { subPt = leadPt; subEta = leadEta; leadPt = jp; leadEta = je; }
      else if (jp > subPt) { subPt = jp; subEta = je; }
      if (jp >= MIN_JET_PT && je >= MIN_ABS_ETA_JET && je <= MAX_ABS_ETA_JET) ++nFwdJets;
    }

    double evtSumPt = 0;
    for (int t : tracks) evtSumPt += branch.trackPt[t];

    // Leading selected lepton's z0: an HS anchor for Z->ll. Kept as a distance
    // rather than membership because Z leptons (~40+ GeV) are cut by
    // MAX_TRACK_PT and so usually never enter a cluster at all.
    double lepZ0 = NaNf, lepPtMax = 0;
    if (branch.trackLeptonID) {
      for (int t = 0; t < (int)branch.trackLeptonID->GetSize(); ++t)
        if (branch.isGoodLepton(t) && branch.trackPt[t] > lepPtMax) {
          lepPtMax = branch.trackPt[t]; lepZ0 = branch.trackZ0[t];
        }
    }

    const bool  hgtdValid = branch.recoVtxValid[0] != 0;
    const float hgtdTime  = hgtdValid ? (float)branch.recoVtxTime[0] : NaNf;
    const float hgtdRes   = hgtdValid ? (float)branch.recoVtxTimeRes[0] : NaNf;

    // ---- Per-cluster features --------------------------------------------
    std::vector<ClusterRow> rows;
    rows.reserve(clusters.size());
    for (size_t ci = 0; ci < clusters.size(); ++ci) {
      const Cluster& cl = clusters[ci];
      if (cl.values.empty() || cl.trackIndices.empty()) continue;
      // NOTE: no nConstituents cut. Dropping small clusters would remove the
      // truth-closest candidate from some events and silently corrupt the
      // ranking label; the candidate set must match what the selector sees.

      ClusterRow R{};
      R.event_num = (float)evtNum; R.cluster_idx = (float)ci;
      R.sample_id = sampleId;      R.weight = *branch.weight;

      const double tClus = cl.values.at(0);
      R.cluster_time       = (float)tClus;
      R.cluster_time_sigma = (float)cl.sigmas.at(0);
      R.n_tracks           = (float)cl.trackIndices.size();

      double znum=0, zden=0, dnum=0, dden=0, qnum=0, qden=0, sumpt=0, sumpt2=0;
      std::vector<float> pts, tres, z0pulls, etas;
      double tchi2=0, zchi2=0, maxTPull=0, maxZPull=0, maxD0Sig=0, sumD0Sig=0;
      int nTPullGt2=0, nZPullGt3=0, nValidT=0, nPtGt2=0, nPtGt5=0, nNhgtdGe2=0;
      double sumptValidT=0, sumNhgtd=0, sumNhgtdPrim=0, sumptWNhgtd=0;
      double minNhgtd=1e9, maxNhgtd=-1e9, sumQual=0, minQual=1e9;
      double nInJet=0, sumptInJet=0, minDrJet=1e9, ptNearJet=NaNf, maxMatchedJetPt=0;
      double nGhost=0, sumptGhost=0, nLep=0, sumptLep=0;
      int nHsTracks=0;

      for (size_t k = 0; k < cl.trackIndices.size(); ++k) {
        const int trk = cl.trackIndices[k];
        const double pt = branch.trackPt[trk];
        const double sz = std::sqrt(branch.trackVarZ0[trk]);
        const double sd = std::sqrt(branch.trackVarD0[trk]);
        const double zp = (branch.trackZ0[trk] - vtxZ) / sz;
        const double d0s= std::abs(branch.trackD0[trk]) / sd;

        znum += branch.trackZ0[trk]/branch.trackVarZ0[trk]; zden += 1.0/branch.trackVarZ0[trk];
        dnum += branch.trackD0[trk]/branch.trackVarD0[trk]; dden += 1.0/branch.trackVarD0[trk];
        qnum += branch.trackQP[trk]/branch.trackVarQp[trk]; qden += 1.0/branch.trackVarQp[trk];

        sumpt += pt; sumpt2 += pt*pt; pts.push_back((float)pt);
        if (pt > 2.0) ++nPtGt2;
        if (pt > 5.0) ++nPtGt5;

        zchi2 += zp*zp; maxZPull = std::max(maxZPull, std::abs(zp));
        if (std::abs(zp) > 3.0) ++nZPullGt3;
        z0pulls.push_back((float)zp);

        sumD0Sig += d0s; maxD0Sig = std::max(maxD0Sig, d0s);
        etas.push_back((float)std::abs(branch.trackEta[trk]));

        const double sigT = branch.trackTimeRes[trk];
        const bool   okT  = branch.trackTimeValid[trk] != 0 && sigT > 0.0;
        if (okT) {
          const double tp = (cl.allTimes[k] - tClus) / sigT;
          tchi2 += tp*tp; maxTPull = std::max(maxTPull, std::abs(tp));
          if (std::abs(tp) > 2.0) ++nTPullGt2;
          tres.push_back((float)sigT);
          ++nValidT; sumptValidT += pt;
        }

        const double nh = branch.trackHgtdHits[trk];
        sumNhgtd += nh; sumNhgtdPrim += branch.trackPrimHits[trk];
        minNhgtd = std::min(minNhgtd, nh); maxNhgtd = std::max(maxNhgtd, nh);
        if (nh >= 2) ++nNhgtdGe2;
        sumptWNhgtd += pt * nh;

        const double q = branch.trackQuality[trk];
        sumQual += q; minQual = std::min(minQual, q);

        const NearJet nj = nearestForwardJet(trk, &branch);
        if (nj.idx >= 0) {
          if (nj.dr < 0.4) { nInJet += 1; sumptInJet += pt; }
          if (nj.dr < minDrJet) { minDrJet = nj.dr; ptNearJet = (float)nj.pt; }
          maxMatchedJetPt = std::max(maxMatchedJetPt, nj.pt);
          const auto& gh = branch.topoJetGhostTrackIdx[nj.idx];
          if (std::find(gh.begin(), gh.end(), trk) != gh.end()) {
            nGhost += 1; sumptGhost += pt;
          }
        }
        if (branch.trackLeptonID && branch.isGoodLepton(trk)) { nLep += 1; sumptLep += pt; }
        if (branch.trackToTruthvtx[trk] == 0) ++nHsTracks;
      }

      const int n = (int)cl.trackIndices.size();
      const double clusZ = znum/zden, clusZSig = 1.0/std::sqrt(zden);
      R.delta_z = (float)(clusZ - vtxZ);
      R.cluster_z_sigma = (float)clusZSig;
      R.delta_z_resunits = (float)(R.delta_z / clusZSig);
      R.cluster_d0 = (float)(dnum/dden);       R.cluster_d0_sigma = (float)(1.0/std::sqrt(dden));
      R.cluster_qOverP = (float)(qnum/qden);   R.cluster_qOverP_sigma = (float)(1.0/std::sqrt(qden));
      R.sumpt = (float)sumpt;

      std::sort(pts.begin(), pts.end(), std::greater<float>());
      R.sumpt2 = (float)sumpt2;
      R.maxpt  = pts.size() > 0 ? pts[0] : NaNf;
      R.pt_2nd = pts.size() > 1 ? pts[1] : 0.f;
      R.pt_3rd = pts.size() > 2 ? pts[2] : 0.f;
      R.lead_pt_frac = sumpt > 0 ? (float)(pts[0]/sumpt) : NaNf;
      R.meanpt = (float)(sumpt/n);  R.medianpt = median(pts);
      R.n_pt_gt2 = (float)nPtGt2;   R.n_pt_gt5 = (float)nPtGt5;

      // chi2/ndf ~ 1 for a genuine single vertex; inflated when the cluster
      // merges two vertices or absorbs mistimed PU -- the discriminant that
      // cluster_time_sigma alone cannot express.
      R.time_chi2_ndf   = nValidT > 1 ? (float)(tchi2/(nValidT-1)) : NaNf;
      R.max_abs_tpull   = nValidT > 0 ? (float)maxTPull : NaNf;
      R.n_tpull_gt2     = (float)nTPullGt2;
      R.time_spread     = (float)const_cast<Cluster&>(cl).timeSpread();
      R.min_timeRes     = tres.empty() ? NaNf : *std::min_element(tres.begin(), tres.end());
      R.max_timeRes     = tres.empty() ? NaNf : *std::max_element(tres.begin(), tres.end());
      R.median_timeRes  = median(tres);
      R.n_valid_time    = (float)nValidT;
      R.frac_valid_time = (float)nValidT/n;
      R.sumpt_valid_time= (float)sumptValidT;

      R.z_chi2_ndf     = (float)(zchi2/n);
      R.max_abs_zpull  = (float)maxZPull;
      R.n_zpull_gt3    = (float)nZPullGt3;
      R.z_spread       = rms(z0pulls) * 0.f + (float)const_cast<Cluster&>(cl).zSpread(&branch);
      R.median_z0_pull = median(z0pulls);

      R.mean_nhgtd = (float)(sumNhgtd/n);  R.min_nhgtd = (float)minNhgtd;
      R.max_nhgtd  = (float)maxNhgtd;      R.mean_nhgtd_primary = (float)(sumNhgtdPrim/n);
      R.frac_nhgtd_ge2 = (float)nNhgtdGe2/n;
      R.sumpt_w_nhgtd  = sumpt > 0 ? (float)(sumptWNhgtd/sumpt) : NaNf;

      R.mean_quality = (float)(sumQual/n);  R.min_quality = (float)minQual;
      R.mean_abs_eta = (float)(std::accumulate(etas.begin(),etas.end(),0.f)/n);
      R.max_abs_eta  = *std::max_element(etas.begin(), etas.end());
      R.min_abs_eta  = *std::min_element(etas.begin(), etas.end());
      R.eta_spread   = rms(etas);
      R.mean_d0_signif = (float)(sumD0Sig/n); R.max_abs_d0_signif = (float)maxD0Sig;

      R.n_tracks_in_fwdjet = (float)nInJet;  R.sumpt_in_fwdjet = (float)sumptInJet;
      R.frac_pt_in_fwdjet  = sumpt > 0 ? (float)(sumptInJet/sumpt) : NaNf;
      R.min_dr_to_fwdjet   = minDrJet < 1e8 ? (float)minDrJet : NaNf;
      R.pt_of_nearest_fwdjet = (float)ptNearJet;
      R.max_matched_jetpt  = (float)maxMatchedJetPt;
      R.n_ghost_tracks = (float)nGhost;      R.sumpt_ghost = (float)sumptGhost;
      R.trkptz_score = cl.scores.count(Score::TRKPTZ.id) ? (float)cl.scores.at(Score::TRKPTZ.id) : NaNf;
      R.waves_score  = cl.scores.count(Score::WAVES.id)  ? (float)cl.scores.at(Score::WAVES.id)  : NaNf;
      R.n_lepton_tracks = (float)nLep;       R.sumpt_lepton = (float)sumptLep;
      R.dz_to_lepton_z0 = std::isnan(lepZ0) ? NaNf : (float)std::abs(clusZ - lepZ0);
      R.dz_to_lepton_signif = std::isnan(lepZ0) ? NaNf : (float)(std::abs(clusZ - lepZ0)/clusZSig);

      R.n_forward_jets = (float)nFwdJets;
      R.lead_jet_pt = (float)leadPt;       R.lead_jet_abseta = (float)leadEta;
      R.sublead_jet_pt = (float)subPt;     R.sublead_jet_abseta = (float)subEta;
      R.n_fwd_tracks_reco = (float)nFTrack; R.event_sumpt = (float)evtSumPt;
      R.n_reco_vertices = (float)branch.recoVtxZ.GetSize();
      R.reco_vtx_z = (float)vtxZ;

      R.hgtd_time = hgtdTime; R.hgtd_valid = hgtdValid ? 1.f : 0.f;
      R.hgtd_time_res = hgtdRes;
      R.dt_cluster_to_hgtd = hgtdValid ? (float)std::abs(tClus - branch.recoVtxTime[0]) : NaNf;

      R.delta_t = (float)(tClus - branch.truthVtxTime[0]);
      R.truth_purity = (float)cl.purity;
      R.truth_n_hs_tracks = (float)nHsTracks;
      R.truth_hs_frac_tracks = (float)nHsTracks/n;

      rows.push_back(R);
    }
    if (rows.empty()) continue;

    // ---- Tier G: within-event context (second pass) -----------------------
    // These make a per-cluster classifier behave like a ranker, which matters
    // because the available GBDT (sklearn HistGradientBoosting) has no listwise
    // objective -- the relative information has to live in the features.
    const int nC = (int)rows.size();
    double maxSumpt = 0, sumSumpt = 0, maxTrkptz = 0;
    for (const auto& r : rows) {
      maxSumpt = std::max(maxSumpt, (double)r.sumpt); sumSumpt += r.sumpt;
      if (!std::isnan(r.trkptz_score)) maxTrkptz = std::max(maxTrkptz, (double)r.trkptz_score);
    }
    for (int i = 0; i < nC; ++i) {
      auto& r = rows[i];
      int rankS = 1, rankT = 1;
      double nearestDt = 1e9, nearestDz = 1e9; int nWithin60 = 0;
      for (int j = 0; j < nC; ++j) {
        if (i == j) continue;
        if (rows[j].sumpt > r.sumpt) ++rankS;
        if (!std::isnan(rows[j].trkptz_score) && !std::isnan(r.trkptz_score) &&
            rows[j].trkptz_score > r.trkptz_score) ++rankT;
        const double dt = std::abs(rows[j].cluster_time - r.cluster_time);
        const double dz = std::abs(rows[j].delta_z - r.delta_z);
        nearestDt = std::min(nearestDt, dt);
        nearestDz = std::min(nearestDz, dz);
        if (dt < 60.0) ++nWithin60;
      }
      r.n_clusters = (float)nC;
      r.sumpt_rank = (float)rankS;
      r.sumpt_frac_of_event = sumSumpt > 0 ? (float)(r.sumpt/sumSumpt) : NaNf;
      r.sumpt_ratio_to_max  = maxSumpt > 0 ? (float)(r.sumpt/maxSumpt) : NaNf;
      r.trkptz_rank = (float)rankT;
      r.trkptz_ratio_to_max = maxTrkptz > 0 ? (float)(r.trkptz_score/maxTrkptz) : NaNf;
      r.is_max_sumpt  = (rankS == 1) ? 1.f : 0.f;
      r.is_max_trkptz = (rankT == 1) ? 1.f : 0.f;
      r.dt_to_nearest_cluster = nC > 1 ? (float)nearestDt : NaNf;
      r.dz_to_nearest_cluster = nC > 1 ? (float)nearestDz : NaNf;
      r.n_clusters_within_60ps = (float)nWithin60;

      C = r;
      clusterTree->Fill();
      ++nClusterRows;
    }

    // ---- Per-track rows ---------------------------------------------------
    for (size_t ci = 0; ci < clusters.size(); ++ci) {
      const Cluster& cl = clusters[ci];
      if (cl.values.empty() || cl.trackIndices.empty()) continue;
      const double tClus = cl.values.at(0);
      const double clusZ = [&]{
        double zn=0, zd=0;
        for (int t : cl.trackIndices) { zn += branch.trackZ0[t]/branch.trackVarZ0[t];
                                        zd += 1.0/branch.trackVarZ0[t]; }
        return zn/zd;
      }();
      for (size_t k = 0; k < cl.trackIndices.size(); ++k) {
        const int trk = cl.trackIndices[k];
        T = TrackRow{};
        T.event_num = (float)evtNum; T.cluster_idx = (float)ci;
        T.track_idx = (float)trk;    T.sample_id = sampleId;
        T.pt = branch.trackPt[trk];  T.eta = branch.trackEta[trk];
        T.phi = branch.trackPhi[trk];T.theta = branch.trackTheta[trk];
        T.z0 = branch.trackZ0[trk];  T.d0 = branch.trackD0[trk];
        T.qOverP = branch.trackQP[trk];
        T.sigma_z0 = std::sqrt(branch.trackVarZ0[trk]);
        T.sigma_d0 = std::sqrt(branch.trackVarD0[trk]);
        T.sigma_qOverP = std::sqrt(branch.trackVarQp[trk]);
        T.time = branch.trackTime[trk]; T.timeRes = branch.trackTimeRes[trk];
        T.time_valid = branch.trackTimeValid[trk] != 0 ? 1.f : 0.f;
        T.quality = branch.trackQuality[trk];
        T.nhgtd_hits = branch.trackHgtdHits[trk];
        T.nhgtd_primary = branch.trackPrimHits[trk];
        T.z0_pull_pv = (float)((branch.trackZ0[trk]-vtxZ)/T.sigma_z0);
        T.t_pull_cluster = (T.time_valid > 0 && T.timeRes > 0)
                           ? (float)((cl.allTimes[k]-tClus)/T.timeRes) : NaNf;
        const NearJet nj = nearestForwardJet(trk, &branch);
        T.dr_nearest_fwdjet = nj.idx >= 0 ? (float)nj.dr : NaNf;
        T.pt_nearest_fwdjet = nj.idx >= 0 ? (float)nj.pt : NaNf;
        T.is_ghost_of_nearest = 0.f;
        if (nj.idx >= 0) {
          const auto& gh = branch.topoJetGhostTrackIdx[nj.idx];
          T.is_ghost_of_nearest =
              (std::find(gh.begin(), gh.end(), trk) != gh.end()) ? 1.f : 0.f;
        }
        T.is_lepton = (branch.trackLeptonID && branch.isGoodLepton(trk)) ? 1.f : 0.f;
        T.cluster_time = (float)tClus;
        T.cluster_delta_z = (float)(clusZ - vtxZ);
        T.truth_is_hs = (branch.trackToTruthvtx[trk] == 0) ? 1.f : 0.f;
        trackTree->Fill();
        ++nTrackRows;
      }
    }
  }

  out->cd();
  clusterTree->Write();
  trackTree->Write();
  out->Close();

  std::cout << "\n\n=== EXPORT SUMMARY ===\n"
            << "  sample          : " << (cfg.sampleName.empty() ? "local" : cfg.sampleName) << '\n'
            << "  events passing  : " << nPassEvent   << '\n'
            << "  cluster rows    : " << nClusterRows << '\n'
            << "  track rows      : " << nTrackRows   << '\n'
            << "  output          : " << outPath      << '\n'
            << "\n  Reminder: truth_* columns and delta_t are NOT model inputs.\n"
            << "  Group key for ranking is (sample_id, event_num).\n";
  return 0;
}
