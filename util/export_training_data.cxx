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
//   Output: one ROOT file per sample containing three TTrees
//     clusters : one row per cluster  (fixed-width features -> GBDT ranking)
//     tracks   : one row per track    (enables a per-track P(HS) model, whose
//                aggregate Sum pT*P(HS) is the learned form of what WAVeS
//                hard-codes as pT_jet/dR)
//     jets     : one row per forward jet. Jets previously reached a model only
//                through per-track association (dr_nearest_fwdjet, ...) and
//                event-level scalars (lead_jet_pt, n_forward_jets); as their own
//                token sequence a set/attention model can attend over them
//                directly, which no per-track summary can express.
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
//   EVENT IDENTITY - the grouping key is (sample_id, file_idx, event_num)
//   ------------------------------------------------------------------------
//   Every model here is a RANKER over the candidate clusters of ONE event, so a
//   grouping key that collides silently merges two events' candidates into one
//   ranking group. That is a training bug with no symptom: the rows are all
//   individually valid.
//
//   `event_num` is therefore the entry number WITHIN ITS FILE, and `file_idx`
//   is that file's position in the sample's sorted file list. Both are exact in
//   float32 (file counts are ~1e3, per-file entries ~1e3-1e4, against float's
//   2^24 integer ceiling). The previous key was the chain-global entry number
//   alone, which is only unique while a sample is exported by a single
//   unsharded process -- under --file-shard each shard builds its own chain and
//   restarts that counter at 0, so merging shards aliased shard 0's event 0
//   onto shard 1's.
//
//   `file_idx` is the index in the FULL sorted list, not in this shard's chain,
//   so it is shard-count invariant: the same event gets the same key whether the
//   sample was exported in 1 job or 12.
//
//   Usage:  ./export_training_data [--sample=<name>] [--max-events=N]
//                                  [--file-shard=<i>/<N>]
//                                  [--vbs-deta=<x>] [--vbs-mjj=<x>]
//
//   Sharding: --file-shard=<i>/<N> restricts the process to files i, i+N, ... of
//   the sorted list and tags the output name, exactly as clustering_hist does.
//   Merge the shards with `hadd`, NOT with util/hist_merge -- see the note at
//   the bottom of this file.
// ---------------------------------------------------------------------------

#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
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

// Nearest jet of ANY eta, unlike nearestForwardJet which is restricted to the
// HGTD forward band. The distinction is the point: in VBF or ttbar the
// hard-scatter tracks sit inside forward jets, which is what the WAVeS-style
// features already capture. In Z+jets the forward jets are usually PILEUP, so
// the hard-scatter content is whatever is left over -- central, or in no jet at
// all -- and no existing feature separates those two populations within a
// cluster. pT > MIN_JET_PT still applies; only the eta restriction is dropped.
NearJet nearestAnyJet(int trk, BranchPointerWrapper* b) {
  NearJet nj;
  const double te = b->trackEta[trk], tp = b->trackPhi[trk];
  for (int j = 0; j < (int)b->topoJetEta.GetSize(); ++j) {
    if (b->isJetRemoved(j)) continue;
    const double jp = b->topoJetPt[j];
    if (jp < MIN_JET_PT) continue;
    const double deta = b->topoJetEta[j] - te;
    const double dphi = TVector2::Phi_mpi_pi(b->topoJetPhi[j] - tp);
    const double dr   = std::hypot(deta, dphi);
    if (dr < nj.dr) { nj.dr = dr; nj.pt = jp; nj.idx = j; }
  }
  return nj;
}

// Relation of a z position to the reconstructed vertex collection. Index 0 is
// the primary by construction; everything else is a reconstructed pileup vertex.
struct VtxRel {
  float dzPU = NaNf;      // to the nearest PU vertex
  float closestIsPV = 0;  // nearest vertex OF ALL is the primary
  float n1 = 0, n2 = 0;   // PU vertices within 1 mm / 2 mm
};

// ---------------------------------------------------------------------------
// Canonical hard-scatter VERTEX-SELECTION scores, per reconstructed vertex.
//
// These are the two scores the experiment actually uses to pick the HS vertex,
// reproduced here as published rather than as this codebase's variants:
//
//   SumPt2(V)  = SUM pT^2                                  over V's tracks
//   WAVeS(V)   = SUM w_jet + SUM_mu w_mu + SUM_e w_e       over V's tracks
//     w_jet = (pT^track)^2 (pT^jet)^2 / dR(track, jet)
//     w_mu  = (pT^mu)^4 / 0.01,   w_e = (pT^e)^4 / 0.01
//
// NOTE how far this is from Score::WAVES in clustering_structs.h, which is a
// DIFFERENT score that happens to share the name. Ours uses pT^1 for both track
// and jet, floors dR at WAVES_DR_FLOOR, multiplies by exp(-1.5|dz|), carries no
// lepton term, and considers only qualifying FORWARD jets. Reporting ours as
// "WAVeS" against the literature would be wrong, which is why this exists.
//
// TRACK-TO-VERTEX ASSIGNMENT HERE IS AN APPROXIMATION, AND IT DOES NOT NEED TO
// BE. An earlier version of this comment claimed the ntuple carried no reco
// track-to-vertex association and that z0-significance assignment was the only
// option. That is FALSE: `Track_recoVtx_idx` holds the vertex fit's own
// assignment, `Track_recoVtx_weight` its per-track fit weight, and
// `RecoVtx_sumPt2` Athena's own SumPt^2 per vertex. None are bound in
// BranchPointerWrapper yet.
//
// The difference is not cosmetic. Measured on the local VBF sample, forward
// region (2.4 < |eta| < 4.0, pt > 1): tracks the fit assigns to vertex 0 are
// 58.0% truly hard-scatter, against ~32% for the incumbent z-significance cut,
// and tracks it assigns to a PILEUP vertex are 98.5% truly pileup -- a
// near-perfect veto. The fit uses full track parameters and covariances, so it
// still separates in the regime where |dz| cannot: at ~110 reconstructed
// vertices per event the nearest pileup vertex sits 0.024 mm from the primary,
// which is far inside the hard-scatter cluster's own z spread.
//
// Until those branches are bound, the scores below are a faithful reproduction
// of the FORMULAS on an approximate track partition, not of the experiment's
// numbers.
//
// The published w_jet has a bare 1/dR and no floor; a track collinear with a
// jet axis would diverge. DR_EPS guards that, and nDrBelowFloor counts how
// often dR lands under WAVES_DR_FLOOR so the guard's influence is measurable
// rather than assumed.
// ---------------------------------------------------------------------------
struct VtxScores {
  std::vector<double> sumpt2, waves;
  int    argmaxSumpt2 = -1, argmaxWaves = -1;
  long   nDrBelowFloor = 0, nTracksUsed = 0;
};

VtxScores computeVertexScores(BranchPointerWrapper* b) {
  constexpr double DR_EPS = 1e-3;
  const int nv = (int)b->recoVtxZ.GetSize();
  const int nt = (int)b->trackZ0.GetSize();
  VtxScores v;
  v.sumpt2.assign(nv, 0.0);
  v.waves.assign(nv, 0.0);
  if (nv == 0) return v;

  for (int t = 0; t < nt; ++t) {
    const double pt = b->trackPt[t];
    if (pt <= 0.0) continue;
    // Assign to the most z0-compatible vertex. sigma_z0 is per track, so this
    // is a significance and not a bare distance -- a badly measured track is
    // correctly allowed to belong to a vertex further away in mm.
    const double sz = std::sqrt(b->trackVarZ0[t]);
    if (!(sz > 0.0)) continue;
    int best = -1; double bestSig = 1e18;
    for (int j = 0; j < nv; ++j) {
      const double sig = std::abs(b->trackZ0[t] - b->recoVtxZ[j]) / sz;
      if (sig < bestSig) { bestSig = sig; best = j; }
    }
    if (best < 0) continue;
    ++v.nTracksUsed;
    v.sumpt2[best] += pt * pt;

    const NearJet nj = nearestAnyJet(t, b);
    if (nj.idx >= 0) {
      if (nj.dr < WAVES_DR_FLOOR) ++v.nDrBelowFloor;
      v.waves[best] += pt * pt * nj.pt * nj.pt / std::max(nj.dr, DR_EPS);
    }
    // Leptons enter as an ADDITIONAL term, not instead of the jet term: the
    // published sum is over tracks PLUS over identified muons and electrons.
    // w_mu and w_e are the same function of pT, so no flavour split is needed
    // numerically -- leptonPdg would give it if that ever changes.
    if (b->trackLeptonID && t < (int)b->trackLeptonID->GetSize() && b->isGoodLepton(t))
      v.waves[best] += std::pow(pt, 4) / 0.01;
  }
  v.argmaxSumpt2 = (int)(std::max_element(v.sumpt2.begin(), v.sumpt2.end()) - v.sumpt2.begin());
  v.argmaxWaves  = (int)(std::max_element(v.waves.begin(),  v.waves.end())  - v.waves.begin());
  return v;
}

// Rank of index `i` in `x` by descending value: 0 means nothing scores higher.
int descRank(const std::vector<double>& x, int i) {
  if (i < 0 || i >= (int)x.size()) return -1;
  int r = 0;
  for (size_t k = 0; k < x.size(); ++k) if (x[k] > x[i]) ++r;
  return r;
}

VtxRel vertexRelation(double z, BranchPointerWrapper* b) {
  VtxRel v;
  const int nv = (int)b->recoVtxZ.GetSize();
  double bestPU = 1e9, bestAll = 1e9; int bestAllIdx = -1;
  for (int i = 0; i < nv; ++i) {
    const double d = std::abs(z - b->recoVtxZ[i]);
    if (d < bestAll) { bestAll = d; bestAllIdx = i; }
    if (i == 0) continue;                       // index 0 is the primary
    if (d < bestPU) bestPU = d;
    if (d < 1.0) v.n1 += 1;
    if (d < 2.0) v.n2 += 1;
  }
  if (bestPU < 1e8) v.dzPU = (float)bestPU;
  v.closestIsPV = (bestAllIdx == 0) ? 1.f : 0.f;
  return v;
}

float median(std::vector<float> v) {
  if (v.empty()) return NaNf;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n % 2) ? v[n / 2] : 0.5f * (v[n / 2 - 1] + v[n / 2]);
}

// Precision-weighted sums over a cluster's TIMED tracks: W = SUM 1/sigma^2 and
// WT = SUM t/sigma^2. The cluster time is WT/W, so leave-one-out follows in
// closed form and needs no refit.
struct TimeSums { double W = 0.0, WT = 0.0; int n = 0; };

TimeSums timeSums(const Cluster& cl, BranchPointerWrapper* b) {
  TimeSums s;
  for (size_t k = 0; k < cl.trackIndices.size(); ++k) {
    const int trk = cl.trackIndices[k];
    const double sig = b->trackTimeRes[trk];
    if (b->trackTimeValid[trk] == 0 || sig <= 0.0) continue;
    const double w = 1.0 / (sig * sig);
    s.W += w; s.WT += w * cl.allTimes[k]; ++s.n;
  }
  return s;
}

// Cluster time with track k removed, minus the full cluster time. NaN when k has
// no usable time, or when it is the ONLY timed track -- removing it would leave
// the mean undefined rather than merely shifted, and reporting 0 there would
// read as "perfectly robust" when it is the opposite.
float looShift(const TimeSums& s, const Cluster& cl, BranchPointerWrapper* b, size_t k) {
  const int trk = cl.trackIndices[k];
  const double sig = b->trackTimeRes[trk];
  if (b->trackTimeValid[trk] == 0 || sig <= 0.0 || s.n < 2) return NaNf;
  const double w = 1.0 / (sig * sig);
  const double Wm = s.W - w;
  if (Wm <= 0.0) return NaNf;
  return (float)(s.WT / s.W - (s.WT - w * cl.allTimes[k]) / Wm);
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
  // (sample_id, file_idx, event_num) is the ranking group key -- see the
  // EVENT IDENTITY note at the top of this file.
  float event_num, file_idx, cluster_idx, sample_id, weight;
  // --- Tier 0: original cluster features ---------------------------------
  float cluster_time, delta_z, delta_z_resunits, cluster_z_sigma;
  float cluster_d0, cluster_d0_sigma, cluster_qOverP, cluster_qOverP_sigma;
  float sumpt, cluster_time_sigma, n_tracks;
  // --- Tier A: pT spectrum / hardness ------------------------------------
  float sumpt2, maxpt, pt_2nd, pt_3rd, lead_pt_frac, meanpt, medianpt;
  float n_pt_gt2, n_pt_gt5;
  // --- Tier B: timing coherence ------------------------------------------
  // LEAVE-ONE-OUT stability. The cluster time is a precision-weighted mean, so
  // dropping one track has a closed form -- no refit needed:
  //     t_-i = (SUM w_j t_j - w_i t_i) / (SUM w_j - w_i),  w = 1/sigma_t^2
  // The shift t - t_-i measures how much one track is dragging the answer. A
  // large max shift means the time rests on a single measurement, which is what
  // an outlier, a bad hit, or a mis-associated HGTD time looks like from the
  // reco side. Distinct from time_chi2_ndf: chi2 asks whether the tracks AGREE,
  // leave-one-out asks whether the answer is ROBUST, and a two-track cluster can
  // have a perfect chi2 while resting entirely on one track.
  float max_abs_loo_shift, mean_abs_loo_shift, loo_shift_lead_pt;
  // Cluster time SPLIT BY JET ASSOCIATION. Two precision-weighted means over
  // disjoint subsets of the cluster's timed tracks: those within dR < 0.4 of any
  // jet, and those in no jet at all. A cluster whose in-jet and out-of-jet halves
  // disagree in time is one whose tracks do not share an origin, which is what a
  // pileup jet absorbing hard-scatter tracks (or vice versa) looks like.
  // t_in_minus_out is that disagreement directly; the counts and sumpt say how
  // well-determined each half is, since a one-track half means nothing.
  float t_in_jet, t_in_jet_sigma, n_in_jet_timed, sumpt_in_jet_timed;
  float t_out_jet, t_out_jet_sigma, n_out_jet_timed, sumpt_out_jet_timed;
  float t_in_minus_out;
  float time_chi2_ndf, max_abs_tpull, n_tpull_gt2, time_spread;
  float min_timeRes, median_timeRes, max_timeRes;
  float n_valid_time, frac_valid_time, sumpt_valid_time;
  // --- Tier C: z compatibility with the reco PV --------------------------
  float z_chi2_ndf, max_abs_zpull, n_zpull_gt3, z_spread, median_z0_pull;
  // --- Tier C2: RELATION TO THE OTHER RECO VERTICES ----------------------
  // RecoVtx_z is a full array and until now only index 0 (the primary) was ever
  // read; every pileup vertex ITk reconstructed was discarded. That is the
  // information these features recover, and it targets the measured failure mode
  // directly: the model picks large pileup clusters, and a pileup cluster sits
  // at a z where ITk has ALREADY reconstructed a pileup vertex. "Is there a
  // reconstructed PU vertex sitting on top of me" is a reco-only question that
  // no existing feature asks.
  //
  // Modelled on the official-ntuple variables cluster_closest_vx_is_pvx,
  // cluster_distance_pu_vx, cluster_n_pu_vertices and local_vx_density.
  float dz_to_nearest_pu_vtx;    // |clusterZ - nearest PU (index>0) reco vertex|
  float closest_vtx_is_pv;       // 1 if the nearest reco vertex of all is index 0
  float n_pu_vtx_within_1mm, n_pu_vtx_within_2mm;
  float local_vtx_density;       // reco vertices per mm within +-2 mm of clusterZ
  // Canonical HS vertex-selection scores (see computeVertexScores). Two
  // families, answering different questions:
  //   clus_*        the published formulas applied to THIS cluster's tracks,
  //                 i.e. "how would the standard scores rank this candidate".
  //                 SumPt2 of the cluster is already exported as `sumpt2`.
  //   nearest_vtx_* the scores of the reconstructed vertex nearest this
  //                 cluster in z, i.e. "does this cluster sit on a vertex that
  //                 the standard selection thinks looks like a hard scatter".
  // The _frac forms are the score divided by the event's maximum, which is what
  // the model gets: raw values are GeV^4 and span many orders of magnitude.
  float clus_waves_canonical;
  float nearest_vtx_sumpt2, nearest_vtx_waves;
  float nearest_vtx_sumpt2_frac, nearest_vtx_waves_frac;
  float nearest_vtx_sumpt2_rank, nearest_vtx_waves_rank;
  float nearest_vtx_is_sumpt2_max, nearest_vtx_is_waves_max;
  // Event-level, repeated on every cluster row: does each canonical score pick
  // the primary? This is the standard algorithms' HS-vertex-selection
  // efficiency on our samples, which is a result in its own right.
  float evt_sumpt2_picks_pv, evt_waves_picks_pv;
  // ---- Vertex-fit assignment (Track_recoVtx_idx) ------------------------
  // What the vertex fit itself thinks these tracks belong to. Not a z proxy:
  // the fit uses full parameters and covariances, so it still separates where
  // |dz| is saturated by pileup density.
  float frac_trk_on_pv, sumpt_frac_on_pv, frac_trk_on_pu_vtx, frac_trk_unassigned;
  float n_distinct_vtx, dominant_vtx_frac, dominant_vtx_is_pv;
  float dominant_vtx_sumpt2, dominant_vtx_sumpt2_frac, dominant_vtx_sumpt2_rank;
  float mean_vtx_weight, sumpt_weighted_vtx_weight;
  float truth_dominant_vtx_is_hs;      // TRUTH (RecoVtx_isHS), diagnostics only
  // ---- Track quality, aggregated over the cluster -----------------------
  float mean_chi2_ndf, max_chi2_ndf, mean_n_pix_hits, mean_n_sct_hits;
  float mean_n_si_holes, frac_trk_with_shared_hits, mean_n_innermost_hits;
  float mean_btag_d0_sig, mean_btag_z0sin_sig;
  float pv_pu_dz_ratio;          // |dz to PV| / |dz to nearest PU vtx|
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
  // The same two times measured against truth. TRUTH-PREFIXED and therefore
  // never model inputs -- they exist to answer "which half was right", which is
  // the diagnostic the split was proposed to settle.
  float truth_dt_in_jet, truth_dt_out_jet;
  // Truth-vertex analogues, diagnostics only -- never inputs.
  float truth_dz_to_nearest_vtx, truth_local_vtx_density;
  float truth_purity;            // HS pT fraction of the cluster
  float truth_n_hs_tracks;       // HS tracks in this cluster
  float truth_hs_frac_tracks;    // HS track fraction
};

struct TrackRow {
  float event_num, file_idx, cluster_idx, track_idx, sample_id;
  float pt, eta, phi, theta, z0, d0, qOverP;
  float sigma_z0, sigma_d0, sigma_qOverP;
  float time, timeRes, time_valid;
  float quality, nhgtd_hits, nhgtd_primary;  // nhgtd_primary -> truth_ branch (MC-only)
  float z0_pull_pv, t_pull_cluster;
  // How far the cluster time moves when THIS track is removed, and that shift in
  // units of the cluster's own time uncertainty. NaN when the track has no valid
  // time, or when it is the only timed track (removing it leaves nothing).
  float loo_shift, loo_pull;
  // Association to ANY jet, not just a forward one. dr_nearest_fwdjet above is
  // NaN for a track with no qualifying FORWARD jet, which in Z+jets is most of
  // them -- so without this the pooling model cannot tell "outside every jet"
  // from "outside the forward jets but inside a central one".
  float dr_nearest_anyjet, pt_nearest_anyjet, is_in_any_jet;
  // Per-track distance to the nearest PILEUP reco vertex, and whether that is
  // closer than the primary. Deep Sets pools tracks, so giving it the relation
  // per track lets it learn "this cluster is built from tracks that each belong
  // to a reconstructed pileup vertex" -- which is what a pileup cluster IS.
  float dz_to_nearest_pu_vtx_trk, closer_to_pu_than_pv;
  float dr_nearest_fwdjet, pt_nearest_fwdjet, is_ghost_of_nearest;
  float is_lepton;
  float cluster_time, cluster_delta_z;   // context, so track rows stand alone
  // Broadcast copies of the cluster-level vertex relation. Deep Sets only ever
  // sees track rows, so a cluster-level feature that is not broadcast here is
  // invisible to it. These two are the ones that survived the conditional test:
  // binned in |delta_z|, the PU-vertex ratio still lifts discrimination
  // 1.78-2.34x in the displaced bins. closest_vtx_is_pv is deliberately NOT
  // broadcast -- conditioned on delta_z its lift collapses to ~1.2 and inverts
  // at large |delta_z|, so it restates delta_z rather than adding to it.
  float cluster_dz_to_nearest_pu_vtx, cluster_pv_pu_dz_ratio;
  // Broadcast canonical vertex-selection context. The RANK is the one that
  // earns its place: conditioned on |delta_z| it lifts 2.01x in the |dz|<0.1mm
  // bin, exactly where pv_pu_dz_ratio is flat (0.99x), so the two cover
  // opposite ends of the displacement range. closest_vtx_is_pv is included by
  // request; conditioned on delta_z it behaves like a restatement of it.
  float cluster_nearest_vtx_waves_rank, cluster_nearest_vtx_waves_frac;
  float cluster_closest_vtx_is_pv;
  // ---- vertex-fit assignment, per track ---------------------------------
  float on_pv, on_pu_vtx, vtx_unassigned, vtx_weight, vtx_sumpt2_frac;
  // ---- track fit quality, per track -------------------------------------
  float chi2_ndf, n_pix_hits, n_sct_hits, n_si_holes, n_shared_hits,
        n_innermost_hits, btag_d0_sig, btag_z0sin_sig;
  // ---- broadcast cluster-level vertex-fit context ------------------------
  float cluster_frac_trk_on_pv, cluster_dominant_vtx_is_pv,
        cluster_dominant_vtx_sumpt2_frac;
  float truth_is_hs;                     // LABEL for the per-track P(HS) model
  // Supervision, never an input: is the track's nearest forward jet matched to
  // a truth HS jet? Distinguishes "this track sits in a real HS jet" from "this
  // track sits in a pileup jet that happens to be forward" -- which the
  // canonical-selection measurement (results/baseline_deta0p0_mjj500p0.md) shows
  // is the majority case on Z+jets, where ~94% of the forward dijet is pileup.
  float truth_nearest_fwdjet_is_hs;
};

// One row per FORWARD jet (the same qualifying band WAVeS uses). Written so a
// set/attention model can treat jets as their own token sequence.
struct JetRow {
  float event_num, file_idx, jet_idx, sample_id;
  float pt, eta, phi, abs_eta;
  float n_ghost_tracks, sumpt_ghost;      // ghost tracks that survived selection
  float n_sel_tracks_dr04, sumpt_dr04;    // selected tracks within dR < 0.4
  float truth_is_hs;                      // matched to a truth HS jet
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
  MyUtl::FILE_SHARD = MyUtl::resolveShard(argc, argv);
  const MyUtl::Shard shard = MyUtl::FILE_SHARD;

  // Numeric sample tag travels with every row so concatenated files still group
  // correctly. The map lives in sample_config.h beside the registry so a newly
  // registered sample cannot silently alias onto an existing id -- see the
  // sampleId() comment there for the bug that motivated moving it.
  const float sampleId = (float)MyUtl::sampleId(cfg.sampleName);

  TChain chain("ntuple");
  // Bind the vertex-fit and track-quality branches. MUST precede the
  // BranchPointerWrapper below: the constructor reads this flag to decide what
  // to bind, so setting it afterwards silently leaves every extended branch
  // unbound and every derived column at its NaN/zero sentinel.
  MyUtl::EXTENDED_BRANCHES = true;

  setupChain(chain, cfg.ntupleDir.c_str(), shard);
  // Record what this sample actually carries, so the wrapper binds only real
  // branches. The productions are NOT uniform: local VBF has 241 branches, the
  // grid samples 183-195, and Track_btagIp_* is missing from every grid sample.
  if (auto* fl = chain.GetListOfFiles(); fl && fl->GetEntries() > 0) {
    std::unique_ptr<TFile> f0(TFile::Open(fl->At(0)->GetTitle()));
    if (f0 && !f0->IsZombie())
      if (auto* t0 = f0->Get<TTree>("ntuple"))
        for (auto* o : *t0->GetListOfBranches())
          MyUtl::AVAILABLE_BRANCHES.insert(o->GetName());
  }
  std::cout << "[branches] sample carries " << MyUtl::AVAILABLE_BRANCHES.size()
            << " branches\n";
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
  auto* jetTree     = new TTree("jets",     "one row per forward jet");
  ClusterRow C{};
  TrackRow   T{};
  JetRow     J{};
  auto BC = [&](const char* n, float* p) { clusterTree->Branch(n, p, Form("%s/F", n)); };
  auto BT = [&](const char* n, float* p) { trackTree  ->Branch(n, p, Form("%s/F", n)); };
  auto BJ = [&](const char* n, float* p) { jetTree    ->Branch(n, p, Form("%s/F", n)); };

  BC("event_num",&C.event_num); BC("file_idx",&C.file_idx);
  BC("cluster_idx",&C.cluster_idx);
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
  BC("max_abs_loo_shift",&C.max_abs_loo_shift);
  BC("mean_abs_loo_shift",&C.mean_abs_loo_shift);
  BC("loo_shift_lead_pt",&C.loo_shift_lead_pt);
  BC("t_in_jet",&C.t_in_jet); BC("t_in_jet_sigma",&C.t_in_jet_sigma);
  BC("n_in_jet_timed",&C.n_in_jet_timed); BC("sumpt_in_jet_timed",&C.sumpt_in_jet_timed);
  BC("t_out_jet",&C.t_out_jet); BC("t_out_jet_sigma",&C.t_out_jet_sigma);
  BC("n_out_jet_timed",&C.n_out_jet_timed); BC("sumpt_out_jet_timed",&C.sumpt_out_jet_timed);
  BC("t_in_minus_out",&C.t_in_minus_out);
  BC("time_chi2_ndf",&C.time_chi2_ndf); BC("max_abs_tpull",&C.max_abs_tpull);
  BC("n_tpull_gt2",&C.n_tpull_gt2); BC("time_spread",&C.time_spread);
  BC("min_timeRes",&C.min_timeRes); BC("median_timeRes",&C.median_timeRes);
  BC("max_timeRes",&C.max_timeRes); BC("n_valid_time",&C.n_valid_time);
  BC("frac_valid_time",&C.frac_valid_time); BC("sumpt_valid_time",&C.sumpt_valid_time);
  BC("z_chi2_ndf",&C.z_chi2_ndf); BC("max_abs_zpull",&C.max_abs_zpull);
  BC("n_zpull_gt3",&C.n_zpull_gt3); BC("z_spread",&C.z_spread);
  BC("median_z0_pull",&C.median_z0_pull);
  BC("dz_to_nearest_pu_vtx",&C.dz_to_nearest_pu_vtx);
  BC("closest_vtx_is_pv",&C.closest_vtx_is_pv);
  BC("n_pu_vtx_within_1mm",&C.n_pu_vtx_within_1mm);
  BC("n_pu_vtx_within_2mm",&C.n_pu_vtx_within_2mm);
  BC("local_vtx_density",&C.local_vtx_density);
  BC("clus_waves_canonical",&C.clus_waves_canonical);
  BC("nearest_vtx_sumpt2",&C.nearest_vtx_sumpt2);
  BC("nearest_vtx_waves",&C.nearest_vtx_waves);
  BC("nearest_vtx_sumpt2_frac",&C.nearest_vtx_sumpt2_frac);
  BC("nearest_vtx_waves_frac",&C.nearest_vtx_waves_frac);
  BC("nearest_vtx_sumpt2_rank",&C.nearest_vtx_sumpt2_rank);
  BC("nearest_vtx_waves_rank",&C.nearest_vtx_waves_rank);
  BC("nearest_vtx_is_sumpt2_max",&C.nearest_vtx_is_sumpt2_max);
  BC("nearest_vtx_is_waves_max",&C.nearest_vtx_is_waves_max);
  BC("evt_sumpt2_picks_pv",&C.evt_sumpt2_picks_pv);
  BC("evt_waves_picks_pv",&C.evt_waves_picks_pv);
  BC("frac_trk_on_pv",&C.frac_trk_on_pv);
  BC("sumpt_frac_on_pv",&C.sumpt_frac_on_pv);
  BC("frac_trk_on_pu_vtx",&C.frac_trk_on_pu_vtx);
  BC("frac_trk_unassigned",&C.frac_trk_unassigned);
  BC("n_distinct_vtx",&C.n_distinct_vtx);
  BC("dominant_vtx_frac",&C.dominant_vtx_frac);
  BC("dominant_vtx_is_pv",&C.dominant_vtx_is_pv);
  BC("dominant_vtx_sumpt2",&C.dominant_vtx_sumpt2);
  BC("dominant_vtx_sumpt2_frac",&C.dominant_vtx_sumpt2_frac);
  BC("dominant_vtx_sumpt2_rank",&C.dominant_vtx_sumpt2_rank);
  BC("mean_vtx_weight",&C.mean_vtx_weight);
  BC("sumpt_weighted_vtx_weight",&C.sumpt_weighted_vtx_weight);
  BC("truth_dominant_vtx_is_hs",&C.truth_dominant_vtx_is_hs);
  BC("mean_chi2_ndf",&C.mean_chi2_ndf); BC("max_chi2_ndf",&C.max_chi2_ndf);
  BC("mean_n_pix_hits",&C.mean_n_pix_hits); BC("mean_n_sct_hits",&C.mean_n_sct_hits);
  BC("mean_n_si_holes",&C.mean_n_si_holes);
  BC("frac_trk_with_shared_hits",&C.frac_trk_with_shared_hits);
  BC("mean_n_innermost_hits",&C.mean_n_innermost_hits);
  BC("mean_btag_d0_sig",&C.mean_btag_d0_sig);
  BC("mean_btag_z0sin_sig",&C.mean_btag_z0sin_sig);
  BC("pv_pu_dz_ratio",&C.pv_pu_dz_ratio);
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
  BC("truth_dt_in_jet",&C.truth_dt_in_jet); BC("truth_dt_out_jet",&C.truth_dt_out_jet);
  BC("truth_dz_to_nearest_vtx",&C.truth_dz_to_nearest_vtx);
  BC("truth_local_vtx_density",&C.truth_local_vtx_density);
  BC("truth_purity",&C.truth_purity); BC("truth_n_hs_tracks",&C.truth_n_hs_tracks);
  BC("truth_hs_frac_tracks",&C.truth_hs_frac_tracks);

  BT("event_num",&T.event_num); BT("file_idx",&T.file_idx);
  BT("cluster_idx",&T.cluster_idx);
  BT("track_idx",&T.track_idx); BT("sample_id",&T.sample_id);
  BT("pt",&T.pt); BT("eta",&T.eta); BT("phi",&T.phi); BT("theta",&T.theta);
  BT("z0",&T.z0); BT("d0",&T.d0); BT("qOverP",&T.qOverP);
  BT("sigma_z0",&T.sigma_z0); BT("sigma_d0",&T.sigma_d0); BT("sigma_qOverP",&T.sigma_qOverP);
  BT("time",&T.time); BT("timeRes",&T.timeRes); BT("time_valid",&T.time_valid);
  BT("quality",&T.quality); BT("nhgtd_hits",&T.nhgtd_hits);
  BT("truth_nhgtd_primary",&T.nhgtd_primary);
  BT("z0_pull_pv",&T.z0_pull_pv); BT("t_pull_cluster",&T.t_pull_cluster);
  BT("loo_shift",&T.loo_shift); BT("loo_pull",&T.loo_pull);
  BT("dz_to_nearest_pu_vtx_trk",&T.dz_to_nearest_pu_vtx_trk);
  BT("closer_to_pu_than_pv",&T.closer_to_pu_than_pv);
  BT("dr_nearest_anyjet",&T.dr_nearest_anyjet);
  BT("pt_nearest_anyjet",&T.pt_nearest_anyjet);
  BT("is_in_any_jet",&T.is_in_any_jet);
  BT("dr_nearest_fwdjet",&T.dr_nearest_fwdjet);
  BT("pt_nearest_fwdjet",&T.pt_nearest_fwdjet);
  BT("is_ghost_of_nearest",&T.is_ghost_of_nearest);
  BT("is_lepton",&T.is_lepton);
  BT("cluster_time",&T.cluster_time); BT("cluster_delta_z",&T.cluster_delta_z);
  BT("cluster_dz_to_nearest_pu_vtx",&T.cluster_dz_to_nearest_pu_vtx);
  BT("cluster_pv_pu_dz_ratio",&T.cluster_pv_pu_dz_ratio);
  BT("cluster_nearest_vtx_waves_rank",&T.cluster_nearest_vtx_waves_rank);
  BT("cluster_nearest_vtx_waves_frac",&T.cluster_nearest_vtx_waves_frac);
  BT("cluster_closest_vtx_is_pv",&T.cluster_closest_vtx_is_pv);
  BT("on_pv",&T.on_pv); BT("on_pu_vtx",&T.on_pu_vtx);
  BT("vtx_unassigned",&T.vtx_unassigned); BT("vtx_weight",&T.vtx_weight);
  BT("vtx_sumpt2_frac",&T.vtx_sumpt2_frac);
  BT("chi2_ndf",&T.chi2_ndf);
  BT("n_pix_hits",&T.n_pix_hits); BT("n_sct_hits",&T.n_sct_hits);
  BT("n_si_holes",&T.n_si_holes); BT("n_shared_hits",&T.n_shared_hits);
  BT("n_innermost_hits",&T.n_innermost_hits);
  BT("btag_d0_sig",&T.btag_d0_sig); BT("btag_z0sin_sig",&T.btag_z0sin_sig);
  BT("cluster_frac_trk_on_pv",&T.cluster_frac_trk_on_pv);
  BT("cluster_dominant_vtx_is_pv",&T.cluster_dominant_vtx_is_pv);
  BT("cluster_dominant_vtx_sumpt2_frac",&T.cluster_dominant_vtx_sumpt2_frac);
  BT("truth_is_hs",&T.truth_is_hs);
  BT("truth_nearest_fwdjet_is_hs",&T.truth_nearest_fwdjet_is_hs);

  BJ("event_num",&J.event_num); BJ("file_idx",&J.file_idx);
  BJ("jet_idx",&J.jet_idx);     BJ("sample_id",&J.sample_id);
  BJ("pt",&J.pt); BJ("eta",&J.eta); BJ("phi",&J.phi); BJ("abs_eta",&J.abs_eta);
  BJ("n_ghost_tracks",&J.n_ghost_tracks); BJ("sumpt_ghost",&J.sumpt_ghost);
  BJ("n_sel_tracks_dr04",&J.n_sel_tracks_dr04); BJ("sumpt_dr04",&J.sumpt_dr04);
  BJ("truth_is_hs",&J.truth_is_hs);

  // NO chain.GetEntries() here. It is not a cheap accessor: it opens EVERY file
  // in the chain to sum tree headers, which costs 0.3-1.2 s per file on the AF's
  // /data and was measured at 5-30 minutes of pre-loop dead time per condor job.
  // See the setupChain comment in src/event_processing.h -- the same call was
  // removed from clustering_hist/rpt_v5_hist for exactly this reason, and this
  // exporter was the last place still paying it. Progress therefore prints a
  // running count with no denominator, as the *_hist executables do.
  Long64_t nPassEvent = 0, nClusterRows = 0, nTrackRows = 0, nJetRows = 0;
  Long64_t nSeen = 0;
  std::cout << "Sample '" << (cfg.sampleName.empty() ? "local" : cfg.sampleName)
            << "' (id " << (int)sampleId << ")"
            << (shard.given ? " shard " + std::to_string(shard.index) + "/" +
                              std::to_string(shard.count) : "")
            << " -- exporting\n" << std::flush;

  while (reader.Next()) {
    // Entry number WITHIN THE CURRENT FILE, paired with that file's index in the
    // sample's full sorted list. A TChain-bound sequential TTreeReader's
    // GetReadEntry() is chain-global, so the per-file entry needs the chain
    // offset subtracted -- and GetTreeNumber() counts position within THIS
    // process's chain, which under sharding holds only files
    // shard.index, shard.index + count, ... of the global list.
    const Long64_t evtNum = chain.GetReadEntry() - chain.GetChainOffset();
    const int fileIdx = shard.active()
        ? shard.index + chain.GetTreeNumber() * shard.count
        : chain.GetTreeNumber();

    // Caps events PER SHARD, not per sample -- this is a local smoke-test knob,
    // so an N-shard run with --max-events=M exports up to N*M events.
    if (maxEvents > 0 && nSeen >= maxEvents) break;
    ++nSeen;
    if (nSeen % 20000 == 0)
      std::cout << "Progress: " << nSeen << " events read\r" << std::flush;

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
    // Once per event: the canonical HS vertex-selection scores for EVERY
    // reconstructed vertex. Hoisted because it is O(n_track x n_vertex) and
    // does not depend on the cluster.
    const VtxScores vsc = computeVertexScores(&branch);
    const double maxSumpt2 = vsc.argmaxSumpt2 >= 0 ? vsc.sumpt2[vsc.argmaxSumpt2] : 0.0;
    const double maxWaves  = vsc.argmaxWaves  >= 0 ? vsc.waves [vsc.argmaxWaves ] : 0.0;
    for (size_t ci = 0; ci < clusters.size(); ++ci) {
      const Cluster& cl = clusters[ci];
      if (cl.values.empty() || cl.trackIndices.empty()) continue;
      // NOTE: no nConstituents cut. Dropping small clusters would remove the
      // truth-closest candidate from some events and silently corrupt the
      // ranking label; the candidate set must match what the selector sees.

      ClusterRow R{};
      R.event_num = (float)evtNum; R.file_idx = (float)fileIdx;
      R.cluster_idx = (float)ci;
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

      // Cluster time split by jet association. Two disjoint precision-weighted
      // means over the cluster's TIMED tracks: inside dR < 0.4 of any jet, and
      // inside no jet. Same weights as the full cluster time, so the three are
      // directly comparable and t_in_minus_out is a like-for-like difference.
      {
        double Wi = 0, WTi = 0, Wo = 0, WTo = 0, ptI = 0, ptO = 0;
        int ni = 0, no = 0;
        for (size_t k = 0; k < cl.trackIndices.size(); ++k) {
          const int trk2 = cl.trackIndices[k];
          const double sg = branch.trackTimeRes[trk2];
          if (branch.trackTimeValid[trk2] == 0 || sg <= 0.0) continue;
          const double w = 1.0 / (sg * sg), tt = cl.allTimes[k];
          const NearJet aj = nearestAnyJet(trk2, &branch);
          if (aj.idx >= 0 && aj.dr < 0.4) {
            Wi += w; WTi += w * tt; ptI += branch.trackPt[trk2]; ++ni;
          } else {
            Wo += w; WTo += w * tt; ptO += branch.trackPt[trk2]; ++no;
          }
        }
        R.t_in_jet        = ni ? (float)(WTi / Wi) : NaNf;
        R.t_in_jet_sigma  = ni ? (float)(1.0 / std::sqrt(Wi)) : NaNf;
        R.n_in_jet_timed  = (float)ni;
        R.sumpt_in_jet_timed = (float)ptI;
        R.t_out_jet       = no ? (float)(WTo / Wo) : NaNf;
        R.t_out_jet_sigma = no ? (float)(1.0 / std::sqrt(Wo)) : NaNf;
        R.n_out_jet_timed = (float)no;
        R.sumpt_out_jet_timed = (float)ptO;
        // NaN unless BOTH halves exist -- a difference against a missing half is
        // not a small disagreement, it is no measurement, and 0 would read as
        // "the halves agree perfectly".
        R.t_in_minus_out  = (ni && no) ? (float)(WTi / Wi - WTo / Wo) : NaNf;
        const double truthT = branch.truthVtxTime[0];
        R.truth_dt_in_jet  = ni ? (float)(WTi / Wi - truthT) : NaNf;
        R.truth_dt_out_jet = no ? (float)(WTo / Wo - truthT) : NaNf;
      }

      // Leave-one-out summary for this cluster.
      {
        const TimeSums ts = timeSums(cl, &branch);
        double mx = 0.0, sum = 0.0; int cnt = 0; double leadPtSeen = -1.0;
        float leadShift = NaNf;
        for (size_t k = 0; k < cl.trackIndices.size(); ++k) {
          const float sh = looShift(ts, cl, &branch, k);
          if (std::isnan(sh)) continue;
          mx = std::max(mx, (double)std::abs(sh)); sum += std::abs(sh); ++cnt;
          const double pt_k = branch.trackPt[cl.trackIndices[k]];
          if (pt_k > leadPtSeen) { leadPtSeen = pt_k; leadShift = sh; }
        }
        R.max_abs_loo_shift  = cnt ? (float)mx : NaNf;
        R.mean_abs_loo_shift = cnt ? (float)(sum / cnt) : NaNf;
        R.loo_shift_lead_pt  = leadShift;
      }

      const int n = (int)cl.trackIndices.size();
      const double clusZ = znum/zden, clusZSig = 1.0/std::sqrt(zden);
      R.delta_z = (float)(clusZ - vtxZ);
      R.cluster_z_sigma = (float)clusZSig;
      {
        const VtxRel vr = vertexRelation(clusZ, &branch);
        R.dz_to_nearest_pu_vtx = vr.dzPU;
        R.closest_vtx_is_pv    = vr.closestIsPV;
        R.n_pu_vtx_within_1mm  = vr.n1;
        R.n_pu_vtx_within_2mm  = vr.n2;
        // vertices per mm over the +-2 mm window (the +1 counts the PV itself
        // when it falls inside, so a lone cluster far from everything reads 0).
        R.local_vtx_density    = (float)((vr.n2 + (std::abs(clusZ - vtxZ) < 2.0 ? 1 : 0)) / 4.0);
        // >1 means the cluster sits closer to a pileup vertex than to the PV,
        // which is the compact statement of "this is a pileup cluster". NaN
        // rather than a huge number when there is no PU vertex to compare to.
        R.pv_pu_dz_ratio = (std::isnan(vr.dzPU) || vr.dzPU <= 0.0)
                           ? NaNf : (float)(std::abs(clusZ - vtxZ) / vr.dzPU);
        double bt = 1e9, dsum = 0;
        for (int tv = 0; tv < (int)branch.truthVtxZ.GetSize(); ++tv) {
          const double d = std::abs(clusZ - branch.truthVtxZ[tv]);
          bt = std::min(bt, d);
          if (d < 2.0) dsum += 1;
        }
        R.truth_dz_to_nearest_vtx  = bt < 1e8 ? (float)bt : NaNf;
        R.truth_local_vtx_density  = (float)(dsum / 4.0);

        // Canonical scores of the reco vertex nearest this cluster in z. The
        // _frac forms are what the model sees; raw GeV^4 spans decades.
        int nv_i = -1; double nv_d = 1e18;
        for (int j = 0; j < (int)branch.recoVtxZ.GetSize(); ++j) {
          const double d = std::abs(clusZ - branch.recoVtxZ[j]);
          if (d < nv_d) { nv_d = d; nv_i = j; }
        }
        if (nv_i >= 0) {
          R.nearest_vtx_sumpt2 = (float)vsc.sumpt2[nv_i];
          R.nearest_vtx_waves  = (float)vsc.waves [nv_i];
          R.nearest_vtx_sumpt2_frac = maxSumpt2 > 0.0
                                      ? (float)(vsc.sumpt2[nv_i] / maxSumpt2) : NaNf;
          R.nearest_vtx_waves_frac  = maxWaves  > 0.0
                                      ? (float)(vsc.waves [nv_i] / maxWaves ) : NaNf;
          R.nearest_vtx_sumpt2_rank = (float)descRank(vsc.sumpt2, nv_i);
          R.nearest_vtx_waves_rank  = (float)descRank(vsc.waves,  nv_i);
          R.nearest_vtx_is_sumpt2_max = (nv_i == vsc.argmaxSumpt2) ? 1.f : 0.f;
          R.nearest_vtx_is_waves_max  = (nv_i == vsc.argmaxWaves ) ? 1.f : 0.f;
        } else {
          R.nearest_vtx_sumpt2 = R.nearest_vtx_waves = NaNf;
          R.nearest_vtx_sumpt2_frac = R.nearest_vtx_waves_frac = NaNf;
          R.nearest_vtx_sumpt2_rank = R.nearest_vtx_waves_rank = NaNf;
          R.nearest_vtx_is_sumpt2_max = R.nearest_vtx_is_waves_max = NaNf;
        }
        // Vertex 0 is the primary by construction, so "picks the PV" is the
        // canonical algorithm getting the HS vertex right on this event.
        R.evt_sumpt2_picks_pv = (vsc.argmaxSumpt2 == 0) ? 1.f : 0.f;
        R.evt_waves_picks_pv  = (vsc.argmaxWaves  == 0) ? 1.f : 0.f;

        // The published WAVeS formula applied to THIS cluster's own tracks.
        double cw = 0.0;
        for (int t : cl.trackIndices) {
          const double ptt = branch.trackPt[t];
          const NearJet nj = nearestAnyJet(t, &branch);
          if (nj.idx >= 0)
            cw += ptt * ptt * nj.pt * nj.pt / std::max(nj.dr, 1e-3);
          if (branch.trackLeptonID && t < (int)branch.trackLeptonID->GetSize()
              && branch.isGoodLepton(t))
            cw += std::pow(ptt, 4) / 0.01;
        }
        R.clus_waves_canonical = (float)cw;

        // ---- vertex-fit assignment over this cluster's tracks -------------
        {
          std::unordered_map<int,double> ptByVtx;   // vertex -> summed pT
          std::unordered_map<int,int>    nByVtx;
          double sPV = 0.0, sAll = 0.0, wSum = 0.0, wPtSum = 0.0;
          int nPV = 0, nPU = 0, nUn = 0, nW = 0;
          for (int t : cl.trackIndices) {
            const int    vi = branch.recoVtxOf(t);
            const double pt = branch.trackPt[t];
            sAll += pt;
            if (vi == 0)      { ++nPV; sPV += pt; }
            else if (vi > 0)  { ++nPU; }
            else              { ++nUn; }
            if (vi >= 0) { ptByVtx[vi] += pt; nByVtx[vi] += 1; }
            const float w = branch.recoVtxWeightOf(t);
            if (!std::isnan(w)) { wSum += w; wPtSum += w * pt; ++nW; }
          }
          const double nT = (double)cl.trackIndices.size();
          R.frac_trk_on_pv       = (float)(nPV / nT);
          R.frac_trk_on_pu_vtx   = (float)(nPU / nT);
          R.frac_trk_unassigned  = (float)(nUn / nT);
          R.sumpt_frac_on_pv     = sAll > 0 ? (float)(sPV / sAll) : NaNf;
          R.n_distinct_vtx       = (float)ptByVtx.size();
          R.mean_vtx_weight      = nW ? (float)(wSum / nW) : NaNf;
          R.sumpt_weighted_vtx_weight = sAll > 0 ? (float)(wPtSum / sAll) : NaNf;
          // The cluster's DOMINANT vertex is the one holding the most of its pT
          // -- a pT-weighted vote, not a track count, so one high-pT track is
          // not outvoted by a handful of soft ones.
          int dom = -1; double domPt = -1.0;
          for (const auto& kv : ptByVtx)
            if (kv.second > domPt) { domPt = kv.second; dom = kv.first; }
          if (dom >= 0) {
            R.dominant_vtx_frac   = sAll > 0 ? (float)(domPt / sAll) : NaNf;
            R.dominant_vtx_is_pv  = (dom == 0) ? 1.f : 0.f;
            R.dominant_vtx_sumpt2 = branch.vtxSumPt2(dom);
            R.truth_dominant_vtx_is_hs = branch.vtxIsHS(dom) ? 1.f : 0.f;
            if (branch.recoVtxSumPt2) {
              const int nv = (int)branch.recoVtxSumPt2->GetSize();
              double mx = 0.0; int rank = 0;
              for (int v = 0; v < nv; ++v) {
                const double sv = (*branch.recoVtxSumPt2)[v];
                mx = std::max(mx, sv);
                if (sv > R.dominant_vtx_sumpt2) ++rank;
              }
              R.dominant_vtx_sumpt2_frac = mx > 0 ? (float)(R.dominant_vtx_sumpt2 / mx) : NaNf;
              R.dominant_vtx_sumpt2_rank = (float)rank;
            }
          } else {
            R.dominant_vtx_frac = R.dominant_vtx_is_pv = NaNf;
            R.dominant_vtx_sumpt2 = R.dominant_vtx_sumpt2_frac = NaNf;
            R.dominant_vtx_sumpt2_rank = R.truth_dominant_vtx_is_hs = NaNf;
          }
        }

        // ---- track fit quality, aggregated ---------------------------------
        {
          double c2=0, c2mx=0, pix=0, sct=0, hol=0, inn=0, d0s=0, z0s=0;
          int nc=0, nh=0, nsh=0, nd=0, nz=0;
          for (int t : cl.trackIndices) {
            if (branch.trackChi2 && branch.trackNdof
                && t < (int)branch.trackNdof->GetSize()
                && (*branch.trackNdof)[t] > 0) {
              const double r = (*branch.trackChi2)[t] / (*branch.trackNdof)[t];
              c2 += r; c2mx = std::max(c2mx, r); ++nc;
            }
            auto U = [&](const std::unique_ptr<TTreeReaderArray<unsigned char>>& a) {
              return (a && t < (int)a->GetSize()) ? (double)(*a)[t] : 0.0; };
            pix += U(branch.trackNPixHits);  sct += U(branch.trackNSctHits);
            hol += U(branch.trackNPixHoles) + U(branch.trackNSctHoles);
            inn += U(branch.trackNInnerHits);
            if (U(branch.trackNPixShared) + U(branch.trackNSctShared) > 0) ++nsh;
            ++nh;
            if (branch.trackBtagD0 && branch.trackBtagD0Unc
                && t < (int)branch.trackBtagD0Unc->GetSize()
                && (*branch.trackBtagD0Unc)[t] > 0) {
              d0s += std::abs((*branch.trackBtagD0)[t]) / (*branch.trackBtagD0Unc)[t]; ++nd;
            }
            if (branch.trackBtagZ0Sin && branch.trackBtagZ0SinUnc
                && t < (int)branch.trackBtagZ0SinUnc->GetSize()
                && (*branch.trackBtagZ0SinUnc)[t] > 0) {
              z0s += std::abs((*branch.trackBtagZ0Sin)[t]) / (*branch.trackBtagZ0SinUnc)[t]; ++nz;
            }
          }
          R.mean_chi2_ndf   = nc ? (float)(c2/nc) : NaNf;
          R.max_chi2_ndf    = nc ? (float)c2mx    : NaNf;
          R.mean_n_pix_hits = nh ? (float)(pix/nh) : NaNf;
          R.mean_n_sct_hits = nh ? (float)(sct/nh) : NaNf;
          R.mean_n_si_holes = nh ? (float)(hol/nh) : NaNf;
          R.mean_n_innermost_hits = nh ? (float)(inn/nh) : NaNf;
          R.frac_trk_with_shared_hits = nh ? (float)((double)nsh/nh) : NaNf;
          R.mean_btag_d0_sig    = nd ? (float)(d0s/nd) : NaNf;
          R.mean_btag_z0sin_sig = nz ? (float)(z0s/nz) : NaNf;
        }
      }
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
      const TimeSums ts = timeSums(cl, &branch);
      const double clusTSig = ts.W > 0.0 ? 1.0 / std::sqrt(ts.W) : -1.0;
      // Hoisted: this is a scan over the event's reco vertices and does not
      // depend on k, so computing it inside the track loop would repeat it once
      // per constituent track.
      const VtxRel cvr = vertexRelation(clusZ, &branch);
      const float cvrRatio = (std::isnan(cvr.dzPU) || cvr.dzPU <= 0.0)
                             ? NaNf : (float)(std::abs(clusZ - vtxZ) / cvr.dzPU);
      float cvNvRank = NaNf, cvNvFrac = NaNf;
      {
        int nvi = -1; double nvd = 1e18;
        for (int j = 0; j < (int)branch.recoVtxZ.GetSize(); ++j) {
          const double d = std::abs(clusZ - branch.recoVtxZ[j]);
          if (d < nvd) { nvd = d; nvi = j; }
        }
        if (nvi >= 0) {
          cvNvRank = (float)descRank(vsc.waves, nvi);
          cvNvFrac = maxWaves > 0.0 ? (float)(vsc.waves[nvi] / maxWaves) : NaNf;
        }
      }
      for (size_t k = 0; k < cl.trackIndices.size(); ++k) {
        const int trk = cl.trackIndices[k];
        T = TrackRow{};
        T.event_num = (float)evtNum; T.file_idx = (float)fileIdx;
        T.cluster_idx = (float)ci;
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
        T.truth_nearest_fwdjet_is_hs =
            nj.idx >= 0 ? (branch.isJetTruthHS(nj.idx) ? 1.f : 0.f) : NaNf;
        T.is_ghost_of_nearest = 0.f;
        if (nj.idx >= 0) {
          const auto& gh = branch.topoJetGhostTrackIdx[nj.idx];
          T.is_ghost_of_nearest =
              (std::find(gh.begin(), gh.end(), trk) != gh.end()) ? 1.f : 0.f;
        }
        T.is_lepton = (branch.trackLeptonID && branch.isGoodLepton(trk)) ? 1.f : 0.f;
        {
          const VtxRel tvr = vertexRelation(branch.trackZ0[trk], &branch);
          T.dz_to_nearest_pu_vtx_trk = tvr.dzPU;
          const double dpv = std::abs(branch.trackZ0[trk] - vtxZ);
          T.closer_to_pu_than_pv = (std::isnan(tvr.dzPU)) ? NaNf
                                   : (tvr.dzPU < dpv ? 1.f : 0.f);
        }
        const NearJet anyj = nearestAnyJet(trk, &branch);
        T.dr_nearest_anyjet = anyj.idx >= 0 ? (float)anyj.dr : NaNf;
        T.pt_nearest_anyjet = anyj.idx >= 0 ? (float)anyj.pt : NaNf;
        T.is_in_any_jet = (anyj.idx >= 0 && anyj.dr < 0.4) ? 1.f : 0.f;
        T.loo_shift = looShift(ts, cl, &branch, k);
        T.loo_pull = (std::isnan(T.loo_shift) || clusTSig <= 0.0)
                     ? NaNf : (float)(T.loo_shift / clusTSig);
        T.cluster_time = (float)tClus;
        T.cluster_delta_z = (float)(clusZ - vtxZ);
        T.cluster_dz_to_nearest_pu_vtx = cvr.dzPU;
        T.cluster_pv_pu_dz_ratio = cvrRatio;
        T.cluster_nearest_vtx_waves_rank = cvNvRank;
        T.cluster_nearest_vtx_waves_frac = cvNvFrac;
        T.cluster_closest_vtx_is_pv = cvr.closestIsPV;
        {
          const int vi = branch.recoVtxOf(trk);
          T.on_pv           = (vi == 0) ? 1.f : 0.f;
          T.on_pu_vtx       = (vi >  0) ? 1.f : 0.f;
          T.vtx_unassigned  = (vi <  0) ? 1.f : 0.f;
          T.vtx_weight      = branch.recoVtxWeightOf(trk);
          T.vtx_sumpt2_frac = NaNf;
          if (vi >= 0 && branch.recoVtxSumPt2) {
            double mx = 0.0;
            for (int v = 0; v < (int)branch.recoVtxSumPt2->GetSize(); ++v)
              mx = std::max(mx, (double)(*branch.recoVtxSumPt2)[v]);
            if (mx > 0) T.vtx_sumpt2_frac = (float)(branch.vtxSumPt2(vi) / mx);
          }
          auto U = [&](const std::unique_ptr<TTreeReaderArray<unsigned char>>& a) {
            return (a && trk < (int)a->GetSize()) ? (float)(*a)[trk] : NaNf; };
          T.chi2_ndf = (branch.trackChi2 && branch.trackNdof
                        && trk < (int)branch.trackNdof->GetSize()
                        && (*branch.trackNdof)[trk] > 0)
                       ? (float)((*branch.trackChi2)[trk] / (*branch.trackNdof)[trk]) : NaNf;
          T.n_pix_hits = U(branch.trackNPixHits);
          T.n_sct_hits = U(branch.trackNSctHits);
          const float ph = U(branch.trackNPixHoles), sh = U(branch.trackNSctHoles);
          T.n_si_holes = (std::isnan(ph)||std::isnan(sh)) ? NaNf : ph + sh;
          const float ps = U(branch.trackNPixShared), ss = U(branch.trackNSctShared);
          T.n_shared_hits = (std::isnan(ps)||std::isnan(ss)) ? NaNf : ps + ss;
          T.n_innermost_hits = U(branch.trackNInnerHits);
          T.btag_d0_sig = (branch.trackBtagD0 && branch.trackBtagD0Unc
                           && trk < (int)branch.trackBtagD0Unc->GetSize()
                           && (*branch.trackBtagD0Unc)[trk] > 0)
                          ? (float)(std::abs((*branch.trackBtagD0)[trk])
                                    / (*branch.trackBtagD0Unc)[trk]) : NaNf;
          T.btag_z0sin_sig = (branch.trackBtagZ0Sin && branch.trackBtagZ0SinUnc
                              && trk < (int)branch.trackBtagZ0SinUnc->GetSize()
                              && (*branch.trackBtagZ0SinUnc)[trk] > 0)
                             ? (float)(std::abs((*branch.trackBtagZ0Sin)[trk])
                                       / (*branch.trackBtagZ0SinUnc)[trk]) : NaNf;
          T.cluster_frac_trk_on_pv           = rows[ci].frac_trk_on_pv;
          T.cluster_dominant_vtx_is_pv       = rows[ci].dominant_vtx_is_pv;
          T.cluster_dominant_vtx_sumpt2_frac = rows[ci].dominant_vtx_sumpt2_frac;
        }
        T.truth_is_hs = (branch.trackToTruthvtx[trk] == 0) ? 1.f : 0.f;
        trackTree->Fill();
        ++nTrackRows;
      }
    }

    // ---- Per-forward-jet rows ---------------------------------------------
    // Same qualifying band as nearestForwardJet / WAVeS, so a jet token and the
    // per-track jet association describe the same object set.
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      if (branch.isJetRemoved(j)) continue;
      const double jp = branch.topoJetPt[j];
      const double je = branch.topoJetEta[j];
      if (jp < MIN_JET_PT) continue;
      if (std::abs(je) < MIN_ABS_ETA_JET || std::abs(je) > MAX_ABS_ETA_JET) continue;

      J = JetRow{};
      J.event_num = (float)evtNum; J.file_idx = (float)fileIdx;
      J.jet_idx = (float)j;        J.sample_id = sampleId;
      J.pt = (float)jp; J.eta = (float)je; J.phi = (float)branch.topoJetPhi[j];
      J.abs_eta = (float)std::abs(je);

      // Ghost tracks are intersected with the event's SELECTED track list: a
      // model can only ever see selected tracks, so an unrestricted ghost count
      // would be a jet feature with no counterpart among the track tokens.
      const auto& gh = branch.topoJetGhostTrackIdx[j];
      double nGh = 0, sumGh = 0, nDr = 0, sumDr = 0;
      for (int t : tracks) {
        if (std::find(gh.begin(), gh.end(), t) != gh.end()) {
          nGh += 1; sumGh += branch.trackPt[t];
        }
        const double deta = je - branch.trackEta[t];
        const double dphi = TVector2::Phi_mpi_pi(branch.topoJetPhi[j] - branch.trackPhi[t]);
        if (std::hypot(deta, dphi) < 0.4) { nDr += 1; sumDr += branch.trackPt[t]; }
      }
      J.n_ghost_tracks = (float)nGh;    J.sumpt_ghost = (float)sumGh;
      J.n_sel_tracks_dr04 = (float)nDr; J.sumpt_dr04 = (float)sumDr;
      J.truth_is_hs = branch.isJetTruthHS(j) ? 1.f : 0.f;
      jetTree->Fill();
      ++nJetRows;
    }
  }

  out->cd();
  clusterTree->Write();
  trackTree->Write();
  jetTree->Write();
  out->Close();

  std::cout << "\n\n=== EXPORT SUMMARY ===\n"
            << "  sample          : " << (cfg.sampleName.empty() ? "local" : cfg.sampleName)
            << "  (id " << (int)sampleId << ")\n"
            << "  shard           : "
            << (shard.given ? std::to_string(shard.index) + "/" + std::to_string(shard.count)
                            : std::string("none (all files)")) << '\n'
            << "  events read     : " << nSeen        << '\n'
            << "  events passing  : " << nPassEvent   << '\n'
            << "  cluster rows    : " << nClusterRows << '\n'
            << "  track rows      : " << nTrackRows   << '\n'
            << "  jet rows        : " << nJetRows     << '\n'
            << "  output          : " << outPath      << '\n'
            << "\n  Reminder: truth_* columns and delta_t are NOT model inputs.\n"
            << "  Group key for ranking is (sample_id, file_idx, event_num).\n"
            << "\n  Merge shards with `hadd`, NOT util/hist_merge: this file holds\n"
            << "  only TTrees, which hadd merges correctly by concatenating rows.\n"
            << "  hist_merge exists for the *_hist outputs, whose TParameter\n"
            << "  scalars hadd would silently take from the first input only --\n"
            << "  there are no such scalars here.\n";
  return 0;
}
