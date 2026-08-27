// ---------------------------------------------------------------------------
// dz_alpha_sweep.cxx — sweep the coefficient alpha in the cluster-selection
// z-penalty, for four arms at once.
//
//   RES_CLUS   (sum pT_i)                       * e^(-b*s_eff)  * e^(-a*|dz_clus|)
//   RES_TRK     sum_i pT_i * e^(-b*dz(eta_i,pT_i))               * e^(-a*|dz_clus|)
//   ... and the same two for the WAVeS jet-proximity weight.
//
// TWO terms, not one, and BOTH are swept: the e^(-a*|dz_clus|) displacement
// penalty is production (alpha = 1.5) but that value was never validated --
// util/scratch/trkptz_alpha_sweep.cxx says so in its own header -- and a
// beta-only scan at alpha = 1.5 measured alpha to be past TRKPTZ's optimum.
//
// (alpha, beta) = (1.5, 0) reproduces the production score EXACTLY, so the
// grid contains its own null hypothesis and any improvement is read directly
// against that cell.
//
// Ranked on core fraction only: core sigma was measured flat to +/-0.1 ps
// across every arm and every coefficient in the 1D scans, so it does not
// discriminate. Counters rather than histograms for the same reason -- the
// grid is 13 x 31 x 4 arms x 3 regions and only needs a pass/total ratio.
//
// What the new factor is for: the association cut has already spent the
// measured displacement -- getAssociatedTracks admits a track when
// |dz|/sigma_z0 < MAX_NSIGMA, so a large-sigma track passes over a wider
// window in z and is that much more likely to be pileup that fell inside it.
// dz(eta,pT) is the width of that window, so e^(-b*dz) down-weights the tracks
// whose presence is the weakest evidence that the cluster is the hard-scatter
// vertex. On its own it would carry no vertex information at all (it depends
// only on eta and pT, not on the event), which is exactly why it multiplies
// the displacement penalty rather than replacing it.
//
// s_eff is the pT-weighted mean of getNewDzpara over the cluster's tracks
// (Cluster::clusterDzPara) -- deliberately multiplicity-independent, unlike
// the cluster's own combined sigma_z which shrinks as 1/sqrt(N).
//
// WHY THIS EXISTS AS A TRACKED EXECUTABLE. util/scratch/trkptz_resunits_sweep.cxx
// already swept almost exactly this, was never committed, and left no output
// PDF, no comment, no doc entry and no memory note -- so the question is still
// open years later. util/scratch/ is gitignored and is where results go to die
// in this repo. This lives in util/ and its findings go into CLAUDE.md.
//
// Related: util/scratch/trkptz_alpha_sweep.cxx states in its own header that
// the production alpha = 1.5 "has never been evaluated". BASE is swept here
// too, so the comparison is best-vs-best rather than best-vs-incumbent.
//
// Clustering runs ONCE per event at DIST_CUT_CONE, exactly as production does.
// Everything an arm needs is then cached per cluster (ClusterInputs below), so
// each additional alpha costs only arithmetic -- no re-clustering, no re-read.
//
//   ./dz_alpha_sweep [--sample=<name>] [--max-events=N]
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "sample_config.h"

#include <iomanip>
#include <iostream>
#include <vector>

using namespace MyUtl;

namespace {

// ── Sweep grid ──────────────────────────────────────────────────────────────
// alpha: displacement coefficient, mm^-1. Production is 1.5; 0 = no
// displacement penalty at all, which the 1D scan measured as costing ~1.3 %pt.
constexpr double AL_MIN = 0.0, AL_MAX = 3.0, AL_STEP = 0.25;
// beta: resolution coefficient, mm^-1. 0 = production (no resolution term).
constexpr double BE_MIN = 0.0, BE_MAX = 1.5, BE_STEP = 0.05;

// Region split on forward HS track count. The association study found gains
// concentrate at low multiplicity, and the failure decomposition says that is
// where the failures live, so the split is worth carrying.
constexpr int HS_SPLIT = 5;
enum Region { R_ALL = 0, R_LOW, R_HIGH, N_REGION };
const char* REGION_NAME[N_REGION] = {"all", "nHS<=5", "nHS>5"};

enum Arm { A_RES_CLUS = 0, A_RES_TRK, W_RES_CLUS, W_RES_TRK, N_ARM };
const char* ARM_NAME[N_ARM] = {"TRKPTZ res/clus", "TRKPTZ res/trk",
                               "WAVES  res/clus", "WAVES  res/trk"};

// Everything a rescoring needs, cached once per cluster so the alpha loop is
// pure arithmetic. The WAVeS jet-proximity weight is the expensive part (a
// nested loop over jets per track), so it is computed here exactly once.
struct ClusterInputs {
  double sumPt      = 0.0;   // Score::TRKPT equivalent
  double sumJetW    = 0.0;   // WAVeS jet-proximity sum
  double dzClus     = 0.0;   // |cluster z - reco vtx z|, mm
  double sigmaEff   = 1.0;   // pT-weighted mean getNewDzpara, mm
  double time       = 0.0;   // cluster time, for the Delta t metric
  std::vector<double> pt, jetW, dzTk, dzPara;   // per-track
};

double armScore(const ClusterInputs& c, Arm arm, double a, double b) {
  const double disp = std::exp(-a * c.dzClus);
  switch (arm) {
    case A_RES_CLUS: return c.sumPt   * std::exp(-b * c.sigmaEff) * disp;
    case W_RES_CLUS: return c.sumJetW * std::exp(-b * c.sigmaEff) * disp;
    case A_RES_TRK: {
      double s = 0.0;
      for (size_t i = 0; i < c.pt.size(); ++i)
        s += c.pt[i] * std::exp(-b * c.dzPara[i]);
      return s * disp;
    }
    case W_RES_TRK: {
      double s = 0.0;
      for (size_t i = 0; i < c.pt.size(); ++i)
        s += c.jetW[i] * std::exp(-b * c.dzPara[i]);
      return s * disp;
    }
    default: return 0.0;
  }
}

}  // namespace

int main(int argc, char** argv) {
  TH1::AddDirectory(kFALSE);
  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::resolveSelection(argc, argv);
  MyUtl::OVERLAP_REMOVAL = sample.overlapRemoval;
  const Long64_t maxEvents = MyUtl::resolveMaxEvents(argc, argv);
  gErrorIgnoreLevel = kFatal;

  std::vector<double> AL, BE;
  for (int i = 0; AL_MIN + i * AL_STEP <= AL_MAX + 1e-9; ++i) AL.push_back(AL_MIN + i * AL_STEP);
  for (int i = 0; BE_MIN + i * BE_STEP <= BE_MAX + 1e-9; ++i) BE.push_back(BE_MIN + i * BE_STEP);
  const int nAL = (int)AL.size(), nBE = (int)BE.size();
  std::cout << "[sweep] alpha " << nAL << " pts [" << AL_MIN << "," << AL_MAX << "]"
            << "  beta " << nBE << " pts [" << BE_MIN << "," << BE_MAX << "]"
            << "  x " << N_ARM << " arms x " << N_REGION << " regions\n";

  // Pass/total counters only -- see the header for why histograms are not used.
  const size_t NC = (size_t)N_ARM * nAL * nBE * N_REGION;
  std::vector<Long64_t> nPass(NC, 0);
  auto IDX = [&](int arm, int ia, int ib, int r) {
    return (((size_t)arm * nAL + ia) * nBE + ib) * N_REGION + r; };
  std::vector<Long64_t> nEvt(N_REGION, 0);

  // ── Migration bookkeeping ────────────────────────────────────────────────
  // A net +X %pt can be a pure gain or a large two-way churn that happens to
  // net positive; only the second is fragile. Tally the 2x2 for the production
  // cell against one challenger cell, per event, so gained and lost are
  // separately visible. Also split the sample in half (even/odd events) so the
  // net can be checked for stability rather than assumed.
  auto findCell = [](const std::vector<double>& v, double x) {
    for (int i = 0; i < (int)v.size(); ++i) if (std::abs(v[i] - x) < 1e-9) return i;
    return -1; };
  const int iaProd = findCell(AL, 1.50), ibProd = findCell(BE, 0.00);
  const int iaNew  = findCell(AL, 0.75), ibNew  = findCell(BE, 0.30);
  if (iaProd < 0 || ibProd < 0 || iaNew < 0 || ibNew < 0) {
    std::cerr << "migration cells are not on the grid\n"; return 1; }
  // mig[arm][region][passedProduction][passedChallenger]
  std::vector<Long64_t> mig((size_t)N_ARM * N_REGION * 4, 0);
  auto MIG = [&](int arm, int r, int p, int q) {
    return (((size_t)arm * N_REGION + r) * 2 + p) * 2 + q; };
  std::vector<Long64_t> half((size_t)N_ARM * 2 * 2, 0);   // [arm][half][prod|new]
  std::vector<Long64_t> halfN(2, 0);

  TChain chain("ntuple");
  setupChain(chain, sample.ntupleDir.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  MyUtl::PhaseTimer phase;
  Long64_t seen = 0;
  while (reader.Next()) {
    ++seen;
    if (maxEvents > 0 && seen > maxEvents) break;
    if (seen % 20000 == 0) phase.mark("processed " + std::to_string(seen));

    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    if (!branch.passBasicCuts())       continue;
    if (!branch.passJetPtCut())        continue;

    auto tracks = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);
    if (tracks.empty()) continue;
    EventCounts ev(&branch, tracks, /*checkValidTimes=*/true);
    auto clusters = clusterTracksInTime(tracks, &branch, DIST_CUT_CONE,
                                        false, true, IDEAL_TRACK_RES,
                                        ClusteringMethod::ITERATIVE, false, false, false);
    if (clusters.empty()) continue;

    std::vector<ClusterInputs> ci;
    ci.reserve(clusters.size());
    for (const auto& c : clusters) {
      ClusterInputs x;
      x.time     = c.values.at(0);
      x.sigmaEff = c.clusterDzPara(&branch);
      double znum = 0.0, zden = 0.0;
      for (int idx : c.trackIndices) {
        double varZ = branch.trackVarZ0[idx];
        if (varZ > 0.0) { znum += branch.trackZ0[idx] / varZ; zden += 1.0 / varZ; }
      }
      x.dzClus = (zden > 0.0) ? std::abs(znum / zden - branch.recoVtxZ[0]) : 0.0;
      for (int idx : c.trackIndices) {
        double pt = branch.trackPt[idx];
        x.sumPt += pt;
        x.pt.push_back(pt);
        x.dzPara.push_back(getNewDzpara(branch.trackEta[idx], pt));
        double minDR = 1e6, nearJetPt = 0.0;
        for (int j = 0; j < (int)branch.topoJetEta.GetSize(); ++j) {
          if (branch.isJetRemoved(j)) continue;
          if (branch.topoJetPt[j] < MIN_JET_PT) continue;
          double ae = std::abs((double)branch.topoJetEta[j]);
          if (ae < MIN_ABS_ETA_JET || ae > MAX_ABS_ETA_JET) continue;
          double dr = std::hypot(branch.topoJetEta[j] - branch.trackEta[idx],
                                 TVector2::Phi_mpi_pi(branch.topoJetPhi[j] - branch.trackPhi[idx]));
          if (dr < minDR) { minDR = dr; nearJetPt = branch.topoJetPt[j]; }
        }
        double w = (nearJetPt > 0.0) ? pt * nearJetPt / std::max(minDR, WAVES_DR_FLOOR) : 0.0;
        x.jetW.push_back(w);
        x.sumJetW += w;
      }
      ci.push_back(std::move(x));
    }

    const double tTruth = branch.truthVtxTime[0];
    const int    region = (ev.nForwardTrackHS <= HS_SPLIT) ? R_LOW : R_HIGH;
    ++nEvt[R_ALL]; ++nEvt[region];

    const int hlf = (int)(seen % 2);
    ++halfN[hlf];
    for (int arm = 0; arm < N_ARM; ++arm) {
      const bool wavesArm = (arm >= W_RES_CLUS);
      bool passProd = false, passNew = false;
      for (int ia = 0; ia < nAL; ++ia)
        for (int ib = 0; ib < nBE; ++ib) {
          int best = 0; double bestS = -1.0;
          for (size_t k = 0; k < ci.size(); ++k) {
            Arm eff = (Arm)arm;
            // No qualifying forward jets anywhere in this cluster: WAVeS has no
            // jet-proximity information, so fall back to its TRKPTZ counterpart
            // exactly as Cluster::updateScores does.
            if (wavesArm && ci[k].sumJetW <= 0.0) eff = (Arm)(arm - W_RES_CLUS);
            double s = armScore(ci[k], eff, AL[ia], BE[ib]);
            if (s > bestS) { bestS = s; best = (int)k; }
          }
          const bool passed = std::abs(ci[best].time - tTruth) < PASS_SIGMA;
          if (passed) {
            ++nPass[IDX(arm, ia, ib, R_ALL)];
            ++nPass[IDX(arm, ia, ib, region)];
          }
          if (ia == iaProd && ib == ibProd) passProd = passed;
          if (ia == iaNew  && ib == ibNew ) passNew  = passed;
        }
      ++mig[MIG(arm, R_ALL,  passProd, passNew)];
      ++mig[MIG(arm, region, passProd, passNew)];
      if (passProd) ++half[((size_t)arm * 2 + hlf) * 2 + 0];
      if (passNew)  ++half[((size_t)arm * 2 + hlf) * 2 + 1];
    }
  }
  phase.mark("event loop done");
  std::cout << "events: all=" << nEvt[R_ALL] << "  low=" << nEvt[R_LOW]
            << "  high=" << nEvt[R_HIGH] << "\n";

  auto CF = [&](int arm, int ia, int ib, int r) {
    return nEvt[r] ? 100.0 * nPass[IDX(arm, ia, ib, r)] / nEvt[r] : 0.0; };
  for (int r = 0; r < N_REGION; ++r) {
    const double prodTrk = CF(A_RES_TRK, iaProd, ibProd, r);
    const double prodWav = CF(W_RES_TRK, iaProd, ibProd, r);
    std::cout << "\n=== REGION " << REGION_NAME[r] << "  (" << nEvt[r] << " events)"
              << "   production(a=1.5,b=0): TRKPTZ " << std::fixed << std::setprecision(2)
              << prodTrk << "%  WAVES " << prodWav << "% ===\n";
    for (int arm = 0; arm < N_ARM; ++arm) {
      int bia = 0, bib = 0; double bf = -1.0;
      for (int ia = 0; ia < nAL; ++ia)
        for (int ib = 0; ib < nBE; ++ib)
          if (CF(arm, ia, ib, r) > bf) { bf = CF(arm, ia, ib, r); bia = ia; bib = ib; }
      const double ref = (arm >= W_RES_CLUS) ? prodWav : prodTrk;
      std::cout << "  " << std::left << std::setw(18) << ARM_NAME[arm]
                << " best alpha=" << std::fixed << std::setprecision(2) << AL[bia]
                << " beta=" << std::setprecision(2) << BE[bib]
                << "  coreFrac=" << std::setprecision(2) << bf
                << "%  (" << std::showpos << std::setprecision(2) << (bf - ref)
                << " %pt vs production)" << std::noshowpos << "\n";
    }
  }

  // Headline arm grid, so the shape of the optimum is visible rather than just
  // its location -- a sharp peak and a broad plateau warrant different trust.
  for (int r = 0; r < N_REGION; ++r) {
    std::cout << "\n--- TRKPTZ res/trk core fraction [%], region " << REGION_NAME[r]
              << "  (rows alpha, cols beta) ---\n" << std::setw(7) << "a\\b";
    for (int ib = 0; ib < nBE; ib += 2)
      std::cout << std::setw(7) << std::fixed << std::setprecision(2) << BE[ib];
    std::cout << "\n";
    for (int ia = 0; ia < nAL; ++ia) {
      std::cout << std::setw(7) << std::fixed << std::setprecision(2) << AL[ia];
      for (int ib = 0; ib < nBE; ib += 2)
        std::cout << std::setw(7) << std::setprecision(2) << CF(A_RES_TRK, ia, ib, r);
      std::cout << "\n";
    }
  }

  std::cout << "\n=== MIGRATION: production (a=1.50, b=0.00) -> challenger (a=0.75, b=0.30) ===\n"
            << "  Is the gain a pure addition, or a two-way churn that nets positive?\n\n"
            << std::left << std::setw(18) << "arm" << std::setw(10) << "region"
            << std::right << std::setw(9) << "gained" << std::setw(9) << "lost"
            << std::setw(9) << "net" << std::setw(11) << "net %pt"
            << std::setw(12) << "churn %pt" << std::setw(10) << "gain:loss" << "\n"
            << std::string(88, '-') << "\n";
  for (int arm = 0; arm < N_ARM; ++arm)
    for (int r = 0; r < N_REGION; ++r) {
      const Long64_t gained = mig[MIG(arm, r, 0, 1)];   // failed prod, passes new
      const Long64_t lost   = mig[MIG(arm, r, 1, 0)];   // passed prod, fails new
      const double   denom  = nEvt[r] ? (double)nEvt[r] : 1.0;
      std::cout << std::left << std::setw(18) << ARM_NAME[arm] << std::setw(10) << REGION_NAME[r]
                << std::right << std::setw(9) << gained << std::setw(9) << lost
                << std::setw(9) << (gained - lost)
                << std::fixed << std::setprecision(2)
                << std::setw(10) << std::showpos << 100.0 * (gained - lost) / denom << std::noshowpos
                << std::setw(12) << 100.0 * (gained + lost) / denom
                << std::setw(10) << (lost ? (double)gained / lost : 0.0) << "\n";
    }
  std::cout << "\n  gained = failed production, passes challenger.  lost = the reverse.\n"
               "  churn  = (gained+lost)/N: how many events the change TOUCHES at all.\n"
               "  A strict improvement would have lost = 0.\n";

  std::cout << "\n=== SPLIT-HALF STABILITY (even vs odd events) ===\n"
            << std::left << std::setw(18) << "arm"
            << std::right << std::setw(14) << "net %pt (even)"
            << std::setw(14) << "net %pt (odd)" << "\n" << std::string(46, '-') << "\n";
  for (int arm = 0; arm < N_ARM; ++arm) {
    std::cout << std::left << std::setw(18) << ARM_NAME[arm] << std::right << std::fixed
              << std::setprecision(2) << std::showpos;
    for (int hh = 0; hh < 2; ++hh) {
      double pr = halfN[hh] ? 100.0 * half[((size_t)arm * 2 + hh) * 2 + 0] / halfN[hh] : 0.0;
      double nw = halfN[hh] ? 100.0 * half[((size_t)arm * 2 + hh) * 2 + 1] / halfN[hh] : 0.0;
      std::cout << std::setw(14) << (nw - pr);
    }
    std::cout << std::noshowpos << "\n";
  }

  std::cout << "\nCSV|region,arm,alpha,beta,core_frac_pct,n_events\n";
  for (int r = 0; r < N_REGION; ++r)
    for (int arm = 0; arm < N_ARM; ++arm)
      for (int ia = 0; ia < nAL; ++ia)
        for (int ib = 0; ib < nBE; ++ib)
          std::cout << "CSV|" << REGION_NAME[r] << "," << ARM_NAME[arm] << ","
                    << std::fixed << std::setprecision(2) << AL[ia] << "," << BE[ib] << ","
                    << CF(arm, ia, ib, r) << "," << nEvt[r] << "\n";
  return 0;
}
