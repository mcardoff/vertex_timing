#ifndef JET_TAGGING_H
#define JET_TAGGING_H

// ---------------------------------------------------------------------------
// jet_tagging.h
//   In-code JVT / fJVT pileup-jet taggers, shared by util/vbs_region_diag.cxx
//   (VBS pair composition) and util/rpt_v5_hist.cxx (region R_pT), so the two
//   remove exactly the same jets before forming the VBS pair.
//
//   The SuperNtuples carry no Jvt/fJvt decoration (241 branches, none match, on
//   either jet collection), so both are computed from the vertex fit's own
//   track assignment (Track_recoVtx_idx -- an EXTENDED branch: set
//   EXTENDED_BRANCHES and call recordAvailableBranches, or recoVtxOf returns a
//   sentinel and every jet fails JVT) and each jet's ghost-associated tracks:
//
//     R_pT       = sum pT(ghost tracks fitted to vertex 0) / pT_jet   [1510.03823]
//     corrJVF    = pT_PV / (pT_PV + pT_PU / (k n_PU^trk)),  k = 0.01  [stored only]
//     fJVT       = max_{i>0} (p_T^miss,i . j_T) / pT_jet              [1705.02211]
//     p_T^miss,i = -1/2 ( sum_{trk fitted to i, |eta|<2.5} p_T
//                       + sum_{central jets whose dominant ghost vertex is i} p_T )
//
//   JVT proper is a k-NN likelihood over (corrJVF, R_pT) that cannot be
//   rebuilt from the ntuple, so the "JVT" cut is R_pT alone, with thresholds
//   set so the paper-HS efficiency on local-VBF central jets inside the JVT
//   window reproduces the published working-point efficiencies (JVT_WPS).
//   fJVT thresholds are the published ones. Windows: JVT |eta| < 2.5 and
//   pT < 60 GeV; fJVT 2.5 <= |eta| < 4.5 and pT < 120 GeV. A jet outside both
//   windows is never removed; a jet with no ghost tracks has R_pT = 0 and
//   fails JVT in its window.
//
//   Note the JVT R_pT here (vertex-FIT assignment, all ghost tracks) is not
//   rpt_v5's R_pT (z0 association, dR < 0.2), which is the discriminant that
//   study measures; they share only the idea.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "clustering_constants.h"
#include "clustering_structs.h"

namespace MyUtl {

  inline constexpr double JVT_ETA_MAX         = 2.5;    // JVT window: |eta| < 2.5 ...
  inline constexpr double JVT_PT_MAX          = 60.0;   // ... and pT < 60 GeV
  inline constexpr double FJVT_ETA_MIN        = 2.5;    // fJVT window: 2.5 <= |eta| < 4.5 ...
  inline constexpr double FJVT_ETA_MAX        = 4.5;
  inline constexpr double FJVT_PT_MAX         = 120.0;  // ... and pT < 120 GeV
  inline constexpr double FJVT_CEN_JET_PT_MIN = 20.0;   // central jets entering p_T^miss,i
  inline constexpr double FJVT_TRK_ETA_MAX    = 2.5;    // tracks entering p_T^miss,i
  inline constexpr double CORRJVF_K           = 0.01;

  struct JvtWP {
    const char* name;     // --jvt= value
    const char* tag;      // output-file tag
    double      rptMin;   // JVT proxy: keep a JVT-window jet iff R_pT >= this
    double      fjvtMax;  // keep an fJVT-window jet iff fJVT <= this
  };
  // R_pT thresholds: calibrated on local-VBF paper-HS jets in the JVT window
  // (|eta| < 2.5, 30 < pT < 60 GeV) to the published JVT working-point HS
  // efficiencies -- loose <-> the default Run-2 "Medium" point (92%), tight <->
  // "Tight" (85%) -- by python/vbs_jvt_calibrate.py on the `jets` tree of a
  // vbs_region_diag --jvt=none run. fJVT: published Loose 0.5 / Tight 0.4.
  inline const JvtWP JVT_WPS[] = {
    {"none",  "",          -1.0, 1e9},
    {"loose", "_jvtLoose", 0.0167, 0.5},   // 92.0% HS eff on local VBF, PU eff 13.2%
    {"tight", "_jvtTight", 0.0767, 0.4},   // 85.0% HS eff on local VBF, PU eff  1.9%
  };

  // The working point called `name`, or nullptr if there is none.
  inline const JvtWP* jvtWPByName(const std::string& name) {
    for (const auto& wp : JVT_WPS) if (name == wp.name) return &wp;
    return nullptr;
  }

  // One jet's discriminants, computed for EVERY jet in the array (not just the
  // pT-passing ones) because the fJVT of a forward jet depends on the central
  // jets' vertex assignment.
  struct JetDisc {
    float rpt = 0.f, corrjvf = -1.f, fjvt = 0.f;
    int   nGhost = 0, nGhostPV = 0;
    bool  jvtWin = false, fjvtWin = false;
  };

  // The removal rule, in one place. The windows are disjoint in |eta|, so a
  // jet can fail at most one of the two.
  inline bool jvtRemoves (const JetDisc& d, const JvtWP& wp) { return d.jvtWin  && d.rpt  < wp.rptMin; }
  inline bool fjvtRemoves(const JetDisc& d, const JvtWP& wp) { return d.fjvtWin && d.fjvt > wp.fjvtMax; }

  inline std::vector<JetDisc> computeJetDiscriminants(const BranchPointerWrapper& b) {
    const int nJ = (int)b.topoJetPt.GetSize();
    const int nV = (int)b.recoVtxZ.GetSize();
    std::vector<JetDisc> D(nJ);

    // Per-vertex transverse momentum of the tracks the fit assigned to it,
    // central tracks only (fJVT's p_T^miss is a central quantity). Index -1 is
    // "fitted to no vertex" (47% of tracks) and contributes nowhere.
    std::vector<double> vpx(nV, 0.0), vpy(nV, 0.0);
    int nPUtrk = 0;
    for (int t = 0; t < (int)b.trackPt.GetSize(); ++t) {
      const int v = b.recoVtxOf(t);
      if (v < 0 || v >= nV) continue;
      if (v > 0) ++nPUtrk;
      if (std::abs((double)b.trackEta[t]) >= FJVT_TRK_ETA_MAX) continue;
      vpx[v] += b.trackPt[t] * std::cos(b.trackPhi[t]);
      vpy[v] += b.trackPt[t] * std::sin(b.trackPhi[t]);
    }

    // Per-jet ghost-track sums split by vertex; a central jet is then assigned
    // to whichever vertex dominates its ghost pT, and if that is a PILEUP vertex
    // the jet's pT enters that vertex's p_T^miss. PV-dominated (i.e. hard
    // scatter) jets enter no PU vertex's balance. Overlap-removed jets are
    // leptons and are skipped for the assignment only.
    std::vector<double> perV(nV);
    for (int j = 0; j < nJ; ++j) {
      std::fill(perV.begin(), perV.end(), 0.0);
      double sumPV = 0.0, sumPU = 0.0;
      JetDisc& d = D[j];
      d.nGhost = (int)b.topoJetGhostTrackIdx[j].size();
      for (int idx : b.topoJetGhostTrackIdx[j]) {
        const int v = b.recoVtxOf(idx);
        if (v < 0 || v >= nV) continue;
        perV[v] += b.trackPt[idx];
        if (v == 0) { sumPV += b.trackPt[idx]; ++d.nGhostPV; }
        else          sumPU += b.trackPt[idx];
      }
      const double pt = b.topoJetPt[j], aeta = std::abs((double)b.topoJetEta[j]);
      d.rpt     = (float)(sumPV / pt);
      d.corrjvf = (sumPV + sumPU > 0.0)
                ? (float)(sumPV / (sumPV + sumPU / (CORRJVF_K * std::max(nPUtrk, 1))))
                : -1.f;
      d.jvtWin  = aeta < JVT_ETA_MAX && pt < JVT_PT_MAX;
      d.fjvtWin = aeta >= FJVT_ETA_MIN && aeta < FJVT_ETA_MAX && pt < FJVT_PT_MAX;
      if (aeta < FJVT_TRK_ETA_MAX && pt > FJVT_CEN_JET_PT_MIN && !b.isJetRemoved(j)) {
        int vBest = -1; double best = 0.0;
        for (int v = 0; v < nV; ++v) if (perV[v] > best) { best = perV[v]; vBest = v; }
        if (vBest > 0) {
          vpx[vBest] += pt * std::cos(b.topoJetPhi[j]);
          vpy[vBest] += pt * std::sin(b.topoJetPhi[j]);
        }
      }
    }

    // fJVT: the pileup vertex whose missing transverse momentum points most
    // along the jet, normalised to the jet pT. Computed for every jet; the
    // window is applied by the caller.
    for (int j = 0; j < nJ; ++j) {
      const double pt = b.topoJetPt[j];
      const double ux = std::cos(b.topoJetPhi[j]), uy = std::sin(b.topoJetPhi[j]);
      double best = -std::numeric_limits<double>::infinity();
      for (int v = 1; v < nV; ++v) {
        const double proj = -0.5 * (vpx[v] * ux + vpy[v] * uy) / pt;
        if (proj > best) best = proj;
      }
      D[j].fjvt = (nV > 1) ? (float)best : 0.f;
    }
    return D;
  }

}  // namespace MyUtl

#endif  // JET_TAGGING_H
