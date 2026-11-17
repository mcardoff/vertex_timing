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
#include <iostream>
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
  // The forward R_pT extension (--jvt=looseext): ITk tracks reach |eta| = 4.0,
  // so R_pT can be measured over the whole fJVT window up to there. Forward
  // jets with 2.5 <= |eta| < 4.0 and pT < FJVT_PT_MAX then need R_pT >= rptMin
  // AS WELL AS fJVT <= fjvtMax. Motivated by a 104 GeV merged-pileup tagging
  // jet at eta 2.80 with R_pT = 0.009 that passed fJVT (0.16): ~90% of the
  // pileup tagging legs surviving fJVT loose have R_pT < 0.0167.
  inline constexpr double RPT_FWD_ETA_MAX     = 4.0;

  struct JvtWP {
    const char* name;     // --jvt= value
    const char* tag;      // output-file tag
    double      rptMin;   // JVT proxy: keep a JVT-window jet iff R_pT >= this
    double      fjvtMax;  // keep an fJVT-window jet iff fJVT <= this
    bool        rptFwd = false;  // also cut R_pT on forward jets (RPT_FWD_ETA_MAX)
    double      fwdPtMax = FJVT_PT_MAX;   // ... below this pT (rptFwd only)
    double      rptMinFwd = -1.0;         // ... at this threshold (< 0: rptMin)
    int         t0Gate = -1;              // >= 0: cut a t0-GATED R_pT instead, t0 source tdr::T0Idx (the caller
                                          // computes it -- util/tdr_common.h's jvtRpt -- and puts it in JetDisc::rpt)
    int         tdrRpt = -2;              // > -2: in the FORWARD window cut the TDR R_pT instead of JVT's (ghost tracks
                                          // z-associated to PV0, util/tdr_common.h's jetRpt): -1 ITk only, k >= 0
                                          // t0-gated with t0 source k. The caller computes it (uncapped track pT)
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
    // loose + the same R_pT cut over 2.5 <= |eta| < 4.0, pT < 120 GeV (NOT
    // recalibrated forward). Appended last: the indices above are stored as
    // meta_vbs_jvt_wp / jvt_wp in existing outputs.
    {"looseext", "_jvtLooseExt", 0.0167, 0.5, true},
    // JVT alone, its R_pT cut extended to |eta| < 4.0 and NO fJVT (fjvtMax 1e9 never removes): the analysis baseline
    // the user was given on 2026-11-06 ("remove fJVT and extend JVT to |eta| < 4.0", "no fJVT at all"). Forward pT
    // range: JVT's own (< 60 GeV) or fJVT's (< 120 GeV), the definition not being fixed yet. Appended last, as above.
    {"jvt4pt60",  "_jvt4pt60",  0.0167, 1e9, true, 60.0},
    {"jvt4pt120", "_jvt4pt120", 0.0167, 1e9, true, 120.0},
    // Two ways to make jvt4pt60 harder on forward pileup (2026-11-09 mini deck; the user: "lets make the time-gate
    // JVT, and tighter forward threshold"). At <mu> = 200 a forward pileup jet passes the loose cut on ONE PV0 track.
    //   t0: the same cut on the t0-gated R_pT (2 sigma, cleaned times; untimed tracks count) -- WAVeS t0, truth t0
    //   fNN: a tighter forward threshold, R_pT >= 0.NN at 2.5 <= |eta| < 4.0 (central unchanged)
    {"jvt4pt60t0",    "_jvt4pt60t0",    0.0167, 1e9, true, 60.0, -1.0, 2},
    {"jvt4pt60t0tru", "_jvt4pt60t0tru", 0.0167, 1e9, true, 60.0, -1.0, 3},
    {"jvt4pt60f03",   "_jvt4pt60f03",   0.0167, 1e9, true, 60.0, 0.03},
    {"jvt4pt60f05",   "_jvt4pt60f05",   0.0167, 1e9, true, 60.0, 0.05},
    {"jvt4pt60f10",   "_jvt4pt60f10",   0.0167, 1e9, true, 60.0, 0.10},
    // Timing in the tagger itself (the user, 2026-11-14: "I want timing to implemented in the first row ... establish a
    // WP at the same efficiency as JVT's 0.0167"). JVT's R_pT cannot carry it: it is 0 or >= ~0.015 (one PV0-fitted
    // track), so 0.0167 means "has a PV0-fitted track" and a gated version cannot be matched back to its efficiency.
    // Instead the forward window (2.5 <= |eta| < 4.0, pT < 60) cuts the TDR R_pT at JVT's forward HS efficiency on grid
    // VBF (84.2%; condor/vbf/vbf_room_rptptmaxnone_tdr_pujet.root, uncapped tracks): ITk only, WAVeS-gated,
    // truth-gated. Central jets keep JVT's 0.0167. Needs --rpt-trk-ptmax=none. results/jvt_timing_wp.md
    {"jvt4pt60tdr",    "_jvt4pt60tdr",    0.0167, 1e9, true, 60.0, 0.2317, -1, -1},
    {"jvt4pt60tdrw",   "_jvt4pt60tdrw",   0.0167, 1e9, true, 60.0, 0.2049, -1, 2},
    {"jvt4pt60tdrtru", "_jvt4pt60tdrtru", 0.0167, 1e9, true, 60.0, 0.1962, -1, 3},
  };

  // The working point called `name`, or nullptr if there is none.
  inline const JvtWP* jvtWPByName(const std::string& name) {
    for (const auto& wp : JVT_WPS) if (name == wp.name) return &wp;
    return nullptr;
  }

  // True for a working point whose cut needs an R_pT computeJetDiscriminants does not compute: a t0-gated one
  // (t0Gate >= 0) or, in the forward window, the TDR one (tdrRpt > -2). The caller has to compute it and put it in
  // JetDisc::rpt, which only util/tdr_fwdpair.cxx does. Anywhere else jvtRemoves would cut JVT's plain R_pT at
  // thresholds calibrated on a different quantity.
  inline bool jvtNeedsCallerRpt(const JvtWP& wp) { return wp.t0Gate >= 0 || wp.tdrRpt > -2; }

  // The working points computeJetDiscriminants alone can evaluate, "none|loose|...", for messages.
  inline std::string jvtSupportedNames() {
    std::string s;
    for (const auto& wp : JVT_WPS)
      if (!jvtNeedsCallerRpt(wp)) s += (s.empty() ? "" : "|") + std::string(wp.name);
    return s;
  }

  // --jvt=<name> for a program that evaluates the taggers with computeJetDiscriminants alone (rpt_v5_hist and
  // rpt_v5_plot, vbs_region_diag, vbs_time_veto): the last --jvt= given, "none" without one. An unknown name, or a
  // working point needing a caller-computed R_pT (jvtNeedsCallerRpt), prints an error naming the supported working
  // points and returns nullptr; the caller exits non-zero. `prog` prefixes the message.
  inline const JvtWP* resolveJvtWP(int argc, char** argv, const std::string& prog) {
    const JvtWP* sel = &JVT_WPS[0];
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      if (a.rfind("--jvt=", 0) != 0) continue;
      const JvtWP* wp = jvtWPByName(a.substr(6));
      if (!wp) {
        std::cerr << "[" << prog << "] unknown " << a << "; supported: " << jvtSupportedNames() << "\n";
        return nullptr;
      }
      if (jvtNeedsCallerRpt(*wp)) {
        std::cerr << "[" << prog << "] " << a << " cuts a t0-gated or TDR R_pT, which this program does not compute "
                     "(only tdr_fwdpair does); supported: " << jvtSupportedNames() << "\n";
        return nullptr;
      }
      sel = wp;
    }
    return sel;
  }

  // One jet's discriminants, computed for EVERY jet in the array (not just the
  // pT-passing ones) because the fJVT of a forward jet depends on the central
  // jets' vertex assignment.
  struct JetDisc {
    float rpt = 0.f, corrjvf = -1.f, fjvt = 0.f;
    int   nGhost = 0, nGhostPV = 0;
    bool  jvtWin = false, fjvtWin = false;
    bool  rptFwdWin = false;   // 2.5 <= |eta| < RPT_FWD_ETA_MAX, pT < FJVT_PT_MAX
    float pt = 0.f;            // for the working point's own forward pT cap (JvtWP::fwdPtMax)
  };

  // The removal rule, in one place. The JVT and fJVT windows are disjoint in
  // |eta|; with wp.rptFwd a forward jet can fail both (callers test jvtRemoves
  // first, so it is counted as a JVT removal).
  inline bool jvtRemoves (const JetDisc& d, const JvtWP& wp) {
    if (d.jvtWin) return d.rpt < wp.rptMin;
    if (wp.rptFwd && d.rptFwdWin && d.pt < wp.fwdPtMax) return d.rpt < (wp.rptMinFwd >= 0.0 ? wp.rptMinFwd : wp.rptMin);
    return false;
  }
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
      d.rptFwdWin = aeta >= FJVT_ETA_MIN && aeta < RPT_FWD_ETA_MAX && pt < FJVT_PT_MAX;
      d.pt = (float)pt;
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
