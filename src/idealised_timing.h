#ifndef IDEALISED_TIMING_H
#define IDEALISED_TIMING_H

// ---------------------------------------------------------------------------
// idealised_timing.h
//   The re-smeared "idealised" timing world behind rpt_v5's truth and
//   waves_ideal scenarios, shared between util/rpt_v5_hist.cxx (which fills
//   the histograms) and python/runHGTD_Clustering.cxx (the event-display
//   macro).
//
//   It lives here rather than inside rpt_v5_hist so that a display of either
//   scenario shows the SAME smeared track times the histograms were filled
//   with, not a fresh draw. That works because the RNG is re-seeded per event
//   from event-stable quantities (idealisedTimingSeed) and the draw order is
//   fixed: every timed track in index order, then the vertex time. Changing
//   that order changes every truth / waves_ideal histogram -- and silently
//   desynchronises the displays from them.
//
//   Only tracks HGTD actually timed are smeared: idealising the time of a
//   track the detector never measured would invent coverage it does not
//   have, which rpt_v5's central baseline catches immediately.
//
//   Mirrors getSmearedTrackTime's priority (particle production time, then
//   truth vertex time, then a pileup draw) but draws from the caller's RNG:
//   that helper uses the GLOBAL gRandom, which is a data race under
//   TTreeProcessorMT -- the same class of bug as the TColor race.
// ---------------------------------------------------------------------------

#include <TRandom3.h>

#include <cmath>
#include <unordered_map>
#include <vector>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"

namespace MyUtl {

  // Reference-study truth-t0 smearing (util/myJet_ana_fr.C).
  inline constexpr double TRUTH_VTX_SMEAR = 10.0;  // ps, on the HS vertex time
  inline constexpr double TRUTH_TRK_SMEAR = 30.0;  // ps, on each track's own truth time

  struct IdealisedTiming {
    std::unordered_map<int, double> times;  // track index -> smeared time (timed tracks only)
    std::unordered_map<int, double> res;    // track index -> TRUTH_TRK_SMEAR
    double tTruthVtx = 0.0;                 // Gaus(t_HS, TRUTH_VTX_SMEAR)
  };

  // Per-event seed from quantities that do not depend on how TTreeProcessorMT
  // distributes entries, or on which file the event sits in -- so a display
  // reading the event from a picked copy reproduces the same draws.
  inline UInt_t idealisedTimingSeed(BranchPointerWrapper* b) {
    UInt_t sd = (UInt_t)std::llround(std::abs(b->truthVtxZ[0])    * 1e4)
              ^ ((UInt_t)std::llround(std::abs(b->truthVtxTime[0]) * 1e3) << 11)
              ^ ((UInt_t)b->trackZ0.GetSize() << 23);
    return sd ? sd : 1u;
  }

  // Re-seeds rng, then draws every timed track's idealised time in index
  // order, then the idealised vertex time -- in that order, always.
  inline IdealisedTiming makeIdealisedTiming(BranchPointerWrapper* b, TRandom3& rng) {
    rng.SetSeed(idealisedTimingSeed(b));
    IdealisedTiming out;
    for (size_t i = 0; i < b->trackZ0.GetSize(); ++i) {
      if (b->trackTimeValid[i] != 1) continue;
      const int pi = b->trackToParticle[i];
      const int vi = b->trackToTruthvtx[i];
      double tPart;
      if (pi != -1)      tPart = b->particleT[pi];
      else if (vi != -1) tPart = b->truthVtxTime[vi];
      else               tPart = rng.Gaus(b->truthVtxTime[0], PILEUP_SMEAR);
      out.times.emplace((int)i, rng.Gaus(tPart, TRUTH_TRK_SMEAR));
      out.res.emplace((int)i, TRUTH_TRK_SMEAR);
    }
    out.tTruthVtx = rng.Gaus(b->truthVtxTime[0], TRUTH_VTX_SMEAR);
    return out;
  }

  // The cluster structure depends on the track times, so the idealised world
  // gets its own collection rather than reusing one built from real times.
  // values[0] of each cluster is the weighted mean of the SMEARED times.
  inline std::vector<Cluster> idealisedClusters(const std::vector<int>& trk,
                                                BranchPointerWrapper* b,
                                                const IdealisedTiming& it) {
    auto cl = makeSimpleClusters(trk, b, /*useSmearedTimes=*/true, it.times, it.res,
                                 /*checkTimeValid=*/true, /*usez0=*/false);
    doIterativeClustering(&cl, DIST_CUT_CONE);
    return cl;
  }

  // waves_ideal's selection: the cluster closest in time to truth. Ranked on
  // |dt|, not purity -- purity measures where a track CAME FROM, not whether
  // its time is right, so a lone HS track with a mis-assigned HGTD hit forms a
  // 100%-pure cluster at a wrong time and would win. Returns -1 if empty.
  inline int closestToTruthCluster(const std::vector<Cluster>& cl, double tTruth) {
    if (cl.empty()) return -1;
    size_t bi = 0; double bestDt = 1e50;
    for (size_t i = 0; i < cl.size(); ++i) {
      const double dt = std::abs(cl[i].values[0] - tTruth);
      if (dt < bestDt) { bestDt = dt; bi = i; }
    }
    return (int)bi;
  }

}  // namespace MyUtl

#endif  // IDEALISED_TIMING_H
