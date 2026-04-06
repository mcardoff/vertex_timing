#ifndef CLUSTERING_CONSTANTS_H
#define CLUSTERING_CONSTANTS_H

// ---------------------------------------------------------------------------
// clustering_constants.h
//   Central repository for all compile-time constants, enums, and small
//   utility functions shared across the analysis.  Everything lives inside
//   the MyUtl namespace so that name collisions with ROOT globals are
//   avoided.  Sections:
//     1. Debug flag
//     2. Plot colour palette
//     3. Event / track selection cuts
//     4. Histogram axis ranges and fold values
//     5. Score struct + SCORE_REGISTRY
//     6. FitParamFields enum + string converter
//     7. folded() helper
// ---------------------------------------------------------------------------

#include "clustering_includes.h"

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // 1. Debug flag
  //   Set to true to enable verbose per-event stdout logging throughout the
  //   analysis.  Compiled out entirely when false (no runtime cost).
  // ---------------------------------------------------------------------------
  const bool DEBUG = false  ;
  // ---------------------------------------------------------------------------
  // 2. Plot colour palette
  //   Eleven colours drawn from ROOT's P10 and P6 palettes.  COLORS is the
  //   ordered vector used when cycling through scores or histogram slices.
  // ---------------------------------------------------------------------------
  const Color_t C01 = kP10Blue  ;
  const Color_t C02 = kP10Red   ;
  const Color_t C03 = kP10Yellow;
  const Color_t C04 = kP10Violet;
  const Color_t C05 = kP10Cyan  ;
  const Color_t C06 = kP10Brown ;
  const Color_t C07 = kP10Orange;
  const Color_t C08 = kP10Green ;
  const Color_t C09 = kP10Ash   ;
  const Color_t C10 = kP10Gray  ;
  const Color_t C11 = kP6Blue   ;
  const std::vector<Color_t> COLORS = {C01, C02, C03, C04, C05, C06,
                                       C07, C08, C09, C10, C11};  

  // ---------------------------------------------------------------------------
  // 3. Event / track selection cuts
  //   All kinematic and association thresholds used by the event selection
  //   and track clustering pipeline.  Changing a cut here propagates
  //   automatically to every function that references these constants.
  // ---------------------------------------------------------------------------
  const int    MIN_JETS           = 2;     // min n jets
  const int    MIN_PASSPT_JETS    = 2;     // min n jets >30 GeV
  const int    MIN_PASSETA_JETS   = 1;     // min n forward jets >30GeV
  const int    MIN_CLUSTER_TRACKS = 0;     // min tracks required to select a cluster
  const int    MIN_NHS_TRACK      = 2;     // testing only
  const int    MAX_NHS_TRACK      = 6;     // testing only
  const double VBS_JET_D_ETA      = 3.0;   // min eta separation for VBS Jets
  const double MIN_JET_PT         = 30.0;  // self explanatory
  const double MAX_VTX_DZ         = 2.0;   // max error for reco HS vertex z
  const double MIN_HGTD_ETA       = 2.38;  // HGTD Min eta
  const double MAX_HGTD_ETA       = 4.0;   // HGTD Max eta
  const double MIN_ABS_ETA_JET    = 2.38;  // "forward" jet min eta
  const double MAX_ABS_ETA_JET    = 4.00;  // "forward" jet max eta
  const double MIN_ABS_ETA_TRACK  = 2.38;  // "forward" track min eta 
  const double MAX_ABS_ETA_TRACK  = 4.00;  // "forward" track max eta
  const double MIN_TRACK_PT       = 1.0;   // clustered track_pt > 1.0 GeV
  const double MAX_TRACK_PT       = 30.0;  // clustered track_pt < 30.0 GeV
  const double MIN_TRACK_PT_COUNT = 1.0;   // track_pt > 1.0 GeV for histgramming
  const double PASS_SIGMA         = 100.0; // Pass threshold for efficiency (ps)
  const double PILEUP_SMEAR       = 175.0; // Pileup track resolution

  const double MAX_TRK_VTX_SIG    = 3.0;   // Pileup removal sigma
  const double MAX_NSIGMA         = 3.0;   // how close a track can be to PV
  const double DIST_CUT_CONE      = 3.0;   // Distance cut for cone clustering
  const double DIST_CUT_SIMUL     = 3.0;   // Distance cut for simul. clustering
  const double DIST_CUT_ITER      = 3.0;   // Distance cut for iterative clustering
  const double DIST_CUT_REFINE    = 2.0;   // Timing distance cut used by ZT_REFINED (legacy)
  const double DIST_CUT_T_REFINED = 2.0;   // Re-clustering distance cut for T_REFINED and ZT_REFINED
  const double TVA_CUT_Z_REFINED  = 2.0;   // z₀ TVA significance cut for Z_REFINED and ZT_REFINED
  const int    CONE_ITER_K        = 3;     // Top cone clusters to refine (REFINED score)
  const double TRUTH_PULL_CUT     = 2.0;   // |pull| < cut keeps track as truth-matched
  const double HS_TIMING_QUALITY_CUT = 3.0;   // |pull| < cut → HS track timing is good (TEST_MISAS gate)
  // Per-track timing resolution used for Ideal Resolution/Efficiency scenarios.
  // Flat per-track value (independent of hit count), representing a hypothetically
  // better detector.  Contrast with real HGTD: ~30 ps/hit → 30/√nHits ≈ 15–21 ps/track.
  const double IDEAL_TRACK_RES   = 1.0;  // ps, flat per-track

  // ---------------------------------------------------------------------------
  // 3b. Clustering method selector
  //   Passed as a single parameter to clusterTracksInTime, replacing the old
  //   bool useCone flag.  Three modes are supported:
  //     SIMULTANEOUS — globally-closest-pair agglomerative merge
  //     CONE         — seed-and-cone, absorbs all in-cone candidates at once
  //     ITERATIVE    — nearest-neighbour iterative, centroid updated per step
  // ---------------------------------------------------------------------------
  enum class ClusteringMethod {
    SIMULTANEOUS, // doSimultaneousClustering — agglomerative minimum-distance
    CONE,         // doConeClustering — seed-and-cone simultaneous absorption
    ITERATIVE,    // doIterativeClustering — nearest-neighbour, centroid-updating
  };

  // ---------------------------------------------------------------------------
  // 3c. Track pre-filter selector
  //   Used by CollSpec to describe which subset of forward tracks a dedicated
  //   collection is built from.  The actual filtering functions are called at
  //   runtime inside processEventData via a switch on this enum.
  //     ALL      — full forward track list (no pre-filtering)
  //     JET      — tracks falling inside a forward jet cone (FILTJET)
  //     Z0_TVA   — tracks passing z₀ TVA at TVA_CUT_Z_REFINED (ZT_REFINED)
  //     HS_ONLY  — truth-HS-linked tracks only (TEST_HS)
  // ---------------------------------------------------------------------------
  enum class TrackFilterType { ALL, JET, Z0_TVA, HS_ONLY };

  // ---------------------------------------------------------------------------
  // 4. Histogram axis ranges and fold values
  //   xMin/xMax define the full histogram axis.  FOLD_* values mark the
  //   point at which overflow is collapsed into the last visible bin so that
  //   sparse high-multiplicity tails don't dominate the plots.
  // ---------------------------------------------------------------------------
  const double DIFF_MIN = -1000.0, DIFF_MAX = 1000.0;
  const double DIFF_WIDTH = 4.0;

  const double PURITY_MIN = 0, PURITY_MAX = 1;
  const double PURITY_WIDTH = 0.05;

  const double FJET_MIN = -0.5, FJET_MAX = 31.5, FOLD_FJET = 5;
  const double FJET_WIDTH = 1.0;

  const double VTX_DZ_MIN = 0, VTX_DZ_MAX = 5.0, FOLD_VTX_DZ = 2.0;
  const double VTX_DZ_WIDTH = 0.1;

  const double TRACK_MIN = -0.5, TRACK_MAX = 100.5, FOLD_TRACK = 20;
  const double TRACK_WIDTH = 1.0;

  const double PU_TRACK_MIN = TRACK_MIN, PU_TRACK_MAX = TRACK_MAX, FOLD_HS_TRACK = FOLD_TRACK;  
  const double PU_TRACK_WIDTH = TRACK_WIDTH;

  const double HS_TRACK_MIN = TRACK_MIN, HS_TRACK_MAX = TRACK_MAX, FOLD_PU_TRACK = FOLD_TRACK;
  const double HS_TRACK_WIDTH = TRACK_WIDTH;

  const double PU_FRAC_WIDTH = 0.1;
  const double PU_FRAC_MIN = 0, PU_FRAC_MAX = 1.0 + PU_FRAC_WIDTH, FOLD_PU_FRAC = 1.0;

  const double Z_MIN = -200, Z_MAX = 200, FOLD_Z = 100;
  const double Z_WIDTH = 10.0;

  const double EFF_YMIN = 0.7, EFF_YMAX = 1.08;
  const double PUR_YMIN = 0.0, PUR_YMAX = 1.5;
  const double RES_YMIN = 0.0, RES_YMAX = 40.0;
  const double BKG_RES_YMIN = 90.0, BKG_RES_YMAX = 500.0;

  // ---------------------------------------------------------------------------
  // 5. Score struct
  //   Single source of truth for every cluster-selection algorithm.
  //   Each named instance carries display metadata, behavioural flags, and —
  //   for scores with a dedicated cluster collection — the full clustering
  //   spec so that processEventData can drive everything from SCORE_REGISTRY
  //   without a separate CollSpec table.
  //
  //   Display / behaviour fields:
  //     id               — stable integer identity (map key, COLORS index)
  //     longName         — ROOT-LaTeX label for plot titles
  //     shortName        — compact identifier for histogram names / tables
  //     usesOwnCollection— true: skipped by the all-scores chooseCluster
  //                        overload; collection built in section E (when
  //                        distCut ≥ 0) or directly in selectClusters
  //                        (HGTD, HGTD_SORT, CONE_BDT where distCut = -1).
  //     requiresPurity   — true: fills gated on cluster purity (TEST_MISCL)
  //     threshold        — score gate for passEfficiency; -1 = no gate
  //
  //   Collection spec fields (only used when usesOwnCollection && distCut≥0):
  //     distCut  — Mahalanobis distance cut for clusterTracksInTime
  //     method   — clustering algorithm (CONE or ITERATIVE)
  //     useZ0    — true: 2D (z₀,t) metric; false: 1D time-only metric
  //     filter   — track pre-filter applied before clustering
  //
  //   Adding a new dedicated-collection score: one declare + one define
  //   in this file.  No other file needs touching.
  // ---------------------------------------------------------------------------
  struct Score {
    int id;
    std::string longName;
    const char* shortName;
    bool             usesOwnCollection = false;
    bool             requiresPurity    = false;
    float            threshold         = -1.f;
    // Collection spec — meaningful only when usesOwnCollection && distCut >= 0
    double           distCut           = -1.0;
    ClusteringMethod method            = ClusteringMethod::ITERATIVE;
    bool             useZ0             = false;
    TrackFilterType  filter            = TrackFilterType::ALL;

    Score() = default;
    Score(int id_, std::string ln, const char* sn,
          bool own=false, bool pur=false, float thr=-1.f,
          double dc=-1.0, ClusteringMethod m=ClusteringMethod::ITERATIVE,
          bool z0=false, TrackFilterType f=TrackFilterType::ALL)
      : id(id_), longName(std::move(ln)), shortName(sn),
        usesOwnCollection(own), requiresPurity(pur), threshold(thr),
        distCut(dc), method(m), useZ0(z0), filter(f) {}

    bool operator<(const Score& o)  const { return id < o.id; }
    bool operator==(const Score& o) const { return id == o.id; }
    bool operator!=(const Score& o) const { return id != o.id; }
    const char* toString()      const { return longName.c_str(); }
    const char* toStringShort() const { return shortName; }
    bool hasThreshold()         const { return threshold >= 0.f; }
    // True when this score's collection is built via SCORE_REGISTRY in section E
    bool buildsCollection()     const { return usesOwnCollection && distCut >= 0.0; }

    static const Score HGTD;
    static const Score TRKPT;
    static const Score TRKPTZ;
    static const Score PASS;
    static const Score T_REFINED;
    static const Score Z_REFINED;
    static const Score ZT_REFINED;
    static const Score CONE;
    static const Score FILTJET;
    static const Score HGTD_SORT;
    static const Score TEST_ML;
    static const Score TEST_MISCL;
    static const Score CONE_BDT;
    static const Score TEST_MISAS;
    static const Score TEST_HS;
    static const Score ZT_ITER;
    static const Score PERF_EVT;
  };

  inline const std::string STR_TRKPTZ = "#Sigma p_{T}e^{-|#Delta z|}";

  // Columns: id  longName  shortName  own    pur    thr     distCut           method                        useZ0  filter
  //          ──  ────────  ─────────  ─────  ─────  ──────  ────────────────  ──────────────────────────    ─────  ──────────────────────
  // Scores without a dedicated collection (distCut omitted → -1, ignored)
  inline const Score Score::HGTD       = {
    0, "HGTD Algorithm", "HGTD",
    true , false, -1.0f
  };
  inline const Score Score::TRKPT      = {
    1, "#Sigma p_{T}", "TRKPT",
    false, false, -1.0f
  };
  inline const Score Score::TRKPTZ     = {  
    2, STR_TRKPTZ + " [Baseline Algorithm]", "TRKPTZ",
    false, false, -1.0f
  };
  inline const Score Score::PASS       = {  
    3, "Pass Cluster", "PASS",
    false, false, -1.0f
  };
  inline const Score Score::Z_REFINED  = {  
    5, "2#sigma z_{0} Refinement","Z_REFINED",
    false, false, -1.0f
  };
  inline const Score Score::CONE_BDT   = {  
    8, "Cone (BDT)", "CONE_BDT",
    true , false,  0.3f
  };
  inline const Score Score::HGTD_SORT  = { 
    10, "HGTD BDT (pT-sorted)", "HGTD_SORT",
    true , false,  0.3f
  };
  inline const Score Score::TEST_ML    = { 
    11, "DNN Selection", "TEST_ML",
    false, false,  0.3f
  };

  inline const Score Score::TEST_MISCL = { 12, STR_TRKPTZ + " [Events with Pure Clusters]"   , "TRKPTZ Pure Clust.",    false, true, -1.f };
  inline const Score Score::TEST_MISAS = { 13, STR_TRKPTZ + " [Events with Perfect Timing]"  , "TRKPTZ Perf. Time",    false, true, -1.f };
  inline const Score Score::PERF_EVT   = { 17, STR_TRKPTZ + " [Pure Clusters + Perf. Timing]", "TRKPTZ Pure Clust./Perf. Time", false, true, -1.f };

  // Scores with a dedicated collection (distCut ≥ 0 → buildsCollection() = true)
  inline const Score Score::CONE       = {  7, "Cone"                       , "CONE",     true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::CONE      };
  inline const Score Score::FILTJET    = {  9, "Filter Tracks in Jets"      , "FILTJET",  true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::CONE,      false, TrackFilterType::JET     };
  inline const Score Score::TEST_HS    = { 14, STR_TRKPTZ + " (HS only)"   , "TEST_HS",  true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::CONE,      false, TrackFilterType::HS_ONLY };
  inline const Score Score::T_REFINED  = {  4, "2#sigma t Re-clustering"    , "T_REFINED",true , false, -1.f, DIST_CUT_T_REFINED, ClusteringMethod::ITERATIVE };
  inline const Score Score::ZT_REFINED = {  6, "2#sigma z_{0}+t Re-clustering", "ZT_REFINED",true, false, -1.f, DIST_CUT_T_REFINED, ClusteringMethod::ITERATIVE, false, TrackFilterType::Z0_TVA  };
  inline const Score Score::ZT_ITER    = { 16, "2D (z_{0},t) Iterative"     , "ZT_ITER",  true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::ITERATIVE, true  };

  inline const std::vector<Score> SCORE_REGISTRY = {
    Score::HGTD,      Score::TRKPT,     Score::TRKPTZ,    Score::PASS,
    Score::T_REFINED, Score::Z_REFINED, Score::ZT_REFINED, Score::CONE,
    Score::FILTJET,   Score::HGTD_SORT,
    Score::TEST_ML,   Score::TEST_MISCL, Score::CONE_BDT,
    Score::TEST_MISAS, Score::TEST_HS,
    Score::ZT_ITER,
    Score::PERF_EVT,
  };

  // ---------------------------------------------------------------------------
  // 6. FitParamFields enum + string converter
  //   Indexes the individual parameters extracted from a double-Gaussian fit:
  //   MEAN, core SIGMA, background BSIGMA, amplitudes CORE / BKG, and their
  //   RATIO.  Used by FitParams::fromEnum() to select which distribution to
  //   draw or fill.
  // ---------------------------------------------------------------------------
  enum FitParamFields {
    MEAN = 0, SIGMA = 1, CORE = 2, BKG = 3, RATIO = 4, BSIGMA = 5, RMS = 6
  };
  const std::vector<FitParamFields> FITPARAM_VEC = {
      FitParamFields::SIGMA, FitParamFields::BSIGMA, FitParamFields::CORE,
      FitParamFields::BKG,   //FitParamFields::RATIO,  FitParamFields::MEAN,
      FitParamFields::RMS
  };  

  auto toString(
    FitParamFields key
  ) -> const char* {
    switch (key) {
    case FitParamFields::SIGMA:  return "SIGMA";
    case FitParamFields::BSIGMA: return "BSIGMA";
    case FitParamFields::CORE:   return "CORE";
    case FitParamFields::BKG:    return "BKG";
    case FitParamFields::RATIO:  return "RATIO";
    case FitParamFields::MEAN:   return "MEAN";
    default:                     return "INVALID";
    }
  }

  // ---------------------------------------------------------------------------
  // 7. folded() helper
  //   Clamps raw to fold when raw >= fold, collapsing overflow into the last
  //   visible histogram bin.  Used when filling efficiency/resolution plots
  //   to avoid sparse tails at high multiplicity.
  // ---------------------------------------------------------------------------
  template <typename T>
  auto folded(T raw, T fold) -> T {
    return (raw >= fold) ? fold : raw;
  }
}
#endif // CLUSTERING_CONSTANTS_H
