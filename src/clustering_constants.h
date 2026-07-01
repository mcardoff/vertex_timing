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
  const double PASS_SIGMA         = 60.0;  // Pass threshold for efficiency (ps)
  const double PILEUP_SMEAR       = 175.0; // Pileup track resolution

  const double MAX_TRK_VTX_SIG    = 3.0;   // Pileup removal sigma
  const double MAX_NSIGMA         = 3.0;   // how close a track can be to PV
  const double DIST_CUT_CONE      = 3.0;   // Distance cut for cone clustering
  const double DIST_CUT_SIMUL     = 3.0;   // Distance cut for simul. clustering
  const double DIST_CUT_ITER      = 3.0;   // Distance cut for iterative clustering
  const double DIST_CUT_T_REFINED = 2.0;   // Re-clustering distance cut for JET_T_REFINED (WAVES_RECLUST)
  const double WAVES_DR_FLOOR    = 0.05;  // Minimum ΔR for WAVeS 1/ΔR weight (prevents divergence)
  const int    CONE_ITER_K        = 3;     // Top cone clusters to refine; used only by util/scratch/cone_iter_k_sweep.cxx
  const double TRUTH_PULL_CUT     = 3.0;   // |pull| < cut keeps track as truth-matched
  // z₀ pull-width inflation factor (measured from z0_pull_diag: σ ≈ 1.15 across all
  // η bins, indicating the covariance matrix underestimates the true z₀ resolution by
  // ~15%).  Inflating var_z0 by this factor before taking the square root matches the
  // observed spread, analogous to the 1.5² = 2.25 inflation used for the timing gate.
  const double Z0_VAR_INFLATION   = 1.15 * 1.15;  // = 1.3225
  // ITERATIVE_SPLIT tunables (inline so diagnostics can sweep them).
  inline double T_PULL_SPLIT_THRESHOLD = 1.5;  // split clusters whose t-pull RMS exceeds this
  inline double DIST_CUT_SPLIT         = 2.0;  // tighter cut used when re-clustering after split
  // Per-track timing resolution used for Ideal Resolution/Efficiency scenarios.
  // Flat per-track value (independent of hit count), representing a hypothetically
  // better detector.  Contrast with real HGTD: ~30 ps/hit → 30/√nHits ≈ 15–21 ps/track.
  const double IDEAL_TRACK_RES   = 1.0;  // ps, flat per-track

  // ---------------------------------------------------------------------------
  // 3b. Clustering method selector
  //   Passed as a single parameter to clusterTracksInTime, replacing the old
  //   bool useCone flag.  Four modes are supported:
  //     SIMULTANEOUS    — globally-closest-pair agglomerative merge
  //     CONE            — seed-and-cone, absorbs all in-cone candidates at once
  //     ITERATIVE       — nearest-neighbour iterative, centroid updated per step
  //     ITERATIVE_SPLIT — ITERATIVE then split clusters with high t-pull RMS
  // ---------------------------------------------------------------------------
  enum class ClusteringMethod {
    SIMULTANEOUS, // doSimultaneousClustering — agglomerative minimum-distance
    CONE,         // doConeClustering — seed-and-cone simultaneous absorption
    ITERATIVE,    // doIterativeClustering — nearest-neighbour, centroid-updating
    ITERATIVE_SPLIT, // ITERATIVE + post-process: split clusters with high t-pull RMS
  };

  // ---------------------------------------------------------------------------
  // 3c. Track pre-filter selector
  //   Used by CollSpec to describe which subset of forward tracks a dedicated
  //   collection is built from.  The actual filtering functions are called at
  //   runtime inside processEventData via a switch on this enum.
  //     ALL      — full forward track list (no pre-filtering)
  //     JET      — tracks falling inside a forward jet cone (FILTJET)
  //     HS_ONLY  — truth-HS-linked tracks only (TEST_HS)
  // ---------------------------------------------------------------------------
  enum class TrackFilterType { ALL, JET, HS_ONLY };

  // ---------------------------------------------------------------------------
  // 4. Histogram axis ranges and fold values
  //   xMin/xMax define the full histogram axis.  FOLD_* values mark the
  //   point at which overflow is collapsed into the last visible bin so that
  //   sparse high-multiplicity tails don't dominate the plots.
  // ---------------------------------------------------------------------------
  const double DIFF_MIN = -1000.0, DIFF_MAX = 1000.0;
  const double DIFF_WIDTH = 5.0;

  const double PULL_MIN = -50.0, PULL_MAX = 50.0;
  const double PULL_WIDTH = 0.25;   // 400 bins booked; display range clipped to ±10 at plot time

  const double PURITY_MIN = 0, PURITY_MAX = 1;
  const double PURITY_WIDTH = 0.05;

  const double FJET_MIN = -0.5, FJET_MAX = 31.5, FOLD_FJET = 5;
  const double FJET_WIDTH = 1.0;

  const double TRACK_MIN = 2.5, TRACK_MAX = 100.5, FOLD_TRACK = 50;
  const double TRACK_WIDTH = 2.0;

  const double PU_TRACK_MIN = 2.5, PU_TRACK_MAX = 100.5, FOLD_HS_TRACK = 20;  
  const double PU_TRACK_WIDTH = 1.0;

  const double HS_TRACK_MIN = 2.5, HS_TRACK_MAX = 100.5, FOLD_PU_TRACK = 20;
  const double HS_TRACK_WIDTH = 1.0;

  const double PU_FRAC_WIDTH = 0.1;
  const double PU_FRAC_MIN = 0, PU_FRAC_MAX = 1.0 + PU_FRAC_WIDTH, FOLD_PU_FRAC = 1.0;

  // Avg nHGTD hits per track in the selected cluster (range 1–4; no folding needed)
  const double NHIT_MIN = 0.75, NHIT_MAX = 4.25;
  const double NHIT_WIDTH = 0.5;

  // Cluster-level PU fraction by track count; mirrors event PU_FRAC_* for direct comparison
  const double CLUS_PU_FRAC_WIDTH = 0.1;
  const double CLUS_PU_FRAC_MIN = 0.0, CLUS_PU_FRAC_MAX = 1.0 + CLUS_PU_FRAC_WIDTH, FOLD_CLUS_PU_FRAC = 1.0;

  // σ_t factor used as the third multiplicative term in cluster quality.
  // Linear roll-off: factor=1 below FLOOR, factor=0 above CEIL, linear between.
  // Floor=5 ps represents a well-determined multi-track multi-hit cluster;
  // Ceil=20 ps marks the boundary where the cluster's own time estimate is unreliable.
  const double CLUS_SIGMA_T_FACTOR_FLOOR = 5.0;
  const double CLUS_SIGMA_T_FACTOR_CEIL  = 20.0;

  const double Z_MIN = -200, Z_MAX = 200, FOLD_Z = 100;
  const double Z_WIDTH = 10.0;

  const double EFF_YMIN = 0.7, EFF_YMAX = 1.08;
  const double PUR_YMIN = 0.5, PUR_YMAX = 1.0;
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
  //     requiresPurity   — true: fills gated on cluster/event purity (TEST_MISAS, WAVES_MISCL, WAVES_MISAS)
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
    static const Score CONE;
    static const Score FILTJET;
    static const Score HGTD_SORT;
    static const Score CONE_BDT;
    static const Score TEST_MISAS;
    static const Score TEST_HS;
    static const Score WAVES;
    static const Score JET_T_REFINED;
    static const Score WAVES_MISCL;
    static const Score WAVES_MISAS;
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
  inline const Score Score::CONE_BDT   = {
    8, "Cone (BDT)", "CONE_BDT",
    true , false,  0.3f
  };
  inline const Score Score::HGTD_SORT  = {
    10, "HGTD BDT (pT-sorted)", "HGTD_SORT",
    true , false,  0.3f
  };

  inline const Score Score::TEST_MISAS = { 13, STR_TRKPTZ + " [Events with Perfect Timing]"  , "TRKPTZ Perf. Time",    false, true, -1.f };
  // WAVES: WAVeS-style selection score — Σ pT·pT_jet/max(ΔR,floor) × exp(−1.5|Δz|).
  // Pure selection: picks the highest-scoring main-collection cluster and reports its
  // standard weighted-mean time and full-cluster purity (no in-jet-only recomputation).
  inline const Score Score::WAVES  = { 18, "WAVeS Score",                      "WAVES",        false, false, -1.f };
  // JET_T_REFINED: clusters only jet-proximate tracks at DIST_CUT_T_REFINED (2σ iterative)
  inline const Score Score::JET_T_REFINED = { 19, "WAVeS 2#sigma t Re-clustering",  "WAVES_RECLUST", true,  false, -1.f,
                                              DIST_CUT_T_REFINED, ClusteringMethod::ITERATIVE,
                                              false, TrackFilterType::JET };
  // WAVeS oracle variants: selection by the WAVeS score, denominator gates applied
  // at fill time — cluster purity (MISCL-style) / HS timing purity (like TEST_MISAS)
  inline const Score Score::WAVES_MISCL  = { 20, "WAVeS [Events with Pure Clusters]" , "WAVES Pure Clust.", false, true, -1.f };
  inline const Score Score::WAVES_MISAS  = { 21, "WAVeS [Events with Perfect Timing]", "WAVES Perf. Time" , false, true, -1.f };

  // Scores with a dedicated collection (distCut ≥ 0 → buildsCollection() = true)
  inline const Score Score::CONE       = {  7, "Cone"                       , "CONE",     true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::CONE      };
  inline const Score Score::FILTJET    = {  9, "Filter Tracks in Jets"      , "FILTJET",  true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::CONE,      false, TrackFilterType::JET     };
  inline const Score Score::TEST_HS    = { 14, STR_TRKPTZ + " (HS only)"   , "TEST_HS",  true , false, -1.f, DIST_CUT_CONE,      ClusteringMethod::CONE,      false, TrackFilterType::HS_ONLY };

  inline const std::vector<Score> SCORE_REGISTRY = {
    Score::HGTD,      Score::TRKPT,     Score::TRKPTZ,    Score::PASS,
    Score::CONE,
    Score::FILTJET,   Score::HGTD_SORT,
    Score::CONE_BDT,
    Score::TEST_MISAS, Score::TEST_HS,
    Score::WAVES,
    Score::JET_T_REFINED,
    Score::WAVES_MISCL,
    Score::WAVES_MISAS,
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
