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
  // VBS candidate-pair |Deta| requirement. RUNTIME-SETTABLE (inline, not const)
  // via --vbs-deta=<x>; see resolveSelection() in sample_config.h.
  //
  // Why it is a knob and not a constant: this is a VBS *topology* requirement, and
  // Z+jets has no VBS topology. Applied to Z+jets it removes 87% of everything
  // surviving the lepton selection (520,126 -> 67,184 out of 3M), and it keeps the
  // events that happen to present two >30 GeV reco jets in opposite hemispheres at
  // |Deta| >= 3 with NO truth-HS matching -- which at this pileup level
  // preferentially selects events where that pair involves a PILEUP jet, since a
  // genuine Z+jets event rarely supplies one. That is precisely the pathology the
  // timing study is trying to measure, so it has to be possible to run with and
  // without the cut rather than assume it is neutral.
  //
  // DEFAULT IS NOW 0 (no requirement), changed 2026-08-26. The cut was never
  // shown to help: vbs_mjj_diag measured it purity-NEUTRAL on VBF (57.1% ->
  // 57.2%) and purity-NEGATIVE on Z+jets (2.69% -> 2.51%) and dijet (40.1% ->
  // 35.0%), i.e. on two samples out of three it actively admits more
  // pileup-legged candidate pairs than it removes. m_jj does the same job
  // better everywhere (see VBS_JET_MJJ below), so the topology requirement is
  // now carried by mass alone.
  //
  // 0 disables the |Deta| test only. passJetPtCut still requires the jet-pT
  // and forward-jet counts, so this is "still a forward-jet event, not forced
  // into a VBS topology" -- not "no topology selection at all".
  const  double VBS_JET_D_ETA_DEFAULT = 0.0;
  inline double VBS_JET_D_ETA         = VBS_JET_D_ETA_DEFAULT;
  // VBS candidate-pair m_jj requirement. RUNTIME-SETTABLE (inline, not const)
  // via --vbs-mjj=<x>; see resolveSelection() in sample_config.h.
  //
  // Measured against |Deta| by vbs_mjj_diag
  // (util/vbs_mjj_diag.cxx) measured both variables against the same
  // candidate-pair truth-HS-match denominator on all three samples. Findings
  // that motivated making this the default topology requirement in place of
  // |Deta|:
  //   - |Deta| >= 3 (the current default) is purity-neutral on VBF (57.1% ->
  //     57.2%) but purity-NEGATIVE on Z+jets (2.69% -> 2.51%) and dijet
  //     (40.1% -> 35.0%) -- it preferentially admits pileup-legged pairs on
  //     both, the pathology documented on VBS_JET_D_ETA above.
  //   - m_jj is at-least-as-good at matched HS efficiency everywhere (VBF
  //     +1.2pp, dijet +8.7pp) but does NOT fix Z+jets (+0.2pp vs the 2.69%
  //     no-cut baseline) -- there is essentially no genuine forward VBS
  //     topology in Z+jets for either kinematic variable to select for
  //     (1,872 genuine pairs out of 2M events processed).
  //   - Neither variable is monotonic in purity: m_jj peaks around 600-800
  //     GeV in VBF/dijet and degrades above ~2 TeV (a far-forward pileup jet
  //     manufactures a large mass); |Deta| degrades past ~4-5 similarly.
  // In short: swapping in m_jj is a real but sample-dependent improvement,
  // not a fix for the underlying issue that neither reco-jet kinematic
  // variable enforces the candidate pair actually being hard-scatter -- it is
  // the better of two imperfect proxies, which is why it is the default and
  // why it stays tunable. See
  // isJetTruthHS (clustering_structs.h) for the one thing that does, at
  // truth level.
  //
  // DEFAULT IS NOW 200 GeV, changed 2026-08-26 -- m_jj replaces |Deta| as the
  // topology requirement, for the reasons tabulated above.
  //
  // 200 is a deliberately CONSERVATIVE baseline, not the purity optimum: the
  // same diagnostic puts the peak nearer 600-800 GeV on VBF/dijet. It is set
  // low so the baseline costs little efficiency, and left tunable via
  // --vbs-mjj= so a study that wants the purer working point can ask for it
  // without a rebuild. Anything above ~2 TeV degrades again (a far-forward
  // pileup jet manufactures a large mass).
  const  double VBS_JET_MJJ_DEFAULT = 200.0;
  inline double VBS_JET_MJJ          = VBS_JET_MJJ_DEFAULT;
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
  const double LEPTON_JET_DR      = 0.5;   // ΔR(lepton, jet) overlap-removal cone (Z+jets only)
  const double LEPTON_MIN_PT       = 20.0;  // GeV — min pT for a selected lepton (Z+jets only)
  // Lepton–jet overlap removal switch. Runtime (not compile-time) so it can be
  // driven from the resolved --sample: set true only for Z+jets in each main()
  // right after resolveSample() (mirrors SAMPLE_NAME/OUTPUT_DIR). When true,
  // BranchPointerWrapper binds Track_leptonID and every jet loop skips reco jets
  // within LEPTON_JET_DR of a flagged lepton track; when false (vbf/dijet/local)
  // the branch is never bound and the removal is a no-op, so behaviour is
  // bit-for-bit unchanged on those samples.
  inline bool  OVERLAP_REMOVAL   = false;
  // Overlap-removal behavior when a Z+jets event has a lepton–jet overlap:
  //   REMOVE_JETS — drop only the overlapping reco jet(s), keep the event.
  //   SKIP_EVENT  — veto the whole event if any lepton–jet overlap exists.
  // Runtime so both can be produced without a recompile. Only consulted when
  // OVERLAP_REMOVAL is true (i.e. Z+jets); a no-op on vbf/dijet/local.
  enum class OverlapMode { REMOVE_JETS, SKIP_EVENT };
  inline OverlapMode OVERLAP_MODE = OverlapMode::REMOVE_JETS;
  const int    CONE_ITER_K        = 3;     // Top cone clusters to refine; used only by util/scratch/cone_iter_k_sweep.cxx
  const double TRUTH_PULL_CUT     = 3.0;   // |pull| < cut keeps track as truth-matched
  // z₀ pull-width inflation factor (measured from z0_pull_diag: σ ≈ 1.15 across all
  // η bins, indicating the covariance matrix underestimates the true z₀ resolution by
  // ~15%).  Inflating var_z0 by this factor before taking the square root matches the
  // observed spread, analogous to the 1.5² = 2.25 inflation used for the timing gate.
  const double Z0_VAR_INFLATION   = 1.15 * 1.15;  // = 1.3225

  // ── Alternative track-to-vertex association: the R_pT parameterization ────
  // getNewDzpara is an empirical parameterization of the observed |z0 - z_vtx|
  // spread (mm) -- a 6th-order polynomial in |eta| with pT-binned coefficients --
  // ported verbatim from the reference R_pT study (util/myJet_ana_fr.C) via
  // util/rpt_v5_hist.cxx, where it already governs FORWARD track-to-jet-vertex
  // association.
  //
  // Using it here is domain-safe in a way it is NOT in rpt_v5's central region:
  // getAssociatedTracks cuts on MIN_HGTD_ETA < |eta| < MAX_HGTD_ETA before
  // anything else, so the clustering only ever sees forward tracks -- exactly
  // the range this fit was derived in. (Centrally it returns sigma_z0 = 31 um at
  // eta = 0, several times tighter than ITk truly delivers; that is why rpt_v5
  // uses a plain z0-significance cut there instead.)
  //
  // It also removes an inconsistency: today the tracks that BUILD t0 and the
  // tracks that COUNT toward R_pT are associated by two different rules.
  //
  // Measured on local VBF forward tracks (pT/eta/quality-selected): the
  // parameterization keeps ~9% fewer tracks than significance < 2.5 (194k vs
  // 214k) for a HS purity of 27.8% vs 25.4%, at a HS recall of 74.9% vs 75.6% --
  // i.e. it sheds pileup at nearly constant hard-scatter efficiency.
  //
  // The pT bins are half-open on the low side and the <=1.5 bin catches
  // everything below, matching the original's if-chain exactly (no else-if, so a
  // value landing on a boundary takes the LAST matching branch). Preserved
  // deliberately rather than "cleaned up" -- changing it would silently shift
  // which coefficients boundary-pT tracks get.
  inline double getNewDzpara(double eta, double pt) {
    eta = std::abs(eta);
    double p[7] = {0, 0, 0, 0, 0, 0, 0};
    auto set = [&p](const double (&src)[7]) { std::copy(src, src + 7, p); };

    if (pt <= 1.5) {
      static const double c[7] = { 0.0314036, 0.790955, -2.65987, 3.62073, -2.18228, 0.614866, -0.0634521 };
      set(c);
    }
    if (pt > 1.5 && pt <= 2.5) {
      static const double c[7] = { 0.0229273, 0.540101, -1.80727, 2.45187, -1.47382, 0.414345, -0.0426769 };
      set(c);
    }
    if (pt > 2.5 && pt <= 5.0) {
      static const double c[7] = { 0.0163773, 0.345112, -1.14474, 1.54382, -0.923523, 0.258617, -0.0265446 };
      set(c);
    }
    if (pt > 5.0 && pt <= 10.0) {
      static const double c[7] = { 0.010919, 0.179329, -0.581971, 0.773186, -0.45679, 0.126608, -0.012875 };
      set(c);
    }
    if (pt > 10.0) {
      static const double c[7] = { 0.00835945, 0.0957783, -0.299255, 0.38722, -0.22351, 0.0607521, -0.00606524 };
      set(c);
    }

    double d = p[0];
    double e = 1.0;
    for (int k = 1; k < 7; ++k) { e *= eta; d += p[k] * e; }
    return std::abs(d);
  }

  // Multiple of getNewDzpara a track's |z0 - z_vtx| must stay within. 1.4 is the
  // reference study's working point, and the value rpt_v5_hist's DZ0_PARA_SCALE
  // already runs at -- kept identical so the two associations agree exactly.
  const double DZ0_PARA_SCALE_CLUSTER = 1.4;

  // Selects which rule passTrackVertexAssociation applies. false (default) is
  // the z0-significance cut this analysis has always used; true switches to the
  // parameterization above, in which case the caller's significanceCut argument
  // is ignored (DZ0_PARA_SCALE_CLUSTER takes its place) -- so the 3.0-then-
  // MAX_NSIGMA two-step in processEventData collapses to a single selection.
  //
  // Runtime rather than compile-time, mirroring OVERLAP_REMOVAL above: set once
  // in main() from --dzpara before any worker thread starts, read-only
  // afterwards, so it is safe under TTreeProcessorMT. Default false means every
  // executable that does not parse the flag is bit-for-bit unchanged.
  inline bool USE_DZ_PARA = false;

  // ── AssocRule: one track-to-vertex association rule, made explicit ────────
  // USE_DZ_PARA / MAX_NSIGMA above are process-wide globals: a run picks one
  // rule and every event uses it. That is fine for production, but it makes
  // comparing rules a matter of re-running the whole event loop once per rule
  // -- 52 GB of local VBF ntuple per pass -- and the passes then disagree on
  // which EVENTS they saw as well as which tracks, so a difference cannot be
  // attributed to the association alone.
  //
  // AssocRule carries the same two pieces of information (which rule, what
  // cut) as an explicit value instead, so several rules can be evaluated
  // side by side against the SAME event in a single loop. See
  // util/assoc_study_hist.cxx, the only current consumer.
  //
  // The globals stay the default everywhere: the AssocRule overloads of
  // passTrackVertexAssociation / getAssociatedTracks and the optional
  // `assoc` argument to processEventData are all additive, so any executable
  // that does not pass one is bit-for-bit unchanged.
  struct AssocRule {
    enum class Kind {
      SIGNIFICANCE,  // |z0 - z_vtx| / sqrt(Z0_VAR_INFLATION * var_z0) < cut
      DZ_PARA        // |z0 - z_vtx| / getNewDzpara(eta, pT)           < cut
    };
    Kind        kind = Kind::SIGNIFICANCE;
    double      cut  = MAX_NSIGMA;
    std::string tag;    // filename-safe, e.g. "signif25" / "dzpara14"
    std::string label;  // human/TLatex label, e.g. "z_{0} signif. < 2.5"

    // Union with ghost association: when true the rule accepts a track if it
    // passes the z-test above OR is ghost-associated to a qualifying forward
    // jet. The premise is that a track can carry hard-scatter information the
    // z-cut misses -- so widen the list rather than tighten it.
    //
    // Unlike `kind`/`cut` this is a property of the SET, not of a single
    // track: deciding it needs the event's jets, so it is applied by
    // getAssociatedTracks (which builds the ghost set once per event) and NOT
    // by passTrackVertexAssociation, which stays a pure per-track z-test.
    // Anything reproducing a rule's track list must therefore go through
    // getAssociatedTracks rather than calling the per-track function in a loop.
    bool        orGhost = false;
  };

  // Counting-scan association, held FIXED across every AssocRule under study.
  // processEventData already separates the two uses of the track list (see its
  // step B): the counting scan feeds EventCounts -- n forward tracks, n HS, n
  // PU, PU fraction -- which become the x-AXES of every efficiency/resolution
  // plot, while the clustering list is what t0 is actually built from. Letting
  // the association move the x-axis as well as the algorithm would mean the
  // "n_hs_tracks = 3" bin holds a different set of events for every rule, and
  // no per-bin difference could be attributed to the association. Pinning the
  // counting scan here keeps the binning common so the comparison is paired.
  const double COUNTING_NSIGMA = 3.0;
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
  // 3d. VBS topology region
  //   The two topologies where forward timing is the deciding information.
  //   Classified from the VBS candidate pair by
  //   BranchPointerWrapper::classifyVbsRegion (clustering_structs.h); lives
  //   here rather than nested in that class so Score can carry one as a
  //   denominator gate (see regionGate below).
  //     R1  both candidate legs forward, one truth-HS + one truth-PU --
  //         timing must say WHICH forward jet is the hard scatter.
  //     R2  forward truth-PU leg + central truth-HS leg -- can timing reject
  //         the forward fake?
  //   Shared with rpt_v5's per-jet RpT regions so "R1" cannot come to mean two
  //   different event sets in the two analyses.
  // ---------------------------------------------------------------------------
  enum class VbsRegion { NONE, R1, R2 };

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
    // Denominator gate on VBS topology region (VBF_R1 / VBF_R2). NONE = no
    // region requirement. Orthogonal to requiresPurity: both defer the
    // denominator fill to the gated block in processEventData's step H, they
    // just test different things.
    VbsRegion        regionGate        = VbsRegion::NONE;

    Score() = default;
    Score(int id_, std::string ln, const char* sn,
          bool own=false, bool pur=false, float thr=-1.f,
          double dc=-1.0, ClusteringMethod m=ClusteringMethod::ITERATIVE,
          bool z0=false, TrackFilterType f=TrackFilterType::ALL,
          VbsRegion rg=VbsRegion::NONE)
      : id(id_), longName(std::move(ln)), shortName(sn),
        usesOwnCollection(own), requiresPurity(pur), threshold(thr),
        distCut(dc), method(m), useZ0(z0), filter(f), regionGate(rg) {}

    bool operator<(const Score& o)  const { return id < o.id; }
    bool operator==(const Score& o) const { return id == o.id; }
    bool operator!=(const Score& o) const { return id != o.id; }
    const char* toString()      const { return longName.c_str(); }
    const char* toStringShort() const { return shortName; }
    bool hasThreshold()         const { return threshold >= 0.f; }
    // True when this score restricts its denominator at fill time (step H)
    // rather than counting every selected event.
    bool gatesDenominator()     const { return requiresPurity || regionGate != VbsRegion::NONE; }
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
    static const Score VBF_R1;
    static const Score VBF_R2;
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
  // VBS topology regions. Same construction as the WAVES_MISAS oracle -- WAVeS
  // selection and WAVeS in-jet-refined time, denominator restricted at fill
  // time -- except the gate is the event's VBS region rather than its timing
  // purity. Reading them against the plain WAVES row isolates what the topology
  // does to WAVeS performance, with the algorithm held fixed.
  //
  // WAVeS is the base (rather than the TRKPTZ baseline) to match WAVES_MISAS
  // and because WAVeS is the algorithm under study; swapping the base is a
  // one-line change in Cluster::updateScores plus calculateTime's score list.
  inline const Score Score::VBF_R1 = { 22, "WAVeS [VBF R1: both tags fwd]"    , "VBF_R1", false, false, -1.f,
                                       -1.0, ClusteringMethod::ITERATIVE, false,
                                       TrackFilterType::ALL, VbsRegion::R1 };
  inline const Score Score::VBF_R2 = { 23, "WAVeS [VBF R2: fwd PU + cen. HS]" , "VBF_R2", false, false, -1.f,
                                       -1.0, ClusteringMethod::ITERATIVE, false,
                                       TrackFilterType::ALL, VbsRegion::R2 };

  // Scores with a dedicated collection (distCut ≥ 0 → buildsCollection() = true)
  inline const Score Score::CONE       = {  7, "Cone"                       , "CONE",     true , false, -1.f, DIST_CUT_CONE, ClusteringMethod::CONE };
  inline const Score Score::FILTJET    = {  9, "Filter Tracks in Jets"      , "FILTJET",  true , false, -1.f, DIST_CUT_CONE, ClusteringMethod::CONE, false, TrackFilterType::JET };
  inline const Score Score::TEST_HS    = { 14, STR_TRKPTZ + " (HS only)"   , "TEST_HS",  true , false, -1.f,  DIST_CUT_CONE, ClusteringMethod::CONE, false, TrackFilterType::HS_ONLY };

  inline const std::vector<Score> SCORE_REGISTRY = {
    Score::HGTD    , Score::TRKPT     , Score::TRKPTZ, Score::PASS,
    Score::CONE    , Score::FILTJET   , Score::HGTD_SORT,
    Score::CONE_BDT, Score::TEST_MISAS, Score::TEST_HS,
    Score::WAVES,    Score::JET_T_REFINED, Score::WAVES_MISCL, Score::WAVES_MISAS,
    Score::VBF_R1,   Score::VBF_R2,
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
