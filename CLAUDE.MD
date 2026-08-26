# CLAUDE.MD - Vertex Timing Analysis Project

## Project Overview

This is a High-Granularity Timing Detector (HGTD) analysis project for ATLAS particle physics research. The project focuses on vertex timing reconstruction and clustering algorithms for analyzing collision events in high-pileup environments.

## Core Functionality

The project performs timing-based vertex reconstruction using:
- HGTD detector timing information
- Track-based timing measurements
- Jet clustering and filtering algorithms
- Multi-variable analysis for vertex identification

## Project Structure

```
vertex_timing/
├── src/                           # Core library headers + main analysis executables
│   ├── clustering_constants.h    # Analysis constants, Score struct + SCORE_REGISTRY, ClusteringMethod enum
│   ├── clustering_functions.h    # Clustering algorithms (clusterTracksInTime, chooseCluster, ...)
│   ├── clustering_includes.h     # Common ROOT/C++ includes
│   ├── clustering_structs.h      # BranchPointerWrapper, Cluster (incl. updateScores/calcFeatures), EventCounts
│   ├── sample_config.h           # --sample=vbf|zjets|dijet + --hist-file= resolution (default: local ../../ntuple-hgtd/)
│   ├── histogram_io.h            # HistWriter/HistReader: TFile-based histogram+scalar I/O shared by the
│   │                             #   *_hist/*_plot split executables (see "Main Executables" below)
│   ├── event_processing.h        # setupChain, getAssociatedTracks, selectClusters, processEventData
│   ├── plotting_utilities.h      # PlotObj / AnalysisObj + all plotting helpers (moneyPlot, inclusivePlot, ...)
│   ├── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
│   ├── clustering_dt_common.h    # Shared by clustering_hist.cxx/clustering_plot.cxx: Scenario enum,
│   │                             #   buildAnalysisMap, makeComparisonPlots
│   ├── clustering_hist.cxx       # Histogramming stage: event loop only, writes <sample-prefixed>_clustering_hist.root
│   ├── clustering_plot.cxx       # Plotting stage: reads that file, produces all PDFs, no ntuple access
│   └── clustering_dt.cxx         # Legacy monolithic executable (event loop + plotting in one run); retained
│                                 #   only as the validation ground truth for the split above, see Main Executables
├── util/                          # Diagnostic / report / ML-data-export executables
│   ├── failure_decomposition.cxx # Per-event principal-cause failure classification (see section below)
│   ├── export_training_data.cxx  # Exports per-cluster CSV rows (8 features + purity/labels) for future DNN retraining
│   ├── generate_rpt.cxx          # RpT-style timing/discrimination report (older version)
│   ├── rpt_v2.cxx .. rpt_v4.cxx  # RpT discrimination study iterations
│   ├── rpt_v5_common.h           # Shared by rpt_v5_hist.cxx/rpt_v5_plot.cxx: Scenario struct, RpT histogram
│   │                             #   binning/booking, saveScenarios/loadScenarios
│   ├── rpt_v5_hist.cxx           # Histogramming stage: event loop only, writes <sample-prefixed>_rpt_v5_hist.root
│   ├── rpt_v5_plot.cxx           # Plotting stage: ROC/console-summary/PDF from that file, no ntuple access
│   ├── rpt_v5.cxx                # Legacy monolithic executable; retained only as validation ground truth
│   ├── assoc_study_common.h      # Shared by assoc_study_hist.cxx/assoc_study_plot.cxx: assocRules() rule
│   │                             #   list, buildAssocMap, scalar-key naming
│   ├── assoc_study_hist.cxx      # Histogramming stage of the z-association study: one event loop, every
│   │                             #   AssocRule evaluated against the same event
│   ├── assoc_study_plot.cxx      # Report stage: console tables + CSV block + PDFs, no ntuple access
│   ├── maxpt_jet_test.cxx        # A/B/C test of the MAX_TRACK_PT cap vs per-event leading-jet pT
│   ├── track_dt.cxx              # Inclusive per-track HGTD-time-vs-truth-particle-time residual
│   └── sweep_utilities.h         # Shared helpers for parameter-sweep executables (no sweep .cxx files currently present)
├── python/                        # Interactive/visualisation scripts (run from within this dir)
│   ├── event_display.py          # Event visualisation — calls runHGTD_Clustering.cxx, produces 2-page PDF
│   ├── runHGTD_Clustering.cxx    # Interactive per-event clustering
│   ├── clustering_animation.py   # Animates the iterative clustering algorithm for one event
│   └── event_tinkering.py        # Quick event inspection and jet printing via uproot
├── share/
│   ├── models/                   # ML artifacts from a prior (currently inactive) DNN/BDT effort — see
│   │                             #   Machine Learning section below
│   ├── scripts/                  # ONNX export/inspection Python utilities (historical; no active consumer in src/)
│   └── docs/                     # Design docs, incl. older ML integration guides
├── condor/                        # HTCondor grid submission templates (see Grid Submission)
│   ├── run_analysis.sh           # Generic job executable (any --sample/--threads-aware target)
│   ├── clustering_dt.sub         # Submit file: one job per sample, for the legacy clustering_dt
│   ├── clustering_hist.sub       # Submit file: one job per sample, for clustering_hist (histogramming stage)
│   ├── rpt_v5.sub                # Submit file: one job per sample, for the legacy rpt_v5
│   ├── rpt_v5_hist.sub           # Submit file: one job per sample, for rpt_v5_hist (histogramming stage)
│   └── logs/                     # condor stdout/stderr/log output
├── CMakeLists.txt
└── build/                        # Build output (created by `cmake ..`)
```

Note: several directories/files described in earlier versions of this doc
(`vertex_time.cxx`, `check_jet_sorting.cxx`, `test_vbf_*.cxx`, `DBSCAN.ipynb`,
`figs/`, `old/`, `src/ml_model.h`, and most of the `util/` diagnostics and
parameter-sweep executables) do not exist in the current tree — this section
reflects what's actually present.

## Main Executables

Since the histogramming/plotting split (see "Histogramming vs. Plotting
Split" below), the two primary analyses each ship as a **pair** of
executables: a `*_hist` stage that runs the `TTreeProcessorMT` event loop
against the ntuples and writes every histogram to a ROOT file, and a
`*_plot` stage that reads that file back and produces the PDFs, touching no
ntuples at all. The original single-executable versions (`clustering_dt`,
`rpt_v5`) are kept building as a reference/fallback but are no longer the
primary workflow.

### Primary Analysis
- **clustering_hist.cxx** / **clustering_plot.cxx** (`src/`): Split of the
  main clustering analysis. `clustering_hist` runs the event loop (only the
  HGTD timing scenario is currently built/processed — `Scenario::IDEAL_RES`/
  `IDEAL_EFF` and `buildAnalysisMap`'s branching for them still exist for
  future reuse, but `main()` no longer instantiates or iterates those maps)
  and writes its histogram file via `MyUtl::histFilePath()` (see
  "Sample-prefixed output paths" below). `clustering_plot`
  reads that file and produces the `theplot_*`, `waves_oracle_*`, and
  `hgtd_trkptz_waves_*` comparison PDF groups (see Score Registry section),
  plus inclusive resolution/pull plots and shape-comparison plots — no
  ntuple access, runs in seconds. Active scores: `TRKPTZ`, `WAVES`,
  `JET_T_REFINED`, `TEST_MISAS`, `WAVES_MISCL`, `WAVES_MISAS`, plus `HGTD`.
  Iterative clustering (`ClusteringMethod::ITERATIVE`, `DIST_CUT_CONE =
  3.0σ`) is the default main collection. Shared booking/plotting logic
  (`Scenario` enum, `buildAnalysisMap`, `makeComparisonPlots`) lives in
  `src/clustering_dt_common.h`.
- **clustering_dt.cxx** (`src/`): Legacy monolithic executable (event loop +
  plotting in one run) that the split above was validated against
  bin-for-bin/pixel-for-pixel; still builds, but `clustering_hist` +
  `clustering_plot` is the primary workflow going forward.

### Diagnostics & Failure Analysis
- **failure_decomposition.cxx** (`util/`): Assigns every TRKPTZ-failing event
  to a single principal cause via a fixed precedence ladder (see Failure
  Decomposition section below).

### RpT / Report Studies
- **rpt_v5_hist.cxx** / **rpt_v5_plot.cxx** (`util/`): Split of the current
  RpT (jet-level pileup-discrimination) study. `rpt_v5_hist` runs the event
  loop (same `TTreeProcessorMT`/`ThreadState` structure as the legacy
  version) and writes its histogram file via `MyUtl::histFilePath()` (see
  "Sample-prefixed output paths" below). `rpt_v5_plot` reads that file and
  does `generate_roc_grid`, `rejAtEff`, the harmonized `[0.80, 1.0]` ROC
  window (see below), the console summary tables, and all PDF drawing — no
  ntuple access. Shared `Scenario`/binning/save-load logic lives in
  `util/rpt_v5_common.h`.

  **Six scenarios**, in `makeScenarios()` order — the index is load-bearing,
  since `fillJets` fills `sv[0..5]` positionally:

  | idx | name | what it is |
  |---|------|------------|
  | 0 | `zonly`       | ITk-only, no time gate |
  | 1 | `hgtd`        | Athena `RecoVtx_time` |
  | 2 | `trkptz`      | TRKPTZ-selected cluster time |
  | 3 | `waves`       | WAVeS-selected cluster time (in-jet refined) |
  | 4 | `waves_ideal` | ORACLE: re-clustered with 30 ps track times AND the cluster closest in time to truth |
  | 5 | `truth`       | vertex = Gaus(truth HS vtx, 10 ps), tracks = Gaus(own truth vtx, 30 ps) |

  **Four jet populations**, each filled for all six scenarios: two forward pT
  slices (`_lo` 30–40 GeV, `_hi` >40 GeV, both 2.4 < |eta| < 3.8), the same
  two at central |eta| < 2.4 (`_lo_cen`/`_hi_cen`), and the two VBS regions
  (`_r1`/`_r2`). The central pair is an ITk-only baseline AND a standing
  correctness check: |eta| < 2.4 is outside HGTD acceptance, so all six rows
  must agree there. That check has caught two real bugs (see Physics
  Findings).

  Both idealised rows only smear tracks with `trackTimeValid == 1` — idealising
  the time of a track HGTD never measured invents coverage it does not have,
  which the central baseline exposes immediately. Their smeared times are
  generated ONCE per track and shared by the re-clustering and both gates, so
  a row never clusters on one draw and gates on another. The per-thread RNG is
  re-seeded per event from event-stable quantities, so results do not depend on
  how `TTreeProcessorMT` distributes entries. `clusterTracksInTime`'s own
  `useSmearedTimes` path is NOT used: it calls `getSmearedTrackTime`, which
  draws from the global `gRandom` and is a data race under MT.
- **generate_rpt.cxx**, **rpt_v2.cxx**..**rpt_v4.cxx** (`util/`): Earlier
  iterations of the same study (z-only, real-HGTD, TRKPTZ-selected timing
  gates), monolithic (event loop + plotting), not split.
- **rpt_v5.cxx** (`util/`): Legacy monolithic executable that
  `rpt_v5_hist`/`rpt_v5_plot` was validated against; still builds, but the
  split pair is the primary workflow going forward.

### ML Data Export
- **export_training_data.cxx** (`util/`): Exports one CSV row per
  cluster/passing event (8 features + purity/labels) mirroring
  `clustering_dt`'s event/track selection exactly. Produces training data for
  a possible future DNN, but no model currently consumes it in `src/` — see
  Machine Learning Model section.

### z-Association Study
- **assoc_study_hist.cxx** / **assoc_study_plot.cxx** (`util/`): Compares the
  two track-to-vertex z-association rules (z0-significance cut vs. the
  `getNewDzpara` parameterization) by measuring the vertex t0 each produces,
  for both `TRKPTZ` and `WAVES`. Reports inclusive core σ and core fraction,
  plus both as a function of n forward HS tracks and n forward jets, plus
  track-level accounting (kept / HS purity / HS recall). Same hist/plot split
  as the primary analyses; the rule list lives in `assocRules()`
  (`util/assoc_study_common.h`). Full write-up and results:
  `share/docs/Z_ASSOCIATION_STUDY.md`; headline findings in Physics Findings
  below.

  All rules are evaluated against the SAME event in one pass, which is what
  makes the comparison paired. Two structural points make that work:
  - The counting scan stays pinned at `COUNTING_NSIGMA` (3.0) for every rule,
    so `EventCounts` — and therefore every plot x-axis — is common. Otherwise
    the "n_hs_tracks = 3" bin would hold a different set of events per rule.
  - Under a rule the clustering list is **re-selected from the full track
    array**, not narrowed from the counting list: a parameterized rule is not
    a subset of a significance cut (`getNewDzpara` is looser than 3σ for some
    (η, pT) and tighter for others), so narrowing would silently intersect
    every parameterized rule with `significance < 3.0`.

### Utilities
- **maxpt_jet_test.cxx** (`util/`): Diagnostic comparing WAVeS/TRKPTZ pass
  rates under different `maxTrkPt` caps (current `MAX_TRACK_PT`, leading-jet
  pT, subleading-jet pT). Pure diagnostic, changes nothing.
- **track_dt.cxx** (`util/`): Inclusive per-track HGTD-time-vs-truth-particle
  residual, no event selection.
- **runHGTD_Clustering.cxx** (`python/`): Interactive per-event clustering.
  Output is parsed by `event_display.py` to produce a 2-page PDF.

## Histogramming vs. Plotting Split

Both primary analyses follow the FastFrames pattern: run the event loop once
against the ntuples, save every histogram needed for plotting into a ROOT
file, then do all plot styling/derivation/output as a separate, fast, local,
iterable step that never touches the ntuples again. This exists because plot
tweaks (axis range, binning, a new working-point table, a style fix) used to
force a full high-stats rerun (zjets alone is 8M events, multi-minute local /
condor-scale for real work) even though the underlying data hadn't changed.

- **`src/histogram_io.h`**: shared `MyUtl::HistWriter`/`MyUtl::HistReader` —
  a thin RAII wrapper around a ROOT `TFile` for writing/reading a named set
  of histograms plus scalar run metadata. Naming convention: every object
  keeps its own existing ROOT name (`GetName()`) as the TFile key — no
  `TDirectory` nesting — since every histogram booked in this codebase
  already embeds enough context (score/variable/scenario/slice) to be
  globally unique within a run. Scalars use `TParameter<Long64_t>`/
  `TParameter<Double_t>` under a `meta_`-prefixed key (event counts,
  accumulator doubles); the energy label is a `TObjString` under
  `meta_energy_label`.
- **Sample-prefixed output paths** (`src/sample_config.h`'s `histFilePath()`/
  `plotFilePath()`): when a run is given `--sample=<name>`, every output file
  name is prefixed with `<name>_` so files stay distinguishable if copied out
  of their per-sample directory into a shared one. For histograms this also
  flattens away the old `hists/` subdirectory — a single sample-prefixed file
  directly under `OUTPUT_DIR` no longer needs it to stay unambiguous — e.g.
  `vbf/vbf_clustering_hist.root`, `vbf/vbf_rpt_v5_hist.root`. Plots keep their
  existing subdirectory layout (`comparisons/`, `inclusive/`, `fullplots/`,
  `rpt_plots/`, ...); only the file name gets prefixed, e.g.
  `vbf/comparisons/vbf_theplot_pu_frac.pdf`,
  `vbf/rpt_plots/vbf_rpt_v5.pdf`. For the local/default run (no `--sample`),
  both helpers fall back to the original unprefixed conventions
  (`<OUTPUT_DIR>/hists/<name>.root`, `<OUTPUT_DIR>/<subdir>/<name>.pdf`) —
  `MyUtl::SAMPLE_NAME` is only set to a non-empty value by the four primary
  split executables (`clustering_hist`/`clustering_plot`/`rpt_v5_hist`/
  `rpt_v5_plot`) when `resolveSample()` matched a real `--sample=` flag; the
  legacy monolithic executables (`clustering_dt.cxx`, `rpt_v5.cxx`) never set
  it and so keep their original unprefixed, un-flattened output paths
  unconditionally, even when run with `--sample=`.
- **Reconstruction pattern**: the plotting stage always "construct-empty-
  then-load" — it re-runs the *same* booking logic the histogramming stage
  used (`buildAnalysisMap()` / `makeScenarios()`) to get correctly-named,
  correctly-binned empty histograms, then calls `loadFrom()`/`loadScenarios()`
  to fill them from the file. This reuses the existing binning/title/name
  construction code instead of duplicating it in a second place, and doubles
  as a binning-mismatch guard (`HistReader::LoadInto` throws if a saved
  histogram's binning doesn't match the freshly-booked destination's).
- **`AnalysisObj`/`PlotObj::saveTo`/`loadFrom`** (`src/plotting_utilities.h`):
  mirror `mergeFrom`'s exact member list (`WriteHist`/`LoadInto` in place of
  `Add`) — deliberately exclude `params`/`effEstimate`/`slicesHists`/
  `slicesFits`/`slicesModels`/`efficiency`, since those are purely derived by
  `plotPostProcessing()` from the saved raw histograms and must be
  recomputed after `loadFrom()`, not persisted.
- **rpt_v5's ordering constraint, preserved structurally**: `generate_roc()`/
  `rejAtEff()` use explicit `Integral(1, bin+1)` bin ranges specifically so
  they're immune to any prior `SetRangeUser()` restriction; the console
  summary's no-arg `Integral()`/`GetMean()` calls are not immune, so they
  must run before `styleScen()`'s cosmetic `SetRangeUser(0.0, 1.5)`. In
  `rpt_v5_plot.cxx` this is enforced structurally, not just by comment:
  `styleScen` is a local lambda whose *definition* is placed lexically after
  the entire derived-numbers block, so calling it earlier would be a compile
  error (undeclared identifier), not a silent runtime bug.
- **Harmonized ROC window**: the ROC x-axis ("Hard Scatter Efficiency") is a
  fixed `[0.80, 1.0]` window for every sample and pT-slice, so plots are
  directly comparable across vbf/zjets/dijet and against the reference study,
  whose ROCs use the same range. 0.80 is also where `rocEffGrid()`'s first
  sampled working point sits — an earlier `[0.85, 1.0]` zoom cropped that
  point off the drawn window. Tunable by editing
  `roc_xmin_default`/`roc_xmax_default` in `rpt_v5_plot.cxx` and rerunning
  `rpt_v5_plot --sample=<sample>` locally against an already-saved histogram
  file — no condor resubmission.
- **Fixed-efficiency ROC sampling**: `generate_roc_grid()` samples each curve
  at the working points in `rocEffGrid()` (0.800 → 0.975 in 0.025 steps)
  rather than once per histogram bin, so every scenario shares identical
  x-points and the ratio panels are exact. Beware when reading the console
  tables: `rejAtEff` returns the first bin with eff <= target, so near a step
  a 0.02% efficiency difference can move the printed rejection by 30%.
  Repeated identical values across working points mean the curve's ceiling
  has been hit and the sampler is clamping, not that the physics is flat.
- **Per-panel y-maxima**: the inclusive slices, the central baseline, and the
  R1/R2 regions each compute their own. They are far apart in scale — with
  the oracle rows, R1 peaks near 880 and R2 near 1180 against an inclusive
  ceiling near 640 — so a shared maximum clips panels.
- **Selection labels wrap** on an explicit `@@` delimiter supplied by the
  caller. The legend is an opaque box whose left edge is at NDC 0.58 and which
  spans the same y range as the label, so a long single line is painted over.
  Splitting on a width estimate does not work (TLatex markup means glyph count
  != string length) and `|` cannot be the delimiter (it appears in every
  `|#eta|` label).
- **Validated** (bin-for-bin histogram identity + pixel-for-pixel PDF
  identity) against the local default VBF sample for both split pairs before
  ever touching a grid sample — see git history for the validation method
  (temporary `saveTo()` instrumentation in the legacy executables, `uproot`
  bin-content diff, `pdftoppm`+ImageMagick `compare` page-by-page PNG diff).

## WAVeS Clustering / Scoring

`Score::WAVES` (jet-proximity weighted score) is the current alternative to
the baseline `Score::TRKPTZ`, implemented in `Cluster::updateScores`
(`src/clustering_structs.h`) and `Cluster::calculateTime`/`calculatePurity`
(`src/event_processing.h`):

- **Score**: `Σ_i pT_i · pT_jet(i) / max(ΔR_i, WAVES_DR_FLOOR)`, summed over
  the cluster's constituent tracks (nearest qualifying forward jet per
  track), then multiplied by `exp(−1.5·|Δz_cluster|)`. Falls back to the
  TRKPTZ score if no qualifying forward jets exist in the event.
- **Time/purity**: for the `WAVES` row specifically (not its oracle
  variants), both the reported time and purity are recomputed using only the
  subset of the selected cluster's tracks that fall within `ΔR < 0.4` of a
  qualifying forward jet — this purifies the time by dropping absorbed PU
  tracks that the jet-agnostic clustering step picked up.
- **`Score::JET_T_REFINED`** (short name `WAVES_RECLUST`, id 19): a
  *dedicated* collection — tracks are first filtered to those within
  `ΔR < 0.4` of a forward jet (`TrackFilterType::JET`), then re-clustered
  with `ClusteringMethod::ITERATIVE` at the tighter `DIST_CUT_T_REFINED =
  2.0σ`. The best cluster in that collection is selected by the TRKPTZ score
  (via the aux-collection path in `selectClusters`).
- **`Score::WAVES_MISCL`** / **`Score::WAVES_MISAS`**: oracle variants that
  reuse the WAVeS selection and time, gated at fill time on cluster purity
  (`> 0.75`) / event-level HS timing purity (`≥ 0.95`) respectively — see
  `requiresPurity` in the Score Registry section. Unlike the `WAVES` row,
  these use the default full-cluster purity (not the in-jet-only
  recomputation), so their purity gate stays comparable to `TEST_MISAS`.

`Score::CONE` (legacy cone clustering, `ClusteringMethod::CONE` at
`DIST_CUT_CONE = 3.0σ`) still exists in the registry as a historical
baseline but is not currently emplaced in `buildAnalysisMap`.

## Failure Decomposition

`util/failure_decomposition.cxx` assigns every event where the selection
score fails the efficiency window to a **single principal cause**, via a
fixed precedence ladder (first match wins) — this superseded an earlier
non-exclusive four-category scheme:

| Precedence | Name | Condition |
|---|------|-----------|
| 1 | Timing misassignment | Ideal Resolution scenario passes, **or** a majority of HS pT is mis-timed vs. truth (`|pull| ≥ 3σ`) — the latter isn't fixed by ideal resolution since the times are offset, not just smeared |
| 2 | Selection failure | Another qualifying cluster in the same real-time collection passes `passEfficiency()` — the right answer existed but wasn't chosen |
| 3 | Misclustering | Grouping *all* truth-HS tracks into one cluster (real times, no PU contamination) passes the window, but the algorithm's chosen cluster does not |
| 4 | Track efficiency loss | Ideal Res. + Ideal Eff. passes but Ideal Res. alone does not — giving currently time-less HS tracks an ideal time is what recovers the event |
| 5 | Irreducible | None of the above — not an algorithm failure |

## Key Analysis Parameters

From `src/clustering_constants.h`:
```cpp
MIN_TRACK_PT          = 1.0   // GeV — track selection floor
MAX_TRACK_PT          = 30.0  // GeV — track selection ceiling
DIST_CUT_CONE         = 3.0   // σ — main iterative clustering distance (also CONE baseline)
DIST_CUT_T_REFINED    = 2.0   // σ — JET_T_REFINED (WAVES_RECLUST) stage-2 iterative distance
MIN_CLUSTER_TRACKS    = 0     // minimum tracks required to select a cluster
IDEAL_TRACK_RES       = 1.0   // ps — flat per-track resolution for ideal scenarios
PASS_SIGMA            = 60.0  // ps — ±window for "Core Fraction" efficiency oracle
TRUTH_PULL_CUT        = 3.0   // σ — |pull| threshold for oracle timing purity gates
WAVES_DR_FLOOR        = 0.05  // minimum ΔR for the WAVeS 1/ΔR weight (prevents divergence)
MAX_NSIGMA            = 3.0   // z0-significance track-to-vertex association (incumbent)
Z0_VAR_INFLATION      = 1.15² // covariance underestimates true z0 resolution by ~15%
DZ0_PARA_SCALE_CLUSTER= 1.4   // multiple of getNewDzpara used when --dzpara is passed
COUNTING_NSIGMA       = 3.0   // association PINNED for the EventCounts scan under an AssocRule
```

**Track-to-vertex association.** Three ways to express it, all in
`clustering_constants.h`:
- `MAX_NSIGMA` — the incumbent z0-significance cut.
- `USE_DZ_PARA` (runtime, from `--dzpara`) — swaps in `getNewDzpara ×
  DZ0_PARA_SCALE_CLUSTER` instead; the caller's significance argument is then
  ignored, since the two rules are not on a common scale.
- `struct AssocRule` — the same information as an explicit *value* (kind +
  cut), so several rules can be evaluated against one event in a single pass.
  `passTrackVertexAssociation` / `getAssociatedTracks` have AssocRule
  overloads that deliberately do NOT consult the globals, and
  `processEventData` takes an optional `assoc` argument. `nullptr` (the
  default) reproduces the global-driven behaviour exactly, so every executable
  that does not pass one is bit-for-bit unchanged. Only
  `util/assoc_study_hist.cxx` uses it today.
- `AssocRule::orGhost` — unions the rule with ghost association (accept a
  track if it passes the z-test OR is ghost-associated to a qualifying forward
  jet, via `forwardJetGhostTracks`). **This is a SET-level property, not a
  per-track one**: it needs the event's jets, so it is applied inside
  `getAssociatedTracks` and is *not* implemented by
  `passTrackVertexAssociation`. Anything reproducing a rule's track list must
  therefore call `getAssociatedTracks` rather than looping the per-track
  function, or it will silently get the z-only list. (Studied and rejected on
  physics grounds — see Physics Findings.)

The kinematic half of the track selection (η window, pT window, quality flag)
lives in `passTrackKinematics` so the two `getAssociatedTracks` overloads
cannot drift apart — only the association test differs between them.

The z-association study (see below) found the incumbent `MAX_NSIGMA = 3.0` is
**not** the best choice for TRKPTZ. It has not been changed, because the study
has so far only been run on the local VBF sample.

RpT-specific, from `util/rpt_v5_hist.cxx`:
```cpp
JET_ETA_MIN / MAX     = 2.4 / 3.8  // forward jet acceptance
CENTRAL_ETA_MAX       = 2.4        // central baseline, |eta| < this
RPT_TRACK_JET_DR      = 0.2        // track must be within this of the jet axis
DZ0_PARA_SCALE        = 1.4        // x getNewDzpara, FORWARD track-to-vertex association
CENTRAL_Z_SIGNIF      = 5.0        // |z0-z_vtx|/sigma_z0, CENTRAL association
GATE_SIGMA            = 2.5        // per-track time-gate half-width, in pulls
TRUTH_VTX_SMEAR       = 10.0       // ps, truth-t0 row's vertex smear
TRUTH_TRK_SMEAR       = 30.0       // ps, idealised per-track smear
```

**The two z-association rules are not interchangeable.** `getNewDzpara` comes
from a FORWARD-region reference study and does not extrapolate inward: it
returns sigma_z0 = 31 um at eta = 0 rising to 4.4 mm at eta = 3.8, and the
central end is several times tighter than ITk z0 resolution at 1 GeV. Using it
centrally roughly halves ITk-only rejection (163 -> 77 at 0.800 HS efficiency).
5.0 is not a guess: it is the value the reference macro itself carries as a
commented alternative, and an independent scan of ours reproduces it as the
best working point.

**Vertex-time inflation is per sample AND per scenario**, resolved by
`inflationFor(SAMPLE_NAME)` in `main()` before any worker thread starts:

| scenario | vbf | zjets | dijet |
|---|---|---|---|
| hgtd | 1.48 | 1.61 | 1.53 |
| trkptz | 1.39 | 1.63 | 1.45 |
| waves | 1.38 | 1.65 | 1.46 |

Measured off `PRINT_PULL_DIAG`'s ratio column; the run echoes `[calibration]`
at startup and prints the freshly measured ratio beside the value in use, so a
stale entry is visible in the same log. The measurement is NOT circular (raw
`var_vtx`, ungated track list, vertex times independent of the gate), so one
pass suffices. Re-measure for any NEW sample. See Physics Findings for why
this matters far less than it looks.
`CONE_ITER_K` also still exists as a constant but its only historical
consumer (a `cone_iter_k_sweep.cxx` sweep executable) is not present in this
tree, so it's currently unused.

## Build System

### Requirements
- CMake ≥ 3.10
- C++17 compiler
- ROOT (RIO, Core, Gpad, Graf, Hist, Tree)
- Boost (filesystem)

TMVA is intentionally excluded (see `CMakeLists.txt:34`) — `HGTD_SORT`/
`CONE_BDT` are registered `Score`s but their selector
(`chooseHGTDSortCluster`) is currently a stub returning an empty `Cluster`,
and neither is emplaced in `buildAnalysisMap`.

### Build Commands
```bash
cd build
cmake ..
make [<target>]   # e.g. make failure_decomposition
```

### Compiler Flags
`-O3 -march=native -g` (no `-DNDEBUG`, no LTO — `CMakeLists.txt:37`).

## Data Format

ROOT ntuples (tree name `"ntuple"`) with the branches bound in
`BranchPointerWrapper` (`src/clustering_structs.h`):

| Group | Branches |
|-------|----------|
| Tracks | `Track_z0/d0/pt/eta/phi/theta/qOverP`, `Track_var_z0/d0/qOverP/theta/phi0`, `Track_time/timeRes/hasValidTime`, `Track_quality`, `Track_truthVtx_idx/truthPart_idx`, `Track_nHGTDHits/nHGTDPrimaryHits` |
| Vertices | `TruthVtx_z/time/isHS`, `RecoVtx_z/time/timeRes/hasValidTime` |
| Reco jets | `AntiKt4EMTopoJets_pt/eta/phi`, `AntiKt4EMTopoJets_truthHSJet_idx`, `AntiKt4EMTopoJets_ghostTrack_idx` |
| Truth jets | `TruthHSJet_pt/eta/phi`, `TruthITPUJet_pt/eta/phi`, `TruthOOTPUJet_pt/eta/phi` |
| Particles | `TruthPart_prodVtx_time` (bound as `particleT`) |
| Event | `weight` |

## Event Display

```bash
cd python && python event_display.py --file_num <NUM> --event_num <ENTRY> --extra_time <dt>
```

Internally calls `runHGTD_Clustering.cxx` and parses its stdout to produce a
3-page PDF. Script lives in **`python/`** (not `share/scripts/`).

The VBS candidate pair is marked with a muted steel-blue cross-hatch drawn
UNDER the jet wedges (fill colour is already spoken for: green = truth HS,
grey = PU), and each tagging jet's label carries `R_pT = <z-only> -> <timed>`.
The timed value applies rpt_v5's gate, so it needs `--infl` set to the
sample's WAVeS inflation (vbf 1.38, zjets 1.65, dijet 1.46). t0 comes from
`--extra_time` directly — NOT from the matched cluster's own weighted-mean
time, which differs from WAVeS's in-jet refined time and made the displayed
R_pT disagree with the analysis.

`--jet_idx`/`--jet_label` paint a target jet orange, which OVERRIDES the
truth-identity colouring. Do not pass them for region displays: which leg is
HS and which is PU is the whole point of those plots.

**The local box cannot run this without help.** Its ROOT needs python3.14
while `~/.local/share/python/lib/python/site-packages` (a version-agnostic
user site dir) holds python3.9's numpy, whose compiled extensions fail to load
under 3.14. Set `PYTHONNOUSERSITE=1` to skip it and use `~/.venv-hgtd`, which
has a matching matplotlib. On the AF use ROOT 6.38.04, not an LCG_106 view:
that ships ROOT 6.32, where `TF1::GetParError` has no `const char*` overload
and `src/plotting_utilities.h` fails to compile.

## Machine Learning Model

**Currently inactive.** There is no `src/ml_model.h`, no ONNX runtime
integration, and no `Score::TEST_ML` in the current `Score` registry — the
WAVeS jet-proximity score (see section above) is the current alternative to
the TRKPTZ baseline. What remains from the prior ML effort:
- `share/models/` still holds artifacts (`model.onnx`, `model_weights.json`,
  `normalization_params.json`, a TMVA `.weights.xml`) from that earlier
  work, but nothing in `src/` currently loads or runs inference against
  them.
- `share/scripts/export_onnx_model.py` (SavedModel → ONNX → JSON weights
  pipeline) and `share/docs/ML_MODEL_INTEGRATION.md` /
  `ONNX_INTEGRATION_GUIDE.md` document that pipeline for reference, but it
  has no active consumer.
- `util/export_training_data.cxx` still exports per-cluster training CSVs
  (8 features, purity/labels) mirroring the current `clustering_dt`
  selection, in case retraining is picked back up later.
- `src/clustering_structs.h`'s `Cluster::calcFeatures()` still computes and
  normalises the same 8-feature vector (used historically as DNN input); its
  `rawDeltaZ`/`rawDeltaZResunits` outputs are reused by `updateScores` for
  TRKPTZ, but the normalised `features` vector itself is otherwise unused
  today.

## Common Workflows

```bash
# Full analysis (two-step: histogram once, then plot — see "Histogramming
# vs. Plotting Split" above)
cd build && ./clustering_hist && ./clustering_plot

# Iterate on plotting only (styling, binning, a new comparison group) without
# re-running the event loop, as long as ../figs/hists/clustering_hist.root
# is already there:
cd build && ./clustering_plot

# Failure decomposition
cd build && ./failure_decomposition

# Event display for a specific event
cd python && python event_display.py --file_num 000009 --event_num 1856 --extra_time 0.0

# RpT study (WAVeS-based), two-step
cd build && ./rpt_v5_hist && ./rpt_v5_plot

# Iterate on rpt_v5's plotting only (e.g. tuning the harmonized ROC window's
# roc_xmin_default/roc_xmax_default)
cd build && ./rpt_v5_plot

# z-association study (11 rules in one pass), two-step. ~35 s + ~15 s on the
# full local VBF sample; add --sample=<name> to both stages for a grid sample.
cd build && ./assoc_study_hist --threads=8 && ./assoc_study_plot

# Iterate on the study's report/plots only (rules unchanged)
cd build && ./assoc_study_plot

# Rebuild after header change
cd build && make
```

## Grid Submission (condor)

Templates in `condor/` run `clustering_hist` / `rpt_v5_hist` on the UChicago AF HTCondor pool, one job per `--sample` value (`vbf`, `zjets`, `dijet`). Condor only runs the **histogramming** stage — the plotting stage is fast/local and runs afterward on whatever machine the `.root` file lands on, with no condor submission needed:

- **`condor/run_analysis.sh`** — generic executable, invoked as `run_analysis.sh <executable> <sample> [<threads>]`. Sources `atlasLocalSetup.sh` + `lsetup "root <version>"` (relaxing `set -e`/`-u` around just that call — those scripts reference unset internal variables like `ALRB_frontlineSite` and aren't strict-mode safe), then runs the transferred `<executable>` **directly from the job's scratch directory** with `--sample=<sample>` and, if given, `--threads=<N>`. No longer moves the binary into a fake `build/` subdir first (see `src/sample_config.h`'s output-path note below) — that indirection existed only to make the old `../<sample>` output convention resolve to the scratch root; the grid-sample output paths are cwd-relative now, so it's unnecessary and was a latent risk (if that step were ever skipped, `../<sample>` would climb out of the scratch directory into its parent, which may not be writable). Unchanged by the histogramming/plotting split otherwise — generic enough for both the legacy and `*_hist` executables.
- **`condor/clustering_hist.sub`** / **`condor/rpt_v5_hist.sub`** — one `.sub` file per histogramming executable, each queuing 3 jobs via `queue sample in (vbf, zjets, dijet)`. Submit each separately (own `condor_submit` call) to get independent job clusters. `arguments` passes `4` as the `<threads>` arg, matching `request_cpus = 4` — both `*_hist` executables parallelize their event loop via `TTreeProcessorMT` (see Multicore in Development Notes), and condor's cgroup throttles any threads spawned beyond `request_cpus`, so the two must track each other.
- **`condor/clustering_dt.sub`** / **`condor/rpt_v5.sub`** — legacy templates for the monolithic executables; still work, kept for reference/fallback, not the primary workflow.

Execute hosts do **not** share a filesystem with the submit host on this pool — confirmed via a `TARGET.FileSystemDomain == MY.FileSystemDomain` match failure when `should_transfer_files = NO` was tried first (0 of 157 slots matched). The `*_hist.sub` files use `transfer_input_files` (stages the prebuilt binary in) and `transfer_output_files = $(sample)/$(sample)_<name>.root` (stages back a single small histogram file — much smaller than the full PDF tree the legacy `.sub` files' `transfer_output_files = $(sample)` transfers). Ntuple inputs are read via the absolute `/data/mcardiff/exotic_superntuples/...` path baked into `--sample=` (see `src/sample_config.h`), so they don't need transferring.

**`transfer_output_remaps` is used** on both `*_hist.sub` files, though its role has narrowed since output file names became sample-prefixed (see "Sample-prefixed output paths" above): condor's default output-transfer behavior flattens any nested `transfer_output_files` path to its basename when staging back, it does not recreate `$(sample)/` in the submit directory on its own. Historically (when every sample wrote a bare `clustering_hist.root`/`rpt_v5_hist.root`), skipping the remap meant all 3 queued jobs silently overwrote each other's output — confirmed by running the real `rpt_v5_hist` jobs once without one: only the last-to-finish sample's histogram survived, and no `vbf/`, `zjets/`, `dijet/` subdirectories were created at all. Now that each job's basename is already unique (`vbf_clustering_hist.root`, `zjets_clustering_hist.root`, `dijet_clustering_hist.root`, ...), that particular collision can no longer happen — the remap is kept purely so the file lands nested under `$(sample)/` on the submit side, matching the local `--sample=` convention (so `clustering_plot --sample=<sample>` finds it without needing `--hist-file=`). (The legacy `.sub` files' `transfer_output_files = $(sample)` isn't affected either way — that argument has no nested subdirectory, so its basename already equals the whole thing they want back.) The obvious alternative, `preserve_relative_paths = True`, was tried first and rejected: it's a submit-wide knob that also governs `transfer_input_files`, and turned `../build/<exe>` (which deliberately climbs *out* of the sandbox to reach the built binary) from silently-flattened into a hard `Attempt to write to illegal sandbox path` hold — confirmed by running the real jobs with it set. `transfer_output_remaps` only touches the output side, so it can't have that side effect.

**`transfer_output_remaps`'s LHS (`name`) must be a bare filename** — per `condor_submit(1)`: "You cannot specify directories to be remapped." A first attempt using a full nested path (`$(sample)/hists/clustering_hist.root`, back when the histogram still lived under a `hists/` subdirectory) on the LHS silently failed to match anything and fell back to the same basename-flattening bug it was meant to fix — confirmed against the real jobs. Directories *are* allowed on the RHS (`newname`) and are created as needed, so the working form is `transfer_output_remaps = "$(sample)_clustering_hist.root=$(sample)/$(sample)_clustering_hist.root"` (bare sample-prefixed basename on the left, full destination path on the right).

```bash
cd build && cmake .. && make clustering_hist rpt_v5_hist   # build before submitting; condor doesn't rebuild
cd condor
condor_submit -dump dryrun.ad clustering_hist.sub          # dry-run: validate without queuing
condor_submit clustering_hist.sub
condor_submit rpt_v5_hist.sub                               # separate job cluster
condor_q
# After jobs complete, the .root files land in
# ../<sample>/<sample>_{clustering,rpt_v5}_hist.root. Plot locally — fast, no
# further condor submission:
cd ../build
./clustering_plot --sample=<sample>
./rpt_v5_plot --sample=<sample>
```

Logs land in `condor/logs/clustering_hist_<sample>.<ClusterId>.<ProcId>.{out,err,log}` (prefixed `rpt_v5_hist_` for the rpt_v5_hist cluster). `stream_output`/`stream_error` are **on**, so `.out`/`.err` are written back to the submit node as the job runs rather than only on exit — without that, `condor_tail` mid-job shows only whatever happened to be flushed.

### Job cost: where the time actually goes

Measured from the three-sample `rpt_v5_hist` clusters of 2026-07-03 (`condor/logs/`, three repeats each):

| sample | files | events | wall | CPU | wall/event |
|---|---:|---:|---:|---:|---:|
| vbf   | 1002 |   499,100 |  24.9 min | 12.2 min | 2.99 ms |
| dijet | 1414 |   699,500 |  25.2 min |  8.4 min | 2.16 ms |
| zjets | 2005 | 8,084,000 | 183.3 min | 108.0 min | 1.36 ms |

**Z+jets is not inefficient — it is 16.2x larger than vbf** (both DSIDs, Zee+Zmumu). Per event it is the *fastest* of the three. File count is not a driver either: dijet has 1.4x vbf's files and 1.4x its events at the same wall time.

**Startup dead time is the real inefficiency, and it is not proportional to anything.** From the condor `.log` image-size timeline (memory sits at ~60 MB with ROOT loaded and nothing allocated, then jumps to ~500 MB when the worker threads book their histograms):

| sample | run 1 | run 2 | run 3 |
|---|---|---|---|
| vbf   | 20.2 min |  5.1 min | 15.2 min |
| dijet | 10.2 min | 20.2 min | 10.1 min |
| zjets | 30.2 min | 25.2 min | 20.2 min |

That is 5–30 min, varying **4x run-to-run on the same sample** — filesystem contention, not a per-sample property. It is 15–20% of a zjets job but **40–80% of a vbf or dijet job**. Sampled at condor's 5-minute granularity, so ±5 min.

Mechanism: `TChain::Add` is lazy (opens nothing), but `TChain::GetEntries()` opens **every file** in the chain to sum entry counts — 0.3–1.2 s per file on `/data` against 4.3 ms on local SSD. That call has been removed from `setupChain` (the "did we find files?" question is answered by `GetListOfFiles()->GetEntries()` with no I/O) and from the `*_hist` progress denominators, leaving only `TTreeProcessorMT`'s own pass. **Do not reintroduce it** — see the note on `setupChain` in `src/event_processing.h`.

How much that actually saved, measured on the local 33-file sample (time from start to the first event reaching a worker):

| | serial startup |
|---|---:|
| with `GetEntries()` | 0.205 s (6.1 ms/file) |
| deferred | 0.103 s (3.1 ms/file) |

**~50%**, not more: `TTreeProcessorMT` still opens every file serially for cluster boundaries before dispatching any task. In the old ordering that second pass was nearly free (0.004 s) only because `GetEntries` had already warmed the files.

**What is left, and what does NOT fix it.** The remaining pass is latency-bound, so overlapping it looks attractive — opening every file concurrently first does collapse the serial scan to 0.004 s. But the prewarm then costs as much as it saves (0.132 s at 8 threads vs a 0.117 s scan): on local SSD the opens are bandwidth-bound and concurrency buys nothing. Whether the AF's shared storage behaves differently is untested and should not be assumed — an earlier plausible-sounding fix here (hoisting `passBasicCuts` above the Z+jets lepton work, so rejected events never touch per-track branches) measured at **1 µs/event**, i.e. 0.1 min over 8M events, and was dropped.

### Sharding a sample across jobs

`--file-shard=<i>/<N>` restricts a `*_hist` process to files i, i+N, i+2N, … of the **sorted** file list, and `hist_merge` combines the per-shard outputs. N jobs each pay 1/N of the per-file startup *and* run 1/N of the event loop concurrently, so it attacks both costs at once.

Both `*_hist.sub` files now take `sample` and `nshards` as submit macros and queue `$(nshards)` jobs for **one** sample, with `$(ProcId)` as the shard index — so each sample gets its own ClusterId and its own shard count:

```bash
cd build && cmake .. && make clustering_hist rpt_v5_hist hist_merge
cd condor
condor_submit -a sample=vbf   -a nshards=4  clustering_hist.sub
condor_submit -a sample=dijet -a nshards=4  clustering_hist.sub
condor_submit -a sample=zjets -a nshards=32 clustering_hist.sub
# after they land:
cd ../build
./hist_merge ../zjets/zjets_clustering_hist.root ../zjets/zjets_clustering_hist.shard*of32.root
./clustering_plot --sample=zjets                                 # unchanged
```

A non-default VBS selection needs **two** matching macros, because
`resolveSelection` stamps a `SELECTION_TAG` onto the output filename (so a
loosened run cannot overwrite the standard one) and `transfer_output_files`
must name that same file:

| `selargs` | `seltag` (note trailing `_`) |
|---|---|
| *(empty)* | *(empty)* |
| `--vbs-deta=0` | `deta0p0_` |
| `--vbs-mjj=200` | `mjj200p0_` |
| `--vbs-deta=0 --vbs-mjj=200` | `deta0p0_mjj200p0_` |

```bash
condor_submit -a sample=zjets -a nshards=32 \
              -a 'selargs=--vbs-deta=0 --vbs-mjj=200' \
              -a seltag=deta0p0_mjj200p0_ clustering_hist.sub
```

A mismatch is not silent — the job completes and then goes on hold failing to
transfer a file that does not exist — but `condor_submit -dump` shows the
expanded name, and running the executable locally with the same `selargs`
prints `[selection] ... (tag: ...)`.

**Shard counts should track EVENT count, not file count** — zjets has 8.08M events against vbf's 0.50M and dijet's 0.70M, so an equal split would leave zjets jobs ~16x longer than the rest. 4/4/32 puts every job in the same ballpark.

`--file-shard` tags the output filename whenever it is passed, **including at N = 1** (`Shard::given`, not `count > 1`). That looks redundant but is what lets one `.sub` template the output name on `$(ProcId)/$(nshards)` for any shard count: tying the tag to `count > 1` would make `nshards=1` write an untagged file the submit file then fails to transfer back — a foot-gun that only fires at the least-tested value.

The sort in `setupChain` is load-bearing, not tidiness: `directory_iterator` order is filesystem-defined and not guaranteed stable, so without it two shards could disagree about which file is "3rd" and silently double-process or skip data. It also makes an unsharded run deterministic, which it previously was not — harmless for full runs, but it changes *which* events a `--max-events` run sees.

**`hadd` will not do this job.** It merges anything with a `Merge()` method and keeps the **first occurrence** of everything else. Every scalar `HistWriter` emits is a `TParameter` — event counts, rejection tallies, the pT accumulators `rpt_v5_plot` turns into untimed-fraction percentages — and `TParameter` has no `Merge()`. `hadd` would produce histograms summed over N shards alongside one shard's event counts, silently. `util/hist_merge.cxx` sums them explicitly, and refuses (rather than dropping) any key type it does not recognise.

**Validation** (local VBF, `util/scratch/compare_hists.py`, unsharded vs sharded-and-merged):

| executable | shards | histograms bit-identical | worst relative bin diff |
|---|---|---|---|
| `rpt_v5_hist`      | 3 | **72 / 72**       | **0** (exactly) |
| `clustering_hist`  | 4 | 461 / 477         | 5.6e-13 |
| `assoc_study_hist` | 5 | 1326 / 1378       | 1.1e-12 |

Every scalar matched exactly in all three, and shard event counts summed to the unsharded total exactly (12353+11785+11635+11762 = 47535). `assoc_study_plot`'s full numeric output was **byte-identical** between the unsharded and merged files.

The only objects that differ are `TProfile2D`s, which accumulate sums of floating-point values whose addition order changes with the partition — machine epsilon, ~12 orders of magnitude below any physics tolerance. Every other histogram is bit-identical because all fills in this codebase are unweighted, so bin contents are integer-valued sums. `rpt_v5_hist` matching exactly also independently confirms that its per-event RNG re-seeding really is independent of how work is distributed (see the Diagnostics section) — sharding is a stronger test of that than varying the thread count.

The remaining dead time is the ATLAS/cvmfs `lsetup` in `run_analysis.sh`, which now prints `[job t+NNNNs]` stamps around itself so setup and analysis are separately attributable; the executables print `[t+NNN.N s]` phase marks (`MyUtl::PhaseTimer`). Both flush immediately — a block-buffered stdout with no flushes is why a working job looked hung under `condor_tail`, which is what sent an earlier investigation after "file loading" in the first place.

**Caveat:** the build uses `-march=native` (see Compiler Flags above), tuned to whichever machine ran `make`. If execute hosts have a different CPU microarchitecture than the build host, jobs can crash with an "Illegal instruction" (SIGILL) rather than a normal error.

## Score Registry

All scores are defined as `static const Score` members of the `Score` struct in `src/clustering_constants.h` and listed in `SCORE_REGISTRY`. Each carries its own collection spec: `distCut`, `method` (`CONE`/`ITERATIVE`/`SIMULTANEOUS`/`ITERATIVE_SPLIT`), `useZ0`, `filter` (`TrackFilterType`: `ALL`, `JET`, `HS_ONLY`), plus `buildsCollection()` to indicate it manages its own cluster collection (only true when `usesOwnCollection && distCut >= 0`).

| id | Short | Long name | Notes |
|----|-------|-----------|-------|
| 0  | HGTD  | HGTD Algorithm | Own collection; selected by proximity to `RecoVtx_time` |
| 1  | TRKPT | Σ pT | Plain pT sum |
| 2  | TRKPTZ | Σ pTe^{−1.5\|Δz\|} | Primary/baseline score |
| 3  | PASS  | Pass Cluster | First cluster passing quality cuts (not currently emplaced in `buildAnalysisMap`) |
| 7  | CONE  | Cone | Legacy cone baseline; own collection at `DIST_CUT_CONE`, `ClusteringMethod::CONE` (not currently emplaced) |
| 8  | CONE_BDT | Cone (BDT) | Threshold 0.3; selector is a stub (TMVA excluded) — currently non-functional |
| 9  | FILTJET | Filter Tracks in Jets | JET filter; own collection at `DIST_CUT_CONE`, `ClusteringMethod::CONE` (not currently emplaced) |
| 10 | HGTD_SORT | HGTD BDT (pT-sorted) | Threshold 0.3; selector is a stub (TMVA excluded) — currently non-functional |
| 13 | TEST_MISAS | ΣpTe^{−1.5\|Δz\|} [Events with Perfect Timing] | Oracle: `requiresPurity`, gated on event-level HS timing purity ≥ 0.95 |
| 14 | TEST_HS | ΣpTe^{−1.5\|Δz\|} (HS only) | Own collection; `HS_ONLY` filter, `ClusteringMethod::CONE` at `DIST_CUT_CONE` (not currently emplaced) |
| 18 | WAVES | WAVeS Score | Jet-proximity weighted score; main collection (see WAVeS section above) |
| 19 | JET_T_REFINED | WAVeS 2σ t Re-clustering | Own collection (`WAVES_RECLUST` short name); `JET` filter, `ClusteringMethod::ITERATIVE` at `DIST_CUT_T_REFINED`; selected by TRKPTZ |
| 20 | WAVES_MISCL | WAVeS [Events with Pure Clusters] | Oracle: `requiresPurity`, gated on cluster purity > 0.75 |
| 21 | WAVES_MISAS | WAVeS [Events with Perfect Timing] | Oracle: `requiresPurity`, gated on event-level HS timing purity ≥ 0.95 |

Oracle scores (`requiresPurity=true`) gate their denominator: events outside the gate are not counted. **Important:** `calcHSTimingPurity` (event-level, used by `TEST_MISAS`/`WAVES_MISAS`) returns `0.0f` — not `1.0f` — when no qualifying HS tracks exist, ensuring that events with zero forward HS tracks are excluded from oracle denominators.

`clustering_hist`/`clustering_plot`'s currently active scores (built by `buildAnalysisMap`, `src/clustering_dt_common.h`) are `TRKPTZ`, `WAVES`, `JET_T_REFINED`, `TEST_MISAS`, `WAVES_MISCL`, `WAVES_MISAS`, plus `HGTD`.

Comparison PDF groups produced by `clustering_plot` from histograms written by `clustering_hist` (`makeComparisonPlots`, one set per `KEYS = {pu_frac, fjet, ftrack, hs_track, truthjets}`):
- `theplot_*.pdf` — HGTD vs TRKPTZ vs TEST_MISAS
- `waves_oracle_*.pdf` — HGTD vs TRKPTZ vs WAVES vs WAVES_MISAS
- `hgtd_trkptz_waves_*.pdf` — HGTD vs TRKPTZ vs WAVES (no oracles)

Plus inclusive resolution/pull plots (`inclusive/inclusivereso_*`,
`inclusive/inclusivepull_*`, log/lin) and shape-comparison plots
(`inclusive/shape_comparison*`) over `{HGTD, TRKPTZ, WAVES, WAVES_MISAS}`
(resolution set) / `{TRKPTZ, WAVES, WAVES_MISAS}` (shape set).

## event_processing.h Internals

Key structures and entry points:

- **`EventResult`**:
  ```cpp
  struct EventResult {
    int    code          = -1;   // -1: rejected by selection; 0: normal; 3: MISAS fail
    double time          = -1.0; // TRKPTZ-selected cluster time
    int    nFwdHS        =  0;   // n forward HS tracks (3σ counting step)
    bool   trkptzPass    = false;
    bool   misasInDenom  = false; // event entered the TEST_MISAS denominator
    bool   misasPass     = false; // in denominator AND passed
    double clusQuality   =  0.0;  // combined quality of the TRKPTZ-selected cluster
  };
  ```
  `processEventData` returns `EventResult`; early-exit returns `{}` (code=-1).

- **Section E** (auxiliary collections): an `applyFilter` lambda dispatches on `TrackFilterType`, then a `SCORE_REGISTRY` loop calls `clusterTracksInTime` for every score where `s.buildsCollection()` is true, populating `unordered_map<int, vector<Cluster>> auxCollections`. No named per-score collection variables.

- **`selectClusters`**: takes the main clusters, tracks, branch, the active `analyses` map, and `auxCollections`; returns `unordered_map<int,Cluster>` keyed on `Score::id`. Handles the main-collection scores via the all-scores `chooseCluster` overload, aux-collection scores via TRKPTZ selection, and `HGTD`/`HGTD_SORT`/`CONE_BDT` via their own dedicated logic.

## plotting_utilities.h Notes

- **`moneyPlot` / `plotWithLegend`**: the `plotWithLegend` lambda accepts optional legend coordinates `(legX1, legY1, legX2, legY2)` so each call can position the legend independently (bottom-right for efficiency plots, top-right for resolution plots).
- **TEfficiency rendering**: `TEfficiency` objects are converted to `TGraphAsymmErrors` with x-errors explicitly set to 0 at construction — this avoids ROOT's repaint-time x-error restoration in `PaintGraph`. All efficiency plots are drawn this way.
- **Fit quality guard**: `plotPostProcessing` skips slices with too few entries and only calls `fillEach` when the relative uncertainty on the core σ is `< 50%`. Bins failing these guards are left empty rather than showing garbage fit results.
- **y-axis labels**: efficiency-style plots use **"Core Fraction"**; the pass window is `PASS_SIGMA = 60 ps`.
- **`AtlasStyle.h`**: calls `atlasStyle->SetErrorX(0)` globally (effective for TH1/TGraphErrors; TEfficiency uses the conversion approach above).
- **Histogram directory registration**: `clustering_hist.cxx`, `rpt_v5_hist.cxx` (and the legacy `clustering_dt.cxx`/`rpt_v5.cxx`) call `TH1::AddDirectory(kFALSE)` early in `main()`, before any histogram construction — required since all four parallelize their event loop (see Multicore below): every worker thread's lazily-built histogram set is identically named, and without this call they'd race on `gDirectory`'s registry. `PlotObj`/`AnalysisObj`'s own constructors still don't manage `gDirectory` themselves; they rely on the caller having set this first. The `*_plot` executables don't need this — they never build histograms concurrently, just load into ones already booked single-threaded.
- **`AnalysisObj`/`PlotObj::saveTo`/`loadFrom`** (adjacent to `mergeFrom`/`plotPostProcessing`): round-trip the raw, event-loop-filled histograms through `src/histogram_io.h`'s `HistWriter`/`HistReader`. Mirror `mergeFrom`'s exact member list by construction — if a new histogram member is added to `AnalysisObj`, add it to `mergeFrom`, `saveTo`, and `loadFrom` together. See "Histogramming vs. Plotting Split" above for the full design.

## Physics Findings

Results as of 2026-08-20, from `rpt_v5` (grid samples, |Deta| >= 0 and
m_jj >= 500 GeV) and `failure_decomposition` (local VBF). Numbers are pileup
rejection (1/mistag) at fixed hard-scatter efficiency unless stated.

### Where the loss actually lives

`failure_decomposition` assigns each failing event one principal cause:

| cause | share |
|---|---|
| Timing misassignment | **67.4%** |
| Selection failure | 26.5% |
| Misclustering | **3.1%** |
| Track efficiency loss | 0.7% |
| Irreducible | 2.3% |

**The clustering is not the bottleneck.** Grouping all truth-HS tracks into
one perfect cluster almost never rescues a lost event, so the good cluster is
essentially always present and simply not chosen. Selection dominates in
sparse events (52.9% at 0–1 forward HS tracks) and that is where nearly all
failures live (78% fail rate at 0–1 tracks, ~1% above 14); timing
misassignment dominates in track-rich ones.

The scenario ladder says the same thing independently. VBF forward >40 GeV at
0.90 efficiency: `zonly 65.0 -> waves 87.1 -> waves_ideal 106.2 -> truth
111.3`. So **idealising times and selection together is worth ~22%, and the
clustering leaves only ~5% more on the table**. dijet agrees (63.9 vs 68.0,
~6%). Two independent measurements plus the 3.1% above.

Caveat: part of the residual `waves_ideal -> truth` gap is not a clustering
error at all — truth is handed a 10 ps vertex time, while even a perfectly
chosen cluster must estimate t0 from a handful of 30 ps tracks. Separating
those is open.

### The gain is concentrated in the topology

ITk -> WAVeS at 0.80 efficiency:

| selection | VBF | dijet | Z+jets |
|---|---|---|---|
| Inclusive forward | x1.3 | x1.2 | x1.2 |
| R1 (both tags forward) | **x1.87** | **x2.53** | x1.19 |
| R2 (fwd PU + central HS) | x1.17 | x1.18 | x1.04 |

The inclusive number understates HGTD roughly two-fold, because it averages
over jets where timing is irrelevant: pileup jets with no PV0 tracks to begin
with, and HS jets with ample track support.

### R2 fails because WAVeS is captured by the fake

The score ordering **inverts** between regions, consistently across all three
samples — WAVeS wins R1, the jet-agnostic TRKPTZ wins R2 (VBF x1.32 vs x1.17;
dijet x1.44 vs x1.18; Z+jets x1.09 vs x1.04).

WAVeS weights tracks by proximity to forward jets. In R2 *the only forward jet
is the pileup one*, so WAVeS up-weights exactly the tracks of the fake vertex,
drags t0 toward it, and the fake then passes its own gate. The fix is either
regional (use TRKPTZ in R2) or structural (compute t0 from tracks outside
jets, so no single jet can capture the score).

### Timing discriminates by severity, not frequency

Per forward jet, VBF, WAVeS gate vs ITk-only:

| jet | n | raised | lowered | -> zero | mean loss |
|---|---|---|---|---|---|
| HS | 56,342 | **0** | 26.5% | 0.6% | 18.4% |
| PU | 15,307 | **0** | 7.3% | 6.0% | **90.4%** |

The gate fires on HS jets more often than PU ones, but a touched PU jet loses
90% of its R_pT against 18% for an HS jet. Most PU jets are untouched because
they had nothing to remove (they pass the z-association at 0.8%). "Raised" is
0 by construction: `applyTimeGate` returns a SUBSET, so R_pT can only fall.

### The efficiency ceiling is the binding constraint

A jet the gate empties goes to R_pT = 0 and becomes indistinguishable from
pileup, so every timed scenario has a LOWER maximum reachable HS efficiency
than ITk-only (VBF 30–40: 93.8% -> 91.0%; Z+jets: 87.6% -> 82.4%). This, not
resolution, is why curves collapse above ~0.93 and why ratio panels fall back
to 1. Proposed fix, NOT implemented: fall back to the jet's ITk-only R_pT when
the gate empties it — the jet-level analogue of TDR self-tagging.

### Calibration is correct but near-irrelevant

The pull denominator is sqrt(f^2 sigma_vtx^2 + sigma_trk^2) with sigma_vtx
~9 ps against sigma_trk ~26 ps, so a **17% change in inflation moves the
effective gate width by 3%**. Worth having right; not worth another iteration.
Offsets are 0.6–1.5 ps (negligible). Non-Gaussian tails outside a 150 ps core
are 14–15% (VBF), 16–18% (dijet), 21–25% (Z+jets) — the tail, not the core, is
what a Gaussian pull cannot handle.

Z+jets is limited by t0 AVAILABILITY, not by the algorithm: Athena produces a
valid vertex time for 80.1% of VBF events, 73.6% of dijet, and **16.5%** of
Z+jets.

### The z-association is mistuned for TRKPTZ, and irrelevant for WAVeS

From `assoc_study_hist`/`assoc_study_plot`, local VBF, 47,535 selected events,
11 association rules evaluated against the same events (2026-08-26). Full
write-up: `share/docs/Z_ASSOCIATION_STUDY.md`.

The incumbent `significance < 3.0` keeps a track list that is **68% pileup**,
and it is the second-worst of the eleven rules for TRKPTZ:

| score | best rule | vs. incumbent |
|---|---|---|
| TRKPTZ | `getNewDzpara × 1.0` | **−1.58 ps core σ** (15.95 → 14.37), **+0.92 %pt core fraction** (90.38 → 91.30) |
| WAVeS  | `getNewDzpara × 1.2` | +0.13 ps, +0.10 %pt — i.e. nothing |

Both families improve as they tighten and then turn over (significance at
~2.0, the parameterization at ~1.0), so there is a real optimum and the
incumbent is not near it. Simply passing the existing `--dzpara` flag
(`getNewDzpara × 1.4`) already captures ~3/4 of the gain with no code change.
The TRKPTZ gain is concentrated at low multiplicity — σ(nHS≤5) improves
29.3 → 25.3 ps (−13%), which is exactly where the failure decomposition says
the failures live.

**The two scores want opposite things, for a mechanical reason.**
`Cluster::calculateTime` recomputes the WAVeS time from only the tracks within
ΔR < 0.4 of a forward jet, and that in-jet filter already does the
pileup-shedding a tighter z-association would do. The two mechanisms are
redundant, so tightening can only cost WAVeS hard-scatter tracks — and it does,
monotonically (`dzpara × 0.6`: +1.86 ps, −1.83 %pt, and worse at *every*
multiplicity, not just sparse events). Expect any future selector that
recomputes its time from a track subset to behave the same way, and to want a
LOOSE association.

**~59% of WAVeS's advantage over TRKPTZ is just the z-association.** The
core-fraction gap is 1.38 %pt at the incumbent rule but 0.56 %pt when each
score gets its own best rule. On core σ the ordering inverts outright: TRKPTZ
at `dzpara × 1.0` (14.37 ps) beats WAVeS at any rule (best 14.95 ps). WAVeS's
residual edge is in core FRACTION, not core width — i.e. it suppresses gross
misassignment, consistent with the failure decomposition, not resolution.

**The parameterization is a better shape, not merely a tighter cut.** Compared
at matched track purity the two families separate, and the gap widens the
tighter the working point: at ~31.5% purity the parameterization is worth +1.0
points of HS recall and the two perform identically; at ~43.7% purity it is
worth **+5.6 points of recall at higher purity**, and TRKPTZ gains 0.5 ps and
0.6 %pt from the shape alone. Likely cause: `sqrt(var_z0)` is a per-track
covariance that does not capture the η dependence of the true z0 spread, and
`Z0_VAR_INFLATION = 1.15²` is independent evidence it is mis-estimated — a
single global factor cannot repair an η-dependent mis-modelling.

**Union with ghost association: tested and rejected.** The idea was that ghost
association, being derived independently of z, might recover HS tracks the
z-cut misses. OR-ing it in costs **20–24 %pt of core fraction** for BOTH
scores (TRKPTZ 91.30 → 68.47, WAVeS 91.80 → 72.98, on a `dzpara × 1.0` base).

The failure mode is worth remembering because it is not the obvious one.
**Core σ barely moves (~1 ps) while core fraction collapses** — so the added
tracks are not dragging the HS cluster's time, they are causing the wrong
cluster to be *selected*. Confirmed differentially: every z-only rule
saturates at ~99.9% core fraction once there are ≥20 forward HS tracks, while
the union plateaus at **84%** and never saturates. It picks wrong one time in
six even in the most track-rich events.

Why: the added pileup is **time-random** (mean |Δt| to the truth HS vertex of
224 ps, against 229 ps for a track drawn at random from the event — measured
in `util/scratch/ghost_assoc_diag.cxx`), but sits inside forward jets and so
carries real pT. It builds *decoy clusters* rather than poisoning the HS one,
feeding the selection channel that the failure decomposition already names as
dominant. WAVeS is hit as hard as TRKPTZ here despite being immune to every
z-only change, because jet-proximate tracks are exactly what its score
up-weights.

Also established while testing this: **ghost association is numerically the
same as `filterTracksInJets` at ΔR < 0.4** (1,142,857 tracks / 24.44% HS vs
1,139,170 / 24.54%), so it is not new information relative to what WAVeS
already uses.

**Not yet acted on.** `MAX_NSIGMA` is unchanged: this has only been run on the
local VBF sample. Z+jets in particular should be checked first, since Athena
produces a valid vertex time for only 16.5% of its events. Rerunning is cheap
(~35 s hist + ~15 s plot on the full local sample; condor via `--sample=`).

### Reference-study replication

Central region at the reference configuration (`util/myJet_ana_fr.C`, present
in the main working tree but not committed — it covers both regions via
commented-out eta lines and carries the z0-significance cut as a commented
alternative): the
**timing improvement factor reproduces** (ours 1.63 vs their ~1.6 at 0.800)
but both absolute rejections sit ~2.3x low. A uniform offset cannot be caused
by the time gate, so the open question is the pileup population or its
normalisation, not the timing model. Jet collection (EMTopo) and the
isHS/isPU definitions are confirmed to match.

### Methodological warnings

- **Never rank a selection oracle by cluster purity.** `calcPurity` is a
  pT-weighted fraction of tracks from truth vertex 0 — it measures where a
  track CAME FROM, not whether its time is right. A lone HS track carrying a
  mis-assigned HGTD hit forms a 100%-pure cluster sitting at a wrong time, so
  max-purity selects exactly the dominant failure mode (HS jets driven to
  R_pT = 0 rise 3.0% -> 11.8%). Rank on `|t_cluster - t_truth|` instead.
- **When selecting display candidates for a failure mode, require the
  failure.** Gating on low HS timing purity alone yields events where the
  tracks are badly mis-timed and the algorithm recovers the right time anyway.
- **A flat 30 ps track smear is WORSE than HGTD's real timing** (Initial
  scenario: 35/25/21 ps for 1/2/3 hits, ~26 ps effective), so the reference's
  30 ps scenario is not an upper bound on what the simulated detector
  delivers.

## Development Notes

- **Namespace**: all utilities live in `namespace MyUtl`
- **Debug mode**: `MyUtl::DEBUG` (`src/clustering_constants.h`) for verbose per-event logging
- **Multicore**: `clustering_hist.cxx`/`rpt_v5_hist.cxx` (and the legacy `clustering_dt.cxx`/`rpt_v5.cxx`) parallelize their per-event loop via `ROOT::TTreeProcessorMT` (`--threads=N`, `MyUtl::resolveThreads()` in `src/sample_config.h`, default `min(hardware_concurrency(), 8)`) — `ROOT::EnableImplicitMT()` alone only parallelises ROOT-internal operations, not the event loop itself. Per-thread accumulation state (`AnalysisObj` map / RpT `ThreadState`) is built lazily once per worker thread and merged into one result after the loop (`AnalysisObj::mergeFrom()` in `src/plotting_utilities.h`). See `util/ttreemt_prototype.cxx` for the standalone prototype this was validated against, and the "TTreeProcessorMT prototype" note below on a real thread-safety bug found via stress-testing. The `*_plot` executables are single-threaded — no event loop, no ntuple access, just histogram load + fit/derive + draw.
- **Constants**: modify analysis cuts in `src/clustering_constants.h`; always rebuild after header changes
- **Plot styling**: ROOT `TColor`/`TCanvas`; ATLAS style via `AtlasStyle.h`/`AtlasLabels.h`
- **Event loop**: `TTreeProcessorMT::Process()` hands each task its own `TTreeReader`; a fresh `BranchPointerWrapper` is bound to it per invocation (cheap — pointer setup only)
- **File identification**: inside the `TTreeProcessorMT` lambda, `reader.GetTree()->GetCurrentFile()->GetName()` gives the current file path, and `reader.GetTree()->GetReadEntry()` gives the local (0-indexed) entry number directly — both needed to reconstruct an event display command. This differs from the old sequential loop's `chain.GetReadEntry() - chain.GetChainOffset()` trick, which needed an outer `TChain` not available inside the lambda; a `TChain`-bound *sequential* `TTreeReader`'s `GetTree()->GetReadEntry()` returns the chain-global entry number, not local, which is exactly why that trick existed.
- **TTreeProcessorMT thread-safety gotcha**: `AnalysisObj`'s constructor calls `SetFillColorAlpha()`, which lazily registers a new `TColor` into ROOT's *global* color table the first time a given (color, alpha) pair is used — not thread-safe. Worker threads racing to build their first `AnalysisObj`/`ThreadState` concurrently corrupted the heap (intermittent `SIGABRT`, found via stress-testing, not a one-off). Fixed by serializing each thread's entire lazy-init construction (not just registry bookkeeping) behind its mutex — see the per-thread registry pattern in both executables' `main()`.
