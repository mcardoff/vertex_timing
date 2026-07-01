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
├── src/                           # Core library headers + main analysis executable
│   ├── clustering_constants.h    # Analysis constants, Score struct + SCORE_REGISTRY, ClusteringMethod enum
│   ├── clustering_functions.h    # Clustering algorithms (clusterTracksInTime, chooseCluster, ...)
│   ├── clustering_includes.h     # Common ROOT/C++ includes
│   ├── clustering_structs.h      # BranchPointerWrapper, Cluster (incl. updateScores/calcFeatures), EventCounts
│   ├── sample_config.h           # --sample=vbf|zjets|dijet resolution (default: local ../../ntuple-hgtd/)
│   ├── event_processing.h        # setupChain, getAssociatedTracks, selectClusters, processEventData
│   ├── plotting_utilities.h      # PlotObj / AnalysisObj + all plotting helpers (moneyPlot, inclusivePlot, ...)
│   ├── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
│   └── clustering_dt.cxx         # Main analysis executable (HGTD scenario only; see Main Executables)
├── util/                          # Diagnostic / report / ML-data-export executables
│   ├── failure_decomposition.cxx # Per-event principal-cause failure classification (see section below)
│   ├── export_training_data.cxx  # Exports per-cluster CSV rows (8 features + purity/labels) for future DNN retraining
│   ├── generate_rpt.cxx          # RpT-style timing/discrimination report (older version)
│   ├── rpt_v2.cxx .. rpt_v5.cxx  # RpT discrimination study iterations; rpt_v5 uses the WAVeS score
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
├── condor/                       # Condor submit templates + run_analysis.sh wrapper
├── CMakeLists.txt
└── build/                        # Build output (created by `cmake ..`)
```

Note: several directories/files described in earlier versions of this doc
(`vertex_time.cxx`, `check_jet_sorting.cxx`, `test_vbf_*.cxx`, `DBSCAN.ipynb`,
`figs/`, `old/`, `src/ml_model.h`, and most of the `util/` diagnostics and
parameter-sweep executables) do not exist in the current tree — this section
reflects what's actually present.

## Main Executables

### Primary Analysis
- **clustering_dt.cxx** (`src/`): Main clustering analysis. Only the HGTD
  timing scenario is currently built/processed (`Scenario::IDEAL_RES`/
  `IDEAL_EFF` and `buildAnalysisMap`'s branching for them still exist for
  future reuse, but `main()` no longer instantiates or iterates those maps).
  Active scores: `TRKPTZ`, `WAVES`, `JET_T_REFINED`, `TEST_MISAS`,
  `WAVES_MISCL`, `WAVES_MISAS`, plus `HGTD`. Iterative clustering
  (`ClusteringMethod::ITERATIVE`, `DIST_CUT_CONE = 3.0σ`) is the default main
  collection. Produces the `theplot_*`, `waves_oracle_*`, and
  `hgtd_trkptz_waves_*` comparison PDF groups (see Score Registry section),
  plus inclusive resolution/pull plots and shape-comparison plots.

### Diagnostics & Failure Analysis
- **failure_decomposition.cxx** (`util/`): Assigns every TRKPTZ-failing event
  to a single principal cause via a fixed precedence ladder (see Failure
  Decomposition section below).

### RpT / Report Studies
- **generate_rpt.cxx**, **rpt_v2.cxx**..**rpt_v5.cxx** (`util/`): Successive
  iterations of an RpT (jet-level pileup-discrimination) study comparing
  z-only, real-HGTD, and algorithm-selected (TRKPTZ or, as of `rpt_v5`,
  WAVeS) timing gates. `rpt_v5` is the most current version.

### ML Data Export
- **export_training_data.cxx** (`util/`): Exports one CSV row per
  cluster/passing event (8 features + purity/labels) mirroring
  `clustering_dt`'s event/track selection exactly. Produces training data for
  a possible future DNN, but no model currently consumes it in `src/` — see
  Machine Learning Model section.

### Utilities
- **maxpt_jet_test.cxx** (`util/`): Diagnostic comparing WAVeS/TRKPTZ pass
  rates under different `maxTrkPt` caps (current `MAX_TRACK_PT`, leading-jet
  pT, subleading-jet pT). Pure diagnostic, changes nothing.
- **track_dt.cxx** (`util/`): Inclusive per-track HGTD-time-vs-truth-particle
  residual, no event selection.
- **runHGTD_Clustering.cxx** (`python/`): Interactive per-event clustering.
  Output is parsed by `event_display.py` to produce a 2-page PDF.

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
```
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

Internally calls `runHGTD_Clustering.cxx` and parses its stdout to produce a 2-page PDF (R-Z geometry view + time histogram). Script lives in **`python/`** (not `share/scripts/`).

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
# Full analysis
cd build && ./clustering_dt

# Failure decomposition
cd build && ./failure_decomposition

# Event display for a specific event
cd python && python event_display.py --file_num 000009 --event_num 1856 --extra_time 0.0

# RpT study (WAVeS-based)
cd build && ./rpt_v5

# Rebuild after header change
cd build && make
```

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

`clustering_dt`'s currently active scores (built by `buildAnalysisMap`) are `TRKPTZ`, `WAVES`, `JET_T_REFINED`, `TEST_MISAS`, `WAVES_MISCL`, `WAVES_MISAS`, plus `HGTD`.

Comparison PDF groups produced by `clustering_dt` (`makeComparisonPlots`, one set per `KEYS = {pu_frac, fjet, ftrack, hs_track, truthjets}`):
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
- **Histogram directory registration**: `PlotObj`/`AnalysisObj` construction does not currently call `TH1::AddDirectory(kFALSE)` or otherwise manage `gDirectory` — fine for the current single-threaded event loop, but anything that constructs multiple `AnalysisObj` instances concurrently (e.g. per-thread duplication) must set `TH1::AddDirectory(kFALSE)` first to avoid racing on ROOT's global directory registry.

## Development Notes

- **Namespace**: all utilities live in `namespace MyUtl`
- **Debug mode**: `MyUtl::DEBUG` (`src/clustering_constants.h`) for verbose per-event logging
- **Multicore**: `ROOT::EnableImplicitMT()` is called in `main()` but only parallelises ROOT-internal operations — the `TTreeReader::Next()` event loop itself remains sequential. See `util/ttreemt_prototype.cxx` for an exploratory (not production-wired) `TTreeProcessorMT` prototype that parallelises the per-event loop itself.
- **Constants**: modify analysis cuts in `src/clustering_constants.h`; always rebuild after header changes
- **Plot styling**: ROOT `TColor`/`TCanvas`; ATLAS style via `AtlasStyle.h`/`AtlasLabels.h`
- **Event loop**: `TTreeReader` + `BranchPointerWrapper` for branch access
- **File identification**: within the sequential event loop, `branch.reader.GetTree()->GetCurrentFile()->GetName()` gives the current file path, and `chain.GetReadEntry() - chain.GetChainOffset()` gives the local (0-indexed) entry number within that file — both needed to reconstruct an event display command
