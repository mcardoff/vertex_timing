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
│   ├── clustering_constants.h    # Analysis constants, Score enum, ClusteringMethod enum
│   ├── clustering_functions.h    # Clustering algorithms; ML model loaded here (line 294)
│   ├── clustering_includes.h     # Common ROOT/C++ includes
│   ├── clustering_structs.h      # Data structures; ML feature extraction (line 328)
│   ├── event_processing.h        # Event loop: setupChain, getAssociatedTracks, processEvent
│   ├── plotting_utilities.h      # Visualization tools (overlay histogram, money plots)
│   ├── ml_model.h                # Standalone C++ neural network
│   ├── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
│   ├── json.hpp                  # JSON parser for model weights
│   └── clustering_dt.cxx         # Main analysis (3 timing scenarios × multiple scores)
├── util/                          # Diagnostic, sweep, and ML utility executables
│   ├── sweep_utilities.h         # Shared helpers for 1-D parameter sweep executables
│   ├── failure_decomposition.cxx # Failure mode categorisation + pie chart
│   ├── rate_diagnostics.cxx      # Selection-agnostic misclustering/misassignment rates
│   ├── hgtd_matching.cxx         # Low-level track/cluster matching diagnostics
│   ├── extract_error_metrics.cxx # Extracts metrics from hgtd_matching output
│   ├── generate_rpt.cxx          # Timing distribution report
│   ├── test_ml_model.cxx         # ML model validation
│   ├── evaluate_ml_lowtrack.cxx  # NN score analysis for low-track-count events
│   ├── export_training_data.cxx  # Exports CSV rows for DNN retraining (8 features)
│   ├── distcut_sweep.cxx         # Sweep: cone distance cut
│   ├── mintrackpt_sweep.cxx      # Sweep: minimum track pT
│   ├── maxtrackpt_sweep.cxx      # Sweep: maximum track pT
│   ├── cone_iter_k_sweep.cxx     # Sweep: CONE_ITER_K (REFINED top-K parameter)
│   ├── dist_refine_sweep.cxx     # Sweep: DIST_CUT_T_REFINED (REFINED iterative distance)
│   └── nsigma_sweep.cxx          # Sweep: distance cut expressed in σ units
├── event_display.py              # Event visualisation — calls runHGTD_Clustering.cxx, produces 2-page PDF
├── runHGTD_Clustering.cxx        # Interactive per-event clustering (root dir)
├── vertex_time.cxx               # Vertex timing distribution analysis with jet selection
├── check_jet_sorting.cxx         # Utility: verify jets are pT-sorted in ntuples
├── test_vbf_integration.cxx      # Tests VBF cut integration
├── test_vbf_selection.cxx        # VBF signal region selection demo
├── share/
│   ├── models/                   # ML model files (model.onnx, model_weights.json,
│   │                             #   normalization_params.json)
│   ├── scripts/                  # Python utilities
│   │   ├── inspect_onnx_model.py
│   │   ├── export_onnx_model.py
│   │   └── model_evaluation_helper.py
│   └── docs/                     # Design docs and ML integration guides
├── CMakeLists.txt
├── build/                        # Build output
├── figs/                         # Output plots (PDF)
├── old/                          # Legacy implementations
└── DBSCAN.ipynb
```

## Main Executables

### Primary Analysis
- **clustering_dt.cxx** (`src/`): Main clustering analysis comparing HGTD times, ideal resolution, and efficiency scenarios. Active scores: `TRKPTZ`, `T_REFINED`, `Z_REFINED`, `ZT_REFINED`, `ZT_ITER`, `TEST_MISCL`, `TEST_MISAS`, `TEST_CTIME`, `PERF_EVT`, `PERF_CLT` (plus `HGTD` in the HGTD scenario). Iterative clustering is the default for all TRKPTZ-based scores. Generates six groups of comparison PDFs (see Score Registry section). At end of run prints up to 20 random pass/fail event display commands for `nFwdHS == 2` events (non-deterministic seed).

- **vertex_time.cxx**: Vertex timing distribution analysis with forward jet selection.

### Diagnostics & Failure Analysis
- **failure_decomposition.cxx** (`util/`): Classifies every TRKPTZ failure into four non-exclusive categories (see section below). Primary diagnostic tool.
- **rate_diagnostics.cxx** (`util/`): Selection-agnostic rates — misclustering (purity > 75%?) and timing misassignment (HGTD vs Ideal Res. PASS oracle). Binned tables matching clustering_dt style.
- **hgtd_matching.cxx** (`util/`): Low-level track/cluster matching diagnostics. Quantifies wrong HGTD matching, wrong cluster selection, and insufficient track count failure modes.
- **extract_error_metrics.cxx** (`util/`): Extracts quantitative metrics from hgtd_matching output.

### Parameter Sweeps
All sweep executables live in `util/` and share helpers from `util/sweep_utilities.h`.
- **distcut_sweep.cxx**: Sweeps `DIST_CUT_CONE` (0.5 σ – 5 σ).
- **mintrackpt_sweep.cxx** / **maxtrackpt_sweep.cxx**: pT threshold sweeps.
- **cone_iter_k_sweep.cxx**: Sweeps `CONE_ITER_K` (top-K clusters for T_REFINED pooling).
- **dist_refine_sweep.cxx**: Sweeps `DIST_CUT_T_REFINED` (iterative distance cut in T_REFINED stage).
- **nsigma_sweep.cxx**: Sweeps the distance cut expressed in σ units.

### ML Tools
- **test_ml_model.cxx**: Validates the loaded neural network model.
- **evaluate_ml_lowtrack.cxx**: Evaluates NN scores specifically on ≤4 HS-track events.
- **export_training_data.cxx**: Exports 8-feature CSV rows + labels for DNN retraining.

### Utilities
- **runHGTD_Clustering.cxx** (root dir): Interactive per-event clustering. Uses `ClusteringMethod::ITERATIVE` at `DIST_CUT_CONE = 3.0 σ` by default. The T_REFINED stage block is present but commented out. Output is parsed by `event_display.py` (project root) to produce a 2-page PDF.
- **generate_rpt.cxx**: Report generation with timing distributions.
- **check_jet_sorting.cxx**: Verifies jet pT ordering in ntuples.

## T_REFINED Clustering Algorithm

`Score::T_REFINED` is a two-stage approach implemented in `src/event_processing.h` and mirrored in `runHGTD_Clustering.cxx`:

1. **Stage 1 (Iterative)**: Iterative nearest-neighbour clustering at `DIST_CUT_CONE = 3.0 σ` on all forward tracks. This is now also the default for all TRKPTZ-based scores.
2. **Stage 2 (Timing refinement)**: Filter clusters to those with `≥ MIN_CLUSTER_TRACKS` tracks, sort by TRKPTZ score descending, pool the track indices of the top `CONE_ITER_K = 3` clusters, then re-cluster the pooled tracks using `ClusteringMethod::ITERATIVE` at the tighter `DIST_CUT_T_REFINED = 2.0 σ`. The ITERATIVE method updates the centroid after each merge, preventing PU track absorption.

`Score::CONE` uses legacy cone clustering (`ClusteringMethod::CONE`) at `DIST_CUT_CONE = 3.0 σ` and is preserved as a comparison baseline.

Relevant constants in `src/clustering_constants.h`:
```cpp
DIST_CUT_CONE      = 3.0   // Main iterative distance (also used for CONE comparison)
DIST_CUT_T_REFINED = 2.0   // T_REFINED / ZT_REFINED stage 2 iterative distance
CONE_ITER_K        = 3     // Top-K clusters to pool for T_REFINED
MIN_CLUSTER_TRACKS = 0     // Minimum tracks to qualify a cluster
```

Score enum: `Score::T_REFINED` (id=4) — selected by TRKPTZ among the refined sub-clusters.

## Failure Decomposition

`src/failure_decomposition.cxx` classifies every TRKPTZ-failing event into four non-exclusive categories (an event can satisfy multiple):

| # | Name | Condition |
|---|------|-----------|
| 1 | Selection failure | Another qualifying HGTD cluster passes the efficiency window |
| 2 | Timing misassignment | Re-running with ideal timing resolution fixes it |
| 3 | Misclustering | No HGTD cluster has purity > 75% **OR** the T_REFINED algorithm fixes it |
| 4 | Track efficiency loss | Ideal resolution + ideal efficiency fixes it (but ideal resolution alone does not) |

Categories 2 and 4 are mutually exclusive by construction.
`cat_none` catches events satisfying none of the four — true fundamental failures.

**Output per run:**
- Console: summary table with per-category counts/fractions, per-bin breakdown vs n Forward HS Tracks, and a sanity check that `|A∪B∪C∪D| + |None| = n_fail`.
- PDF page 1: TEfficiency curves per category, normalised to **failing events** in each bin (`h_fail` denominator). The overall failure rate (`n_fail / n_total`) is overlaid as a grey dashed reference.
- PDF page 2: Pie chart of exclusive Venn decomposition (12 slices summing to `n_fail`); a `printf` sanity check verifies the sum before drawing.

**cat_none debug output**: When `cat_none` fires, a ready-to-run event display command is printed:
```
python3 event_display.py --file_num <NUM> --event_num <ENTRY> --extra_time 0.00
```
followed by the full cluster dump (time, score, purity, tracks, passes?).

## Key Analysis Parameters

From `src/clustering_constants.h`:
```cpp
MIN_TRACK_PT          = 1.0   // GeV — track selection floor
MAX_TRACK_PT          = 30.0  // GeV — track selection ceiling
DIST_CUT_CONE         = 3.0   // σ — main iterative clustering distance (also CONE baseline)
DIST_CUT_T_REFINED    = 2.0   // σ — T_REFINED / ZT_REFINED stage 2 iterative distance
CONE_ITER_K           = 3     // top-K clusters for T_REFINED pooling
IDEAL_TRACK_RES       = 1.0   // ps — flat per-track resolution for ideal scenarios
PASS_SIGMA            = 100.0 // ps — ±window for "Core Fraction" efficiency oracle
HS_TIMING_QUALITY_CUT = 3.0   // σ — |pull| threshold for oracle timing purity gates
```

## Build System

### Requirements
- CMake ≥ 3.10
- C++17 compiler
- ROOT (RIO, Core, Gpad, Graf, Hist, Tree)
- Boost (filesystem)

### Build Commands
```bash
cd build
cmake ..
make [<target>]   # e.g. make failure_decomposition
```

### Compiler Flags
`-O3 -march=native -DNDEBUG` plus thin LTO.

## Data Format

ROOT ntuples with the following key branches:

| Group | Branches |
|-------|----------|
| Jets | `AntiKt4EMTopoJets_pt`, `_eta`, `_track_idx` |
| Vertices | `TruthVtx_z/time`, `RecoVtx_z/time/timeRes/hasValidTime` |
| Tracks | `Track_z0/pt/eta/time/timeRes/hasValidTime/var_z0/truthVtx_idx/truthPart_idx` |
| Particles | `TruthPart_prodVtx_z` |

## Event Display

```bash
python event_display.py --file_num <NUM> --event_num <ENTRY> --extra_time <dt>
```

Internally calls `runHGTD_Clustering.cxx` and parses its stdout to produce a 2-page PDF (R-Z geometry view + time histogram). Script lives at the **project root** (not `share/scripts/`).

## Machine Learning Model

- **Architecture**: 8 → 128 → 64 → 32 → 1 (ReLU hidden, Sigmoid output; 11 521 parameters)
- **Features**: delta_z, delta_z_resunits, cluster_z_sigma, cluster_d0, cluster_d0_sigma, cluster_qOverP, cluster_qOverP_sigma, cluster_sumpt
- **Loading**: Lazy static init in `src/clustering_functions.h:294` — once per run, ~0.01 ms per inference
- **Score**: `Score::TEST_ML`
- **Files**: `share/models/model_weights.json`, `model.onnx`, `normalization_params.json`

### Retraining Workflow
1. Train in Python, export to ONNX: `tf2onnx.convert.from_keras()`
2. Save normalization params (`scaler.mean_`, `scaler.scale_`)
3. `python share/scripts/inspect_onnx_model.py --model share/models/model.onnx --export-weights`
4. Copy new files to `share/models/`
5. `cd build && make`

## Common Workflows

```bash
# Full analysis
cd build && ./clustering_dt

# Failure decomposition
cd build && ./failure_decomposition

# Event display for a specific event
python event_display.py --file_num 000009 --event_num 1856 --extra_time 0.0

# Parameter sweep
cd build && ./dist_refine_sweep

# Rebuild after header change
cd build && make
```

## Score Registry

All scores are defined as `static const Score` members of the `Score` struct in `src/clustering_constants.h` and listed in `SCORE_REGISTRY`. Each carries its own collection spec: `distCut`, `method` (`CONE`/`ITERATIVE`), `useZ0`, `filter` (`TrackFilterType`: `ALL`, `JET`, `Z0_TVA`, `HS_ONLY`), plus `buildsCollection()` to indicate it manages its own cluster collection.

| id | Short | Long name | Notes |
|----|-------|-----------|-------|
| 0  | HGTD  | HGTD Algorithm | Own collection; HGTD-seeded |
| 1  | TRKPT | Σ pT | Plain pT sum |
| 2  | TRKPTZ | Σ pTe^{−\|Δz\|} | Primary score |
| 3  | PASS  | Pass Cluster | First cluster passing quality cuts |
| 4  | T_REFINED | 2σ t Reclustering | Two-stage; ITERATIVE at DIST_CUT_T_REFINED |
| 5  | Z_REFINED | 2σ z Refinement | z-based refinement |
| 6  | ZT_REFINED | 2σ z+t Reclustering | Two-stage with Z0_TVA filter |
| 7  | CONE  | Cone | Legacy cone baseline |
| 8  | CONE_BDT | Cone (BDT) | Cone + BDT threshold 0.3 |
| 9  | FILTJET | Filter Tracks in Jets | JET filter; tracks within ΔR < 0.4 of jet |
| 10 | HGTD_SORT | HGTD BDT (pT-sorted) | pT-sorted + BDT threshold |
| 11 | TEST_ML | DNN Selection | Neural net score threshold 0.3 |
| 12 | MISCL | TRKPTZ (pure clust) | Oracle: cluster purity > 75% |
| 13 | MISAS | TRKPTZ (no t misassign.) | Oracle: event-level HS timing purity = 100% |
| 14 | TEST_HS | TRKPTZ (HS only) | HS_ONLY filter; truth-matched tracks only |
| 15 | CTIME | TRKPTZ (clust t pure) | Oracle: cluster-level HS timing purity = 100% |
| 16 | ZT_ITER | 2D (z₀,t) Iterative | `useZ0=true`; 2D Mahalanobis distance √((Δt/σt)²+(Δz₀/σz₀)²) |
| 17 | PERF_EVT | TRKPTZ (Pure Clust+Evt Times) | Combined oracle: MISCL ∧ MISAS |
| 18 | PERF_CLT | TRKPTZ (Pure Clust.+Clust. Times) | Combined oracle: MISCL ∧ CTIME |

Oracle scores (`requiresPurity=true`) gate their denominator: events outside the gate are not counted. **Important:** both `calcHSTimingPurity` (event-level, used by MISAS) and `calcClusterHSTimingPurity` (cluster-level, used by CTIME) return `0.0f` — not `1.0f` — when no qualifying HS tracks exist, ensuring that events with zero forward HS tracks are excluded from all oracle denominators.

Money plot groups produced by `clustering_dt`:
- `theplot_*.pdf` — main score comparison
- `timing_oracle_*.pdf` — MISAS vs CTIME
- `perfect_timing_oracle_*.pdf` — PERF_EVT vs PERF_CLT
- `hgtd_trkptz_*.pdf` — HGTD vs TRKPTZ
- `z_refined_*.pdf` — z/zt refinement variants
- `2d_clust_*.pdf` — ZT_ITER comparison

## event_processing.h Internals

Key structures and refactored entry points:

- **`EventResult`** (replaces `std::pair<int,double>`):
  ```cpp
  struct EventResult { int code; double time; int nFwdHS; bool trkptzPass; };
  ```
  `processEventData` returns `EventResult`; early-exit returns `{}` (code=-1).

- **Section E** (auxiliary collections): an `applyFilter` lambda dispatches on `TrackFilterType`, then a `SCORE_REGISTRY` loop calls `clusterTracksInTime` for every score where `s.buildsCollection()` is true, populating `unordered_map<int, vector<Cluster>> auxCollections`. No more named per-score collection variables.

- **`selectClusters`** signature: takes `const unordered_map<int,vector<Cluster>>& auxCollections` (previously 10 named collection parameters).

## plotting_utilities.h Notes

- **`moneyPlot` / `plotWithLegend`**: the `plotWithLegend` lambda accepts optional legend coordinates `(legX1, legY1, legX2, legY2)` so each call can position the legend independently (bottom-right for efficiency plots, top-right for resolution plots).
- **TEfficiency rendering**: `TEfficiency` objects are converted to `TGraphAsymmErrors` with x-errors explicitly set to 0 at construction — this avoids ROOT's repaint-time x-error restoration in `PaintGraph`. All efficiency plots are drawn this way.
- **Fit quality guard**: `plotPostProcessing` skips slices with `< 50` entries (was 6) and only calls `fillEach` when `percErrCore < 0.5` (relative uncertainty on the core σ). Bins failing these guards are left empty rather than showing garbage fit results.
- **y-axis labels**: efficiency-style plots use **"Core Fraction"**; the pass window is `PASS_SIGMA = 100 ps`.
- **`AtlasStyle.h`**: calls `atlasStyle->SetErrorX(0)` globally (effective for TH1/TGraphErrors; TEfficiency uses the conversion approach above).

## Development Notes

- **Namespace**: all utilities live in `namespace MyUtl`
- **Debug mode**: `#define debug true` at top of executable for verbose logging
- **Multicore**: `ROOT::EnableImplicitMT()` parallelises histogram filling; `TTreeReader::Next()` loop remains sequential
- **Constants**: modify analysis cuts in `src/clustering_constants.h`; always rebuild after header changes
- **Plot styling**: ROOT `TColor`/`TCanvas`; ATLAS style via `AtlasStyle.h`/`AtlasLabels.h`
- **Event loop**: `TTreeReader` + `BranchPointerWrapper` for branch access
- **File identification**: within a TChain loop, `chain.GetCurrentFile()->GetName()` gives the file path and `chain.GetTree()->GetReadEntry()` gives the local (0-indexed) entry number — both needed to reconstruct an event display command
