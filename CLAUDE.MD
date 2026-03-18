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
├── src/                           # Analysis executables and core library headers
│   ├── clustering_constants.h    # Analysis constants, Score enum, ClusteringMethod enum
│   ├── clustering_functions.h    # Clustering algorithms; ML model loaded here (line 294)
│   ├── clustering_includes.h     # Common ROOT/C++ includes
│   ├── clustering_structs.h      # Data structures; ML feature extraction (line 328)
│   ├── event_processing.h        # Event loop: setupChain, getAssociatedTracks, processEvent
│   ├── plotting_utilities.h      # Visualization tools (overlay histogram, money plots)
│   ├── sweep_utilities.h         # Shared helpers for 1-D parameter sweep executables
│   ├── ml_model.h                # Standalone C++ neural network
│   ├── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
│   ├── json.hpp                  # JSON parser for model weights
│   │
│   ├── clustering_dt.cxx         # Main analysis (3 timing scenarios × multiple scores)
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
│   └── dist_refine_sweep.cxx     # Sweep: DIST_CUT_REFINE (REFINED iterative distance)
├── runHGTD_Clustering.cxx        # Interactive event-display clustering script (root dir)
├── vertex_time.cxx               # Vertex timing distribution analysis with jet selection
├── check_jet_sorting.cxx         # Utility: verify jets are pT-sorted in ntuples
├── test_vbf_integration.cxx      # Tests VBF cut integration
├── test_vbf_selection.cxx        # VBF signal region selection demo
├── share/
│   ├── models/                   # ML model files (model.onnx, model_weights.json,
│   │                             #   normalization_params.json)
│   ├── scripts/                  # Python utilities
│   │   ├── event_display.py      # Event visualisation (calls runHGTD_Clustering.cxx)
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
- **clustering_dt.cxx** (`src/`): Main clustering analysis comparing HGTD times, ideal resolution, and efficiency scenarios. Evaluates `Score::HGTD`, `TRKPTZ`, `REFINED`, and optionally `TESTML`.
- **vertex_time.cxx**: Vertex timing distribution analysis with forward jet selection.

### Diagnostics & Failure Analysis
- **failure_decomposition.cxx** (`src/`): Classifies every TRKPTZ failure into four non-exclusive categories (see section below). Primary diagnostic tool.
- **rate_diagnostics.cxx** (`src/`): Selection-agnostic rates — misclustering (purity > 75%?) and timing misassignment (HGTD vs Ideal Res. PASS oracle). Binned tables matching clustering_dt style.
- **hgtd_matching.cxx** (`src/`): Low-level track/cluster matching diagnostics. Quantifies wrong HGTD matching, wrong cluster selection, and insufficient track count failure modes.
- **extract_error_metrics.cxx** (`src/`): Extracts quantitative metrics from hgtd_matching output.

### Parameter Sweeps
All sweep executables live in `src/` and share helpers from `sweep_utilities.h`.
- **distcut_sweep.cxx**: Sweeps `DIST_CUT_CONE` (0.5 σ – 5 σ).
- **mintrackpt_sweep.cxx** / **maxtrackpt_sweep.cxx**: pT threshold sweeps.
- **cone_iter_k_sweep.cxx**: Sweeps `CONE_ITER_K` (top-K cone clusters for REFINED pooling).
- **dist_refine_sweep.cxx**: Sweeps `DIST_CUT_REFINE` (iterative distance cut in REFINED stage).

### ML Tools
- **test_ml_model.cxx**: Validates the loaded neural network model.
- **evaluate_ml_lowtrack.cxx**: Evaluates NN scores specifically on ≤4 HS-track events.
- **export_training_data.cxx**: Exports 8-feature CSV rows + labels for DNN retraining.

### Utilities
- **runHGTD_Clustering.cxx** (root dir): Interactive per-event clustering. Runs CONE at 3.0 σ then applies the REFINED stage (top-`CONE_ITER_K` cone clusters re-clustered with ITERATIVE at `DIST_CUT_REFINE`). To revert to plain cone output comment out the `===== REFINED STAGE =====` block. Output format is consumed by `event_display.py`.
- **generate_rpt.cxx**: Report generation with timing distributions.
- **check_jet_sorting.cxx**: Verifies jet pT ordering in ntuples.

## REFINED Clustering Algorithm

REFINED is a two-stage approach implemented in `src/event_processing.h` and mirrored in `runHGTD_Clustering.cxx`:

1. **Stage 1 (Cone)**: Standard cone clustering at `DIST_CUT_CONE = 3.0 σ` on all forward tracks.
2. **Stage 2 (Iterative refinement)**: Filter cone clusters to those with `≥ MIN_CLUSTER_TRACKS` tracks, sort by TRKPTZ score descending, pool the track indices of the top `CONE_ITER_K = 3` clusters, then re-cluster the pooled tracks using `ClusteringMethod::ITERATIVE` at the tighter `DIST_CUT_REFINE = 2.0 σ`. The ITERATIVE method updates the centroid after each merge, preventing PU track absorption that plagues CONE in spread events.

Relevant constants in `src/clustering_constants.h`:
```cpp
DIST_CUT_CONE   = 3.0   // Stage 1 cone distance
DIST_CUT_REFINE = 2.0   // Stage 2 iterative distance
CONE_ITER_K     = 3     // Top-K cone clusters to pool
MIN_CLUSTER_TRACKS = 2  // Minimum tracks to qualify a cluster
```

Score enum: `Score::REFINED` — selected by TRKPTZ among the refined sub-clusters.

## Failure Decomposition

`src/failure_decomposition.cxx` classifies every TRKPTZ-failing event into four non-exclusive categories (an event can satisfy multiple):

| # | Name | Condition |
|---|------|-----------|
| 1 | Selection failure | Another qualifying HGTD cluster passes the efficiency window |
| 2 | Timing misassignment | Re-running with ideal timing resolution fixes it |
| 3 | Misclustering | No HGTD cluster has purity > 75% **OR** the REFINED algorithm fixes it |
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
MIN_TRACK_PT     = 0.5   // GeV — track selection floor
MAX_TRACK_PT     = 30.0  // GeV — track selection ceiling
DIST_CUT_CONE    = 3.0   // σ — cone clustering distance
DIST_CUT_REFINE  = 2.0   // σ — REFINED iterative distance
CONE_ITER_K      = 3     // top-K cone clusters for REFINED pooling
IDEAL_TRACK_RES  = 10.0  // ps — flat per-track resolution for ideal scenarios
PASS_SIGMA       = 20.0  // pass threshold for efficiency oracle
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
python share/scripts/event_display.py --file_num <NUM> --event_num <ENTRY> --extra_time <dt>
```

Internally calls `runHGTD_Clustering.cxx` and parses its stdout to produce a 2-page PDF (R-Z geometry view + time histogram).

## Machine Learning Model

- **Architecture**: 8 → 128 → 64 → 32 → 1 (ReLU hidden, Sigmoid output; 11 521 parameters)
- **Features**: delta_z, delta_z_resunits, cluster_z_sigma, cluster_d0, cluster_d0_sigma, cluster_qOverP, cluster_qOverP_sigma, cluster_sumpt
- **Loading**: Lazy static init in `src/clustering_functions.h:294` — once per run, ~0.01 ms per inference
- **Score**: `Score::TESTML`
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
python share/scripts/event_display.py --file_num 000009 --event_num 1856 --extra_time 0.0

# Parameter sweep
cd build && ./dist_refine_sweep

# Rebuild after header change
cd build && make
```

## Development Notes

- **Namespace**: all utilities live in `namespace MyUtl`
- **Debug mode**: `#define debug true` at top of executable for verbose logging
- **Multicore**: `ROOT::EnableImplicitMT()` parallelises histogram filling; `TTreeReader::Next()` loop remains sequential
- **Constants**: modify analysis cuts in `src/clustering_constants.h`; always rebuild after header changes
- **Plot styling**: ROOT `TColor`/`TCanvas`; ATLAS style via `AtlasStyle.h`/`AtlasLabels.h`
- **Event loop**: `TTreeReader` + `BranchPointerWrapper` for branch access
- **File identification**: within a TChain loop, `chain.GetCurrentFile()->GetName()` gives the file path and `chain.GetTree()->GetReadEntry()` gives the local (0-indexed) entry number — both needed to reconstruct an event display command
