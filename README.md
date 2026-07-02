# Vertex Timing Analysis

HGTD vertex timing reconstruction and clustering analysis for ATLAS.

## Directory Structure

```
vertex_timing/
├── src/                          # Core library headers + main analysis executable
│   ├── clustering_dt.cxx         # Main clustering analysis (primary entry point)
│   ├── clustering_constants.h    # Analysis constants, cuts, Score enum
│   ├── clustering_functions.h    # Core clustering and scoring algorithms
│   ├── clustering_includes.h     # Aggregated ROOT/Boost/STL includes
│   ├── clustering_structs.h      # BranchPointerWrapper and Cluster data structures
│   ├── event_processing.h        # High-level event loop orchestration
│   ├── plotting_utilities.h      # Gaussian fitting and plot generation utilities
│   └── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
├── util/                         # Diagnostic and utility executables
│   ├── sweep_utilities.h         # Shared helpers for 1-D parameter sweep executables
│   ├── generate_rpt.cxx          # ROC curve / timing distribution report generator
│   ├── export_training_data.cxx  # Exports cluster features to CSV for DNN retraining
│   └── failure_decomposition.cxx # Failure mode categorisation + pie chart
├── share/
│   ├── models/
│   │   ├── model.onnx                  # Original ONNX model (8 input features)
│   │   ├── model_weights.json          # Exported weights for C++ inference
│   │   └── normalization_params.json   # Feature mean/scale for input normalization
│   ├── scripts/                        # Legacy/utility scripts (see share/scripts below)
│   └── docs/                           # Reference documentation (ML integration, VBF, etc.)
├── python/                       # Interactive/visualisation scripts (run from within this dir)
│   ├── runHGTD_Clustering.cxx    # Single-event clustering wrapper / debug tool
│   ├── event_display.py          # Interactive event visualization (vertices, tracks, jets)
│   ├── event_tinkering.py        # Quick event inspection and jet printing via uproot
│   └── clustering_animation.py   # Animates the iterative clustering algorithm for one event
├── condor/                       # HTCondor grid submission templates
│   ├── run_analysis.sh           # Generic job executable (any --sample-aware target)
│   ├── clustering_dt.sub         # Submit file: one job per sample, for clustering_dt
│   ├── rpt_v5.sub                # Submit file: one job per sample, for rpt_v5
│   └── logs/                     # condor stdout/stderr/log output (created empty)
├── CMakeLists.txt                # Build configuration
├── build/                        # CMake build output directory
├── figs/                         # Output plots (PDF)
└── CLAUDE.md                     # Full project documentation for Claude Code
```

---

## Build System Executables

The following targets are defined in `CMakeLists.txt` and built in the `build/` directory:

| Executable | Source | Description |
|---|---|---|
| `clustering_dt` | `src/clustering_dt.cxx` | Main analysis. Runs timing-based vertex reconstruction across three scenarios (real HGTD, ideal resolution, ideal efficiency) and several scoring algorithms. Supports `--sample=vbf\|zjets\|dijet` to switch input ntuple, energy-label annotation, and output directory. Generates comparison PDFs in `figs/` (or `vbf/`, `zjets/`, `dijet/`). |
| `failure_decomposition` | `util/failure_decomposition.cxx` | Classifies every TRKPTZ-failing event into four failure categories (selection, timing misassignment, misclustering, track efficiency). Produces efficiency curves + pie chart. |
| `generate_rpt` | `util/generate_rpt.cxx` | Generates timing distribution plots and ROC curves comparing pT-weighted and pT+time-weighted track timing for hard-scatter vs. pileup separation. |
| `rpt_v2` | `util/rpt_v2.cxx` | RpT discrimination study using the main-analysis event selection (z-only / HGTD / TRKPTZ scenarios, with purity and HS-timing oracle gates). |
| `rpt_v3` | `util/rpt_v3.cxx` | Jet-level RpT analysis with no event selection, matching the paper's approach. |
| `rpt_v4` | `util/rpt_v4.cxx` | Jet-level RpT analysis with a 30–40 GeV / >40 GeV jet-pT split. |
| `rpt_v5` | `util/rpt_v5.cxx` | Jet-level RpT analysis using the WAVeS-selected cluster time; supports `--sample=vbf\|zjets\|dijet`. |
| `export_training_data` | `util/export_training_data.cxx` | Exports per-cluster feature vectors to CSV for DNN retraining. |
| `maxpt_jet_test` | `util/maxpt_jet_test.cxx` | Compares WAVeS/TRKPTZ efficiency under different track pT upper-cut choices. |
| `track_dt` | `util/track_dt.cxx` | Inclusive per-track timing residual (HGTD track time − truth particle production-vertex time), no event/track selection. |

---

## Root-Level Scripts

### `python/` Scripts (run from within `python/`)

| File | Description |
|---|---|
| `runHGTD_Clustering.cxx` | Single-event clustering wrapper for interactive debugging. Loads one event, runs iterative-time clustering, and prints cluster details (time, z, score, track list, pass/fail). Output is consumed by `event_display.py`. |
| `event_display.py` | Generates event display plots showing reconstructed and truth vertices, tracks, timing information, cluster assignments, and jets. Accepts `--file_num`, `--event_num`, and `--extra_time` arguments. |
| `event_tinkering.py` | Quick inspection tool that loads ROOT ntuples via uproot and prints jet information (truth HS jets, in-time/out-of-time pileup jets, reconstructed topo jets) for a given event. |
| `clustering_animation.py` | Animates the iterative HGTD clustering algorithm step-by-step for one event. |

---

## `share/scripts/` Contents

These are legacy or utility scripts not built by CMake. The ONNX/DNN tooling
below is orphaned — the `TEST_ML` score and its C++ inference engine
(`src/ml_model.h`) have been removed from the active codebase, so these
scripts no longer feed into anything `clustering_dt` reads:

| File | Description |
|---|---|
| `generate_timeplot.cxx` | ROOT macro for generating timing resolution plots from pre-processed histogram files. |
| `generate_onefile_timeplot.cxx` | Variant of `generate_timeplot.cxx` operating on a single input ROOT file. |
| `testing_separation.cxx` | Tests vertex separation algorithms using common clustering utilities. |
| `test_onnx_model.C` | ROOT macro that verifies ONNX model loading via ROOT's TMVA SOFIE interface. (Orphaned — see note above.) |
| `inspect_onnx_model.py` | Inspects the ONNX model graph, prints layer shapes, and can export weights to `model_weights.json`. (Orphaned — see note above.) |
| `export_onnx_model.py` | All-in-one pipeline: converts a TF SavedModel → ONNX → `model_weights.json`. (Orphaned — see note above.) |
| `model_evaluation_helper.py` | Guide and reference code for evaluating the ML model in Python and comparing against C++ inference. (Orphaned — see note above.) |

---

## Quick Start

```bash
# Build all executables
cd build
cmake ..
make

# Run main analysis (outputs to figs/, or vbf|zjets|dijet with --sample)
./clustering_dt

# Run failure decomposition
./failure_decomposition

# Visualize an event (script lives in python/)
cd python && python event_display.py --file_num 000009 --event_num 1856 --extra_time 0.0
```

---

## Grid Submission (condor)

Templates for running `clustering_dt` / `rpt_v5` on the UChicago AF HTCondor pool live in `condor/`:

| File | Purpose |
|---|---|
| `run_analysis.sh` | Generic job executable — `run_analysis.sh <executable> <sample>`. Sources the ATLAS/LCG environment (`atlasLocalSetup.sh` + `lsetup "root <version>"`) then runs the transferred `<executable>` with `--sample=<sample>` (or no flag if `<sample>` is `default`). |
| `clustering_dt.sub` | Submit file: one job per `--sample` value (`vbf`, `zjets`, `dijet`) for `clustering_dt`. |
| `rpt_v5.sub` | Same, for `rpt_v5`. Submitted separately (its own `condor_submit` call) so it runs as an independent job cluster from `clustering_dt.sub`. |

Execute hosts don't share a filesystem with the submit host on this pool (confirmed via a `FileSystemDomain` match failure when `should_transfer_files = NO` was tried), so both `.sub` files stage the prebuilt executable in via `transfer_input_files` and stage the output directory (`vbf/`, `zjets/`, or `dijet/`) back out via `transfer_output_files`. The ntuple inputs are read via the absolute `/data/mcardiff/exotic_superntuples/...` path baked into `--sample=`, so they don't need transferring — only the tiny compiled binary does.

```bash
# Build first — condor does not rebuild on the execute host
cd build && cmake .. && make clustering_dt rpt_v5

# Submit (from inside condor/, so relative paths resolve correctly)
cd ../condor
condor_submit -dump dryrun.ad clustering_dt.sub   # dry-run: validate without queuing
condor_submit clustering_dt.sub
condor_submit rpt_v5.sub                          # separate job cluster

# Monitor
condor_q
condor_q -better-analyze <ClusterId>.<ProcId>     # if a job is idle/held
```

Logs land in `condor/logs/<sample>.<ClusterId>.<ProcId>.{out,err,log}` (prefixed `rpt_v5_` for the rpt_v5 cluster).

**Caveat:** `CMakeLists.txt` compiles with `-march=native`, tuned to whichever machine ran `make`. If execute hosts have a different CPU microarchitecture than the build host, jobs can crash with an "Illegal instruction" (SIGILL) rather than a normal error — not a bug in the job scripts.

## Scoring Algorithms

Each run of `clustering_dt` evaluates several scoring algorithms simultaneously. Scores are defined as `struct Score` instances in `src/clustering_constants.h` and carry their own metadata (name, threshold, flags). The active set is:

| Short name | Long name | Description |
|---|---|---|
| `HGTD` | HGTD Algorithm | Selects the cluster chosen by the HGTD detector's own timing algorithm. Own collection; only active in the real-HGTD scenario. |
| `TRKPT` | ΣpT | Scalar sum of track pT. |
| `TRKPTZ` | ΣpTe^{−\|Δz\|} | Primary score. Weights pT by exp(−\|Δz\|/σ_z); tracks closer in z contribute more. |
| `PASS` | Pass Cluster | First cluster passing quality cuts. |
| `CONE` | Cone | Legacy cone-clustering baseline. |
| `CONE_BDT` | Cone (BDT) | Cone clustering with a TMVA BDT selector (currently disabled — stub only). |
| `HGTD_SORT` | HGTD BDT (pT-sorted) | pT-sorted simultaneous clustering with a TMVA BDT selector (currently disabled — stub only). |
| `FILTJET` | Filter Tracks in Jets | Re-clusters using only tracks within ΔR < 0.4 of a jet, then TRKPTZ. |
| `TEST_HS` | TRKPTZ (HS only) | Cone clustering restricted to truth-HS-linked tracks only. |
| `WAVES` | WAVeS Score | Σ pT·pT_jet/max(ΔR,floor) weighted by exp(−1.5\|Δz\|); in-jet timing refinement drops absorbed PU tracks. |
| `JET_T_REFINED` | WAVeS 2σ t Re-clustering | Own collection: clusters only jet-proximate tracks at the tighter 2σ iterative distance, selected by TRKPTZ. |
| `TEST_MISAS` (MISAS) | TRKPTZ (no t misassign.) | **Oracle**: fills only when 100% of all fwd HS tracks have \|pull\| < 3σ (no timing misassignment). |
| `WAVES_MISCL` | WAVeS (pure clust) | **Oracle**: WAVeS-selected cluster, gated on cluster purity. |
| `WAVES_MISAS` | WAVeS (no t misassign.) | **Oracle**: WAVeS-selected cluster, gated on event-level HS timing purity. |

Scores with `requiresPurity = true` (oracle scores) gate their denominator — events outside the gate are excluded entirely. All oracle functions return `0.0f` when no qualifying HS tracks exist, so zero-HS-track events are excluded from all oracle denominators.

---

## Removed Scoring Algorithms

Several scores have been removed from the codebase over time as they were superseded or found unreliable:

- **`T_REFINED`, `Z_REFINED`, `ZT_REFINED`, `ZT_ITER`** — earlier two-stage re-clustering/refinement variants, superseded by the WAVeS-based scores (`WAVES`, `JET_T_REFINED`) above.
- **`TEST_MISCL`, `PERF_EVT`** — TRKPTZ purity/combined oracles, superseded by the `WAVES_MISCL`/`WAVES_MISAS` oracle pair.
- **`TEST_ML`** — a DNN selection score (8→128→64→32→1 neural network). The C++ inference engine (`src/ml_model.h`) and its JSON weight loader (`src/json.hpp`) have been deleted along with the score.
- **`CALO60` / `CALO90`, `JUST60` / `JUST90`, `FILT60` / `FILT90`** — experimental scores based on an independent calorimeter-derived timing estimate to gate, replace, or pre-filter the HGTD track timing. The calorimeter time was a Gaussian-weighted average of track times within a fixed cone around the calorimeter cluster centroid; in practice the estimator proved unreliable and all three approaches degraded efficiency without a compensating improvement in purity or resolution. The supporting infrastructure (`filterCaloTracks`, the calorimeter time variables in `chooseCluster`, and the JUST synthetic-cluster branches) has been deleted along with the score definitions.

---

## Event Selection

Every event must pass `branch->passBasicCuts()` (minimum jet count, reco/truth HS vertex within `MAX_VTX_DZ`) and `branch->passJetPtCut()` (both in `src/clustering_structs.h`), the latter requiring:
- At least `MIN_PASSPT_JETS` reco jets above `MIN_JET_PT`, at least `MIN_PASSETA_JETS` of which are in the forward HGTD acceptance — no truth-HS jet matching required.
- The VBS candidate pair — the opposite-hemisphere, pT-passing jet pair with the largest invariant mass `m_jj` (`calcBestVbsDeltaEta()`) — must have an η separation ≥ `VBS_JET_D_ETA`.

`branch->passForwardHsTracks(nForwardHSTrack)` (`src/clustering_structs.h`) is an available but currently-unused testing-only helper that gates on a minimum forward hard-scatter track count; not applied in the primary analysis.

---
