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
│   ├── ml_model.h                # Standalone C++ neural network (8→128→64→32→1)
│   ├── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
│   └── json.hpp                  # Header-only JSON parser for model weight loading
├── util/                         # Diagnostic, sweep, and ML utility executables
│   ├── sweep_utilities.h         # Shared helpers for 1-D parameter sweep executables
│   ├── hgtd_matching.cxx         # Error analysis: quantifies failure modes
│   ├── extract_error_metrics.cxx # Post-processor: generates failure mode summary table
│   ├── test_ml_model.cxx         # ML model smoke test
│   ├── generate_rpt.cxx          # ROC curve / timing distribution report generator
│   ├── evaluate_ml_lowtrack.cxx  # ML performance in low-track-count events
│   ├── export_training_data.cxx  # Exports cluster features to CSV for DNN retraining
│   ├── failure_decomposition.cxx # Failure mode categorisation + pie chart
│   ├── rate_diagnostics.cxx      # Selection-agnostic misclustering/misassignment rates
│   ├── distcut_sweep.cxx         # Sweep: cone distance cut
│   ├── mintrackpt_sweep.cxx      # Sweep: minimum track pT
│   ├── maxtrackpt_sweep.cxx      # Sweep: maximum track pT
│   ├── cone_iter_k_sweep.cxx     # Sweep: CONE_ITER_K
│   └── dist_refine_sweep.cxx     # Sweep: DIST_CUT_REFINE
├── share/
│   ├── models/
│   │   ├── model.onnx                  # Original ONNX model (8 input features)
│   │   ├── model_weights.json          # Exported weights for C++ inference
│   │   └── normalization_params.json   # Feature mean/scale for input normalization
│   ├── scripts/                        # Legacy/utility scripts (see share/scripts below)
│   └── docs/                           # Reference documentation (ML integration, VBF, etc.)
├── vertex_time.cxx               # Standalone ROOT script: vertex timing distributions
├── runHGTD_Clustering.cxx        # Single-event clustering wrapper / debug tool
├── check_jet_sorting.cxx         # Diagnostic: verifies jets are pT-sorted in ntuples
├── test_vbf_selection.cxx        # Validates VBF cut sequence and prints pass rates
├── test_vbf_integration.cxx      # Integration test for VBF cuts in event processing
├── event_display.py              # Interactive event visualization (vertices, tracks, jets)
├── event_tinkering.py            # Quick event inspection and jet printing via uproot
├── analyze_pileup_removal.py     # Python analysis of pileup removal cut effects
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
| `clustering_dt` | `src/clustering_dt.cxx` | Main analysis. Runs timing-based vertex reconstruction across three scenarios (real HGTD, ideal resolution, ideal efficiency) and multiple scoring algorithms (HGTD, TRKPTZ, TESTML). Generates comparison plots in `figs/`. |
| `hgtd_matching` | `util/hgtd_matching.cxx` | Error analysis tool. Quantifies three failure modes: wrong HGTD track-time matching, wrong hard-scatter cluster selection, and insufficient forward tracks. Produces 11 diagnostic plots and a ROOT histogram file. |
| `extract_error_metrics` | `util/extract_error_metrics.cxx` | Reads histograms from `hgtd_matching_analysis.root` and prints a quantitative breakdown of each failure mode's contribution to overall inefficiency. |
| `test_ml_model` | `util/test_ml_model.cxx` | Loads the neural network from `share/models/model_weights.json` and runs inference on dummy inputs to verify the C++ ML implementation is functioning correctly. |
| `generate_rpt` | `util/generate_rpt.cxx` | Generates timing distribution plots and two-sample ROC curves comparing pT-weighted and pT+time-weighted track timing for hard-scatter vs. pileup separation. |
| `test_vbf_selection` | `test_vbf_selection.cxx` | Applies the full VBF signal region cut sequence (jet multiplicity, pT, η, Δφ, Δη, dijet mass, forward jets) and reports per-step pass rates. |
| `evaluate_ml_lowtrack` | `util/evaluate_ml_lowtrack.cxx` | Evaluates ML model performance specifically on low-track-multiplicity events (N_HS ≤ 4), producing score distributions and purity/efficiency curves in this challenging regime. |
| `export_training_data` | `util/export_training_data.cxx` | Exports per-cluster feature vectors (delta_z, uncertainties, sumpt, label) to a CSV file for use in retraining the neural network. Mirrors the event selection of `clustering_dt`. |
| `distcut_sweep` | `util/distcut_sweep.cxx` | Sweeps over cone clustering distance cut values and plots efficiency, timing resolution (from a double-Gaussian fit), and mean cluster count as a function of cut threshold. |

---

## Root-Level Scripts

### C++ Scripts (run with ROOT or compile separately)

| File | Description |
|---|---|
| `vertex_time.cxx` | Vertex timing distribution analysis. Selects forward jets and the hard-scatter vertex, fits timing residuals with single- and double-Gaussian models, and produces resolution vs. forward jet multiplicity plots. |
| `runHGTD_Clustering.cxx` | Single-event clustering wrapper for interactive debugging. Loads one event, runs cone-time clustering with configurable parameters, and prints cluster details (time, z, score, track list, pass/fail). |
| `check_jet_sorting.cxx` | Diagnostic tool. Checks whether jets in input ROOT files are stored in descending-pT order and prints a per-event summary plus aggregate pass rate. |
| `test_vbf_selection.cxx` | Validates the VBF selection cut chain against the ntuple data and reports per-cut pass rates. (Also a CMake build target — see above.) |
| `test_vbf_integration.cxx` | Integration test checking that `processEventData()` correctly rejects non-VBF events, with printed pass-rate statistics. |

### Python Scripts

| File | Description |
|---|---|
| `event_display.py` | Generates event display plots showing reconstructed and truth vertices, tracks, timing information, cluster assignments, and jets. Accepts `--file_num`, `--event_num`, and `--extra_time` arguments. |
| `event_tinkering.py` | Quick inspection tool that loads ROOT ntuples via uproot and prints jet information (truth HS jets, in-time/out-of-time pileup jets, reconstructed topo jets) for a given event. |
| `analyze_pileup_removal.py` | Analyzes the effect of the 3σ pileup removal cut on track composition, measuring how many hard-scatter and pileup tracks are removed at each stage. |

---

## `share/scripts/` Contents

These are legacy or utility scripts not built by CMake:

| File | Description |
|---|---|
| `generate_timeplot.cxx` | ROOT macro for generating timing resolution plots from pre-processed histogram files. |
| `generate_onefile_timeplot.cxx` | Variant of `generate_timeplot.cxx` operating on a single input ROOT file. |
| `testing_separation.cxx` | Tests vertex separation algorithms using common clustering utilities. |
| `test_onnx_model.C` | ROOT macro that verifies ONNX model loading via ROOT's TMVA SOFIE interface. |
| `inspect_onnx_model.py` | Inspects the ONNX model graph, prints layer shapes, and can export weights to `model_weights.json`. |
| `export_onnx_model.py` | Helper for converting a trained Keras/TensorFlow model to ONNX format. |
| `model_evaluation_helper.py` | Guide and reference code for evaluating the ML model in Python and comparing against C++ inference. |

---

## Quick Start

```bash
# Build all executables
cd build
cmake ..
make

# Run main analysis (outputs to figs/)
./clustering_dt

# Run error analysis
./hgtd_matching
./extract_error_metrics   # prints quantitative breakdown

# Visualize an event
python ../event_display.py --file_num 0 --event_num 42
```

## Scoring Algorithms

Each run of `clustering_dt` evaluates several scoring algorithms simultaneously. Scores are defined as `struct Score` instances in `src/clustering_constants.h` and carry their own metadata (name, threshold, flags). The active set is:

| Short name | Long name | Description |
|---|---|---|
| `HGTD` | HGTD Algorithm | Selects the cluster chosen by the HGTD detector's own timing algorithm. Uses the HGTD-seeded cluster collection; only active in the real-HGTD scenario. |
| `TRKPT` | ΣpT | Scores each cluster by the scalar sum of track pT. Selects the highest-scoring cluster. |
| `TRKPTZ` | ΣpT·exp(−\|Δz\|) | Like `TRKPT` but weights each track's pT by exp(−\|Δz\|/σ_z) so tracks closer to the vertex in z contribute more. |
| `PASS` | Pass Cluster | Selects the first cluster that passes all quality requirements, without using a score. |
| `FILTJET` | Filter Tracks in Jets | Re-clusters using only tracks within ΔR < 0.4 of a reconstructed jet before scoring with TRKPTZ. Only active in the real-HGTD scenario. |
| `TESTML` | DNN Selection | Scores clusters with a neural network (8→128→64→32→1, ReLU/sigmoid). Applies a score threshold of 0.3. Weights loaded from `share/models/model_weights.json`. |
| `MISCL` | TRKPTZ (pure clusters) | Runs TRKPTZ but only fills histograms for events where the selected cluster is "pure" (≥80% hard-scatter tracks). Isolates the intrinsic timing resolution from cluster-selection errors. |
| `HGTD_SORT` | HGTD BDT (pT-sorted) | Like `HGTD` but sorts tracks by pT before clustering (order-sensitive cone algorithm) and applies a TMVA BDT post-selection with threshold 0.3. |

Scores that set `usesOwnCollection = true` (`HGTD`, `HGTD_SORT`) bypass the main cone-clustering output and select from their own pre-built collections. Scores with `requiresPurity = true` (`MISCL`) only fill histograms when the selected cluster passes a purity check.

---

## Removed Scoring Algorithms

The following experimental scores were removed from the codebase. They were all based on using an independent calorimeter-derived timing estimate to gate or replace the HGTD track timing. The calorimeter time was computed as a Gaussian-weighted average of track times within a fixed cone around the calorimeter cluster centroid. In practice the estimator proved unreliable and all three approaches degraded efficiency without a compensating improvement in purity or resolution.

| Short name | How it worked |
|---|---|
| `CALO60` / `CALO90` | Applied the standard `TRKPTZ` score, but vetoed cluster selection if the cluster's time differed from the calorimeter time estimate by more than 60 ps / 90 ps. The cluster collection and z reconstruction were unchanged; only the selection decision was affected by the calorimeter gate. |
| `JUST60` / `JUST90` | Constructed a synthetic cluster whose assigned time came entirely from the calorimeter estimate. Only tracks whose times fell within 60 ps / 90 ps of the calorimeter time were included in the synthetic cluster. Unlike the FILT variants, the underlying cone-clustering and z reconstruction used the full track set; only the reported cluster time was replaced. |
| `FILT60` / `FILT90` | Called `filterCaloTracks` before cone clustering to remove any track whose time lay more than 60 ps / 90 ps from the calorimeter estimate. Because the filtering happened before clustering, both the z position and the track composition of the resulting clusters were affected, not just the timing. |

The supporting infrastructure (`filterCaloTracks`, the calorimeter time variables in `chooseCluster`, and the JUST synthetic-cluster branches) has been deleted along with the score definitions.

---

## Disabled Event Selection Cuts

The following cuts are commented out in `src/event_processing.h` and can be re-enabled for specific studies:

- **VBF signal region** (`branch->passVBFSignalRegion()`): Applies VBF H→invisible selection cuts. Disabled to run over the full inclusive sample; re-enable for VBF-specific efficiency studies.
- **Forward HS track requirement** (`branch->pass_forward_hs_tracks(nForwardTrack_HS)`): Requires a minimum number of hard-scatter tracks in the forward region. Useful for studies of track multiplicity dependence; disabled by default to retain low-track events in the denominator.

---
