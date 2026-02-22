# Vertex Timing Analysis

HGTD vertex timing reconstruction and clustering analysis for ATLAS.

## Directory Structure

```
vertex_timing/
├── src/                          # Analysis executables and core library headers
│   ├── clustering_dt.cxx         # Main clustering analysis (primary entry point)
│   ├── hgtd_matching.cxx         # Error analysis: quantifies failure modes
│   ├── extract_error_metrics.cxx # Post-processor: generates failure mode summary table
│   ├── test_ml_model.cxx         # ML model smoke test
│   ├── generate_rpt.cxx          # ROC curve / timing distribution report generator
│   ├── evaluate_ml_lowtrack.cxx  # ML performance in low-track-count events
│   ├── export_training_data.cxx  # Exports cluster features to CSV for DNN retraining
│   ├── distcut_sweep.cxx         # Parameter sweep over cone clustering distance cuts
│   ├── clustering_constants.h    # Analysis constants, cuts, Score enum
│   ├── clustering_functions.h    # Core clustering and scoring algorithms
│   ├── clustering_includes.h     # Aggregated ROOT/Boost/STL includes
│   ├── clustering_structs.h      # BranchPointerWrapper and Cluster data structures
│   ├── event_processing.h        # High-level event loop orchestration
│   ├── plotting_utilities.h      # Gaussian fitting and plot generation utilities
│   ├── ml_model.h                # Standalone C++ neural network (8→128→64→32→1)
│   ├── suppress_stdout.h         # RAII stdout suppressor for noisy ROOT messages
│   └── json.hpp                  # Header-only JSON parser for model weight loading
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
| `hgtd_matching` | `src/hgtd_matching.cxx` | Error analysis tool. Quantifies three failure modes: wrong HGTD track-time matching, wrong hard-scatter cluster selection, and insufficient forward tracks. Produces 11 diagnostic plots and a ROOT histogram file. |
| `extract_error_metrics` | `src/extract_error_metrics.cxx` | Reads histograms from `hgtd_matching_analysis.root` and prints a quantitative breakdown of each failure mode's contribution to overall inefficiency. |
| `test_ml_model` | `src/test_ml_model.cxx` | Loads the neural network from `share/models/model_weights.json` and runs inference on dummy inputs to verify the C++ ML implementation is functioning correctly. |
| `generate_rpt` | `src/generate_rpt.cxx` | Generates timing distribution plots and two-sample ROC curves comparing pT-weighted and pT+time-weighted track timing for hard-scatter vs. pileup separation. |
| `test_vbf_selection` | `test_vbf_selection.cxx` | Applies the full VBF signal region cut sequence (jet multiplicity, pT, η, Δφ, Δη, dijet mass, forward jets) and reports per-step pass rates. |
| `evaluate_ml_lowtrack` | `src/evaluate_ml_lowtrack.cxx` | Evaluates ML model performance specifically on low-track-multiplicity events (N_HS ≤ 4), producing score distributions and purity/efficiency curves in this challenging regime. |
| `export_training_data` | `src/export_training_data.cxx` | Exports per-cluster feature vectors (delta_z, uncertainties, sumpt, label) to a CSV file for use in retraining the neural network. Mirrors the event selection of `clustering_dt`. |
| `distcut_sweep` | `src/distcut_sweep.cxx` | Sweeps over cone clustering distance cut values and plots efficiency, timing resolution (from a double-Gaussian fit), and mean cluster count as a function of cut threshold. |

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

## Disabled Event Selection Cuts

The following cuts are commented out in `src/event_processing.h` and can be re-enabled for specific studies:

- **VBF signal region** (`branch->passVBFSignalRegion()`): Applies VBF H→invisible selection cuts. Disabled to run over the full inclusive sample; re-enable for VBF-specific efficiency studies.
- **Forward HS track requirement** (`branch->pass_forward_hs_tracks(nForwardTrack_HS)`): Requires a minimum number of hard-scatter tracks in the forward region. Useful for studies of track multiplicity dependence; disabled by default to retain low-track events in the denominator.

---
