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
├── src/                          # Core library headers
│   ├── clustering_constants.h   # Analysis constants and configuration
│   ├── clustering_functions.h   # Clustering algorithms (ML model loaded here)
│   ├── clustering_includes.h    # Common ROOT/C++ includes
│   ├── clustering_structs.h     # Data structures (ML feature extraction)
│   ├── event_processing.h       # Event loop processing
│   ├── plotting_utilities.h     # Visualization tools
│   ├── ml_model.h               # Standalone C++ neural network
│   └── json.hpp                 # JSON parser for model weights
├── share/                        # Shared resources (NEW!)
│   ├── models/                  # ML model files
│   │   ├── model.onnx          # ONNX model (8 features)
│   │   ├── model_weights.json  # Exported weights for C++ (330 KB)
│   │   └── normalization_params.json  # Feature normalization parameters
│   ├── scripts/                 # Utility scripts
│   │   ├── event_display.py    # Python event visualization
│   │   ├── inspect_onnx_model.py  # Model inspection tool
│   │   ├── export_onnx_model.py   # ONNX export helper
│   │   └── model_evaluation_helper.py  # Model integration guide
│   ├── docs/                    # Documentation
│   │   ├── SUMMARY.md           # Project status overview
│   │   ├── ML_QUICK_START.md    # Quick ML integration guide
│   │   ├── ML_FEATURE_EXTRACTION.md  # Feature calculation guide
│   │   ├── ML_MODEL_INTEGRATION.md   # Complete ML integration
│   │   └── ERROR_ANALYSIS_PLAN.md    # Error analysis methodology
│   └── README.md                # Share directory documentation
├── *.cxx                         # Analysis executables (see below)
├── CMakeLists.txt               # Build configuration
├── build/                       # Build directory
│   ├── clustering_dt            # Main analysis executable
│   ├── hgtd_matching            # Error analysis executable
│   ├── extract_error_metrics    # Metrics extraction
│   ├── test_ml_model            # ML model test
│   ├── error_metrics_summary.txt  # Analysis results
│   └── example_events.txt       # Event display candidates
├── figs/                        # Output plots
├── old/                         # Legacy implementations
└── DBSCAN.ipynb                 # Jupyter notebook for DBSCAN clustering
```

## Main Executables

### Primary Analysis
- **clustering_dt.cxx**: Main clustering analysis comparing HGTD times, ideal resolution, and efficiency scenarios
- **vertex_time.cxx**: Vertex timing distribution analysis with forward jet selection
- **hgtd_matching.cxx**: Comprehensive error analysis quantifying three major failure modes:
  1. Wrong HGTD matching (incorrect track time assignments)
  2. Wrong cluster selection (incorrect hard scatter identification)
  3. Insufficient tracks (too few tracks for reliable clustering)

### Plotting Utilities
- **generate_timeplot.cxx**: Generate timing resolution plots
- **generate_onefile_timeplot.cxx**: Single-file timing plots
- **generate_rpt.cxx**: Report generation with timing distributions

### Support Scripts
- **runHGTD_Clustering.cxx**: Wrapper for running HGTD clustering
- **testing_separation.cxx**: Test vertex separation algorithms

## Key Analysis Parameters

From vertex_time.cxx:
```cpp
min_jets          = 2      // Minimum number of jets per event
min_jetpt         = 30.0   // Minimum jet pT (GeV)
min_abs_eta       = 2.0    // Minimum |eta| for forward vertices
min_abs_track_eta = 2.4    // Minimum |eta| for forward tracks
min_track_pt      = 1.0    // Minimum track pT (GeV)
max_vtx_dz        = 1.0e6  // Maximum vertex z error
max_nsigma        = 3.0    // Track-to-PV association threshold
```

## Build System

### Requirements
- CMake >= 3.10
- C++17 compiler
- ROOT (with RIO, Core, Gpad, Graf, Hist, Tree components)
- Boost (filesystem)

### Build Commands
```bash
cd build
cmake ..
make
```

### Compiler Flags
- `-O3`: Maximum optimization
- `-march=native`: CPU-specific optimization
- `-DNDEBUG`: Disable debug assertions

## Data Format

The analysis reads ROOT ntuples with the following branch structure:

### Jets
- `AntiKt4EMTopoJets_pt`, `AntiKt4EMTopoJets_eta`
- `AntiKt4EMTopoJets_track_idx`: Track association

### Vertices
- `TruthVtx_z`, `TruthVtx_time`: Truth-level vertices
- `RecoVtx_z`, `RecoVtx_time`, `RecoVtx_timeRes`: Reconstructed vertices
- `RecoVtx_hasValidTime`: Timing validity flag

### Tracks
- `Track_z0`, `Track_pt`, `Track_eta`: Kinematics
- `Track_time`, `Track_timeRes`, `Track_hasValidTime`: Timing info
- `Track_var_z0`: z0 variance
- `Track_truthVtx_idx`, `Track_truthPart_idx`: Truth matching

### Particles
- `TruthPart_prodVtx_z`: Production vertex z-position

## Analysis Modes

The clustering_dt.cxx implements three analysis scenarios:

1. **HGTD Times**: Real HGTD detector timing performance
2. **Ideal Resolution**: Perfect timing resolution, real efficiency
3. **Ideal Efficiency**: Perfect timing resolution and efficiency

Each mode evaluates multiple scoring algorithms:
- `HGTD`: Basic HGTD timing
- `TRKPTZ`: Track pT-weighted timing
- `PASS`: Passing clusters
- (Commented: CALO60, CALO90, JUST60, JUST90, FILT60, FILT90, FILTJET)

## Output

### Plots
Generated plots are saved to `figs/` directory in PDF format:
- Vertex timing distributions
- Timing resolution vs. forward jet multiplicity
- Efficiency and purity curves
- 2D correlation plots

### Event Display
Events of interest can be visualized using:
```bash
python share/scripts/event_display.py --file_num <file> --event_num <event> --extra_time <dt>
```

## Machine Learning Model

### Overview
The project integrates a neural network for cluster selection, providing improved hard scatter identification compared to traditional TRKPT/TRKPTZ scoring methods.

### Model Architecture
- **Input**: 8 features
  1. delta_z - Cluster z position relative to reference vertex
  2. delta_z_resunits - delta_z in resolution units
  3. cluster_z_sigma - Cluster z uncertainty
  4. cluster_d0 - Transverse impact parameter
  5. cluster_d0_sigma - d0 uncertainty
  6. cluster_qOverP - Charge over momentum
  7. cluster_qOverP_sigma - q/p uncertainty
  8. cluster_sumpt - Sum of track pT
- **Architecture**: 8 → 128 → 64 → 32 → 1
- **Activations**: ReLU (hidden), Sigmoid (output)
- **Parameters**: 11,521
- **Performance**: ~0.01-0.02 ms per inference

### Integration
- **Model loading**: Lazy static initialization in `src/clustering_functions.h:294`
- **Feature extraction**: `src/clustering_structs.h:328` (calcFeatures function)
- **Scoring**: Part of `updateScores()` called after clustering
- **Usage**: Score assigned to `cluster.scores[TESTML]` or `cluster.scores[Score::TESTML]`

### Files
- `share/models/model_weights.json` - Model weights loaded by C++
- `share/models/model.onnx` - Original ONNX model
- `share/models/normalization_params.json` - Feature normalization parameters
- `src/ml_model.h` - Standalone C++ neural network implementation

### Performance Optimization
Model is loaded **once** on first event (static initialization), not per-event. Expect ~100x speedup vs naive per-event loading.

### Documentation
- Quick start: `share/docs/ML_QUICK_START.md`
- Feature extraction: `share/docs/ML_FEATURE_EXTRACTION.md`
- Full integration: `share/docs/ML_MODEL_INTEGRATION.md`

### Retraining Workflow
1. Train model in Python with proper normalization
2. Export to ONNX: `tf2onnx.convert.from_keras()`
3. Save normalization params: `scaler.mean_` and `scaler.scale_`
4. Run: `python share/scripts/inspect_onnx_model.py --model share/models/model.onnx --export-weights`
5. Copy new files to `share/models/`
6. Update normalization in C++ if needed
7. Rebuild: `cd build && make`

## Development Notes

### Debug Mode
Set `#define debug true` at the top of executables to enable verbose logging.

### Multicore Processing
`ROOT::EnableImplicitMT()` is used to parallelize event processing across all CPU cores.

### Color Scheme
Standard colors defined in vertex_time.cxx:
- c1: #3f90da (blue)
- c2: #ffa90e (orange)
- c3: #bd1f01 (red)
- c4: #94a4a2 (gray)
- c5: #832db6 (purple)

## Git Status (Last Snapshot)

Current branch: main

Untracked files:
- .DS_Store
- hgtd_matching.cxx
- hgtd_matching_analysis.root

Recent commits:
- 78bff80: Up to date before claude implementation
- 261a1a8: add cmakelists.txt
- 5a61760: Update scripts to latest

## Common Workflows

### Running Full Analysis
```bash
cd build
./clustering_dt
# Output: Various plots in ../figs/ and event lists
```

### Analyzing Vertex Times
```bash
cd build
cmake .. && make
./clustering_dt
# Check figs/ for timing distribution plots
```

### Single File Analysis
Use generate_onefile_timeplot.cxx for quick checks on individual ROOT files.

## Error Analysis Tool

The `hgtd_matching.cxx` executable provides comprehensive diagnostics for algorithm failure modes:

### Usage
```bash
# Run the error analysis
root -l -b -q hgtd_matching.cxx

# Extract quantitative metrics and generate summary table
root -l -b -q extract_error_metrics.cxx
```

### Output
- **ROOT file**: `hgtd_matching_analysis.root` (all histograms)
- **PDF plots**: `error_analysis_plots/plot*.pdf` (11 diagnostic plots)
- **Summary table**: `error_metrics_summary.txt` (quantitative breakdown)

### Generated Plots

**Error Source 1: Wrong HGTD Matching**
- Plot 1: Track time residual distributions (all, HS, PU, matched, mismatched, cumulative)
- Plot 2: Matching efficiency vs kinematics (η, pT, nHGTDHits)
- Plot 3: Resolution vs matching quality (RMS vs nHits, resolution vs avg residual)

**Error Source 2: Wrong Cluster Selection**
- Plot 4: Cluster score distributions (HS vs PU for TRKPT and TRKPTZ)
- Plot 5: Selection efficiency vs event properties (n_clusters, n_jets, n_HS, PU fraction)
- Plot 6: Resolution impact (correct vs wrong selection, extra error)

**Error Source 3: Insufficient Tracks**
- Plot 7: Track multiplicity (HS valid distribution, success rate, 2D HS vs PU, track loss waterfall)
- Plot 8: Cluster sizes (HS vs PU, HS tracks in HS cluster)
- Plot 9: HS track fragmentation (time spread, z spread, success vs spread)
- Plot 10: Resolution vs track count (resolution scaling with N_HS)
- Plot 11: Additional diagnostics (selection outcomes, HS cluster rank, score ratios)

### Key Metrics Provided
- Track-level time residuals and matching quality
- Event-level matching efficiency by kinematics
- Cluster selection confusion matrix (TP, FP, FN, TN)
- Track loss through selection pipeline
- Success rate vs track multiplicity with efficiency thresholds
- HS track fragmentation in time and space

## Notes for Claude Code

- The main analysis entry point is `src/clustering_dt.cxx`
- Error analysis tool is `hgtd_matching.cxx` (comprehensive failure mode diagnostics)
- Core algorithms are in header files under `src/`
- When modifying analysis cuts, update constants in `src/clustering_constants.h`
- Always rebuild after header changes: `cd build && make`
- Plot styling uses ROOT's TColor and TCanvas
- Event loop uses TTreeReader for efficient branch access
- The project uses namespace `MyUtl` for utilities
