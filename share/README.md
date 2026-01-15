# Shared Resources Directory

This directory contains common resources, documentation, and utilities for the HGTD vertex timing analysis project.

## Directory Structure

```
share/
├── models/                  # ML model files
│   ├── model.onnx          # ONNX model (8 features)
│   ├── model_weights.json  # Exported weights for C++ (330 KB)
│   └── normalization_params.json  # Feature normalization parameters
│
├── scripts/                 # Utility scripts
│   ├── event_display.py           # Event visualization
│   ├── inspect_onnx_model.py      # Model inspection tool
│   ├── export_onnx_model.py       # ONNX export helper
│   └── model_evaluation_helper.py # Model integration guide
│
└── docs/                    # Documentation
    ├── SUMMARY.md                  # Project status overview
    ├── ML_QUICK_START.md           # Quick ML integration guide
    ├── ML_FEATURE_EXTRACTION.md    # Feature calculation guide
    ├── ML_MODEL_INTEGRATION.md     # Complete ML integration
    ├── ONNX_INTEGRATION_GUIDE.md   # Alternative ONNX approaches
    └── ERROR_ANALYSIS_PLAN.md      # Error analysis methodology

## ML Model Files

### model_weights.json
- **Size**: 330 KB
- **Architecture**: 8 → 128 → 64 → 32 → 1
- **Parameters**: 11,521
- **Usage**: Loaded by C++ code at `src/clustering_functions.h:296`
- **Path in code**: `../share/models/model_weights.json` (from build directory)

### normalization_params.json
- **Content**: Feature means and standard deviations
- **Format**:
  ```json
  {
    "means": [mean_feat0, mean_feat1, ...],
    "stds": [std_feat0, std_feat1, ...],
    "feature_names": ["delta_z", "delta_z_resunits", ...]
  }
  ```
- **Usage**: Apply normalization in C++ before model inference

### model.onnx
- **Original model**: Complete ONNX format
- **Input**: 8 features (batch_size, 8)
- **Output**: Sigmoid probability (batch_size, 1)

## Scripts

### event_display.py
Event visualization for error analysis.

**Usage**:
```bash
python share/scripts/event_display.py --file_num X --event_num Y --extra_time 0
```

### inspect_onnx_model.py
Inspect ONNX model structure and export weights.

**Usage**:
```bash
python share/scripts/inspect_onnx_model.py --model share/models/model.onnx
python share/scripts/inspect_onnx_model.py --model share/models/model.onnx --export-weights
python share/scripts/inspect_onnx_model.py --model share/models/model.onnx --test
```

## Documentation

### Quick Start
Start with `docs/SUMMARY.md` for project overview, then `docs/ML_QUICK_START.md` for ML integration.

### Feature Extraction
See `docs/ML_FEATURE_EXTRACTION.md` for detailed guide on calculating the 8 ML features.

### Error Analysis
See `docs/ERROR_ANALYSIS_PLAN.md` for methodology and `build/error_metrics_summary.txt` for results.

## Integration Notes

### C++ Code References
The ML model is loaded from:
- `src/clustering_functions.h`: Line 296
- `src/test_ml_model.cxx`: Line 19

### Python Training
When retraining the model, update these files:
1. Export new weights: `python share/scripts/inspect_onnx_model.py --export-weights`
2. Copy new `model.onnx` to `share/models/`
3. Copy new `model_weights.json` to `share/models/`
4. Copy new `normalization_params.json` to `share/models/`
5. Update normalization values in `src/clustering_structs.h` (if needed)

## File Sizes
```
models/model_weights.json       330 KB
models/model.onnx                49 KB
models/normalization_params.json  1 KB
scripts/                         ~30 KB total
docs/                            ~30 KB total
```

## Access from Code

All paths are relative to the `build/` directory where executables run:

```cpp
// From build/ directory
ml_model.load_weights("../share/models/model_weights.json");
```

```bash
# From project root
python share/scripts/event_display.py --file_num 0 --event_num 123
```
