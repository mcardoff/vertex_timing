# ONNX Model Integration Guide

This guide shows how to integrate your `model.onnx` into the vertex timing analysis.

## Prerequisites

- ROOT 6.26+ (check with `root --version`)
- The model.onnx file in your project directory

## Option 1: Using ROOT's SOFIE (Recommended for ROOT 6.26+)

SOFIE (System for Optimized Fast Inference code Emit) compiles ONNX models to C++ code for fast inference.

### Step 1: Generate C++ Code from ONNX

Create a script to convert your ONNX model to C++ header:

```cpp
// generate_model_code.C
void generate_model_code() {
    using namespace TMVA::Experimental;

    // Parse the ONNX model
    SOFIE::RModelParser_ONNX parser;
    SOFIE::RModel model = parser.Parse("model.onnx");

    // Generate C++ code
    model.Generate();

    // Save to header file
    model.OutputGenerated("model_inference.hxx");

    std::cout << "Generated model_inference.hxx from model.onnx" << std::endl;
}
```

Run it:
```bash
root -l -b -q generate_model_code.C
```

### Step 2: Integrate into Your Analysis

In your analysis code (e.g., `clustering_dt.cxx`):

```cpp
// At the top of your file
#include "model_inference.hxx"

// One-time setup (before event loop)
TMVA_SOFIE_model::Session session("model.dat");

// In event loop, for each cluster
std::vector<float> features = {
    cluster.score_trkpt,
    cluster.score_trkptz,
    cluster.n_tracks,
    cluster.sumpt,
    // ... add all your features in the correct order
};

// Get prediction
std::vector<float> output = session.infer(features.data());
float ml_score = output[0];  // Probability or score

// Use score for selection
if (ml_score > threshold) {
    // Select this cluster
}
```

## Option 2: Using External ONNX Runtime (For older ROOT versions)

If ROOT < 6.26, use the external ONNX Runtime library.

### Installation:
```bash
# macOS
brew install onnxruntime

# or download from https://github.com/microsoft/onnxruntime/releases
```

### CMakeLists.txt addition:
```cmake
find_package(onnxruntime REQUIRED)

target_include_directories(clustering_dt PRIVATE ${ONNXRUNTIME_INCLUDE_DIRS})
target_link_libraries(clustering_dt PRIVATE ${ONNXRUNTIME_LIBRARIES})
```

### C++ Code:
```cpp
#include <onnxruntime/core/session/onnxruntime_cxx_api.h>

// Setup (once)
Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "vertex_timing");
Ort::SessionOptions session_options;
Ort::Session session(env, "model.onnx", session_options);

// Get input/output info
auto input_name = session.GetInputName(0, allocator);
auto input_shape = session.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();

// In event loop
std::vector<float> input_tensor_values = {feat1, feat2, feat3, ...};
auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
    memory_info,
    input_tensor_values.data(),
    input_tensor_values.size(),
    input_shape.data(),
    input_shape.size()
);

// Run inference
auto output_tensors = session.Run(
    Ort::RunOptions{nullptr},
    &input_name, &input_tensor, 1,
    &output_name, 1
);

float* output = output_tensors.front().GetTensorMutableData<float>();
float ml_score = output[0];
```

## Critical: Feature Order

**The features MUST be in the exact same order as during training!**

Create a helper function to ensure consistency:

```cpp
std::vector<float> extract_features(const Cluster& cluster, const Branch& branch) {
    std::vector<float> features;

    // Add features in EXACT training order
    features.push_back(cluster.score_trkpt);      // Feature 0
    features.push_back(cluster.score_trkptz);     // Feature 1
    features.push_back(cluster.n_tracks);          // Feature 2
    features.push_back(cluster.sumpt);             // Feature 3
    features.push_back(cluster.avg_eta);           // Feature 4
    // ... continue for all features

    return features;
}
```

## Testing Your Integration

Create a simple test script:

```cpp
// test_onnx_model.C
void test_onnx_model() {
    #include "model_inference.hxx"

    TMVA_SOFIE_model::Session session("model.dat");

    // Test with dummy input
    std::vector<float> test_input = {1.0, 2.0, 3.0, 4.0, 5.0};  // Adjust size
    std::vector<float> output = session.infer(test_input.data());

    std::cout << "Test prediction: " << output[0] << std::endl;

    // Test with different inputs
    test_input = {-1.0, 0.0, 1.0, 2.0, 3.0};
    output = session.infer(test_input.data());
    std::cout << "Test prediction 2: " << output[0] << std::endl;
}
```

## Feature Extraction Example

Based on typical cluster scoring features:

```cpp
std::vector<float> extract_cluster_features(
    const Cluster& cluster,
    const Branch& branch,
    int n_jets
) {
    std::vector<float> features;

    // Cluster properties
    features.push_back(cluster.trackIndices.size());  // n_tracks
    features.push_back(cluster.score_trkpt);           // score_trkpt
    features.push_back(cluster.score_trkptz);          // score_trkptz

    // Calculate additional features
    double sumpt = 0, sumpt2 = 0, avg_eta = 0;
    int n_hs = 0, n_pu = 0;

    for (auto trk_idx : cluster.trackIndices) {
        double pt = branch.trackPt[trk_idx];
        double eta = branch.trackEta[trk_idx];

        sumpt += pt;
        sumpt2 += pt * pt;
        avg_eta += eta;

        if (branch.trackToTruthvtx[trk_idx] == 0) n_hs++;
        else n_pu++;
    }

    avg_eta /= cluster.trackIndices.size();

    features.push_back(sumpt);
    features.push_back(sumpt2);
    features.push_back(avg_eta);
    features.push_back(n_hs);
    features.push_back(n_pu);
    features.push_back(n_jets);

    // Add cluster time if available
    if (cluster.values.size() > 0) {
        features.push_back(cluster.values[0]);
    }

    return features;
}
```

## Integration into Main Analysis

In `clustering_dt.cxx` or your main analysis:

```cpp
// After clustering, before selection
for (size_t ic = 0; ic < clusters.size(); ic++) {
    auto& cluster = clusters[ic];

    // Extract features
    std::vector<float> features = extract_cluster_features(cluster, branch, n_jets);

    // Get ML score
    std::vector<float> output = session.infer(features.data());
    float ml_score = output[0];

    // Store score (add new field to Cluster struct if needed)
    cluster.ml_score = ml_score;
}

// Select cluster with highest ML score
auto best_cluster = std::max_element(
    clusters.begin(),
    clusters.end(),
    [](const Cluster& a, const Cluster& b) {
        return a.ml_score < b.ml_score;
    }
);
```

## Next Steps

1. Check your ROOT version: `root --version`
2. If ROOT ≥ 6.26: Use Option 1 (SOFIE)
3. If ROOT < 6.26: Use Option 2 (ONNX Runtime)
4. Document your feature order in a comment or separate file
5. Create test script to verify inference works
6. Integrate into main analysis
7. Compare ML-based selection vs current scoring methods

## Troubleshooting

**"SOFIE not found"**: ROOT too old, use Option 2
**"Wrong input shape"**: Check feature vector size matches model input
**"Segmentation fault"**: Verify feature order and types (float vs double)
**"Different predictions than Python"**: Check feature normalization/scaling

## Performance Tips

- Create session once, reuse for all events
- Pre-allocate feature vectors to avoid repeated allocation
- Consider batching if processing many clusters per event
- Profile to ensure ML inference doesn't dominate runtime
