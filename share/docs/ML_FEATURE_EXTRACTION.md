# ML Model Feature Extraction Guide

## Model Status: ✓ Ready

The ML model has been successfully integrated and tested with the updated 8-feature model.

## Feature Order (CRITICAL)

The model expects **exactly 8 features in this order**:

1. **delta_z** - Difference between cluster z position and reference vertex z
2. **delta_z_resunits** - delta_z normalized by uncertainty (resolution units)
3. **cluster_z_sigma** - Uncertainty in cluster z position
4. **cluster_d0** - Transverse impact parameter of cluster
5. **cluster_d0_sigma** - Uncertainty in d0
6. **cluster_qOverP** - Charge over momentum (q/p) for cluster
7. **cluster_qOverP_sigma** - Uncertainty in q/p
8. **cluster_sumpt** - Sum of transverse momentum of tracks in cluster

## C++ Implementation Example

### Step 1: Include Header

```cpp
#include "ml_model.h"
```

### Step 2: Load Model (Before Event Loop)

```cpp
// One-time initialization
MLModel ml_model;
try {
    ml_model.load_weights("model_weights.json");
    std::cout << "✓ ML model loaded successfully" << std::endl;
} catch (const std::exception& e) {
    std::cerr << "ERROR loading ML model: " << e.what() << std::endl;
    return 1;
}
```

### Step 3: Extract Features for Each Cluster

Here's how to calculate each feature from your cluster data:

```cpp
std::vector<float> extract_ml_features(
    const Cluster& cluster,
    const Branch& branch,
    double reference_vtx_z
) {
    std::vector<float> features(8);

    // Calculate cluster properties
    double cluster_z = 0.0;
    double cluster_z_sigma = 0.0;
    double cluster_d0 = 0.0;
    double cluster_d0_sigma = 0.0;
    double cluster_qOverP = 0.0;
    double cluster_qOverP_sigma = 0.0;
    double cluster_sumpt = 0.0;

    int n_tracks = cluster.trackIndices.size();

    // Aggregate track properties
    for (auto trk_idx : cluster.trackIndices) {
        // z position (weighted by pT)
        double pt = branch.trackPt[trk_idx];
        double z0 = branch.trackZ0[trk_idx];
        cluster_z += z0 * pt;
        cluster_sumpt += pt;

        // d0
        double d0 = branch.trackD0[trk_idx];  // If you have this branch
        cluster_d0 += d0 * pt;

        // q/p
        double qOverP = branch.trackQOverP[trk_idx];  // If you have this branch
        cluster_qOverP += qOverP * pt;

        // Uncertainties (add in quadrature)
        double z0_var = branch.trackVar_z0[trk_idx];
        cluster_z_sigma += z0_var;

        // Add d0_sigma and qOverP_sigma if you have those branches
    }

    // Normalize by sumpt
    if (cluster_sumpt > 0) {
        cluster_z /= cluster_sumpt;
        cluster_d0 /= cluster_sumpt;
        cluster_qOverP /= cluster_sumpt;
    }

    // Take sqrt for uncertainties
    cluster_z_sigma = std::sqrt(cluster_z_sigma / n_tracks);

    // Calculate delta_z
    double delta_z = cluster_z - reference_vtx_z;
    double delta_z_resunits = delta_z / cluster_z_sigma;

    // Fill feature vector IN ORDER
    features[0] = delta_z;              // Feature 0
    features[1] = delta_z_resunits;     // Feature 1
    features[2] = cluster_z_sigma;      // Feature 2
    features[3] = cluster_d0;           // Feature 3
    features[4] = cluster_d0_sigma;     // Feature 4
    features[5] = cluster_qOverP;       // Feature 5
    features[6] = cluster_qOverP_sigma; // Feature 6
    features[7] = cluster_sumpt;        // Feature 7

    return features;
}
```

### Step 4: Score Clusters in Event Loop

```cpp
// In your event loop, after clustering
for (auto& cluster : clusters) {
    // Extract features
    std::vector<float> features = extract_ml_features(
        cluster,
        branch,
        reference_vtx_z  // Could be RecoVtx_z[0] or truth vertex
    );

    // Get ML prediction (sigmoid output: 0-1, higher = more likely HS)
    float ml_score = ml_model.predict(features);

    // Store in cluster
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

## Feature Normalization

**Important**: Check if your Python training script applied any feature normalization (StandardScaler, MinMaxScaler, etc.). If so, you MUST apply the same normalization in C++ using the same mean/std values from training.

Example with StandardScaler:

```cpp
// If you used StandardScaler in Python, you need these values
const std::vector<float> feature_means = {
    /* means from your scaler.mean_ */
};
const std::vector<float> feature_stds = {
    /* stds from your scaler.scale_ */
};

// Normalize features before prediction
for (size_t i = 0; i < features.size(); i++) {
    features[i] = (features[i] - feature_means[i]) / feature_stds[i];
}

float ml_score = ml_model.predict(features);
```

## Branch Requirements

To calculate all features, you need these branches in your ROOT file:

**Required**:
- `Track_z0` or `Track_z0_wrtPV` ✓
- `Track_pt` ✓
- `Track_var_z0` ✓

**May Need to Add**:
- `Track_d0` or `Track_d0_wrtPV`
- `Track_var_d0`
- `Track_qOverP`
- `Track_var_qOverP`

If these branches don't exist, you may need to:
1. Compute them from other track parameters
2. Use proxy variables (e.g., approximate d0_sigma from z0_sigma)
3. Set them to reasonable default values

## Testing Your Implementation

```cpp
// Test with known values
std::vector<float> test_features = {
    0.5f,   // delta_z
    1.0f,   // delta_z_resunits
    0.2f,   // cluster_z_sigma
    0.1f,   // cluster_d0
    0.05f,  // cluster_d0_sigma
    0.001f, // cluster_qOverP
    0.0001f,// cluster_qOverP_sigma
    50.0f   // cluster_sumpt
};

float score = ml_model.predict(test_features);
std::cout << "Test score: " << score << std::endl;
```

## Performance Expectations

Based on your training results:
- Model has sigmoid output: 0 = pileup, 1 = hard scatter
- Training accuracy: ~60.8% of clusters are hard scatter in your dataset
- The model should significantly outperform simple TRKPT/TRKPTZ scoring

Threshold tuning:
- Start with threshold = 0.5 (standard for binary classification)
- Tune based on efficiency vs purity trade-off
- May want lower threshold (~0.3-0.4) to maximize efficiency

## Debugging

If predictions seem wrong:

1. **Check feature order** - Most common error!
2. **Print feature values** - Are they in reasonable ranges?
3. **Check for NaN/Inf** - Handle missing values properly
4. **Verify normalization** - Must match training exactly
5. **Compare with Python** - Run same features through Python model

## Integration Checklist

- [ ] Added `#include "ml_model.h"` to analysis file
- [ ] Load model before event loop with error handling
- [ ] Implement feature extraction function
- [ ] Verify all 8 features are calculated correctly
- [ ] Check if normalization is needed (compare with Python training)
- [ ] Add ml_score field to Cluster struct
- [ ] Score all clusters in event loop
- [ ] Select best cluster by ml_score
- [ ] Compare ML results with TRKPT/TRKPTZ methods
- [ ] Plot ML score distributions (HS vs PU)
- [ ] Tune threshold for optimal efficiency
