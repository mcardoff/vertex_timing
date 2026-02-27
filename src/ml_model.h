#ifndef ML_MODEL_H
#define ML_MODEL_H

// ---------------------------------------------------------------------------
// ml_model.h
//   Standalone C++ implementation of the trained DNN used for cluster
//   selection.  No external inference framework is required; weights are
//   loaded from a JSON file exported by temp-conversion-dir/convert_model.py.
//
//   Architecture: 11 → 128 → 64 → 32 → 1
//   Activations:  ReLU on hidden layers, Sigmoid on the output layer.
//   Parameters:   11,841
//
//   Performance design:
//     • Weight matrices are stored as flat row-major 1D vectors for cache
//       locality — W[i][j] is at w[i * N_OUT_L + j].
//     • Intermediate activations are fixed-size std::array on the stack;
//       the predict() hot path performs zero heap allocation per call.
//     • The model is loaded once via static initialisation in
//       clusterTracksInTime (clustering_functions.h); all subsequent calls
//       use the already-loaded weights at ~0.01–0.02 ms per inference.
//
//   To update the model after retraining:
//     1. cd temp-conversion-dir && ./convert_model.py
//     2. Rebuild: cd ../vertex_timing/build && make
//
//   Public interface:
//     load_weights(path) — parse JSON and populate weight vectors
//     predict(features)  — single-sample forward pass → [0, 1]
//     predict_batch(...)  — convenience wrapper over predict()
// ---------------------------------------------------------------------------

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <fstream>
#include <stdexcept>
#include "json.hpp"

// ---------------------------------------------------------------------------
// MLModel
//   Encapsulates the weight storage, forward-pass computation, and weight
//   loading for the cluster-selection DNN.
// ---------------------------------------------------------------------------
class MLModel {
private:
    // -----------------------------------------------------------------------
    // Layer dimensions
    //   Compile-time constants; must match the exported model exactly.
    //   Changing these requires retraining and re-exporting weights.
    // -----------------------------------------------------------------------
    static constexpr int N_IN  = 11;
    static constexpr int N_H1  = 128;
    static constexpr int N_H2  = 64;
    static constexpr int N_H3  = 32;
    static constexpr int N_OUT = 1;

    // -----------------------------------------------------------------------
    // Weight storage
    //   Flat row-major vectors: W[i][j] == w[i * N_OUT_L + j].  Keeping
    //   each row contiguous means the inner dot-product loop streams
    //   sequentially through cache lines.
    // -----------------------------------------------------------------------
    std::vector<float> w1;  // N_IN  * N_H1
    std::vector<float> b1;  // N_H1
    std::vector<float> w2;  // N_H1 * N_H2
    std::vector<float> b2;  // N_H2
    std::vector<float> w3;  // N_H2 * N_H3
    std::vector<float> b3;  // N_H3
    std::vector<float> w4;  // N_H3 * N_OUT
    std::vector<float> b4;  // N_OUT

    bool weights_loaded = false;
    int n_input_features = N_IN;

    inline float relu   (float x) const { return x > 0.0f ? x : 0.0f; }
    inline float sigmoid(float x) const { return 1.0f / (1.0f + std::exp(-x)); }

    // -----------------------------------------------------------------------
    // dense_relu / dense_sigmoid
    //   Templated dense-layer helpers.  Each computes:
    //     output[j] = activation( bias[j] + Σ_i input[i] * W[i*N_OUT_L+j] )
    //   Output is written into a caller-owned buffer; no heap allocation.
    //   Template parameters fix the loop bounds at compile time so the
    //   compiler can unroll and vectorise aggressively.
    // -----------------------------------------------------------------------
    template<int N_IN_L, int N_OUT_L>
    void dense_relu(
        const float* __restrict__ input,
        const float* __restrict__ W,   // flat row-major [N_IN_L * N_OUT_L]
        const float* __restrict__ b,
        float* __restrict__ output
    ) const {
        for (int j = 0; j < N_OUT_L; ++j) {
            float sum = b[j];
            for (int i = 0; i < N_IN_L; ++i)
                sum += input[i] * W[i * N_OUT_L + j];
            output[j] = relu(sum);
        }
    }

    // Dense layer with Sigmoid — same pattern, different activation.
    template<int N_IN_L, int N_OUT_L>
    void dense_sigmoid(
        const float* __restrict__ input,
        const float* __restrict__ W,
        const float* __restrict__ b,
        float* __restrict__ output
    ) const {
        for (int j = 0; j < N_OUT_L; ++j) {
            float sum = b[j];
            for (int i = 0; i < N_IN_L; ++i)
                sum += input[i] * W[i * N_OUT_L + j];
            output[j] = sigmoid(sum);
        }
    }

public:
    // -----------------------------------------------------------------------
    // load_weights
    //   Parses the JSON weight file at json_path and populates the six
    //   weight/bias vectors.  Keys are stable names written by
    //   convert_model.py (w1/b1/w2/b2/w3/b3); layer order is determined
    //   by tensor shape in the script, not by TF node names.
    // -----------------------------------------------------------------------
    void load_weights(const std::string& json_path) {
        std::ifstream file(json_path);
        if (!file.is_open())
            throw std::runtime_error("Cannot open weights file: " + json_path);

        nlohmann::json j;
        file >> j;

        w1 = j["w1"]["values"].get<std::vector<float>>();
        b1 = j["b1"]["values"].get<std::vector<float>>();
        w2 = j["w2"]["values"].get<std::vector<float>>();
        b2 = j["b2"]["values"].get<std::vector<float>>();
        w3 = j["w3"]["values"].get<std::vector<float>>();
        b3 = j["b3"]["values"].get<std::vector<float>>();
        w4 = j["w4"]["values"].get<std::vector<float>>();
        b4 = j["b4"]["values"].get<std::vector<float>>();

        // Detect input size in case model was retrained with different feature count
        n_input_features = static_cast<int>(w1.size()) / N_H1;

        weights_loaded = true;
    }

    // -----------------------------------------------------------------------
    // predict
    //   Single-sample forward pass.  All intermediate activation buffers
    //   are stack-allocated fixed-size arrays so no heap allocation occurs
    //   in the hot path.  Returns the sigmoid output in [0, 1]; values > 0.5
    //   indicate a predicted hard-scatter cluster.
    //   Throws if weights have not been loaded via load_weights().
    // -----------------------------------------------------------------------
    float predict(const std::vector<float>& features) const {
        if (!weights_loaded)
            throw std::runtime_error("Model weights not loaded. Call load_weights() first.");

        // Stack-allocated activation buffers — no heap involvement
        std::array<float, N_H1>  h1;
        std::array<float, N_H2>  h2;
        std::array<float, N_H3>  h3;
        std::array<float, N_OUT> out;

        dense_relu   <N_IN, N_H1> (features.data(), w1.data(), b1.data(), h1.data());
        dense_relu   <N_H1, N_H2> (h1.data(),       w2.data(), b2.data(), h2.data());
        dense_relu   <N_H2, N_H3> (h2.data(),       w3.data(), b3.data(), h3.data());
        dense_sigmoid<N_H3, N_OUT>(h3.data(),        w4.data(), b4.data(), out.data());

        return out[0];
    }

    // -----------------------------------------------------------------------
    // predict_batch
    //   Convenience wrapper that calls predict() for each sample in batch
    //   and returns the results as a vector.
    // -----------------------------------------------------------------------
    std::vector<float> predict_batch(const std::vector<std::vector<float>>& batch) const {
        std::vector<float> predictions;
        predictions.reserve(batch.size());
        for (const auto& f : batch)
            predictions.push_back(predict(f));
        return predictions;
    }
};

#endif // ML_MODEL_H
