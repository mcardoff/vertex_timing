#ifndef ML_MODEL_H
#define ML_MODEL_H

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <fstream>
#include <stdexcept>
#include "json.hpp"

// Standalone Neural Network implementation
// Architecture: 8 -> 128 -> 64 -> 32 -> 1
// Activation: ReLU for hidden layers, Sigmoid for output
//
// Performance notes:
//   - Weight matrices are stored as flat 1D arrays (row-major) for cache locality
//   - Intermediate activations use fixed-size std::array to avoid heap allocation
//   - The predict() hot path does zero heap allocation per call
class MLModel {
private:
    // Layer dimensions (compile-time constants match the trained model)
    static constexpr int N_IN  = 8;
    static constexpr int N_H1  = 128;
    static constexpr int N_H2  = 64;
    static constexpr int N_H3  = 32;
    static constexpr int N_OUT = 1;

    // Weights stored as flat row-major arrays: W[i][j] == w1[i*N_H1 + j]
    // This keeps each row contiguous in memory so the inner dot-product loop
    // streams sequentially through cache lines.
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

    // Dense layer with ReLU — writes into a pre-allocated output array.
    // No heap allocation; caller owns the output buffer.
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
    void load_weights(const std::string& json_path) {
        std::ifstream file(json_path);
        if (!file.is_open())
            throw std::runtime_error("Cannot open weights file: " + json_path);

        nlohmann::json j;
        file >> j;

        // Key names — try current naming scheme first, fall back to legacy
        auto pick = [&](const std::string& primary, const std::string& fallback) -> std::string {
            return j.contains(primary) ? primary : fallback;
        };

        std::string w1k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_24_1/Cast/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_1/Cast/ReadVariableOp:0");
        std::string b1k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_24_1/BiasAdd/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_1/BiasAdd/ReadVariableOp:0");
        std::string w2k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_25_1/Cast/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_1_2/Cast/ReadVariableOp:0");
        std::string b2k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_25_1/BiasAdd/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_1_2/BiasAdd/ReadVariableOp:0");
        std::string w3k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_26_1/Cast/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_2_1/Cast/ReadVariableOp:0");
        std::string b3k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_26_1/BiasAdd/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_2_1/BiasAdd/ReadVariableOp:0");
        std::string w4k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_27_1/Cast/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_3_1/Cast/ReadVariableOp:0");
        std::string b4k = pick(
            "StatefulPartitionedCall/sequential_6_1/dense_27_1/Add/ReadVariableOp:0",
            "StatefulPartitionedCall/sequential_1/dense_3_1/Add/ReadVariableOp:0");

        // Load flat weight vectors directly — no reshape needed
        w1 = j[w1k]["values"].get<std::vector<float>>();
        b1 = j[b1k]["values"].get<std::vector<float>>();
        w2 = j[w2k]["values"].get<std::vector<float>>();
        b2 = j[b2k]["values"].get<std::vector<float>>();
        w3 = j[w3k]["values"].get<std::vector<float>>();
        b3 = j[b3k]["values"].get<std::vector<float>>();
        w4 = j[w4k]["values"].get<std::vector<float>>();
        b4 = j[b4k]["values"].get<std::vector<float>>();

        // Detect input size in case model was retrained with different feature count
        n_input_features = static_cast<int>(w1.size()) / N_H1;

        weights_loaded = true;
    }

    // Forward pass: zero heap allocation in the hot path.
    // Intermediate activations live on the stack as fixed-size arrays.
    float predict(const std::vector<float>& features) const {
        if (!weights_loaded)
            throw std::runtime_error("Model weights not loaded. Call load_weights() first.");

        // Stack-allocated activation buffers — no heap involvement
        std::array<float, N_H1>  h1;
        std::array<float, N_H2>  h2;
        std::array<float, N_H3>  h3;
        std::array<float, N_OUT> out;

        dense_relu   <N_IN, N_H1>(features.data(), w1.data(), b1.data(), h1.data() );
        dense_relu   <N_H1, N_H2>(h1.data(),       w2.data(), b2.data(), h2.data() );
        dense_relu   <N_H2, N_H3>(h2.data(),       w3.data(), b3.data(), h3.data() );
        dense_sigmoid<N_H3, N_OUT>(h3.data(),      w4.data(), b4.data(), out.data());

        return out[0];
    }

    // Batch prediction (convenience wrapper)
    std::vector<float> predict_batch(const std::vector<std::vector<float>>& batch) const {
        std::vector<float> predictions;
        predictions.reserve(batch.size());
        for (const auto& f : batch)
            predictions.push_back(predict(f));
        return predictions;
    }
};

#endif // ML_MODEL_H
