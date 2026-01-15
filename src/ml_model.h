#ifndef ML_MODEL_H
#define ML_MODEL_H

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <stdexcept>
#include "json.hpp"  // Will need nlohmann/json or simple parser

// Standalone Neural Network implementation
// Architecture: 9 -> 128 -> 64 -> 32 -> 1
// Activation: ReLU for hidden layers, Sigmoid for output
class MLModel {
private:
    // Layer weights and biases
    std::vector<std::vector<float>> W1;  // 8 x 128 (or 9 x 128 for old model)
    std::vector<float> b1;                // 128
    std::vector<std::vector<float>> W2;  // 128 x 64
    std::vector<float> b2;                // 64
    std::vector<std::vector<float>> W3;  // 64 x 32
    std::vector<float> b3;                // 32
    std::vector<std::vector<float>> W4;  // 32 x 1
    std::vector<float> b4;                // 1

    bool weights_loaded = false;
    int n_input_features = 8;  // Expected number of input features

    // ReLU activation
    inline float relu(float x) const {
        return x > 0.0f ? x : 0.0f;
    }

    // Sigmoid activation
    inline float sigmoid(float x) const {
        return 1.0f / (1.0f + std::exp(-x));
    }

    // Dense layer: output = activation(input * W + b)
    std::vector<float> dense_layer(
        const std::vector<float>& input,
        const std::vector<std::vector<float>>& W,
        const std::vector<float>& b,
        bool use_relu = true
    ) const {
        size_t n_out = b.size();
        std::vector<float> output(n_out, 0.0f);

        // Matrix multiplication: output = input * W + b
        for (size_t i = 0; i < n_out; i++) {
            float sum = b[i];
            for (size_t j = 0; j < input.size(); j++) {
                sum += input[j] * W[j][i];
            }
            output[i] = use_relu ? relu(sum) : sigmoid(sum);
        }

        return output;
    }

public:
    // Load weights from JSON file
    void load_weights(const std::string& json_path) {
        std::ifstream file(json_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open weights file: " + json_path);
        }

        nlohmann::json weights_json;
        file >> weights_json;

        // Extract weights - NOTE: The layer naming from TensorFlow/ONNX
        // Try new naming scheme first (sequential_1_1/dense_4_1), fallback to old (sequential_1/dense_1)
        std::string w1_key = "StatefulPartitionedCall/sequential_6_1/dense_24_1/Cast/ReadVariableOp:0";
	std::string b1_key = "StatefulPartitionedCall/sequential_6_1/dense_24_1/BiasAdd/ReadVariableOp:0";

        if (!weights_json.contains(w1_key)) {
            // Fallback to old naming
            w1_key = "StatefulPartitionedCall/sequential_1/dense_1/Cast/ReadVariableOp:0";
            b1_key = "StatefulPartitionedCall/sequential_1/dense_1/BiasAdd/ReadVariableOp:0";
        }

        // Layer 1 (8 -> 128)
        auto w1_flat = weights_json[w1_key]["values"].get<std::vector<float>>();
        b1 = weights_json[b1_key]["values"].get<std::vector<float>>();

        // Detect input size from weight matrix
        int w1_size = w1_flat.size();
        int detected_inputs = w1_size / 128;
        n_input_features = detected_inputs;

        // Reshape W1: [n_input_features, 128]
        W1.resize(n_input_features, std::vector<float>(128));
        for (int i = 0; i < n_input_features; i++) {
            for (int j = 0; j < 128; j++) {
                W1[i][j] = w1_flat[i * 128 + j];
            }
        }

        // Layer 2 (128 -> 64)
	std::string w2_key = "StatefulPartitionedCall/sequential_6_1/dense_25_1/Cast/ReadVariableOp:0";
	std::string b2_key = "StatefulPartitionedCall/sequential_6_1/dense_25_1/BiasAdd/ReadVariableOp:0";
        if (!weights_json.contains(w2_key)) {
            w2_key = "StatefulPartitionedCall/sequential_1/dense_1_2/Cast/ReadVariableOp:0";
            b2_key = "StatefulPartitionedCall/sequential_1/dense_1_2/BiasAdd/ReadVariableOp:0";
        }

        auto w2_flat = weights_json[w2_key]["values"].get<std::vector<float>>();
        b2 = weights_json[b2_key]["values"].get<std::vector<float>>();

        // Reshape W2: [128, 64]
        W2.resize(128, std::vector<float>(64));
        for (int i = 0; i < 128; i++) {
            for (int j = 0; j < 64; j++) {
                W2[i][j] = w2_flat[i * 64 + j];
            }
        }

        // Layer 3 (64 -> 32)
        std::string w3_key = "StatefulPartitionedCall/sequential_6_1/dense_26_1/Cast/ReadVariableOp:0";
	std::string b3_key = "StatefulPartitionedCall/sequential_6_1/dense_26_1/BiasAdd/ReadVariableOp:0";
        if (!weights_json.contains(w3_key)) {
            w3_key = "StatefulPartitionedCall/sequential_1/dense_2_1/Cast/ReadVariableOp:0";
            b3_key = "StatefulPartitionedCall/sequential_1/dense_2_1/BiasAdd/ReadVariableOp:0";
        }

        auto w3_flat = weights_json[w3_key]["values"].get<std::vector<float>>();
        b3 = weights_json[b3_key]["values"].get<std::vector<float>>();

        // Reshape W3: [64, 32]
        W3.resize(64, std::vector<float>(32));
        for (int i = 0; i < 64; i++) {
            for (int j = 0; j < 32; j++) {
                W3[i][j] = w3_flat[i * 32 + j];
            }
        }

        // Layer 4 (32 -> 1)
        std::string w4_key = "StatefulPartitionedCall/sequential_6_1/dense_27_1/Cast/ReadVariableOp:0";
	std::string b4_key = "StatefulPartitionedCall/sequential_6_1/dense_27_1/Add/ReadVariableOp:0";
        if (!weights_json.contains(w4_key)) {
            w4_key = "StatefulPartitionedCall/sequential_1/dense_3_1/Cast/ReadVariableOp:0";
            b4_key = "StatefulPartitionedCall/sequential_1/dense_3_1/Add/ReadVariableOp:0";
        }

        auto w4_flat = weights_json[w4_key]["values"].get<std::vector<float>>();
        b4 = weights_json[b4_key]["values"].get<std::vector<float>>();

        // Reshape W4: [32, 1]
        W4.resize(32, std::vector<float>(1));
        for (int i = 0; i < 32; i++) {
            W4[i][0] = w4_flat[i];
        }

        weights_loaded = true;
    }

    // Forward pass: predict(features) -> score
    float predict(const std::vector<float>& features) const {
        if (!weights_loaded) {
            throw std::runtime_error("Model weights not loaded. Call load_weights() first.");
        }

        if (features.size() != static_cast<size_t>(n_input_features)) {
            throw std::runtime_error("Expected " + std::to_string(n_input_features) +
                                   " features, got " + std::to_string(features.size()));
        }

        // Forward pass through network
        auto h1 = dense_layer(features, W1, b1, true);   // 9 -> 128, ReLU
        auto h2 = dense_layer(h1, W2, b2, true);         // 128 -> 64, ReLU
        auto h3 = dense_layer(h2, W3, b3, true);         // 64 -> 32, ReLU
        auto output = dense_layer(h3, W4, b4, false);    // 32 -> 1, Sigmoid

        return output[0];
    }

    // Batch prediction
    std::vector<float> predict_batch(const std::vector<std::vector<float>>& batch) const {
        std::vector<float> predictions;
        predictions.reserve(batch.size());

        for (const auto& features : batch) {
            predictions.push_back(predict(features));
        }

        return predictions;
    }
};

#endif // ML_MODEL_H
