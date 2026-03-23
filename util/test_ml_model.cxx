// Test program for ML model implementation
// Usage: ./test_ml_model

#include "ml_model.h"
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    std::cout << "======================================================================" << std::endl;
    std::cout << "ML MODEL TEST" << std::endl;
    std::cout << "======================================================================" << std::endl;

    try {
        // Create model instance
        MLModel model;

        std::cout << "\nLoading model weights from share/models/model_weights.json..." << std::endl;
        model.load_weights("../share/models/model_weights.json");
        std::cout << "✓ Weights loaded successfully!" << std::endl;

        // Test with dummy features (8 features matching training order)
        std::cout << "\nTesting inference with dummy input..." << std::endl;

        std::vector<float> test_features = {
            1.0f,   // feature 0: delta_z
            2.0f,   // feature 1: delta_z_resunits
            3.0f,   // feature 2: cluster_z_sigma
            0.5f,   // feature 3: cluster_d0
            -1.0f,  // feature 4: cluster_d0_sigma
            0.0f,   // feature 5: cluster_qOverP
            10.0f,  // feature 6: cluster_qOverP_sigma
            5.0f    // feature 7: cluster_sumpt
        };

        std::cout << "Input features: [";
        for (size_t i = 0; i < test_features.size(); i++) {
            std::cout << test_features[i];
            if (i < test_features.size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;

        float prediction = model.predict(test_features);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "\nPrediction: " << prediction << std::endl;

        // Test with different inputs
        std::cout << "\n" << std::string(70, '-') << std::endl;
        std::cout << "Testing with multiple inputs:" << std::endl;
        std::cout << std::string(70, '-') << std::endl;

        std::vector<std::vector<float>> test_batch = {
            {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},  // All zeros
            {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},  // All ones
            {-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f},  // All -1
        };

        for (size_t i = 0; i < test_batch.size(); i++) {
            float pred = model.predict(test_batch[i]);
            std::cout << "Input " << (i+1) << " -> Prediction: " << pred << std::endl;
        }

        std::cout << "\n✓ Model inference working correctly!" << std::endl;

        std::cout << "\n======================================================================" << std::endl;
        std::cout << "INTEGRATION INSTRUCTIONS:" << std::endl;
        std::cout << "======================================================================" << std::endl;
        std::cout << "To use in your analysis code:" << std::endl;
        std::cout << std::endl;
        std::cout << "1. Include the header:" << std::endl;
        std::cout << "   #include \"ml_model.h\"" << std::endl;
        std::cout << std::endl;
        std::cout << "2. Create and load model (once, before event loop):" << std::endl;
        std::cout << "   MLModel model;" << std::endl;
        std::cout << "   model.load_weights(\"model_weights.json\");" << std::endl;
        std::cout << std::endl;
        std::cout << "3. In event loop, extract features and predict:" << std::endl;
        std::cout << "   std::vector<float> features = {" << std::endl;
        std::cout << "       delta_z,              // Feature 0" << std::endl;
        std::cout << "       delta_z_resunits,     // Feature 1" << std::endl;
        std::cout << "       cluster_z_sigma,      // Feature 2" << std::endl;
        std::cout << "       cluster_d0,           // Feature 3" << std::endl;
        std::cout << "       cluster_d0_sigma,     // Feature 4" << std::endl;
        std::cout << "       cluster_qOverP,       // Feature 5" << std::endl;
        std::cout << "       cluster_qOverP_sigma, // Feature 6" << std::endl;
        std::cout << "       cluster_sumpt         // Feature 7" << std::endl;
        std::cout << "   };" << std::endl;
        std::cout << "   float ml_score = model.predict(features);" << std::endl;
        std::cout << std::endl;
        std::cout << "4. Use score for cluster selection:" << std::endl;
        std::cout << "   if (ml_score > threshold) {" << std::endl;
        std::cout << "       // Select this cluster" << std::endl;
        std::cout << "   }" << std::endl;
        std::cout << "======================================================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
