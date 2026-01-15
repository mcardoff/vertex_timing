// Test script to verify ONNX model can be loaded with ROOT SOFIE
// Usage: root -l -b -q test_onnx_model.C

#include "TMVA/SOFIE_common.hxx"
#include "TMVA/RModelParser_ONNX.hxx"
#include <iostream>

void test_onnx_model() {
    using namespace TMVA::Experimental;

    std::cout << "=== ONNX Model Test ===" << std::endl;
    std::cout << "Loading model.onnx..." << std::endl;

    try {
        // Parse the ONNX model
        SOFIE::RModelParser_ONNX parser;
        SOFIE::RModel model = parser.Parse("model.onnx");

        std::cout << "\n✓ Model loaded successfully!" << std::endl;

        // Print model information
        std::cout << "\nModel Information:" << std::endl;
        std::cout << "  Input tensors: " << model.GetInputTensorNames().size() << std::endl;
        for (const auto& name : model.GetInputTensorNames()) {
            std::cout << "    - " << name << std::endl;
            auto shape = model.GetInputTensorShape(name);
            std::cout << "      Shape: [";
            for (size_t i = 0; i < shape.size(); i++) {
                std::cout << shape[i];
                if (i < shape.size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }

        std::cout << "\n  Output tensors: " << model.GetOutputTensorNames().size() << std::endl;
        for (const auto& name : model.GetOutputTensorNames()) {
            std::cout << "    - " << name << std::endl;
        }

        // Generate C++ code
        std::cout << "\nGenerating C++ inference code..." << std::endl;
        model.Generate();
        model.OutputGenerated("model_inference.hxx");

        std::cout << "✓ Generated model_inference.hxx" << std::endl;
        std::cout << "\nNext steps:" << std::endl;
        std::cout << "1. Include the generated header in your analysis code:" << std::endl;
        std::cout << "   #include \"model_inference.hxx\"" << std::endl;
        std::cout << "2. Create a session and run inference (see ONNX_INTEGRATION_GUIDE.md)" << std::endl;

    } catch (const std::exception& e) {
        std::cout << "✗ Error loading model: " << e.what() << std::endl;
        std::cout << "\nTroubleshooting:" << std::endl;
        std::cout << "- Check that model.onnx exists in the current directory" << std::endl;
        std::cout << "- Verify ROOT was compiled with SOFIE support (ROOT >= 6.26)" << std::endl;
        std::cout << "- Try using ONNX Runtime instead (see ONNX_INTEGRATION_GUIDE.md)" << std::endl;
    }
}
