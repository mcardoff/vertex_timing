#!/usr/bin/env python3
"""
Inspect ONNX model structure and optionally evaluate it or export parameters.

Usage:
    python inspect_onnx_model.py --model model.onnx [--export-weights] [--test]
"""

import argparse
import sys

def inspect_model(model_path):
    """Inspect ONNX model structure"""
    try:
        import onnx
        from onnx import numpy_helper
    except ImportError:
        print("ERROR: onnx package not installed")
        print("Install with: pip install onnx")
        return False

    print("=" * 70)
    print(f"ONNX MODEL INSPECTION: {model_path}")
    print("=" * 70)

    model = onnx.load(model_path)

    print("\nMODEL INFORMATION:")
    print(f"  IR Version: {model.ir_version}")
    print(f"  Producer: {model.producer_name} {model.producer_version}")
    print(f"  Domain: {model.domain}")

    # Get the graph
    graph = model.graph

    # Input information
    print("\nINPUT TENSORS:")
    for inp in graph.input:
        print(f"  Name: {inp.name}")
        shape = [dim.dim_value if dim.dim_value > 0 else 'dynamic'
                 for dim in inp.type.tensor_type.shape.dim]
        print(f"    Shape: {shape}")
        print(f"    Type: {inp.type.tensor_type.elem_type}")

    # Output information
    print("\nOUTPUT TENSORS:")
    for out in graph.output:
        print(f"  Name: {out.name}")
        shape = [dim.dim_value if dim.dim_value > 0 else 'dynamic'
                 for dim in out.type.tensor_type.shape.dim]
        print(f"    Shape: {shape}")
        print(f"    Type: {out.type.tensor_type.elem_type}")

    # Layer information
    print(f"\nMODEL ARCHITECTURE:")
    print(f"  Total operators: {len(graph.node)}")

    op_types = {}
    for node in graph.node:
        op_types[node.op_type] = op_types.get(node.op_type, 0) + 1

    print("  Operator counts:")
    for op_type, count in sorted(op_types.items()):
        print(f"    {op_type}: {count}")

    # Initializers (weights)
    print(f"\nMODEL PARAMETERS:")
    print(f"  Total initializers: {len(graph.initializer)}")
    total_params = 0
    for init in graph.initializer:
        tensor = numpy_helper.to_array(init)
        total_params += tensor.size
        print(f"    {init.name}: shape={tensor.shape}, size={tensor.size}")

    print(f"\n  Total parameters: {total_params:,}")

    return True

def export_weights(model_path, output_path="model_weights.json"):
    """Export model weights to JSON for manual C++ implementation"""
    try:
        import onnx
        from onnx import numpy_helper
        import json
        import numpy as np
    except ImportError:
        print("ERROR: Required packages not installed")
        print("Install with: pip install onnx numpy")
        return False

    print("\n" + "=" * 70)
    print("EXPORTING MODEL WEIGHTS")
    print("=" * 70)

    model = onnx.load(model_path)
    graph = model.graph

    weights = {}

    # Extract initializers (weights and biases)
    for init in graph.initializer:
        tensor = numpy_helper.to_array(init)
        weights[init.name] = {
            'shape': list(tensor.shape),
            'values': tensor.flatten().tolist()
        }

    # Export to JSON
    with open(output_path, 'w') as f:
        json.dump(weights, f, indent=2)

    print(f"\n✓ Exported weights to {output_path}")
    print(f"  File size: {len(json.dumps(weights)) / 1024:.1f} KB")

    return True

def test_inference(model_path):
    """Test model inference with dummy data"""
    try:
        import onnxruntime as ort
        import numpy as np
    except ImportError:
        print("ERROR: onnxruntime not installed")
        print("Install with: pip install onnxruntime")
        return False

    print("\n" + "=" * 70)
    print("TESTING MODEL INFERENCE")
    print("=" * 70)

    # Create inference session
    session = ort.InferenceSession(model_path)

    # Get input details
    input_name = session.get_inputs()[0].name
    input_shape = session.get_inputs()[0].shape

    print(f"\nInput: {input_name}")
    print(f"  Expected shape: {input_shape}")

    # Create dummy input (batch_size=1, n_features)
    n_features = input_shape[1] if len(input_shape) > 1 else input_shape[0]
    dummy_input = np.random.randn(1, n_features).astype(np.float32)

    print(f"  Using dummy input shape: {dummy_input.shape}")

    # Run inference
    outputs = session.run(None, {input_name: dummy_input})

    print(f"\nOutput shape: {outputs[0].shape}")
    print(f"Output values: {outputs[0]}")

    print("\n✓ Model inference successful!")

    return True

def create_cpp_inference_wrapper(model_path):
    """Create a Python script that can be called from C++ for inference"""
    wrapper_path = "evaluate_model.py"

    wrapper_code = '''#!/usr/bin/env python3
"""
Model evaluation wrapper for C++ integration.
Reads features from stdin, outputs prediction to stdout.

Usage from C++:
    echo "feat1,feat2,feat3,..." | python evaluate_model.py
"""

import sys
import numpy as np
import onnxruntime as ort

def main():
    # Load model
    session = ort.InferenceSession("MODEL_PATH")

    # Read features from stdin
    line = sys.stdin.readline().strip()
    features = [float(x) for x in line.split(',')]

    # Prepare input
    input_name = session.get_inputs()[0].name
    input_array = np.array([features], dtype=np.float32)

    # Run inference
    outputs = session.run(None, {input_name: input_array})

    # Output prediction
    print(outputs[0][0])

if __name__ == "__main__":
    main()
'''.replace("MODEL_PATH", model_path)

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_code)

    import os
    os.chmod(wrapper_path, 0o755)

    print(f"\n✓ Created inference wrapper: {wrapper_path}")
    print("\nC++ integration example:")
    print("""
    // In your C++ code:
    #include <cstdio>
    #include <string>

    float evaluate_model(const std::vector<float>& features) {
        // Build feature string
        std::string feat_str;
        for (size_t i = 0; i < features.size(); i++) {
            feat_str += std::to_string(features[i]);
            if (i < features.size() - 1) feat_str += ",";
        }

        // Call Python script
        std::string cmd = "echo '" + feat_str + "' | python evaluate_model.py";
        FILE* pipe = popen(cmd.c_str(), "r");
        float score;
        fscanf(pipe, "%f", &score);
        pclose(pipe);

        return score;
    }
    """)

def main():
    parser = argparse.ArgumentParser(
        description="Inspect and work with ONNX models",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--model", default="model.onnx", help="Path to ONNX model")
    parser.add_argument("--export-weights", action="store_true",
                       help="Export weights to JSON")
    parser.add_argument("--test", action="store_true",
                       help="Test inference with dummy data")
    parser.add_argument("--create-wrapper", action="store_true",
                       help="Create Python wrapper for C++ integration")

    args = parser.parse_args()

    # Always inspect
    if not inspect_model(args.model):
        return 1

    # Optional exports
    if args.export_weights:
        export_weights(args.model)

    if args.test:
        test_inference(args.model)

    if args.create_wrapper:
        create_cpp_inference_wrapper(args.model)

    print("\n" + "=" * 70)
    print("RECOMMENDATIONS FOR C++ INTEGRATION:")
    print("=" * 70)
    print("Since ROOT doesn't have SOFIE support, you have 3 options:")
    print()
    print("1. INSTALL ONNX RUNTIME (Recommended):")
    print("   brew install onnxruntime")
    print("   Then update CMakeLists.txt (see ONNX_INTEGRATION_GUIDE.md)")
    print()
    print("2. USE PYTHON WRAPPER:")
    print("   Run with --create-wrapper to generate evaluate_model.py")
    print("   Call it from C++ using popen() (slower but simple)")
    print()
    print("3. MANUAL IMPLEMENTATION:")
    print("   Run with --export-weights to get model_weights.json")
    print("   Implement the network manually in C++ (for simple models)")
    print("=" * 70)

if __name__ == "__main__":
    sys.exit(main() or 0)
