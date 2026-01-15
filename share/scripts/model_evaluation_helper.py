#!/usr/bin/env python3
"""
Helper script to export Python ML models for evaluation in C++/ROOT

Usage:
    python model_evaluation_helper.py --model-type [onnx|tmva|params] --model-file model.pkl

This script provides code snippets and exports for evaluating your Python-trained
model in the C++ vertex timing analysis.
"""

import argparse
import json
import sys

def export_to_onnx(model_file, output_file="model.onnx"):
    """Export model to ONNX format for ROOT SOFIE evaluation"""
    print("=" * 70)
    print("OPTION 1: ONNX Export (Recommended for Neural Networks)")
    print("=" * 70)

    print("\n1. First, install required packages:")
    print("   pip install onnx skl2onnx  # For scikit-learn")
    print("   pip install torch onnx      # For PyTorch")

    print("\n2. Python export code:")
    print("""
# For PyTorch model:
import torch
model.eval()
dummy_input = torch.randn(1, n_features)
torch.onnx.export(model, dummy_input, "{output}")

# For scikit-learn model:
import pickle
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType

with open("{model}", "rb") as f:
    model = pickle.load(f)

# Define input shape
initial_type = [('float_input', FloatTensorType([None, n_features]))]
onx = convert_sklearn(model, initial_types=initial_type)

with open("{output}", "wb") as f:
    f.write(onx.SerializeToString())
    """.format(model=model_file, output=output_file))

    print("\n3. C++ evaluation code (requires ROOT 6.26+):")
    print("""
// In your analysis code:
#include "TMVA/SOFIE_common.hxx"
#include "TMVA/RModelParser_ONNX.hxx"

// One-time setup
TMVA::Experimental::SOFIE::RModelParser_ONNX parser;
TMVA::Experimental::SOFIE::RModel model = parser.Parse("model.onnx");
model.Generate();
model.OutputGenerated("model.hxx");

// Include generated header
#include "model.hxx"

// In event loop:
TMVA_SOFIE_model::Session s("model.dat");
std::vector<float> input = {{feature1, feature2, feature3}};
std::vector<float> output = s.infer(input.data());
float score = output[0];
    """)

def export_tmva_compatible(model_file):
    """Show how to use TMVA Reader"""
    print("=" * 70)
    print("OPTION 2: TMVA Reader (For models trained in ROOT)")
    print("=" * 70)

    print("\nThis is the method already used in your codebase!")
    print("Check existing usage in clustering_functions.h")

    print("\nC++ evaluation code:")
    print("""
// Setup (once per analysis):
TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

// Add variables (must match training)
Float_t var1, var2, var3;
reader->AddVariable("var1_name", &var1);
reader->AddVariable("var2_name", &var2);
reader->AddVariable("var3_name", &var3);

// Book MVA method
reader->BookMVA("MethodName", "weights/TMVA.weights.xml");

// In event loop:
var1 = value1;
var2 = value2;
var3 = value3;
Float_t score = reader->EvaluateMVA("MethodName");
    """)

def export_parameters(model_file, output_file="model_params.json"):
    """Export model parameters for manual C++ implementation"""
    print("=" * 70)
    print("OPTION 3: Parameter Export (For simple models)")
    print("=" * 70)

    print("\n1. Python export code:")
    print("""
import pickle
import json
import numpy as np

with open("{model}", "rb") as f:
    model = pickle.load(f)

# For LogisticRegression, LinearRegression, etc:
params = {{
    'type': type(model).__name__,
    'coefficients': model.coef_.tolist(),
    'intercept': float(model.intercept_),
    'features': feature_names  # Add your feature names
}}

with open("{output}", "w") as f:
    json.dump(params, f, indent=2)

# For tree-based models (more complex):
# Consider using export_text or tree structure export
    """.format(model=model_file, output=output_file))

    print("\n2. C++ evaluation code:")
    print("""
#include <fstream>
#include <nlohmann/json.hpp>  // or use ROOT's TBufferJSON

// Load once
std::ifstream f("model_params.json");
nlohmann::json params = nlohmann::json::parse(f);

// For logistic regression:
double predict(const std::vector<double>& features) {{
    double score = params["intercept"];
    auto coeffs = params["coefficients"];
    for (size_t i = 0; i < features.size(); i++) {{
        score += coeffs[i].get<double>() * features[i];
    }}
    return 1.0 / (1.0 + std::exp(-score));  // sigmoid
}}

// In event loop:
std::vector<double> features = {{var1, var2, var3}};
double score = predict(features);
    """)

def check_existing_tmva():
    """Check if TMVA weights already exist"""
    import glob
    import os

    weights = glob.glob("TMVA*.xml")
    if weights:
        print("\n" + "=" * 70)
        print("Found existing TMVA weights:")
        print("=" * 70)
        for w in weights:
            print(f"  - {w} ({os.path.getsize(w) / 1024:.1f} KB)")
        print("\nYou're already using TMVA! See OPTION 2 above.")
        print("Check clustering_functions.h for existing evaluation code.")

def main():
    parser = argparse.ArgumentParser(
        description="Export Python ML model for C++/ROOT evaluation",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--model-type",
        choices=["onnx", "tmva", "params", "all"],
        default="all",
        help="Export method to use"
    )
    parser.add_argument(
        "--model-file",
        default="model.pkl",
        help="Path to Python model file (pickle)"
    )
    parser.add_argument(
        "--output",
        help="Output file name (default: auto-generated)"
    )

    args = parser.parse_args()

    print("\nMODEL EVALUATION HELPER")
    print("=" * 70)

    # Check for existing TMVA weights
    check_existing_tmva()

    if args.model_type == "all":
        print("\n")
        export_to_onnx(args.model_file)
        print("\n")
        export_tmva_compatible(args.model_file)
        print("\n")
        export_parameters(args.model_file)
    elif args.model_type == "onnx":
        export_to_onnx(args.model_file, args.output or "model.onnx")
    elif args.model_type == "tmva":
        export_tmva_compatible(args.model_file)
    elif args.model_type == "params":
        export_parameters(args.model_file, args.output or "model_params.json")

    print("\n" + "=" * 70)
    print("RECOMMENDATIONS:")
    print("=" * 70)
    print("1. If you already have TMVA XML weights → Use OPTION 2 (TMVA Reader)")
    print("2. For new neural network models → Use OPTION 1 (ONNX)")
    print("3. For simple linear/logistic models → Use OPTION 3 (Parameter export)")
    print("\nNeed help? Check ROOT TMVA documentation:")
    print("https://root.cern.ch/doc/master/group__tutorial__tmva.html")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    main()
