#!/usr/bin/env python3
"""
Export trained Keras model to ONNX format.
Run this after training your model with the 8 features.
"""

import tensorflow as tf
import tf2onnx
import onnx

# Assuming you have your trained model in 'neural_net'
# If you need to load it:
# neural_net = tf.keras.models.load_model('path_to_saved_model')

# Export to ONNX
model_proto, _ = tf2onnx.convert.from_keras(
    neural_net,
    input_signature=[tf.TensorSpec(shape=(None, 8), dtype=tf.float32, name='keras_tensor')],
    opset=13
)

# Save the ONNX model
onnx.save(model_proto, "model.onnx")

print("✓ Model exported to model.onnx")
print("  Input shape: (batch_size, 8)")
print("  Output shape: (batch_size, 1)")

# Verify the exported model
onnx_model = onnx.load("model.onnx")
onnx.checker.check_model(onnx_model)
print("✓ ONNX model verification passed")
