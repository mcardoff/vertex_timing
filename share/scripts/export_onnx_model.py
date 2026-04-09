#!/usr/bin/env python3
"""
export_onnx_model.py — ML model conversion pipeline for vertex_timing.

Converts a TensorFlow SavedModel into the three artifacts consumed by the
C++ inference engine (ml_model.h):

  1. saved_model_dir  →  model.onnx            (via tf2onnx)
  2. model.onnx       →  model_weights.json    (stable w1/b1/w2/b2 keys)
  3. normalization_params.json  →  clustering_structs.h MEANS/STDS (in-place)
                                →  share/models/ (copy)
  4. model_weights.json         →  share/models/ (install)

Usage:
    python export_onnx_model.py [--saved-model-dir DIR] [--dest-dir DIR]

Defaults:
    --saved-model-dir  ~/project/clusterclassification/saved_model_dir
    --dest-dir         ~/project/vertex_timing/share/models

Dependencies (install into a venv if needed):
    pip install tensorflow tf2onnx onnx numpy

After running, rebuild vertex_timing:
    cd ~/project/vertex_timing/build && make

Notes:
    - Weight keys are STABLE (w1/b1/w2/b2/…) regardless of what tf2onnx
      names the internal graph nodes.  Layer order is inferred from tensor
      shapes: 2-D matrices are chained by shape[0] → shape[1]; 1-D bias
      vectors are matched by size.
    - clustering_structs.h is edited in-place to update MEANS[N] / STDS[N].
    - The ONNX file is written next to the SavedModel, not into share/models,
      because it is large (~50 KB) and only needed at conversion time.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths (overridable via CLI)
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent          # share/scripts/
_PROJECT_ROOT = _HERE.parent.parent              # vertex_timing/
_DEFAULT_SRC = Path.home() / "project" / "clusterclassification" / "saved_model_dir"
_DEFAULT_DST = _PROJECT_ROOT / "share" / "models"
_STRUCTS_H   = _PROJECT_ROOT / "src" / "clustering_structs.h"


# ---------------------------------------------------------------------------
# Step 1: SavedModel → ONNX (uses the tf2onnx CLI to avoid TF import issues)
# ---------------------------------------------------------------------------
def step1_convert_to_onnx(saved_model_dir: Path, onnx_out: Path) -> None:
    print("=" * 60)
    print("Step 1: SavedModel → ONNX  (tf2onnx)")
    print("=" * 60)

    if not saved_model_dir.is_dir():
        raise FileNotFoundError(f"SavedModel directory not found: {saved_model_dir}")

    cmd = [
        sys.executable, "-m", "tf2onnx.convert",
        "--saved-model", str(saved_model_dir),
        "--output",      str(onnx_out),
        "--opset",       "13",
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("\n--- stdout ---")
        print(result.stdout)
        print("--- stderr ---")
        print(result.stderr)
        raise RuntimeError("tf2onnx conversion failed (see output above)")

    size_kb = onnx_out.stat().st_size / 1024
    print(f"✓  Wrote {onnx_out}  ({size_kb:.1f} KB)")


# ---------------------------------------------------------------------------
# Step 2: ONNX → model_weights.json with stable w1/b1/w2/b2/… keys.
#
#   Layer ordering: chain 2-D weight matrices by shape (first dim of each
#   layer matches the last dim of the previous layer, starting from n_input).
#   Biases are matched by output size (= corresponding weight's shape[1]).
#   Values are flattened row-major, matching ml_model.h: W[i*N_OUT + j].
# ---------------------------------------------------------------------------
def step2_export_weights(onnx_out: Path, json_out: Path) -> None:
    print("\n" + "=" * 60)
    print("Step 2: ONNX → model_weights.json")
    print("=" * 60)

    try:
        import onnx
        from onnx import numpy_helper
    except ImportError:
        raise ImportError("onnx not found — pip install onnx")

    model = onnx.load(str(onnx_out))
    graph = model.graph

    # Input feature count
    n_input = graph.input[0].type.tensor_type.shape.dim[-1].dim_value
    if n_input <= 0:
        raise RuntimeError(f"Could not read input feature count from ONNX graph (got {n_input})")

    # Separate 2-D weights from 1-D biases
    weights_2d, biases_1d = {}, {}
    for init in graph.initializer:
        t = numpy_helper.to_array(init)
        if t.ndim == 2:
            weights_2d[init.name] = t
        elif t.ndim == 1:
            biases_1d[init.name] = t
        else:
            print(f"  WARNING: skipping initializer {init.name} shape={t.shape}")

    # Chain weight matrices by shape
    ordered_weights, remaining, current_size = [], dict(weights_2d), n_input
    while remaining:
        matches = {n: t for n, t in remaining.items() if t.shape[0] == current_size}
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected exactly 1 weight matrix with {current_size} input rows, "
                f"found {len(matches)}: {list(matches.keys())}"
            )
        name, tensor = next(iter(matches.items()))
        ordered_weights.append((name, tensor))
        current_size = tensor.shape[1]
        del remaining[name]

    # Match biases by output size
    ordered_biases, remaining_b = [], dict(biases_1d)
    for _, w in ordered_weights:
        out_size = w.shape[1]
        matches = {n: t for n, t in remaining_b.items() if t.shape[0] == out_size}
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected exactly 1 bias vector of size {out_size}, "
                f"found {len(matches)}: {list(matches.keys())}"
            )
        name, tensor = next(iter(matches.items()))
        ordered_biases.append((name, tensor))
        del remaining_b[name]

    # Emit with stable names
    print(f"\n  Input features: {n_input}")
    print(f"  Detected {len(ordered_weights)} dense layers:")
    output, total_vals = {}, 0
    for i, ((wn, wt), (bn, bt)) in enumerate(zip(ordered_weights, ordered_biases), 1):
        wkey, bkey = f"w{i}", f"b{i}"
        output[wkey] = {"shape": list(wt.shape), "values": wt.flatten().tolist()}
        output[bkey] = {"shape": list(bt.shape), "values": bt.flatten().tolist()}
        total_vals += wt.size + bt.size
        print(f"    Layer {i}: {wt.shape[0]} → {wt.shape[1]}")
        print(f"      {wkey} ← {wn}  |  {bkey} ← {bn}")

    print(f"\n  Total parameters: {total_vals:,}")

    with open(json_out, "w") as f:
        json.dump(output, f, indent=2)

    size_kb = json_out.stat().st_size / 1024
    print(f"✓  Wrote {json_out}  ({size_kb:.1f} KB)")


# ---------------------------------------------------------------------------
# Step 3: Update MEANS/STDS arrays in clustering_structs.h and copy JSON
# ---------------------------------------------------------------------------
def step3_update_normalization(norm_json: Path, dest_dir: Path) -> None:
    print("\n" + "=" * 60)
    print("Step 3: Update normalization in clustering_structs.h")
    print("=" * 60)

    if not norm_json.exists():
        raise FileNotFoundError(f"normalization_params.json not found: {norm_json}")

    params = json.loads(norm_json.read_text())
    means, stds, names = params["means"], params["stds"], params["feature_names"]
    n = len(means)

    if len(stds) != n or len(names) != n:
        raise ValueError(f"normalization_params.json inconsistency: "
                         f"means={len(means)}, stds={len(stds)}, names={len(names)}")

    print(f"\n  {n} features:")
    for nm, m, s in zip(names, means, stds):
        print(f"    {nm:<28}  mean={m:>20.13g}  std={s:>20.13g}")

    def build_block(arr_name: str, values: list, kind: str) -> str:
        lines = [f"      static const float {arr_name}[{n}] = {{"]
        for i, (v, nm) in enumerate(zip(values, names)):
            comma = "," if i < n - 1 else ""
            lines.append(f"        {repr(v)}{comma}  // {nm} {kind}")
        lines.append("      };")
        return "\n".join(lines)

    means_block = build_block("MEANS", means, "mean")
    stds_block  = build_block("STDS",  stds,  "std")

    src = _STRUCTS_H.read_text()

    pat_means = re.compile(r'[ \t]*static const float MEANS\[\d+\] = \{[^}]*\};', re.DOTALL)
    pat_stds  = re.compile(r'[ \t]*static const float STDS\[\d+\]\s*= \{[^}]*\};', re.DOTALL)

    for pat, name in [(pat_means, "MEANS"), (pat_stds, "STDS")]:
        if len(pat.findall(src)) != 1:
            raise RuntimeError(f"Expected exactly 1 {name} array in {_STRUCTS_H}")

    src = pat_means.sub(means_block, src)
    src = pat_stds.sub(stds_block,   src)

    _STRUCTS_H.write_text(src)
    print(f"\n✓  Updated {_STRUCTS_H}")

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_norm = dest_dir / "normalization_params.json"
    shutil.copy2(str(norm_json), str(dest_norm))
    print(f"✓  Copied → {dest_norm}")


# ---------------------------------------------------------------------------
# Step 4: Install model_weights.json → share/models/
# ---------------------------------------------------------------------------
def step4_install(json_out: Path, dest_dir: Path) -> None:
    print("\n" + "=" * 60)
    print(f"Step 4: Install model_weights.json → {dest_dir}")
    print("=" * 60)

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / "model_weights.json"
    shutil.copy2(str(json_out), str(dest))
    print(f"✓  Copied → {dest}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--saved-model-dir", type=Path, default=_DEFAULT_SRC,
        help=f"TF SavedModel directory (default: {_DEFAULT_SRC})"
    )
    parser.add_argument(
        "--dest-dir", type=Path, default=_DEFAULT_DST,
        help=f"Destination for model artifacts (default: {_DEFAULT_DST})"
    )
    args = parser.parse_args()

    saved_model_dir: Path = args.saved_model_dir
    dest_dir: Path        = args.dest_dir

    # ONNX file lives alongside the SavedModel (large, not checked in)
    onnx_out  = saved_model_dir.parent / "model.onnx"
    json_out  = saved_model_dir.parent / "model_weights.json"
    norm_json = saved_model_dir / "normalization_params.json"

    step1_convert_to_onnx(saved_model_dir, onnx_out)
    step2_export_weights(onnx_out, json_out)
    step3_update_normalization(norm_json, dest_dir)
    step4_install(json_out, dest_dir)

    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print(f"  Artifacts installed to: {dest_dir}")
    print("\nNext step:")
    print(f"  cd {_PROJECT_ROOT / 'build'} && make")
    print("=" * 60)


if __name__ == "__main__":
    main()
