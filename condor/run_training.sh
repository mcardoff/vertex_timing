#!/bin/bash
# ---------------------------------------------------------------------------
# run_training.sh — condor job executable for python/train_deepsets.py.
#
#   run_training.sh <input-dir> <out-dir> [extra args for train_deepsets.py...]
#
# Torch is located in this order, so the job works whether or not $HOME is
# mounted on the execute node (see probe_env.sub -- run that first):
#   1. $VTX_PYTHON if set          (explicit override from the .sub file)
#   2. ~/venvs/vtx/bin/python      (works when $HOME is shared)
#   3. python3 on PATH             (works if the base image already has torch)
# It fails loudly with the reason rather than falling back to a python without
# torch and dying halfway through the event loop.
#
# The script is transferred in; input ROOT files are read from an ABSOLUTE path
# (/data/...) rather than transferred, because they are GB-scale and /data is
# mounted on the workers -- the same reason the ntuples are read that way.
# ---------------------------------------------------------------------------
set -euo pipefail

# The AF login shell sources an ATLAS/LCG setup that exports PYTHONHOME and
# PYTHONPATH pointing at LCG's Python 3.13 tree. Any OTHER interpreter launched
# with those set looks for its stdlib in the wrong prefix and dies before it can
# run a line of code:
#   Fatal Python error: init_fs_encoding: failed to get the Python codec of the
#   filesystem encoding / ModuleNotFoundError: No module named 'encodings'
# The .sub files use `getenv = True`, so those variables reach the execute node
# too. Clear them before touching any venv python.
unset PYTHONHOME PYTHONPATH

INPUT_DIR=$1
OUT_DIR=$2
shift 2

echo "host: $(hostname)   date: $(date)"
echo "input-dir: ${INPUT_DIR}"
echo "out-dir:   ${OUT_DIR}"

PY=""
if [ -n "${VTX_PYTHON:-}" ] && [ -x "${VTX_PYTHON}" ]; then
  PY="${VTX_PYTHON}"
elif [ -x "${HOME}/venvs/vtx/bin/python" ]; then
  PY="${HOME}/venvs/vtx/bin/python"
elif command -v python3 >/dev/null 2>&1; then
  PY="$(command -v python3)"
fi
[ -n "${PY}" ] || { echo "FATAL: no python found"; exit 1; }
echo "python: ${PY} ($(${PY} -V 2>&1))"

# Verify torch BEFORE the long run, so a missing dependency costs seconds.
"${PY}" - <<'PY' || { echo "FATAL: torch/uproot not importable with ${PY} -- see probe_env.sub"; exit 1; }
import torch, uproot, numpy, pandas
print(f"  torch {torch.__version__}  cuda={torch.cuda.is_available()}  uproot {uproot.__version__}")
PY

mkdir -p "${OUT_DIR}"
echo "running: ${PY} train_deepsets.py --input-dir ${INPUT_DIR} --out ${OUT_DIR} $*"
"${PY}" train_deepsets.py --input-dir "${INPUT_DIR}" --out "${OUT_DIR}" "$@"
