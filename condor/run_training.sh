#!/bin/bash
# ---------------------------------------------------------------------------
# run_training.sh — condor job executable for python/train_deepsets.py.
#
#   run_training.sh <script.py> <input-dir> <out-dir> [extra args for the script...]
#
# <script.py> is transferred in by the .sub (train_deepsets.py or
# train_transformer.py) and run from the job's scratch directory.
#
# TORCH COMES FROM AN LCG VIEW ON CVMFS, not from a venv in $HOME. The LCG views
# already ship everything this needs, version-consistent with the ROOT the rest
# of the jobs lsetup:
#     LCG_110 x86_64-el9-gcc15-opt : python 3.13.11, torch 2.11.0 (CPU),
#                                    uproot 5.7.1, numpy 2.4.4, pandas 2.2.3
#     LCG_110_cuda x86_64-el9-gcc13-opt : the same but torch built against CUDA
#                                    12.5, for a GPU slot
# That avoids a ~1 GB venv in $HOME, and avoids the trap that a venv built from
# the container python (~/venvs/vtx) is a dangling symlink outside the Jupyter
# container -- which is exactly how it failed on the login node and on workers.
#
# Override with LCG_VIEW=<path> (or VTX_PYTHON=<python> to bypass views entirely).
# ---------------------------------------------------------------------------
set -euo pipefail

# The AF login shell sources an ATLAS/LCG setup exporting PYTHONHOME/PYTHONPATH
# for ITS python. `getenv = True` carries those to the execute node, where they
# would send any other interpreter looking for its stdlib in the wrong prefix:
#   Fatal Python error: init_fs_encoding ... No module named 'encodings'
# Clear them BEFORE sourcing the view, which then sets its own correctly.
unset PYTHONHOME PYTHONPATH

SCRIPT=$1
INPUT_DIR=$2
OUT_DIR=$3
shift 3
[ -f "${SCRIPT}" ] || { echo "FATAL: ${SCRIPT} not transferred in"; exit 1; }

echo "host: $(hostname)   date: $(date)"
echo "script:    ${SCRIPT}"
echo "input-dir: ${INPUT_DIR}"
echo "out-dir:   ${OUT_DIR}"

# Prefer a CUDA view when the slot actually has a GPU, so the same submit file
# works for both without editing.
if [ -z "${LCG_VIEW:-}" ]; then
  if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_110_cuda/x86_64-el9-gcc13-opt
    echo "GPU detected -> CUDA view"
  else
    LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc15-opt
    echo "no GPU -> CPU view"
  fi
fi

PY="${VTX_PYTHON:-}"
if [ -z "${PY}" ]; then
  [ -r "${LCG_VIEW}/setup.sh" ] || { echo "FATAL: no LCG view at ${LCG_VIEW}"; exit 1; }
  echo "sourcing ${LCG_VIEW}/setup.sh"
  # LCG setup scripts reference unset internals and are not strict-mode safe
  set +eu
  source "${LCG_VIEW}/setup.sh"
  set -euo pipefail
  PY="$(command -v python3)"
fi
[ -n "${PY}" ] || { echo "FATAL: no python found"; exit 1; }
echo "python: ${PY} ($(${PY} -V 2>&1))"

# Verify the stack BEFORE the long run, so a bad view costs seconds not hours.
"${PY}" - <<'PY' || { echo "FATAL: torch/uproot not importable with ${PY}"; exit 1; }
import torch, uproot, numpy, pandas
print(f"  torch {torch.__version__} (cuda build {torch.version.cuda}, "
      f"available {torch.cuda.is_available()})  uproot {uproot.__version__}")
PY

mkdir -p "${OUT_DIR}"
echo "running: ${PY} ${SCRIPT} --input-dir ${INPUT_DIR} --out ${OUT_DIR} $*"
"${PY}" "${SCRIPT}" --input-dir "${INPUT_DIR}" --out "${OUT_DIR}" "$@"
