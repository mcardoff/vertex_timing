#!/bin/bash
# ---------------------------------------------------------------------------
# probe_env.sh — 30-second preflight for the Deep Sets training jobs.
#
# The pool's execute hosts do NOT share a filesystem with the submit host (see
# CLAUDE.md: a `TARGET.FileSystemDomain == MY.FileSystemDomain` requirement
# matched 0 of 157 slots). But the ntuples ARE read from an absolute
# /data/mcardiff/... path by the existing jobs, so /data is clearly mounted --
# which leaves it genuinely unknown whether $HOME (and therefore ~/venvs/vtx)
# is visible from a worker.
#
# Rather than guess and burn a long job, this reports exactly what a worker can
# see. Run it first; its .out tells you which torch strategy the real submit
# file should use.
# ---------------------------------------------------------------------------
echo "=== host ==="
hostname; date; echo "cwd: $(pwd)"; echo "whoami: $(whoami)"; echo "HOME=$HOME"

echo; echo "=== cpu / memory ==="
nproc 2>/dev/null || echo "nproc unavailable"
grep -E "^MemTotal" /proc/meminfo 2>/dev/null || true

echo; echo "=== GPU ==="
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi -L
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
else
  echo "no nvidia-smi -- CPU-only slot"
fi

echo; echo "=== is \$HOME visible (venv + repo)? ==="
for p in "$HOME" "$HOME/venvs/vtx/bin/python" "$HOME/project/vertex_timing"; do
  [ -e "$p" ] && echo "  OK      $p" || echo "  MISSING $p"
done

echo; echo "=== is /data visible (ntuples + staged training files)? ==="
for p in /data /data/mcardiff /data/mcardiff/exotic_superntuples /data/mcardiff/training; do
  [ -e "$p" ] && echo "  OK      $p" || echo "  MISSING $p"
done
echo "  training ROOT files found:"
ls -lh /data/mcardiff/training/*_training.root 2>/dev/null | sed 's/^/    /' || echo "    none"

echo; echo "=== python / torch ==="
echo "system python3: $(command -v python3) $(python3 -V 2>&1)"
if [ -x "$HOME/venvs/vtx/bin/python" ]; then
  "$HOME/venvs/vtx/bin/python" - <<'PY' 2>&1 | sed 's/^/  venv: /'
import sys
print("python", sys.version.split()[0])
for m in ("torch", "uproot", "numpy", "pandas"):
    try:
        mod = __import__(m); print(f"{m} {getattr(mod,'__version__','?')}")
    except Exception as e:
        print(f"{m} MISSING ({type(e).__name__})")
try:
    import torch; print("cuda available:", torch.cuda.is_available())
    if torch.cuda.is_available(): print("device:", torch.cuda.get_device_name(0))
except Exception:
    pass
PY
else
  echo "  venv python not reachable from this node"
fi

echo; echo "=== apptainer/singularity available? (fallback if HOME is not shared) ==="
command -v apptainer >/dev/null 2>&1 && apptainer --version || echo "  no apptainer"
command -v singularity >/dev/null 2>&1 && singularity --version || echo "  no singularity"

echo; echo "=== PROBE COMPLETE ==="
