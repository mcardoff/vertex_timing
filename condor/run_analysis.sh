#!/bin/bash
# ---------------------------------------------------------------------------
# run_analysis.sh — condor job executable template for vertex_timing.
#
# Invoked by clustering_dt.sub as:
#   run_analysis.sh <executable> <sample>
#     <executable>  clustering_dt | rpt_v5   (any --sample-aware build/ target)
#     <sample>      vbf | zjets | dijet | default   ("default" = no --sample flag,
#                    i.e. the local ../../ntuple-hgtd/ + ../figs/ behavior)
#
# Assumes a shared filesystem between the condor submit and execute hosts
# (true on the UChicago AF pool) — PROJECT_DIR must already contain a build/
# directory built with `cd build && cmake .. && make` *before* submitting.
# This script does not rebuild, so concurrently-queued jobs don't race on
# the same build/ directory.
# ---------------------------------------------------------------------------
set -euo pipefail

# EDIT ME: absolute path to the vertex_timing checkout on the shared filesystem.
PROJECT_DIR=/home/mcardiff/project/vertex_timing

EXECUTABLE=$1
SAMPLE=$2

# ATLAS/LCG environment (provides ROOT + Boost via cvmfs, matching the
# `lsetup root` assumption baked into CMakeLists.txt's cvmfs discovery).
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"
lsetup "root 6.38.04-x86_64-el9-gcc15-opt"

cd "${PROJECT_DIR}/build"

if [ "${SAMPLE}" = "default" ]; then
  ./"${EXECUTABLE}"
else
  ./"${EXECUTABLE}" --sample="${SAMPLE}"
fi
