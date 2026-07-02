#!/bin/bash
# ---------------------------------------------------------------------------
# run_analysis.sh — condor job executable template for vertex_timing.
#
# Invoked by clustering_dt.sub as:
#   run_analysis.sh <executable> <sample> [<threads>]
#     <executable>  clustering_dt | rpt_v5   (any --sample/--threads-aware build/ target)
#     <sample>      vbf | zjets | dijet | default   ("default" = no --sample flag,
#                    i.e. the local ../../ntuple-hgtd/ + ../figs/ behavior)
#     <threads>     optional; forwarded as --threads=<N>. Should match this
#                    job's request_cpus in the .sub file -- condor's cgroup
#                    throttles any threads spawned beyond what was requested,
#                    which would silently eat the whole parallelization
#                    benefit without any visible error. Omit to fall back to
#                    the executable's own default (min(hardware_concurrency(), 8)).
#
# No shared filesystem between submit and execute hosts: <executable>
# arrives in this job's scratch directory via transfer_input_files. It's
# moved into a build/ subdir so the program's own "../<sample>" output path
# convention (mirroring `cd build && ./clustering_dt`) lands at the scratch
# root, where transfer_output_files can stage it back.
# ---------------------------------------------------------------------------
set -euo pipefail

EXECUTABLE=$1
SAMPLE=$2
THREADS=${3:-}

# ATLAS/LCG environment (provides ROOT + Boost via cvmfs). atlasLocalSetup.sh
# / lsetup reference unset variables internally (e.g. ALRB_frontlineSite) and
# aren't `set -e`/`set -u` safe, so relax those flags around them.
set +eu
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"
lsetup "root 6.38.04-x86_64-el9-gcc15-opt"   # match the version used to build
set -euo pipefail

mkdir -p build
mv "${EXECUTABLE}" build/
cd build

SAMPLE_ARG=""
if [ "${SAMPLE}" != "default" ]; then
  SAMPLE_ARG="--sample=${SAMPLE}"
fi

THREADS_ARG=""
if [ -n "${THREADS}" ]; then
  THREADS_ARG="--threads=${THREADS}"
fi

./"${EXECUTABLE}" ${SAMPLE_ARG} ${THREADS_ARG}
