#!/bin/bash
# ---------------------------------------------------------------------------
# run_analysis.sh — condor job executable template for vertex_timing.
#
# Invoked by clustering_hist.sub as:
#   run_analysis.sh <executable> <sample> [<threads>] [<max-events>]
#     <executable>  clustering_hist | rpt_v5_hist   (any --sample/--threads-aware target)
#     <sample>      vbf | zjets | dijet | default   ("default" = no --sample flag,
#                    i.e. the local ../../ntuple-hgtd/ + ../figs/ behavior --
#                    not actually queued by any .sub file today, since condor
#                    jobs always run against a named grid sample)
#     <threads>     optional; forwarded as --threads=<N>. Should match this
#                    job's request_cpus in the .sub file -- condor's cgroup
#                    throttles any threads spawned beyond what was requested,
#                    which would silently eat the whole parallelization
#                    benefit without any visible error. Omit to fall back to
#                    the executable's own default (min(hardware_concurrency(), 8)).
#     <max-events>  optional; forwarded as --max-events=<N> (see resolveMaxEvents
#                    in src/sample_config.h). For a quick sanity-check job over
#                    a small event prefix instead of the full sample -- omit
#                    for a normal production run (unlimited).
#
# No shared filesystem between submit and execute hosts: <executable>
# arrives in this job's scratch directory via transfer_input_files and is run
# directly from there. The grid-sample output paths (vbf/zjets/dijet in
# src/sample_config.h) are resolved relative to the executable's own cwd, not
# its parent, specifically so this never needs to escape the scratch
# directory it's confined to -- see the comment in sample_config.h. (The
# "default" sample above still uses the "../figs" local-dev convention, but
# since no .sub file ever queues it, that's never exercised on condor.)
# ---------------------------------------------------------------------------
set -euo pipefail

EXECUTABLE=$1
SAMPLE=$2
THREADS=${3:-}
MAX_EVENTS=${4:-}

# ATLAS/LCG environment (provides ROOT + Boost via cvmfs). atlasLocalSetup.sh
# / lsetup reference unset variables internally (e.g. ALRB_frontlineSite) and
# aren't `set -e`/`set -u` safe, so relax those flags around them.
set +eu
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"
lsetup "root 6.38.04-x86_64-el9-gcc15-opt"   # match the version used to build
set -euo pipefail

SAMPLE_ARG=""
if [ "${SAMPLE}" != "default" ]; then
  SAMPLE_ARG="--sample=${SAMPLE}"
fi

THREADS_ARG=""
if [ -n "${THREADS}" ]; then
  THREADS_ARG="--threads=${THREADS}"
fi

MAX_EVENTS_ARG=""
if [ -n "${MAX_EVENTS}" ]; then
  MAX_EVENTS_ARG="--max-events=${MAX_EVENTS}"
fi

./"${EXECUTABLE}" ${SAMPLE_ARG} ${THREADS_ARG} ${MAX_EVENTS_ARG}
