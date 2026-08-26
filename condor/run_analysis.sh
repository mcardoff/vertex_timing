#!/bin/bash
# ---------------------------------------------------------------------------
# run_analysis.sh — condor job executable template for vertex_timing.
#
# Invoked by clustering_hist.sub as:
#   run_analysis.sh <executable> <sample> [<threads>] [<max-events>] [--flag ...]
#
# Trailing arguments are classified by SHAPE, not position: anything starting
# with "--" is forwarded to the executable verbatim; anything else fills
# <threads> then <max-events>. So a flag-only job needs no empty placeholders.
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
#     --flag         optional; any argument starting with "--" is forwarded
#                    verbatim. Used for --vbs-deta=<x> (loosened VBS topology
#                    selection, see resolveSelection in src/sample_config.h),
#                    which also tags the output file name so a loosened run
#                    cannot overwrite the standard one.
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
shift 2

# Remaining args are classified by shape rather than by position: anything
# starting with "--" is a flag forwarded verbatim, anything else fills <threads>
# then <max-events> in order. This means a job that wants only a flag can write
#   arguments = export_training_data $(sample) --vbs-deta=0
# instead of needing empty positional placeholders -- condor's argument parser
# rejects bare "" ("Found illegal unescaped double-quote"), and its new-style
# quoting ('' inside an enclosing "...") is easy to get wrong. Existing
# invocations like `... clustering_hist vbf 4` are unaffected.
THREADS=""
MAX_EVENTS=""
EXTRA_ARGS=()
for a in "$@"; do
  case "$a" in
    --*) EXTRA_ARGS+=("$a") ;;
    "")  ;;                                  # tolerate stray empty placeholders
    *)   if   [ -z "${THREADS}" ];    then THREADS="$a"
         elif [ -z "${MAX_EVENTS}" ]; then MAX_EVENTS="$a"
         else EXTRA_ARGS+=("$a"); fi ;;
  esac
done

# Clear the positional parameters now that they're captured in THREADS/
# MAX_EVENTS/EXTRA_ARGS/EXECUTABLE/SAMPLE: `source file` with no explicit
# trailing args leaves the CALLING script's "$@" unchanged, so without this,
# atlasLocalSetup.sh below inherits this script's own leftover args (e.g.
# "4 --vbs-deta=0 --vbs-mjj=200") as if they were its own CLI options. It
# errors "unrecognized option '--vbs-deta=0'" on any "--" flag and returns
# without defining `lsetup`, which then cascades into "lsetup: command not
# found", no ROOT environment, and the transferred executable failing to
# start with "error while loading shared libraries: libtbb.so.12" -- no
# output file is ever written, so the job holds on output-transfer failure
# with no indication the real cause was upstream. Silent for any job whose
# leftover args happened to be flag-free (a bare "4"), which is why this went
# unnoticed until the first job that combined a numeric <threads> AND a
# --vbs-* flag in the same submission.
set --

# Elapsed-seconds stamps around each phase. The ATLAS/cvmfs setup below and
# the analysis itself are charged to the same job, and separating them is the
# only way to attribute a slow job: startup dead time has been measured at
# 5-30 minutes and varies ~4x run-to-run on the SAME sample, which no
# per-sample property explains. Without these stamps the only evidence
# available was condor's 5-minute image-size sampling in the .log file.
T0=$SECONDS
stamp() { printf '[job t+%4ds] %s\n' "$((SECONDS - T0))" "$1"; }
stamp "job start on $(hostname -s)"

# ATLAS/LCG environment (provides ROOT + Boost via cvmfs). atlasLocalSetup.sh
# / lsetup reference unset variables internally (e.g. ALRB_frontlineSite) and
# aren't `set -e`/`set -u` safe, so relax those flags around them.
set +eu
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"
lsetup "root 6.38.04-x86_64-el9-gcc15-opt"   # match the version used to build
set -euo pipefail
stamp "ATLAS/cvmfs environment ready" 

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

stamp "running: ./${EXECUTABLE} ${SAMPLE_ARG} ${THREADS_ARG} ${MAX_EVENTS_ARG} ${EXTRA_ARGS[*]:-}"
./"${EXECUTABLE}" ${SAMPLE_ARG} ${THREADS_ARG} ${MAX_EVENTS_ARG} ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
stamp "analysis exited" 
