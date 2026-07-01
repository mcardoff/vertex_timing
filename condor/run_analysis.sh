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
# No shared filesystem between submit and execute hosts: <executable>
# arrives in this job's scratch directory via transfer_input_files. It's
# moved into a build/ subdir so the program's own "../<sample>" output path
# convention (mirroring `cd build && ./clustering_dt`) lands at the scratch
# root, where transfer_output_files can stage it back.
# ---------------------------------------------------------------------------
set -euo pipefail

EXECUTABLE=$1
SAMPLE=$2

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

if [ "${SAMPLE}" = "default" ]; then
  ./"${EXECUTABLE}"
else
  ./"${EXECUTABLE}" --sample="${SAMPLE}"
fi
