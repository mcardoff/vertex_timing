#!/bin/bash
# ---------------------------------------------------------------------------
# merge_exports.sh — hadd each sample's export shards into one file per
# (sample, selection).
#
#   ./condor/merge_exports.sh [<condor-dir>]      # default: condor/
#
# WHY hadd AND NOT util/hist_merge. This is the inverse of the rule for the
# *_hist outputs, and the reason is what the files contain. hist_merge exists
# because those carry TParameter scalars (event counts, pT accumulators) which
# have no Merge() method, so hadd would silently keep the FIRST input's copy and
# pair N-shard-summed histograms with one shard's event count. A training export
# holds only TTrees, which hadd concatenates correctly -- and hist_merge would
# refuse them outright, since it recognises histograms and scalars only.
#
# GLOB ON THE SHARD COUNT, NOT ON shard*. A sample directory can hold two
# different shardings of the same data (condor/zjets/ already did once, at 32
# and 12); "shard*" matches both sets and double-counts every event. This script
# therefore discovers the (tag, count) pairs present and merges each separately,
# rather than globbing loosely.
#
# Re-running is safe: hadd -f overwrites, and the merged file's own name carries
# no shard tag so it can never be picked up as an input.
# ---------------------------------------------------------------------------
set -euo pipefail

DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
command -v hadd >/dev/null || { echo "hadd not on PATH -- lsetup root first"; exit 1; }

shopt -s nullglob
merged=0

for sampledir in "$DIR"/*/; do
  sample=$(basename "$sampledir")
  [ "$sample" = "logs" ] && continue

  # Discover every distinct "<tag>training.shard<i>of<N>" family in this
  # directory: key on the part before ".shard" plus the "of<N>" count, so two
  # shardings of the same selection stay separate.
  declare -A fam=()
  for f in "$sampledir"/*_training.shard*of*.root; do
    base=$(basename "$f")
    prefix="${base%%.shard*}"                 # <sample>_<tag>training
    count="${base##*of}"; count="${count%.root}"
    fam["$prefix|$count"]=1
  done
  [ ${#fam[@]} -eq 0 ] && continue

  for key in "${!fam[@]}"; do
    prefix="${key%|*}"; count="${key#*|}"
    out="$sampledir/$prefix.root"
    inputs=("$sampledir/$prefix".shard*of"$count".root)
    if [ ${#inputs[@]} -ne "$count" ]; then
      echo "  !! $prefix: ${#inputs[@]} of $count shards present -- SKIPPING"
      echo "     (a merge of an incomplete set silently loses those events)"
      continue
    fi
    hadd -f "$out" "${inputs[@]}" > /dev/null 2>&1
    ev=$(python3 - "$out" <<'PY' 2>/dev/null || echo "?"
import sys, uproot
f = uproot.open(sys.argv[1])
c = f["clusters"].arrays(["sample_id","file_idx","event_num"], library="pd")
print(f'{len(c.groupby(["sample_id","file_idx","event_num"])):,}')
PY
)
    printf "  %-46s %2s shards -> %s events\n" "$prefix" "$count" "$ev"
    merged=$((merged+1))
  done
  unset fam
done

echo
echo "merged $merged file set(s)."
echo "Point the trainer at this directory:"
echo "  python python/train_deepsets.py --input-dir $DIR --selection canonical --out runs/canonical"
