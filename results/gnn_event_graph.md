# The event as a graph: message passing helps the TIME, not the SELECTION

`python/train_gnn.py`, VBF, local export, 2026-09-09. First run of the graph idea
against the cluster-selection problem, with the `--layers 0` ablation and a
5-seed paired sweep that settles what message passing is actually worth.

**Headline: message passing is worth +0.80 +- 0.14 on the timing head and
+0.26 +- 0.10 on selection, both 5/5 seeds. The graph's value is in deciding how
much to trust each track's time, not in choosing the cluster.**

## Why a graph at all

Not capacity -- that is a measured null (5.7x parameters scored 0.2 LOWER, see
`train_deepsets.py`). The motivation is that ONE hand-built round of message
passing was already the second most productive thing in the ML study.
`final_method_eval.py:148-165` joins all track pairs within an event, gates the
neighbour's HS probability by hard boxes, and sums:

    nb_phzt = SUM_{j != i} ph_j * (|dz_ij| < 2.0) * (|dt_ij| < 60)

That is one GNN layer with a fixed, non-learned edge kernel. Its six columns take
37% of the round-2 tagger's gain (`nb_phzt` alone 23.76 against 25.98 for the top
feature), and `nb_ph60` is the #2 feature in the entire e2e t0 model. This asks what
the hand version leaves behind: raw-|dt| boxes that discard per-track resolution,
one hop, degree-dominated sum aggregation, and a frozen round-1 -> round-2 bootstrap.

## Design: three things this does NOT do

1. **It does not re-partition.** Nodes are tracks in the EVENT, edges are all-pairs
   with two types (same-cluster / cross-cluster), so the partition is a prior the
   model can override, never the topology. Building edges only inside clusters
   would disconnect the graph and blind it to selection, which is the dominant
   failure channel. Re-partitioning is separately and repeatedly measured to lose:
   `WAVES_RECLUST` -13.5 on zjets (`baseline_deta0p0_mjj500p0.md`), split
   clustering -3.14 / 25.34% at -2.61 real complementarity
   (`split_clustering.md`), misclustering 3.1% of principal causes. The mechanism
   is statistics: at ~29 timed tracks and a 60 ps window, a soft weight can push a
   contaminant to 0.05 and keep the other 28, while a partition boundary throws it
   out and takes the statistics with it.
2. **It does not learn the averaging.** `t_hat = SUM w_i t_i / SUM w_i` stays
   structural. The no-clustering transformer collapsed on VBF (77%, plateaued by
   epoch 4) and the diagnosis was that it had to LEARN the averaging cluster
   methods get for free; the untested fix named there was "emit per-track weights,
   output SUM w_i t_i". That is what this is.
3. **It does not start from scratch.** Both heads are residual from the incumbent:

       score_c = tau * log(trkptz_score_c) + f(pool_c)     f last layer zero-init
       w_i     = (1/sigma_t,i^2) * exp(g(h_i))             g last layer zero-init

   so at epoch 0 the model IS TRKPTZ selection with inverse-variance timing,
   exactly. **This is the methodology contribution and it should be reused**: the
   TRKPTZ control column that `train_deepsets.py` prints as a diagnostic becomes
   structural, and "did it help?" is read against a verified zero rather than a
   remembered baseline. It caught a real error on the first run (see Methodology).

## The ablation (seed 0, 11,000 test events, local VBF export)

| run | feats | layers | kernel | sel | + learned t0 |
|---|---:|---:|---|---:|---:|
| TRKPTZ (these events) | -- | -- | -- | 91.08 | -- |
| **A** reference | 48 | 2 | learned | **94.22** | **95.20** |
| **B** old-export feature set | 22 | 2 | learned | 93.98 | 94.90 |
| **C** no message passing | 48 | 0 | -- | 93.99 | 94.38 |
| **D** frozen hand-built gate | 48 | 2 | box | 94.15 | 94.90 |
| oracle ceiling | | | | 98.39 | |

| isolated by | selection | timing |
|---|---:|---:|
| message passing (A-C) | +0.23 | +0.82 |
| the 26 restored features (A-B) | +0.24 | +0.30 |
| learned vs hand-built gate (A-D) | +0.06 | +0.30 |

**The edge kernel's SHAPE barely matters; having edges does.** D reproduces A to
+0.06 on selection. An earlier control that compared learned vs box alone was
therefore the wrong test and nearly produced the wrong conclusion -- box is still a
full GNN with only the gate frozen. `--layers 0` is the control that answers it.

## The paired 5-seed sweep

A and C at MATCHED seeds: `--seed` fixes the event split and `--init-seed` defaults
to it, so each replicate shares split, events and batch order and the difference is
paired.

| seed | TRKPTZ | A sel | C sel | d sel | A full | C full | d full |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 91.08 | 94.22 | 93.99 | +0.23 | 95.20 | 94.38 | +0.82 |
| 1 | 90.65 | 94.00 | 93.60 | +0.40 | 94.99 | 94.06 | +0.93 |
| 2 | 90.91 | 94.20 | 93.87 | +0.33 | 95.16 | 94.29 | +0.87 |
| 3 | 90.35 | 93.63 | 93.45 | +0.17 | 94.82 | 93.98 | +0.84 |
| 4 | 90.75 | 93.66 | 93.49 | +0.17 | 94.95 | 94.39 | +0.56 |

    message passing, selection : +0.26 +- 0.10   5/5 positive
    message passing, timing    : +0.80 +- 0.14   5/5 positive

**Pairing was load-bearing and the data proves it.** TRKPTZ's OWN core fraction
varies by sd 0.28 across these seeds -- pure event-split variance in a fixed
algorithm. The selection effect (+0.26) is smaller than the split noise it sits in,
so unpaired it is unmeasurable. This is the local, VBF-sized instance of the
standing rule that run-to-run variance is split-dominated (zjets +-1.9 split+init
vs +-0.5 init-only): **pair on the split, or do not compare.**

## The mechanism

The learned per-track weight head is worth **+0.54 +- 0.21 alone and +1.08 +- 0.15
with the graph**. Message passing roughly DOUBLES what learned weights are worth.

That is the whole result, and it is mechanistically forced: `timeRes` separates HS
from PU at **AUC 0.508, a coin flip**, yet it is the only weight the
inverse-variance mean uses. A track's own features cannot say whether its time is
trustworthy; its neighbourhood can. The graph is a trust estimator, not a selector.

Decomposition of the total, means over the 5 seeds:

| step | delta | share |
|---|---:|---:|
| TRKPTZ 90.75 -> 0-layer model 93.68 (pooling + listwise + residual head) | +2.93 | 68% |
| -> message passing on selection 93.94 | +0.26 | 6% |
| -> graph-informed per-track weights 95.03 | +1.08 | 25% |
| **total** | **+4.28** | |

Two thirds of the gain is NOT the graph. It is the listwise objective and the
pooled selector -- consistent with the study's finding that the training objective
was the only lever that ever moved the number.

## What is NOT established

- **VBF only.** Ceiling 98.29, TRKPTZ 90.75 -- the narrowest headroom of the four
  samples, and the one where the transformer collapsed rather than where the
  problem lives. zjets (61.87) decides anything.
- **The timing gain may be a VBF fact, not a physics fact.** VBF carries 47% of HS
  timing weight inside jets against zjets' 10% (`injet_decomposition.md`), and that
  is exactly the quantity governing whether jet-neighbourhood information means
  anything. Expect the +0.80 to SHRINK on zjets. That is the falsifying test.
- **Small sample.** 23,092 fit events (the local export is 39,681 events total,
  vs 519,115 in the grid `highstats_vbf` export). Absolute numbers here are not
  comparable with `ml_overview.md`'s ladder.
- **`--no-residual` never run.** Whether the incumbent initialisation helps the
  final number, or only makes it readable, is untested.

## Methodology notes

- **The residual init caught a real bug on its first use.** Epoch 0 read 89.3%
  against a 90.4% TRKPTZ reference and the check fired -- but the fault was the
  CONTROL: the reference covered the whole 155k-event test fold while the tensor
  held 1,000 events, where binomial spread at 90% is ~1 point. Fixed by computing
  TRKPTZ on the tensor's own events (`tensor_ref()`). Now exact to -0.000 on all
  10 sweep runs. A control compared against the wrong denominator is worse than no
  control.
- **The time loss must be scaled.** In raw picoseconds the Huber term outweighed
  the listwise CE ~1000:1 (loss 1525 vs O(1)) and the selection head was barely
  training. It is now expressed in units of PASS_PS with delta=0.5 (== 30 ps), the
  same physical transition point.
- **Two exports of "VBF" are different productions.** The local 33-file sample
  (39,681 events) and `/data/mcardiff/exotic_superntuples/highstats_vbf/`
  (519,115) are not the same data, so a feature ablation across the two files
  confounds the feature set with the sample. Hence `--drop-features`, which does
  the A/B on one file at one seed.

## Reproduction

```bash
# export, local VBF, canonical selection (4 shards then hadd -- NOT hist_merge)
cd build
for i in 0 1 2 3; do ./export_training_data --sample=vbf \
   --ntuple-dir=/Users/mcard/project/ntuple-hgtd/ --vbs-mjj=500 --file-shard=$i/4 & done; wait
cd vbf && hadd -f vbf_mjj500p0_training.root vbf_mjj500p0_training.shard*.root

# the ablation (local box: PYTHONNOUSERSITE=1 and ~/.venv-hgtd, see CLAUDE.md)
PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/train_gnn.py --input-dir build \
   --samples vbf --train-events 27000 --test-events 11000 --epochs 12 \
   --batch-events 256 --torch-threads 4 --out runs/A            # A
#  ... --layers 0                                               # C
#  ... --edge-kernel box                                        # D
#  ... --drop-features "<the 26>"                               # B
```

Runs live in `runs/` (gitignored: the checkpoints are binary blobs).
