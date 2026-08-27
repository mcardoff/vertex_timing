# Cluster-selection ML: the agreed plan, and how to run it

Nine-item priority order agreed 2026-08-27. This file is the runbook and the
status board; `results/README.md` holds the previous round's numbers (and the
banner explaining why they are not directly comparable).

## Status

| # | Item | State |
|---|------|-------|
| 1 | Re-export / retrain at the canonical selection | exported + merged + staged; training **running** (cluster 2936092) |
| 2 | Freeze the Deep Sets architecture | **done** — it is the default in `train_deepsets.py` |
| 3 | 60 ps timing criterion as the selection target | **done** — was already the label |
| 4 | Core fraction at 60 ps as the primary metric | **done** — was already the metric |
| 5 | Z+jets learning curve | code ready, `condor/train_learning_curve.sub` |
| 6 | Add μ=200 ttbar, test topology transfer | export running; `condor/train_transfer.sub` |
| 7 | Deep Sets vs XGBoost error correlation | code ready, `python/compare_models.py` |
| 8 | μ=0 samples as control/floor, not training data | **done** — excluded from `SAMPLES` |
| 9 | Reopen architecture only on a concrete reason | **done** — encoded in `parse_args` |

## First canonical-selection baselines (2026-08-27)

Reference selectors on the test fold of the re-exported canonical samples. These
are the numbers every later result is measured against, and the first ones taken
at `--vbs-mjj=500` on data carrying the `file_idx` group key.

| sample | events | TRKPTZ | WAVeS | oracle | headroom |
|---|---:|---:|---:|---:|---:|
| vbf   | 519,115 | 90.4% | 91.1% | 98.3% | 7.9 |
| **zjets** | 36,531 | **61.2%** | 59.4% | **92.6%** | **31.4** |
| dijet | 80,383 | 86.6% | 86.2% | 97.7% | 11.1 |
| ttbar | 561,745 | 87.0% | 86.3% | 98.0% | 11.0 |

**Z+jets is still the headroom sample, by a wide margin.** ttbar was a candidate
to displace it — it has 15x the events — but at 11.0 points of headroom it
behaves like dijet, not like Z+jets. So the learning curve (item 5) stays on
Z+jets, and ttbar's value is as a transfer topology (item 6), not as the sample
with room to improve.

**WAVeS trails TRKPTZ on every sample except VBF**, now including ttbar
(86.3 vs 87.0). ttbar is a third independent confirmation of the topology-bound
behaviour, on a sample whose forward jets are scarce for a different reason than
Z+jets'.

## Run-to-run variance, decomposed (2026-08-27)

Every remaining study compares runs, so this is the band any claimed effect must
clear. Measured on canonical, 15 epochs, four runs per row.

| variance source | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| infrastructure — 4x identical seed | **0.000** | **0.000** | **0.000** | **0.000** |
| init only — one split, 4 init seeds | 0.3 | 0.5 | 0.5 | 0.2 |
| split + init — 4 seeds | 1.2 | **1.9** | 0.6 | 0.6 |

**Condor is bit-deterministic**: four runs at one seed agree to every printed
digit. So there is no infrastructure noise to budget for at all.

**The variance is dominated by the SPLIT, not the init** — which events land in
train versus test matters several times more than where the optimiser starts.
That has two consequences:

- The learning curve (item 5) must run at several SPLIT seeds per fraction. A
  single-seed curve on Z+jets would be reading steps of ~1 point against a
  1.9-point band. Three seeds put the error on the mean near 0.55, so a ~1.1
  point step becomes resolvable.
- **Seed ensembling is structurally unable to help here, and measured null.**
  Members must share a split or member B trains on member A's test events, so an
  ensemble can only exploit init variance — the small term. Four members gave
  -0.0 / +0.1 / +0.0 / +0.1 against the best single member. The premise (large
  usable variance) was real, but it lives in the axis ensembling cannot touch
  without leaking. Exploiting it would need split-bagging with a held-out fold
  none of the members saw, which costs test statistics; not pursued.

## Epoch selection: fixed, and worth real points

Validation core fraction was selecting the epoch and is unreliable — on VBF it
collapsed to 17.6% and DEGRADED with training while the same weights scored
91.7% on test. The slice is sound: its oracle matches test's to 0.2 points and
TRKPTZ pushed through the identical tensor scores 90.0%. The fault is in the
argmax-of-model-scores step on that slice, and the ROOT CAUSE REMAINS OPEN — a
control column now prints TRKPTZ through every tensor so a recurrence is caught
immediately.

Selection now uses validation LOSS, which is the training objective itself,
reads the whole score distribution rather than only its argmax, and moves
smoothly. Same seed, same config, same machine:

| selection | epoch | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|---:|
| val macro (old) | 1 | 91.7 | 65.4 | 89.2 | 89.8 |
| **val loss (new)** | 5 | **92.8** | **67.0** | **89.9** | **90.2** |

Better on every sample. Core fraction is still reported per epoch, because the
two disagreeing is the signature of the open bug.

Items 3, 4 and 8 needed no change: the label was already `|Δt| < 60 ps`, the
metric was already core fraction, and the μ=0 exclusion went in with the
group-key fix. They are listed because *verifying* them was part of the plan,
and because item 3 is the one place the plan under review got it wrong — see
"Why the label is not purity" below.

## 1. Re-export at the canonical selection

Seven samples × two selections. Shard counts track events **read**, not events
accepted: the loop runs over every event and only then applies the selection.

```bash
cd build && cmake .. && make export_training_data     # build BEFORE submitting
cd ../condor
for s in "vbf 10" "zjets 13" "dijet 5" "ttbar 20" "vbf_mu0 12" "zeejets_mu0 4" "ttbar_mu0 12"; do
  set -- $s
  condor_submit -a sample=$1 -a nshards=$2 -a 'selargs=--vbs-mjj=500' -a seltag=mjj500p0_ export_training_data_sharded.sub
  condor_submit -a sample=$1 -a nshards=$2 -a 'selargs=--vbs-deta=-1'  -a seltag=novbs_    export_training_data_sharded.sub
done
```

Then merge and stage:

```bash
./condor/merge_exports.sh                    # hadd, NOT hist_merge — see the script
mkdir -p /data/mcardiff/training
cp condor/*/*_training.root /data/mcardiff/training/
```

`merge_exports.sh` globs on the shard **count**, not on `shard*`: a sample
directory can hold two different shardings of the same data (`condor/zjets/`
already did, at 32 and 12), and the loose glob double-counts every event.

## 2 & 9. The architecture is frozen

The defaults in `train_deepsets.py`'s `parse_args` **are** the frozen
configuration — run 1's 24k-parameter model, `pool=sum`, `hidden=96`,
`embed=64`, auxiliary head off. Four levers were scanned and every one came
back null or negative:

| lever | result |
|---|---|
| capacity (5.7× parameters) | **0.2 lower** on all three samples, 44 s/epoch vs 25 |
| pooling (gated vs sum) | +0.1 on Z+jets — inside run-to-run noise |
| auxiliary per-track head | within noise at every λ |
| cluster features / tagger features | null |

The only change that ever moved the number was the **training objective**
(multi-positive listwise), which is not a flag — it is the loss.

Reopening a lever means overriding a default, which is deliberately possible.
Item 9 is that you do it because a study gave you a reason, not on spec.

## 3. Why the label is not purity

The reviewed plan proposed labelling a cluster HS when enough of its tracks come
from the true HS vertex. That is cluster **purity**, and it is the one
construction this project has a standing warning against.

Purity measures where a track *came from*, not whether its time is *right*. A
lone HS track carrying a mis-assigned HGTD hit forms a 100%-pure cluster sitting
at the wrong time, so ranking on purity selects **for** timing misassignment —
which the failure decomposition puts at 67.4% of all losses. Measured directly:
max-purity selection drives HS jets to R_pT = 0 in 11.8% of cases against 3.0%
for the shipped score.

The label is `|delta_t| < 60 ps`. The exporter emits `delta_t` as the target for
exactly this reason and keeps `truth_purity` as a diagnostic only.

## 4. Why core fraction and not AUC

AUC is **anti-correlated** with the figure of merit here. Across five tagger
variants it saturated at 1.000 while core fraction still had 6.4 points to
climb, and the ordering between the two inverted.

Two rules that follow:

- **Quote per sample, never pooled.** VBF starts near 90% with ~8 points of
  headroom and hides everything; a truth feature worth ~14 points on Z+jets
  ablated to exactly 0.0 on VBF, and the pooled 88.9% hides the entire Z+jets
  story.
- **A scan must clear its own noise.** Picking a winner by bare argmax once
  reported "the auxiliary head helped" on a +0.006 margin against a 0.05
  within-run sd.

## 5. Z+jets learning curve

```bash
cd condor && condor_submit train_learning_curve.sub
```

Five fractions of the Z+jets training quota, everything else held fixed.

**Read the curve off `train_events_fitted` in `results.json`, not off the
requested fraction.** Z+jets contributed every available train-fold event last
round, so a requested 1.0 can silently mean "all there was"; plotted against the
request that shows a flat top which is a quota artefact, not saturation — the
exact thing the study is trying to measure.

| outcome | what it means |
|---|---|
| still climbing at 1.0 | more Z+jets MC is the highest-value next step; the `novbs` export is the cheap way to get some without generating anything |
| flat well before 1.0 | Z+jets is not data-limited; the remaining ~24 points need new information, not more events |

## 6. Topology transfer, and what ttbar adds

```bash
cd condor && condor_submit train_transfer.sub
```

Each job fits on one sample and evaluates on all, so one job is one **row** of
the matrix. The train budget is divided over the samples being fit on, not over
all samples loaded, so a single-sample row trains on the full `--train-events`
rather than a fifth of it — otherwise a bad off-diagonal could just mean less
data.

This is the test WAVeS fails: it assumes forward jets are hard-scatter objects,
and on Z+jets the forward dijet is ~94% pileup-supplied, so it gains +0.7 on VBF
and loses 1.9 on Z+jets. **ttbar at μ=200 stresses the same assumption through a
different mechanism** — real heavy-flavour tracks in a busy final state rather
than a nearly-empty forward region — so it is a third independent leg rather
than a repeat. Its first smoke run showed a 20.5% acceptance against Z+jets'
1.5%, so it will contribute real statistics.

## 7. Do Deep Sets and the GBDT fail on the same events?

```bash
python python/compare_models.py --input-dir /data/mcardiff/training \
    --model runs/canonical --selection canonical --out runs/canonical/compare.json
```

The decisive number is the **union ceiling**: the fraction of events where at
least one model is right. A combiner can only pick between what the two offer,
so if the union sits barely above the better model alone, item 7 is closed no
matter how any stack performs.

Reported next to a concrete combiner (per-event rank-average) and φ of the two
error indicators. φ rather than raw agreement, because agreement is dominated by
the easy events both get right — on VBF that is ~90% of them, so two models
could agree 90% of the time while failing on disjoint events.

Uses the **shipped checkpoint**, not a reimplementation, and reproduces the
event split from the same `--seed`, so both models are read on identical test
events through the same argmax-within-event rule.

## 8. μ=0 samples are controls, not training data

`SAMPLES` in `train_deepsets.py` deliberately excludes `vbf_mu0`,
`zeejets_mu0` and `ttbar_mu0`. At μ=0 every method scores ~99.9% core fraction
on every process: one vertex means one cluster and a trivial label, so they
contribute essentially no gradient while diluting the mix with free wins.

The first landed export confirms it independently — `vbf_mu0` produced
**1.1 clusters per passing event**.

This is where the problem differs from the HGTD track–vertex-association BDT,
which *does* train on μ=0 usefully. Do not copy that sample mix over.

They are still worth exporting: they measure the intrinsic timing floor, which
is what turns "the residual needs new detector information" from an elimination
argument into a measurement.

## Things not to re-learn

- **The group key is `(sample_id, file_idx, event_num)`.** Each shard restarts
  `event_num` at 0, so the older two-column key merges unrelated events into one
  listwise group once shards are merged — and every row still looks valid.
- **`hadd` for the exports, `hist_merge` for the `*_hist` files.** Opposite
  rules, same reason: `TParameter` scalars have no `Merge()`, and the exports
  have none while the histogram files do.
- **A negative `--vbs-deta` is a third mode**, not a smaller threshold. No
  combination of the two thresholds drops the VBS *pair* requirement, because
  the pair's `dEta`/`mjj` keep a −1 sentinel when no pair exists.
- **Every number in `results/README.md` predates the canonical selection.**
