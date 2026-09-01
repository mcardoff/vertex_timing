# Classical scoring: where the value is, and why central information keeps failing

Z+jets (`training_new`, novbs union, 93,977 events; canonical subset 36,531).
Everything here is model-free and truth-free at inference.

## The structural constraint

Cluster selection is an **argmax within an event**. Any purely event-level
quantity multiplies every candidate equally and therefore **cancels exactly**.
This is not a soft effect -- it is why the vertex-WAVeS proposal was dead before
it was measured. Only three things can help:

1. quantities that vary **across the clusters of one event**;
2. quantities inside the **per-track** weight (they reshape, not rescale);
3. a **per-event decision** that is not an argmax -- i.e. a chooser.

## What works classically

| | Z+jets |
|---|---:|
| TRKPTZ | 62.12% |
| **TRKPTZ_TZ** — per-track `exp(-0.7·dz)` | **63.89%** |
| classical aggregate t0 (weighted median + double Winsorisation) | **62.89%** |
| — same, plain weighted median only | 60.02% |
| — same, inverse-variance mean over all tracks | **41.43%** |
| **union(TRKPTZ_TZ cluster, classical t0)** | **72.31%** |

Two things stand out.

**The estimator is the single largest classical lever measured all session.**
Median + double Winsorisation beats the inverse-variance mean by **21.5 points**
on identical tracks (62.89 vs 41.43). The inverse-variance mean is the estimator
most exposed to the ~13.5% of genuine HS tracks sitting >200 ps from truth, and
it is what the code currently uses everywhere.

**A model-free aggregate t0 already beats TRKPTZ** (62.89 vs 62.12) while
choosing no cluster at all.

## The classical union is within 0.6 of the ML union

Canonical subset, like for like:

| | cluster answer | t0 answer | **union** | delivered |
|---|---:|---:|---:|---:|
| classical | 63.63 | 62.55 | **72.25** | 64.85 |
| ML pipeline | 69.43 | 69.01 | **72.87** | 70.00 |

The ML's two answers are each ~6 points better, yet the union is only 0.6 better
-- **the ML components are far more correlated with each other than the
classical ones are.** The classical pair is individually weaker but covers
almost the same events between them.

So the ceiling is essentially reachable classically. What the ML bought was the
ability to *realise* more of it, not a higher ceiling.

## The chooser is the classical prize, and it is hard

+8.4 points sit between the best single classical answer (63.89) and the
classical union (72.31). Single-variable threshold rules, threshold chosen
out-of-fold on half the events and scored on the other half:

| feature | best | vs 63.89 |
|---|---:|---:|
| `abs(t_cluster - t0)` | 64.85 | **+0.96** |
| `cluster_time_sigma` | 64.74 | +0.85 |
| `time_chi2_ndf` | 64.60 | +0.71 |
| `n_tracks` | 64.59 | +0.70 |
| `sumpt` | 64.45 | +0.56 |
| `local_vtx_density` | 64.01 | +0.12 |
| **`n_forward_jets`** | 63.89 | **+0.00** |
| `n_reco_vertices` | 63.84 | -0.05 |

Only ~11% of the available headroom. For scale the ML meta captured ~17% of its
own (smaller) headroom, so **the chooser is hard for everyone** -- this is not a
classical-vs-learned gap.

Note *which* features work: they are all **cluster-quality** measures plus the
separation between the two answers. Every event-topology variable is ~0.

## Central information has now failed four ways

| route | result |
|---|---|
| vertex-WAVeS as a score multiplier | -0.54 to -5.69 |
| `closer_to_pu_than_pv` veto / down-weight | -1.27 / superseded by TZ |
| event topology as chooser input | +0.00 (`n_forward_jets`), -0.05 (`n_reco_vertices`) |
| **event-adaptive exponent** `a = a(event)` | **+0.08 in-sample, ~0 honest** |

The adaptive-exponent test is the cleanest of the four because it is the one
route that survives the cancellation argument. Best `a` per bin, in-sample:

| split | bin | n | core frac | best a | vs a=0.7 |
|---|---|---:|---:|---:|---:|
| `local_vtx_density` | 0.00-1.00 | 54,726 | 66.74 | 0.7 | +0.00 |
| | 1.25 | 22,907 | 61.89 | 0.9 | +0.15 |
| | 1.50-2.25 | 16,344 | 57.58 | 1.6 | +0.24 |
| `n_reco_vertices` | 65-102 | 31,872 | 65.39 | **0.7** | +0.00 |
| | 103-111 | 33,790 | 63.88 | **0.7** | +0.00 |
| | 112-148 | 28,315 | 62.20 | **0.7** | +0.00 |

On `n_reco_vertices` the optimum is 0.7 in **every** bin. The direction on
`local_vtx_density` is at least physical (a denser vertex neighbourhood wants a
tighter weight) and the density strongly predicts *difficulty* -- 66.74% down to
57.34% -- but it does not call for a different weighting. **a = 0.7 is genuinely
global.**

## Ranked list of what is still worth trying, classically

1. **Re-time the chosen cluster with the same weight used to select it.**
   `Cluster::calculateTime` still returns the plain inverse-variance mean. The
   classical analogue of the ML pipeline's step 3 is
   `t = SUM(w_t·t_t/sigma^2) / SUM(w_t/sigma^2)` with
   `w_t = pT·exp(-0.7·dz)`, Winsorised twice. In the ML pipeline the equivalent
   move was worth **+0.9 to +1.7 on every sample** and changed no selection at
   all. Given the 21.5-point estimator effect above, this is the best-motivated
   untested idea here.
2. **A 2-3 variable classical chooser** over `abs(t_cluster-t0)`,
   `cluster_time_sigma` and `n_tracks`. Single variables give +0.96 of +8.4;
   these three are plausibly not redundant.
3. **A better classical t0.** The second Winsorisation pass was worth +2.87
   (60.02 -> 62.89). A proper Huber/Tukey M-estimator, or more passes, is cheap.
4. **Cluster in (t, z) jointly rather than t alone.** The clustering has been
   untouched all session; `useZ0` and the distance-cut machinery already exist.
5. Stop spending effort on central information via the routes above.
