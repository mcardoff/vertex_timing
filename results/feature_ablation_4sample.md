# Feature ablation, 4-sample canonical (2026-08-28)

Paired A/B: identical `--seed` (same split) AND identical `--init-seed`, so the
only difference between arms is which columns reach the model. 5 seeds per arm,
20 runs. Every run asserted to have trained on the same sample set --
an earlier attempt was voided because dijet arrived mid-experiment and 5 runs
saw 2 samples while 11 saw 3.

Data: the 2026-08-28 re-export, canonical selection (`--vbs-mjj=500`).
vbf 519,115 ev / zjets 36,531 / dijet 80,383 / ttbar 561,745.

## Core fraction

| arm | dijet | ttbar | vbf | zjets |
|---|---|---|---|---|
| base (22 feat)   | 90.00 | 90.30 | 93.30 | 66.94 |
| + vertex (29)    | -0.02 | -0.06 (0/5) | -0.10 (0/5) | **+0.34 +/- 0.11 (5/5)** |
| + jet/loo (25)   | -0.04 | -0.08 (0/5) | **+0.12 (5/5)** | **-0.16 (0/5)** |
| full (32)        | -0.08 | -0.04 | +0.10 | **+0.40 +/- 0.16 (4/5)** |
| TRKPTZ           | 86.60 | 87.00 | 90.40 | 61.20 |
| oracle           | 97.70 | 98.00 | 98.30 | 92.60 |

`base` reproduces the published 67.1 +/- 1.0 Z+jets baseline (66.94), so the
harness is consistent before anything is layered on.

## The finding: the two groups are sample-ANTAGONISTIC

- **vertex** (dz_to_nearest_pu_vtx, pv_pu_dz_ratio, nearest_vtx_waves_rank/frac,
  closest_vtx_is_pv, and the two per-track forms) helps ONLY Z+jets, and is
  mildly negative everywhere else -- 0/5 seeds up on both vbf and ttbar.
- **jet/loo** (is_in_any_jet, dr_nearest_anyjet, loo_pull) does the reverse:
  +0.12 on vbf with 5/5 up, -0.16 on Z+jets with 0/5 up.

Both are consistent across every seed in the direction stated, so neither is
noise. The mechanism explains it: vertex relations matter only where the failure
is picking a displaced pure-pileup cluster (65.9% of Z+jets recoverable
failures, rare elsewhere), and jet association matters only where the hard
scatter IS the forward jets, which is VBF and not Z+jets.

This is why the local-VBF ablation was misleading. There, jet/loo carried
everything (+0.38) and vertex was a clean null on top (-0.02, 0/5) -- the
correct reading of a sample that does not have the failure the vertex features
target.

## Size

+0.40 on a 25.7-point Z+jets gap closes ~1.6% of it. Real, reproducible, small.
The remaining gap is still dominated by the 13.2-point time-quality term (which
no selector reaches) plus a selection problem the rank distribution characterises
as top-3 discrimination, not search.

## Open, suggested by the antagonism

Each group helps exactly one sample and costs the others. A single feature set
is therefore a compromise. Worth testing whether a per-sample or gated model
beats the pooled one -- but note the transfer result already showed the learned
selector generalises across topologies where WAVeS does not, so specialising
trades that away.
