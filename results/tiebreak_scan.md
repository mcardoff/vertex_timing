# What else could separate the failure mode -- population-first scan

`python/tiebreak_scan.py` + `python/plot_candidate_vars.py`
(`figs/candidate_vars.pdf`). After the jet-fraction lesson, this scan never
conditions on losing: every number is computed over the whole population of both
samples, and a candidate must clear three bars --

1. **population AUC** (identify the HS-dominant cluster among all clusters);
2. **tie-break accuracy** on the reco-defined ambiguous pairs (top-2 by our
   score) on BOTH samples, against the score's own accuracy on the same pairs
   (VBF 95.32%, Z+jets 81.04%);
3. **switch accounting**: re-rank the ambiguous pair by the candidate with ONE
   global direction and threshold, count better vs worse.

## The ranking (top non-score variables)

| variable | pop VBF | pop zjets | tie VBF | tie zjets |
|---|---:|---:|---:|---:|
| **dt_cluster_to_hgtd** | 0.977 | **0.901** | 0.957 | **0.853** |
| cluster_d0_sigma | 0.975 | 0.884 | 0.949 | 0.809 |
| cluster_z_sigma | 0.970 | 0.879 | 0.939 | 0.805 |
| dt_to_agg_cl (classical aggregate) | 0.969 | 0.857 | 0.949 | 0.800 |
| (score's own accuracy on the pairs) | -- | -- | 0.953 | 0.810 |

Everything else -- time/qOverP sigmas, nhgtd fractions, waves_score -- sits at or
below the score's own accuracy on at least one sample.

## Switch accounting: exactly one survivor

| variable | f | VBF net (bet:wor) | zjets net (bet:wor) |
|---|---|---|---|
| **dt_cluster_to_hgtd** | 0.5 | **+0.48 (464:275)** | **+0.43 (1704:1304)** |
| **dt_cluster_to_hgtd** | 0.8 | **+0.29 (213:96)** | **+0.31 (811:521)** |
| cluster_d0_sigma | 0.5 | +0.04 (510:493) | -0.07 (3277:3339) |
| cluster_z_sigma | 0.5 | -0.32 (471:599) | +0.24 (3512:3287) |

`cluster_d0_sigma` and `cluster_z_sigma` ranked well on AUC and are coin flips
when actually switched -- the same pattern that killed the jet term, caught here
before implementation.

**`dt_cluster_to_hgtd` -- agreement with Athena's own HGTD vertex time
(`RecoVtx_time`) -- is the one variable that survives on both samples**, with
better:worse ratios of 1.7-2.2 on VBF and 1.3-1.6 on Z+jets. The exporter sets it
NaN when the HGTD time is invalid (checked at `export_training_data.cxx:1196`),
so the Z+jets gain is earned inside the 16.5% of events where Athena produced a
time -- it is real corroboration, not a garbage-zero artifact; where the time is
absent the rule simply never fires, which is what makes it sample-independent.

Physical reading: `RecoVtx_time` is an INDEPENDENT reconstruction of the vertex
time. A candidate cluster that agrees with it is corroborated by an algorithm
that shares none of our clustering decisions -- the same "agreement with an
independent estimate" mechanism that made `dt_tagw_agg` the ML selector's top
feature, but with stronger independence and no model.

## Plots

`figs/candidate_vars.pdf`: nine panels, four curves each ({VBF, Z+jets} x
{HS-dominant, other}), unit-normalised, per-sample AUC and direction annotated.
The rejected `frac_pt_in_fwdjet` is included as the contrast case: visibly
separating on VBF and nearly overlapping on Z+jets.
