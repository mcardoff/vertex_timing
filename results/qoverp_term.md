# The momentum-precision term: sigma_qP passes where the jet term failed

User-proposed candidate after rejecting dt_cluster_to_hgtd (agreement with the
Athena time is borrowed-answer information, not independent selection power).

`cluster_qOverP_sigma = 1/sqrt(SUM_t 1/var_qOverP_t)` (export line 1131) -- the
inverse-variance-combined q/P uncertainty. `1/sigma_qP^2` is "amount of
well-measured momentum": it rises with track count AND per-track momentum
quality.

## The battery (`python/qoverp_term.py`), run population-first

**Redundancy.** corr(log sigma_qP, log sumpt) = -0.95 / -0.92, and on pairs with
EQUAL n_tracks its zjets tie-break accuracy collapses to 0.550 -- so it is mostly
a repackaging of "more, better-measured tracks". Honest framing: this is a mild
reshaping of the score toward well-measured clusters, not a new physics handle.

**Switch accounting** (hard pair re-rank): coin flip (-0.05/-0.44 at f=0.5).
As a TIE-BREAKER it is useless -- but unlike the jet term it was not proposed as
one.

**Multiplicative** `score' = tz * (1/sigma_qP)^beta`, one global beta, vs the
controls it might be repackaging:

| best working point | worst-sample delta |
|---|---:|
| **(1/sigma_qP)^0.8** | **+0.21** |
| sumpt^0.4 | +0.18 |
| n_tracks^0.2 (== n_valid_time) | +0.14 |

It narrowly beats both controls because it combines their information in one
number. The three curves are nearly identical, as the -0.95 correlation demands.

## Four-sample validation (novbs, ~2.0M events, one global beta)

| beta | vbf | zjets | dijet | ttbar | worst |
|---|---:|---:|---:|---:|---:|
| 0.40 | +0.19 | +0.26 | +0.20 | +0.19 | **+0.19** |
| 0.80 | +0.26 | +0.24 | +0.25 | +0.21 | **+0.21** |

**Positive on every sample at both working points** -- the first additive term
since TRKPTZ_TZ itself to pass the full consistency battery. Contrast the jet
term (+0.32 vbf / -0.01 zjets) and the hard veto (-1.27 vbf).

## In the main program: `Score::TRKPTZ_TZQ` (id 33)

`SUM_t pT e^{-0.7|dz_t|} * e^{-1.5|dz_cl|} * (SUM_t 1/var_qP)^{0.4}`, guarded
in-jet time. `QP_PRECISION_WEIGHT = 0.8`. Local VBF (mjj>=200), all 48,314
selected events:

| | core frac | vs TRKPTZ |
|---|---:|---:|
| **TRKPTZ_TZQ** | **91.83%** | **+1.47** |
| WAVES (deployed) | 91.71% | +1.36 |
| TRKPTZ_TZ_GIJ | 91.58% | +1.22 |
| TRKPTZ_TZ (raw t) | 91.07% | +0.71 |
| TRKPTZ | 90.36% | -- |

The qP term adds +0.24 in the C++ against the python scan's +0.22-0.26, and
WAVES reproduces at 91.71% exactly.

**This flips the deployable-candidate comparison: TRKPTZ_TZQ now beats deployed
WAVeS on VBF as well** (novbs python numbers: vbf +0.03, zjets +5.33, dijet
+1.84, ttbar +1.97) -- positive on all four samples with no per-sample
configuration anywhere in the chain.

## Incidental ROOT footgun, recorded in the header comment

The score's long name originally said "q/P precision". Histogram names embed the
long name as the TFile key, and `TFile::Get` (and uproot's `file[key]`) treat
'/' in a key path as a DIRECTORY separator -- the write succeeds and every
readback fails with "key not found" (clustering_plot aborted rc=134). Renamed to
"qP precision"; do not put a slash back into any Score long name.

---

# Follow-up: the z0 and d0 siblings (`python/zd0_precision_terms.py`)

The same construction applied to the other two per-track variances. All three
precision sums are near-duplicates -- mutual correlations of log(1/sigma):
z0-d0 **+0.93 to +0.95**, d0-qP **+0.92 to +0.95**, z0-qP +0.81 to +0.88 --
so this was a test of which packaging of ONE signal is best, and whether a
second copy adds anything.

## Alone on TZ + guarded timing (worst-sample delta at each optimum)

| term | optimum | worst-sample | vbf | zjets | dijet | ttbar |
|---|---|---:|---:|---:|---:|---:|
| **(1/sigma_d0)^b** | **1.2-1.6 plateau** | **+0.29/+0.30** | +0.29 | +0.42 | +0.33 | +0.30 |
| (1/sigma_qP)^b | 0.8 | +0.22 | +0.26 | +0.24 | +0.25 | +0.21 |
| (1/sigma_z0)^b | 0.4 | +0.09 | +0.09 | +0.43 | +0.22 | +0.23 |

d0 is the best on every sample; z0 is zjets-loving but VBF-capped. The d0 curve
is a clean plateau at 1.2-1.6 turning over by 3.2.

## Stacked -- they are one signal, carried once

| stack | worst-sample |
|---|---:|
| z0 on top of qP@0.8 | +0.03 |
| d0 on top of qP@0.8 | +0.02 |
| z0 on top of d0@1.6 | **-0.04 / -0.09** |

Nothing survives stacking. Carry exactly one precision term.

## Consequence: the shipped term is now d0, not qP

`Score` id 33 renamed `TRKPTZ_TZP` ("[per-track #Deltaz + d0 precision]"),
computed as `(SUM_t 1/var_d0)^{D0_PRECISION_WEIGHT/2}` with
`D0_PRECISION_WEIGHT = 1.2` (plateau centre). d0@1.6 vs the superseded qP@0.8:
vbf +0.04, zjets +0.15, dijet +0.05, ttbar +0.09 -- uniformly >=.

C++ validation (local VBF, mjj>=200, all 48,314 events): TRKPTZ_TZP = 91.83%
(+1.47 over TRKPTZ, +0.12 over deployed WAVeS), d0 term worth +0.25 over
TZ_GIJ against the scan's +0.29 (novbs); WAVES reproduces at 91.71% exactly.
QP_PRECISION_WEIGHT retained in the header as a record but unused.
