# Making the in-jet re-timing sample-independent

## The constraint

Z+jets is the background to VBF H->inv. A per-SAMPLE configuration therefore
gives signal and background different algorithms, and any resulting separation is
manufactured rather than physical. Earlier notes in
`results/injet_decomposition.md` proposed exactly that ("behind a
sample-appropriate gate", with the in-jet track fraction estimated per sample).
**That proposal is withdrawn.**

The distinction that matters: **a per-EVENT function of reconstructed quantities
is legitimate even when its distribution differs by sample** -- the algorithm is
fixed, the events differ. What is not allowed is the configuration depending on
which sample it is applied to.

## One rule, same threshold everywhere

TRKPTZ_TZ selection throughout; only the reported time varies.

| | VBF | Z+jets | fires on |
|---|---:|---:|---|
| raw (full cluster) | 91.49% | 63.89% | -- |
| in-jet ALWAYS *(not implementable)* | 92.31% (+0.82) | 62.46% (**-1.43**) | 100% |
| **in-jet if >= 3 in-jet tracks in the cluster** | **91.98% (+0.50)** | **63.90% (+0.02)** | 70% / 9% |
| in-jet if >= 2 | 92.28% (+0.80) | 63.81% (-0.08) | 81% / 16% |
| in-jet if in-jet timing-weight frac >= 0.2 | 91.99% (+0.50) | 63.89% (+0.00) | 73% / 12% |

The rule never inspects the sample. It fires on 70% of VBF events and 9% of
Z+jets ones purely because Z+jets clusters rarely contain three in-jet tracks --
**the sample dependence emerges from the data rather than being configured in.**
At >= 3 it captures 61% of the unattainable sample-dependent gain at no Z+jets
cost. At >= 2 it captures 98% of it for a -0.08 Z+jets cost, which is a judgement
call rather than a measurement.

Full scan: `python/unified_timesource.py`.

## Implementation

`Score::minSubsetTracks` (default **0**) -- the minimum tracks an IN_JET/OUT_JET
subset must contain before `calculateTime` uses its time; below it, the
full-cluster time is returned. 0 reproduces the historical "any surviving track
will do" behaviour, which the WAVES family keeps, so **WAVES is unchanged at
94.11% exactly**. `MIN_INJET_TRACKS_FOR_TIME = 3`.

`Score::TRKPTZ_TZ_GIJ` (id 30) is TRKPTZ_TZ's selection with the guarded in-jet
time. Local VBF, `clustering_hist` (m_jj >= 200):

| | core frac | vs TRKPTZ |
|---|---:|---:|
| TRKPTZ_TZ_IJ, unguarded | 94.35% | +1.32 |
| **TRKPTZ_TZ_GIJ, guarded** | **94.18%** | **+1.15** |
| WAVES (deployed) | 94.11% | +1.08 |
| TRKPTZ_TZ, raw t | 93.68% | +0.64 |
| TRKPTZ | 93.04% | -- |

The guard costs 0.17 on VBF and buys sample-independence. **94.18% still beats
deployed WAVeS (94.11%)** -- with a single algorithm that needs no knowledge of
which sample it is running on.

Cross-check: the export-based scan (m_jj >= 500) put the guarded rule +0.50 over
TRKPTZ_TZ; the C++ run (m_jj >= 200) gives 94.18 - 93.68 = +0.50. Two selections,
two implementations, same answer.

## What else in this work needs the same audit

- `TRACK_DZ_WEIGHT = 0.7` -- **clean**. Scanned independently on VBF and Z+jets,
  both peak at 0.7, so one value serves both.
- `PU_SIDE_PT_WEIGHT = 0.4` -- scanned on VBF only, and superseded by TRKPTZ_TZ
  anyway. Should not be used without a Z+jets scan.
- The `|t_in - t_out| < 60 ps` confidence flag -- **clean**, already a per-event
  reco function.
- Anything keyed on the in-jet TRACK FRACTION of a sample -- **withdrawn**.

Open: dijet and ttbar have not been checked against the >= 3 threshold. A rule
that is safe on two samples is not yet demonstrated safe on four.
