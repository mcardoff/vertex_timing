# Selection: weighting tracks by how reliable their TIME is (2026-10-05)

> Numbers here are on the nominal track list. `wide_track_list.md` (same day)
> widens the list and takes the fails removed to 29.0 / 17.8 / 22.6 / 27.6%.

Follows `kde_mean_shift.md` (clustering at its ceiling, selection is what is
left) and `time_quality.md` (the flag Q). Scripts: `python/kde_study/sel.py`
(candidates, top-2 structure, aggregate), `sel2.py` / `sel9.py` (tie-break
scans), `sel3..7.py` (score terms), `sel8.py` (validation), `sel10.py` (Q
working points).

## Result

Two factors on the per-track weight of `TZP_KDE` / `TZP_KDE_TZ`:

    h   = 0.5 for a track with ONE HGTD hit, 1 otherwise        (selection density only)
    r_z = 1 / (1 + sigma_z0 / 2 mm)                             (selection density and time weight)

    w_sel  = sqrt(pT) e^{-|z0 - z_PV|} * r_z * h
    w_time = sqrt(pT) e^{-|z0 - z_PV|} * r_z / sigma_t^2

Grid C++ (clusters 3262766-73, 92 shards, none held; m_jj >= 200):

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| events | 630,720 | 64,950 | 105,137 | 863,552 |
| TZP | 91.69 | 64.48 | 88.30 | 88.50 |
| TZP_KDE_TZ before (2026-10-03) | 93.03 | 67.78 | 89.94 | 90.29 |
| **TZP_KDE_TZ with reliability** | **93.33** | **68.43** | **90.23** | **90.66** |
| TZP fails removed, before | 16.2% | 9.3% | 14.0% | 15.6% |
| **TZP fails removed, now** | **19.8%** | **11.1%** | **16.5%** | **18.8%** |

All four samples are now past the 10% bar. Local VBF C++: 3,877 -> 3,122 fails
(19.5%), offline prediction 18.5-20.6%.

## Mechanism: a wrong time is a random time

Fraction of truth-HS tracks whose HGTD time is within 3 sigma of the truth HS
time (novbs exports):

| | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| 1 HGTD hit (31% of tracks, sigma_t 35 ps) | 62.5 | 68.5 | 69.7 | 67.7 |
| 2 hits (33%, 25 ps) | 89.1 | 92.2 | 92.8 | 91.9 |
| 3 hits (26%, 25 ps) | 88.8 | 91.6 | 92.4 | 91.3 |
| 4 hits (10%, 20 ps) | 88.1 | 90.5 | 91.5 | 90.3 |

It is a step at one hit, flat above. Within each class it falls with the
track's extrapolation precision and rises with pT (VBF, 1 hit / >= 2 hits):

| sigma_z0 [mm] | < 0.3 | 0.3-0.6 | 0.6-1.2 | 1.2-2.4 |
|---|---|---|---|---|
| time right | 87 / 98 | 74 / 96 | 65 / 92 | 57 / 86 |

| pT [GeV] | 1-1.5 | 2-3 | 5-10 | 10-30 |
|---|---|---|---|---|
| time right | 60 / 87 | 71 / 94 | 84 / 97 | 89 / 98 |

It does not depend on event multiplicity. A third of single-hit hard-scatter
tracks carry a time unrelated to the vertex, and those tracks look perfect in z.
They are what the competing fake candidates are built from. The time weight
already halves a single-hit track through 1/sigma_t^2 (35^2 vs 25^2); the
selection density had no such term.

This is also why "HS purity by n hits" was flat in the earlier weight study
(`exp3.py`): hit count says nothing about where a track came from, only about
whether its time can be believed.

## How it was found

The top two candidates (winner and best distinct-time runner-up) contain the
right time far more often than the winner alone:

| | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| pick right | 67.96 | 90.54 | 92.93 | 90.70 |
| pick or runner-up right | 80.43 | 94.87 | 96.54 | 95.26 |
| pick or aggregate t0 right | 72.41 | 92.03 | 94.30 | 92.10 |
| any candidate right | 89.53 | 97.62 | 98.46 | 97.90 |

So the game is a tie-break between two candidates. On the pairs where exactly
one is right and Q < 1, the score itself is right 60-63% of the time. A scan of
32 per-candidate variables ("take the candidate with the larger X") found
nothing above that, but three weak ones that are not functions of the score:
mean time resolution of the tracks at the mode (58-59%), mean hit count
(56-58%) and |eta| inverted (55-57%). All three are the same thing, and as a
per-track weight it is worth more than as a tie-break. After adding it the
same scan shows mean time resolution at 53-54%: absorbed.

## Scan record (percent of TZP fails removed, zjets / dijet / vbf / ttbar)

| selection density weight | |
|---|---|
| reference (TZP_KDE_TZ, 2026-10-03) | 9.2 / 14.3 / 15.1 / 16.0 |
| x (25 ps / sigma_t)^1 | 9.9 / 15.3 / 16.6 / 17.0 |
| x (25 ps / sigma_t)^2 | 10.0 / 15.2 / 17.1 / 16.9 |
| x (25 ps / sigma_t)^4 | 9.0 / 12.4 / 13.9 / 13.6 |
| single-hit x 0.65 | 10.0 / 15.5 / 16.9 / 17.2 |
| **single-hit x 0.5** | **10.3 / 15.7 / 17.0 / 17.5** |
| single-hit x 0.35 | 10.5 / 15.7 / 16.9 / 17.6 |
| single-hit x 0.25 | 10.5 / 15.5 / 16.2 / 17.0 |
| factors (1, 2, 3, 4) by hit count | 9.7 / 14.5 / 15.5 / 15.7 |
| single-hit x 0.5 and 1/(1 + sigma_z0/1 mm) | 10.7 / 16.2 / 17.5 / 19.2 |
| **single-hit x 0.5 and 1/(1 + sigma_z0/2 mm)** | **10.9 / 16.2 / 18.2 / 19.0** |
| single-hit x 0.5 and 1/(1 + sigma_z0/4 mm) | 10.7 / 16.2 / 18.1 / 18.6 |
| + the sigma_z0 factor in the TIME weight too (adopted) | 11.1 / 16.7 / 18.5 / 19.6 |
| + the hit factor in the time weight | 10.4 / 15.4 / 17.3 / 18.7 |
| + reliability in the clustering seeds | 10.9 / 16.1 / 18.3 / 18.9 |
| + hit factor inside the d0 information sum | 10.4 / 15.6 / 17.2 / 17.7 |

Graded factors by hit count lose to the step, as the mechanism table says they
should.

## Out of sample, mu = 0, by track count

| set | events | TZP | before | now | fails removed |
|---|---:|---:|---:|---:|---:|
| vbf novbs, files 120-299 | 87,187 | 91.65 | 93.09 | 93.38 | 20.6% |
| ttbar novbs, files 60-179 | 130,330 | 89.11 | 90.85 | 91.18 | 19.0% |
| zjets mjj500 | 36,531 | 64.31 | 67.81 | 68.47 | 11.7% |
| dijet mjj500 | 80,383 | 88.58 | 90.19 | 90.49 | 16.7% |
| vbf mjj500, files 380-499 | 42,548 | 91.96 | 93.28 | 93.57 | 20.0% |
| ttbar mjj500, files 300-399 | 69,305 | 88.60 | 90.40 | 90.77 | 19.1% |

mu = 0, passing events (TZP / before / now): vbf 57,592 / 57,606 / 57,599,
ttbar 42,692 / 42,697 / 42,694, zeejets 34,227 / 34,228 / 34,226. Neutral.

Core fraction by number of timed HS tracks, TZP -> now:

| n HS timed | zjets | dijet | vbf | ttbar |
|---|---|---|---|---|
| 0 | 19.5 -> 19.0 | 20.4 -> 20.7 | 23.9 -> 21.4 | 20.1 -> 19.9 |
| 1 | 34.7 -> 38.3 | 37.9 -> 41.6 | 43.6 -> 49.5 | 38.7 -> 44.7 |
| 2-3 | 51.5 -> 57.1 | 60.8 -> 66.2 | 71.3 -> 75.4 | 62.5 -> 68.2 |
| 4-7 | 73.3 -> 78.1 | 85.4 -> 88.5 | 90.3 -> 92.6 | 85.3 -> 88.7 |
| >= 8 | 90.8 -> 94.0 | 97.4 -> 98.3 | 98.1 -> 98.8 | 97.1 -> 98.2 |

## The quality flag after the change

The reliability factors scale every score down, and Q with it, so the threshold
moves from 2.0 to **1.5** (`KDE_QUALITY_MIN`), which keeps the acceptance.
Grid C++:

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| time provided at Q >= 1.5 | 89.7% | 56.9% | 85.7% | 84.8% |
| right among provided | 97.41 | 83.44 | 95.88 | 96.45 |
| core over all events (withheld = fail) | 87.36 | 47.49 | 82.20 | 81.80 |
| same, before (Q >= 2, 2026-10-03) | 86.07 | 45.42 | 80.95 | 80.19 |
| Athena, core over all events | 84.53 | 41.60 | 77.88 | 77.24 |
| mean cluster purity, no cut -> cut | 72.6 -> 76.7 | 40.2 -> 51.9 | 68.4 -> 74.3 | 68.2 -> 74.1 |

At Athena's own acceptance the provided time is right 96.9 / 83.9 / 95.8 /
96.4% of the time (vbf / zjets / dijet / ttbar; was 96.6 / 83.4 / 95.5 / 96.1).

## RpT

Forward jets > 40 GeV, rejection at matched 0.80 efficiency:

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| ITk-only | 121 | 64.7 | 71.1 | 62.4 |
| Athena time | 154 | 68.8 | 82.8 | 79.2 |
| tzp | 156 | 77.7 | 84.6 | 87.8 |
| kde before / now | 157 / 156 | 78.0 / 80.2 | 85.3 / 85.9 | 89.1 / 89.6 |
| kde, Q >= 1.5 | 146 | 69.4 | 79.4 | 77.3 |

`kde` is now level with or above `tzp` on every sample (zjets +3%). The quality
cut still costs rejection, as before.

## Closed in this round

- **Low-Q fallback to the aggregate t0.** The aggregate is worse than the pick
  in every Q bin, including the lowest (zjets Q < 0.5: pick 36.5%, aggregate
  29.7%). Third time this has been measured, third negative.
- **Low-Q switch to the runner-up.** The runner-up is right 29-35% of the time
  in the lowest Q bin against 36-45% for the pick.
- **A beam-time prior on the candidate** (`S x exp(-k t^2 / 2 sigma^2)`): +0.4
  on zjets at k = 0.1, negative on the others by k = 0.25.
- **|eta| weights**: positive on zjets, strongly negative on VBF.
- **Peakedness, sideband contrast, sigma-scaled kernels, kernel chi2, z
  coherence of the tracks at the mode**: at or below coin-flip as tie-breaks.
- In-jet count is a strong tie-break on VBF (69% on low-Q pairs) and useless on
  Z+jets (52%), the known sample dependence.

## What is left

Pick-or-runner-up is still 80.9% on zjets against a pick of 68.6% (VBF 96.6 vs
93.2). On low-Q decisive pairs the score is right 60-63% and no reco variable
tested says which of the two it is.
