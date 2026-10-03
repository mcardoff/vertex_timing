# When NOT to use the time: the quality flag Q (2026-10-03)

Athena withholds its HGTD vertex time when a BDT score is below 0.3. This is the
equivalent for `TZP_KDE_TZ` (results/kde_mean_shift.md), with no training and
no truth: `Cluster::kdeQuality`, filled by `applyKernelDensityScore`.

    Q = (S1 - S2) / sqrt(S1 + S2)

S1 is the selected candidate's kernel-density score. S2 is the best score among
candidates whose time is more than 60 ps away, so the best COMPETING time (0 if
there is none). Q is the significance by which the chosen time beats the
runner-up. No time is provided below `KDE_QUALITY_MIN = 2`.

"Purity" below means the purity of the provided time: of the events where a
time is handed out, the fraction within 60 ps of truth. Mean cluster purity
(pT-weighted hard-scatter fraction of the selected cluster) is quoted beside it.

Scripts in `python/kde_study/`: `qual.py` (variables), `qual2.py`, `qual4.py`
(combinations), `qual3.py` (calibration, working points), `quality_curve.py`
(figure), `quality_report.py` (the before/after tables from the C++ outputs).

## Which variable

22 reco-only candidates ranked on the four novbs exports by how pure the best
90 / 70 / 50% of events are. AUC for separating right from wrong times:

| variable | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| **Q = (S1 - S2)/sqrt(S1 + S2)** | **0.758** | **0.897** | **0.922** | **0.900** |
| S1 - S2 | 0.749 | 0.893 | 0.919 | 0.895 |
| S1 (winning score alone) | 0.705 | 0.862 | 0.888 | 0.860 |
| margin 1 - S2/S1 | 0.703 | 0.833 | 0.856 | 0.836 |
| envelope of the picked cluster | 0.715 | 0.862 | 0.879 | 0.862 |
| kernel density / event total | 0.693 | 0.838 | 0.862 | 0.839 |
| n in-jet tracks at the mode | 0.572 | 0.841 | 0.897 | 0.817 |
| \|z_C - z_PV\| | 0.665 | 0.751 | 0.750 | 0.748 |
| kernel time RMS | 0.584 | 0.686 | 0.690 | 0.684 |
| time chi2/ndf of the cluster | 0.512 | 0.560 | 0.558 | 0.566 |
| n clusters, n tracks | 0.50-0.55 | | | |

The information is the GAP to the competing time, not any property of the
winner. Everything combined with the gap is flat: dividing by the kernel RMS,
multiplying by the density fraction, by an in/out-of-jet agreement penalty or by
exp(-|dz|) moves the AUC by under 0.01. In/out-of-jet agreement on its own is a
strong flag where defined (VBF: 98.1% right when the two agree within 30 ps,
71.9% beyond 60 ps) but is undefined for 46% of Z+jets events and adds nothing
on top of Q.

## Q means nearly the same thing on every sample

P(|dt| < 60 ps) in bins of Q, novbs exports:

| Q | zjets | dijet | vbf | ttbar | vbf mu=0 |
|---|---:|---:|---:|---:|---:|
| < 0.25 | 34 | 38 | 40 | 40 | 94 |
| 0.5-0.75 | 43 | 50 | 53 | 53 | 96 |
| 1.0-1.5 | 55 | 64 | 65 | 65 | 97 |
| 1.5-2.0 | 64 | 73 | 73 | 75 | 99 |
| 2.0-3.0 | 74 | 83 | 86 | 85 | 99.6 |
| 3.0-4.0 | 82 | 90 | 93 | 93 | 99.8 |
| 4.0-6.0 | 89 | 96 | 97 | 97 | 99.9 |
| > 6.0 | 96 | 99.4 | 99.7 | 99.6 | 100 |

dijet, vbf and ttbar agree within about 3 points in every bin, so one threshold
means one reliability on those. Z+jets sits 8-10 points lower at the same Q: a
tenth of its events have no timed hard-scatter track, and a clean winner among
pileup candidates looks the same to Q. What differs between samples is how many
events reach a given Q, which is the right behaviour for a quality flag.

## Working points against Athena (novbs exports)

Time provided / right among provided:

| | zjets | dijet | vbf | ttbar |
|---|---|---|---|---|
| no cut | 100 / 68.0 | 100 / 90.5 | 100 / 92.9 | 100 / 90.7 |
| Q >= 1 | 74.9 / 77.0 | 92.4 / 94.0 | 93.8 / 95.8 | 91.8 / 94.3 |
| **Q >= 2** | **54.2 / 83.8** | **85.2 / 96.2** | **88.0 / 97.6** | **83.6 / 96.7** |
| Q >= 3 | 37.2 / 88.3 | 77.5 / 97.5 | 81.6 / 98.5 | 74.9 / 98.0 |
| Athena (BDT >= 0.3) | 55.6 / 74.6 | 87.8 / 89.9 | 91.5 / 92.2 | 86.1 / 90.4 |
| ours at Athena's acceptance | 55.6 / 83.4 | 87.8 / 95.5 | 91.5 / 96.6 | 86.1 / 96.1 |

At the same fraction of events served, the provided time is wrong about half as
often as Athena's on dijet, vbf and ttbar (4-5% against 8-10%), and a third less
often on Z+jets (17% against 25%). Q >= 2 lands close to Athena's own
acceptance on every sample without being tuned to. Figure:
`figs/quality/quality_curve.png` (table in the `.csv` beside it).

## Before and after the cut, C++ on the grid

Clusters 3248510-17, 92 shards, none held; default selection m_jj >= 200.
`TZP_KDE_TZ_Q` (Score id 36) is `TZP_KDE_TZ`'s pick with the time withheld where
Q < 2. **A withheld time stays in the denominator and counts as a failure**,
exactly as an invalid Athena time does in the HGTD row, so the differential
core-fraction curves of the three rows share one denominator. (The first
version gated the denominator instead; clusters 3250027-30 are the rerun.)

| | | vbf | zjets | dijet | ttbar |
|---|---|---:|---:|---:|---:|
| events | | 630,720 | 64,950 | 105,137 | 863,552 |
| time provided | before | 100% | 100% | 100% | 100% |
| | after | 88.3% | 54.2% | 84.3% | 83.1% |
| core fraction of provided times | before | 93.03 | 67.78 | 89.94 | 90.29 |
| | after | **97.48** | **83.85** | **96.00** | **96.55** |
| mean cluster purity | TZP | 69.9 | 37.3 | 65.4 | 65.0 |
| | before | 72.5 | 40.0 | 68.4 | 68.1 |
| | after | **77.0** | **52.4** | **74.6** | **74.5** |
| core counted over ALL events (no time = fail) | before | 93.03 | 67.78 | 89.94 | 90.29 |
| | after | 86.07 | 45.42 | 80.95 | 80.19 |
| Athena, core over all events | | 84.53 | 41.60 | 77.88 | 77.24 |

The last block is the price: the cut removes 4.5-7 wrong times per 100 events
on vbf / dijet / ttbar and gives up 7-10 right ones (Z+jets: 23 wrong, 22
right). Counted that way the gated time still beats Athena's on every sample.

Differential plots: `<s>/comparisons/<s>_quality_<key>.pdf` (HGTD / TZP_KDE_TZ /
TZP_KDE_TZ_Q), bundled per sample in `figs/quality/<s>_quality.pdf`. All three
curves are core fraction over every selected event. Read that way the gated
curve sits below the ungated one everywhere (it gives up right times too) and
above Athena's: on VBF it tracks Athena's below ~8 forward HS tracks and
reaches 99% above 13 where Athena's plateaus near 95%; on Z+jets it is 2-4
points above Athena's in every populated bin. The purity of the provided time
is the table row above, not these curves.

Athena's row from the same files, for scale (provided / right among provided /
right over all): vbf 91.6 / 92.2 / 84.5, zjets 56.0 / 74.3 / 41.6, dijet 86.9 /
89.6 / 77.9, ttbar 85.7 / 90.1 / 77.2.

## RpT: the cut does not help the jet tagger

`rpt_v5_hist` scenarios 7-9: `kde` (time used in every event), `kde_q` (time
used only where Q >= 2, else the jet keeps its ITk-only R_pT), `kde_q1`
(Q >= 1). The flag passes in 77.4 / 50.0 / 71.5 / 71.5% of vbf / zjets / dijet /
ttbar events in the RpT selection (no VBS cut there, so lower than above).

Pileup rejection at MATCHED hard-scatter efficiency, interpolated on the traced
curve (`quality_report.py`), forward jets > 40 GeV:

| | | zonly | hgtd | tzp | kde | kde, Q >= 1 | kde, Q >= 2 |
|---|---|---:|---:|---:|---:|---:|---:|
| vbf | eff ceiling | 97.5 | 96.6 | 96.6 | 96.6 | 96.9 | 97.1 |
| | rej @ 0.80 | 121 | 154 | 156 | 157 | 149 | 144 |
| | rej @ 0.90 | 58.7 | 75.1 | 80.4 | 79.7 | 76.1 | 73.3 |
| zjets | eff ceiling | 92.8 | 91.0 | 90.8 | 90.9 | 91.5 | 91.9 |
| | rej @ 0.80 | 64.7 | 68.8 | 77.7 | 78.0 | 72.2 | 68.7 |
| | rej @ 0.90 | 12.4 | 9.6 | 11.3 | 11.3 | 11.4 | 11.8 |
| dijet | eff ceiling | 97.6 | 96.8 | 96.8 | 96.8 | 97.0 | 97.1 |
| | rej @ 0.80 | 71.1 | 82.8 | 84.6 | 85.3 | 81.9 | 79.0 |
| | rej @ 0.90 | 40.9 | 47.1 | 51.7 | 52.2 | 49.2 | 46.5 |
| ttbar | eff ceiling | 92.6 | 91.1 | 91.2 | 91.2 | 91.6 | 91.8 |
| | rej @ 0.80 | 62.4 | 79.2 | 87.8 | 89.1 | 81.4 | 76.6 |
| | rej @ 0.90 | 10.0 | 9.6 | 11.0 | 11.1 | 10.8 | 10.5 |

- **Without the cut, `kde` is level with `tzp`** (within 1-2% everywhere, R2 on
  VBF +3%). The 15% fewer wrong t0s do not show up in the tagger, which is the
  core-fraction / RpT dissociation already recorded in CLAUDE.MD.
- **The cut lowers the rejection at 0.80 by 7-14%** and Q >= 1 by about half of
  that, monotonically. A low-Q time is still right 40-75% of the time, and the
  3 sigma gate with a right t0 removes pileup tracks; withholding the time gives
  that up. With a wrong t0 the gate mostly empties the timed component of every
  jet, pileup ones included, so wrong times cost the tagger less than they cost
  the core fraction.
- **What the cut buys is efficiency reach.** The ceiling rises by 0.3-1.0
  points, and on the two samples whose ceiling is near 0.90 the rejection at
  0.90 is level or slightly up (zjets 11.3 -> 11.8, 30-40 GeV 2.2 -> 2.5).
- Regions: VBF R1 at 0.80 is 164 -> 157, R2 128 -> 111; ttbar R1 is unchanged
  (69.6 -> 70.5) and R2 56 -> 49. The Z+jets regions have 749 HS legs and move
  inside their noise.

So the flag is for consumers that need the TIME to be right (a vertex-time
measurement, a jet-vs-t0 veto, anything using t0 as a number). For the R_pT
gate, apply the time in every event.

## Calibration note

`kde`'s vertex-time inflation is seeded from tzp's except VBF (measured 1.33,
in use). Measured on this run: zjets 1.52 (in use 1.65), ttbar 1.38 (1.46),
vbf 1.32. A 17% change in inflation moves the gate width by 3%, so the numbers
above stand; update `inflationFor` before the next run.

## Not done

- A purity-targeted variant, Q x (mean e^{-|dz|} at the mode) / (kernel time
  RMS), gives +1 to +3 points of mean cluster purity at fixed acceptance and
  the same time purity. Not adopted: it adds two factors for a small gain.
- Using Q as a continuous gate width in R_pT (widen the time gate as Q falls)
  instead of an on/off switch. The table says on/off is the wrong shape there.
- The exporter and the event display do not carry Q yet.
