# TZP in the RpT study (rpt_v5, local VBF) — first consumer-side result

New scenario `tzp` (index 6, appended so the load-bearing 0-5 are untouched):
TZP-selected cluster + guarded in-jet time, gated with the standard
`GATE_SIGMA = 2.5` machinery. Inflation seeded from the waves value (1.38);
the calibration table measured **1.37** on the first run, so no iteration.

## Calibration row — TZP's picks are also the best-timed

| scenario | n core | core sig | tail>150 | ratio |
|---|---:|---:|---:|---:|
| hgtd | 542,652 | 41.7 ps | 13.9% | 1.48 |
| trkptz | 579,567 | 38.6 ps | 15.2% | 1.39 |
| waves | 579,389 | 38.5 ps | 15.2% | 1.38 |
| **tzp** | **583,339** | **38.0 ps** | **14.6%** | 1.37 |

Most core tracks, smallest core width, lowest tail of any real scenario.

## Central baseline (the standing correctness check): PASSES

The tzp row agrees with every other scenario in both central slices to within
the expected few-jet wiggle. No leak.

## Forward >40 GeV — TZP is the best real scenario at every WP to 0.95

| scenario | 0.85 | 0.90 | 0.93 | 0.95 | 0.97 |
|---|---:|---:|---:|---:|---:|
| zonly | 93.2 | 57.1 | 37.5 | 21.3 | 6.3 |
| hgtd | 126.8 | 71.2 | 51.1 | 16.2 | 8.8 |
| trkptz | 117.4 | 75.5 | 38.0 | 17.6 | 9.7 |
| waves | 117.4 | 78.2 | 46.3 | 25.4 | 9.9 |
| **tzp** | **126.8** | **83.4** | **51.5** | **26.9** | 9.8 |
| waves_ideal | 162.5 | 107.4 | 56.1 | 27.1 | 9.3 |
| truth | 176.1 | 115.2 | 62.1 | 29.1 | 10.6 |

At 0.90 HS efficiency: **83.4 vs waves' 78.2 (+7%) and trkptz's 75.5 (+10%)**.
At 0.95, tzp (26.9) is within 0.2 of the waves_ideal oracle (27.1). It closes
~18% of the waves -> waves_ideal gap at 0.90.

## Forward 30-40 GeV — parity with trkptz, waves ahead at loose WPs

tzp 101.9 / 27.4 / 10.9... vs waves 104.3 / 33.6 / 10.8. Consistent with the
core-fraction differentials: WAVeS's remaining edge is sparse/low-pT jets, and
the soft slice is where those live. rejAtEff quantization caveat applies (the
10.9 repeats are the sampler clamping at the efficiency ceiling).

## Caveats

- Local VBF only; the zjets/dijet RpT runs need condor (`rpt_v5_hist.sub`) and
  the tzp inflation seeds (zjets 1.65, dijet 1.46) checked against their
  measured ratios.
- `rejAtEff` steps: near a bin edge a 0.02% eff difference moves printed
  rejection by ~30% (the shared 126.8 at 0.85 with hgtd is same-bin clamping).
  Read the ordering across WPs, not single cells.
