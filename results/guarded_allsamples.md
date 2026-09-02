# The guarded in-jet re-timing on all four mu=200 samples

`python/guarded_allsamples.py`, novbs selection. Selection and timing varied
independently, so every cell is comparable. Deployed WAVeS is the
(WAVeS, in-jet always) cell.

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| events | 711,901 | 93,977 | 132,138 | 1,049,983 |
| in-jet track fraction | 18.5% | **5.0%** | 15.8% | 13.2% |

## The grid

| selection | timing | vbf | zjets | dijet | ttbar |
|---|---|---:|---:|---:|---:|
| TRKPTZ | raw | 90.01 | 62.12 | 87.37 | 87.33 |
| TRKPTZ | in-jet always | 90.83 **+0.82** | 60.75 **-1.37** | 87.32 -0.04 | 87.21 -0.12 |
| TRKPTZ | **in-jet if >=3** | 90.49 +0.48 | 62.15 +0.02 | 87.59 +0.23 | 87.54 +0.21 |
| **TRKPTZ_TZ** | **raw** | **90.74** | **63.89** | **88.27** | **88.42** |
| TRKPTZ_TZ | in-jet always | 91.53 +0.79 | 62.46 **-1.43** | 88.20 -0.07 | 88.27 -0.15 |
| **TRKPTZ_TZ** | **in-jet if >=3** | **91.21** +0.47 | **63.90** +0.02 | **88.49** +0.23 | **88.63** +0.21 |
| WAVeS | raw | 90.62 | 60.21 | 86.91 | 86.98 |
| WAVeS | in-jet always *(deployed)* | 91.44 +0.82 | 58.81 **-1.40** | 86.90 -0.02 | 86.87 -0.11 |
| WAVeS | in-jet if >=3 | 91.09 +0.47 | 60.23 +0.02 | 87.13 +0.22 | 87.18 +0.21 |

## Three findings

**1. The >= 3 guard is safe on all four samples and positive on all four.** Worst
case +0.02. It is not a two-sample coincidence.

**2. Unconditional in-jet re-timing -- what deployed WAVeS does -- is a net
NEGATIVE on three of the four samples.** zjets -1.40, ttbar -0.11, dijet -0.02;
only VBF gains. The re-timing has been carried as an unqualified improvement and
it is not one outside VBF.

**3. The re-timing remains a pure additive term orthogonal to selection.** The
guard is worth +0.47/+0.48/+0.47 on VBF, +0.02 on zjets, +0.22/+0.23/+0.22 on
dijet and +0.21 everywhere on ttbar -- identical to two decimal places across
three different selectors. This is the fourth independent confirmation.

## TRKPTZ_TZ is the best selector on every sample

At matched (raw) timing: VBF 90.74 vs WAVeS 90.62 vs TRKPTZ 90.01; zjets 63.89
vs 60.21 vs 62.12; dijet 88.27 vs 86.91 vs 87.37; ttbar 88.42 vs 86.98 vs 87.33.
**WAVeS is the worst of the three selectors on zjets, dijet and ttbar**, and
wins only on VBF -- and even there only by 0.61 over TRKPTZ.

## The best sample-independent combination

`TRKPTZ_TZ` selection + in-jet time when the cluster holds >= 3 in-jet tracks,
against deployed WAVeS:

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| deployed WAVeS | **91.44** | 58.81 | 86.90 | 86.87 |
| TRKPTZ_TZ + guard | 91.21 | **63.90** | **88.49** | **88.63** |
| **difference** | **-0.23** | **+5.09** | **+1.59** | **+1.76** |

It loses 0.23 on VBF and wins 1.6-5.1 points on the other three, with **no
per-sample configuration anywhere**. Given Z+jets is the background to VBF
H->inv, +5.09 on the background against -0.23 on the signal is the trade that
matters.

## Caveats

- novbs selection throughout; the C++ numbers elsewhere use m_jj >= 200 and are
  not directly comparable (same direction and magnitude: C++ VBF gives WAVeS
  91.71 against TRKPTZ_TZ_GIJ 91.58, i.e. -0.13).
- Only the >= 3 working point was scanned on dijet/ttbar. >= 2 was not.
- Core fraction only. Nothing here says how the choice propagates to the RpT
  jet-rejection side, which is where the timing is actually consumed.
