# Why the Z+jets VBS pair almost never lands in R2

`util/vbs_region_diag.cxx`, full grid samples read via the skimmed copies
(`/data/mcardiff/skimmed_ntuples/`), window m_jj >= 500 GeV and |dEta| >= 2.5 —
the same window as `results/vbs_region_mjj_fourband.md`.

The question: Z+jets should be **dominated** by R2 — a forward pileup fake
paired with the central hard-scatter jet the Z recoils against — and the
four-band composition put only 3.46% of its pairs there, with 92.5% in "Other".
"Other" had never been decomposed, so three explanations were live: the
topology is absent; it is present but `calcBestVbsPair`'s max-m_jj ranking
outranks it; or the winning pair is not an analysis pair at all.

**All three are true, in that order of importance.**

## 1. The central hard-scatter leg usually does not exist

Event content, Z+jets, 36,588 events in the window:

| | Z+jets | dijet | VBF |
|---|---|---|---|
| >= 1 **central** paper-HS jet | **26.84%** | 93.08% | 69.58% |
| >= 1 forward (2.38–4.0) paper-PU jet | 88.17% | 43.17% | 28.35% |
| >= 1 forward paper-HS jet | 13.55% | 71.87% | 89.34% |
| **NO paper-HS jet anywhere** | **63.80%** | 1.62% | 0.04% |
| >= 1 paper-PU jet beyond \|eta\| 4 | 59.79% | 28.85% | 20.06% |
| mean jets > 30 GeV (of which \|eta\| > 4) | 5.51 (1.18) | 6.26 (0.58) | 4.47 (0.44) |

R2 requires a central paper-HS leg, and **only 26.84% of selected Z+jets events
have one** — 63.80% have no hard-scatter jet anywhere above 30 GeV. The Z's
visible energy goes into the two leptons; the recoil is frequently not a 30 GeV
jet inside a paper-HS cone. That caps R2 at ~27% before any pairing question is
asked, and it is the single largest term.

## 2. When the topology IS there, the m_jj picker usually discards it

| | R2 pair EXISTS | picker CHOSE it | missed |
|---|---|---|---|
| Z+jets | 16.22% | 3.46% | **79%** of those available |
| dijet | 33.90% | 9.29% | 73% |
| VBF | 14.98% | 3.53% | 76% |

m_jj grows like cosh(dEta), so the ranking systematically prefers a wide pair
over a genuine forward-fake + central-hard-scatter one. On Z+jets the pair that
wins instead is overwhelmingly pileup + pileup reaching past the acceptance —
`far-PU + fwd-PU` 25%, `far-PU + far-PU` 15%, `fwd-PU + fwd-PU` 15% (55%
pileup-pileup) — with its more-forward leg at |eta| median **4.17** (q75 4.31)
and median m_jj 1223, against the R2 pair's 393.

The R2 pair itself is exactly the expected object:

| Z+jets R2-shaped pair | \|eta\| (q25/med/q75) | pT (q25/med/q75) |
|---|---|---|
| forward PU leg | 2.79 / **3.09** / 3.70 | 34 / 40 / 48 GeV |
| central HS leg | 0.68 / **1.22** / 1.86 | 37 / 48 / 69 GeV |

## 3. The definitions discard it twice over

Every quantity computed three times over the same events, differing only in
what counts as a taggable jet and where "forward" ends:

| convention | | R1 | R2 | both-HS | Other | R1+R2 |
|---|---|---|---|---|---|---|
| `all_` every jet, fwd = 2.38–4.0 (status quo) | Z+jets | 2.39% | 3.46% | 1.63% | 92.51% | 5.85% |
| `acc_` \|eta\| < 4.0 jets only, fwd = 2.38–4.0 | Z+jets | 5.91% | **7.80%** | 2.75% | 83.55% | 13.70% |
| `wide_` every jet, fwd = anything > 2.38 | Z+jets | 5.94% | 5.38% | 1.63% | 87.05% | 11.32% |
| `all_` | dijet | 6.00% | 9.29% | 38.29% | 46.42% | 15.30% |
| `acc_` | dijet | 10.22% | **15.52%** | 47.33% | 26.92% | 25.74% |
| `wide_` | dijet | 14.64% | 15.31% | 38.29% | 31.76% | 29.96% |
| `all_` | VBF | 5.82% | 3.53% | 58.83% | 31.82% | 9.35% |
| `acc_` | VBF | 9.60% | 6.14% | 66.02% | 18.25% | 15.73% |
| `wide_` | VBF | 15.44% | 6.57% | 58.83% | 19.16% | 22.01% |

**Restricting the pair search to analysis-quality jets (|eta| < 4.0) more than
doubles Z+jets R2**, 3.46% -> 7.80%, and cuts "Other" from 92.5% to 83.6%. It
is the single most effective change, because a jet past |eta| 4 is neither
forward nor central under the status-quo zones and so can only ever land the
pair in Other — while being exactly the jet that maximises m_jj.

**Under `acc_`, Z+jets IS R2-dominated as expected**: R2 (7.80%) > R1 (5.91%),
the same ordering as dijet (15.52 vs 10.22) and the opposite of VBF (6.14 vs
9.60, signal-like). That ordering is invisible under the status quo only
because half the Z+jets pairs are won by |eta| > 4 pileup jets.

## Caveat on `wide_`

Extending "forward" with no upper edge more than doubles R1 (2.39 -> 5.94% on
Z+jets, 5.82 -> 15.44% on VBF), but a leg past |eta| 4 is outside HGTD and
cannot be timed. `wide_` R1 is a TOPOLOGY statement, not a "timing can act
here" statement — up to now those have been the same claim, and they should not
be conflated when quoting a reachable fraction.

## Validation

- `all_` reproduces the published four-band numbers: Z+jets R1+R2 5.85% and
  dijet 15.30% match `results/vbs_region_mjj_fourband.md` exactly. VBF reads
  9.35% against the published 9.51% because that figure was measured on the
  LOCAL 33-file VBF sample and this is the grid `highstats_vbf`.
- The local R1/R2 reconstruction is cross-checked against the shared
  `classifyVbsRegion` per event: **0 disagreements** on every event of every
  sample.
- Read via the skimmed copies, whose per-sample counters reproduce the full
  samples exactly — see `skim_n_*` in CLAUDE.md's "Skimmed + slimmed copies".

## Not answered here

Why the R_pT gain in R1/R2 is small on Z+jets is a separate question with a
separate answer: the truth-t0 row only reaches x1.33 in Z+jets R1 (against
x1.89 VBF, x1.74 dijet), so the CEILING is low rather than the algorithm being
at fault. Z+jets forward hard-scatter jets are track-poor (median R_pT 0.21 vs
0.40 on VBF), which forces the 80%-efficiency working point down to R_pT ~ 0.06
where the pileup distribution is dense.
