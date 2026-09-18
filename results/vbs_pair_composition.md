# The VBS pair, fully decomposed

Every selected event's chosen tagging pair, classified by the **truth identity**
of each leg (paper HS / paper PU / neither) and by the **detector region** each
leg falls in, split at `|eta| = 2.4` with **no upper edge** — so a tagging jet at
`|eta| 4.3` is still forward; HGTD simply cannot time it.

Six leg types (2 regions x 3 identities) give 21 unordered pair compositions.
R1 and R2 are two of those cells; everything else is what the four-band plot
calls "Other".

    R1 = F-HS + F-PU        R2 = F-PU + C-HS

Source: `util/vbs_region_diag.cxx` `wide_` block (`figs/diagnostics/vbs_region_diag/*.root`),
re-binned offline at 2.4 by `figs/diagnostics/vbs_region_diag/decompose.py`.
Overlap removal ON (the analysis configuration). Interactive version:
https://claude.ai/code/artifact/798dc470-a78a-41a8-84aa-36ba12a9328c

## Headline

**Z+jets has no hard-scatter leg to pair with.** At `m_jj >= 200 GeV`, 79.7% of
Z+jets pairs carry no paper-HS jet on either leg, and the two pileup-only cells
hold 65.6% between them. R2's scarcity is not primarily the m_jj picker
outranking a real R2 pair — the central hard-scatter leg mostly is not there.

VBF inverts this on the same code and the same jets: 55.3% of its pairs put a
paper-HS jet on BOTH legs (53.9% of them with at least one leg forward), i.e.
the picker finds the real tagging pair. So the Z+jets result is a
statement about the sample, not the algorithm.

## Full table, m_jj >= 200 GeV

| pair composition | Z+jets % | n | VBF % | n |
|---|---:|---:|---:|---:|
| F-PU + C-PU        | 37.79 | 24,543 |  2.06 |  12,991 |
| F-PU + F-PU        | 27.84 | 18,079 |  3.55 |  22,367 |
| **F-PU + C-HS** (R2) | **7.80** | 5,064 | **6.91** | 43,613 |
| F-PU + C-X         |  5.08 |  3,298 |  0.28 |   1,739 |
| F-HS + C-PU        |  4.77 |  3,099 | 13.38 |  84,377 |
| F-PU + F-X         |  4.34 |  2,822 |  0.71 |   4,489 |
| **F-HS + F-PU** (R1) | **3.33** | 2,165 | **12.53** | 79,054 |
| F-X + C-PU         |  3.20 |  2,080 |  0.24 |   1,495 |
| F-HS + C-HS        |  1.68 |  1,094 | 37.30 | 235,276 |
| F-X + C-HS         |  0.78 |    509 |  0.88 |   5,542 |
| F-HS + C-X         |  0.74 |    478 |  1.94 |  12,253 |
| F-X + C-X          |  0.73 |    477 |  0.04 |     281 |
| C-PU + C-PU        |  0.49 |    318 |  0.01 |      87 |
| C-HS + C-PU        |  0.37 |    238 |  0.72 |   4,532 |
| F-HS + F-X         |  0.32 |    205 |  1.32 |   8,317 |
| F-X + F-X          |  0.25 |    165 |  0.04 |     270 |
| F-HS + F-HS        |  0.19 |    126 | 16.61 | 104,779 |
| C-PU + C-X         |  0.13 |     86 |  0.00 |      31 |
| C-HS + C-HS        |  0.09 |     59 |  1.37 |   8,640 |
| C-HS + C-X         |  0.06 |     36 |  0.09 |     587 |
| C-X + C-X          |  0.01 |      9 |  0.00 |       0 |
| **total** | | **64,950** | | **630,720** |

## Marginals

| | no HS leg | one HS | two HS | both fwd | fwd+cen | both cen | R1 | R2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| zjets, no m_jj cut | 79.7 | 18.3 |  2.0 | 33.8 | 64.7 | 1.5 |  3.11 | 7.80 |
| zjets, m_jj >= 500 | 83.4 | 14.9 |  1.6 | 63.4 | 36.4 | 0.1 |  5.81 | 5.46 |
| vbf,   no m_jj cut |  7.0 | 37.9 | 55.1 | 34.6 | 63.2 | 2.3 | 12.46 | 6.92 |
| vbf,   m_jj >= 500 |  6.8 | 34.4 | 58.8 | 42.2 | 56.4 | 1.4 | 15.20 | 6.69 |

## Where the forward leg sits (m_jj >= 200 GeV)

| cell | sample | % | fwd leg median \|eta\| | past \|eta\| 4.0 | median m_jj |
|---|---|---:|---:|---:|---:|
| R1  F-HS + F-PU | zjets |  3.33 | 3.44 | 34% | 1373 |
| R1  F-HS + F-PU | vbf   | 12.53 | 3.57 | 34% | 1996 |
| R2  F-PU + C-HS | zjets |  7.80 | 3.14 | **18%** |  436 |
| R2  F-PU + C-HS | vbf   |  6.91 | 3.71 | 39% |  804 |
| F-PU + F-PU     | zjets | 27.84 | 4.08 | **54%** | 1757 |
| F-PU + F-PU     | vbf   |  3.55 | 4.18 | **69%** | 2387 |
| F-PU + C-PU     | zjets | 37.79 | 3.23 | 20% |  404 |
| F-PU + C-PU     | vbf   |  2.06 | 3.92 | 46% |  501 |

Two things follow:

- **The pileup-only pairs are the far-forward ones.** F-PU + F-PU sits at median
  `|eta| 4.08` (zjets) with 54% of its forward legs past 4.0. m_jj grows like
  cosh(dEta), so the ranking reaches for whatever is widest, and what is widest
  is uninstrumented pileup. This is the (c) explanation from `vbs_region_diag`'s
  header, quantified.
- **R2 is the cell HGTD is best placed to act on** — its forward leg is the most
  central of the four (median 3.14 in zjets, only 18% past 4.0) — and it is
  simply rare.

## This is the overlap-removal result seen from the other end

In Z(ee)+jets, 47% of truth HS jets (pT > 10) sit within `dR < 0.1` of a truth
electron with pT > 20, and the fraction barely moves out to `dR < 0.5` (47.6%) —
i.e. these jets **are** the electrons, not jets near them. Overlap removal
(`LEPTON_JET_DR = 0.5`, `REMOVE_JETS`) then strips the reco jets that match them,
and what is left to tag on is pileup.

Rerunning the diagnostic with `--no-or` moves Z+jets R2 from 7.80% to ~20.7%
(acceptance zone) and roughly triples it in every zone convention, while R1
barely moves. **That is an upper bound, not a fix**: those jets are genuinely
lepton-contaminated and an analysis would remove them.

Open: the truth HS jet collection in the ntuple appears to be a bare
`AntiKt4TruthJets`-style collection with no lepton removal, so the paper-HS
*label* inherits the same contamination. If a dressed/lepton-removed collection
(`AntiKt4TruthDressedWZJets`) is available, every HS-labelled number above
should be re-measured against it.

## Caveats

- Paper HS and PU are **not complements** (HS: `dR < 0.3` of a truth HS jet
  pT > 10; PU: `dR > 0.6` from every truth HS jet pT > 4). The `X` columns are
  the gap between them, not an error.
- No clustering and no time gates run here — jet content only. Per-region
  rejection comes from `rpt_v5`, which keeps its own 2.4–3.8 forward window
  because it fills R_pT *for* each leg and an untimeable leg would dilute R1.
- Splitting at 2.4 rather than the clustering side's `MIN_ABS_ETA_JET = 2.38`
  moves everything by <= 0.25 pts: zjets R1+R2 at `m_jj >= 500` reads 11.28%
  here against 11.32% at 2.38.
- dijet was not included in this pass; its tree is present and
  `decompose.py` will pick it up by adding it to `SAMPLES`.

---

## The composition stack, redrawn

`util/scratch/vbs_stack_decomp.C` redraws `vbs_region_mjj`'s page-1 stack with
this labelling — same denominator, same `MJJ_EDGES`, same per-column
normalisation, same `|d-eta| >= 2.5` on the chosen pair — collapsing the 21
cells into seven bands ordered "timing can act on this" first:

    figs/diagnostics/vbs_region_diag/vbs_stack_decomp.pdf   (zjets, vbf)
                                     <sample>_vbs_stack_decomp.png

It reads the diag trees, not the ntuples, so it reruns in seconds.

### m_jj >= 500 GeV, |d-eta| > 2.5

| band | Z+jets | VBF |
|---|---:|---:|
| R1  fwd HS + fwd PU        |  5.81% | 15.21% |
| R2  fwd PU + central HS    |  5.46% |  6.69% |
| HS + PU, other eta config  |  2.31% |  8.76% |
| HS + HS (genuine pair)     |  1.63% | 58.83% |
| **PU + PU, both forward**  | **48.73%** | 4.29% |
| PU + PU, other eta config  | 21.96% |  1.25% |
| a leg carries neither label| 14.09% |  4.96% |
| plotted | 36,588 | 518,950 |

Events dropped below the axis floor: zjets 199 (`|d-eta|`) + 32,914 (m_jj);
vbf 1,794 + 113,821.

**The `|d-eta| > 2.5` cut is nearly redundant with `m_jj >= 500`.** Applied on
top it removes 4 further Z+jets pairs and 238 VBF pairs (36,592 -> 36,588 and
519,188 -> 518,950). With pT > 30 GeV legs, 500 GeV already forces a wide pair,
so the m_jj cut is doing essentially all the topology selection.

### What the stack shows

- **Z+jets is a PU+PU plot.** 70.7% of the selection is a pileup pair, and the
  both-forward half alone is 48.7% and *rises steeply* with m_jj: 12% in the
  500-750 column, 54% by 1 TeV, 80-85% above 2 TeV. The genuine-pair band is a
  1.6% sliver.
- **R2 dies with m_jj in BOTH samples; R1 grows in both.** R2 runs
  11% -> 0% (zjets) and 9% -> 0% (vbf) across the columns, while R1 runs
  2% -> 4% (zjets, peaking at 11% near 1.8 TeV) and 2% -> 36% (vbf). A wider
  pair is more likely to have both legs forward, which is R1's shape and not
  R2's — so the m_jj cut actively selects against the region HGTD is best
  placed to help. Worth weighing before raising the working point.
- **What differs between the samples is what fills the space R2 vacates**:
  genuine R1 in VBF, far-forward pileup pairs in Z+jets.
- **The 14.1% unlabelled band in Z+jets is not noise.** It is the paper labels'
  own gap (HS cone 0.3, PU cone 0.6), and at this size it is large enough that
  "Other" in the four-band plot was partly a labelling artefact.

---

## How much of the truth HS jet content is inside the lepton OR cone

Direct truth-level measurement (`~/skimcheck/truth_or.C` on the AF, full
skimmed Z+jets, 530,962 events, 100% of which carry >=1 selected lepton): a
truth HS jet counted as "inside the cone" if it lies within dR of a selected
lepton -- `Track_leptonID`-flagged AND pt > `LEPTON_MIN_PT` (20 GeV), i.e.
exactly `isGoodLepton`, with the truth jet substituted for the reco jet in
`computeOverlapRemoval`.

| truth HS jet | n | dR<0.2 | dR<0.3 | dR<0.4 | **dR<0.5** |
|---|---:|---:|---:|---:|---:|
| pt > 10        | 1,048,172 | 33.98% | 34.38% | 34.86% | **35.55%** |
| " \|eta\|<2.4  |   861,466 | 40.86% | 41.34% | 41.90% | **42.70%** |
| pt > 20        |   566,637 | 61.51% | 61.74% | 62.01% | **62.43%** |
| " \|eta\|<2.4  |   523,331 | 65.83% | 66.08% | 66.36% | **66.79%** |
| pt > 30        |   439,332 | 73.58% | 73.72% | 73.88% | **74.14%** |
| " \|eta\|<2.4  |   419,144 | 76.26% | 76.41% | 76.57% | **76.83%** |

`dR < 0.5` is `LEPTON_JET_DR`, the cone the analysis removes with.

**At the analysis jet floor (pt > 30 GeV) it is 74% -- 77% centrally.** The
number is nearly flat in cone size (73.6% at 0.2 vs 74.1% at 0.5), so these are
jets sitting essentially ON the lepton, not jets clipped by a generous cone.
It rises steeply with pT (35.6% -> 62.4% -> 74.1%) because the harder the truth
"jet", the more likely it simply IS one of the two ~45 GeV Z decay electrons.

Reco side, same events, same lepton definition, OR as actually applied:
**13.90%** of all reco jets and **26.67%** of reco jets with pt > 30 are removed
(486,257 / 3,497,108 and 395,591 / 1,483,020). The reco figure is far lower
because its denominator is dominated by pileup jets, which no lepton is near.

Relation to the earlier 47% figure: that measured truth HS jets (pt>10) within
dR<0.1 of a TRUTH ELECTRON (pt>20) over 20 raw unselected files. This measures
against RECONSTRUCTED lepton tracks over the Z-selected skim. Different lepton
definition and different event population, so the two are not the same
quantity; both say the same thing about the mechanism.
