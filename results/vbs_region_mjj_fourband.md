# VBS pair composition vs m_jj -- four-band scheme, m_jj >= 500 GeV, |dEta| >= 2.5

`util/vbs_region_mjj.cxx` at af3b390: regions classified from the two VBS
candidate-pair legs (reverted from the event-level definition at user request),
axis restricted to the analysis working point. Bands: R1, R2, both-VBS-jets-
truth-HS (any eta), Other (old R3 + both-PU + out-of-acceptance + neither-label).
R1/R2 cross-checked per event against classifyVbsRegion: agrees on every event,
all three samples.

The pair |dEta| >= 2.5 requirement and the 500 GeV floor are **counted
exclusions inside the executable**, deliberately not routed through the
`VBS_JET_*` globals -- `classifyVbsRegion` consults those internally, so moving
the cut there would have silently changed what R1/R2 mean and broken the
cross-check that validates the local classification.

## Composition (VBF local 33-file; Z+jets and dijet full grid samples)

| band | VBF | Z+jets | dijet |
|---|---:|---:|---:|
| plotted events (pair, mjj>=500, dEta>=2.5) | 39,666 | 36,588 | 75,676 |
| R1: both tags fwd, HS + PU | 5.92% | 2.39% | 6.00% |
| R2: fwd PU + central HS | 3.59% | **3.46%** | 9.29% |
| Both VBS jets truth HS | **59.80%** | 1.63% | 38.29% |
| Other | 30.69% | **92.51%** | 46.42% |
| **HGTD can contribute (R1+R2)** | **9.51%** | **5.85%** | **15.30%** |

Readings:

- **VBF**: above 500 GeV the selected pair is the genuine VBS pair 60% of the
  time; the timing-actionable topologies are 9.5%.
- **Z+jets**: the pair is almost never a clean truth configuration -- 92.5%
  lands in Other (mixed pairs in other eta configs, out-of-acceptance legs, and
  neither-label jets; the paper HS/PU cones are not complements). Only 1.6% of
  background pairs are both-truth-HS, and R2 > R1 (3.46 vs 2.39): when timing
  does have a defined job on the background pair, it is more often rejecting a
  forward fake against a central genuine jet than picking between two forward
  candidates.
- **dijet**: the largest reachable fraction (15.3%), R2-dominated like zjets.

## What the |dEta| >= 2.5 cut changed

Same runs without the pair |dEta| requirement (commit e8025ab), for comparison:

| band | VBF | Z+jets | dijet |
|---|---:|---:|---:|
| R1 | 5.92 -> **5.92** | 2.39 -> **2.39** | 5.65 -> **6.00** |
| R2 | 3.59 -> **3.59** | 3.46 -> **3.46** | 8.75 -> **9.29** |
| Both truth HS | 59.82 -> **59.80** | 1.64 -> **1.63** | 41.86 -> **38.29** |
| Other | 30.68 -> **30.69** | 92.50 -> **92.51** | 43.75 -> **46.42** |
| R1+R2 | 9.50 -> **9.51** | 5.85 -> **5.85** | 14.40 -> **15.30** |

**The cut is nearly inert on VBF and Z+jets, and only mildly active on dijet** --
because m_jj >= 500 already implies a large rapidity gap for the jet pT spectra
these samples produce. Events excluded by |dEta| < 2.5, as a fraction of those
passing the pair requirement:

| | VBF | Z+jets | dijet |
|---|---:|---:|---:|
| pairs before the two cuts | 48,640 | 69,701 | 106,926 |
| excluded by \|dEta\| < 2.5 | 136 (0.3%) | 199 (0.3%) | 6,338 (5.9%) |
| excluded by m_jj < 500 | 8,838 (18.2%) | 32,914 (47.2%) | 24,912 (23.3%) |

**m_jj carries the topology selection essentially by itself on VBF and Z+jets**
(a 60x and 165x larger exclusion than |dEta| respectively), which is independent
confirmation of the 2026-08-26 finding that swapped the default selection from
|dEta| to m_jj. Only dijet has a meaningful population of high-mass narrow-gap
pairs -- 5.9% -- and removing them shifts the mix toward the HGTD-relevant
regions (R1+R2 14.40 -> 15.30) at the expense of clean two-HS-jet pairs
(41.86 -> 38.29). Those narrow-gap dijet pairs are disproportionately genuine
both-HS configurations, which is the expected direction: a wide-gap forward pair
in a QCD dijet sample is more likely to have picked up a pileup leg.

Definitional note: this is a deliberately STRICTER statement than the old
event-level CAN_HELP figure (which asked whether the EVENT contains anything
timeable and put ~60% of zjets in reach). Here the question is whether the
SELECTED PAIR presents a topology where the gate has a defined job. The two
numbers answer different questions; this one is the per-pair statement the
tagging actually acts on.

Sample-size note: Z+jets and dijet are the full grid samples (2.0M and 699.5k
events read); VBF is the local 33-file sample (112,400 read), so its
percentages carry more statistical noise than the other two columns.
