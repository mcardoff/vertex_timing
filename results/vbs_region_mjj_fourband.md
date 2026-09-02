# VBS pair composition vs m_jj -- four-band scheme, m_jj >= 500 GeV

`util/vbs_region_mjj.cxx` at e8025ab: regions classified from the two VBS
candidate-pair legs (reverted from the event-level definition at user request),
axis restricted to the analysis working point. Bands: R1, R2, both-VBS-jets-
truth-HS (any eta), Other (old R3 + both-PU + out-of-acceptance + neither-label).
R1/R2 cross-checked per event against classifyVbsRegion: agrees on every event,
all three samples.

| band | VBF | Z+jets | dijet |
|---|---:|---:|---:|
| plotted events (pair, mjj>=500) | 39,687 | 36,592 | 80,402 |
| R1: both tags fwd, HS + PU | 5.92% | 2.39% | 5.65% |
| R2: fwd PU + central HS | 3.59% | **3.46%** | 8.75% |
| Both VBS jets truth HS | **59.82%** | 1.64% | 41.86% |
| Other | 30.68% | **92.50%** | 43.75% |
| **HGTD can contribute (R1+R2)** | **9.50%** | **5.85%** | **14.40%** |

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
- **dijet**: the largest reachable fraction (14.4%), R2-dominated like zjets.

Definitional note: this is a deliberately STRICTER statement than the old
event-level CAN_HELP figure (which asked whether the EVENT contains anything
timeable and put ~60% of zjets in reach). Here the question is whether the
SELECTED PAIR presents a topology where the gate has a defined job. The two
numbers answer different questions; this one is the per-pair statement the
tagging actually acts on.
