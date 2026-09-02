# The oracle overstates the Z+jets headroom by ~2.5x (2026-08-28)

Prompted by a direct challenge to the claim that a large fraction of Z+jets is
helpable by HGTD. The challenge was right. Three independent tests, all on the
canonical Z+jets export (36,531 events, 209,451 clusters, 5.7 per event).

## The oracle is largely combinatorial

The "oracle" is best-|dt|-over-all-candidates, i.e. a BEST-OF-N statistic, and N
is 5.7. Null test: keep each event's cluster COUNT but draw its cluster times
from the pooled distribution over all events, destroying every bit of event
structure and leaving only "best of N draws".

    NULL oracle  53.4%
    REAL oracle  92.6%      excess over chance +39.2 points

**58% of the ceiling is reachable with no information at all.** Two corroborating
signatures:

- the oracle rate RISES with candidate count -- 72.2% at 1 cluster, 86.1% at
  2-3, 94.5% at 6-8, 96.8% at 9+. Physics does not improve with more pileup
  clusters; a best-of-N statistic does.
- **13.6%** of oracle-selected clusters have truth purity exactly 0. A pure
  pileup cluster at the right time is a coincidence, not a hard-scatter time
  a selector could learn to find.

## A ceiling that does not have this problem

Restrict the oracle to clusters that actually carry hard scatter:

| sample | TRKPTZ | naive oracle | oracle, purity>0.5 | oracle, any HS track |
|---|---|---|---|---|
| zjets | 61.9 | 92.6 | **51.2** | 83.9 |
| vbf   | 90.5 | 98.3 | 86.0 | 97.1 |
| dijet | 86.9 | 97.9 | 81.5 | 95.9 |
| ttbar | 87.1 | 98.0 | 81.9 | 96.5 |

On Z+jets the majority-HS ceiling (51.2) is BELOW what TRKPTZ already delivers
(61.9), and the same ordering holds on all four samples. That is only possible
if a large share of current performance already comes from clusters that are not
majority hard-scatter but land at the right time regardless.

## What is genuinely recoverable

Of the 11,225 Z+jets "recoverable failures" (30.7% of events):

    target cluster has truth purity exactly 0   22.1%   <- unlearnable
    target cluster has NO HS tracks at all      22.1%
    target cluster is majority hard-scatter     40.6%

    GENUINELY recoverable                       4,555 events = 12.5%

**The learnable headroom over TRKPTZ is ~12.5 points, not 30.7**, and the
realistic Z+jets ceiling is ~74.3%, not 92.6%. The remainder asks a selector to
pick a pileup cluster that happens to sit near the truth time, which nothing
reconstructable distinguishes -- unlearnable by construction, not merely hard.

## What this does to the ML result

Deep Sets reaches 67.1-67.3 against TRKPTZ's 61.9.

    old framing: +5.4 of 30.7 available   = 18% of the headroom, 25.3 left
    new framing: +5.4 of ~12.4 available  = 44% of the headroom, ~7.0 left

The work looks BETTER under the honest ceiling, not worse, and the residual gap
is a third of what was reported. Every "25.5-point gap" statement in the status
page and the failure gallery was measured against an inflated target.

## Methodological rule this establishes

**A best-of-N oracle must be null-tested before it is quoted as headroom.** The
cost of not doing it here: a target 2.5x too large, which made six successive
null results look like failures to close a gap that was mostly not there.
