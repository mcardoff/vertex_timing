# Integrating the aggregate t0 into TZP: measured, and it does not work

`python/t0_integration.py`. Both sample-independent entry points for the
classical aggregate t0, scanned on all four novbs samples against the
joint-tuned TZP baseline.

## Baselines

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| TZP pick | 91.66% | 64.71% | 88.97% | 89.08% |
| t0 alone | 91.24% | 62.89% | 87.81% | 87.87% |
| **union** | **94.12%** | **71.31%** | **91.64%** | **92.03%** |

The union headroom is real (2.5-6.6 points). Note TZP now beats the aggregate
standalone on every sample, VBF included -- the joint tune closed the one place
the aggregate used to win inclusively.

## A. As a score agreement term -- monotonically harmful

`S' = S * f(|t_C - t0|)`, Cauchy and exponential forms:

| tau [ps] | worst-sample (Cauchy) | worst-sample (exp) |
|---|---:|---:|
| 60 | -2.19 | -2.40 |
| 250 | -0.65 | -0.92 |
| 1000 | -0.01 (inert) | -0.10 |

Strictly worse the stronger it is; "safe" only when it does nothing. An earlier
zjets-only scan on the WEAKER TZ base had found tau=250 mildly positive (+0.27)
-- on the tuned TZP base the same term is -0.65. **The improved score already
carries what agreement-with-the-aggregate was proxying**, so the term now only
adds the aggregate's own mistakes.

## B. As a timing fallback -- worse than harmful

Report t0 instead of the cluster time when they disagree by more than X:

| X [ps] | fires (zjets) | worst-sample |
|---|---:|---:|
| 60 | 19% | -2.52 |
| 100 | 15% | **-2.81** |
| 250 | 5% | -0.76 |
| 400 | 1% | -0.04 (inert) |

Costs 2-3 points wherever it actually fires. The mechanism was already measured
from the ML side: **conditional on the two answers disagreeing, the aggregate is
the right one only ~27% of the time** (the ML meta's t0-wins fraction). A
wholesale switch on disagreement therefore replaces a right answer with a wrong
one ~3x more often than the reverse. Disagreement flags hard events; it does not
say to trust the aggregate.

## Standing conclusion

The aggregate's value survives as (1) a CONFIDENCE flag -- |t_pick - t0| still
separates 79%-core events from 58%-core events -- and (2) the union ceiling a
learned chooser could in principle chase (the ML meta captured ~+0.9 of it).
There is no simple global rule that converts it into core fraction on top of
TZP: two entry points, two functional forms each, all negative wherever active.
This closes the classical t0-integration question.
