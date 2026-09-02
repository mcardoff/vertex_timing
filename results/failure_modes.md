# Why the guarded score loses the events it loses — local VBF

`python/failure_modes.py`. 39,681 events; ours 91.98%, WAVeS 92.18%. Raw
disagreements 1,001 ours-only / 1,077 WAVeS-only; after the provenance screen
(a right answer existed, the two picked different clusters, the loser's cluster
is genuinely pileup-like) **543 clean WAVeS-wins and 488 clean ours-wins**, i.e.
**1.37% vs 1.23% of all events**. Half of each raw set survives the screen, so
the screen is symmetric.

## Where we lose: our score is maximised by the WRONG cluster

Medians over the 543 events, our (losing) pick against WAVeS's (winning) pick:

| quantity | **our** cluster | WAVeS' cluster |
|---|---:|---:|
| n tracks | **11** | 9 |
| sum pT [GeV] | **18.76** | 16.48 |
| **\|delta z\| [mm]** | **0.233** | 0.548 |
| n in-jet tracks | **0** | 2 |
| frac pT in fwd jet | **0.000** | 0.393 |
| min dR to fwd jet | 0.626 | 0.076 |
| truth purity | **0.000** | 0.594 |
| n truth HS tracks | **0** | 4 |

**The cluster we pick has more pT and sits closer to the primary vertex in z than
the real hard-scatter cluster.** Our score is
`SUM_t pT_t exp(-0.7|dz_t|) * exp(-1.5|dz_cl|)` -- pT and z only -- so on these
events it is not mis-tuned, it is *correctly maximised by a pileup cluster*.
Both of its inputs point the wrong way.

Failure modes of our pick (not exclusive):

| | share |
|---|---:|
| **closer to the PV in z than the winner's cluster** | **77.0%** |
| **winner's cluster has more in-jet tracks** | **79.7%** |
| **our cluster has NO in-jet tracks at all** | **71.6%** |
| higher sum pT than the winner's cluster | 54.5% |
| both bigger AND closer in z (fully degenerate on our observables) | 34.1% |

No amount of retuning `TRACK_DZ_WEIGHT` fixes a case where the wrong cluster is
closer in z. **The only discriminating information present is jet association,
and our score does not use it** -- 71.6% of our wrong picks contain no
jet-proximate track whatsoever.

## Where WAVeS loses: it is captured by a jet that has no hard scatter in it

Same table for the 488 events we win:

| quantity | WAVeS' (losing) cluster | **our** cluster |
|---|---:|---:|
| n tracks | 8 | 8 |
| sum pT [GeV] | 12.66 | **15.11** |
| \|delta z\| [mm] | 0.375 | 0.305 |
| n in-jet tracks | **1** | 0 |
| frac pT in fwd jet | **0.106** | 0.000 |
| truth purity | **0.000** | 0.635 |
| n truth HS tracks | **0** | 4 |

WAVeS picks a cluster with a jet-proximate track and **zero** hard-scatter
content, over ours which has no jet association at all and purity 0.635. Only
**7.2%** of its losses are cases where our cluster had more in-jet tracks -- so
this is not the mirror of our failure, it is the documented "WAVeS is captured by
the fake" mode: a pileup cluster sitting inside a forward jet.

## The two failure modes are complementary, not competing

We fail when z and pT both favour pileup and only jet association could save us.
WAVeS fails when jet association favours pileup and z/pT would have saved it.
Neither is a tuning problem in the other's direction.

That is an argument for combining them rather than choosing between them: a score
carrying pT, z AND jet association has, on this evidence, a different and larger
reachable set than either. It is also a warning -- WAVeS's jet term is exactly
what makes it fail on Z+jets (forward jets there are usually pileup), so the
combination has to weight jet association by something that knows when the jets
are trustworthy.

## Displays

`_000001` ev 798, ev 3096: bigger + closer in z + no jet content (129 of the 543).
`_000001` ev 1293: bigger in pT but further in z (111 of the 543) -- 20 tracks at
purity 0.00 beating a 4-track all-hard-scatter cluster.
