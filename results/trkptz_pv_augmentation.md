# Augmenting TRKPTZ with vertex-proximity information — local VBF + Z+jets, 2026-09-01

Four related proposals tested against the TRKPTZ baseline. All are implemented
in the main program as registry scores (`src/clustering_constants.h`, ids 24–27),
computed in `Cluster::updateScores` beside TRKPTZ and sharing its identical
`exp(-1.5|dz_cluster|)` factor, so ONLY the pT sum differs. None appears in
`Cluster::calculateTime`'s score list, so all report the standard
inverse-variance mean and every comparison is **pure selection**.

## Headline

| variant | VBF | Z+jets |
|---|---:|---:|
| WAVES (reference) | +1.08 | — |
| **TRKPTZ_TZ** — per-track `exp(-0.7·\|z0_t − z_PV\|)` | **+0.64** | **+1.76** |
| TRKPTZ_PUW — PU-side pT down-weighted, w = 0.4 | +0.28 | +1.04 |
| TRKPTZ_PV — PU-side pT vetoed outright | −1.27 | — |
| TRKPTZ_PU — orientation control | −9.03 | — |
| TRKPTZ baseline | 93.04% | 62.12% |

**The per-track Δz weighting is the winner** and is the one to carry forward.
VBF and Z+jets peak at the **same exponent a = 0.7** from very different
baselines, which is the main reason to believe it rather than a per-sample fit;
the curve is broad, with anything in 0.4–1.2 within 0.1 of the peak on both.

Details of each below.

---

## 1. The `closer_to_pu_than_pv` hard veto — costs 1.27

Proposal: `score = exp(-1.5|dz_cl|) · SUM_t { pT if <flag matches>, else 0 }`,
with `RecoVtx_z[0]` as the primary. Flag helper
`BranchPointerWrapper::closerToPuThanPv`, mirroring `vertexRelation()` in
`util/export_training_data.cxx` — reco-only, no truth.

The orientation control settles that this is not a sign error: keeping the
PU-side tracks instead is **9 points worse**, so the flag is informative and
`_PV` reads it the right way round. It is simply too blunt (local VBF, 15.25M
timed tracks):

| quantity | value |
|---|---:|
| P(flag = 1 \| genuine HS track) | **46.8%** |
| P(flag = 1 \| PU track) | 88.0% |
| clusters vetoed to exactly zero | **56.0%** |
| events where every cluster scores zero (argmax arbitrary) | 1.02% |

It removes 88% of pileup but **also 47% of the hard scatter**. TRKPTZ ranks
clusters by hard-scatter pT content, so destroying half that signal costs more
than the pileup suppression buys.

## 2. Softening the veto to a down-weight — gains 0.28

`SUM_PV pT + w·SUM_PU pT`, so w = 1 is TRKPTZ and w = 0 is the veto. Scanned on
exported VBF (TRKPTZ reproduced from `sumpt`/`delta_z` to max relative error
9.6e-05 over 2.93M clusters, and `sumpt` from the track table to 1.2e-07,
before trusting any of it):

| w | 0.0 | 0.2 | 0.3 | **0.4** | 0.5 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| core frac | 89.22 | 90.65 | 90.74 | **90.76** | 90.73 | 90.70 | 90.58 | 90.45 |

The python w = 0 point (−1.23) reproduces the C++ `TRKPTZ_PV` row (−1.27) across
a different selection and a separate implementation.

## 3. Per-track Δz weighting — gains 0.64 / 1.76, and supersedes 2

`score = exp(-1.5|dz_cluster|) · SUM_t pT_t · exp(-a·|z0_t − z_PV|)`. The
cluster-level term is **kept, not replaced**: replacing it is 0.1–0.2 worse.

| a | 0.1 | 0.3 | 0.5 | **0.7** | 0.9 | 1.2 | 2.0 |
|---|---|---|---|---|---|---|---|
| Z+jets | +0.64 | +1.40 | +1.66 | **+1.76** | +1.73 | +1.67 | +1.22 |
| VBF | +0.33 | +0.59 | +0.67 | **+0.68** | +0.63 | +0.54 | +0.21 |

Two functional-form findings, both consistent with prior work:

- **Raw \|Δz\| in mm beats the z0 significance** (best +1.76 vs +1.21 on Z+jets).
  Same conclusion the z-association study reached about `sqrt(var_z0)`
  mis-modelling the true z0 spread.
- **`exp` beats `1/Δz`**, which is worse than the baseline at every floor tried.
  A `1/Δz` weight diverges and over-rewards a track that happens to sit on the
  vertex.

This contradicts the standing comment on `Score::WAVES` in
`src/clustering_structs.h` ("the cluster-level z-term is more effective than
per-track z-pull weighting because it averages over track-by-track z noise").
That holds for a *steep* per-track weight — `exp(-6·dz)` is −1.77 — but a gentle
one applied ON TOP of the cluster term is worth more than either alone.

**Not complementary with 2.** TZ +1.76 and PUW +1.04 combine to **+1.58**, i.e.
worse than TZ alone. `closer_to_pu_than_pv` is a coarse binary version of the
same track-to-PV distance that `exp(-0.7·dz)` encodes continuously. Use TZ, drop
PUW.

## 4. Weighting by the canonical vertex-WAVeS score — refuted

Proposal: `w_t = SUM_V WAVeS(V) / dz(V, t)`, so a track close in z to a vertex
that looks hard-scatter-like is boosted — the appeal being that
`computeVertexScores` uses `nearestAnyJet` (ALL jets, not just forward), so it
would inject **central** information into an otherwise forward-only selector.

Sound motivation, but it does not survive contact with the numbers.

**The canonical WAVeS score is overwhelmingly concentrated on one vertex, and
that vertex is the PV.**

| | VBF | Z+jets |
|---|---:|---:|
| canonical WAVeS picks the PV | 85.9% | **96.4%** |
| canonical sumpt2 picks the PV | 65.0% | 77.3% |
| median n reco vertices | 103 | 106 |

For clusters whose nearest vertex is NOT the WAVeS-best one, `waves(nearest) /
waves(best)` has median **9.9e-05** on Z+jets (VBF 9.0e-04); only 1.29% exceed
0.10. So on Z+jets the sum is dominated by a single term:

    w_t  ~=  WAVeS(V*) / dz(V*, t)      with V* = the PV, 96.4% of the time

and `WAVeS(V*)` is a **per-event constant**. Cluster selection is an argmax
*within* an event, so that constant cancels exactly and the score reduces to
`SUM_t pT / dz(t, PV)` — which is proposal 3, with the worse `1/dz` functional
form. **The central information enters as an event-level constant and therefore
carries no ranking information at all.**

Why Z+jets is the extreme case: the lepton term `pT^4 / 0.01` dominates
everything (a 45 GeV muon contributes 4.1e8), so WAVeS there is close to "which
vertex has the Z leptons" — a near-binary signal. VBF, with no leptons and a
jet-driven score, is far less concentrated (9.1% of non-best-nearest clusters
above 0.10), so this is the one place the multi-vertex term could have distinct
content — but VBF has little headroom to begin with.

Measured directly on top of TZ (Z+jets), using the exported cluster-level
`nearest_vtx_waves_frac`:

| multiplier on TZ | frac^0.05 | frac^0.1 | frac^0.25 | frac^0.5 | frac^1.0 |
|---|---|---|---|---|---|
| vs TZ | −0.54 | −1.47 | −3.94 | −5.13 | −5.69 |

Monotonically harmful, and a binary "halve any cluster whose nearest vertex is
not WAVeS-best" is −0.98. The failure mode is the same one as proposal 1:
`nearest_vtx_waves_frac` is effectively a step function (~1 near the PV, ~1e-4
elsewhere), so multiplying by it is a near-veto that destroys more hard scatter
than it removes pileup.

## The lesson common to 1 and 4

Both failures are steep, near-binary z-based penalties, and both lose for the
same reason: at this pileup level the hard-scatter cluster is often NOT the one
nearest the primary vertex in z, so anything that sharply penalises distance
throws away real signal. The gentle continuous weight of proposal 3 is worth
+1.76 where the sharp versions of the same information are worth −1.27 and −5.69.

## Standing caveats

- **VBF and Z+jets only.** `TRACK_DZ_WEIGHT = 0.7` has not been scanned on dijet
  or ttbar. Do that before treating it as general.
- **TRKPTZ_TZ's +0.64 on VBF is still below WAVeS's +1.08**, though on Z+jets
  its +1.76 has no WAVeS comparison in this table yet — worth running.
- Z+jets numbers here come from the exported training data (`training_new/`,
  the m_jj≥500 + novbs union, 93,977 events), not from a `clustering_hist` run;
  the local ntuples are VBF only. The C++ `TRKPTZ_TZ` row was confirmed against
  VBF (+0.64 vs the python scan's +0.68 on a different selection).

---

## TRKPTZ_TZ in the main program — differential core fraction (local VBF)

`clustering_hist` + `clustering_plot`, 45,335 events. `makeComparisonPlots` now
emits a `trkptz_tz_<key>.pdf` group (TRKPTZ vs TRKPTZ_TZ vs WAVeS) for each of
the five KEYs; the two TRKPTZ curves share clusters, cluster-level envelope and
time estimator, so their separation is exactly the per-track term.

| | core fraction |
|---|---:|
| TRKPTZ | 93.04% |
| **TRKPTZ_TZ** | **93.68%  (+0.64, 283 events)** |
| WAVES | 94.11% |

**vs n forward HS tracks** — event-level binning, so this IS a paired comparison
(per-bin totals differ by at most 2 events):

| n fwd HS trk | events | TRKPTZ | TRKPTZ_TZ | change |
|---|---:|---:|---:|---:|
| 3 | 2319 | 73.74 | 74.69 | +0.95 |
| **4** | 3003 | 80.55 | 82.91 | **+2.35** |
| 5 | 3604 | 85.99 | 87.29 | +1.30 |
| 6 | 3906 | 89.40 | 90.52 | +1.12 |
| 7 | 4112 | 92.34 | 93.48 | +1.14 |
| 8 | 4094 | 94.36 | 95.24 | +0.88 |
| 9 | 3823 | 95.68 | 95.87 | +0.18 |
| 10+ | — | — | — | ~0 |

**The entire gain is at low multiplicity and it is gone by ~10 tracks.** That is
exactly the population `failure_decomposition` names as dominant (selection
failure is 52.9% of failures at 0-1 forward HS tracks, and ~78% of all events
that fail do so below 2 tracks). Above 10 tracks the baseline is already >96%
and there is nothing left to win.

**vs n forward jets** — the gain is largest where there is least jet information
to exploit, which is also where most events are:

| n fwd jets | events | TRKPTZ | TRKPTZ_TZ | change |
|---|---:|---:|---:|---:|
| 1 | 30896 (68%) | 88.86 | 89.77 | +0.91 |
| 2 | 13062 | 92.82 | 93.25 | +0.43 |
| 3 | 3162 | 94.05 | 94.43 | +0.38 |

**vs cluster PU fraction** — read with care: this axis is a property of the
SELECTED cluster, so the binning is score-dependent and up to 312 events migrate
between bins. It is NOT a paired comparison, unlike the two above.

| PU frac | TRKPTZ | TRKPTZ_TZ | change |
|---|---:|---:|---:|
| 0.0-0.1 | 97.82 | 95.74 | -2.08 |
| 0.5-0.6 | 97.08 | 96.58 | -0.49 |
| 0.8-0.9 | 58.65 | 60.53 | +1.89 |
| **0.9-1.0** | 26.57 | **31.07** | **+4.50** |

Suggestive rather than conclusive: the large gains sit in the dirtiest selected
clusters and there is an apparent loss in the cleanest bin, but some of both is
population migration. Worth a paired re-measurement before being quoted.
