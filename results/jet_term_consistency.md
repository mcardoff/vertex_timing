# The bounded jet term: derived from a subpopulation, rejected on the population

## What motivated it

`missing_information.py`, run on the **543 events TRKPTZ_TZ loses to WAVeS**,
found that our own inputs carry nothing there -- cluster `|delta_z|` separates
the wrong pick from the true cluster at **AUC 0.502** and `sumpt` at **0.466** --
while every jet-association variable separates at **0.86-0.90** against a 0.978
truth ceiling. That looked decisive.

## The population check says otherwise

**1. Population AUC** -- identify the HS-dominant cluster among ALL clusters in
multi-cluster events, not just the losses:

| variable | VBF (3,989,385 clusters) | Z+jets (525,600) |
|---|---:|---:|
| `base` (our score) | **0.974** | **0.882** |
| `sumpt` | 0.971 | 0.871 |
| `n_tracks` | 0.937 | 0.853 |
| `\|delta_z\|` | 0.904 | 0.805 |
| **`frac_pt_in_fwdjet`** | **0.882** | **0.612** |

The jet fraction is the **weakest** variable available on both samples, and on
Z+jets it is barely above chance. Its 0.90 on the loss subset was a property of
that subset.

**2. Switch accounting** -- of the events whose pick actually changes, how many
improve?

| alpha | VBF switched | better:worse | net | Z+jets switched | better:worse | net |
|---|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.32% | **2.54** | +0.12 | 0.35% | 1.03 | +0.00 |
| 1.00 | 0.97% | **2.20** | +0.32 | 1.28% | **0.97** | -0.01 |
| 2.00 | 1.51% | 1.94 | +0.41 | 2.35% | 0.94 | -0.05 |
| 4.00 | 2.14% | 1.70 | +0.46 | 3.98% | **0.85** | -0.22 |

**On Z+jets the term is a coin flip.** At alpha = 1 it changes 1.28% of picks and
gets 417 right against 430 wrong -- 49.2%. Its earlier "safe on every sample"
verdict was safe only because it does *nothing* there, not because it helps.

**3. Differential.** On VBF the gain is real and sits where it should (+1.08 at
2-3 HS tracks, +0.76 at 4-5, decaying to 0 by 9 and slightly negative above 13).
On Z+jets every bin is within +-0.22 with no consistent sign.

## Verdict: not worth carrying

- It is **signal-only**: +0.32 on VBF, nil on Z+jets. Z+jets is the background to
  VBF H->inv, so a signal-only gain is the less valuable half of the trade.
- Even on VBF it is the weakest input available (0.882 against the score's own
  0.974), so the +0.32 is residual correlation rather than new power.
- It is the same information WAVeS uses, in a weaker form.

`JET_FRAC_WEIGHT` defaults to **0** and `Score::TRKPTZ_TZJ` is retained as a
reproducible diagnostic, not a recommendation.

## The transferable lesson

**A variable ranked on the events a method loses will overstate itself.** The
loss set is by construction the set where the incumbent score failed, so anything
correlated with the winner looks strong there. The AUC gap between the two
framings is large -- 0.90 on the subset against 0.612 population-wide on
Z+jets -- and only the population number predicts what happens when the term is
actually switched on. Rank on the subset to generate hypotheses; test on the
population before believing them.
