#!/usr/bin/env python3
"""Derive the R_pT thresholds for vbs_region_diag's JVT working points, and
report what each working point does to the jet population.

Reads the `jets` tree of a `--jvt=none` run (every pT-passing jet BEFORE any
cut, with its paper HS/PU truth labels and discriminants) and finds the R_pT
cut whose paper-HS efficiency in the JVT window matches each target. The
targets are the published JVT working-point efficiencies: the default Run-2
"Medium" point (92%) for loose, "Tight" (85%) for tight. fJVT thresholds are
the published ones (0.5 / 0.4) and are only REPORTED here, not derived.

    python/vbs_jvt_calibrate.py figs/local_vbs_region_diag.root

Paste the printed R_pT thresholds into JVT_WPS in util/vbs_region_diag.cxx.
"""
import sys
import numpy as np
import uproot

TARGET_EFF = {"loose": 0.92, "tight": 0.85}
FJVT_MAX   = {"loose": 0.50, "tight": 0.40}

fn = sys.argv[1] if len(sys.argv) > 1 else "figs/local_vbs_region_diag.root"
t  = uproot.open(fn)["jets"]
a  = t.arrays(["pt", "abseta", "hs", "pu", "rpt", "corrjvf", "fjvt",
               "n_ghost", "jvt_win", "fjvt_win"], library="np")
hs, pu = a["hs"] > 0.5, a["pu"] > 0.5
print(f"{fn}: {len(hs):,} pT-passing jets, {hs.sum():,} paper-HS, {pu.sum():,} paper-PU, "
      f"{(~hs & ~pu).sum():,} neither")

def eff(mask_pass, mask_pop):
    n = mask_pop.sum()
    return (mask_pass & mask_pop).sum() / n if n else float("nan")

# ── JVT window: |eta| < 2.5, pT < 60 ───────────────────────────────────────
w = a["jvt_win"] > 0.5
print(f"\nJVT window (|eta|<2.5, pT<60): {w.sum():,} jets, HS {(w&hs).sum():,}, PU {(w&pu).sum():,}")
print(f"  jets with no ghost tracks: {(w & (a['n_ghost']==0)).sum():,} "
      f"(HS {(w&hs&(a['n_ghost']==0)).sum():,}, PU {(w&pu&(a['n_ghost']==0)).sum():,})")
rpt_hs = np.sort(a["rpt"][w & hs])
print(f"  R_pT of HS jets: median {np.median(rpt_hs):.3f}, "
      f"5th pct {np.quantile(rpt_hs, .05):.3f}, frac == 0: {(rpt_hs==0).mean():.3%}")
print(f"  R_pT of PU jets: median {np.median(a['rpt'][w&pu]):.3f}, "
      f"frac == 0: {(a['rpt'][w&pu]==0).mean():.3%}")
print(f"\n  {'WP':6s} {'target':>7s} {'R_pT >=':>8s} {'HS eff':>7s} {'PU eff':>7s} {'PU rej':>7s}")
thresholds = {}
for wp, target in TARGET_EFF.items():
    # keep the fraction `target` of HS jets: threshold at the (1-target) quantile
    thr = float(np.quantile(rpt_hs, 1.0 - target))
    keep = a["rpt"] >= thr
    e_hs, e_pu = eff(keep, w & hs), eff(keep, w & pu)
    thresholds[wp] = thr
    print(f"  {wp:6s} {target:7.2f} {thr:8.4f} {e_hs:7.3%} {e_pu:7.3%} {1/e_pu if e_pu else float('inf'):7.1f}")

# ── fJVT window: 2.5 <= |eta| < 4.5, pT < 120 ───────────────────────────────
f = a["fjvt_win"] > 0.5
print(f"\nfJVT window (2.5<=|eta|<4.5, pT<120): {f.sum():,} jets, HS {(f&hs).sum():,}, PU {(f&pu).sum():,}")
print(f"  fJVT of HS jets: median {np.median(a['fjvt'][f&hs]):.3f}, 90th pct {np.quantile(a['fjvt'][f&hs], .9):.3f}")
print(f"  fJVT of PU jets: median {np.median(a['fjvt'][f&pu]):.3f}, 10th pct {np.quantile(a['fjvt'][f&pu], .1):.3f}")
print(f"\n  {'WP':6s} {'fJVT <=':>8s} {'HS eff':>7s} {'PU eff':>7s} {'PU rej':>7s}")
for wp, fmax in FJVT_MAX.items():
    keep = a["fjvt"] <= fmax
    e_hs, e_pu = eff(keep, f & hs), eff(keep, f & pu)
    print(f"  {wp:6s} {fmax:8.2f} {e_hs:7.3%} {e_pu:7.3%} {1/e_pu if e_pu else float('inf'):7.1f}")

# direction check: a correct fJVT gives PU jets HIGHER values than HS jets
rng = np.random.default_rng(0)
x_hs, x_pu = a["fjvt"][f & hs], a["fjvt"][f & pu]
i, j = rng.integers(0, len(x_hs), 200000), rng.integers(0, len(x_pu), 200000)
print(f"\n  P(fJVT_PU > fJVT_HS) over random pairs: {(x_pu[j] > x_hs[i]).mean():.3f}  (should be well above 0.5)")
x_hs, x_pu = a["rpt"][w & hs], a["rpt"][w & pu]
i, j = rng.integers(0, len(x_hs), 200000), rng.integers(0, len(x_pu), 200000)
print(f"  P(R_pT_HS > R_pT_PU)  over random pairs: {(x_hs[i] > x_pu[j]).mean():.3f}  (should be well above 0.5)")

print("\nJVT_WPS thresholds to paste:")
for wp, thr in thresholds.items():
    print(f'  {{"{wp}", "_jvt{wp.capitalize()}", {thr:.4f}, {FJVT_MAX[wp]:.1f}}},')
