#!/usr/bin/env python3
"""Hard-scatter jet-time quality against the number of timed tracks: why the jet tests
improve as the 'in HGTD acceptance' requirement tightens.

For every paper-HS tagging leg with a jet time in the signal file (all selected events, as
the jet-time calibration), (t_jet - t_HS^truth) is binned in
  left   n tracks in the jet time (the highest-sum-pT time cluster within dR < 0.2)
  right  n timed tracks in HGTD (the acceptance count, hgtd_ntimed)
and three quantities drawn per bin: the core resolution 1.4826 x median |dt| (ps), the core
fraction |dt| < 60 ps (the project's PASS_SIGMA), and the fraction passing the test's own
3 sigma against the true t0 (sigma x f(n)). The share of HS jets per bin is printed above.

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_jet_time_vs_ntrk.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root --out figs/time_veto/jet_time_vs_ntrk
"""
import argparse, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
plt.rcParams.update({"font.size": 13})

cols = ["t0_truth"] + [f"{L}_{k}" for L in "ab" for k in ["hs", "core_n", "core_t", "core_sig", "hgtd_ntimed"]]
e = uproot.open(args.sig)["events"].arrays(cols, library="np")
cal = C.measure(e, True, "signal HS legs")
dt, n, nh, pull = [], [], [], []
for L in "ab":
    m = e[L + "_hs"] & (e[L + "_core_n"] > 0)
    d = e[L + "_core_t"][m] - e["t0_truth"][m]
    dt.append(d); n.append(e[L + "_core_n"][m]); nh.append(e[L + "_hgtd_ntimed"][m])
    pull.append(np.abs(d) / (cal(e[L + "_core_n"][m]) * e[L + "_core_sig"][m]))
dt, n, nh, pull = map(np.concatenate, (dt, n, nh, pull))

BINS_CORE = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 8), (9, 10**6)]
BINS_HGTD = [(1, 2), (3, 4), (5, 6), (7, 9), (10, 14), (15, 19), (20, 10**6)]


def lab(lo, hi):
    return f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 10**6 else f"$\\geq${lo}")


fig, axes = plt.subplots(2, 2, figsize=(13, 7.4), sharex="col", gridspec_kw=dict(height_ratios=[1, 1.1], hspace=0.08, wspace=0.28))
for col, (x, bins, xl) in enumerate([(n, BINS_CORE, "tracks in the jet time (time cluster, $\\Delta R < 0.2$)"),
                                     (nh, BINS_HGTD, "timed tracks in HGTD (the acceptance count)")]):
    res, cf, pas, share = [], [], [], []
    for lo, hi in bins:
        m = (x >= lo) & (x <= hi)
        res.append(C.mad(dt[m])); cf.append(np.mean(np.abs(dt[m]) < 60)); pas.append(np.mean(pull[m] < 3))
        share.append(m.mean())
    xs = np.arange(len(bins))
    a0, a1 = axes[0, col], axes[1, col]
    a0.plot(xs, res, "o-", color="#3B6FB6", lw=2, ms=7)
    a0.set_ylabel("core resolution [ps]\n1.4826 median $|t_{\\rm jet} - t_{\\rm HS}|$")
    a0.set_ylim(0, max(res) * 1.35)
    for xi, (r, sh) in enumerate(zip(res, share)):
        a0.text(xi, r + max(res) * 0.06, f"{r:.0f}", ha="center", fontsize=11, color="#3B6FB6")
        a0.text(xi, max(res) * 1.27, f"{100 * sh:.0f}%", ha="center", va="top", fontsize=10, color="#777777")
    a0.set_title("grey: share of HS jets in each bin", fontsize=10, color="#777777", loc="right")
    a1.plot(xs, cf, "s-", color="#2E9E4F", lw=2, ms=7, label="core fraction: $|t_{\\rm jet} - t_{\\rm HS}| < 60$ ps")
    a1.plot(xs, pas, "^-", color="#C0392B", lw=2, ms=7, label="passes the test vs the true $t_0$ ($< 3\\sigma$)")
    a1.set_ylim(0.4, 1.02); a1.set_ylabel("fraction of HS jets")
    a1.set_xticks(xs); a1.set_xticklabels([lab(*b) for b in bins]); a1.set_xlabel(xl)
    a1.grid(alpha=0.3); a0.grid(alpha=0.3)
    if col == 0:
        a1.legend(loc="lower right", fontsize=11, frameon=False)
fig.suptitle("VBF $H\\to$inv., truth-matched hard-scatter tagging jets with a jet time "
             f"({len(dt):,} jets)   ATLAS Simulation Internal", fontsize=12, x=0.08, ha="left")
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(f"{args.out}.{ext}", dpi=150, bbox_inches="tight")
print("wrote", args.out + ".pdf / .png")
for b, r, c, p_ in zip(BINS_CORE, *[[C.mad(dt[(n >= lo) & (n <= hi)]) for lo, hi in BINS_CORE],
                                     [np.mean(np.abs(dt[(n >= lo) & (n <= hi)]) < 60) for lo, hi in BINS_CORE],
                                     [np.mean(pull[(n >= lo) & (n <= hi)] < 3) for lo, hi in BINS_CORE]]):
    print(f"  core n {lab(*b):>6}: {r:5.1f} ps  core frac {100 * c:5.1f}%  pass {100 * p_:5.1f}%")
