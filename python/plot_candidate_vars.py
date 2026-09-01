#!/usr/bin/env python3
"""
plot_candidate_vars.py -- population-wide distributions of variables NOT in the
current score, for the HS-dominant cluster vs every other cluster, on VBF and
Z+jets together.

One panel per variable, four curves: {VBF, Z+jets} x {HS-dominant, other},
all normalised to unit area. A variable can only be worth developing if the
solid and dashed curves separate in the SAME direction on BOTH samples -- the
per-sample population AUC is printed in each panel (folded so >0.5 always means
"separates", with the arrow giving the direction).

Population-wide by construction: every cluster of every multi-cluster event
enters, with no conditioning on any method winning or losing.
"""
import argparse, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tiebreak_scan import load, auc, EVT

# (column, pretty label, log-x?)  -- none of these enter the current score
VARS = [
    ("dt_cluster_to_hgtd", r"$|t_{clus}-t_{HGTD\,vtx}|$ [ps]  (valid only)", True),
    ("dt_to_agg_cl",       r"$|t_{clus}-t_0^{agg}|$ [ps]  (classical aggregate)", True),
    ("cluster_time_sigma", r"$\sigma_t$ of cluster [ps]", True),
    ("cluster_z_sigma",    r"$\sigma_z$ of cluster [mm]", True),
    ("cluster_d0_sigma",   r"$\sigma_{d0}$ of cluster [mm]", True),
    ("cluster_qOverP_sigma", r"$\sigma_{q/P}$ of cluster", True),
    ("frac_nhgtd_ge2",     r"frac. tracks with $\geq$2 HGTD hits", False),
    ("waves_score",        r"WAVeS score  (jet info, reference)", True),
    ("frac_pt_in_fwdjet",  r"frac. $p_T$ in fwd jet  (REJECTED, contrast)", False),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vbf", default="figs/hists/vbf_mjj500p0_training.root")
    ap.add_argument("--zjets",
                    default="/Users/mcard/project/vertex_timing/training_new/zjets_novbs_training.root")
    ap.add_argument("--out", default="figs/candidate_vars.pdf")
    a = ap.parse_args()

    D = {}
    for tag, path in (("VBF", a.vbf), ("Z+jets", a.zjets)):
        print(f"  loading {tag} ...", flush=True)
        c, _ = load(path)
        multi = c.groupby(EVT, sort=False)["cluster_idx"].transform("size") > 1
        D[tag] = c[multi]

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    col = {"VBF": "#1f77b4", "Z+jets": "#d62728"}
    for ax, (v, label, logx) in zip(axes.ravel(), VARS):
        # shared binning from the pooled 1-99% range
        pool = np.concatenate([D[t][v].to_numpy(float) for t in D])
        pool = pool[np.isfinite(pool)]
        lo, hi = np.percentile(pool, [1, 99])
        if logx:
            pos = pool[pool > 0]
            lo = np.percentile(pos, 1) if len(pos) else 1e-3
            hi = max(hi, lo * 10.0)          # degenerate range guard
            bins = np.geomspace(lo, hi, 41)
            ax.set_xscale("log")
        else:
            bins = np.linspace(lo, hi, 41)
        txt = []
        for tag in ("VBF", "Z+jets"):
            c = D[tag]
            x = c[v].to_numpy(float); dom = c["is_dom"].to_numpy()
            m = np.isfinite(x)
            pa = auc(x, dom)
            dirn = np.sign(pa - 0.5) if np.isfinite(pa) else 1.0
            fold = 0.5 + abs(pa - 0.5)
            arrow = r"$\uparrow$HS" if dirn > 0 else r"$\downarrow$HS"
            txt.append(f"{tag}: AUC {fold:.3f} {arrow}")
            for dm, ls, lw in ((True, "-", 1.8), (False, "--", 1.2)):
                sel = m & (dom == dm)
                h, e = np.histogram(np.clip(x[sel], bins[0], bins[-1]),
                                    bins=bins, density=True)
                ax.stairs(h, e, color=col[tag], linestyle=ls, linewidth=lw,
                          label=f"{tag} {'HS-dom' if dm else 'other'}")
        ax.set_yscale("log")
        ax.set_xlabel(label, fontsize=9)
        ax.text(0.03, 0.97, "\n".join(txt), transform=ax.transAxes, fontsize=8,
                va="top", ha="left",
                bbox=dict(fc="white", alpha=0.75, ec="none"))
        ax.tick_params(labelsize=8)
    axes[0, 0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Variables NOT in the current score: HS-dominant cluster (solid) "
                 "vs other clusters (dashed), whole population", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out)
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
