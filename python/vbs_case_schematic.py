#!/usr/bin/env python3
"""Schematic R-z views of the two cases the timing study splits the VBS tagging pair into,
for the talk's section dividers (python/vbs_case_schematic.py --out figs/time_veto/schematic).

  forward_forward   both tagging jets reach HGTD (2.38 < |eta| < 4.0): each has a time,
                    so they can be tested against each other AND against the vertex t0
  central_forward   one tagging jet is central (no HGTD time), one reaches HGTD: only the
                    forward jet can be tested, and only against t0

Both jets are drawn from the hard-scatter vertex: JVT/fJVT have already found them
compatible with it in z, which is why timing is the handle left. Geometry is approximate
(HGTD at |z| = 3.5 m, 120-640 mm in radius; ITk to |z| = 3 m, R = 1 m), not to scale in the
vertex region.
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--out", required=True, help="output stem: <out>_{forward_forward,central_forward}.(pdf|png)")
args = ap.parse_args()
plt.rcParams.update({"font.size": 13})

Z_HGTD, R_IN, R_OUT, DZ = 3.5, 0.12, 0.64, 0.12   # m
GREEN, HGTD_C, GRID = "#2E9E4F", "#3B6FB6", "#B9C3CF"


def r_at(z, eta):
    return abs(z) / np.sinh(abs(eta))


def jet(ax, eta, length, label, timed):
    """a jet cone from the origin, R = 0.4 wide in eta, to `length` m along its axis"""
    th = 2 * np.arctan(np.exp(-abs(eta)))
    sgn = np.sign(eta)
    half = 0.4 * np.sin(th) / 2 * 1.6   # widened a little so it reads at slide size
    pts = [(0, 0)]
    for d in (-half, half):
        pts.append((sgn * length * np.cos(th + d), length * np.sin(th + d)))
    ax.add_patch(Polygon(pts, closed=True, facecolor=GREEN, alpha=0.55, edgecolor=GREEN, lw=1.5, zorder=4))
    zm, rm = sgn * 0.62 * length * np.cos(th), 0.62 * length * np.sin(th)   # the cone's middle
    tz, tr = (sgn * 2.1, 1.2) if abs(eta) > 2 else (2.75, 1.05)
    ax.annotate(label, (zm, rm), xytext=(tz, tr), ha="center", fontsize=12, color=GREEN, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=1.2), zorder=7,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    if timed:   # where it crosses HGTD
        rc = r_at(Z_HGTD, eta)
        ax.plot([sgn * Z_HGTD], [rc], marker="*", ms=17, color="#E69F00", mec="black", mew=0.8, zorder=6)


def detector(ax):
    ax.add_patch(Rectangle((-3.0, 0), 6.0, 1.0, facecolor="#F2F4F7", edgecolor=GRID, lw=1, zorder=0))
    ax.text(0, 0.86, "ITk", ha="center", fontsize=12, color="#6B7785")
    ax.add_patch(Rectangle((-3.2, 1.4), 6.4, 0.6, facecolor="#EFE6D8", edgecolor="#C9B79A", lw=1, zorder=0))
    ax.text(-1.7, 1.62, "calorimeter (barrel)", ha="center", fontsize=11, color="#8A7454")
    for sgn in (-1, 1):
        ax.add_patch(Rectangle((min(sgn * 3.65, sgn * 4.4), 0.1), 0.75, 1.9, facecolor="#EFE6D8",
                               edgecolor="#C9B79A", lw=1, zorder=0))
        ax.add_patch(Rectangle((sgn * Z_HGTD - DZ / 2, R_IN), DZ, R_OUT - R_IN, facecolor=HGTD_C,
                               edgecolor="black", lw=0.8, zorder=3))
        ax.text(sgn * Z_HGTD, R_OUT + 0.1, "HGTD", ha="center", fontsize=12, color=HGTD_C, fontweight="bold")
        for eta, ls in ((2.38, "--"), (4.0, ":")):
            zz = 4.6
            ax.plot([0, sgn * zz], [0, r_at(zz, eta)], ls=ls, color="#7C8794", lw=1, zorder=1)
        ax.text(sgn * 4.62, r_at(4.6, 2.38), "$|\\eta| = 2.38$", fontsize=10, color="#7C8794",
                ha="left" if sgn > 0 else "right", va="center")
        ax.text(sgn * 4.62, r_at(4.6, 4.0), "$|\\eta| = 4.0$", fontsize=10, color="#7C8794",
                ha="left" if sgn > 0 else "right", va="center")
    ax.plot([-4.7, 4.7], [0, 0], color="black", lw=1.2, zorder=2)
    # the hard-scatter vertex and a few pileup vertices along the beam (spread exaggerated)
    for z in (-0.35, -0.2, 0.15, 0.3, 0.42):
        ax.plot([z], [0], "o", ms=5, color="#8C8C8C", zorder=5)
    ax.plot([0], [0], "o", ms=9, color=GREEN, mec="black", zorder=6)
    ax.text(0, -0.13, "hard-scatter vertex ($z$, $t_0$)", ha="center", va="top", fontsize=11)
    ax.set_xlim(-5.3, 5.3); ax.set_ylim(-0.3, 2.2); ax.set_aspect("equal")
    ax.set_xlabel("$z$ [m]"); ax.set_ylabel("$R$ [m]")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


for key, jets, note in [
    ("forward_forward", [(3.1, 3.6, "tagging jet 1: $t_{\\rm jet,1}$", True), (-3.4, 3.6, "tagging jet 2: $t_{\\rm jet,2}$", True)],
     "both jets reach HGTD: $t_{\\rm jet,1}$ vs $t_{\\rm jet,2}$ (jet time compatibility), and each vs $t_0$"),
    ("central_forward", [(0.6, 1.5, "central tagging jet: no HGTD time", False), (-3.0, 3.6, "forward tagging jet: $t_{\\rm jet}$", True)],
     "only the forward jet reaches HGTD: $t_{\\rm jet}$ vs $t_0$ alone"),
]:
    fig, ax = plt.subplots(figsize=(12.5, 3.9))
    detector(ax)
    for eta, length, lab, timed in jets:
        jet(ax, eta, length, lab, timed)
    ax.set_title(note + "\nJVT/fJVT: both jets already compatible with the hard-scatter vertex in $z$",
                 fontsize=12.5, loc="left")
    ax.plot([], [], marker="*", ls="none", ms=13, color="#E69F00", mec="black", label="jet timed by HGTD")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}_{key}.{ext}", dpi=170)
    plt.close(fig)
    print("wrote", f"{args.out}_{key}.pdf / .png")
