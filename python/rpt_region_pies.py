#!/usr/bin/env python3
"""Why the R_pT time gate gets a jet's tracks wrong, as pies: the HS pT it
removes from HS jets, and the PU pT it lets survive in PU jets. One pie per
real-t0 case (TRKPTZ, WAVeS, HGTD Athena), from rpt_v5_hist's region tree.

Every attribution is judged against the truth HS vertex time ONLY, through one
counterfactual -- the same gate (same sigma_t0, inflation, 3 sigma) with t0
moved to t_HS (see RegionRow in util/rpt_v5_hist.cxx):

  HS pT removed from an HS jet (R1 forward HS legs, HS-vertex tracks):
    incorrect vertex t0      kept with the true t0; split into a wrong cluster
                             (|t0 - t_HS| >= 60 ps) and a small offset
    assigned track time      removed even with the true t0: the track's own
      incorrect              time is incompatible with the HS time
  PU pT surviving in a PU jet (R1 + R2 forward PU legs, non-HS tracks):
    no HGTD time             never gated
    no vertex t0             this case had no t0 in the event: no gate at all
    compatible with the      kept even with the true t0 -- timing cannot
      HS time                remove it at this gate width
    incorrect vertex t0      removed with the true t0 (wrong cluster / offset)

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/rpt_region_pies.py \\
        condor/vbf/vbf_jvtLoose_rpt_v5_regions.root --sample vbf_jvtLoose \\
        --label "VBF H#rightarrowinv." --out-dir figs/rpt_regions
"""
import argparse, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CASES = [("trkptz", "TRKPTZ $t_0$"), ("waves", "WAVeS $t_0$"), ("hgtd", "HGTD $t_0$ (Athena)")]

# Colour follows the cause, identically in both figures (validated with the
# dataviz skill's validate_palette.js: every pair that shares a figure clears
# the CVD and normal-vision floors). The two neutrals mean "no timing
# information"; the lighter one is also hatched so it never rests on colour.
T0_FAR, T0_NEAR = "#1c5cab", "#6da7ec"
HS_SLICES = [  # (tree stem, label, colour, hatch)
    ("hs_rm_t0far",  "Incorrect vertex $t_0$: wrong cluster ($|t_0-t_{HS}| \\geq$ 60 ps)", T0_FAR,  None),
    ("hs_rm_t0near", "Incorrect vertex $t_0$: small offset ($|t_0-t_{HS}|$ < 60 ps)",       T0_NEAR, None),
    ("hs_rm_trk",    "Assigned track time is incorrect (removed even with $t_0 = t_{HS}$)", "#eb6834", None),
]
PU_SLICES = [
    ("pu_keep_untimed", "Track has no HGTD time (never gated)",                              "#898781", None),
    ("pu_keep_not0",    "No vertex $t_0$ in the event (no gate)",                            "#c3c2b7", "////"),
    ("pu_keep_intime",  "Track time compatible with $t_{HS}$ (kept even with $t_0 = t_{HS}$)", "#eda100", None),
    ("pu_keep_t0far",   "Incorrect vertex $t_0$: wrong cluster ($|t_0-t_{HS}| \\geq$ 60 ps)", T0_FAR,    None),
    ("pu_keep_t0near",  "Incorrect vertex $t_0$: small offset ($|t_0-t_{HS}|$ < 60 ps)",       T0_NEAR,   None),
]
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("file", help="merged <prefix>rpt_v5_regions.root")
ap.add_argument("--sample", required=True, help="output-name prefix, e.g. vbf_jvtLoose")
ap.add_argument("--label", default="VBF H$\\rightarrow$inv.", help="sample name for the subtitle")
ap.add_argument("--out-dir", default="figs/rpt_regions")
args = ap.parse_args()
args.label = args.label.replace("#rightarrow", "$\\rightarrow$")
os.makedirs(args.out_dir, exist_ok=True)

t = uproot.open(args.file)["regions"]
a = t.arrays(library="np")
if "hs_rm_trk_trkptz" not in a:
    raise SystemExit(f"{args.file} predates the attribution columns (rpt_v5_hist c57a46f); rerun it")
reg = a["region"]
jvt = ["none", "loose", "tight"][int(a["jvt_wp"][0])] if "jvt_wp" in a else "none"
sel_txt = (f"max-$m_{{jj}}$ pair, $m_{{jj}}$ > {float(a['vbs_mjj_cut'][0]):g} GeV"
           + (f", JVT + fJVT {jvt} before pairing" if jvt != "none" else ""))


def slice_values(stems, case, mask):
    vals = []
    for stem, *_ in stems:
        col = stem if stem in a else f"{stem}_{case}"
        vals.append(float(a[col][mask].sum()))
    return np.array(vals)


def draw(stems, title, mask, total_col, what, out_stem, n_jets_txt):
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.6))
    fig.patch.set_facecolor(SURFACE)
    total_cone = float(a[total_col][mask].sum())
    rows = []
    for ax, (case, case_lbl) in zip(axes, CASES):
        v = slice_values(stems, case, mask)
        tot = v.sum()
        ax.set_facecolor(SURFACE)
        if tot <= 0:
            ax.text(0.5, 0.5, "nothing to attribute", ha="center", va="center", color=INK2)
            ax.axis("off"); continue
        frac = v / tot
        keep = v > 0
        wedges, _ = ax.pie(v[keep], colors=[s[2] for s, k in zip(stems, keep) if k],
                           startangle=90, counterclock=False,
                           wedgeprops=dict(edgecolor=SURFACE, linewidth=2))
        for w, (stem, _, _, hatch) in zip(wedges, [s for s, k in zip(stems, keep) if k]):
            if hatch:
                w.set_hatch(hatch)
                w.set_edgecolor(SURFACE)
        # Direct labels: percentages inside the large wedges, outside with a
        # short leader for the thin ones, always in ink rather than the fill.
        for w, f in zip(wedges, frac[keep]):
            ang = np.deg2rad((w.theta1 + w.theta2) / 2)
            x, y = np.cos(ang), np.sin(ang)
            txt = f"{100 * f:.0f}%" if f >= 0.01 else f"{100 * f:.1f}%"
            if f >= 0.07:
                ax.text(0.66 * x, 0.66 * y, txt, ha="center", va="center", fontsize=13,
                        fontweight="bold", color="white" if w.get_facecolor()[0] < 0.5 else INK)
            else:
                ax.annotate(txt, xy=(0.98 * x, 0.98 * y), xytext=(1.22 * x, 1.22 * y),
                            ha="center", va="center", fontsize=11, color=INK,
                            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
        ax.set_title(case_lbl, fontsize=15, color=INK, pad=14)
        ax.text(0, -1.42, f"{what}: {tot:,.0f} GeV\n= {100 * tot / total_cone:.1f}% of the cones' {total_col_txt[total_col]}",
                ha="center", va="top", fontsize=10.5, color=INK2)
        ax.set_aspect("equal")
        rows.append((case, tot, frac))
    fig.suptitle(title, fontsize=19, color=INK, y=0.985)
    fig.text(0.5, 0.905, f"{args.label}, {n_jets_txt}; {sel_txt}; judged against the truth HS vertex time $t_{{HS}}$",
             ha="center", fontsize=11.5, color=INK2)
    fig.text(0.012, 0.975, "ATLAS", fontsize=15, fontweight="bold", fontstyle="italic", color=INK, va="top")
    fig.text(0.075, 0.975, "Simulation Internal", fontsize=15, color=INK, va="top")
    handles = [Patch(facecolor=c, edgecolor=INK2 if h else SURFACE, hatch=h, label=l) for _, l, c, h in stems]
    # two columns: three long cause labels in one row overrun the figure width
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=11.5, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.03, right=0.97, top=0.80, bottom=0.25, wspace=0.25)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.out_dir, f"{out_stem}.{ext}"), dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return rows


total_col_txt = {"hs_hspt": "HS-track pT", "pu_pupt": "PU-track pT"}
r1, pu_legs = reg == 1, (reg == 1) | (reg == 2)
hs_rows = draw(HS_SLICES, "Cases where HS $p_T$ gets removed from a HS jet", r1, "hs_hspt",
               "HS pT removed", f"{args.sample}_hs_pt_removed_pie",
               f"R1 forward HS legs ({r1.sum():,} jets)")
pu_rows = draw(PU_SLICES, "Cases where PU $p_T$ survives in a PU jet", pu_legs, "pu_pupt",
               "PU pT kept", f"{args.sample}_pu_pt_survives_pie",
               f"R1 + R2 forward PU legs ({pu_legs.sum():,} jets)")
for name, stems, rows in (("HS pT removed from HS jets", HS_SLICES, hs_rows),
                          ("PU pT surviving in PU jets", PU_SLICES, pu_rows)):
    print(f"\n{args.sample}: {name}")
    for case, tot, frac in rows:
        print(f"  {case:7s} total {tot:10,.0f} GeV: " + ", ".join(
            f"{s[0].replace('hs_rm_', '').replace('pu_keep_', '')} {100 * f:.1f}%" for s, f in zip(stems, frac)))
print(f"\nwrote {args.out_dir}/{args.sample}_hs_pt_removed_pie.(pdf|png) and {args.sample}_pu_pt_survives_pie.(pdf|png)")
