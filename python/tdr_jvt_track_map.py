#!/usr/bin/env python3
"""Why a forward pileup jet passes JVT, in one figure: the jet's ghost tracks in (z0 - z_PV0, t - t0), the track(s) the
vertex fit gave PV0 -- the only ones JVT's R_pT counts -- marked, and the WAVeS t0 gate (2 sigma) as a band. For the
2026-11-09 mini deck on a time-gated / tighter forward JVT (slides/slides_jvt_variants.tex).

Events: rows of a tdr_display_events.py candidates.csv, read from the picked RAW events (the cleaning needs
Track_hitInLastLayer). Times: cleaned (the TDR's cuts) as the analysis uses them; a track without a usable time sits in
the "no time" strip. t0: the row's WAVeS t0; the gate band uses sigma = sqrt((infl sigma_t0)^2 + (30 ps)^2), a typical
track's -- each track's own pull is in the printout.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_jvt_track_map.py --csv ../figs/tdr/displays_pu60/candidates.csv \\
      --picked ../condor/tdr_display_events_pu60 --out-dir ../figs/tdr/displays_pu60 --only pu60_caught_2,pu60_kept_4
"""
import argparse
import csv
import itertools
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True); ap.add_argument("--picked", required=True); ap.add_argument("--out-dir", required=True)
ap.add_argument("--only", required=True, help="comma list of names (case_rank)")
ap.add_argument("--gate", type=float, default=2.0)
args = ap.parse_args()

_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "tdr_display_events.py")).read()
_ns = {"np": np, "itertools": itertools, "sys": sys}
exec(_src[_src.index("def jet_index"):_src.index("def cmd_render")], _ns)   # jet_index, core_tracks: the display's own
BR = ["AntiKt4EMTopoJets_pt", "AntiKt4EMTopoJets_eta", "AntiKt4EMTopoJets_phi", "AntiKt4EMTopoJets_ghostTrack_idx",
      "Track_pt", "Track_eta", "Track_phi", "Track_quality", "Track_time", "Track_timeRes", "Track_hasValidTime",
      "Track_hitInLastLayer", "Track_nHGTDHits", "Track_truthVtx_idx", "Track_z0", "Track_recoVtx_idx", "RecoVtx_z"]
rows = {f"{r['case']}_{r['rank']}": r for r in csv.DictReader(open(args.csv))}
mp = {(m["sample"], int(m["file_idx"]), int(m["entry"])): int(m["picked_entry"]) for m in csv.DictReader(open(args.picked + "_map.csv"))}
plt.rcParams.update({"font.size": 13})
TOP = 640.0   # |t - t0| beyond this is drawn at the frame edge; the "no time" strip sits above it

for name in args.only.split(","):
    r = rows[name]
    pe = mp[(r["sample"], int(r["file_idx"]), int(r["entry"]))]
    g = uproot.open(f"{args.picked}_{r['sample']}.root")["ntuple"].arrays(BR, entry_start=pe, entry_stop=pe + 1)
    g = {k: g[k][0].to_list() for k in BR}
    k = 1 if float(r["j1_pu"]) == 1 else 2
    pt_j, eta_j = float(r[f"j{k}_pt"]), float(r[f"j{k}_eta"])
    j = _ns["jet_index"](g, pt_j, eta_j)
    core = set(_ns["core_tracks"](g, j, float(r[f"j{k}_core_t"]), int(float(r[f"j{k}_core_n"]))))
    t0, sg0 = float(r["t0_waves"]), float(r["sig0_waves"]) * float(r["infl_waves"])
    zpv = g["RecoVtx_z"][0]
    gh = g["AntiKt4EMTopoJets_ghostTrack_idx"][j]

    def timed(i):
        if g["Track_hasValidTime"][i] != 1 or not g["Track_hitInLastLayer"][i]: return False
        ae = abs(g["Track_eta"][i])
        return not (3.5 < ae < 3.9 and g["Track_nHGTDHits"][i] < 2)

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    band = args.gate * np.hypot(sg0, 30.0)
    ax.axhspan(-band, band, color="#E69F00", alpha=0.18, lw=0)
    ax.axhline(0.0, color="#E69F00", lw=1.2)
    ax.axvline(0.0, color="#444444", lw=0.8, ls=":")
    pv_sum, lines = 0.0, []
    for i in gh:
        dz = g["Track_z0"][i] - zpv
        tm = timed(i)
        y = np.clip(g["Track_time"][i] - t0, -TOP, TOP) if tm else TOP + 80
        v = g["Track_recoVtx_idx"][i]
        s = 25 + 22 * g["Track_pt"][i]
        if v == 0:
            pv_sum += g["Track_pt"][i]
            ax.scatter([dz], [y], s=4 * s, marker="*", color="#C0392B", edgecolor="black", lw=0.6, zorder=6)
            pull = (g["Track_time"][i] - t0) / np.hypot(sg0, g["Track_timeRes"][i]) if tm else np.nan
            lab = (f"fitted to PV0: {g['Track_pt'][i]:.1f} GeV, {dz:+.1f} mm"
                   + (f", {g['Track_time'][i] - t0:+.0f} ps ({abs(pull):.1f}$\\sigma$)" if tm else ", no time")
                   + f"  [truth: {'hard scatter' if g['Track_truthVtx_idx'][i] == 0 else 'pileup'}]")
            lines.append(lab)
        else:
            ax.scatter([dz], [y], s=s, marker="o", color="#7F7F7F" if v > 0 else "white",
                       edgecolor="#555555", lw=0.8, alpha=0.85, zorder=3)
        if i in core:
            ax.scatter([dz], [y], s=s + 160, marker="D", facecolor="none", edgecolor="#008B8B", lw=1.6, zorder=5)
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlim(-150, 150)
    ax.set_xticks([-100, -10, -1, 0, 1, 10, 100]); ax.set_xticklabels(["$-$100", "$-$10", "$-$1", "0", "1", "10", "100"])
    ax.set_ylim(-TOP - 40, TOP + 140)
    ax.axhspan(TOP + 40, TOP + 140, color="#EEEEEE", lw=0, zorder=0)
    ax.text(140, TOP + 90, "no time", va="center", ha="right", fontsize=11, color="#666666")
    ax.set_xlabel("track $z_0 - z_{\\rm PV0}$ [mm]   (linear within $\\pm$1 mm)")
    ax.set_ylabel("track time $-$ WAVeS $t_0$ [ps]")
    ax.set_yticks([-600, -400, -200, 0, 200, 400, 600])
    rpt = pv_sum / pt_j
    head = (f"{'Z+jets' if r['sample'] == 'zjets' else 'VBF'}: the pileup VBS jet, {pt_j:.0f} GeV at $\\eta$ = {eta_j:+.2f}, "
            f"{len(gh)} ghost tracks\nJVT $R_{{p_T}}$ = {pv_sum:.1f} / {pt_j:.0f} GeV = {rpt:.3f} "
            f"{'$\\geq$' if rpt >= 0.0167 else '<'} 0.0167: {'passes' if rpt >= 0.0167 else 'fails'}")
    ax.set_title(head, fontsize=12.5, loc="left")
    from matplotlib.lines import Line2D
    hand = [Line2D([], [], marker="*", ls="", color="#C0392B", mec="black", ms=15, label="fitted to PV0 (JVT counts these)"),
            Line2D([], [], marker="o", ls="", color="#7F7F7F", mec="#555555", ms=8, label="fitted to a pileup vertex"),
            Line2D([], [], marker="o", ls="", color="white", mec="#555555", ms=8, label="fitted to no vertex"),
            Line2D([], [], marker="D", ls="", mfc="none", mec="#008B8B", ms=10, mew=1.6, label="the jet's time (core tracks)"),
            plt.Rectangle((0, 0), 1, 1, color="#E69F00", alpha=0.3, label=f"$\\pm${args.gate:g}$\\sigma$ t$_0$ gate (30 ps track)")]
    ax.legend(handles=hand, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=10.5, frameon=False)
    fig.text(0.125, -0.02, "\n".join(lines), fontsize=11.5, color="#C0392B", va="top")
    out = os.path.join(args.out_dir, f"jvtmap_{name}")
    fig.savefig(out + ".png", dpi=170, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"{name}: R_pT {rpt:.3f} from {pv_sum:.1f} GeV; " + " | ".join(lines) + f"  -> {out}.png")
