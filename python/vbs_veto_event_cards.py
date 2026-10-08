#!/usr/bin/env python3
"""Event cards for the jet-vs-t0 timing veto: one figure per event, drawn from the
util/vbs_time_veto trees alone (no ntuple access).

  left    the event's jets in eta-phi (radius 0.4), hard-scatter green, pileup grey,
          the tagging pair outlined; HGTD's 2.38 < |eta| < 4.0 shaded
  right   the time axis: the truth hard-scatter time, each reco t0 (TRKPTZ, WAVeS,
          HGTD) and each tagging jet's time, with the +-1 sigma the test uses (t0:
          quoted x inflation; jet: quoted x f(n), f measured on the signal's HS legs
          as in vbs_time_veto_plot.py), and each jet-vs-t0 pull with its verdict

Events are chosen by fixed rules, one per finding of the study (see CASES), after
the same recoil cuts as the slides; the rule and the event are printed.

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_veto_event_cards.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root \\
        --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root --out figs/time_veto/cards/card
"""
import argparse, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 13})   # cards are shown one per slide at ~half size
from matplotlib.patches import Circle, Rectangle
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--met-min", type=float, default=100.0)
ap.add_argument("--zpt-min", type=float, default=100.0)
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--out", required=True, help="output stem: <out>_<case>.(pdf|png)")
args = ap.parse_args()

SRC = [("trkptz", "TRKPTZ $t_0$", "#C0392B"), ("waves", "WAVeS $t_0$", "#E69F00"), ("hgtd", "HGTD $t_0$ (Athena)", "#3B8ED0")]


def load(path, col, cut):
    t = uproot.open(path)["events"]
    keys = [k for k in t.keys() if not k.startswith("jet_")]
    e_all = t.arrays(keys, library="np")
    J = t.arrays(["jet_idx", "jet_pt", "jet_eta", "jet_phi", "jet_hs", "jet_pu", "jet_t", "jet_sig", "jet_n"], library="np")
    k = e_all[col] > cut
    e = {c: v[k] for c, v in e_all.items()}
    Jk = {c: v[k] for c, v in J.items()}
    return e, Jk, e_all


S, SJ, S_all = load(args.sig, "met_truth", args.met_min)
B, BJ, _ = load(args.bkg, "z_pt", args.zpt_min)
cal = C.measure(S_all, True, "signal HS legs")   # one jet-time calibration for S and B, as the slides


def derived(e):
    d = {}
    for L in "ab":
        d[L + "_timed"] = e[L + "_core_n"] >= 1
        d[L + "_s"] = cal(e[L + "_core_n"]) * e[L + "_core_sig"]
        d[L + "_in"] = e[L + "_hgtd_ntimed"] >= 1
    for s, _, _ in SRC + [("truth", "", "")]:
        s0 = e[f"infl_{s}"] * e[f"sig0_{s}"]
        veto = np.zeros(len(s0), bool)
        for L in "ab":
            pull = np.abs(e[L + "_core_t"] - e[f"t0_{s}"]) / np.sqrt(d[L + "_s"] ** 2 + s0 ** 2)
            d[f"{L}_pull_{s}"] = np.where(d[L + "_timed"] & e[f"ok_{s}"], pull, np.nan)
            veto |= e[f"ok_{s}"] & d[L + "_timed"] & (pull >= args.nsigma)
        d["veto_" + s] = veto
    d["nin"] = d["a_in"].astype(int) + d["b_in"].astype(int)
    return d


SD, BD = derived(S), derived(B)


def pick(e, d, mask, score, rule):
    """the event in mask with the largest score (deterministic)"""
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        raise SystemExit(f"no event for rule: {rule}")
    i = idx[np.argmax(score[idx])]
    print(f"  {len(idx):,} candidates; picked entry {e['entry'][i]} of {os.path.basename(str(e['file_path'][i]))}")
    return i


def short_file(path):
    """zjets_skim_shard3of20_0.root stays; the raw VBF user.mcardiff.<task>.Output._002625.SuperNtuple.root -> raw _002625"""
    b = os.path.basename(str(path))
    if b.startswith("user.") and "._" in b:
        return "raw ntuple _" + b.split("._")[-1].split(".")[0]
    return b


def leg_role(e, L, i):
    return "HS" if e[L + "_hs"][i] else ("PU" if e[L + "_pu"][i] else "unlabelled")


# ── the cases: (key, sample, title, rule, mask, score) ───────────────────────
cases = []
both, one = BD["nin"] == 2, BD["nin"] == 1
hspu = (B["a_hs"] & B["b_pu"]) | (B["a_pu"] & B["b_hs"])
cases.append(("A_both_pu_caught", "bkg",
              "$Z$+jets, both tagging jets in HGTD: the pileup jet is caught by the reco $t_0$s; Athena has no time",
              "Z+jets, both legs in HGTD, HS+PU pair, both legs timed; vetoed by TRKPTZ and WAVeS, no valid HGTD t0 "
              "(every such event a valid HGTD t0 keeps has none: 5 of 5); highest PU-leg pT",
              both & hspu & BD["a_timed"] & BD["b_timed"] & BD["veto_trkptz"] & BD["veto_waves"] & ~B["ok_hgtd"],
              np.where(B["a_pu"], B["a_pt"], B["b_pt"])))
pu_in = np.where(BD["a_in"], B["a_pu"], B["b_pu"])
oth_hs_cen = np.where(BD["a_in"], B["b_hs"] & (np.abs(B["b_eta"]) < 2.4), B["a_hs"] & (np.abs(B["a_eta"]) < 2.4))
cases.append(("B_one_r2_caught", "bkg",
              "$Z$+jets, one tagging jet in HGTD: a forward pileup jet beside a central hard-scatter jet, caught",
              "Z+jets, exactly one leg in HGTD and it is PU, the other HS with |eta| < 2.4; vetoed by TRKPTZ, WAVeS and HGTD; "
              "highest PU-leg pT",
              one & pu_in & oth_hs_cen & BD["veto_trkptz"] & BD["veto_waves"] & BD["veto_hgtd"],
              np.where(BD["a_in"], B["a_pt"], B["b_pt"])))
sboth = (SD["nin"] == 2) & S["a_hs"] & S["b_hs"] & SD["a_timed"] & SD["b_timed"]
fail1 = np.where(SD["a_pull_truth"] >= args.nsigma, S["a_core_n"], np.where(SD["b_pull_truth"] >= args.nsigma, S["b_core_n"], 0))
worst = np.fmax(SD["a_pull_truth"], SD["b_pull_truth"])
cases.append(("C_sig_jet_tail", "sig",
              "VBF $H\\to$inv., both jets hard-scatter and in HGTD: lost to the jet-time tail, even with the true $t_0$",
              "VBF, both legs HS, in HGTD and timed; vetoed by WAVeS AND by the truth t0, the failing leg single-track; "
              "pull vs truth closest to 4",
              sboth & SD["veto_waves"] & SD["veto_truth"] & (fail1 == 1),
              -np.abs(worst - 4.0)))
hshs = B["a_hs"] & B["b_hs"] & (BD["a_timed"] | BD["b_timed"])
with np.errstate(divide="ignore", invalid="ignore"):
    t0off = np.abs(B["t0_waves"] - B["t0_truth"]) / (B["infl_waves"] * B["sig0_waves"])
cases.append(("D_bkg_wrong_t0", "bkg",
              "$Z$+jets, both tagging jets hard-scatter: vetoed only because the reco $t_0$ is wrong",
              "Z+jets, both legs HS, >= 1 timed; vetoed by WAVeS and TRKPTZ, kept by the truth t0, WAVeS t0 > 3 sigma from truth; "
              "highest m_jj",
              hshs & BD["veto_waves"] & BD["veto_trkptz"] & ~BD["veto_truth"] & (t0off > 3.0),
              B["mjj"]))
edgeA = B["a_pu"] & (np.abs(B["a_eta"]) > 4.0) & BD["a_in"] & BD["a_timed"]
edgeB = B["b_pu"] & (np.abs(B["b_eta"]) > 4.0) & BD["b_in"] & BD["b_timed"]
edge_pull = np.where(edgeA, BD["a_pull_waves"], BD["b_pull_waves"])
cases.append(("E_edge_pu_caught", "bkg",
              "$Z$+jets: a pileup tagging jet with its axis past $|\\eta| = 4.0$, timed by HGTD and caught",
              "Z+jets, a PU leg with |eta_axis| > 4.0 and a timed HGTD track; vetoed by TRKPTZ and WAVeS on THAT leg; "
              "highest edge-leg |eta|",
              (edgeA | edgeB) & BD["veto_waves"] & BD["veto_trkptz"] & (edge_pull >= args.nsigma),
              np.where(edgeA, np.abs(B["a_eta"]), np.abs(B["b_eta"]))))


def draw_times(at, e, d, i, names):
    """the time axis: one row per t0 (truth, TRKPTZ, WAVeS, HGTD), then one per
    tagging jet, each with the +-1 sigma the test uses, and each t0's verdict"""
    tt = e["t0_truth"][i]
    rows = [("truth $t_{\\rm HS}$", tt, 0.0, "#5B3A29", "truth")]
    for s, lab, col in SRC:
        ok = e[f"ok_{s}"][i]
        rows.append((lab + ("" if ok else "\n(no time)"), e[f"t0_{s}"][i] if ok else np.nan,
                     e[f"infl_{s}"][i] * e[f"sig0_{s}"][i] if ok else np.nan, col, s))
    jets = []
    for L, name in zip("ab", names):
        role_ = leg_role(e, L, i)
        txt = f"{name}: {role_}, $\\eta$ = {e[L + '_eta'][i]:+.2f}\n{e[L + '_hgtd_ntimed'][i]} timed HGTD tracks"
        if d[L + "_timed"][i]:
            jets.append((txt + f", time from {e[L + '_core_n'][i]}", e[L + "_core_t"][i], d[L + "_s"][i],
                         "#2E9E4F" if role_ == "HS" else "#555555", L))
        else:
            jets.append((txt + "\nno jet time", np.nan, np.nan, "#999999", L))
    allrows = rows + jets
    ys = np.arange(len(allrows))[::-1]
    vals = [v for _, v, _, _, _ in allrows if np.isfinite(v)]
    sigs = [x for _, v, x, _, _ in allrows if np.isfinite(v) and np.isfinite(x)]
    lo, hi = min(vals) - 3 * max(sigs + [20]), max(vals) + 3 * max(sigs + [20])
    pad = 0.08 * (hi - lo)
    for y, (lab, v, sg, col, tag) in zip(ys, allrows):
        if np.isfinite(v):
            at.errorbar([v], [y], xerr=[[sg], [sg]] if sg > 0 else None, fmt="s" if tag in [x[0] for x in SRC] else "o",
                        color=col, ms=8, capsize=4, lw=2)
        if tag in [x[0] for x in SRC] + ["truth"] and tag != "truth" or tag == "truth":
            if not e[f"ok_{tag}"][i]:
                verdict, vcol = "no $t_0$: kept", "#777777"
            else:
                p = [d[f"{L}_pull_{tag}"][i] for L in "ab"]
                ps = ", ".join("–" if not np.isfinite(x) else f"{x:.1f}$\\sigma$" for x in p)
                verdict = ("VETOED" if d["veto_" + tag][i] else "kept") + f"  ({ps})"
                vcol = "#B03A2E" if d["veto_" + tag][i] else "#2E7D32"
            at.text(1.015, y, verdict, transform=at.get_yaxis_transform(), ha="left", va="center", fontsize=12,
                    color=vcol, fontweight="bold")
    at.set_yticks(ys)
    at.set_yticklabels([r[0] for r in allrows], fontsize=12)
    for tl, r in zip(at.get_yticklabels(), allrows): tl.set_color(r[3])
    at.axhline(len(jets) - 0.5, color="#BBBBBB", lw=0.8)
    at.axvline(tt, color="#5B3A29", lw=0.8, ls=":")
    at.set_xlim(lo - pad, hi + pad); at.set_ylim(-0.6, len(allrows) - 0.4)
    at.set_xlabel("time [ps]   (bars: $\\pm 1\\sigma$ as used by the test)")
    at.text(1.015, len(allrows) - 0.25, f"jet vs $t_0$ at {args.nsigma:g}$\\sigma$\n({names[0]}, {names[1]} pulls)",
            transform=at.get_yaxis_transform(), ha="left", va="bottom", fontsize=11, color="#555555")


def card(key, role, title, rule, mask, score):
    e, d, J = (S, SD, SJ) if role == "sig" else (B, BD, BJ)
    print(f"[{key}] {rule}")
    i = pick(e, d, mask, score, rule)
    fig = plt.figure(figsize=(13.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.75, left=0.05, right=0.83, top=0.82, bottom=0.13)
    # ── eta-phi ──
    ax = fig.add_subplot(gs[0])
    for sgn in (-1, 1):
        ax.add_patch(Rectangle((min(sgn * 2.38, sgn * 4.0), -np.pi), 4.0 - 2.38, 2 * np.pi, color="#DCE6F2", zorder=0))
        ax.text(sgn * 3.19, -np.pi * 0.95, "HGTD", ha="center", va="bottom", fontsize=12, color="#4A6FA5")
    legs = {int(e["a_idx"][i]): "a", int(e["b_idx"][i]): "b"}
    for k in range(len(J["jet_pt"][i])):
        eta, phi, pt = J["jet_eta"][i][k], J["jet_phi"][i][k], J["jet_pt"][i][k]
        hs, pu = J["jet_hs"][i][k] > 0, J["jet_pu"][i][k] > 0
        col = "#2E9E4F" if hs else ("#8C8C8C" if pu else "#C9B458")
        leg = legs.get(int(J["jet_idx"][i][k]))
        for dphi in (0.0, 2 * np.pi, -2 * np.pi):   # wrap round phi = +-pi
            ax.add_patch(Circle((eta, phi + dphi), 0.4, facecolor=col, alpha=0.55 if leg else 0.25,
                                edgecolor="black" if leg else col, lw=2.2 if leg else 0.8, zorder=3))
        lab = (("jet 1" if leg == "a" else "jet 2") + f"\n{pt:.0f} GeV") if leg else f"{pt:.0f}"
        ty = min(max(phi, -np.pi + 0.45), np.pi - 0.45)
        ax.text(eta, ty, lab, ha="center", va="center", fontsize=11 if leg else 9, zorder=4,
                fontweight="bold" if leg else "normal")
    ax.set_xlim(-4.6, 4.6); ax.set_ylim(-np.pi, np.pi); ax.set_aspect("equal")
    ax.set_xlabel("jet $\\eta$"); ax.set_ylabel("jet $\\phi$")
    ax.set_title("green: hard scatter, grey: pileup\noutlined: tagging pair", fontsize=12)
    # ── time axis: one row per t0, then one per tagging jet ──
    draw_times(fig.add_subplot(gs[1]), e, d, i, ("jet 1", "jet 2"))
    rec = (f"truth $E_T^{{\\rm miss}}$ = {e['met_truth'][i]:.0f} GeV" if role == "sig" else f"$p_T^{{\\ell\\ell}}$ = {e['z_pt'][i]:.0f} GeV")
    # the slide carries the title; the card keeps only the event's identity
    fig.text(0.02, 0.955, f"{'VBF H→inv.' if role == 'sig' else 'Z(→ℓℓ)+jets'}  ·  {short_file(e['file_path'][i])}  "
             f"entry {e['entry'][i]}  ·  $m_{{jj}}$ = {e['mjj'][i]:.0f} GeV  ·  {rec}  ·  ATLAS Simulation Internal",
             fontsize=12, color="#555555")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}_{key}.{ext}", dpi=150)
    plt.close(fig)
    print(f"  wrote {args.out}_{key}.pdf / .png")
    # the time panel alone, beside an event display: its jets are named as event_display.py
    # numbers them (ntuple jet index + 1)
    f2 = plt.figure(figsize=(8.6, 4.4))
    at = f2.add_axes([0.35, 0.15, 0.40, 0.72])
    draw_times(at, e, d, i, (f"Jet {int(e['a_idx'][i]) + 1}", f"Jet {int(e['b_idx'][i]) + 1}"))
    for ext in ("pdf", "png"):
        f2.savefig(f"{args.out}_{key}_times.{ext}", dpi=170, bbox_inches="tight", pad_inches=0.08)
    plt.close(f2)
    print(f"  wrote {args.out}_{key}_times.pdf / .png")
    return dict(key=key, rule=rule, file=os.path.basename(str(e["file_path"][i])), entry=int(e["entry"][i]))


picked = [card(*c) for c in cases]
with open(args.out + "_events.md", "w") as f:
    f.write("# Timing-veto event cards\n\n| card | rule | file | entry |\n|---|---|---|---|\n")
    for p in picked:
        f.write(f"| {p['key']} | {p['rule']} | {p['file']} | {p['entry']} |\n")
print("wrote", args.out + "_events.md")
