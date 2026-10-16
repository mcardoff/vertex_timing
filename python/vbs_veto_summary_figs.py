#!/usr/bin/env python3
"""Summary figures for the timing-veto talk, in place of number tables. Computed directly from
the util/vbs_time_veto trees with the slides' veto (vbs_time_veto_common.t0_veto; jet time
compatibility as vbs_time_veto_plot.decisions()), recoil > 100 GeV, generator-weighted, one
jet-time calibration from the signal's HS legs for both samples.

  cases   per case (Forward-Forward / Central-Forward / all events): the S/sqrt(B) change vs
          no veto per method (top), signal and background efficiency (bottom)
  cf      Central-Forward events by the one jet HGTD times (HS / pileup) and where the other
          jet sits: the fraction the jet-vs-t0 veto keeps, per t0, VBF and Z+jets
  shares  how many tagging jets HGTD sees: the share of VBF and Z+jets in each case
  cone    the jet-time cone, --sig/--bkg against --sig2/--bkg2 (e.g. dR < 0.2 vs 0.4), all
          events: S/sqrt(B) change and efficiencies per method at each cone

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_veto_summary_figs.py cases \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root \\
        --out figs/time_veto_k2/summary_cases
"""
import argparse
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("what", choices=["cases", "cf", "shares", "cone"])
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--sig2", help="cone: the second cone's signal tree")
ap.add_argument("--bkg2", help="cone: the second cone's background tree")
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--tag", default="JVT + fJVT loose")
ap.add_argument("--out", required=True)
args = ap.parse_args()
plt.rcParams.update({"font.size": 12})

METHODS = [("jj", "Jet time compatibility", "black", "o"), ("trkptz", "TRKPTZ $t_0$", "#C0392B", "s"),
           ("waves", "WAVeS $t_0$", "#E69F00", "^"), ("hgtd", "HGTD $t_0$ (Athena)", "#3B8ED0", "v"),
           ("truth", "truth $t_0$ (perfect)", "#8C5A3C", "P")]
CASE_COL = {"both": "#3F8FFF", "one": "#B3231C", "neither": "#CFCFCF"}
FOOT = (f"VBF $H\\to$inv. against $Z(\\to\\ell\\ell)$+jets, {args.tag}, recoil > 100 GeV, "
        f"jet time from $\\geq${args.min_trk} tracks   ATLAS Simulation Internal")


def load(path, rc):
    t = uproot.open(path)["events"]
    e = t.arrays([k for k in t.keys() if not k.startswith("jet_")], library="np")
    k = e[rc] > 100
    return {c: v[k] for c, v in e.items()}, e


def vetoes(e, cal):
    out = {}
    sA, sB = cal(e["a_core_n"]) * e["a_core_sig"], cal(e["b_core_n"]) * e["b_core_sig"]
    both = (e["a_core_n"] >= args.min_trk) & (e["b_core_n"] >= args.min_trk)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["jj"] = both & (np.abs(e["a_core_t"] - e["b_core_t"]) / np.sqrt(sA ** 2 + sB ** 2) >= args.nsigma)
    for s in ["trkptz", "waves", "hgtd", "truth"]:
        out[s] = C.t0_veto(e, cal, s, args.min_trk, args.nsigma)[0]
    return out


def sample(sig, bkg):
    S, S_all = load(sig, "met_truth")
    B, _ = load(bkg, "z_pt")
    cal = C.measure(S_all, True, "signal HS legs")
    return dict(S=S, B=B, vS=vetoes(S, cal), vB=vetoes(B, cal))


def nin(e):
    return (e["a_hgtd_ntimed"] >= 1).astype(int) + (e["b_hgtd_ntimed"] >= 1).astype(int)


def figures(d, ms, mb, m):
    """eps_S, eps_B (with errors) and the S/sqrt(B) change (with error) in the events selected"""
    ws, wb = d["S"]["weight"][ms].astype(float), d["B"]["weight"][mb].astype(float)
    es, ses, _, _ = C.eff(ws, ~d["vS"][m][ms])
    eb, seb, _, _ = C.eff(wb, ~d["vB"][m][mb])
    r = es / np.sqrt(eb)
    return es, ses, eb, seb, r - 1, r * np.hypot(ses / es, 0.5 * seb / eb)


def dodge(i, n, width=0.6):
    return width * ((i + 0.5) / n - 0.5)


def save(fig):
    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}.{ext}", dpi=150, bbox_inches="tight")
    print("wrote", args.out + ".pdf / .png")


d = sample(args.sig, args.bkg)
S, B = d["S"], d["B"]
ws, wb = S["weight"].astype(float), B["weight"].astype(float)

if args.what == "cases":
    cases = [("Forward-Forward\nboth jets in HGTD", nin(S) == 2, nin(B) == 2, "both"),
             ("Central-Forward\none jet in HGTD", nin(S) == 1, nin(B) == 1, "one"),
             ("All selected events", np.ones(len(ws), bool), np.ones(len(wb), bool), None)]
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(10, 7.4), sharex=True, gridspec_kw=dict(height_ratios=[1.15, 1], hspace=0.06))
    for i, (key, lab, col, mk) in enumerate(METHODS):
        for x, (cl, ms, mb, cc) in enumerate(cases):
            if key == "jj" and cc == "one":
                continue      # needs two timed jets: nothing to do in Central-Forward
            es, ses, eb, seb, ch, sch = figures(d, ms, mb, key)
            xx = x + dodge(i, len(METHODS))
            a0.errorbar(xx, 100 * ch, yerr=100 * sch, fmt=mk, color=col, ms=8, capsize=3, lw=1.6,
                        label=lab if x == 0 else None)
            a1.errorbar(xx, 100 * eb, yerr=100 * seb, fmt=mk, color=col, ms=8, capsize=3, lw=1.6)
            a1.plot(xx, 100 * es, mk, color=col, ms=8, mfc="white", mew=1.6)
            print(f"{cl.splitlines()[0]:16s} {key:7s} eps_S {100 * es:5.1f}  eps_B {100 * eb:5.1f} +- {100 * seb:3.1f}  "
                  f"S/sqrt(B) {100 * ch:+5.1f} +- {100 * sch:3.1f}%")
    a0.axhline(0, color="#999999", lw=0.8, ls=":")
    a0.set_ylabel("$S/\\sqrt{B}$ change vs no veto [%]")
    a0.legend(fontsize=10, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.17))
    a1.set_ylabel("efficiency [%]\nfilled: $Z$+jets, open: VBF")
    a1.set_ylim(55, 102)
    for a in (a0, a1):
        a.grid(alpha=0.3, axis="y")
        for x in (0.5, 1.5): a.axvline(x, color="#BBBBBB", lw=0.8)
    ticks = []
    for cl, ms, mb, cc in cases:
        fs, fb = ws[ms].sum() / ws.sum(), wb[mb].sum() / wb.sum()
        ticks.append(f"{cl}\n" + (f"{100 * fs:.1f}% of $S$, {100 * fb:.1f}% of $B$\n" if cc else "\n")
                     + f"({int(mb.sum()):,} $Z$+jets MC)")
    a1.set_xticks(range(len(cases))); a1.set_xticklabels(ticks, fontsize=10)
    a1.set_xlim(-0.5, len(cases) - 0.5)
    fig.text(0.5, -0.07, FOOT + ";  'neither' (2.4% of $S$, 5.1% of $B$) is untouched", ha="center", fontsize=9.5, color="#555555")
    save(fig)

elif args.what == "cf":
    one_s, one_b = nin(S) == 1, nin(B) == 1
    def groups(e):
        ain = e["a_hgtd_ntimed"] >= 1
        pick = lambda col: np.where(ain, e["a_" + col], e["b_" + col])
        hs, pu = pick("hs"), pick("pu")
        oth = np.abs(np.where(ain, e["b_eta"], e["a_eta"]))
        return [("HS in HGTD,\nother central", hs & (oth < 2.4)),
                ("HS in HGTD,\nother $|\\eta| > 2.4$", hs & (oth >= 2.4)),
                ("PU in HGTD,\nother central", pu & (oth < 2.4)),
                ("PU in HGTD,\nother $|\\eta| > 2.4$", pu & (oth >= 2.4))]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True, gridspec_kw=dict(wspace=0.06))
    for ax, (nm, e, one, vk) in zip(axes, [("VBF $H\\to$inv.", S, one_s, "vS"), ("$Z(\\to\\ell\\ell)$+jets", B, one_b, "vB")]):
        w = e["weight"].astype(float)
        G = groups(e)
        for i, (key, lab, col, mk) in enumerate(METHODS[1:]):
            for x, (gl, gm) in enumerate(G):
                m = one & gm
                k, sk, _, _ = C.eff(w[m], ~d[vk][key][m])
                ax.errorbar(x + dodge(i, 4, 0.5), k, yerr=sk, fmt=mk, color=col, ms=8, capsize=3, lw=1.6,
                            label=lab if x == 0 else None)
        ticks = [f"{gl}\n{100 * w[one & gm].sum() / w[one].sum():.1f}% of the case\n({int((one & gm).sum()):,} MC)" for gl, gm in G]
        ax.set_xticks(range(4)); ax.set_xticklabels(ticks, fontsize=9.5)
        ax.set_xlim(-0.5, 3.5); ax.set_ylim(0.4, 1.02); ax.grid(alpha=0.3, axis="y")
        ax.axvline(1.5, color="#999999", lw=1.0)
        ax.set_title(f"{nm}: Central-Forward events", fontsize=13, loc="left", fontweight="bold")
    axes[0].set_ylabel("fraction the jet-vs-$t_0$ veto keeps")
    axes[0].legend(fontsize=10, frameon=False, loc="lower left")
    fig.text(0.5, -0.12, FOOT + ";  events whose HGTD jet has no truth label not shown", ha="center", fontsize=9.5, color="#555555")
    save(fig)

elif args.what == "shares":
    fig, ax = plt.subplots(figsize=(10, 3.0))
    rows = [("VBF $H\\to$inv.", S), ("$Z(\\to\\ell\\ell)$+jets", B)]
    labs = [("both", "Forward-Forward: both in HGTD"), ("one", "Central-Forward: one in HGTD"),
            ("neither", "neither: untouched by any timing veto")]
    for y, (nm, e) in enumerate(rows):
        w = e["weight"].astype(float); n = nin(e); left = 0.0
        for k, (cc, cl) in enumerate(labs):
            f = w[n == {"both": 2, "one": 1, "neither": 0}[cc]].sum() / w.sum()
            ax.barh(y, f, left=left, color=CASE_COL[cc], edgecolor="black", lw=1, label=cl if y == 0 else None)
            if f > 0.08:
                ax.text(left + f / 2, y, f"{100 * f:.1f}%", ha="center", va="center", fontsize=13, color="white", fontweight="bold")
            else:   # too narrow to hold its label: write it just past the bar
                ax.text(1.01, y, f"{100 * f:.1f}%", ha="left", va="center", fontsize=12, color="#444444")
            left += f
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=13)
    ax.invert_yaxis(); ax.set_xlim(0, 1.08); ax.set_xticks([])
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
    ax.legend(fontsize=10.5, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.45, 1.32))
    ax.set_xlabel("share of the selected events (generator-weighted); a jet is in HGTD with $\\geq$ 1 timed ghost track", fontsize=10.5)
    save(fig)

else:   # cone
    if not (args.sig2 and args.bkg2): raise SystemExit("cone needs --sig2 and --bkg2")
    d2 = sample(args.sig2, args.bkg2)
    drs = []
    for path in (args.sig, args.sig2):
        drs.append(float(np.unique(uproot.open(path)["meta"].arrays(["jet_time_dr"], library="np")["jet_time_dr"])[0]))
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(8.6, 6.6), sharex=True, gridspec_kw=dict(height_ratios=[1.15, 1], hspace=0.06))
    alls, allb = np.ones(len(ws), bool), np.ones(len(wb), bool)
    for i, (key, lab, col, mk) in enumerate(METHODS):
        pts = [figures(dd, alls, allb, key) for dd in (d, d2)]
        xs = np.arange(2) + dodge(i, len(METHODS), 0.35)
        a0.errorbar(xs, [100 * p[4] for p in pts], yerr=[100 * p[5] for p in pts], fmt=mk + "-", color=col, ms=8, capsize=3, lw=1.6, label=lab)
        a1.errorbar(xs, [100 * p[2] for p in pts], yerr=[100 * p[3] for p in pts], fmt=mk + "-", color=col, ms=8, capsize=3, lw=1.6)
        a1.plot(xs, [100 * p[0] for p in pts], mk + "--", color=col, ms=8, mfc="white", mew=1.6, lw=1.2)
        print(f"{key:7s} " + "  ".join(f"dR<{dr}: eps_S {100 * p[0]:5.1f} eps_B {100 * p[2]:5.1f} S/sqrt(B) {100 * p[4]:+5.1f}"
                                       for dr, p in zip(drs, pts)))
    a0.axhline(0, color="#999999", lw=0.8, ls=":")
    a0.set_ylabel("$S/\\sqrt{B}$ change vs no veto [%]")
    a0.legend(fontsize=10, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.24))
    a1.set_ylabel("efficiency [%]\nfilled: $Z$+jets, open: VBF")
    a1.set_ylim(55, 102)
    a1.set_xticks([0, 1]); a1.set_xticklabels([f"jet-time cone\n$\\Delta R < {dr:.1f}$" for dr in drs])
    a1.set_xlim(-0.5, 1.5)
    for a in (a0, a1): a.grid(alpha=0.3, axis="y")
    fig.text(0.5, -0.03, f"all selected events, jet time from $\\geq${args.min_trk} tracks   ATLAS Simulation Internal",
             ha="center", fontsize=9.5, color="#555555")
    save(fig)
