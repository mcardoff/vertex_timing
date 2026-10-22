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
  budget  the S/sqrt(B) accounting (printed, no figure): an oracle that vetoes every truth-pileup
          tagging jet, then -- applied in this order -- only those HGTD reaches, those with a jet
          time, from >= --min-trk tracks, >= nsigma from the true t0; then the real test against
          the true t0 (adds hard-scatter / unlabelled false vetoes) and against TRKPTZ / WAVeS
  money   the money plot, per m_jj column and for all events: the S/sqrt(B) change from no veto up
          to an oracle that vetoes every truth-pileup tagging jet, split into the gain achieved
          (WAVeS), our failures (WAVeS vetoes an ideal veto would not make), the loss reducible by
          a perfect t0 and by a perfect jet time, and the irreducible loss. The ideal veto: the
          true t0 and pure jet times (util/vbs_time_veto --pure-times: a hard-scatter jet timed by
          its own tracks, a pileup jet by its own collision's; unlabelled jets untested). Needs the
          _pure trees.
  puretime  the perfect-assignment jet times (--pure-times trees) against the true t_HS: hard-scatter
          tagging jets (their hard-scatter tracks) and pileup tagging jets (their own, non-HS
          tracks), all selected events -- what the "+ perfect jet-time assignment" rung tests
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
ap.add_argument("what", choices=["cases", "cf", "shares", "cone", "budget", "money", "puretime"])
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

elif args.what == "budget":
    cal = C.measure(load(args.sig, "met_truth")[1], True, "signal HS legs")
    def chain(e):
        s0 = e["infl_truth"] * e["sig0_truth"]
        lg = {}
        for L in "ab":
            sg = cal(e[f"{L}_core_n"]) * e[f"{L}_core_sig"]
            with np.errstate(divide="ignore", invalid="ignore"):
                pt = np.abs(e[f"{L}_core_t"] - e["t0_truth"]) / np.sqrt(sg ** 2 + s0 ** 2)
            tm = e[f"{L}_core_n"] >= args.min_trk
            lg[L] = dict(pu=e[f"{L}_pu"], reach=e[f"{L}_hgtd_ntimed"] >= 1, time=e[f"{L}_core_ncand"] >= 1,
                         timed=tm, out=tm & (pt >= args.nsigma))
        anyl = lambda f: f(lg["a"]) | f(lg["b"])
        return [("no veto", np.zeros(len(e["weight"]), bool)),
                ("oracle: veto every truth-pileup tagging jet", anyl(lambda x: x["pu"])),
                ("... only those HGTD reaches (>= 1 timed track)", anyl(lambda x: x["pu"] & x["reach"])),
                ("... with a jet time (a timed track within the cone)", anyl(lambda x: x["pu"] & x["time"])),
                (f"... from >= {args.min_trk} tracks", anyl(lambda x: x["pu"] & x["timed"])),
                (f"... out of time (>= {args.nsigma:g} sigma from the true t0)", anyl(lambda x: x["pu"] & x["out"])),
                ("the test against the true t0 (adds false vetoes)", d["vS" if e is S else "vB"]["truth"]),
                ("the test against the WAVeS t0", d["vS" if e is S else "vB"]["waves"]),
                ("the test against the TRKPTZ t0", d["vS" if e is S else "vB"]["trkptz"])]
    # each step against the one before it, from the ROUNDED cumulative values so the steps sum exactly;
    # the WAVeS and TRKPTZ lines both against the true-t0 line
    prev, truth = None, None
    print(f"{'step':58s} {'eps_S':>7s} {'eps_B':>7s} {'S/sqrt(B)':>10s} {'step':>7s}")
    for (lab, vs), (_, vb) in zip(chain(S), chain(B)):
        es = ws[~vs].sum() / ws.sum(); eb = wb[~vb].sum() / wb.sum(); g = round(100 * (es / np.sqrt(eb) - 1), 1)
        ref = truth if ("WAVeS" in lab or "TRKPTZ" in lab) else prev
        print(f"{lab:58s} {100 * es:6.1f}% {100 * eb:6.1f}% {g:+9.1f}% " + (f"{g - ref:+6.1f}" if ref is not None else ""))
        prev = g
        if "true t0" in lab: truth = g

elif args.what == "money":
    if "a_pure_hs_n" not in S: raise SystemExit("money needs the --pure-times trees (*_pure_vbs_time_veto.root)")
    S_all = load(args.sig, "met_truth")[1]
    # the pure times' sigma, calibrated like f(n): on the signal's HS legs, against the true t0
    pure = {"t0_truth": S_all["t0_truth"]}
    for L in "ab":
        pure.update({f"{L}_hs": S_all[f"{L}_hs"], f"{L}_core_n": S_all[f"{L}_pure_hs_n"],
                     f"{L}_core_t": S_all[f"{L}_pure_hs_t"], f"{L}_core_sig": S_all[f"{L}_pure_hs_sig"]})
    pcal = C.measure(pure, True, "pure HS legs")
    def ideal(e):
        """the true t0 and perfect-assignment jet times; unlabelled legs untested"""
        v = np.zeros(len(e["weight"]), bool)
        for L in "ab":
            for lab, key in (("hs", "pure_hs"), ("pu", "pure_pu")):
                n = e[f"{L}_{key}_n"]
                sg = pcal(n) * e[f"{L}_{key}_sig"]
                with np.errstate(divide="ignore", invalid="ignore"):
                    pull = np.abs(e[f"{L}_{key}_t"] - e["t0_truth"]) / sg
                v |= e[f"{L}_{lab}"] & (n >= 1) & (pull >= args.nsigma)
        return v
    oracle = lambda e: e["a_pu"] | e["b_pu"]
    V = {r: dict(waves=d["v" + r]["waves"], t0=d["v" + r]["truth"], ideal=ideal(d[r]), oracle=oracle(d[r])) for r in "SB"}
    # between the ideal veto and the oracle: perfect time resolution (every pileup tagging jet with a timed
    # track in HGTD caught), then perfect HGTD timing efficiency (every one with a track in HGTD's range,
    # timed or not); what is left is HGTD's acceptance (pileup jets with no track in 2.38 < |eta| < 4.0)
    for r in "SB":
        e = d[r]
        pu_with = lambda col: (e["a_pu"] & (e["a_" + col] >= 1)) | (e["b_pu"] & (e["b_" + col] >= 1))
        V[r]["res"] = V[r]["ideal"] | pu_with("pure_pu_n")
        V[r]["eff"] = V[r]["res"] | pu_with("hgtd_n")
    edges = np.array([200, 1000, 2000, 4000.])
    cs, cb = C.mjj_index(S["mjj"], edges), C.mjj_index(B["mjj"], edges)
    cols = [(f"{edges[k]:.0f}-{edges[k + 1]:.0f} GeV" if k < 2 else f"> {edges[k]:.0f} GeV", cs == k, cb == k) for k in range(3)]
    cols.append(("all events", np.ones(len(ws), bool), np.ones(len(wb), bool)))
    def Z(vs, vb, ms, mb):
        es = ws[ms & ~vs].sum() / ws[ms].sum(); eb = wb[mb & ~vb].sum() / wb[mb].sum()
        return 100 * (es / np.sqrt(eb) - 1)
    rows = []
    for lab, ms, mb in cols:
        A, Bv = V["S"], V["B"]
        net = Z(A["waves"], Bv["waves"], ms, mb)
        cor = Z(A["waves"] & A["ideal"], Bv["waves"] & Bv["ideal"], ms, mb)              # our correct vetoes only
        r_t0 = lambda X: (X["waves"] & X["ideal"]) | (X["ideal"] & ~X["waves"] & X["t0"])  # + what a perfect t0 adds
        zt0 = Z(r_t0(A), r_t0(Bv), ms, mb)
        zid = Z(A["ideal"], Bv["ideal"], ms, mb)
        zor = Z(A["oracle"], Bv["oracle"], ms, mb)
        rung_t0 = Z(A["t0"], Bv["t0"], ms, mb)
        rows.append(dict(lab=lab, net=net, cor=cor, zt0=zt0, zid=zid, zor=zor, rung_t0=rung_t0, nB=int(mb.sum())))
        print(f"{lab:16s} where we are {net:+5.1f} | correct only {cor:+5.1f} (failures {cor - net:+4.1f}) | "
              f"+perfect t0 {zt0:+5.1f} ({zt0 - cor:+4.1f}) | +perfect jet time = ideal {zid:+5.1f} ({zid - zt0:+4.1f}) | "
              f"oracle {zor:+5.1f} (irreducible {zor - zid:+4.1f}) | rung: truth t0 + real jet times {rung_t0:+5.1f} | {mb.sum()} B MC")
    # ── the figures, in the style of the deck's m_jj-column plots (vbs_time_veto_common.column_figure):
    # absolute S/sqrt(B) at 3000 fb^-1 per column (left) and its change vs no veto (right), one dodged
    # marker per rung of the ladder, the case-stacked S + B MC count underneath
    R = C.setup_root(args.out)
    EDG = np.array([200, 1000, 2000, 4000.])
    ks = 3000e3 * C.NORM["vbf"]["xs_pb"] / C.NORM["vbf"]["sumw"]
    kb = 3000e3 * C.NORM["zjets"]["xs_pb"] / C.NORM["zjets"]["sumw"]
    nofail = lambda X: X["waves"] & X["ideal"]
    RUNGS = {"none": (np.zeros(len(ws), bool), np.zeros(len(wb), bool)),
             "waves": (V["S"]["waves"], V["B"]["waves"]),
             "nofail": (nofail(V["S"]), nofail(V["B"])),
             "t0": (V["S"]["t0"], V["B"]["t0"]),
             "ideal": (V["S"]["ideal"], V["B"]["ideal"]),
             "res": (V["S"]["res"], V["B"]["res"]),
             "eff": (V["S"]["eff"], V["B"]["eff"]),
             "oracle": (V["S"]["oracle"], V["B"]["oracle"])}
    STY = {"none":   ("No veto: where we were",          "#7F7F7F",    25),
           "waves":  ("WAVeS veto: where we are",        "kP10Yellow", 22),
           "nofail": ("WAVeS, no failures",              "kP10Yellow", 26),
           "t0":     ("+ perfect t_{0}",                 "kP10Brown",  34),
           "ideal":  ("+ perfect jet-time assignment",   "#1B7837",    20),
           "res":    ("+ perfect time resolution",       "kP10Blue",   23),
           "eff":    ("+ perfect HGTD efficiency",       "#7B3294",    32),
           "oracle": ("+ full acceptance: ceiling",      "#000000",    21)}
    absol, change, ex_abs, ex_chg = {}, {}, {}, {}
    for m, (vs, vb) in RUNGS.items():
        absol[m], change[m] = [], []
        for k in range(len(EDG) - 1):
            ms, mb = cs == k, cb == k
            S_ = ks * ws[ms & ~vs].sum(); dS = ks * np.sqrt((ws[ms & ~vs] ** 2).sum())
            B_ = kb * wb[mb & ~vb].sum(); dB = kb * np.sqrt((wb[mb & ~vb] ** 2).sum())
            z = S_ / np.sqrt(B_)
            absol[m].append((z, z * np.hypot(dS / S_, 0.5 * dB / B_), z * np.hypot(dS / S_, 0.5 * dB / B_), 0))
            es, ses, _, _ = C.eff(ws[ms], ~vs[ms]); eb, seb, _, _ = C.eff(wb[mb], ~vb[mb])
            r = es / np.sqrt(eb); er = r * np.hypot(ses / es, 0.5 * seb / eb)
            change[m].append((100 * (r - 1), 100 * er, 100 * er, 0))
        Sall, Ball = ks * ws[~vs].sum(), kb * wb[~vb].sum()
        ex_abs[m] = f"  ({Sall / np.sqrt(Ball):.0f})"
        ex_chg[m] = f"  ({100 * ((ws[~vs].sum() / ws.sum()) / np.sqrt(wb[~vb].sum() / wb.sum()) - 1):+.1f}%)"
    ceil = "ceiling (veto every pileup tagging jet): " + ", ".join(f"{v[0]:+.0f}%" for v in change["oracle"][:3]) \
        + f";  all events {ex_chg['oracle'].strip()[1:-1]}"
    # the S + B MC panel: each column's MC events by case
    def parts(e, col):
        n = nin(e)
        return [(lab, cc, [int(((n == nn) & (col == k)).sum()) for k in range(len(EDG) - 1)])
                for lab, cc, nn in [("Forward-Forward", "#3F8FFF", 2), ("Central-Forward", "#B3231C", 1), ("Other", "#CFCFCF", 0)]]
    sp, bp = parts(S, cs), parts(B, cb)
    scount = [sum(x[2][k] for x in sp) for k in range(3)]; bcount = [sum(x[2][k] for x in bp) for k in range(3)]
    label = "S: VBF H#rightarrowinv.,  B: Z(#rightarrowll)+jets"
    head = ["**All selected events: where we were, are, could be"]
    C.column_figure(args.out + "_yields", "S / #sqrt{B}  per m_{jj} column", absol, ex_abs, label + ",  3000 fb^{-1}", head,
                    bcount, "B MC", EDG, STY, ["none", "waves", "nofail", "t0", "ideal", "res", "eff", "oracle"], ref=None,
                    compact=True, count_parts=bp, sig_counts=(scount, sp))
    C.column_figure(args.out, "S/#sqrt{B} change vs no veto [%]", change, ex_chg, label, head + [ceil],
                    bcount, "B MC", EDG, STY, ["waves", "nofail", "t0", "ideal", "res", "eff"], ref=0.0, pct=True,
                    compact=True, count_parts=bp, sig_counts=(scount, sp))

elif args.what == "puretime":
    S_all = load(args.sig, "met_truth")[1]
    B_all = load(args.bkg, "z_pt")[1]
    pure = {"t0_truth": S_all["t0_truth"]}
    for L in "ab":
        pure.update({f"{L}_hs": S_all[f"{L}_hs"], f"{L}_core_n": S_all[f"{L}_pure_hs_n"],
                     f"{L}_core_t": S_all[f"{L}_pure_hs_t"], f"{L}_core_sig": S_all[f"{L}_pure_hs_sig"]})
    pcal = C.measure(pure, True, "pure HS legs")
    def series(e, lab, key):
        dt, pull, n, inh = [], [], [], []
        for L in "ab":
            m = e[f"{L}_{lab}"]
            nn = e[f"{L}_{key}_n"][m]
            d = (e[f"{L}_{key}_t"] - e["t0_truth"])[m]
            sg = pcal(nn) * e[f"{L}_{key}_sig"][m]
            dt.append(d); pull.append(np.where(nn >= 1, np.abs(d) / np.where(sg > 0, sg, np.inf), np.nan))
            n.append(nn); inh.append(e[f"{L}_hgtd_ntimed"][m] >= 1)
        return [np.concatenate(x) for x in (dt, pull, n, inh)]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    bins = np.linspace(-600, 600, 121)
    for e, lab, key, nm, col, fill in [(S_all, "hs", "pure_hs", "hard-scatter tagging jets (VBF): their hard-scatter tracks", "#2E9E4F", True),
                                       (B_all, "pu", "pure_pu", "pileup tagging jets ($Z$+jets): their own tracks", "#555555", False),
                                       (S_all, "pu", "pure_pu", "pileup tagging jets (VBF): their own tracks", "#999999", False)]:
        dt, pull, n, inh = series(e, lab, key)
        ok = n >= 1
        h, _ = np.histogram(np.clip(dt[ok], -599, 599), bins=bins)
        h = h / h.sum()
        ax.step(0.5 * (bins[1:] + bins[:-1]), h, where="mid", color=col, lw=2,
                label=f"{nm}\n  {100 * ok[inh].mean():.0f}% of those in HGTD timed; "
                      f"{100 * np.mean(pull[ok] >= args.nsigma):.0f}% beyond {args.nsigma:g}$\\sigma$ of $t_{{\\rm HS}}$")
        if fill: ax.fill_between(0.5 * (bins[1:] + bins[:-1]), h, step="mid", color=col, alpha=0.15)
        print(f"{nm}: timed (>= 1 own track) {100 * ok[inh].mean():.1f}% of those in HGTD; beyond {args.nsigma:g} sigma of t_HS {100 * np.mean(pull[ok] >= args.nsigma):.1f}%; "
              f"|dt| < 60 ps {100 * np.mean(np.abs(dt[ok]) < 60):.1f}%; median sigma {np.median((pcal(n) * 1.0)[ok]):.2f} x quoted")
    ax.set_yscale("log"); ax.set_ylim(1e-4, 40)
    ax.set_xlabel("perfect-assignment jet time $-\\ t_{\\rm HS}^{\\rm truth}$ [ps]   (overflow in the edge bins)")
    ax.set_ylabel("fraction of jets")
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    ax.grid(alpha=0.25)
    ax.text(1.0, 1.01, "all selected events   ATLAS Simulation Internal", transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color="#555555")
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
