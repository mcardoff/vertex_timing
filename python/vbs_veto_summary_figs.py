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
          its hard-scatter tracks, a pileup jet by its non-hard-scatter tracks; unlabelled jets untested). Needs the
          _pure trees.
  steps   the money plot's ladder taken apart (all events; same trees and definitions as money):
          <out>_waterfall  where the S/sqrt(B) gain comes from, step by step, with the ceiling
                           split into the hard-scatter false vetoes it undoes and the acceptance;
                           signal and Z+jets kept underneath
          <out>_pu_legs    every pileup tagging jet by the step that catches it
          <out>_failures   WAVeS's failures (vetoes the perfect veto would not make) by cause and pair
  chainevents  one typical event per step of the money plot's chain that the step before gets wrong
          (vetoes when it should keep, or keeps when it should veto) and this step fixes; writes
          <out>.csv in the candidates format of vbs_veto_event_cards.py / vbs_veto_displays.py
          (rpt_region_displays.py pick copies them out of the ntuples on the AF) and prints each
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
ap.add_argument("what", choices=["cases", "cf", "shares", "cone", "budget", "money", "steps", "chainevents", "puretime"])
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--sig2", help="cone: the second cone's signal tree")
ap.add_argument("--bkg2", help="cone: the second cone's background tree")
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--tag", default="JVT + fJVT loose")
ap.add_argument("--out", required=True)
ap.add_argument("--min-pu-pt", type=float, default=0.0, help="chainevents: ignore nothing by default")
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

elif args.what in ("money", "steps", "chainevents"):
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
    def jj_ideal(e):
        """jet time compatibility (jet 1 vs jet 2, no t0) with the perfect-assignment jet times: both legs
        labelled and timed (>= 1 track), vetoed when they differ by >= nsigma"""
        t, sg, ok = {}, {}, {}
        for L in "ab":
            hs, pu = e[f"{L}_hs"], e[f"{L}_pu"]
            n = np.where(hs, e[f"{L}_pure_hs_n"], np.where(pu, e[f"{L}_pure_pu_n"], 0))
            t[L] = np.where(hs, e[f"{L}_pure_hs_t"], e[f"{L}_pure_pu_t"])
            sg[L] = pcal(n) * np.where(hs, e[f"{L}_pure_hs_sig"], e[f"{L}_pure_pu_sig"])
            ok[L] = (hs | pu) & (n >= 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            pull = np.abs(t["a"] - t["b"]) / np.sqrt(sg["a"] ** 2 + sg["b"] ** 2)
        return ok["a"] & ok["b"] & (pull >= args.nsigma)
    oracle = lambda e: e["a_pu"] | e["b_pu"]
    V = {r: dict(waves=d["v" + r]["waves"], t0=d["v" + r]["truth"], ideal=ideal(d[r]), oracle=oracle(d[r]),
                 jjideal=jj_ideal(d[r])) for r in "SB"}
    # between the ideal veto and the oracle: perfect time resolution (every pileup tagging jet with a timed
    # track in HGTD caught), then perfect HGTD timing efficiency (every one with a track in HGTD's range,
    # timed or not); what is left is HGTD's acceptance (pileup jets with no track in 2.38 < |eta| < 4.0)
    for r in "SB":
        e = d[r]
        pu_with = lambda col: (e["a_pu"] & (e["a_" + col] >= 1)) | (e["b_pu"] & (e["b_" + col] >= 1))
        V[r]["res"] = V[r]["ideal"] | pu_with("pure_pu_n")
        V[r]["eff"] = V[r]["res"] | pu_with("hgtd_n")
    # THE CHAIN (the money plot since 2026-10-23): each step fixes one kind of mistake and builds on the
    # one before -- WAVeS -> the true t0 (real jet times) -> the out-of-time pileup tagging jets caught
    # (perfect-assignment time >= nsigma from t_HS: what better jet times could reach) -> the in-time ones
    # (timed, < nsigma: what needs better resolution) -> the untimed ones with a track in HGTD -> the false
    # vetoes undone (only events with a pileup tagging jet stay vetoed) -> every pileup tagging jet (the
    # ceiling: the rest is pileup HGTD cannot reach). Unlabelled jets keep their real-time test throughout.
    def pu_classes(e):
        o_, i_, u_ = (np.zeros(len(e["weight"]), bool) for _ in range(3))
        for L in "ab":
            n = e[f"{L}_pure_pu_n"]; sg = pcal(n) * e[f"{L}_pure_pu_sig"]
            with np.errstate(divide="ignore", invalid="ignore"):
                pull = np.abs(e[f"{L}_pure_pu_t"] - e["t0_truth"]) / sg
            pu = e[f"{L}_pu"]
            o_ |= pu & (n >= 1) & (pull >= args.nsigma)
            i_ |= pu & (n >= 1) & (pull < args.nsigma)
            u_ |= pu & (n < 1) & (e[f"{L}_hgtd_n"] >= 1)
        return o_, i_, u_
    for r, e in (("S", d["S"]), ("B", d["B"])):
        o_, i_, u_ = pu_classes(e)
        V[r]["oot"] = V[r]["t0"] | o_
        V[r]["itpu"] = V[r]["oot"] | i_
        V[r]["untimed"] = V[r]["itpu"] | u_
        V[r]["nofalse"] = V[r]["untimed"] & (e["a_pu"] | e["b_pu"])
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
        zjj = Z(A["jjideal"], Bv["jjideal"], ms, mb)
        chain = " -> ".join(f"{Z(A[k], Bv[k], ms, mb):+5.1f}" for k in ("waves", "t0", "oot", "itpu", "untimed", "nofalse", "oracle"))
        rows.append(dict(lab=lab, net=net, cor=cor, zt0=zt0, zid=zid, zor=zor, rung_t0=rung_t0, nB=int(mb.sum())))
        print(f"{lab:16s} where we are {net:+5.1f} | correct only {cor:+5.1f} (failures {cor - net:+4.1f}) | "
              f"+perfect t0 {zt0:+5.1f} ({zt0 - cor:+4.1f}) | +perfect jet time = ideal {zid:+5.1f} ({zid - zt0:+4.1f}) | "
              f"oracle {zor:+5.1f} (irreducible {zor - zid:+4.1f}) | rung: truth t0 + real jet times {rung_t0:+5.1f} | "
              f"jet 1 vs jet 2, perfect assignment {zjj:+5.1f} | chain {chain} | {mb.sum()} B MC")
    if args.what == "chainevents":
        import csv
        cal = C.measure(S_all, True, "signal HS legs")   # the veto's jet-time calibration, for the printout
        def legs_of(e, i):
            return [L for L in "ab"]
        def pu_info(e, L):
            n = e[f"{L}_pure_pu_n"]; sg = pcal(n) * e[f"{L}_pure_pu_sig"]
            with np.errstate(divide="ignore", invalid="ignore"):
                pull = np.abs(e[f"{L}_pure_pu_t"] - e["t0_truth"]) / sg
            return n, pull
        def typical(e, cand, keys):
            """the candidate nearest the middle of the candidates in each key (rank distance)"""
            idx = np.flatnonzero(cand)
            if len(idx) == 0: raise SystemExit("no candidate")
            dist = np.zeros(len(idx))
            for k in keys:
                v = k[idx]; r = np.argsort(np.argsort(v)) / max(len(v) - 1, 1)
                dist += (r - 0.5) ** 2
            return idx[np.argmin(dist)], len(idx)
        hs_pu = lambda e, L, O: e[f"{L}_pu"] & e[f"{O}_hs"]
        picks = []
        # + perfect t0: a hard-scatter pair WAVeS vetoes because its t0 is wrong, kept with the true t0
        e = S; X = V["S"]
        dt0 = np.abs(e["t0_waves"] - e["t0_truth"]); s0 = e["infl_waves"] * e["sig0_waves"]
        cand = X["waves"] & ~X["t0"] & e["a_hs"] & e["b_hs"] & (dt0 > args.nsigma * s0)
        i, n = typical(e, cand, [dt0, e["mjj"]]); picks.append(("chain_t0", "S", i, n, "WAVeS vetoes a hard-scatter pair: its t0 is wrong"))
        # the pileup steps (Z+jets, a hard-scatter + pileup pair): the step before keeps the event
        e = B; X = V["B"]
        for key, prev, sel, what in [("chain_oot", "t0", "oot", "an out-of-time pileup jet the jet time misses"),
                                     ("chain_itpu", "oot", "itpu", "an in-time pileup jet"),
                                     ("chain_untimed", "itpu", "untimed", "a pileup jet whose HGTD tracks got no time")]:
            best = None
            for L, O in (("a", "b"), ("b", "a")):
                n_, pull = pu_info(e, L)
                inside = (np.abs(e[f"{L}_eta"]) > 2.4) & (np.abs(e[f"{L}_eta"]) < 4.0)
                # the in-time example: a pileup jet inside HGTD timed by >= 2 of its tracks -- in time, not
                # merely within a wide one-track sigma
                c = {"oot": (n_ >= 1) & (pull >= args.nsigma), "itpu": (n_ >= 2) & (pull < args.nsigma) & inside,
                     "untimed": (n_ < 1) & (e[f"{L}_hgtd_n"] >= 1)}[sel]
                m = ~X[prev] & X[sel] & hs_pu(e, L, O) & c
                best = m if best is None else (best | m)
            # the pileup leg's pT and its perfect-time pull, for the typical pick
            pt_pu = np.where(e["a_pu"], e["a_pt"], e["b_pt"])
            i, n = typical(e, best, [pt_pu, e["mjj"]]); picks.append((key, "B", i, n, f"the step before keeps {what}"))
        # + false vetoes undone: a hard-scatter pair vetoed (the true t0, the real jet times): a jet timed by pileup tracks
        e = S; X = V["S"]
        pu_frac = lambda L: (e[f"{L}_core_hsfrac"] >= 0) & (e[f"{L}_core_hsfrac"] < 0.5) & (e[f"{L}_core_n"] >= args.min_trk)
        cand = X["untimed"] & ~X["nofalse"] & e["a_hs"] & e["b_hs"] & (pu_frac("a") | pu_frac("b"))
        print(f"VBF false vetoes at the untimed step: {100 * ws[X['untimed'] & ~X['nofalse']].sum() / ws.sum():.1f}% of events; "
              f"hard-scatter pairs with a pileup-timed jet {100 * ws[cand].sum() / ws.sum():.1f}%")
        i, n = typical(e, cand, [e["mjj"]]); picks.append(("chain_false", "S", i, n, "a hard-scatter pair vetoed: a jet timed by pileup tracks"))
        # the ceiling: a pileup jet past |eta| = 4, no track in HGTD's range, kept by every timing step
        e = B; X = V["B"]
        best = None
        for L, O in (("a", "b"), ("b", "a")):
            m = ~X["nofalse"] & X["oracle"] & hs_pu(e, L, O) & (e[f"{L}_hgtd_n"] < 1) & (np.abs(e[f"{L}_eta"]) > 4.0)
            best = m if best is None else (best | m)
        pt_pu = np.where(e["a_pu"], e["a_pt"], e["b_pt"])
        i, n = typical(e, best, [pt_pu, e["mjj"]]); picks.append(("chain_ceiling", "B", i, n, "a pileup jet past HGTD, kept by every timing step"))
        COLS = ["file_path", "entry", "a_idx", "b_idx", "a_hs", "a_pu", "b_hs", "b_pu", "t0_waves", "sig0_waves", "infl_waves",
                "t0_trkptz", "t0_hgtd", "ok_hgtd", "t0_truth", "mjj"]
        with open(args.out + ".csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["card", "sample"] + COLS); w.writeheader()
            for key, r, i, n, what in picks:
                e = S if r == "S" else B
                w.writerow(dict(card=key, sample="vbf" if r == "S" else "zjets", **{c: (e[c][i].decode() if isinstance(e[c][i], bytes) else e[c][i]) for c in COLS}))
                print(f"\n{key} ({'VBF' if r == 'S' else 'Z+jets'}, {n:,} candidates): {what}")
                print(f"  {e['file_path'][i]} #{e['entry'][i]}  m_jj {e['mjj'][i]:.0f} GeV  t_HS {e['t0_truth'][i]:+.0f} ps  WAVeS t0 {e['t0_waves'][i]:+.0f} ps "
                      f"(sigma {e['infl_waves'][i] * e['sig0_waves'][i]:.0f})")
                for L in "ab":
                    role = "HS" if e[f"{L}_hs"][i] else ("PU" if e[f"{L}_pu"][i] else "unlabelled")
                    sg = cal(np.array([e[f"{L}_core_n"][i]]))[0] * e[f"{L}_core_sig"][i]
                    pn, ppull = pu_info(e, L)
                    print(f"  jet {L} (idx {e[f'{L}_idx'][i]}): {role}, pT {e[f'{L}_pt'][i]:.0f} GeV, eta {e[f'{L}_eta'][i]:+.2f}; HGTD tracks {e[f'{L}_hgtd_n'][i]} "
                          f"({e[f'{L}_hgtd_ntimed'][i]} timed); jet time {e[f'{L}_core_t'][i]:+.0f} ps from {e[f'{L}_core_n'][i]} tracks (sigma {sg:.0f}, "
                          f"{100 * (1 - e[f'{L}_core_hsfrac'][i]):.0f}% PU)" + (f"; perfect PU time {e[f'{L}_pure_pu_t'][i]:+.0f} ps from {pn[i]} tracks, {ppull[i]:.1f} sigma from t_HS"
                          if e[f"{L}_pu"][i] and pn[i] >= 1 else ""))
        print("wrote", args.out + ".csv")
        raise SystemExit(0)
    if args.what == "steps":
        P = lambda w, m: 100 * w[m].sum() / w.sum()
        Zc = lambda vs, vb: 100 * ((ws[~vs].sum() / ws.sum()) / np.sqrt(wb[~vb].sum() / wb.sum()) - 1)
        cal = C.measure(S_all, True, "signal HS legs")
        def leg_veto(e, L, lab_mask, pure_key=None):
            """leg L vetoed against the true t0: by its perfect-assignment time (pure_key) or as used"""
            if pure_key:
                n = e[f"{L}_{pure_key}_n"]; sg = pcal(n) * e[f"{L}_{pure_key}_sig"]; t = e[f"{L}_{pure_key}_t"]; ok = n >= 1
            else:
                n = e[f"{L}_core_n"]; sg = cal(n) * e[f"{L}_core_sig"]; t = e[f"{L}_core_t"]; ok = n >= args.min_trk
            with np.errstate(divide="ignore", invalid="ignore"):
                return lab_mask & ok & (np.abs(t - e["t0_truth"]) / sg >= args.nsigma)
        def variant(e, hs_pure, pu_pure, unl):
            v = np.zeros(len(e["weight"]), bool)
            for L in "ab":
                hs, pu = e[f"{L}_hs"], e[f"{L}_pu"]
                v |= leg_veto(e, L, hs, "pure_hs" if hs_pure else None)
                v |= leg_veto(e, L, pu, "pure_pu" if pu_pure else None)
                if unl: v |= leg_veto(e, L, ~hs & ~pu)
            return v
        has_pu = lambda e: e["a_pu"] | e["b_pu"]
        # the money plot's chain: (label, colour, rung it ends, S vetoes, B vetoes)
        CHAIN = [("WAVeS\nveto", "#FFA90E", "waves", V["S"]["waves"], V["B"]["waves"]),
                 ("+ perfect\n$t_0$", "#A96B59", "t0", V["S"]["t0"], V["B"]["t0"]),
                 ("+ out-of-\ntime pileup\ncaught", "#1B7837", "oot", V["S"]["oot"], V["B"]["oot"]),
                 ("+ in-time\npileup\ncaught", "#3F90DA", "itpu", V["S"]["itpu"], V["B"]["itpu"]),
                 ("+ untimed\npileup in\nHGTD caught", "#7B3294", "untimed", V["S"]["untimed"], V["B"]["untimed"]),
                 ("+ false\nvetoes\nundone", "#717581", "nofalse", V["S"]["nofalse"], V["B"]["nofalse"]),
                 ("+ pileup\nHGTD cannot\nreach:\nceiling", "#000000", "oracle", V["S"]["oracle"], V["B"]["oracle"])]
        z = [0.0] + [Zc(vs, vb) for *_, vs, vb in CHAIN]
        kS = [100.0] + [100 - P(ws, vs) for *_, vs, vb in CHAIN]; kB = [100.0] + [100 - P(wb, vb) for *_, vs, vb in CHAIN]
        for (lab, *_), z0, z1, s_, b_ in zip(CHAIN, z[:-1], z[1:], kS[1:], kB[1:]):
            print(f"{lab.replace(chr(10), ' '):34s} {z0:+5.1f} -> {z1:+5.1f}  ({z1 - z0:+5.1f})   VBF kept {s_:5.1f}%  Z+jets kept {b_:5.1f}%")
        # figure 1: the waterfall
        fig, (a0, a1) = plt.subplots(2, 1, figsize=(12, 7.6), sharex=True, gridspec_kw=dict(height_ratios=[1.5, 1], hspace=0.06))
        xs = np.arange(len(CHAIN))
        for i, ((lab, col, rung, *_), z0, z1) in enumerate(zip(CHAIN, z[:-1], z[1:])):
            a0.bar(i, z1 - z0, bottom=z0, color=col, width=0.62, edgecolor="black" if col != "#000000" else "#000000", lw=0.6,
                   alpha=0.95 if rung else 0.7, hatch=None if rung else "//")
            dz = round(z1, 1) - round(z0, 1)   # the difference of the rounded totals: the labels add up to the last one
            a0.text(i, max(z0, z1) + 0.6, f"{dz:+.1f}" + (f"\n= {z1:+.1f}%" if rung else ""), ha="center", va="bottom",
                    fontsize=12, fontweight="bold", linespacing=1.15)
            if rung: a0.plot([i - 0.42, i + 0.42], [z1, z1], color="black", lw=1.4)
            if i + 1 < len(CHAIN): a0.plot([i + 0.31, i + 0.69], [z1, z1], color="#888888", lw=0.8, ls=":")
        a0.axhline(0, color="#888888", lw=0.8)
        a0.set_ylabel("$S/\\sqrt{B}$ change vs no veto [%]"); a0.set_ylim(min(z) - 2, max(z) + 7)
        a0.text(0.01, 0.97, "each step builds on the one before; the line: where the money plot's marker sits\n"
                "the side branch (WAVeS, no wrong-reason vetoes) is not in this chain", transform=a0.transAxes, va="top", fontsize=9.5, color="#444444")
        a0.grid(alpha=0.3, axis="y")
        a1.plot(np.arange(-1, len(CHAIN)), kB, "s-", color="#555555", ms=7, lw=1.6, label="$Z$+jets kept")
        a1.plot(np.arange(-1, len(CHAIN)), kS, "o-", color="#2E9E4F", mfc="white", ms=7, lw=1.6, label="VBF kept")
        for i in range(len(CHAIN)):
            a1.text(i, kS[i + 1] + 2.5, f"{kS[i + 1]:.0f}", ha="center", fontsize=10.5, color="#2E9E4F")
            a1.text(i, kB[i + 1] - 8.5, f"{kB[i + 1]:.0f}", ha="center", fontsize=10.5, color="#555555")
        a1.set_ylabel("kept [%]"); a1.set_ylim(20, 108); a1.legend(fontsize=10, frameon=False, loc="lower left"); a1.grid(alpha=0.3, axis="y")
        a1.set_xticks(np.arange(-1, len(CHAIN))); a1.set_xticklabels(["no\nveto"] + [c[0] for c in CHAIN], fontsize=11)
        a1.set_xlim(-1.4, len(CHAIN) - 0.5)
        fig.text(0.5, -0.07, "all selected events, recoil > 100 GeV, JVT + fJVT loose; from the perfect $t_0$ on, $t_0 = t_{\\rm HS}$; a pileup tagging jet's time: "
                 "its non-hard-scatter tracks (perfect assignment)   ATLAS Simulation Internal", ha="center", fontsize=9, color="#555555")
        args.out, base = args.out + "_waterfall", args.out
        save(fig); args.out = base
        # figure 2: every pileup tagging jet by the step that catches it
        fig, ax = plt.subplots(figsize=(11, 3.4))
        CATS = [("out of time ($\\geq 3\\sigma$):\nthe out-of-time step", "#1B7837"), ("in time ($< 3\\sigma$):\nthe in-time step", "#3F90DA"),
                ("tracks in HGTD, none\ntimed: the untimed step", "#7B3294"), ("no track in HGTD,\naxis $|\\eta| < 2.4$", "#555555"),
                ("no track in HGTD,\naxis $|\\eta| > 2.4$", "#AAAAAA")]
        for row, (r, e, nm) in enumerate((("B", B, "$Z$+jets"), ("S", S, "VBF $H\\to$inv."))):
            w = np.concatenate([e["weight"][e[f"{L}_pu"]] for L in "ab"]).astype(float)
            g = lambda k: np.concatenate([e[f"{L}_{k}"][e[f"{L}_pu"]] for L in "ab"])
            n, hn, eta = g("pure_pu_n"), g("hgtd_n"), np.abs(g("eta"))
            with np.errstate(divide="ignore", invalid="ignore"):
                pull = np.abs(g("pure_pu_t") - np.concatenate([e["t0_truth"][e[f"{L}_pu"]] for L in "ab"])) / (pcal(n) * g("pure_pu_sig"))
            masks = [(n >= 1) & (pull >= 3), (n >= 1) & (pull < 3), (n < 1) & (hn >= 1), (n < 1) & (hn < 1) & (eta < 2.4), (n < 1) & (hn < 1) & (eta >= 2.4)]
            left = 0.0
            for (lab, col), m in zip(CATS, masks):
                f = 100 * w[m].sum() / w.sum()
                ax.barh(row, f, left=left, color=col, height=0.6, label=lab if row == 0 else None, edgecolor="white")
                if f > 3: ax.text(left + f / 2, row, f"{f:.0f}%", ha="center", va="center", color="white", fontsize=11, fontweight="bold")
                left += f
            beyond = (n < 1) & (hn < 1) & (eta >= 2.4)
            print(f"{nm}: pileup legs " + ", ".join(f"{100 * w[m].sum() / w.sum():.1f}%" for m in masks)
                  + f";  of 'no track, |eta| > 2.4': axis |eta| > 4.0 {100 * w[beyond & (eta > 4.0)].sum() / max(w[beyond].sum(), 1e-12):.0f}%")
        ax.set_yticks([0, 1]); ax.set_yticklabels(["$Z(\\to\\ell\\ell)$+jets", "VBF $H\\to$inv."]); ax.set_xlim(0, 100)
        ax.set_xlabel("pileup tagging jets [%]")
        ax.legend(fontsize=9, frameon=False, ncol=5, loc="lower center", bbox_to_anchor=(0.5, 1.0))
        fig.text(0.5, -0.08, "every truth-pileup tagging jet, all selected events: what its perfect-assignment time (non-hard-scatter tracks) does against $t_{\\rm HS}$   "
                 "ATLAS Simulation Internal", ha="center", fontsize=9, color="#555555")
        args.out = base + "_pu_legs"; save(fig); args.out = base
        # figure 2b: in-time pileup -- the pileup jets' perfect-assignment times against t_HS (generator-weighted),
        # each with a Gaussian fit; the expectation for two independent collisions, from a Gaussian fit to the
        # true t_HS spread; and the 3-sigma test window, sigma from the double-Gaussian core width of the
        # hard-scatter jets' perfect-assignment times (the project's core-resolution recipe)
        import ROOT as _R
        _R.gROOT.SetBatch(True); _R.gErrorIgnoreLevel = _R.kWarning; C.ROOT = _R
        def gaus_fit(x, w, lo, hi, nb):
            h = _R.TH1D(f"h_gfit_{len(C.KEEP)}", "", nb, lo, hi); h.Sumw2()
            for xi, wi in zip(x, w): h.Fill(xi, wi)
            f = _R.TF1(f"f_gfit_{len(C.KEEP)}", "gaus", lo, hi); f.SetParameters(h.GetMaximum(), 0.0, h.GetStdDev())
            h.Fit(f, "RQN0"); C.KEEP.extend([h, f])
            return f.GetParameter(1), abs(f.GetParameter(2)), f.GetParError(2)
        hsdt = np.concatenate([(S[f"{L}_pure_hs_t"] - S["t0_truth"])[S[f"{L}_hs"] & (S[f"{L}_pure_hs_n"] >= 1)] for L in "ab"])
        s_res, s_res_e = C.dgaus_core(hsdt)
        _, s_hs, s_hs_e = gaus_fit(S["t0_truth"], ws, -800, 800, 160)
        fig, ax = plt.subplots(figsize=(8.6, 5.2))
        bins = np.linspace(-800, 800, 41); xc = 0.5 * (bins[1:] + bins[:-1])
        for r, e, w_ev, nm, col in (("B", B, wb, "pileup tagging jets ($Z$+jets)", "#555555"), ("S", S, ws, "pileup tagging jets (VBF)", "#AAAAAA")):
            g = lambda k: np.concatenate([e[f"{L}_{k}"][e[f"{L}_pu"]] for L in "ab"])
            w = np.concatenate([w_ev[e[f"{L}_pu"]] for L in "ab"]).astype(float)
            n = g("pure_pu_n"); m = n >= 1
            dt = (g("pure_pu_t") - np.concatenate([e["t0_truth"][e[f"{L}_pu"]] for L in "ab"]))[m]
            sg = (pcal(n) * g("pure_pu_sig"))[m]; wm = w[m]
            in_test = 100 * wm[np.abs(dt) < 3 * sg].sum() / wm.sum()
            in_band = 100 * wm[np.abs(dt) < 3 * s_res].sum() / wm.sum()
            neff = wm.sum() ** 2 / (wm ** 2).sum()
            mu, sfit, sfit_e = gaus_fit(dt, wm, -800, 800, 80)
            h, _ = np.histogram(dt, bins=bins, weights=wm); h = h / wm.sum() / (bins[1] - bins[0]) * 100
            ax.step(xc, h, where="mid", color=col, lw=2,
                    label=f"{nm}: {m.sum():,} MC jets ($N_{{\\rm eff}}$ = {neff:,.0f}); Gaussian fit $\\sigma$ = {sfit:.0f} $\\pm$ {sfit_e:.0f} ps")
            ax.plot(xc, 100 * np.exp(-0.5 * ((xc - mu) / sfit) ** 2) / (np.sqrt(2 * np.pi) * sfit), color=col, lw=1.0, ls=":")
            print(f"{nm}: {m.sum()} MC jets, N_eff {neff:.0f}; Gaussian fit mean {mu:.1f} ps, sigma {sfit:.1f} +- {sfit_e:.1f} ps; "
                  f"within the per-jet 3-sigma test {in_test:.1f}%, within +-3 x {s_res:.1f} ps {in_band:.1f}%")
        sp = np.sqrt(2) * s_hs
        ax.plot(xc, 100 * np.exp(-0.5 * (xc / sp) ** 2) / (np.sqrt(2 * np.pi) * sp), color="#B3231C", lw=1.4, ls="--",
                label=f"two independent collisions: $\\sqrt{{2}} \\times$ {s_hs:.0f} ps = {sp:.0f} ps ($t_{{\\rm HS}}$ spread, Gaussian fit)")
        ax.axvspan(-3 * s_res, 3 * s_res, color="#3F90DA", alpha=0.18, lw=0,
                   label=f"$\\pm 3\\sigma$, $\\sigma$ = {s_res:.1f} ps (hard-scatter perfect-assignment times, double-Gaussian core)")
        print(f"t_HS spread (VBF, Gaussian fit) {s_hs:.1f} +- {s_hs_e:.1f} ps -> sqrt(2) x = {sp:.1f} ps; HS perfect-assignment core sigma {s_res:.1f} +- {s_res_e:.1f} ps")
        ax.set_xlim(-800, 800); ax.set_ylim(0, 0.34)
        ax.set_xlabel("perfect-assignment pileup-jet time $-\\ t_{\\rm HS}^{\\rm truth}$ [ps]")
        ax.set_ylabel("pileup jets [% per ps]")
        ax.legend(fontsize=8.6, frameon=False, loc="upper left")
        ax.grid(alpha=0.25)
        ax.text(1.0, 1.01, "all selected events, generator-weighted; dotted: the Gaussian fits   ATLAS Simulation Internal",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color="#555555")
        args.out = base + "_intime"; save(fig); args.out = base
        # figure 3: WAVeS's failures
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.2), sharey=True, gridspec_kw=dict(wspace=0.08))
        for row, (r, e, w, nm) in enumerate((("B", B, wb, "$Z(\\to\\ell\\ell)$+jets"), ("S", S, ws, "VBF $H\\to$inv."))):
            X = V[r]; f = X["waves"] & ~X["ideal"]
            hs2 = e["a_hs"] & e["b_hs"]; pu = has_pu(e); un = ~hs2 & ~pu
            for ax, parts in ((axes[0], [("wrong $t_0$ (the true $t_0$ keeps it)", "#A96B59", f & ~X["t0"]),
                                         ("jet time (the true $t_0$ still vetoes it)", "#1B7837", f & X["t0"])]),
                              (axes[1], [("both jets hard scatter", "#2E9E4F", f & hs2), ("a pileup jet, in time or untimed", "#555555", f & pu),
                                         ("an unlabelled jet", "#BBBBBB", f & un)])):
                left = 0.0
                for lab, col, m in parts:
                    v = P(w, m)
                    ax.barh(row, v, left=left, color=col, height=0.6, label=lab if row == 0 else None, edgecolor="white")
                    if v > 0.85: ax.text(left + v / 2, row, f"{v:.1f}", ha="center", va="center", color="white", fontsize=10.5, fontweight="bold")
                    left += v
            share = 100 * w[f].sum() / w[X["waves"]].sum()
            axes[1].text(left + 0.2, row, f"{P(w, f):.1f}% of the events\n= {share:.0f}% of WAVeS's vetoes", va="center", fontsize=9.5)
        axes[0].set_yticks([0, 1]); axes[0].set_yticklabels(["$Z(\\to\\ell\\ell)$+jets", "VBF $H\\to$inv."])
        axes[0].set_title("by cause", fontsize=11); axes[1].set_title("by the tagging pair", fontsize=11)
        for ax in axes:
            ax.set_xlabel("wrong-reason vetoes [% of all events]"); ax.set_xlim(0, 10.5)
            ax.legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=1)
        fig.text(0.5, 1.02, "a wrong-reason veto: WAVeS vetoes an event that the perfect veto (true $t_{\\rm HS}$, perfect jet-time assignment) keeps   "
                 "ATLAS Simulation Internal", ha="center", fontsize=9.5, color="#555555")
        args.out = base + "_failures"; save(fig); args.out = base
        raise SystemExit(0)
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
             "oot": (V["S"]["oot"], V["B"]["oot"]),
             "itpu": (V["S"]["itpu"], V["B"]["itpu"]),
             "untimed": (V["S"]["untimed"], V["B"]["untimed"]),
             "nofalse": (V["S"]["nofalse"], V["B"]["nofalse"]),
             "oracle": (V["S"]["oracle"], V["B"]["oracle"])}
    STY = {"none":   ("No veto: where we were",          "#7F7F7F",    25),
           "waves":  ("WAVeS veto: where we are",        "kP10Yellow", 22),
           "nofail": ("WAVeS, no wrong-reason vetoes", "kP10Yellow", 26),
           "t0":     ("+ perfect t_{0}",                 "kP10Brown",  34),
           "oot":    ("+ out-of-time pileup caught",     "#1B7837",    20),
           "itpu":   ("+ in-time pileup caught",         "kP10Blue",   23),
           "untimed": ("+ untimed pileup in HGTD caught", "#7B3294",   32),
           "nofalse": ("+ false vetoes undone",          "#717581",    33),
           "oracle": ("+ pileup HGTD cannot reach: ceiling", "#000000", 21)}
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
                    bcount, "B MC", EDG, STY, ["none", "waves", "nofail", "t0", "oot", "itpu", "untimed", "nofalse", "oracle"], ref=None,
                    compact=True, count_parts=bp, sig_counts=(scount, sp))
    C.column_figure(args.out, "S/#sqrt{B} change vs no veto [%]", change, ex_chg, label, head + [ceil],
                    bcount, "B MC", EDG, STY, ["waves", "nofail", "t0", "oot", "itpu", "untimed", "nofalse"], ref=0.0, pct=True,
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
    offs = []
    for e, lab, key, nm, col, fill in [(S_all, "hs", "pure_hs", "hard-scatter tagging jets (VBF): their hard-scatter tracks", "#2E9E4F", True),
                                       (B_all, "pu", "pure_pu", "pileup tagging jets ($Z$+jets): their non-hard-scatter tracks", "#555555", False),
                                       (S_all, "pu", "pure_pu", "pileup tagging jets (VBF): their non-hard-scatter tracks", "#999999", False)]:
        dt, pull, n, inh = series(e, lab, key)
        ok = n >= 1
        h, _ = np.histogram(dt[ok], bins=bins)   # no overflow bins: entries beyond the range go off the plot
        h = h / ok.sum()                         # ... but stay in the normalisation (a fraction of the timed jets)
        offs.append(100 * np.mean(np.abs(dt[ok]) >= bins[-1]))
        ax.step(0.5 * (bins[1:] + bins[:-1]), h, where="mid", color=col, lw=2,
                label=f"{nm}\n  {100 * ok[inh].mean():.0f}% of those in HGTD timed; "
                      f"{100 * np.mean(pull[ok] >= args.nsigma):.0f}% beyond {args.nsigma:g}$\\sigma$ of $t_{{\\rm HS}}$")
        if fill: ax.fill_between(0.5 * (bins[1:] + bins[:-1]), h, step="mid", color=col, alpha=0.15)
        print(f"{nm}: timed (>= 1 own track) {100 * ok[inh].mean():.1f}% of those in HGTD; beyond {args.nsigma:g} sigma of t_HS {100 * np.mean(pull[ok] >= args.nsigma):.1f}%; "
              f"|dt| < 60 ps {100 * np.mean(np.abs(dt[ok]) < 60):.1f}%; median sigma {np.median((pcal(n) * 1.0)[ok]):.2f} x quoted")
    ax.set_yscale("log"); ax.set_ylim(1e-4, 40)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("perfect-assignment jet time $-\\ t_{\\rm HS}^{\\rm truth}$ [ps]   (beyond $\\pm$600 ps, off scale: "
                  + ", ".join(f"{o:.1f}%" for o in offs) + " of the jets)")
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
