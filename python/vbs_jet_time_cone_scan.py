#!/usr/bin/env python3
"""The jet-time cone scanned: util/vbs_time_veto run at --jet-time-dr = 0.1 ... 1.0, everything
else at the deck's settings (JVT + fJVT loose, jet time from >= --min-trk tracks, 2-sigma time
clustering, 3 sigma veto). The cone only restricts which of the jet's GHOST-associated tracks
enter its time, so beyond the extent of an R = 0.4 jet's ghost tracks the scan saturates.

Two figures:
  <out>_jets  the jet times (truth-matched tagging jets with >= 1 timed track in HGTD, all
              selected events, as python/vbs_jet_time_diag.py):
                coverage       HS jets (VBF) and pileup jets (Z+jets) with a jet time; the HS
                               core fraction (|t_jet - t_HS| < 60 ps, untimed = failure)
                contamination  of the timed HS jets: a mostly-pileup time cluster (< 50% of its
                               pT from the hard scatter), and |t_jet - t_HS| > 60 ps
                resolution     the HS core sigma (narrow width of a double-Gaussian fit)
  <out>_veto  the veto (recoil > 100 GeV, all events): the S/sqrt(B) change vs no veto per
              method; the signal and Z+jets kept (WAVeS); and the change against the default
              dR < 0.2, with a paired bootstrap error -- every cone sees the SAME events (the
              cone does not touch the selection), so the difference is far better measured
              than either point.
Each cone's jet-time sigma calibration f(n) is measured on its own signal HS legs, as in the deck.

    cd python && PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python vbs_jet_time_cone_scan.py --out ../figs/time_veto_k2/cone_scan
"""
import argparse
import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "../condor"))
ap.add_argument("--cones", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0")
ap.add_argument("--ref", type=float, default=0.2, help="the default cone (untagged trees); differences are against it")
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--nboot", type=int, default=200)
ap.add_argument("--out", required=True)
args = ap.parse_args()
C.ROOT = ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning
plt.rcParams.update({"font.size": 12})
CONES = [float(x) for x in args.cones.split(",")]
assert args.ref in CONES
METHODS = [("jj", "Jet time compatibility", "black", "o"), ("trkptz", "TRKPTZ $t_0$", "#C0392B", "s"),
           ("waves", "WAVeS $t_0$", "#E69F00", "^"), ("hgtd", "HGTD $t_0$ (Athena)", "#3B8ED0", "v"),
           ("truth", "truth $t_0$ (perfect)", "#8C5A3C", "P")]


def path(sample, dr):
    tag = "" if abs(dr - args.ref) < 1e-9 else f"_jtdr{dr:.1f}".replace(".", "p")
    return f"{args.dir}/{sample}/{sample}_jvtLoose{tag}_vbs_time_veto.root"


def load(p):
    t = uproot.open(p)
    dr = float(np.unique(t["meta"].arrays(["jet_time_dr"], library="np")["jet_time_dr"])[0])
    ev = t["events"]
    e = ev.arrays([k for k in ev.keys() if not k.startswith("jet_")], library="np")
    order = np.lexsort((e["entry"], e["file_path"]))      # one event order for every cone
    return {k: v[order] for k, v in e.items()}, dr


def vetoes(e, cal):
    out = {}
    sA, sB = cal(e["a_core_n"]) * e["a_core_sig"], cal(e["b_core_n"]) * e["b_core_sig"]
    both = (e["a_core_n"] >= args.min_trk) & (e["b_core_n"] >= args.min_trk)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["jj"] = both & (np.abs(e["a_core_t"] - e["b_core_t"]) / np.sqrt(sA ** 2 + sB ** 2) >= args.nsigma)
    for s, *_ in METHODS[1:]:
        out[s] = C.t0_veto(e, cal, s, args.min_trk, args.nsigma)[0]
    return out


def legs(e, lab):
    w = np.concatenate([e["weight"][e[f"{L}_{lab}"]] for L in "ab"]).astype(float)
    g = lambda k: np.concatenate([e[f"{L}_{k}"][e[f"{L}_{lab}"]] for L in "ab"])
    return dict(w=w, inh=g("hgtd_ntimed") >= 1, n=g("core_n"), hsf=g("core_hsfrac"),
                dt=g("core_t") - np.concatenate([e["t0_truth"][e[f"{L}_{lab}"]] for L in "ab"]))


def frac(w, num, den):
    e, _, lo, hi = C.eff(w[den], num[den])
    return 100 * e, 100 * lo, 100 * hi   # C.eff's lo / hi are distances below / above e


# ── per cone: the jet-time figures and the vetoes ──
J, V, KEYS = {}, {}, None
for dr in CONES:
    ps, pb = path("vbf", dr), path("zjets", dr)
    if not (os.path.exists(ps) and os.path.exists(pb)):
        print(f"dR < {dr:.1f}: missing {ps if not os.path.exists(ps) else pb}, skipped"); continue
    S_all, drs = load(ps); B_all, drb = load(pb)
    assert abs(drs - dr) < 1e-6 and abs(drb - dr) < 1e-6, f"{ps}: meta jet_time_dr {drs} / {drb}, expected {dr}"
    key = (S_all["file_path"], S_all["entry"], B_all["file_path"], B_all["entry"])
    if KEYS is None: KEYS = key
    assert all(len(a) == len(b) and np.array_equal(a, b) for a, b in zip(KEYS, key)), f"dR < {dr}: event lists differ"
    cal = C.measure(S_all, True, f"signal HS legs, dR < {dr:.1f}")
    # jet times (all selected events, as vbs_jet_time_diag)
    hs, pu = legs(S_all, "hs"), legs(B_all, "pu")
    tm = hs["n"] >= args.min_trk
    good = tm & (np.abs(hs["dt"]) < 60)
    J[dr] = dict(hs_timed=frac(hs["w"], tm, hs["inh"]), pu_timed=frac(pu["w"], pu["n"] >= args.min_trk, pu["inh"]),
                 hs_core=frac(hs["w"], good, hs["inh"]), hs_mostlypu=frac(hs["w"], hs["hsf"] < 0.5, hs["inh"] & tm),
                 hs_tail=frac(hs["w"], np.abs(hs["dt"]) > 60, hs["inh"] & tm),
                 hs_sigma=C.dgaus_core(hs["dt"][hs["inh"] & tm]))
    # vetoes (recoil > 100 GeV)
    ks, kb = S_all["met_truth"] > 100, B_all["z_pt"] > 100
    S = {k: v[ks] for k, v in S_all.items()}; B = {k: v[kb] for k, v in B_all.items()}
    V[dr] = dict(vS=vetoes(S, cal), vB=vetoes(B, cal))
    j = J[dr]
    print(f"dR < {dr:.1f}: HS timed {j['hs_timed'][0]:5.1f}%  PU timed {j['pu_timed'][0]:5.1f}%  HS core fraction {j['hs_core'][0]:5.1f}%  "
          f"timed HS mostly-pileup {j['hs_mostlypu'][0]:5.1f}%  beyond 60 ps {j['hs_tail'][0]:5.1f}%  core sigma {j['hs_sigma'][0]:5.1f} ps")
DRS = sorted(J)
assert args.ref in J, "the reference cone's trees are needed"
nin = lambda e: (e["a_hgtd_ntimed"] >= 1).astype(int) + (e["b_hgtd_ntimed"] >= 1).astype(int)
ws, wb = S["weight"].astype(float), B["weight"].astype(float)
CASES = {"all events": (np.ones(len(ws), bool), np.ones(len(wb), bool)),
         "Forward-Forward": (nin(S) == 2, nin(B) == 2), "Central-Forward": (nin(S) == 1, nin(B) == 1)}

# paired bootstrap: one set of Poisson-like resamplings of the VBF and Z+jets MC events, shared by every cone
rng = np.random.default_rng(20261022)
BS = [np.bincount(rng.integers(0, len(ws), len(ws)), minlength=len(ws)).astype(float) for _ in range(args.nboot)]
BB = [np.bincount(rng.integers(0, len(wb), len(wb)), minlength=len(wb)).astype(float) for _ in range(args.nboot)]


def zchg(vs, vb, ms, mb, bs=1.0, bb=1.0):
    a, b = ws * bs, wb * bb
    return 100 * ((a[ms & ~vs].sum() / a[ms].sum()) / np.sqrt(b[mb & ~vb].sum() / b[mb].sum()) - 1)


R = {}   # (case, method, dr) -> (change, error, diff vs ref, paired error, eps_S, eps_B, their errors)
for cname, (ms, mb) in CASES.items():
    for m, *_ in METHODS:
        ref_v = V[args.ref]
        ref_b = np.array([zchg(ref_v["vS"][m], ref_v["vB"][m], ms, mb, bs, bb) for bs, bb in zip(BS, BB)])
        for dr in DRS:
            vs, vb = V[dr]["vS"][m], V[dr]["vB"][m]
            z = zchg(vs, vb, ms, mb)
            zb = np.array([zchg(vs, vb, ms, mb, bs, bb) for bs, bb in zip(BS, BB)])
            es, ses, _, _ = C.eff(ws[ms], ~vs[ms]); eb, seb, _, _ = C.eff(wb[mb], ~vb[mb])
            R[cname, m, dr] = (z, np.std(zb), z - zchg(ref_v["vS"][m], ref_v["vB"][m], ms, mb), np.std(zb - ref_b),
                               100 * es, 100 * ses, 100 * eb, 100 * seb)
for cname in CASES:
    print(f"\n{cname}: S/sqrt(B) change vs no veto [%] (difference vs dR < {args.ref:.1f}, paired error)")
    for m, lab, *_ in METHODS:
        print(f"  {m:7s} " + "  ".join(f"{dr:.1f}: {R[cname, m, dr][0]:+5.1f} ({R[cname, m, dr][2]:+4.1f}+-{R[cname, m, dr][3]:.1f})" for dr in DRS))
    print("  WAVeS eps_S / eps_B: " + "  ".join(f"{dr:.1f}: {R[cname, 'waves', dr][4]:.1f}/{R[cname, 'waves', dr][6]:.1f}" for dr in DRS))

FOOT = "ATLAS Simulation Internal"
x = np.array(DRS)
SUB = f"jet time from $\\geq${args.min_trk} timed ghost tracks within $\\Delta R$ of the axis, 2$\\sigma$ time clustering"


def finish(fig, axes, foot):
    for a in axes:
        a.grid(alpha=0.3); a.axvline(args.ref, color="#888888", lw=0.9, ls=":")
    axes[-1].set_xlabel("jet-time cone $\\Delta R$")
    axes[-1].set_xticks(x); axes[-1].set_xlim(x[0] - 0.05, x[-1] + 0.05)
    fig.text(0.5, 0.035, foot, ha="center", va="top", fontsize=9.5, color="#555555")


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=150, bbox_inches="tight")
    print("wrote", stem + ".pdf / .png")


# ── figure 1: the jet times ──
fig, (a0, a1, a2) = plt.subplots(3, 1, figsize=(7.6, 9.4), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.25, 1, 0.8], hspace=0.07))
eb = lambda key: np.array([[J[d][key][1] for d in DRS], [J[d][key][2] for d in DRS]])
val = lambda key: np.array([J[d][key][0] for d in DRS])
for key, lab, col, mk, mfc in [("hs_timed", "hard-scatter jets (VBF) with a jet time", "#2E9E4F", "o", None),
                                ("pu_timed", "pileup jets ($Z$+jets) with a jet time", "#555555", "s", None),
                                ("hs_core", "hard-scatter core fraction (untimed = failure)", "#2E9E4F", "o", "white")]:
    a0.errorbar(x, val(key), yerr=eb(key), fmt=mk + "-", color=col, mfc=mfc or col, ms=6.5, capsize=2, lw=1.6, label=lab)
a0.set_ylabel("of the jets in HGTD [%]"); a0.set_ylim(0, 105)
a0.legend(fontsize=9.5, frameon=False, loc="lower right")
a0.set_title(f"truth-matched tagging jets with $\\geq$1 timed track in HGTD, all selected events\n{SUB}",
             fontsize=9.5, color="#444444", loc="left")
for key, lab, col, mk in [("hs_mostlypu", "time cluster mostly pileup ($<$50% of its $p_T$ hard scatter)", "#B3231C", "D"),
                          ("hs_tail", "$|t_{\\rm jet} - t_{\\rm HS}| > 60$ ps", "#E69F00", "v")]:
    a1.errorbar(x, val(key), yerr=eb(key), fmt=mk + "-", color=col, ms=6.5, capsize=2, lw=1.6, label=lab)
a1.set_ylabel("of the timed\nhard-scatter jets [%]"); a1.set_ylim(0, None)
a1.legend(fontsize=9.5, frameon=False, loc="lower right")
sg = np.array([J[d]["hs_sigma"][0] for d in DRS]); sge = np.array([J[d]["hs_sigma"][1] for d in DRS])
a2.errorbar(x, sg, yerr=sge, fmt="o-", color="#2E9E4F", ms=6.5, capsize=2, lw=1.6)
a2.set_ylabel("hard-scatter\ncore $\\sigma$ [ps]")
lo, hi = np.nanmin(sg - sge), np.nanmax(sg + sge); a2.set_ylim(lo - 0.25 * (hi - lo) - 0.5, hi + 0.25 * (hi - lo) + 0.5)
finish(fig, (a0, a1, a2), FOOT + "   (dotted: the default cone)")
save(fig, args.out + "_jets")

# ── figure 2: the veto ──
fig, (a0, a1, a2) = plt.subplots(3, 1, figsize=(7.6, 9.4), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.25, 0.85, 0.9], hspace=0.07))
n = len(METHODS)
for i, (m, lab, col, mk) in enumerate(METHODS):
    xx = x + 0.012 * (i - (n - 1) / 2)
    a0.errorbar(xx, [R["all events", m, d][0] for d in DRS], yerr=[R["all events", m, d][1] for d in DRS],
                fmt=mk + "-", color=col, ms=6.5, capsize=2, lw=1.5, label=lab)
a0.axhline(0, color="#999999", lw=0.8)
a0.set_ylabel("$S/\\sqrt{B}$ change\nvs no veto [%]")
a0.legend(fontsize=9.5, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), columnspacing=1.0, handletextpad=0.4)
for col_, lab, idx, mfc in [("#555555", "$Z$+jets kept", 6, None), ("#2E9E4F", "VBF kept", 4, "white")]:
    a1.errorbar(x, [R["all events", "waves", d][idx] for d in DRS], yerr=[R["all events", "waves", d][idx + 1] for d in DRS],
                fmt="o-", color=col_, mfc=mfc or col_, ms=6.5, capsize=2, lw=1.5, label=lab)
a1.set_ylabel("WAVeS veto:\nkept [%]"); a1.set_ylim(52, 100)
a1.legend(fontsize=9.5, frameon=False, loc="lower right", ncol=2)
for cname, ls, mfc, off in [("all events", "-", None, 0.0), ("Forward-Forward", "--", "white", -0.012),
                            ("Central-Forward", ":", "white", 0.012)]:
    d_ = np.array([R[cname, "waves", d][2] for d in DRS]); e_ = np.array([R[cname, "waves", d][3] for d in DRS])
    a2.errorbar(x + off, d_, yerr=e_, fmt="^", ls=ls, color="#E69F00", mfc=mfc or "#E69F00", ms=6.5, capsize=2, lw=1.5,
                label=f"WAVeS, {cname}")
a2.axhline(0, color="#999999", lw=0.8)
a2.set_ylabel(f"WAVeS: change vs\n$\\Delta R < {args.ref:.1f}$ [% points]")
a2.legend(fontsize=9, frameon=False, loc="lower right", ncol=1)
a2.text(0.02, 0.96, "paired bootstrap: the same events\nat every cone; each case's own events", transform=a2.transAxes,
        ha="left", va="top", fontsize=8.5, color="#555555")
finish(fig, (a0, a1, a2), "VBF $H\\to$inv. against $Z(\\to\\ell\\ell)$+jets, JVT + fJVT loose, recoil > 100 GeV, 3$\\sigma$, all events unless labelled\n"
       + SUB + "   " + FOOT)
save(fig, args.out + "_veto")
