#!/usr/bin/env python3
"""How good are the tagging-jet times? An efficiency / assignment / resolution diagnostic of the
jet time the timing veto uses (util/vbs_time_veto core_: timed ghost tracks within dR < 0.2,
2-sigma time clustering, highest-sum-pT cluster), against the jet's axis |eta|, its pT and the
number of its ghost tracks HGTD times. Truth-matched hard-scatter tagging jets of VBF H->inv.
(wrong time -> false veto) and truth-matched pileup tagging jets of Z+jets and VBF (no time, or
the hard-scatter time -> missed), all selected events. Needs the --pure-times trees.

  top     efficiency: of the jets with >= 1 timed track in HGTD, the share with a jet time from
          >= --min-trk tracks (the veto's 'timed')
  middle  assignment purity: of the timed jets, the share whose time cluster is mostly of the
          right origin (HS jets: >= 50% of its pT from the hard scatter; pileup jets: < 50%)
  bottom  resolution of the HS jets: core sigma (narrow width of a double-Gaussian fit,
          FIT_DBLGAUS's recipe) of t_jet - t_HS, as used (filled) and with a perfect track
          assignment (open: the pure_hs time, HS tracks only)

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_jet_time_diag.py \\
        --sig condor/vbf/vbf_jvtLoose_pure_vbs_time_veto.root --bkg condor/zjets/zjets_jvtLoose_pure_vbs_time_veto.root \\
        --out figs/time_veto_k2/jet_time_diag
"""
import argparse
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--sig2", help="the same signal tree at another --jet-time-dr (e.g. 0.4): its as-used curves added, dashed")
ap.add_argument("--bkg2", help="the same background tree at that cone")
ap.add_argument("--out", required=True)
args = ap.parse_args()
C.ROOT = ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning
plt.rcParams.update({"font.size": 12})

COLS = ["weight", "t0_truth"] + [f"{L}_{k}" for L in "ab" for k in
        ["eta", "pt", "hs", "pu", "hgtd_ntimed", "core_n", "core_t", "core_hsfrac", "pure_hs_n", "pure_hs_t"]]


def legs(path, lab):
    """the tagging legs carrying truth label `lab` ('hs' / 'pu'), both legs stacked"""
    t = uproot.open(path)["events"]
    cols = [c for c in COLS if c in t.keys()]
    e = t.arrays(cols, library="np")
    out = {}
    for k in ["eta", "pt", "hgtd_ntimed", "core_n", "core_t", "core_hsfrac", "pure_hs_n", "pure_hs_t"]:
        if f"a_{k}" not in e: continue
        out[k] = np.concatenate([e[f"{L}_{k}"][e[f"{L}_{lab}"]] for L in "ab"])
    out["t0"] = np.concatenate([e["t0_truth"][e[f"{L}_{lab}"]] for L in "ab"])
    out["w"] = np.concatenate([e["weight"][e[f"{L}_{lab}"]] for L in "ab"]).astype(float)
    out["aeta"] = np.abs(out["eta"])
    return out


def cone(path):
    return float(np.unique(uproot.open(path)["meta"].arrays(["jet_time_dr"], library="np")["jet_time_dr"])[0])
HS = legs(args.sig, "hs")
PU_Z = legs(args.bkg, "pu")
DR1 = cone(args.sig)
ALT = None
if args.sig2:
    ALT = (legs(args.sig2, "hs"), legs(args.bkg2, "pu"), cone(args.sig2))
XV = [("aeta", "tagging-jet axis $|\\eta|$", np.array([2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4])),
      ("pt", "tagging-jet $p_T$ [GeV]", np.array([30, 40, 60, 100, 200, 1000.])),
      ("hgtd_ntimed", "timed ghost tracks in HGTD", np.array([1, 2, 3, 5, 9, 15, 100]))]


def frac(w, num, den):
    e, err, lo, hi = C.eff(w[den], num[den]) if den.sum() else (np.nan,) * 4
    return e, lo, hi


fig, axes = plt.subplots(3, 3, figsize=(15, 11), sharex="col", gridspec_kw=dict(hspace=0.07, wspace=0.18))
CAT = {"pt", "hgtd_ntimed"}   # categorical x (bins as equal-width columns)
for c, (key, xlab, edges) in enumerate(XV):
    xc = np.arange(len(edges) - 1) if key in CAT else 0.5 * (edges[:-1] + edges[1:])
    def inbin(pop, b):
        return (pop[key] >= edges[b]) & (pop[key] < edges[b + 1])
    # row 0: efficiency, HS and pileup jets
    for pop, nm, col, mk in [(HS, "hard-scatter jets (VBF)", "#2E9E4F", "o"), (PU_Z, "pileup jets ($Z$+jets)", "#555555", "s")]:
        inh, timed = pop["hgtd_ntimed"] >= 1, pop["core_n"] >= args.min_trk
        v = np.array([frac(pop["w"], timed, inbin(pop, b) & inh) for b in range(len(edges) - 1)])
        axes[0, c].errorbar(xc, v[:, 0], yerr=[v[:, 1], v[:, 2]], fmt=mk + "-", color=col, ms=7, lw=1.6, capsize=3,
                            label=nm if c == 0 else None)
    if ALT:
        for pop, col, mk in [(ALT[0], "#2E9E4F", "o"), (ALT[1], "#555555", "s")]:
            inh2, timed2 = pop["hgtd_ntimed"] >= 1, pop["core_n"] >= args.min_trk
            v = np.array([frac(pop["w"], timed2, inbin(pop, b) & inh2) for b in range(len(edges) - 1)])
            axes[0, c].errorbar(xc, v[:, 0], yerr=[v[:, 1], v[:, 2]], fmt=mk + "--", color=col, ms=5, lw=1.2, capsize=2, alpha=0.8,
                                label=(f"  ... jet-time cone $\\Delta R < {ALT[2]:.1f}$" if pop is ALT[0] else None) if c == 0 else None)
    # the HS jets with a perfect track assignment: >= min-trk timed HS tracks (any ghost track in HGTD)
    inh, timed = HS["hgtd_ntimed"] >= 1, HS["pure_hs_n"] >= args.min_trk
    v = np.array([frac(HS["w"], timed, inbin(HS, b) & inh) for b in range(len(edges) - 1)])
    axes[0, c].errorbar(xc, v[:, 0], yerr=[v[:, 1], v[:, 2]], fmt="o-", color="#2E9E4F", ms=7, lw=1.6, capsize=3, mfc="white",
                        label="hard-scatter jets, perfect track assignment" if c == 0 else None)
    # row 1: HS core fraction (|t_jet - t_HS| < 60 ps of the jets in HGTD; an untimed jet counts as a failure),
    # as used and with a perfect track assignment; the assignment purity of the timed ones dashed
    inh = HS["hgtd_ntimed"] >= 1
    for tcol, ncol, fill, lab in [("core_t", "core_n", True, "as used"), ("pure_hs_t", "pure_hs_n", False, "perfect track assignment")]:
        ok = (HS[ncol] >= args.min_trk) & (np.abs(HS[tcol] - HS["t0"]) < 60)
        v = np.array([frac(HS["w"], ok, inbin(HS, b) & inh) for b in range(len(edges) - 1)])
        axes[1, c].errorbar(xc, v[:, 0], yerr=[v[:, 1], v[:, 2]], fmt="o-", color="#2E9E4F", ms=7, lw=1.6, capsize=3,
                            mfc="#2E9E4F" if fill else "white", label=("core fraction, " + lab) if c == 0 else None)
    if ALT:
        H2 = ALT[0]; inh2 = H2["hgtd_ntimed"] >= 1
        ok = (H2["core_n"] >= args.min_trk) & (np.abs(H2["core_t"] - H2["t0"]) < 60)
        v = np.array([frac(H2["w"], ok, (H2[key] >= edges[b]) & (H2[key] < edges[b + 1]) & inh2) for b in range(len(edges) - 1)])
        axes[1, c].errorbar(xc, v[:, 0], yerr=[v[:, 1], v[:, 2]], fmt="o--", color="#2E9E4F", ms=5, lw=1.2, capsize=2, alpha=0.8,
                            label=f"core fraction, jet-time cone $\\Delta R < {ALT[2]:.1f}$" if c == 0 else None)
    timed = HS["core_n"] >= args.min_trk
    v = np.array([frac(HS["w"], HS["core_hsfrac"] >= 0.5, inbin(HS, b) & timed) for b in range(len(edges) - 1)])
    axes[1, c].errorbar(xc, v[:, 0], yerr=[v[:, 1], v[:, 2]], fmt="d--", color="#7B3294", ms=6, lw=1.3, capsize=3,
                        label="assignment purity: timed by mostly HS tracks" if c == 0 else None)
    # row 2: HS resolution, as used and with a perfect track assignment
    for tcol, ncol, fill, lab in [("core_t", "core_n", True, "as used (the veto's jet time)"),
                                  ("pure_hs_t", "pure_hs_n", False, "perfect track assignment (HS tracks only)")]:
        sig = []
        for b in range(len(edges) - 1):
            m = inbin(HS, b) & (HS[ncol] >= args.min_trk)
            sig.append(C.dgaus_core(HS[tcol][m] - HS["t0"][m])[0] if m.sum() >= 300 else np.nan)
        axes[2, c].plot(xc, sig, "o-", color="#2E9E4F", ms=7, lw=1.6, mfc="#2E9E4F" if fill else "white",
                        label=lab if c == 0 else None)
    if ALT:
        H2 = ALT[0]; sig = []
        for b in range(len(edges) - 1):
            m = (H2[key] >= edges[b]) & (H2[key] < edges[b + 1]) & (H2["core_n"] >= args.min_trk)
            sig.append(C.dgaus_core(H2["core_t"][m] - H2["t0"][m])[0] if m.sum() >= 300 else np.nan)
        axes[2, c].plot(xc, sig, "o--", color="#2E9E4F", ms=5, lw=1.2, alpha=0.8,
                        label=f"as used, jet-time cone $\\Delta R < {ALT[2]:.1f}$" if c == 0 else None)
    axes[2, c].set_xlabel(xlab)
    if key in CAT:
        if key == "pt":
            labs = [f"{edges[b]:.0f}-{edges[b + 1]:.0f}" if edges[b + 1] < 1000 else f"> {edges[b]:.0f}" for b in range(len(edges) - 1)]
        else:
            labs = [f"{edges[b]}" if edges[b + 1] - edges[b] == 1 else (f"{edges[b]}-{edges[b + 1] - 1}" if edges[b + 1] < 100 else f"$\\geq${edges[b]}")
                    for b in range(len(edges) - 1)]
        axes[2, c].set_xticks(xc); axes[2, c].set_xticklabels(labs)
        for a in axes[:, c]: a.set_xlim(-0.5, len(edges) - 1.5)
    else:
        for a in axes[:, c]: a.axvspan(2.38, 4.0, color="#EEEEEE", zorder=0)
    for a in axes[:, c]: a.grid(alpha=0.3)
axes[0, 0].set_ylabel(f"efficiency: jet time from\n$\\geq${args.min_trk} tracks, of jets in HGTD")
axes[1, 0].set_ylabel("HS jets: core fraction\n$|t_{\\rm jet} - t_{\\rm HS}| < 60$ ps")
axes[2, 0].set_ylabel("HS jets: core $\\sigma$ of\n$t_{\\rm jet} - t_{\\rm HS}$ [ps]")
for a in axes[0]: a.set_ylim(0, 1.02)
for a in axes[1]: a.set_ylim(0, 1.02)
for a in axes[2]: a.set_ylim(0, 35)
axes[0, 0].legend(fontsize=10, frameon=False, loc="lower center")
axes[1, 0].legend(fontsize=9.5, frameon=False, loc="lower center")
axes[2, 0].legend(fontsize=10, frameon=False, loc="upper left")
fig.suptitle("Tagging-jet times (HGTD): VBF $H\\to$inv. and $Z(\\to\\ell\\ell)$+jets, JVT + fJVT loose, all selected events;  "
             "shaded: the jet axis inside HGTD;  core fraction: an untimed jet counts as a failure   ATLAS Simulation Internal",
             fontsize=10.5, y=0.91)
for ext in ("pdf", "png"):
    fig.savefig(f"{args.out}.{ext}", dpi=150, bbox_inches="tight")
print("wrote", args.out + ".pdf / .png")

# the inclusive numbers, for the caption
for pop, nm in [(HS, "HS (VBF)"), (PU_Z, "PU (Z+jets)")]:
    inh = pop["hgtd_ntimed"] >= 1; timed = pop["core_n"] >= args.min_trk
    good = (pop["core_hsfrac"] >= 0.5) if pop is HS else (pop["core_hsfrac"] < 0.5)
    print(f"{nm:18s}: {int(inh.sum()):,} jets in HGTD; timed (>= {args.min_trk} tracks) {100 * frac(pop['w'], timed, inh)[0]:.1f}%; "
          f"of the timed, right-origin cluster {100 * frac(pop['w'], good, timed)[0]:.1f}%")
inh = HS["hgtd_ntimed"] >= 1
for tcol, ncol, nm in [("core_t", "core_n", "as used"), ("pure_hs_t", "pure_hs_n", "perfect assignment")]:
    ok = (HS[ncol] >= args.min_trk) & (np.abs(HS[tcol] - HS["t0"]) < 60)
    print(f"HS core fraction (of jets in HGTD, untimed = failure), {nm}: {100 * frac(HS['w'], ok, inh)[0]:.1f}%")
print(f"HS timed with a perfect assignment (>= {args.min_trk} timed HS tracks), of jets in HGTD: {100 * frac(HS['w'], HS['pure_hs_n'] >= args.min_trk, inh)[0]:.1f}%")
if ALT:
    H2, P2, dr2 = ALT; inh2 = H2["hgtd_ntimed"] >= 1
    ok = (H2["core_n"] >= args.min_trk) & (np.abs(H2["core_t"] - H2["t0"]) < 60)
    m2 = H2["core_n"] >= args.min_trk
    print(f"cone dR < {dr2}: HS timed {100 * frac(H2['w'], m2, inh2)[0]:.1f}%, core fraction {100 * frac(H2['w'], ok, inh2)[0]:.1f}%, "
          f"core sigma {C.dgaus_core(H2['core_t'][m2] - H2['t0'][m2])[0]:.1f} ps ({100 * np.mean(np.abs(H2['core_t'][m2] - H2['t0'][m2]) > 60):.1f}% beyond 60 ps); "
          f"PU timed {100 * frac(P2['w'], P2['core_n'] >= args.min_trk, P2['hgtd_ntimed'] >= 1)[0]:.1f}%; "
          f"timed HS jets with a mostly-HS cluster {100 * frac(H2['w'], H2['core_hsfrac'] >= 0.5, m2)[0]:.1f}%")
mt = HS["core_n"] >= args.min_trk; mp = HS["pure_hs_n"] >= args.min_trk
print(f"HS core sigma: as used {C.dgaus_core(HS['core_t'][mt] - HS['t0'][mt])[0]:.1f} ps ({100 * np.mean(np.abs(HS['core_t'][mt] - HS['t0'][mt]) > 60):.1f}% beyond 60 ps); "
      f"perfect assignment {C.dgaus_core(HS['pure_hs_t'][mp] - HS['t0'][mp])[0]:.1f} ps ({100 * np.mean(np.abs(HS['pure_hs_t'][mp] - HS['t0'][mp]) > 60):.1f}% beyond 60 ps)")
