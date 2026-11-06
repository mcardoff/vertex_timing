#!/usr/bin/env python3
"""Signal and background efficiency of each timing veto against m_jj (the user, 2026-11-06: "putting signal efficiency
and background efficiency as a function of mjj after the event displays") -- the two halves of the S/sqrt(B) slides:
jet vs jet keeps about as much signal as the WAVeS veto but much more background.

Events: tdr_fwdpair --pair=one --jvt=loose trees (section 3), all pairs or forward-forward pairs only (--ff), with no
BDT or after the kinematic BDT at --bdt-eff. Vetoes as tdr_fwdpair_kinplots.py's compact S/sqrt(B) figures: the 3 sigma
WAVeS and truth t0 vetoes, jet vs jet, and the perfect HGTD veto (truth labels). Efficiencies are weighted, relative to
the events before the veto in each m_jj bin, with weighted-binomial errors; a bin with < --min-mc MC events of a
sample is left out. Upper pad: VBF H->inv. (signal); lower pad: Z+jets (background).

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_effs.py --sig ../condor/vbf/vbf_onefwd_jvtLoose_tdr_fwdpair.root \\
      --bkg ../condor/zjets/zjets_onefwd_jvtLoose_tdr_fwdpair.root --bdt ../condor/zjets/vbf_vs_zjets_onefwd_jvtLoose_kinbdt.npz \\
      --out ../figs/tdr/grid/onefwd_jvt/grid_onefwd_jvt_effs
writes <out>_{all,ff}_{nobdt,bdt80}.png / .pdf
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC
import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True); ap.add_argument("--bkg", required=True); ap.add_argument("--bdt", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--bdt-eff", type=float, default=0.80)
ap.add_argument("--edges", default="0,500,1000,1500,2000,4000", help="m_jj bin edges [GeV]; the last bin includes overflow")
ap.add_argument("--min-mc", type=int, default=5)
ap.add_argument("--nsigma", type=float, default=3.0)
args = ap.parse_args()
R = S.setup()
EDGES = np.array([float(x) for x in args.edges.split(",")])
E = {"sig": uproot.open(args.sig)["events"].arrays(library="np"), "bkg": uproot.open(args.bkg)["events"].arrays(library="np")}
D = np.load(args.bdt); SC = {"sig": D["sig_score"], "bkg": D["bkg_score"]}
w0 = E["sig"]["weight"].astype(float); o = np.argsort(-SC["sig"]); cs = np.cumsum(w0[o]) / w0.sum()
THR = SC["sig"][o][np.searchsorted(cs, args.bdt_eff)]
cal = FC.calibration(E["sig"], verbose=False)
K = {}
for k, e in E.items():
    K[k] = FC.vetoes(e, cal, args.nsigma)
    fwd = lambda j: np.abs(e[f"j{j}_eta"]) > 2.4
    K[k]["perfect"] = ~(((e["j1_pu"] == 1) & fwd(1)) | ((e["j2_pu"] == 1) & fwd(2)))
ROWS = [("waves", "WAVeS t_{0}", S.COL["waves"], 20, 1), ("truth", "Truth t_{0}", S.COL["truth"], 25, 3),
        ("jj", "Jet vs jet (no t_{0})", "kP10Violet", 33, 10), ("perfect", "Perfect HGTD veto (truth)", "#000000", 34, 2)]
JVT = "JVT + fJVT loose"


def effs(key, base, v):
    e = E[key]; w = e["weight"].astype(float); x = np.clip(e["mjj"], EDGES[0], EDGES[-1] - 1e-6)
    out = []
    for c in range(len(EDGES) - 1):
        m = base & (x >= EDGES[c]) & (x < EDGES[c + 1])
        if m.sum() < args.min_mc: out.append((np.nan, np.nan)); continue
        kk = K[key][v][m]; wm = w[m]
        ef = wm[kk].sum() / wm.sum()
        var = ((wm[kk] ** 2).sum() * (1 - ef) ** 2 + (wm[~kk] ** 2).sum() * ef ** 2) / wm.sum() ** 2
        out.append((ef, np.sqrt(var)))
    allm = base
    kk = K[key][v][allm]; wm = E[key]["weight"].astype(float)[allm]
    return out, wm[kk].sum() / wm.sum()


for ff, ftag in ((False, "all"), (True, "ff")):
    for emj, btag in ((None, "nobdt"), (args.bdt_eff, f"bdt{round(100 * args.bdt_eff)}")):
        base = {}
        for k, e in E.items():
            b = np.ones(len(e["weight"]), bool) if emj is None else SC[k] >= THR
            if ff: b &= e["cat"] == 2
            base[k] = b
        stem = f"{args.out}_{ftag}_{btag}"
        cv = S.keep(R.TCanvas(stem, "", 800, 860))
        TOP, MID = 0.80, 0.47                                   # header band above 0.80; signal 0.47-0.80; background below
        cv.cd(); S.pad(stem + "_h", 0.0, TOP, 1.0, 1.0, top=0.0, bottom=0.0)
        S.header(0.16, 0.80, ["#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200", ("Forward-forward pairs, " if ff else "#geq 1 VBS jet in HGTD, ") + JVT,
                              "No BDT: the veto alone" if emj is None else f"BDT at #varepsilon_{{S}} = {100 * emj:.0f}%, then the veto"],
                 size=0.11, step=0.15)
        leg = S.legend(0.60, 0.06, 0.95, 0.92, size=0.09, margin=0.14)
        leg.SetHeader(f"{args.nsigma:g}#sigma timing veto (all m_{{jj}}: #varepsilon_{{S}} / #varepsilon_{{B}}):")
        print(f"\n{stem.split('/')[-1]}: efficiency per m_jj bin " + ", ".join(f"[{EDGES[c]:.0f},{EDGES[c + 1]:.0f})" for c in range(len(EDGES) - 1)))
        graphs = []
        for side, (key, ylab, ylo, yhi, y0, y1) in enumerate((("sig", "VBF H#rightarrowinv. kept", 0.55, 1.05, MID, TOP),
                                                               ("bkg", "Z+jets kept", 0.0, 1.05, 0.0, MID))):
            cv.cd(); S.pad(f"{stem}_{key}", 0.0, y0, 1.0, y1, top=0.04 if side == 0 else 0.0, bottom=0.03 if side == 0 else 0.25)
            ph = y1 - y0
            ts = 0.05 * 0.33 / ph * 1.6
            xlo = 500.0 if ff else EDGES[0]                      # forward-forward pairs start at ~500 GeV
            fr = S.keep(R.gPad.DrawFrame(xlo, ylo, EDGES[-1], yhi))
            S.axes(fr, "m_{jj} [GeV]   (last bin includes overflow)" if side == 1 else "", ylab, size=ts, ytoff=1.0 * 0.05 / ts * 1.25,
                   xtoff=1.1, ndiv_y=505)
            if side == 0: fr.GetXaxis().SetLabelSize(0)
            one = S.keep(R.TLine(xlo, 1.0, EDGES[-1], 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
            for c in range(1, len(EDGES) - 1):
                if EDGES[c] <= xlo: continue
                ln = S.keep(R.TLine(EDGES[c], ylo, EDGES[c], yhi)); ln.SetLineStyle(3); ln.SetLineColor(R.kGray + 1); ln.Draw()
            for k_, (v, name, colour, mk, ls) in enumerate(ROWS):
                vals, allv = effs(key, base[key], v)
                # around the bin centre, a fixed small step per veto (the widest bin would scatter them otherwise)
                xs = np.array([0.5 * (EDGES[c] + EDGES[c + 1]) + (k_ - 1.5) * 45.0 for c in range(len(EDGES) - 1)])
                g = S.style(S.graph(xs, np.array([a for a, _ in vals]), np.array([b for _, b in vals])), colour, mk, ls, 2, 1.2)
                g.Draw("PL SAME")
                if side == 0: graphs.append((g, name, allv))
                else: graphs[k_] = (graphs[k_][0], graphs[k_][1], (graphs[k_][2], allv))
                print(f"  {key} {name.replace('#', '').replace('_{0}', '0'):26s} " + " ".join("  --  " if not np.isfinite(a) else f"{a:.3f}" for a, _ in vals) + f"   all {allv:.3f}")
            R.gPad.RedrawAxis()
        for g, name, (es, eb) in graphs:
            leg.AddEntry(g, f"{name}: {es:.2f} / {eb:.2f}", "PL")
        cv.cd(); cv.GetListOfPrimitives().First().cd(); leg.Draw()
        S.save(cv, stem)
