#!/usr/bin/env python3
"""What the timing veto adds on top of each pileup-tagging baseline -- the figure of slides_tdr_timing's backup
slide "If the forward jets were tagged with tracks instead of fJVT" (the user's request, 2026-11-05).

Events: tdr_fwdpair --pair=one trees, without (--jvt=none) and with JVT + fJVT loose before pairing (--jvt=loose);
both from the same generated events with the same weights, so their weighted yields compare directly. Strategies:
  JVT + fJVT loose            the --jvt=loose trees (the analysis baseline, section 3)
  ITk R_pT                    the --jvt=none trees, both VBS jets passing the TDR's ITk-only R_pT tagger at its
                              85% working points (python/tdr_pujet_plots.py's wp json)
  ... + WAVeS veto            and the 3-sigma jet-vs-WAVeS-t0 veto (tdr_fwdpair_common.vetoes)
Each strategy's own kinematic BDT (tdr_fwdpair_kinbdt.py on its trees) is cut at the signal efficiency on the x axis,
relative to that strategy's events before the BDT. Upper panel: S/sqrt(B) relative to no pileup tagging at the same
BDT signal efficiency (the --jvt=none trees, all events). Lower panel: the veto's own gain, with / without, per
baseline. Errors: weighted binomial on the background efficiency only.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_baselines.py \\
      --sig ../condor/vbf/vbf_onefwd_tdr_fwdpair.root --bkg ../condor/zjets/zjets_onefwd_tdr_fwdpair.root \\
      --bdt ../condor/zjets/vbf_vs_zjets_onefwd_tdr_fwdpair_kinbdt.npz \\
      --sig-jvt ../condor/vbf/vbf_onefwd_jvtLoose_tdr_fwdpair.root --bkg-jvt ../condor/zjets/zjets_onefwd_jvtLoose_tdr_fwdpair.root \\
      --bdt-jvt ../condor/zjets/vbf_vs_zjets_onefwd_jvtLoose_kinbdt.npz \\
      --wp ../figs/tdr/grid/grid_tdr_pujet_wp85.json --out ../figs/tdr/grid/onefwd_jvt/grid_onefwd_baselines
"""
import argparse
import json

import numpy as np
import uproot

import tdr_fwdpair_common as FC
import tdr_root_style as S

ap = argparse.ArgumentParser()
for k in ("sig", "bkg", "bdt", "sig_jvt", "bkg_jvt", "bdt_jvt", "wp", "out"):
    ap.add_argument("--" + k.replace("_", "-"), required=True)
args = ap.parse_args()
R = S.setup()
EFFS = np.round(np.arange(1.0, 0.249, -0.05), 2)


def load(sp, bp, npz):
    D = np.load(npz)
    Es, Eb = (uproot.open(p)["events"].arrays(library="np") for p in (sp, bp))
    cal = FC.calibration(Es, verbose=False)
    return dict(Es=Es, Eb=Eb, ss=D["sig_score"], sb=D["bkg_score"], ws=Es["weight"].astype(float), wb=Eb["weight"].astype(float),
                Ks=FC.vetoes(Es, cal)["waves"], Kb=FC.vetoes(Eb, cal)["waves"])


def itk_pass(E, wp):
    th = {k: wp[k]["threshold"] for k in ("central_30_50", "central_50_inf", "forward_30_50", "forward_50_inf")}
    ok = np.ones(len(E["weight"]), bool)
    for k in (1, 2):
        ae, pt, r = np.abs(E[f"j{k}_eta"]), E[f"j{k}_pt"], E[f"j{k}_rpt_itk"]
        ok &= r >= np.where(ae < 2.4, np.where(pt < 50, th["central_30_50"], th["central_50_inf"]),
                            np.where(pt < 50, th["forward_30_50"], th["forward_50_inf"]))
    return ok


N = load(args.sig, args.bkg, args.bdt)                   # no JVT: the reference and the ITk strategy
J = load(args.sig_jvt, args.bkg_jvt, args.bdt_jvt)       # the analysis baseline (JVT + fJVT loose, or JVT to |eta| < 4.0)
_m = uproot.open(args.sig_jvt)["meta"]
BASE = FC.jvt_short(str(_m.arrays(library="np")["jvt_wp"][0]) if "jvt_wp" in _m.keys() else "loose")
WP = json.load(open(args.wp))["working_points"]["rpt_itk"]
pS, pB = itk_pass(N["Es"], WP), itk_pass(N["Eb"], WP)


def yields(T, bs, bb, veto, e):
    """S, B, var(B) after the strategy's events (bs / bb), its BDT at signal efficiency e, and optionally the veto"""
    ss, ws = T["ss"][bs], T["ws"][bs]
    thr = -np.inf if e >= 1.0 else ss[np.argsort(-ss)][np.searchsorted(np.cumsum(ws[np.argsort(-ss)]) / ws.sum(), e)]
    ms, mb = bs & (T["ss"] >= thr), bb & (T["sb"] >= thr)
    if veto: ms, mb = ms & T["Ks"], mb & T["Kb"]
    return T["ws"][ms].sum(), T["wb"][mb].sum(), np.sum(T["wb"][mb] ** 2)


allS, allB = np.ones(len(N["ws"]), bool), np.ones(len(N["wb"]), bool)
jS, jB = np.ones(len(J["ws"]), bool), np.ones(len(J["wb"]), bool)
STRATS = [(BASE, J, jS, jB, False, "#000000", 20, 1),
          (BASE + " + WAVeS veto", J, jS, jB, True, S.COL["waves"], 20, 1),
          ("ITk R_{p_{T}} (85%)", N, pS, pB, False, "#000000", 24, 2),
          ("ITk R_{p_{T}} (85%) + WAVeS veto", N, pS, pB, True, S.COL["waves"], 24, 2)]
ref = [yields(N, allS, allB, False, e) for e in EFFS]
rows = []
print("S/sqrt(B) relative to no pileup tagging, at the BDT signal efficiency " + " ".join(f"{e:.2f}" for e in EFFS))
for name, T, bs, bb, veto, colour, mk, ls in STRATS:
    y, ey = [], []
    for (s0, b0, _), e in zip(ref, EFFS):
        s, b, vb = yields(T, bs, bb, veto, e)
        f = (s / s0) / np.sqrt(b / b0)
        y.append(f); ey.append(f * 0.5 * np.sqrt(vb) / b)     # Poisson-like, weighted: sqrt(sum w^2) / sum w
    rows.append((name, colour, mk, ls, np.array(y), np.array(ey)))
    print(f"  {name.replace('#', '').replace('_{p_{T}}', 'pT'):34s} " + " ".join(f"{a:.2f}" for a in y))
gain = {}
for base, (a, b) in ((BASE, (0, 1)), ("ITk R_{p_{T}} (85%)", (2, 3))):
    gain[base] = rows[b][4] / rows[a][4]
print("the WAVeS veto's own gain:")
for k, v in gain.items(): print(f"  on {k.replace('#', '').replace('_{p_{T}}', 'pT'):24s} " + " ".join(f"{a:.3f}" for a in v))

MAIN, RAT = 0.68, 0.32
TS, TSR = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)
cv = S.keep(R.TCanvas("baselines", "", 800, 800))
cv.cd(); S.pad("bl_m", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
top = max(np.nanmax(r[4] + r[5]) for r in rows)
fr = S.keep(R.gPad.DrawFrame(0.225, 0.0, 1.025, 1.9 * top))
S.axes(fr, "", "S/#sqrt{B} relative to no pileup tagging", size=TS, ytoff=1.35)
fr.GetXaxis().SetLabelSize(0)
one = S.keep(R.TLine(0.225, 1.0, 1.025, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
leg = S.legend(0.19, 0.40, 0.66, 0.64, size=0.62 * TS, margin=0.15)   # upper left: the curves stay below ~1.3 there
for name, colour, mk, ls, y, ey in rows:
    g = S.style(S.graph(EFFS, y, ey), colour, mk, ls, 2, 1.2); g.Draw("LP SAME")
    leg.AddEntry(g, name, "LP")
leg.Draw()
S.header(0.19, 0.88, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200", "S: VBF H#rightarrowinv., B: Z+jets, #geq 1 VBS jet in HGTD",
                      "Kinematic BDT, then the 3#sigma WAVeS t_{0} veto"], size=0.8 * TS)
R.gPad.RedrawAxis()
cv.cd(); S.pad("bl_r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
fr2 = S.keep(R.gPad.DrawFrame(0.225, 0.85, 1.025, 1.45))
S.axes(fr2, "Signal efficiency of the kinematic BDT cut (1 = no BDT)", "Veto: with / without", size=TSR, ytoff=1.35 * TS / TSR, xtoff=1.2, ndiv_y=503)
one2 = S.keep(R.TLine(0.225, 1.0, 1.025, 1.0)); one2.SetLineStyle(2); one2.SetLineColor(R.kGray + 2); one2.Draw()
leg2 = S.legend(0.19, 0.70, 0.58, 0.94, size=0.62 * TSR, margin=0.2)   # left: the solid curve is below 1.2 there
for (base, v), ls, mk in zip(gain.items(), (1, 2), (20, 24)):
    g2 = S.style(S.graph(EFFS, v), S.COL["waves"], mk, ls, 2, 1.0); g2.Draw("LP SAME")
    leg2.AddEntry(g2, f"on {base}", "LP")
leg2.Draw()
S.save(cv, args.out)
