#!/usr/bin/env python3
"""ATL-HGTD-PUB-2022-001 Sec. 6 figures from util/tdr_fwdpair.cxx's trees (+ tdr_fwdpair_bdt.py's scores).

Signal: VBF events whose two jets are both hard-scatter; background: events with at least one pileup jet
(the note's Sec. 3). Rejection = 1 / (weighted background efficiency) against the weighted signal
efficiency, sampled on a common efficiency grid (linear interpolation between thresholds), binomial
errors on the effective background count.

  <out>_jvfcut[_jz<k>]  the note's Fig 2: the same minimum JVF on both jets, per t0 (no timing / Athena's /
                        TRKPTZ / WAVeS / truth); a jet without associated tracks fails every cut, so these
                        curves end below 100% (the note: "the events where one jet does not have any
                        associated track are lost")
  <out>_bdt[_jz<k>]     the note's Fig 3: the BDTs of its Table 1, every event scored by the fold that did not train on it -- JVF only, JVF with
                        each t0, and JVF plus the two jets' time differences (no t0). With --bdt-extra, also
                        the PV0-track jet time (beyond the note)
--bkg = --sig: the VBF sample's own events with a pileup jet as the background (34,760 events on the grid,
10x the Z+jets ones; the TDR's and the note's per-jet figures take their pileup jets from VBF events too).
Points with fewer than 3 effective background events surviving are not drawn.
Dijet (--bkg-label Dijet): one pair of figures per JZ slice -- the slices have no cross-sections, and a
rejection within one slice needs none.

--set veto (with tdr_fwdpair_bdt.py --set veto; trees with the core jet time): our timing study's vetoes instead
of the note's timing --
  <out>_jvfveto  the same JVF cut on both jets, alone and with each 3-sigma veto applied on top (jet vs jet;
                 jet vs HGTD / TRKPTZ / WAVeS / truth t0; tdr_fwdpair_common.vetoes): a vetoed event never passes,
                 so a vetoed curve ends at the veto's signal efficiency
  <out>_bdtveto  BDT on the two JVFs, and with the veto pulls as extra inputs
The selection line follows the trees' pair mode (the note's, or at least one jet in HGTD).

Run with ~/.venv-hgtd (ROOT):
  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_plots.py --sig <vbf>.root --bkg <zjets>.root \\
      --bdt <x>.npz --bkg-label Z+jets --out ../figs/tdr/grid/grid_fwdpair_zjets
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S
import tdr_fwdpair_common as FC

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--bdt", default="")
ap.add_argument("--bkg-label", default="Z+jets")
ap.add_argument("--bdt-extra", action="store_true", help="also draw the PV0-track jet-time BDT")
ap.add_argument("--min-bkg", type=int, default=200, help="dijet: skip slices with fewer background events")
ap.add_argument("--set", default="note", choices=["note", "veto"])
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()

ENERGY = "#sqrt{s} = 14 TeV, #LT#mu#GT = 200"
SEL = "p_{T}^{j1} > 50 GeV, p_{T}^{j2} > 30 GeV, 2.4 < |#eta^{jet}| < 3.8"
EFF = np.round(np.arange(0.60, 0.9801, 0.02), 4)
T0S = [("hgtd", "HGTD t_{0}", S.COL["hgtd"], 23, 1), ("trkptz", "TRKPTZ t_{0}", S.COL["trkptz"], 21, 1),
       ("waves", "WAVeS t_{0}", S.COL["waves"], 20, 1), ("truth", "Truth t_{0}", S.COL["truth"], 25, 3)]
DT = ("dtime", "#Deltat of the jets (no t_{0})", "kP10Violet", 33, 10)
DT_HS = ("dtime_hs", "#Deltat, PV0 tracks (no t_{0})", "kP10Green", 27, 7)
MAIN, RAT = 0.68, 0.32
TS_M, TS_R = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)


def plain(s):
    return s.replace("#", "").replace("{", "").replace("}", "").replace("_", "")


def roc(ds, ws, db, wb):
    """background efficiency (and its error) at each signal efficiency of EFF; discriminant: higher = more
    signal-like; NaN never passes (it still counts in the totals)"""
    Ws, Wb = ws.sum(), wb.sum()
    neff = Wb ** 2 / np.sum(wb ** 2)
    vs, vb = np.isfinite(ds), np.isfinite(db)
    v = np.concatenate([ds[vs], db[vb]])
    w = np.concatenate([ws[vs], wb[vb]])
    lab = np.concatenate([np.ones(vs.sum(), bool), np.zeros(vb.sum(), bool)])
    o = np.argsort(-v, kind="stable")
    v, w, lab = v[o], w[o], lab[o]
    cs = np.cumsum(np.where(lab, w, 0.0)) / Ws
    cb = np.cumsum(np.where(lab, 0.0, w)) / Wb
    last = np.r_[v[1:] != v[:-1], True]                  # one point per distinct threshold
    es, eb = np.maximum.accumulate(np.r_[0.0, cs[last]]), np.maximum.accumulate(np.r_[0.0, cb[last]])
    pe = np.interp(EFF, es, eb, right=np.nan)            # beyond the reachable efficiency: undefined
    pe = np.where(pe * neff < 3.0, np.nan, pe)           # under 3 effective background events: no measurement
    return pe, np.sqrt(np.clip(pe * (1 - pe), 0, None) / neff), es[-1]


def draw(rows, title, lines, out, xlab="Efficiency for signal events (two HS jets)"):
    """rows: (name, colour, marker, line style, background efficiency, error); the first is the reference"""
    ref = None
    print(f"\n{title}: rejection (ratio to {plain(rows[0][0])}) at signal efficiency 0.70 / 0.80 / 0.90")
    curves = []
    for name, colour, mk, ls, pe, perr, emax in rows:
        with np.errstate(divide="ignore", invalid="ignore"):
            rej = np.where(pe > 0, 1.0 / pe, np.nan)
            rerr = np.where(pe > 0, perr / pe ** 2, np.nan)
        if ref is None: ref = rej
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = rej / ref
        curves.append((name, colour, mk, ls, rej, rerr, ratio))
        at = [int(np.argmin(np.abs(EFF - e))) for e in (0.70, 0.80, 0.90)]
        print(f"  {plain(name):34s}" + "".join(f"{rej[i]:9.1f} ({ratio[i]:4.2f})" for i in at) + f"   max eff {emax:.3f}")
    ymax = 10 ** np.ceil(np.log10(100.0 * np.nanmax([np.nanmax(c[4] + np.nan_to_num(c[5])) for c in curves])))
    cv = S.keep(R.TCanvas(out, "", 800, 800))
    cv.cd(); p1 = S.pad(out + "_m", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
    p1.SetLogy()
    fr = S.keep(p1.DrawFrame(EFF[0], 1.0, 1.0, ymax))
    S.axes(fr, "", f"Rejection of {args.bkg_label} with #geq 1 PU jet", size=TS_M, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0)
    leg = S.legend(0.19, 0.05, 0.62, 0.05 + 0.05 * len(curves), size=0.62 * TS_M, margin=0.2)   # bottom-left: empty on a falling ROC
    for name, colour, mk, ls, rej, rerr, _ in curves:
        ok = np.isfinite(rej)
        g = S.style(S.graph(EFF[ok], rej[ok]), colour, mk, ls, 2, 1.1); g.Draw("L SAME")
        S.style(S.graph(EFF[ok], rej[ok], rerr[ok]), colour, mk, 1, 2, 1.1).Draw("P SAME")
        leg.AddEntry(g, name, "LP")
    leg.Draw()
    vs = (f"VBF H#rightarrowinv.: 2 HS jets vs #geq 1 PU jet" if args.bkg == args.sig else f"VBF H#rightarrowinv. vs {args.bkg_label}")
    S.header(0.20, 0.88, [ENERGY, vs, SEL] + lines, size=0.85 * TS_M)
    cv.cd(); p2 = S.pad(out + "_r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
    rmax = np.ceil(1.15 * np.nanmax([np.nanmax(c[6][np.isfinite(c[6])]) if np.isfinite(c[6]).any() else 1 for c in curves]) * 2) / 2
    fr2 = S.keep(p2.DrawFrame(EFF[0], 0.5, 1.0, max(rmax, 1.5)))
    S.axes(fr2, xlab, f"Ratio to {curves[0][0].split(' (')[0]}", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2, ndiv_y=505)
    one = S.keep(R.TLine(EFF[0], 1.0, 1.0, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
    for name, colour, mk, ls, rej, rerr, ratio in curves[1:]:
        ok = np.isfinite(ratio)
        S.style(S.graph(EFF[ok], ratio[ok]), colour, mk, ls, 2, 1.0).Draw("LP SAME")
    S.save(cv, out)


Sg = uproot.open(args.sig)["events"].arrays(library="np")
Bg = uproot.open(args.bkg)["events"].arrays(library="np")
_meta = uproot.open(args.sig)["meta"]
if "pair_mode" in _meta.keys() and str(_meta.arrays(library="np")["pair_mode"][0]) == "one":
    SEL = "j_{1}, j_{2}: |#eta| < 3.8, #geq 1 in 2.4#minus3.8, opposite #eta, > 50 / 30 GeV"
sig = (Sg["j1_hs"] == 1) & (Sg["j2_hs"] == 1)
bkg = (Bg["j1_pu"] == 1) | (Bg["j2_pu"] == 1)
dijet = bool(np.any(Bg["jz"] >= 0))
slices = [None] if not dijet else [z for z in np.unique(Bg["jz"]) if np.sum(bkg & (Bg["jz"] == z)) >= args.min_bkg]
if dijet:
    print("dijet events with >= 1 PU jet per JZ slice: " + ", ".join(f"JZ{z} {int(np.sum(bkg & (Bg['jz'] == z)))}" for z in np.unique(Bg["jz"])))
print(f"signal {int(sig.sum())} events, background {int(bkg.sum())}")
D = np.load(args.bdt) if args.bdt else None


def mincut(E, key):
    """the same cut on both jets = a cut on the smaller JVF; a jet without associated tracks (JVF < 0) never passes"""
    d = np.minimum(E["j1_" + key], E["j2_" + key])
    return np.where((E["j1_" + key] >= 0) & (E["j2_" + key] >= 0), d, np.nan)


VETO_ROWS = [("jj", "jet vs jet", "kP10Violet", 33, 10)] + [(s_, n_.replace(" t_{0}", "") + " t_{0}", c_, m_, l_) for s_, n_, c_, m_, l_ in T0S]
if args.set == "veto":
    cal = FC.calibration(Sg)
    KS, KB = FC.vetoes(Sg, cal), FC.vetoes(Bg, cal)
    ws, wb = Sg["weight"][sig].astype(float), Bg["weight"][bkg].astype(float)
    js, jb = mincut(Sg, "jvf_itk"), mincut(Bg, "jvf_itk")
    rows = [("ITk (JVF)", "#000000", 24, 1, *roc(js[sig], ws, jb[bkg], wb))]
    for key, name, colour, mk, ls in VETO_ROWS:
        rows.append((f"+ {name} veto", colour, mk, ls,
                     *roc(np.where(KS[key], js, np.nan)[sig], ws, np.where(KB[key], jb, np.nan)[bkg], wb)))
    draw(rows, "JVF cut + 3-sigma vetoes", ["JVF, the same cut on both jets, + 3#sigma timing veto"], f"{args.out}_jvfveto")
    if D is not None:
        rows = [("BDT ITk (JVF)", "#000000", 24, 1, *roc(D["sig_itk"], D["sig_weight"], D["bkg_itk"], D["bkg_weight"]))]
        for key, name, colour, mk, ls in VETO_ROWS:
            k = "veto_jj" if key == "jj" else f"veto_{key}"
            rows.append((f"BDT + {name} pulls", colour, mk, ls, *roc(D[f"sig_{k}"], D["sig_weight"], D[f"bkg_{k}"], D["bkg_weight"])))
        draw(rows, "BDTs with the veto pulls", ["BDT: JVFs + timing-veto pulls (two-fold)"], f"{args.out}_bdtveto")
    raise SystemExit(0)

for z in slices:
    mb = bkg if z is None else bkg & (Bg["jz"] == z)
    tag = "" if z is None else f"_jz{z}"
    zline = [] if z is None else [f"Dijet JZ{z} slice"]
    ws, wb = Sg["weight"][sig].astype(float), Bg["weight"][mb].astype(float)
    rows = [("ITk", "#000000", 24, 1, *roc(mincut(Sg, "jvf_itk")[sig], ws, mincut(Bg, "jvf_itk")[mb], wb))]
    for s, name, colour, mk, ls in T0S:
        rows.append((name, colour, mk, ls, *roc(mincut(Sg, f"jvf_t0_{s}")[sig], ws, mincut(Bg, f"jvf_t0_{s}")[mb], wb)))
    draw(rows, f"JVF, same cut on both jets{tag}", ["JVF, the same cut on both jets"] + zline, f"{args.out}_jvfcut{tag}")
    if D is None: continue
    mz = np.ones(len(D["bkg_weight"]), bool) if z is None else D["bkg_jz"] == z
    ws_, wb_ = D["sig_weight"], D["bkg_weight"][mz]
    rows = [("BDT ITk", "#000000", 24, 1, *roc(D["sig_itk"], ws_, D["bkg_itk"][mz], wb_))]
    for s, name, colour, mk, ls in T0S:
        rows.append((f"BDT {name}", colour, mk, ls, *roc(D[f"sig_t0_{s}"], ws_, D[f"bkg_t0_{s}"][mz], wb_)))
    for key, name, colour, mk, ls in [DT] + ([DT_HS] if args.bdt_extra else []):
        rows.append((f"BDT {name}", colour, mk, ls, *roc(D[f"sig_{key}"], ws_, D[f"bkg_{key}"][mz], wb_)))
    draw(rows, f"BDTs{tag}", ["BDT on the two jets (two-fold, held out)"] + zline, f"{args.out}_bdt{tag}")
