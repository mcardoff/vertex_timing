#!/usr/bin/env python3
"""TDR Sec. 3.3.3 figures from util/tdr_pujet's per-jet tree, in the ATLAS style (python/tdr_root_style.py).

  <out>_roc      TDR Fig 3.20: pileup jet rejection (1 / pileup efficiency, linear axis) at the sampled
                 hard-scatter efficiencies 0.800-0.975, 2.4 < |eta| < 4.0, for 30-50 GeV (left) and
                 > 50 GeV (right) jets; the t0-gated taggers (HGTD, TRKPTZ, WAVeS, truth) against
                 ITk only, with the ratio to ITk underneath. The tree's self-tagging and combined
                 columns are not drawn (user, 2026-11-03: too many series; combined = HGTD wherever
                 Athena has a t0)
  <out>_eta      TDR Fig 3.21: pileup jet rate relative to ITk only, against |eta|, at 85% hard-scatter
                 efficiency (the working point set per |eta| bin), 30-50 GeV
  <out>_wp85.json  every tagger's R_pT threshold at 85% hard-scatter efficiency per pT bin (30-50,
                 > 50 GeV) and region (central |eta| < 2.4, forward 2.4-4.0): the TDR preselection's
                 "every jet passes the ITk-only R_pT tagger at 85%", and the timing taggers that
                 replace it, read by tdr_vbf_plots.py

Jets are counted unweighted, as the TDR's per-jet performance. A jet passes a tagger when its R_pT is
>= the threshold; the curve is traced over every threshold, and rejection is read at matched
hard-scatter efficiency, never by snapping to bins. Between two thresholds the pileup efficiency is
interpolated linearly: that is the tagger accepting the jets between them at random. Above the last
threshold before R_pT = 0 (every tagger has a spike of hard-scatter jets at R_pT = 0, larger with a
wrong t0) that interpolation is the only way to reach the efficiency, so the curve there is a
mixture, not a cut. The working points in the json are real cuts: the highest threshold keeping >= 85%.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_pujet_plots.py \\
        ../figs/hists/tdr_pujet.root --label "VBF H#rightarrowinv. (local sample)" --out ../figs/tdr/local_tdr_pujet
"""
import argparse, json, os
import numpy as np
import uproot
import tdr_root_style as S

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("tree", help="tdr_pujet output")
ap.add_argument("--label", default="VBF H#rightarrowinv.", help="sample, in ROOT TLatex")
ap.add_argument("--out", required=True)
args = ap.parse_args()
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
R = S.setup()

# taggers: (column, legend, short legend, colour, ROOT marker, ROOT line style); ITk first (the reference)
TAGGERS = [("rpt_itk", "ITk only", "ITk only", S.COL["itk"], 24, 2),
           ("rpt_t0_hgtd", "HGTD (Athena) t_{0}", "HGTD t_{0}", S.COL["hgtd"], 23, 1),
           ("rpt_t0_trkptz", "TRKPTZ t_{0}", "TRKPTZ t_{0}", S.COL["trkptz"], 21, 1),
           ("rpt_t0_waves", "WAVeS t_{0}", "WAVeS t_{0}", S.COL["waves"], 20, 1),
           ("rpt_t0_truth", "Truth t_{0} (reference)", "Truth t_{0}", S.COL["truth"], 25, 3)]
EFF = np.round(np.arange(0.800, 0.9751, 0.025), 3)    # the sampled working points
plain = lambda s: s.replace("_{0}", "0").replace("#", "")

t = uproot.open(args.tree)
J = t["jets"].arrays(library="np")
meta = t["meta"].arrays(library="np")
gate = float(np.unique(meta["gate_sigma"])[0])
clean = str(meta["time_cleaning"][0]) if "time_cleaning" in meta else "off"   # trees before 2026-11-03: uncleaned
CLEAN_TXT = {"all": "Cleaned track times (TDR)", "t0": "Cleaned t_{0} only", "off": "Track times as in the ntuple"}

ae = np.abs(J["eta"])
ENERGY = "#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200"


def curve(hs, pu):
    """the ROC: for each distinct hard-scatter efficiency (ascending), the lowest pileup efficiency
    reaching it and the threshold (the highest one) that does; pass = R_pT >= threshold"""
    th = np.unique(np.concatenate([hs, pu]))
    e_hs = 1.0 - np.searchsorted(np.sort(hs), th, side="left") / len(hs)
    e_pu = 1.0 - np.searchsorted(np.sort(pu), th, side="left") / len(pu)
    u, inv = np.unique(e_hs, return_inverse=True)
    pmin = np.full(len(u), np.inf); np.minimum.at(pmin, inv, e_pu)
    tmax = np.full(len(u), -np.inf); np.maximum.at(tmax, inv, th)
    return u, pmin, tmax


def at_eff(hs, pu, effs):
    """pileup efficiency and its binomial error at each hard-scatter efficiency (interpolated)"""
    e_hs, e_pu, _ = curve(hs, pu)
    pe = np.interp(effs, e_hs, e_pu)
    return pe, np.sqrt(np.clip(pe * (1 - pe), 0, None) / len(pu))


def working_point(hs, pu, eff):
    """a real cut: the highest threshold keeping >= eff of the hard-scatter jets"""
    e_hs, e_pu, th = curve(hs, pu)
    i = np.flatnonzero(e_hs >= eff)[0]
    return float(th[i]), float(e_hs[i]), float(e_pu[i])


def sel(pt_lo, pt_hi, eta_lo, eta_hi):
    return (J["pt"] > pt_lo) & (J["pt"] < pt_hi) & (ae >= eta_lo) & (ae < eta_hi)


SLICES = [("30 < p_{T} < 50 GeV", 30.0, 50.0), ("p_{T} > 50 GeV", 50.0, 1e9)]

# ── Fig 3.20 analogue: two 800 x 800 panels side by side, each a main pad over a ratio pad ──
cv = S.keep(R.TCanvas("roc", "", 1600, 800))
MAIN, RAT = 0.68, 0.32                     # pad heights: text sizes are scaled to match 0.05 of 600 px
TS_M, TS_R = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)
print(f"{args.tree}: gate {gate:g} sigma, track-time cleaning {clean}\n\nrejection at matched hard-scatter efficiency, 2.4 < |eta| < 4.0 (ratio to ITk in brackets)")
for c, (lab, lo, hi) in enumerate(SLICES):
    m = sel(lo, hi, 2.4, 4.0)
    nhs, npu = int(np.sum(m & (J["hs"] == 1))), int(np.sum(m & (J["pu"] == 1)))
    print(f"\n  {plain(lab)}: {nhs} hard-scatter / {npu} pileup jets")
    print("  " + f"{'tagger':40s}" + "".join(f"{e:>14.3f}" for e in EFF))
    rows, ref = [], None
    for col, name, short, colour, mk, ls in TAGGERS:
        hs, pu = J[col][m & (J["hs"] == 1)], J[col][m & (J["pu"] == 1)]
        pe, perr = at_eff(hs, pu, EFF)
        rej = np.where(pe > 0, 1.0 / pe, np.nan)
        rerr = np.where(pe > 0, perr / pe ** 2, np.nan)
        if ref is None: ref = rej
        rows.append((rej, rerr, rej / ref))
        print("  " + f"{plain(name):40s}" + "".join(f"{r:7.1f} ({q:4.2f})" for r, q in zip(rej, rows[-1][2])))
    ymax = 1.7 * np.nanmax([np.nanmax(r[0] + r[1]) for r in rows])
    rmax = np.ceil(1.1 * np.nanmax([np.nanmax(r[2]) for r in rows]) / 0.5) * 0.5

    cv.cd()
    p1 = S.pad(f"main{c}", 0.5 * c, RAT, 0.5 * c + 0.5, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
    fr = S.keep(p1.DrawFrame(0.80, 0.0, 1.0, ymax))
    S.axes(fr, "", "Pileup jet rejection", size=TS_M, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0); fr.GetXaxis().SetNdivisions(505)
    leg = S.legend(0.52, 0.32, 0.94, 0.56, size=0.75 * TS_M)
    for (col, name, short, colour, mk, ls), (rej, rerr, _) in zip(TAGGERS, rows):
        g = S.style(S.graph(EFF, rej), colour, mk, ls, 2, 1.3); g.Draw("L SAME")      # the line in its style
        S.style(S.graph(EFF, rej, rerr), colour, mk, 1, 2, 1.3).Draw("P SAME")         # solid error bars
        leg.AddEntry(g, name, "LP")
    leg.Draw()
    S.header(0.20, 0.87, [ENERGY, f"{args.label}, t_{{0}} gate {gate:g}#sigma", f"{lab}, 2.4 < |#eta| < 4.0", CLEAN_TXT[clean]], size=TS_M)

    cv.cd()
    p2 = S.pad(f"ratio{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, RAT, top=0.03, bottom=0.36)
    fr2 = S.keep(p2.DrawFrame(0.80, 0.0, 1.0, rmax))
    S.axes(fr2, "Hard-scatter jet efficiency", "Ratio to ITk", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2, ndiv_y=505)
    fr2.GetXaxis().SetNdivisions(505)
    fr2.GetYaxis().ChangeLabel(1, -1, 0)          # the 0 would sit on the x axis's first label
    one = S.keep(R.TLine(0.80, 1.0, 1.0, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
    for (col, name, short, colour, mk, ls), (_, _, ratio) in zip(TAGGERS[1:], rows[1:]):
        S.style(S.graph(EFF, ratio), colour, mk, ls, 2, 1.1).Draw("LP SAME")
S.save(cv, f"{args.out}_roc")

# ── Fig 3.21 analogue ─────────────────────────────────────────────────────────
EB = np.array([2.4, 2.8, 3.2, 3.6, 4.0])
rel, nbin = {}, []
for col, *_ in TAGGERS:
    v = []
    for b in range(len(EB) - 1):
        m = sel(30.0, 50.0, EB[b], EB[b + 1])
        hs, pu = J[col][m & (J["hs"] == 1)], J[col][m & (J["pu"] == 1)]
        hi, pi = J["rpt_itk"][m & (J["hs"] == 1)], J["rpt_itk"][m & (J["pu"] == 1)]
        pe = at_eff(hs, pu, np.array([0.85]))[0][0]
        pr = at_eff(hi, pi, np.array([0.85]))[0][0]
        v.append(pe / pr if pr > 0 else np.nan)
        if col == "rpt_itk": nbin.append((len(hs), len(pu)))
    rel[col] = v
cv = S.keep(R.TCanvas("eta", "", 800, 600))
ymax = max(2.8, 1.15 * np.nanmax([x for v in rel.values() for x in v]))
fr = S.keep(cv.DrawFrame(2.4, 0.0, 4.0, ymax))
S.axes(fr, "Jet |#eta|", "Pileup jet rate relative to ITk only")
fr.GetYaxis().ChangeLabel(1, -1, 0)
leg = S.legend(0.18, 0.51, 0.74, 0.66, size=0.036, ncol=2)
for col, name, short, colour, mk, ls in TAGGERS:
    h = S.keep(R.TH1D(f"h_eta_{col}", "", len(EB) - 1, EB))
    for b, x in enumerate(rel[col]): h.SetBinContent(b + 1, x)
    S.style(h, colour, 0, ls, 3 if col in ("rpt_itk", "rpt_t0_truth") else 2); h.Draw("HIST SAME")
    leg.AddEntry(h, short, "L")
leg.Draw()
S.header(0.20, 0.87, [ENERGY + ", " + args.label, f"30 < p_{{T}} < 50 GeV, 85% hard-scatter efficiency, t_{{0}} gate {gate:g}#sigma", CLEAN_TXT[clean]])
cv.RedrawAxis()
S.save(cv, f"{args.out}_eta")
print("\nrelative pileup-jet rate at 85%, 30-50 GeV, |eta| bins " + " ".join(f"{EB[b]:.1f}-{EB[b + 1]:.1f}" for b in range(len(EB) - 1)))
for col, name, *_ in TAGGERS:
    print(f"  {plain(name):40s}" + "".join(f"{x:8.2f}" for x in rel[col]))
print(f"  {'(hard-scatter / pileup jets per bin)':40s}" + "".join(f"  {h}/{p}" for h, p in nbin))

# ── the 85% working points, for the TDR preselection ──────────────────────────
wp = {}
print("\nR_pT working points at 85% hard-scatter efficiency (threshold / achieved HS eff / pileup eff):")
for col, name, *_ in TAGGERS:
    wp[col] = {}
    line = f"  {plain(name):40s}"
    for region, (elo, ehi) in {"central": (0.0, 2.4), "forward": (2.4, 4.0)}.items():
        for lab, lo, hi in SLICES:
            m = sel(lo, hi, elo, ehi)
            hs, pu = J[col][m & (J["hs"] == 1)], J[col][m & (J["pu"] == 1)]
            thr, ehs, epu = working_point(hs, pu, 0.85)
            key = f"{region}_{int(lo)}_{'inf' if hi > 1e8 else int(hi)}"
            wp[col][key] = dict(threshold=thr, hs_eff=ehs, pu_eff=epu, n_hs=int(len(hs)), n_pu=int(len(pu)))
            line += f"  {key} {thr:.3f}/{ehs:.3f}/{epu:.4f}"
    print(line)
json.dump(dict(tree=os.path.abspath(args.tree), hs_eff=0.85, gate_sigma=gate, time_cleaning=clean, working_points=wp), open(args.out + "_wp85.json", "w"), indent=1)
print("wrote", args.out + "_wp85.json")
