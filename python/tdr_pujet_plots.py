#!/usr/bin/env python3
"""TDR Sec. 3.3.3 figures from util/tdr_pujet's per-jet tree, in the ATLAS style (python/tdr_root_style.py).

  <out>_roc      TDR Fig 3.20: pileup jet rejection (1 / pileup efficiency, linear axis) at the sampled
                 hard-scatter efficiencies 0.800-0.975, 2.4 < |eta| < 4.0, for 30-50 GeV (left) and
                 > 50 GeV (right) jets; the t0-gated taggers (HGTD, TRKPTZ, WAVeS, truth) against
                 ITk only, with the ratio to ITk underneath. The tree's self-tagging and combined
                 columns are not drawn (user, 2026-11-03: too many series; combined = HGTD wherever
                 Athena has a t0)
  <out>_eta      TDR Fig 3.21: pileup jet rate relative to ITk only, against |eta|, 30-50 GeV, each tagger
                 at ONE 85% hard-scatter working point set over 2.4 < |eta| < 4.0 (the TDR's legend; the
                 user's choice, 2026-11-03) and applied in every |eta| bin -- not a working point per bin,
                 which forces the 3.6-4.0 bin to 85% right at the timed taggers' efficiency ceiling
  <out>_eta_tdr  the TDR's own Fig 3.21 method as a reference point: its "combined" tagger (t0 gate where
                 the vertex t0 exists -- here Athena's RecoVtx_time -- self-tagging otherwise), on the same
                 single working point, with statistical errors, against the TDR's published values (read
                 off its Fig 3.21 by eye) and ITk only
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
           ("rpt_t0_truth", "Truth t_{0}", "Truth t_{0}", S.COL["truth"], 25, 3)]
EFF = np.round(np.arange(0.800, 0.9751, 0.025), 3)    # the sampled working points
plain = lambda s: s.replace("_{0}", "0").replace("#", "")

t = uproot.open(args.tree)
J = t["jets"].arrays(library="np")
meta = t["meta"].arrays(library="np")
gate = float(np.unique(meta["gate_sigma"])[0])
clean = str(meta["time_cleaning"][0]) if "time_cleaning" in meta else "off"   # trees before 2026-11-03: uncleaned
match = str(meta["rpt_match"][0]) if "rpt_match" in meta else "dr02"         # trees before the switch: the cone
CLEAN_TXT = {"all": None, "t0": "Cleaned t_{0} only", "off": "Track times as in the ntuple"}   # the default (all) unlabelled
MATCH_TXT = {"ghost": None, "dr02": "R_{p_{T}} tracks within #DeltaR < 0.2"}   # the default (ghost) unlabelled
clean_line = [t for t in (CLEAN_TXT[clean], MATCH_TXT[match]) if t]   # the non-default modes, one line each
gate_txt = "" if gate == 2.0 else f", t_{{0}} gate {gate:g}#sigma"     # the TDR's 2 sigma is the standard: unlabelled

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
print(f"{args.tree}: gate {gate:g} sigma, track-time cleaning {clean}, R_pT matching {match}\n\nrejection at matched hard-scatter efficiency, 2.4 < |eta| < 4.0 (ratio to ITk in brackets)")
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
    # y range: the tallest curve just under the header (whose depth depends on its line count), a round
    # number, never clipping an error bar
    n_lines = 2 + len(clean_line)
    frac = (0.87 - n_lines * 1.2 * TS_M - 0.012 - 0.02) / 0.936
    need = np.nanmax([np.nanmax(r[0] + r[1]) for r in rows]) / frac   # error bars included: they stay under the header
    step = 50.0 if need > 200 else 25.0 if need > 100 else 10.0
    ymax = step * np.ceil(need / step)
    rmax = np.ceil(1.1 * np.nanmax([np.nanmax(r[2]) for r in rows]) / 0.5) * 0.5

    cv.cd()
    p1 = S.pad(f"main{c}", 0.5 * c, RAT, 0.5 * c + 0.5, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
    fr = S.keep(p1.DrawFrame(0.80, 0.0, 1.0, ymax))
    S.axes(fr, "", "Pileup jet rejection", size=TS_M, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0); fr.GetXaxis().SetNdivisions(505)
    leg = S.legend(0.72, 0.66, 0.94, 0.93, size=0.7 * TS_M, margin=0.25)   # top-right corner, short names
    for (col, name, short, colour, mk, ls), (rej, rerr, _) in zip(TAGGERS, rows):
        g = S.style(S.graph(EFF, rej), colour, mk, ls, 2, 1.3); g.Draw("L SAME")      # the line in its style
        S.style(S.graph(EFF, rej, rerr), colour, mk, 1, 2, 1.3).Draw("P SAME")         # solid error bars
        leg.AddEntry(g, short, "LP")
    leg.Draw()
    S.header(0.20, 0.87, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}{gate_txt}", f"{lab}, 2.4 < |#eta| < 4.0"] + clean_line, size=TS_M)

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
rel, effhs, nbin = {}, {}, []
mf = sel(30.0, 50.0, 2.4, 4.0)
thr = {col: working_point(J[col][mf & (J["hs"] == 1)], J[col][mf & (J["pu"] == 1)], 0.85)[0] for col, *_ in TAGGERS}
for col, *_ in TAGGERS:
    v, e = [], []
    for b in range(len(EB) - 1):
        m = sel(30.0, 50.0, EB[b], EB[b + 1])
        H, P = m & (J["hs"] == 1), m & (J["pu"] == 1)
        pe = np.mean(J[col][P] >= thr[col])
        pr = np.mean(J["rpt_itk"][P] >= thr["rpt_itk"])
        v.append(pe / pr if pr > 0 else np.nan)
        e.append(np.mean(J[col][H] >= thr[col]))
        if col == "rpt_itk": nbin.append((int(H.sum()), int(P.sum())))
    rel[col], effhs[col] = v, e
cv = S.keep(R.TCanvas("eta", "", 800, 600))
fr = S.keep(cv.DrawFrame(2.4, 0.0, 4.0, 1.6))                         # the TDR's y range
S.axes(fr, "Jet |#eta|", "Pileup jet rate relative to ITk only")
fr.GetYaxis().ChangeLabel(1, -1, 0)
leg = S.legend(0.20, 0.19, 0.93, 0.33 + (0.06 if clean_line else 0), size=0.04, ncol=3)   # along the bottom
if clean_line: leg.SetHeader(", ".join(clean_line))   # the non-default modes: no room in the header on the TDR's y range
for col, name, short, colour, mk, ls in TAGGERS:
    h = S.keep(R.TH1D(f"h_eta_{col}", "", len(EB) - 1, EB))
    for b, x in enumerate(rel[col]): h.SetBinContent(b + 1, x)
    S.style(h, colour, 0, ls, 3 if col in ("rpt_itk", "rpt_t0_truth") else 2); h.Draw("HIST SAME")
    leg.AddEntry(h, short, "L")
leg.Draw()
S.header(0.20, 0.87, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}{gate_txt}",
                      "2.4 < |#eta| < 4.0, 30 < p_{T} < 50 GeV, 85% HS efficiency"])
cv.RedrawAxis()
S.save(cv, f"{args.out}_eta")

# ── Fig 3.21 the TDR's way: its combined tagger against its published values ──
# ATLAS-TDR-031 Fig 3.21 (PDF p. 69), the ITk+HGTD (combined) curve, read off the figure by eye (+- 0.01)
TDR_FIG321 = [0.87, 0.70, 0.72, 0.57]
tC = working_point(J["rpt_comb"][mf & (J["hs"] == 1)], J["rpt_comb"][mf & (J["pu"] == 1)], 0.85)[0]
rc, rce, nc = [], [], []
for b in range(len(EB) - 1):
    m = sel(30.0, 50.0, EB[b], EB[b + 1])
    P = m & (J["pu"] == 1)
    kc, ki = int(np.sum(J["rpt_comb"][P] >= tC)), int(np.sum(J["rpt_itk"][P] >= thr["rpt_itk"]))
    r = kc / ki if ki else np.nan
    rc.append(r); nc.append((kc, ki))
    rce.append(r * np.sqrt(1.0 / max(kc, 1) + 1.0 / max(ki, 1)))   # Poisson, the two counts taken as independent: conservative
cv = S.keep(R.TCanvas("eta_tdr", "", 800, 600))
fr = S.keep(cv.DrawFrame(2.4, 0.0, 4.0, 1.6))
S.axes(fr, "Jet |#eta|", "Pileup jet rate relative to ITk only")
fr.GetYaxis().ChangeLabel(1, -1, 0)
h_itk = S.keep(R.TH1D("h_eta_tdr_itk", "", len(EB) - 1, EB))
for b in range(len(EB) - 1): h_itk.SetBinContent(b + 1, 1.0)
S.style(h_itk, S.COL["itk"], 0, 2, 3); h_itk.Draw("HIST SAME")
h_tdr = S.keep(R.TH1D("h_eta_tdr_pub", "", len(EB) - 1, EB))
for b, x in enumerate(TDR_FIG321): h_tdr.SetBinContent(b + 1, x)
S.style(h_tdr, "kP10Ash", 0, 1, 3); h_tdr.Draw("HIST SAME")
h_c = S.keep(R.TH1D("h_eta_tdr_comb", "", len(EB) - 1, EB))
for b, x in enumerate(rc): h_c.SetBinContent(b + 1, x)
S.style(h_c, "kP10Violet", 0, 1, 3); h_c.Draw("HIST SAME")
g_c = S.style(S.graph(0.5 * (EB[:-1] + EB[1:]), rc, rce), "kP10Violet", 20, 1, 2, 1.2); g_c.Draw("P SAME")
leg = S.legend(0.20, 0.19, 0.93, 0.35 + (0.055 if clean_line else 0), size=0.038)
if clean_line: leg.SetHeader(", ".join(clean_line))
leg.AddEntry(h_itk, "ITk only", "L")
leg.AddEntry(h_tdr, "ATLAS-TDR-031, Fig. 3.21 (combined)", "L")
leg.AddEntry(g_c, "This study, TDR method: Athena t_{0}, else self-tagging", "LP")
leg.Draw()
# the header starts right of the first bin, whose error bar reaches ~1.3
S.header(0.40, 0.87, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}{gate_txt}",
                      "2.4 < |#eta| < 4.0, 30 < p_{T} < 50 GeV", "85% HS efficiency"], size=0.045)
cv.RedrawAxis()
S.save(cv, f"{args.out}_eta_tdr")
okH = np.mean(J["ok_hgtd"][mf])
print(f"\nFig 3.21 the TDR's way (combined, one 85% working point over 2.4-4.0): threshold {tC:.3f}; Athena t0 exists for {100 * okH:.1f}% of these jets")
print("  TDR published  " + "  ".join(f"{x:5.2f}" for x in TDR_FIG321))
print("  this study     " + "  ".join(f"{x:5.2f}+-{e:4.2f}" for x, e in zip(rc, rce)) + "   (pileup jets passing, combined / ITk: " + " ".join(f"{a}/{b}" for a, b in nc) + ")")
print("\nrelative pileup-jet rate, 30-50 GeV, ONE 85% working point over 2.4-4.0, |eta| bins " + " ".join(f"{EB[b]:.1f}-{EB[b + 1]:.1f}" for b in range(len(EB) - 1)) + "   (HS efficiency per bin)")
for col, name, *_ in TAGGERS:
    print(f"  {plain(name):40s}" + "".join(f"{x:8.2f}" for x in rel[col]) + "   (" + " ".join(f"{e:.3f}" for e in effhs[col]) + ")")
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
json.dump(dict(tree=os.path.abspath(args.tree), hs_eff=0.85, gate_sigma=gate, time_cleaning=clean, rpt_match=match, working_points=wp), open(args.out + "_wp85.json", "w"), indent=1)
print("wrote", args.out + "_wp85.json")
