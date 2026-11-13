#!/usr/bin/env python3
"""R_pT of the VBS jets HGTD can see, without and with a t0 gate, split by the pair's eta configuration -- the
lead-in from the pileup-jet section to the analysis section of slides/slides_tdr_timing.tex (the user's request,
2026-11-04: like the VBS-region R_pT plot, python/rpt_region_dists.py, but divided by eta configuration and NOT
by truth identity, so every panel is well filled).

Events: util/tdr_fwdpair.cxx --pair=one trees (the two leading jets with |eta| < 3.8, >= 1 of them in
2.4 < |eta| < 3.8, opposite eta, > 50 / 30 GeV). Jets counted:
  central-forward pairs (cat 1)  the forward jet only (the central one has no HGTD time)
  forward-forward pairs (cat 2)  both jets
R_pT: tdr_fwdpair's TDR values -- ghost-associated tracks from PV0 / jet pT; with a t0, timed tracks only within
2 sigma of it (the TDR gate, cleaned track times). One 2x2 canvas: rows = configuration, columns = VBF H->inv /
Z+jets; each panel the ITk-only distribution (filled) against the HGTD / TRKPTZ / WAVeS / truth t0 gates, unit
area (weighted), with the share of jets at R_pT = 0 in the legend. Truth labels only annotate (the pileup share).

Also <out>_roc (the user's choice for the slide, 2026-11-04: ROC curves instead of the distributions): the same
2x2, each panel the pileup-jet rejection against the hard-scatter-jet efficiency of those jets (truth labels per
jet; the event split stays by eta configuration), as section 2's ROC (python/tdr_pujet_plots.py: unweighted,
sampled at 0.800 ... 0.975, linear interpolation between thresholds, binomial errors), ratio to ITk below.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_rpt.py --sig <vbf>_onefwd_tdr_fwdpair.root \\
      --bkg <zjets>_onefwd_tdr_fwdpair.root --out ../figs/tdr/grid/onefwd/grid_onefwd_rpt
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--t0s", default=None, help="comma list of the t0s to draw (hgtd, trkptz, waves, truth); default all. "
                "The Global Performance talk shows WAVeS against Athena only (user, 2026-11-11): --t0s hgtd,waves,truth")
args = ap.parse_args()
R = S.setup()

CASES = [("rpt_itk", "ITk only", S.COL["itk"], 1), ("rpt_t0_hgtd", "HGTD t_{0}", S.COL["hgtd"], 1),
         ("rpt_t0_trkptz", "TRKPTZ t_{0}", S.COL["trkptz"], 1), ("rpt_t0_waves", "WAVeS t_{0}", S.COL["waves"], 1),
         ("rpt_t0_truth", "Truth t_{0}", S.COL["truth"], 3)]
if args.t0s:   # ITk (the reference) always stays
    CASES = [c for c in CASES if c[0] == "rpt_itk" or c[0].replace("rpt_t0_", "") in args.t0s.split(",")]
CONF = [(1, "Central-forward pairs: the forward jet"), (2, "Forward-forward pairs: both jets")]
SAMPLES = [("VBF H#rightarrowinv.", args.sig), ("Z+jets", args.bkg)]
EDGES = np.linspace(0.0, 1.2, 61)


def jets(E, cat):
    """the counted jets of events in configuration `cat`: {column: values}, plus weights and PU / HS flags"""
    out = {k: [] for k, *_ in CASES}
    w, pu, hs = [], [], []
    for k in (1, 2):
        m = (E["cat"] == cat) & ((np.abs(E[f"j{k}_eta"]) > 2.4) if cat == 1 else True)
        for c, *_ in CASES: out[c].append(E[f"j{k}_{c}"][m])
        w.append(E["weight"][m].astype(float)); pu.append(E[f"j{k}_pu"][m] == 1); hs.append(E[f"j{k}_hs"][m] == 1)
    return {c: np.concatenate(v) for c, v in out.items()}, np.concatenate(w), np.concatenate(pu), np.concatenate(hs)


cv = S.keep(R.TCanvas("rpt2x2", "", 1600, 1200))
print("R_pT = 0 share of the counted jets (weighted)")
for row, (cat, ctxt) in enumerate(CONF):
    for col, (slab, path) in enumerate(SAMPLES):
        E = uproot.open(path)["events"].arrays(library="np")
        V, w, pu, hs = jets(E, cat)
        cv.cd()
        p = S.pad(f"p{row}{col}", 0.5 * col, 0.5 * (1 - row), 0.5 * col + 0.5, 0.5 * (2 - row), top=0.05, bottom=0.15, left=0.15, right=0.04)
        p.SetLogy()
        hs_ = []
        for c, name, colour, ls in CASES:
            h = S.keep(R.TH1D(f"h_{row}{col}_{c}", "", len(EDGES) - 1, EDGES))
            for v_, ww in zip(np.clip(V[c], 0.0, EDGES[-1] - 1e-6), w): h.Fill(v_, ww)
            h.Scale(1.0 / w.sum())
            hs_.append((c, name, colour, ls, h, w[V[c] <= 0].sum() / w.sum()))
        top = max(h.GetMaximum() for *_, h, _ in hs_)
        fr = S.keep(p.DrawFrame(0.0, 2e-5, 1.2, 5000.0 * top))   # headroom: the curves stay under the header and legend
        S.axes(fr, "R_{p_{T}}   (last bin: overflow)", "Fraction of jets / 0.02", size=0.05, ytoff=1.4, xtoff=1.25)
        for c, name, colour, ls, h, z in hs_:
            if c == "rpt_itk":
                h.SetFillColor(R.TColor.GetColor("#D9D9D9")); h.SetLineColor(R.kGray + 2); h.SetLineWidth(1); h.Draw("HIST SAME")
            else:
                h.SetLineColor(S.colour(colour)); h.SetLineWidth(3 if ls == 1 else 2); h.SetLineStyle(ls); h.Draw("HIST SAME")
        leg = S.legend(0.52, 0.48, 0.94, 0.76, size=0.042, margin=0.22)
        leg.SetHeader("R_{p_{T}} = 0 for")
        for c, name, colour, ls, h, z in hs_:
            leg.AddEntry(h, f"{name}: {100 * z:.0f}%", "F" if c == "rpt_itk" else "L")
        leg.Draw()
        S.header(0.19, 0.88, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {slab}", ctxt,
                              f"{len(w):,} jets, {100 * w[pu].sum() / w.sum():.0f}% pileup"], size=0.045)
        R.gPad.RedrawAxis()
        print(f"  {slab.replace('#rightarrow', '->'):14s} {ctxt:40s} {len(w):7d} jets ({100 * w[pu].sum() / w.sum():4.1f}% PU, "
              f"{100 * w[hs].sum() / w.sum():4.1f}% HS): " + ", ".join(f"{name.replace('_{0}', '0').replace('#', '')} {100 * z:.1f}%" for c, name, *_, z in hs_))
S.save(cv, args.out)


# ── the ROC version ──────────────────────────────────────────────────────────
EFF = np.round(np.arange(0.800, 0.9751, 0.025), 3)
SHORT = {"rpt_itk": "ITk only", "rpt_t0_hgtd": "HGTD t_{0}", "rpt_t0_trkptz": "TRKPTZ t_{0}",
         "rpt_t0_waves": "WAVeS t_{0}", "rpt_t0_truth": "Truth t_{0}"}
MARK = {"rpt_itk": 24, "rpt_t0_hgtd": 23, "rpt_t0_trkptz": 21, "rpt_t0_waves": 20, "rpt_t0_truth": 25}


def curve(hs, pu):
    """python/tdr_pujet_plots.py's: per distinct HS efficiency, the lowest PU efficiency reaching it"""
    th = np.unique(np.concatenate([hs, pu]))
    e_hs = 1.0 - np.searchsorted(np.sort(hs), th, side="left") / len(hs)
    e_pu = 1.0 - np.searchsorted(np.sort(pu), th, side="left") / len(pu)
    u, inv = np.unique(e_hs, return_inverse=True)
    pmin = np.full(len(u), np.inf); np.minimum.at(pmin, inv, e_pu)
    return u, pmin


def at_eff(hs, pu):
    e_hs, e_pu = curve(hs, pu)
    pe = np.interp(EFF, e_hs, e_pu)
    return pe, np.sqrt(np.clip(pe * (1 - pe), 0, None) / len(pu))


cr = S.keep(R.TCanvas("roc2x2", "", 1600, 1300))
MAIN_, RAT_ = 0.68, 0.32
print("\nROC: pileup-jet rejection at hard-scatter efficiency 0.85 / 0.90 (ratio to ITk)")
for row, (cat, ctxt) in enumerate(CONF):
    for col, (slab, path) in enumerate(SAMPLES):
        E = uproot.open(path)["events"].arrays(library="np")
        V, w, pu, hs = jets(E, cat)
        rows_ = []
        for c, name, colour, ls in CASES:
            pe, perr = at_eff(V[c][hs], V[c][pu])
            with np.errstate(divide="ignore", invalid="ignore"):
                rej, rerr = np.where(pe > 0, 1.0 / pe, np.nan), np.where(pe > 0, perr / pe ** 2, np.nan)
            rows_.append((c, colour, ls, rej, rerr))
        ref = rows_[0][3]
        x0, y0 = 0.5 * col, 0.5 * (1 - row)
        cr.cd()
        p1 = S.pad(f"r{row}{col}m", x0, y0 + 0.5 * RAT_, x0 + 0.5, y0 + 0.5, top=0.06, bottom=0.02, left=0.15, right=0.04)
        ymax = max(np.nanmax(r[3] + np.nan_to_num(r[4])) for r in rows_) * 1.9
        ymax = float(np.ceil(ymax / (50 if ymax > 200 else 10 if ymax > 40 else 5)) * (50 if ymax > 200 else 10 if ymax > 40 else 5))
        fr = S.keep(p1.DrawFrame(0.80, 0.0, 1.0, ymax))
        S.axes(fr, "", "Pileup jet rejection", size=0.065, ytoff=1.1)
        fr.GetXaxis().SetLabelSize(0); fr.GetXaxis().SetNdivisions(505); fr.GetYaxis().ChangeLabel(1, -1, 0)
        leg = S.legend(0.70, 0.40, 0.95, 0.90, size=0.055, margin=0.25)
        for c, colour, ls, rej, rerr in rows_:
            ok = np.isfinite(rej)
            g = S.style(S.graph(EFF[ok], rej[ok]), colour, MARK[c], ls, 2, 1.2); g.Draw("L SAME")
            S.style(S.graph(EFF[ok], rej[ok], rerr[ok]), colour, MARK[c], 1, 2, 1.2).Draw("P SAME")
            leg.AddEntry(g, SHORT[c], "LP")
        leg.Draw()
        S.header(0.20, 0.85, [f"{slab}, #LT#mu#GT = 200", ctxt, f"{int(hs.sum()):,} HS / {int(pu.sum()):,} PU jets"], size=0.06)
        S.end_pad()
        cr.cd()
        p2 = S.pad(f"r{row}{col}r", x0, y0, x0 + 0.5, y0 + 0.5 * RAT_, top=0.03, bottom=0.36, left=0.15, right=0.04)
        ratios = [r[3] / ref for r in rows_[1:]]
        rmax = float(np.ceil(1.1 * np.nanmax([np.nanmax(q) for q in ratios]) / 0.5) * 0.5)
        fr2 = S.keep(p2.DrawFrame(0.80, 0.5, 1.0, max(rmax, 1.5)))
        S.axes(fr2, "Hard-scatter jet efficiency", "Ratio to ITk", size=0.13, ytoff=0.55, xtoff=1.1, ndiv_y=505)
        fr2.GetXaxis().SetNdivisions(505)
        fr2.GetYaxis().ChangeLabel(1, -1, 0)          # the 0.5 would sit on the x axis's first label
        one = S.keep(R.TLine(0.80, 1.0, 1.0, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
        for (c, colour, ls, rej, rerr), q in zip(rows_[1:], ratios):
            ok = np.isfinite(q)
            S.style(S.graph(EFF[ok], q[ok]), colour, MARK[c], ls, 2, 1.0).Draw("LP SAME")
        S.end_pad()
        i85, i90 = int(np.argmin(np.abs(EFF - 0.85))), int(np.argmin(np.abs(EFF - 0.90)))
        print(f"  {slab.replace('#rightarrow', '->'):14s} {ctxt:40s} ({int(hs.sum())} HS / {int(pu.sum())} PU): " +
              ", ".join(f"{SHORT[c].replace('_{0}', '0')} {rej[i85]:.1f} ({rej[i85] / ref[i85]:.2f}) / {rej[i90]:.1f} ({rej[i90] / ref[i90]:.2f})" for c, colour, ls, rej, rerr in rows_))
S.save(cr, args.out + "_roc")

