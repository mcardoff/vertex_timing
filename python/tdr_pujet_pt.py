#!/usr/bin/env python3
"""Where JVT's p_T range should end at <mu> = 200: the p_T of pileup jets against hard-scatter ones.

Run 2/3 JVT stops at 60 GeV because pileup jets above it are rare (Ariel Schwartzman, 2026-11-08: "it is extremely
unlikely to have a pileup jet above 60 GeV ... For run 4 this may go a little bit up but not much. You could make a
plot of PU jet pt and put a cut based on that"). This is that plot, on our samples.

Jets: util/tdr_pujet's per-jet trees (EMTopo, pT > 30 GeV, |eta| < 4.0, every event passing the vertex selection and,
for Z+jets, Z->ll with lepton-jet overlap removal -- the population JVT acts on), TDR truth labels (HS: a truth jet
> 10 GeV within dR 0.3; PU: none > 4 GeV within 0.6). Optionally also the pileup VBS jets of section 3's pair before
any pileup tagging (tdr_fwdpair --pair=one, --jvt=none trees): the pair takes the two leading jets, so its pileup legs
are the hardest pileup jets of the event.
Per sample, central (|eta| < 2.5, JVT's Run 2 range) and forward (2.5 <= |eta| < 4.0) columns:
  upper  the p_T shape of HS and PU jets (each normalised to unit area), log scale
  lower  the pileup fraction of all jets in each p_T bin
Dashed lines at 60 GeV (JVT's p_T range) and 120 GeV (fJVT's). The last bin holds the overflow.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_pujet_pt.py \\
      --jets ../condor/vbf/vbf_tdr_pujet.root ../condor/zjets/zjets_tdr_pujet.root \\
      --pairs ../condor/vbf/vbf_onefwd_tdr_fwdpair.root ../condor/zjets/zjets_onefwd_tdr_fwdpair.root \\
      --out ../figs/tdr/grid/grid_tdr_pujet_pt
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--jets", nargs=2, required=True, metavar=("VBF", "ZJETS"), help="tdr_pujet trees")
ap.add_argument("--pairs", nargs=2, metavar=("VBF", "ZJETS"), help="tdr_fwdpair --pair=one --jvt=none trees (optional)")
ap.add_argument("--out", required=True)
args = ap.parse_args()

EDGES = np.array([30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 200, 250, 300], float)
REGIONS = [("central", "|#eta| < 2.5", 0.0, 2.5), ("forward", "2.5 #leq |#eta| < 4.0", 2.5, 4.0)]
SAMPLES = [("vbf", "VBF H#rightarrowinv."), ("zjets", "Z+jets")]


def hist(pt, w):
    """weighted counts per bin, the last bin holding the overflow"""
    return np.histogram(np.minimum(pt, EDGES[-1] - 1e-6), EDGES, weights=w)[0]


def quantile(pt, w, q):
    o = np.argsort(pt); c = np.cumsum(w[o]) / w.sum()
    return pt[o][np.searchsorted(c, q)]


R = S.setup()
for (key, lab), jpath, ppath in zip(SAMPLES, args.jets, args.pairs or (None, None)):
    J = uproot.open(jpath)["jets"].arrays(["weight", "pt", "eta", "hs", "pu"], library="np")
    w, pt, ae = J["weight"].astype(float), J["pt"].astype(float), np.abs(J["eta"])
    P = None
    if ppath:
        E = uproot.open(ppath)["events"].arrays(library="np")
        legs = [(E[f"j{k}_pt"], np.abs(E[f"j{k}_eta"]), E[f"j{k}_pu"] == 1, E["weight"].astype(float)) for k in (1, 2)]
        P = {k: np.concatenate([l[i] for l in legs]) for i, k in enumerate(("pt", "ae", "pu", "w"))}
    print(f"\n{lab.replace('#rightarrow', '->')}: {len(pt):,} jets")
    cv = S.keep(R.TCanvas(f"pt_{key}", "", 1400, 800))
    for c, (rk, rlab, lo, hi) in enumerate(REGIONS):
        r = (ae >= lo) & (ae < hi)
        hs, pu = r & (J["hs"] == 1), r & (J["pu"] == 1)
        n_all, n_hs, n_pu = hist(pt[r], w[r]), hist(pt[hs], w[hs]), hist(pt[pu], w[pu])
        frac = n_pu / np.maximum(n_all, 1e-12)
        neff = n_all ** 2 / np.maximum(hist(pt[r], w[r] ** 2), 1e-12)
        efrac = np.sqrt(frac * (1 - frac) / np.maximum(neff, 1))
        above = lambda m, x: 100 * w[m & (pt > x)].sum() / w[m].sum()
        print(f"  {rk:8s} {r.sum():8,} jets, PU {100 * w[pu].sum() / w[r].sum():4.1f}%: PU jets above 60 GeV {above(pu, 60):5.1f}%, "
              f"above 120 GeV {above(pu, 120):4.1f}%; PU-jet pT quantiles 90/95/99%: "
              + " / ".join(f"{quantile(pt[pu], w[pu], q):.0f}" for q in (0.90, 0.95, 0.99)) + " GeV")
        print(f"           PU fraction of all jets per pT bin: "
              + ", ".join(f"{EDGES[i]:.0f}-{EDGES[i + 1]:.0f}{'+' if i == len(EDGES) - 2 else ''}: {100 * frac[i]:.1f}%" for i in range(len(frac))))
        sh = {}
        x1 = c * 0.5
        cv.cd(); S.pad(f"u_{key}_{rk}", x1, 0.36, x1 + 0.5, 1.0, top=0.06, bottom=0.02, left=0.15, right=0.04)
        fr = S.keep(R.gPad.DrawFrame(EDGES[0], 2e-5, EDGES[-1], 30))
        R.gPad.SetLogy()
        S.axes(fr, "", "Fraction of jets / bin", size=0.055, ytoff=1.25)
        fr.GetXaxis().SetLabelSize(0)
        x = 0.5 * (EDGES[:-1] + EDGES[1:])
        leg = S.legend(0.50, 0.56, 0.95, 0.74, size=0.042)
        curves = [("HS jets", n_hs, "kP10Green", 1), ("PU jets", n_pu, "kP10Red", 1)]
        if P is not None and rk == "forward":
            m = P["pu"] & (P["ae"] > 2.4)
            curves.append(("PU VBS jets of the pair", hist(P["pt"][m], P["w"][m]), "kP10Red", 2))
            print(f"           PU VBS legs at 2.4 < |eta| < 3.8 (section-3 pair, no JVT): {m.sum():,} MC; above 60 GeV "
                  f"{100 * P['w'][m & (P['pt'] > 60)].sum() / P['w'][m].sum():.1f}%, above 120 GeV "
                  f"{100 * P['w'][m & (P['pt'] > 120)].sum() / P['w'][m].sum():.1f}%; quantiles 50/90/95%: "
                  + " / ".join(f"{quantile(P['pt'][m], P['w'][m], q):.0f}" for q in (0.5, 0.9, 0.95)) + " GeV")
        for name, n, col, ls in curves:
            h = S.keep(R.TH1D(f"h_{key}_{rk}_{name}", "", len(EDGES) - 1, EDGES))
            for i, v in enumerate(n / n.sum()): h.SetBinContent(i + 1, v)
            h.SetLineColor(S.colour(col)); h.SetLineStyle(ls); h.SetLineWidth(3); h.Draw("HIST SAME")
            leg.AddEntry(h, name, "L")
        for v in (60, 120):
            ln = S.keep(R.TLine(v, 2e-5, v, 0.5)); ln.SetLineStyle(2); ln.SetLineColor(R.kGray + 2); ln.Draw()   # stops below the header
        leg.Draw()
        S.header(0.19, 0.86, ["#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200", f"{lab}, jets p_{{T}} > 30 GeV, {rlab}"], size=0.045)
        R.gPad.RedrawAxis()
        cv.cd(); S.pad(f"l_{key}_{rk}", x1, 0.0, x1 + 0.5, 0.36, top=0.03, bottom=0.30, left=0.15, right=0.04)
        top = min(1.0, 1.25 * np.max(frac + efrac))
        fr2 = S.keep(R.gPad.DrawFrame(EDGES[0], 0.0, EDGES[-1], top))
        S.axes(fr2, "Jet p_{T} [GeV]   (last bin: overflow)", "Pileup fraction", size=0.095, ytoff=0.72, xtoff=1.15, ndiv_y=504)
        g = S.style(S.graph(x, frac, efrac), "kP10Red", 20, 1, 2, 1.0); g.Draw("P SAME")
        for v in (60, 120):
            ln = S.keep(R.TLine(v, 0.0, v, top)); ln.SetLineStyle(2); ln.SetLineColor(R.kGray + 2); ln.Draw()
        R.gPad.RedrawAxis()
    S.save(cv, f"{args.out}_{key}")
