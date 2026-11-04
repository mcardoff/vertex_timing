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
args = ap.parse_args()
R = S.setup()

CASES = [("rpt_itk", "ITk only", S.COL["itk"], 1), ("rpt_t0_hgtd", "HGTD t_{0}", S.COL["hgtd"], 1),
         ("rpt_t0_trkptz", "TRKPTZ t_{0}", S.COL["trkptz"], 1), ("rpt_t0_waves", "WAVeS t_{0}", S.COL["waves"], 1),
         ("rpt_t0_truth", "Truth t_{0}", S.COL["truth"], 3)]
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
