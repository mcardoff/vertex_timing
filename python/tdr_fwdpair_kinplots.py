#!/usr/bin/env python3
"""What remains after a VBF H->inv vs Z+jets BDT on the 2022 PUB note's loose selection
(tdr_fwdpair_kinbdt.py): the m_jj stack by the truth identity of the pair (tdr_stack.py, the same look as the
TDR-preselection stacks) at the selection and after BDT cuts at fixed signal efficiencies, plus a composition
table. Both jets are in 2.4 < |eta| < 3.8 by construction, so every pair is forward-forward and every pileup
leg is inside HGTD: the category and "where is the pileup leg" stacks would be trivial and are not drawn.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_kinplots.py --sig <vbf>_tdr_fwdpair.root \\
      --bkg <zjets>_tdr_fwdpair.root --bdt <x>.npz --out ../figs/tdr/grid/fwdpair/grid_fwdpair_kin
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S
import tdr_stack as ST

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--bdt", required=True)
ap.add_argument("--effs", default="0.8,0.5,0.2")
ap.add_argument("--stack-edges", default="0,1000,2000,4000")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R_ = S.setup()
EDGES = np.array([float(x) for x in args.stack_edges.split(",")])
ENERGY = "#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200"
SEL = "j_{1}, j_{2}: 2.4 < |#eta| < 3.8, opposite #eta, > 50 / 30 GeV"

D = np.load(args.bdt)
samples = [("sig", "VBF H#rightarrowinv.", uproot.open(args.sig)["events"].arrays(library="np"), D["sig_score"]),
           ("bkg", "Z+jets", uproot.open(args.bkg)["events"].arrays(library="np"), D["bkg_score"])]
for key, lab, E, sc in samples:
    assert len(sc) == len(E["weight"]) and np.all(np.isfinite(sc)), f"{args.bdt}: {key} scores do not match the tree"
E0, s0 = samples[0][2], samples[0][3]
o = np.argsort(-s0); cs = np.cumsum(E0["weight"][o]) / E0["weight"].sum()
cuts = [(None, None)] + [(e, s0[o][np.searchsorted(cs, e)]) for e in (float(x) for x in args.effs.split(","))]

print(f"{'':28s} {'events':>7s} {'HS + HS':>15s} {'HS + PU':>15s} {'PU + PU':>15s} {'neither':>15s}")
for key, lab, E, sc in samples:
    w = E["weight"].astype(float)
    nhs, npu = E["j1_hs"] + E["j2_hs"], E["j1_pu"] + E["j2_pu"]
    groups = [nhs == 2, (nhs == 1) & (npu == 1), npu == 2, (nhs + npu) < 2]
    for eff, thr in cuts:
        ok = np.ones(len(w), bool) if thr is None else sc >= thr
        tag = "" if eff is None else f"_bdt{round(100 * eff)}"
        line = "the note's selection" if eff is None else f"+ BDT (j_{{1}}, j_{{2}} kinematics) at #varepsilon_{{S}} = {100 * eff:.0f}%"
        name = "selection" if eff is None else f"BDT eps_S = {eff:.0%}"
        print(f"{lab.replace('#rightarrow', '->') + ', ' + name:28s} {int(ok.sum()):7d} "
              + " ".join(f"{100 * w[ok & g].sum() / w[ok].sum():6.1f}% ({int((ok & g).sum()):5d})" for g in groups))
        if ok.sum() == 0: continue
        ST.mjj_stack(R_, f"{args.out}_truth_mjj_stack{tag}_{key}", f"kin_truth{tag}_{key}", w, ok, E["mjj"],
                     [(t, m, c) for (t, c), m in zip(ST.IDENT, groups)], EDGES,
                     [f"{ENERGY}, {lab}", SEL, line], f"{lab.replace('#rightarrow', '->')} truth stack{tag}")
