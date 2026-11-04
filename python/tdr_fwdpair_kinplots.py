#!/usr/bin/env python3
"""What remains after a VBF H->inv vs Z+jets BDT on the 2022 PUB note's loose selection
(tdr_fwdpair_kinbdt.py): the m_jj stack by the truth identity of the pair (tdr_stack.py, the same look as the
TDR-preselection stacks) at the selection and after BDT cuts at fixed signal efficiencies, plus a composition
table. On the note's selection (--pair=ff) both jets are in 2.4 < |eta| < 3.8, so every pair is forward-forward
and every pileup leg is inside HGTD: only the truth stack is drawn. On tdr_fwdpair --pair=one (at least one
jet in HGTD, the other anywhere inside |eta| < 3.8) the category stack (central-forward / forward-forward)
and the "where is the pileup leg" stack (central / inside HGTD) are drawn too. The mode is read from the
trees' meta.

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

D = np.load(args.bdt)
samples = [("sig", "VBF H#rightarrowinv.", uproot.open(args.sig)["events"].arrays(library="np"), D["sig_score"]),
           ("bkg", "Z+jets", uproot.open(args.bkg)["events"].arrays(library="np"), D["bkg_score"])]
modes = {str(uproot.open(f)["meta"].arrays(library="np")["pair_mode"][0]) if "pair_mode" in uproot.open(f)["meta"].keys() else "ff"
         for f in (args.sig, args.bkg)}
if len(modes) != 1: raise SystemExit(f"the two samples were made with different --pair modes: {modes}")
ONE = modes.pop() == "one"
SEL = ("j_{1}, j_{2}: |#eta| < 3.8, #geq 1 in 2.4#minus3.8, opposite #eta, > 50 / 30 GeV" if ONE
       else "j_{1}, j_{2}: 2.4 < |#eta| < 3.8, opposite #eta, > 50 / 30 GeV")
CATS = [("Central-forward", 1, S.COL["cf"]), ("Forward-forward", 2, S.COL["ff"])]
REACH = [("HS + HS", "#7BCB6B"), ("PU leg at |#eta| < 2.4", "#B3231C"), ("PU leg at 2.4#minus3.8 (HGTD)", "#3F90DA"),
         ("A leg has neither label", "#C79BF5")]
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
    fwd = lambda k: np.abs(E[f"j{k}_eta"]) > 2.4
    pu_hgtd = ((E["j1_pu"] == 1) & fwd(1)) | ((E["j2_pu"] == 1) & fwd(2))       # a pileup leg HGTD sees counts first
    pu_cen = ~pu_hgtd & (((E["j1_pu"] == 1) & ~fwd(1)) | ((E["j2_pu"] == 1) & ~fwd(2)))
    reach = [(npu == 0) & (nhs == 2), pu_cen, pu_hgtd, (npu == 0) & (nhs < 2)]
    for eff, thr in cuts:
        ok = np.ones(len(w), bool) if thr is None else sc >= thr
        tag = "" if eff is None else f"_bdt{round(100 * eff)}"
        line = "before the BDT" if eff is None else f"+ BDT (j_{{1}}, j_{{2}} kinematics) at #varepsilon_{{S}} = {100 * eff:.0f}%"
        name = "selection" if eff is None else f"BDT eps_S = {eff:.0%}"
        print(f"{lab.replace('#rightarrow', '->') + ', ' + name:28s} {int(ok.sum()):7d} "
              + " ".join(f"{100 * w[ok & g].sum() / w[ok].sum():6.1f}% ({int((ok & g).sum()):5d})" for g in groups))
        if ok.sum() == 0: continue
        ST.mjj_stack(R_, f"{args.out}_truth_mjj_stack{tag}_{key}", f"kin_truth{tag}_{key}", w, ok, E["mjj"],
                     [(t, m, c) for (t, c), m in zip(ST.IDENT, groups)], EDGES,
                     [f"{ENERGY}, {lab}", SEL, line], f"{lab.replace('#rightarrow', '->')} truth stack{tag}")
        if not ONE: continue
        print(f"{'':28s}         CF {100 * w[ok & (E['cat'] == 1)].sum() / w[ok].sum():5.1f}%; PU leg central "
              f"{100 * w[ok & pu_cen].sum() / w[ok].sum():5.1f}% ({int((ok & pu_cen).sum())}), in HGTD "
              f"{100 * w[ok & pu_hgtd].sum() / w[ok].sum():5.1f}% ({int((ok & pu_hgtd).sum())})")
        ST.mjj_stack(R_, f"{args.out}_cat_mjj_stack{tag}_{key}", f"kin_cat{tag}_{key}", w, ok, E["mjj"],
                     [(t, E["cat"] == c, col) for t, c, col in CATS], EDGES,
                     [f"{ENERGY}, {lab}", SEL, line], f"{lab.replace('#rightarrow', '->')} category stack{tag}")
        ST.mjj_stack(R_, f"{args.out}_reach_mjj_stack{tag}_{key}", f"kin_reach{tag}_{key}", w, ok, E["mjj"],
                     [(t, m, c) for (t, c), m in zip(REACH, reach)], EDGES,
                     [f"{ENERGY}, {lab}", SEL, line], f"{lab.replace('#rightarrow', '->')} reach stack{tag}")
