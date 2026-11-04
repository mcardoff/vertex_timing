#!/usr/bin/env python3
"""What remains after a VBF H->inv vs Z+jets BDT on the 2022 PUB note's loose selection
(tdr_fwdpair_kinbdt.py): the m_jj stack by the truth identity of the pair (tdr_stack.py, the same look as the
TDR-preselection stacks) at the selection and after BDT cuts at fixed signal efficiencies, plus a composition
table. On the note's selection (--pair=ff) both jets are in 2.4 < |eta| < 3.8, so every pair is forward-forward
and every pileup leg is inside HGTD: only the truth stack is drawn. On tdr_fwdpair --pair=one (at least one
jet in HGTD, the other anywhere inside |eta| < 3.8) the category stack (central-forward / forward-forward)
and the "where is the pileup leg" stack (central / inside HGTD) are drawn too. The mode is read from the
trees' meta.

Timing veto (--veto): after each BDT cut, the event vetoes of the earlier timing study (util/vbs_time_veto.cxx,
python/vbs_time_veto_plot.py decisions(), reproduced here), on tdr_fwdpair's core jet time (needs trees with
core_t / core_sig / core_n, commit 53d16c3 on):
  jet vs jet  both legs timed: veto when |t_1 - t_2| / sqrt(s_1^2 + s_2^2) >= 3
  jet vs t0   per t0 source: veto when a timed leg has |t_jet - t0| / sqrt(s_jet^2 + (infl sigma_t0)^2) >= 3
s_jet = f(n) x sigma_jet, f(n) measured from the VBF file's hard-scatter forward legs (the MAD of
(t_core - t_truth) / sigma per track-count bin, vbs_time_veto_common's NBINS; the global MAD below 20 legs);
a leg is timed with >= 1 core track. Prints signal and Z+jets efficiencies of every veto at every BDT cut, and
draws the stacks after BDT + each veto named in --veto-stacks.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_kinplots.py --sig <vbf>_tdr_fwdpair.root \\
      --bkg <zjets>_tdr_fwdpair.root --bdt <x>.npz --out ../figs/tdr/grid/fwdpair/grid_fwdpair_kin
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S
import tdr_stack as ST
import tdr_fwdpair_common as FC

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--bdt", required=True)
ap.add_argument("--effs", default="0.8,0.5,0.2")
ap.add_argument("--stack-edges", default="0,1000,2000,4000")
ap.add_argument("--veto", action="store_true", help="apply the timing vetoes after the BDT cuts")
ap.add_argument("--veto-stacks", default="jj,waves", help="vetoes whose post-veto stacks are drawn")
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--summary", action="store_true", help="with --veto: <out>_veto_gain, the S/sqrt(B) factor of each veto "
                "against the BDT cut (signal efficiency 1.0 = no BDT), the Z+jets MC events left below")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R_ = S.setup()
EDGES = np.array([float(x) for x in args.stack_edges.split(",")])
ENERGY = "#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200"

plain = lambda x: x.replace('#rightarrow', '->').replace('_{0}', '0').replace('#', '')
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

VETOES = [("jj", "jet vs jet"), ("hgtd", "jet vs HGTD t_{0}"), ("trkptz", "jet vs TRKPTZ t_{0}"),
          ("waves", "jet vs WAVeS t_{0}"), ("truth", "jet vs truth t_{0}")]
calibration = FC.calibration


def vetoes(E, f, nsigma):
    return FC.vetoes(E, f, nsigma)


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

if args.veto:
    if "j1_core_n" not in samples[0][2]:
        raise SystemExit("--veto needs trees with the core jet time (tdr_fwdpair from commit 53d16c3)")
    f = calibration(samples[0][2])
    K = {key: vetoes(E, f, args.nsigma) for key, lab, E, sc in samples}
    draw = [v for v in args.veto_stacks.split(",") if v]
    print(f"\nTiming vetoes at {args.nsigma:g} sigma after each BDT cut: efficiency relative to the cut (MC events kept / before); "
          "S/sqrt(B) factor = eps_S / sqrt(eps_B)")
    print(f"  {'cut':18s} {'veto':22s} {'VBF H->inv':>22s} {'Z+jets':>22s} {'S/sqrt(B)':>10s}")
    for eff, thr in cuts:
        base = {key: (np.ones(len(E["weight"]), bool) if thr is None else sc >= thr) for key, lab, E, sc in samples}
        name = "selection" if eff is None else f"BDT eps_S = {eff:.0%}"
        for v, vlab in VETOES:
            r = {}
            for key, lab, E, sc in samples:
                w = E["weight"].astype(float); b = base[key]; k = b & K[key][v]
                r[key] = (w[k].sum() / w[b].sum(), int(k.sum()), int(b.sum()))
            fac = r["sig"][0] / np.sqrt(r["bkg"][0]) if r["bkg"][0] > 0 else np.inf
            print(f"  {name:18s} {plain(vlab):22s} {r['sig'][0]:8.3f} ({r['sig'][1]:6d}/{r['sig'][2]:6d}) "
                  f"{r['bkg'][0]:8.3f} ({r['bkg'][1]:5d}/{r['bkg'][2]:5d}) {fac:9.3f}")
        tag = "" if eff is None else f"_bdt{round(100 * eff)}"
        for v in draw:
            vlab = dict(VETOES)[v]
            for key, lab, E, sc in samples:
                w = E["weight"].astype(float)
                ok = base[key] & K[key][v]
                if ok.sum() == 0: continue
                nhs, npu = E["j1_hs"] + E["j2_hs"], E["j1_pu"] + E["j2_pu"]
                groups = [nhs == 2, (nhs == 1) & (npu == 1), npu == 2, (nhs + npu) < 2]
                line = ("before the BDT" if eff is None else f"BDT at #varepsilon_{{S}} = {100 * eff:.0f}%") + f", {vlab} veto"
                ST.mjj_stack(R_, f"{args.out}_truth_mjj_stack{tag}_veto{v}_{key}", f"kin_truth{tag}_veto{v}_{key}", w, ok, E["mjj"],
                             [(t_, m, c) for (t_, c), m in zip(ST.IDENT, groups)], EDGES,
                             [f"{ENERGY}, {lab}", SEL, line], f"{plain(lab)} truth stack{tag} + {v} veto")
                if ONE:
                    fwd = lambda k: np.abs(E[f"j{k}_eta"]) > 2.4
                    pu_hgtd = ((E["j1_pu"] == 1) & fwd(1)) | ((E["j2_pu"] == 1) & fwd(2))
                    pu_cen = ~pu_hgtd & (((E["j1_pu"] == 1) & ~fwd(1)) | ((E["j2_pu"] == 1) & ~fwd(2)))
                    reach = [(npu == 0) & (nhs == 2), pu_cen, pu_hgtd, (npu == 0) & (nhs < 2)]
                    ST.mjj_stack(R_, f"{args.out}_reach_mjj_stack{tag}_veto{v}_{key}", f"kin_reach{tag}_veto{v}_{key}", w, ok, E["mjj"],
                                 [(t_, m, c) for (t_, c), m in zip(REACH, reach)], EDGES,
                                 [f"{ENERGY}, {lab}", SEL, line], f"{plain(lab)} reach stack{tag} + {v} veto")

if args.veto and args.summary:
    # S/sqrt(B) factor of each veto, eps_S / sqrt(eps_B) relative to the BDT cut, against that cut
    EFFS = np.round(np.arange(1.0, 0.249, -0.05), 2)
    Es, ss = samples[0][2], samples[0][3]; Eb, sb = samples[1][2], samples[1][3]
    ws_, wb_ = Es["weight"].astype(float), Eb["weight"].astype(float)
    o_ = np.argsort(-ss); cs_ = np.cumsum(ws_[o_]) / ws_.sum()
    ROWS = [("waves", "WAVeS t_{0}", S.COL["waves"], 20, 1), ("trkptz", "TRKPTZ t_{0}", S.COL["trkptz"], 21, 1),
            ("hgtd", "HGTD t_{0}", S.COL["hgtd"], 23, 1), ("truth", "Truth t_{0}", S.COL["truth"], 25, 3),
            ("jj", "Jet vs jet (no t_{0})", "kP10Violet", 33, 10)]
    G, nB = {}, []
    print("\nS/sqrt(B) factor of each veto against the BDT cut (signal efficiency): " + " ".join(f"{e:.2f}" for e in EFFS))
    for v, name, colour, mk, ls in ROWS:
        y, ey = [], []
        for e in EFFS:
            thr = -np.inf if e >= 1.0 else ss[o_][np.searchsorted(cs_, e)]
            bS, bB = ss >= thr, sb >= thr
            es_ = ws_[bS & K["sig"][v]].sum() / ws_[bS].sum()
            wB = wb_[bB]; kB = K["bkg"][v][bB]
            eb_ = wB[kB].sum() / wB.sum()
            var = ((wB[kB] ** 2).sum() * (1 - eb_) ** 2 + (wB[~kB] ** 2).sum() * eb_ ** 2) / wB.sum() ** 2   # weighted binomial
            fac = es_ / np.sqrt(eb_)
            y.append(fac); ey.append(fac * 0.5 * np.sqrt(var) / eb_)
            if v == "waves": nB.append(int(bB.sum()))
        G[v] = (np.array(y), np.array(ey))
        print(f"  {plain(name):24s} " + " ".join(f"{a:.3f}" for a in y))
    MAIN_, RAT_ = 0.70, 0.30
    cv = S.keep(R_.TCanvas("veto_gain", "", 800, 800))
    cv.cd(); S.pad("vg_m", 0.0, RAT_, 1.0, 1.0, top=0.05 / MAIN_ * 0.6, bottom=0.02)
    lo, hi = 0.7, 1.55
    fr = S.keep(R_.gPad.DrawFrame(0.225, lo, 1.025, hi))
    TS = 0.05 * 600 / (800 * MAIN_)
    S.axes(fr, "", "S/#sqrt{B} with the veto / without", size=TS, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0)
    one = S.keep(R_.TLine(0.225, 1.0, 1.025, 1.0)); one.SetLineStyle(2); one.SetLineColor(R_.kGray + 2); one.Draw()
    leg = S.legend(0.19, 0.05, 0.94, 0.23, size=0.68 * TS, ncol=2, margin=0.15)   # along the bottom: empty below the curves
    for v, name, colour, mk, ls in ROWS:
        y, ey = G[v]
        g = S.style(S.graph(EFFS, y, ey), colour, mk, ls, 2, 1.2); g.Draw("LP SAME")
        leg.AddEntry(g, name, "LP")
    leg.SetHeader(f"{args.nsigma:g}#sigma timing veto:")
    leg.Draw()
    S.header(0.19, 0.88, [f"{ENERGY}", "VBF H#rightarrowinv. vs Z+jets", SEL, "After a kinematic VBF vs Z+jets BDT (j_{1}, j_{2})"], size=0.8 * TS)
    R_.gPad.RedrawAxis()
    cv.cd(); S.pad("vg_r", 0.0, 0.0, 1.0, RAT_, top=0.03, bottom=0.36)
    R_.gPad.SetLogy()
    TSR = 0.05 * 600 / (800 * RAT_)
    hn = S.keep(R_.TGraph(len(EFFS)))
    for i, (x, n_) in enumerate(zip(EFFS, nB)): hn.SetPoint(i, x, n_)
    frn = S.keep(R_.gPad.DrawFrame(0.225, 10.0, 1.025, 3e4))
    S.axes(frn, "Signal efficiency of the kinematic BDT cut (1 = no BDT)", "Z+jets MC", size=TSR, ytoff=1.35 * TS / TSR, xtoff=1.2)
    S.style(hn, "#4D4D4D", 20, 1, 2, 1.0); hn.Draw("LP SAME")
    S.save(cv, f"{args.out}_veto_gain")

if args.summary:
    # the kinematic BDT's held-out score distributions, unit area (weighted), with the cuts of the analysis
    nb = 40
    cvs = S.keep(R_.TCanvas("bdt_score", "", 800, 600))
    cvs.SetLogy()
    hs_ = {}
    for key, lab, E, sc in samples:
        h = S.keep(R_.TH1D(f"score_{key}", "", nb, 0.0, 1.0))
        for v_, w_ in zip(sc, E["weight"].astype(float)): h.Fill(min(v_, 1 - 1e-9), w_)
        h.Scale(1.0 / h.Integral())
        hs_[key] = h
    ymax_ = 30.0 * max(h.GetMaximum() for h in hs_.values())
    ymin_ = 0.3 * min(min(h.GetBinContent(b) for b in range(1, nb + 1) if h.GetBinContent(b) > 0) for h in hs_.values())
    frs = S.keep(cvs.DrawFrame(0.0, ymin_, 1.0, ymax_))
    S.axes(frs, "Kinematic BDT score (held out)", "Fraction of events")
    leg = S.legend(0.20, 0.17, 0.56, 0.31, size=0.04, margin=0.2)   # lower left: empty on both curves
    for key, colour, ls in (("sig", S.COL["hgtd"], 1), ("bkg", S.COL["trkptz"], 1)):
        h = hs_[key]; h.SetLineColor(S.colour(colour)); h.SetLineWidth(3); h.SetLineStyle(ls); h.Draw("HIST SAME")
        leg.AddEntry(h, dict((k_, l_) for k_, l_, _, _ in samples)[key], "L")
    leg.Draw()
    t_ = R_.TLatex(); t_.SetTextSize(0.032); t_.SetTextAlign(11); t_.SetTextAngle(90); t_.SetTextColor(R_.kGray + 2)
    for eff, thr in cuts[1:]:   # the cuts of the analysis: labels vertical, left of each line, so close cuts do not collide
        ln = S.keep(R_.TLine(thr, ymin_, thr, ymax_ / 300.0)); ln.SetLineStyle(2); ln.SetLineColor(R_.kGray + 2); ln.SetLineWidth(2); ln.Draw()
        S.keep(t_.DrawLatex(thr - 0.006, ymin_ * 1.5, f"#varepsilon_{{S}} = {100 * eff:.0f}%"))
    S.header(0.19, 0.87, [ENERGY, "VBF H#rightarrowinv. vs Z+jets, kinematic BDT (j_{1}, j_{2})", SEL], size=0.04)
    cvs.RedrawAxis()
    S.save(cvs, f"{args.out}_bdt_score")
    print("BDT cut thresholds: " + ", ".join(f"{100 * e:.0f}% -> score > {t:.3f}" for e, t in cuts[1:]))

