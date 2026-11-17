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
ap.add_argument("--mjj2-leg", type=lambda t: tuple(float(x) for x in t.split(",")), default=(0.61, 0.665, 0.94, 0.855),
                help="x1,y1,x2,y2 (NDC) of the compact m_jj figures' legend (default: the top-right corner, below the strip "
                     "where off-scale points are marked)")
ap.add_argument("--mjj2-ymax", type=float, default=150.0, help="the top of the compact S/sqrt(B) figures' y axis [%%]")
ap.add_argument("--mjj2-inframe", action="store_true", help="the compact S/sqrt(B) figures (veto_mjj2_{nobdt,bdt<x>}) with the "
                "ATLAS label and legend INSIDE the frame, the y axis raised until every point sits below them (the Global "
                "Performance talk, 2026-11-12); default: both in a band above the frame, the top fixed at --mjj2-ymax")
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--fixes", action="store_true", help="with --veto: also the vetoes on the truth-fixed WAVeS t0s (tdr_fwdpair "
                "trees from commit dcf6a30: perfect HS track times, perfect selection, both; todo vt-23j), and the compact "
                "m_jj pair shows WAVeS, the three fixes, truth and the perfect veto instead of jet vs jet")
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
# the JVT / fJVT working point applied before pairing (tdr_fwdpair --jvt), from the trees' meta: in the selection line
_jm = uproot.open(args.sig)["meta"]
JVT_WP = str(_jm.arrays(library="np")["jvt_wp"][0]) if "jvt_wp" in _jm.keys() else "none"
if JVT_WP != "none":
    SEL = f"j_{{1}}, j_{{2}}: |#eta| < 3.8, #geq 1 in 2.4#minus3.8, {FC.jvt_label(JVT_WP)}"

REACH = [("HS + HS", "#7BCB6B"), ("PU leg at |#eta| < 2.4", "#B3231C"), ("PU leg at 2.4#minus3.8 (HGTD)", "#3F90DA"),
         ("A leg has neither label", "#C79BF5")]
for key, lab, E, sc in samples:
    assert len(sc) == len(E["weight"]) and np.all(np.isfinite(sc)), f"{args.bdt}: {key} scores do not match the tree"
E0, s0 = samples[0][2], samples[0][3]
o = np.argsort(-s0); cs = np.cumsum(E0["weight"].astype(float)[o]) / E0["weight"].astype(float).sum()   # float64: a float32 cumsum moves the cut
cuts = [(None, None)] + [(e, s0[o][np.searchsorted(cs, e)]) for e in (float(x) for x in args.effs.split(","))]

VETOES = [("jj", "jet vs jet"), ("hgtd", "jet vs HGTD t_{0}"), ("trkptz", "jet vs TRKPTZ t_{0}"),
          ("waves", "jet vs WAVeS t_{0}"), ("truth", "jet vs truth t_{0}"), ("perfect", "perfect HGTD")]
FIX_ROWS = [("wfix_time", "Perfect timing", "kP10Cyan", 22, 1), ("wfix_sel", "Perfect selection", "kP10Green", 33, 1),
            ("wfix_both", "Perfect timing + selection", "kP10Red", 34, 1)]   # the WAVeS t0 with those failures fixed by truth
if args.fixes:
    if not all(f"t0_{k}" in E for k, *_ in FIX_ROWS for _, _, E, _ in samples):
        raise SystemExit("--fixes needs tdr_fwdpair trees with the truth-fixed t0s (commit dcf6a30)")
    VETOES[4:4] = [(k, f"jet vs WAVeS t_{{0}}, {lab.lower()}") for k, lab, *_ in FIX_ROWS]


def pu_in_hgtd(E):
    """a VBS jet is truth-labelled pileup and inside HGTD (2.4 < |eta|; the selection stops at 3.8)"""
    fwd = lambda k: np.abs(E[f"j{k}_eta"]) > 2.4
    return ((E["j1_pu"] == 1) & fwd(1)) | ((E["j2_pu"] == 1) & fwd(2))
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
    # the ceiling (the user's request, 2026-11-05): a veto that removes exactly the events with a pileup VBS jet
    # HGTD can see, by the truth labels, and no other -- every pileup jet in HGTD rejected, no signal cost
    for key, lab, E, sc in samples: K[key]["perfect"] = ~pu_in_hgtd(E)
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
    # each veto's efficiencies relative to the BDT cut, against that cut: drawn as the S/sqrt(B) factor
    # (eps_S / sqrt(eps_B), <out>_veto_gain) and the S/B factor (eps_S / eps_B, <out>_veto_sb)
    EFFS = np.round(np.arange(1.0, 0.249, -0.05), 2)
    Es, ss = samples[0][2], samples[0][3]; Eb, sb = samples[1][2], samples[1][3]
    ws_, wb_ = Es["weight"].astype(float), Eb["weight"].astype(float)
    o_ = np.argsort(-ss); cs_ = np.cumsum(ws_[o_]) / ws_.sum()
    ROWS = [("waves", "WAVeS t_{0}", S.COL["waves"], 20, 1), ("trkptz", "TRKPTZ t_{0}", S.COL["trkptz"], 21, 1),
            ("hgtd", "HGTD t_{0}", S.COL["hgtd"], 23, 1), ("truth", "Truth t_{0}", S.COL["truth"], 25, 3),
            ("jj", "Jet vs jet (no t_{0})", "kP10Violet", 33, 10),
            ("perfect", "Perfect HGTD veto (truth)", "#000000", 34, 2)]
    if args.fixes: ROWS[1:1] = FIX_ROWS
    thr_of = lambda e: -np.inf if e >= 1.0 else ss[o_][np.searchsorted(cs_, e)]

    def veto_effs(v, bS, bB):
        """eps_S, eps_B of veto v on the events bS / bB, the weighted-binomial variance of eps_B, and that of eps_S"""
        wS = ws_[bS]; kS = K["sig"][v][bS]
        es_ = wS[kS].sum() / wS.sum()
        wB = wb_[bB]; kB = K["bkg"][v][bB]
        eb_ = wB[kB].sum() / wB.sum()
        var = ((wB[kB] ** 2).sum() * (1 - eb_) ** 2 + (wB[~kB] ** 2).sum() * eb_ ** 2) / wB.sum() ** 2
        var_s = ((wS[kS] ** 2).sum() * (1 - es_) ** 2 + (wS[~kS] ** 2).sum() * es_ ** 2) / wS.sum() ** 2
        return es_, eb_, var, var_s

    EV, nB = {}, []
    for v, *_ in ROWS:
        EV[v] = [veto_effs(v, ss >= thr_of(e), sb >= thr_of(e)) for e in EFFS]
    nB = [int(np.sum(sb >= thr_of(e))) for e in EFFS]
    MAIN_, RAT_ = 0.70, 0.30
    TS = 0.05 * 600 / (800 * MAIN_)
    TSR = 0.05 * 600 / (800 * RAT_)

    def draw_gain(stem, ytitle, fac, lo, hi, what):
        print(f"\n{what} factor of each veto against the BDT cut (signal efficiency): " + " ".join(f"{e:.2f}" for e in EFFS))
        cv = S.keep(R_.TCanvas(stem, "", 800, 800))
        cv.cd(); S.pad(stem + "_m", 0.0, RAT_, 1.0, 1.0, top=0.05 / MAIN_ * 0.6, bottom=0.02)
        fr = S.keep(R_.gPad.DrawFrame(0.225, lo, 1.025, hi))
        S.axes(fr, "", ytitle, size=TS, ytoff=1.35)
        fr.GetXaxis().SetLabelSize(0)
        one = S.keep(R_.TLine(0.225, 1.0, 1.025, 1.0)); one.SetLineStyle(2); one.SetLineColor(R_.kGray + 2); one.Draw()
        leg = S.legend(0.19, 0.45, 0.75, 0.645, size=0.68 * TS, ncol=2, margin=0.15)   # under the header, left of the ceiling's rise
        for v, name, colour, mk, ls in ROWS:
            y, ey = np.array([fac(*t)[0] for t in EV[v]]), np.array([fac(*t)[1] for t in EV[v]])
            g = S.style(S.graph(EFFS, y, ey), colour, mk, ls, 2, 1.2); g.Draw("LP SAME")
            leg.AddEntry(g, name, "LP")
            print(f"  {plain(name):24s} " + " ".join(f"{a:.3f}" for a in y))
        leg.SetHeader(f"{args.nsigma:g}#sigma timing veto:")
        leg.Draw()
        S.header(0.19, 0.88, [f"{ENERGY}", "VBF H#rightarrowinv. vs Z+jets", SEL, "After a kinematic VBF vs Z+jets BDT (j_{1}, j_{2})"], size=0.8 * TS)
        R_.gPad.RedrawAxis()
        cv.cd(); S.pad(stem + "_r", 0.0, 0.0, 1.0, RAT_, top=0.03, bottom=0.36)
        R_.gPad.SetLogy()
        hn = S.keep(R_.TGraph(len(EFFS)))
        for i, (x, n_) in enumerate(zip(EFFS, nB)): hn.SetPoint(i, x, n_)
        frn = S.keep(R_.gPad.DrawFrame(0.225, 10.0, 1.025, 3e4))
        S.axes(frn, "Signal efficiency of the kinematic BDT cut (1 = no BDT)", "Z+jets MC", size=TSR, ytoff=1.35 * TS / TSR, xtoff=1.2)
        S.style(hn, "#4D4D4D", 20, 1, 2, 1.0); hn.Draw("LP SAME")
        S.save(cv, f"{args.out}_{stem}")

    draw_gain("veto_gain", "S/#sqrt{B} with the veto / without",
              lambda es_, eb_, var, var_s=0.0: (es_ / np.sqrt(eb_), es_ / np.sqrt(eb_) * 0.5 * np.sqrt(var) / eb_), 0.7, 2.6, "S/sqrt(B)")
    draw_gain("veto_sb", "S/B with the veto / without",
              lambda es_, eb_, var, var_s=0.0: (es_ / eb_, es_ / eb_ * np.sqrt(var) / eb_), 0.0, 5.5, "S/B")   # from 0: the legend sits below every curve

    # the money-plot style: the S/sqrt(B) (<out>_veto_mjj) and S/B (<out>_veto_mjj_sb, the user's request 2026-11-05:
    # "a percent just like" the S/sqrt(B) one) change of each veto per m_jj column, at one BDT cut
    EMJ = 0.8                                            # the loosest cut, as the "what the veto removes" slides
    bS, bB = ss >= thr_of(EMJ), sb >= thr_of(EMJ)
    xs_, xb_ = np.clip(Es["mjj"], EDGES[0], EDGES[-1] - 1e-6), np.clip(Eb["mjj"], EDGES[0], EDGES[-1] - 1e-6)
    ncol = len(EDGES) - 1

    def mc_panel(cvm, stem, bS, bB, xs_, xb_, EDGES, ncol):
        """the S + B MC per m_jj column under the m_jj figures, by pair category"""
        cvm.cd(); S.pad(stem + "_r", 0.0, 0.0, 1.0, RAT_, top=0.03, bottom=0.36)
        stk = S.keep(R_.THStack(stem + "_stack", ""))
        hmax = 0.0
        for ctxt, cat_, colour in (CATS[1], CATS[0]) if ONE else ():   # (name, cat, colour): FF stacked under CF
            h = S.keep(R_.TH1D(f"{stem}_{cat_}", "", ncol, EDGES))
            for c in range(ncol):
                nS_ = np.sum(bS & (Es["cat"] == cat_) & (xs_ >= EDGES[c]) & (xs_ < EDGES[c + 1]))
                nB_ = np.sum(bB & (Eb["cat"] == cat_) & (xb_ >= EDGES[c]) & (xb_ < EDGES[c + 1]))
                h.SetBinContent(c + 1, (nS_ + nB_) / 1000.0)
            h.SetFillColor(S.colour(colour)); h.SetLineColor(R_.kBlack); stk.Add(h)
        tot = [(np.sum(bS & (xs_ >= EDGES[c]) & (xs_ < EDGES[c + 1])), np.sum(bB & (xb_ >= EDGES[c]) & (xb_ < EDGES[c + 1]))) for c in range(ncol)]
        hmax = max((a + b) / 1000.0 for a, b in tot)
        frr = S.keep(R_.gPad.DrawFrame(EDGES[0], 0.0, EDGES[-1], 1.9 * hmax))
        S.axes(frr, "m_{jj} [GeV]   (last bin includes overflow)", "S + B MC [10^{3}]", size=TSR, ytoff=1.35 * TS / TSR, xtoff=1.2, ndiv_y=503)
        stk.Draw("HIST SAME")
        legr = S.legend(0.19, 0.74, 0.94, 0.90, size=0.7 * TSR, ncol=3, margin=0.2)
        for h in reversed(list(stk.GetHists())):
            if h.Integral() <= 0: continue                       # e.g. no central-forward pairs in a forward-forward figure
            legr.AddEntry(h, "Central-forward" if h.GetName() == f"{stem}_1" else "Forward-forward", "F")
        legr.AddEntry(S.keep(R_.TH1D(stem + "_dummy", "", 1, 0, 1)), "numbers: S / B MC", "")
        legr.Draw()
        tl = R_.TLatex(); tl.SetTextSize(0.6 * TSR); tl.SetTextAlign(21)
        for c, (a, b) in enumerate(tot):
            S.keep(tl.DrawLatex(0.5 * (EDGES[c] + EDGES[c + 1]), (a + b) / 1000.0 + 0.08 * hmax, f"{a / 1000:.1f}k / {b}"))
        R_.gPad.RedrawAxis()

    def mjj_extent(fac, rows, edges, emj, cat=None):
        """the lowest and highest point (change in %, +- its error) draw_mjj would draw"""
        E_ = np.asarray(edges, float)
        bS_, bB_ = (np.ones(len(ss), bool), np.ones(len(sb), bool)) if emj is None else (ss >= thr_of(emj), sb >= thr_of(emj))
        if cat is not None: bS_, bB_ = bS_ & (Es["cat"] == cat), bB_ & (Eb["cat"] == cat)
        xs__, xb__ = np.clip(Es["mjj"], E_[0], E_[-1] - 1e-6), np.clip(Eb["mjj"], E_[0], E_[-1] - 1e-6)
        pts = []
        for v, *_ in rows:
            for c in range(len(E_) - 1):
                f_, ef_ = fac(*veto_effs(v, bS_ & (xs__ >= E_[c]) & (xs__ < E_[c + 1]), bB_ & (xb__ >= E_[c]) & (xb__ < E_[c + 1])))
                pts += [100 * (f_ - 1) - 100 * ef_, 100 * (f_ - 1) + 100 * ef_]
        return min(pts), max(pts)

    def draw_mjj(stem, what, ytitle, fac, cap, legpos=(0.19, 0.39, 0.50, 0.645), legsize=0.6, rows=ROWS, edges=EDGES, emj=EMJ,
                 clear_legend=False, ylim=None, offscale_arrows=False, cat=None, header_above=False, as_kept=False):
        """fac(eps_S, eps_B, var(eps_B)) -> (factor, its error); drawn as a change in %. rows: the vetoes drawn; edges: the
        m_jj columns (the last includes overflow); emj: the BDT cut's signal efficiency, None for no BDT"""
        ROWS, EDGES, EMJ = rows, np.asarray(edges, float), emj
        bS, bB = (np.ones(len(ss), bool), np.ones(len(sb), bool)) if emj is None else (ss >= thr_of(emj), sb >= thr_of(emj))
        if cat is not None: bS, bB = bS & (Es["cat"] == cat), bB & (Eb["cat"] == cat)   # one pair type (1 CF, 2 FF)
        xs_, xb_ = np.clip(Es["mjj"], EDGES[0], EDGES[-1] - 1e-6), np.clip(Eb["mjj"], EDGES[0], EDGES[-1] - 1e-6)
        ncol = len(EDGES) - 1
        cvm = S.keep(R_.TCanvas(stem, "", 800, 800))
        # header_above: the header (top left) and the legend (top right) in a band above the frame, so the whole frame
        # is free for the points and the off-scale arrows (the user's 150% cap put them under an in-frame header)
        TOPM = 0.30 if header_above else 0.05 / MAIN_ * 0.6
        cvm.cd(); S.pad(stem + "_m", 0.0, RAT_, 1.0, 1.0, top=TOPM, bottom=0.02)
        vals = {}
        for v, *_ in ROWS:
            out = []
            for c in range(ncol):
                inS = bS & (xs_ >= EDGES[c]) & (xs_ < EDGES[c + 1]); inB = bB & (xb_ >= EDGES[c]) & (xb_ < EDGES[c + 1])
                f_, ef_ = fac(*veto_effs(v, inS, inB))
                out.append((100 * (f_ - (0 if as_kept else 1)), 100 * ef_))
            vals[v] = (out, 100 * (fac(*veto_effs(v, bS, bB))[0] - (0 if as_kept else 1)))
        top = max(max(a + e for a, e in vals[v][0]) for v in vals)
        bot = min(min(a - e for a, e in vals[v][0]) for v in vals)
        # room above the first two columns for the legend; the last column (few Z+jets MC events) may be clipped
        lo_, hi_ = min(-20.0, 10 * np.floor(bot / 10) - 10), max(100.0, min(cap, 10 * np.ceil((top + 30) / 10)))
        if ylim is not None:                             # a range shared with a companion figure
            lo_, hi_ = ylim
        elif clear_legend:   # raise the frame until every point left of the legend's right edge sits below it
            xpt = lambda k, c: EDGES[c] + (EDGES[c + 1] - EDGES[c]) * (k + 1) / (len(ROWS) + 1)
            x2 = EDGES[0] + (legpos[2] - 0.16) / (1 - 0.16 - 0.05) * (EDGES[-1] - EDGES[0])   # the pad's margins
            tl_ = max(vals[v][0][c][0] + vals[v][0][c][1] for k, (v, *_) in enumerate(ROWS) for c in range(ncol) if xpt(k, c) < x2)
            fr_ = (legpos[1] - 0.02) / (1 - 0.02 - 0.05 / MAIN_ * 0.6)                          # legend bottom, in frame units
            hi_ = max(hi_, 10 * np.ceil((lo_ + (tl_ - lo_) / (0.92 * fr_)) / 10))
        frm = S.keep(R_.gPad.DrawFrame(EDGES[0], lo_, EDGES[-1], hi_))
        S.axes(frm, "", ytitle, size=TS, ytoff=1.35)
        frm.GetXaxis().SetLabelSize(0)
        yref = 100.0 if as_kept else 0.0                 # everything kept / no change
        zero = S.keep(R_.TLine(EDGES[0], yref, EDGES[-1], yref)); zero.SetLineStyle(2); zero.SetLineColor(R_.kGray + 2); zero.Draw()
        for c in range(1, ncol):
            ln = S.keep(R_.TLine(EDGES[c], lo_, EDGES[c], hi_)); ln.SetLineStyle(3); ln.SetLineColor(R_.kGray + 1); ln.Draw()
        legm = S.legend(*legpos, size=legsize * TS, margin=0.18)   # default: above column 1, below the header, left of the ceiling's column-2 point
        print(f"\n{what} {'[%]' if as_kept else 'change [%]'} per m_jj column " + ("with no BDT" if EMJ is None else f"at the BDT cut eps_S = {EMJ:.0%}") + f" ({stem}): " +
              ", ".join(f"[{EDGES[c]:.0f},{EDGES[c + 1]:.0f})" for c in range(ncol)) + "; all")
        for k, (v, name, colour, mk, ls) in enumerate(ROWS):
            xs = [EDGES[c] + (EDGES[c + 1] - EDGES[c]) * (k + 1) / (len(ROWS) + 1) for c in range(ncol)]
            g = S.style(S.graph(np.array(xs), np.array([a for a, _ in vals[v][0]]), np.array([e for _, e in vals[v][0]])), colour, mk, 1, 2, 1.3)
            g.Draw("P SAME")                                 # points with errors, no lines: the columns are independent
            legm.AddEntry(g, f"{name} ({vals[v][1]:.0f}%)" if as_kept else f"{name} ({vals[v][1]:+.0f}%)", "PE")
            print(f"  {plain(name):24s} " + " ".join(f"{a:+6.1f} +- {e:4.1f}" for a, e in vals[v][0]) + f"   all {vals[v][1]:+.1f}")
        legm.SetHeader(f"{args.nsigma:g}#sigma timing veto (all m_{{jj}}):")
        legm.Draw()
        if offscale_arrows:                              # a point above the frame: an arrow at its own place, with its value
            for k, (v, name, colour, mk, ls) in enumerate(ROWS):
                for c in range(ncol):
                    a_ = vals[v][0][c][0]
                    if a_ <= hi_: continue
                    x_ = EDGES[c] + (EDGES[c + 1] - EDGES[c]) * (k + 1) / (len(ROWS) + 1)
                    ar = S.keep(R_.TArrow(x_, hi_ - 0.10 * (hi_ - lo_), x_, hi_ - 0.01 * (hi_ - lo_), 0.012, "|>"))
                    ar.SetLineColor(S.colour(colour)); ar.SetFillColor(S.colour(colour)); ar.SetLineWidth(2); ar.Draw()
                    tv = R_.TLatex(); tv.SetTextSize(0.5 * TS); tv.SetTextAlign(32); tv.SetTextColor(S.colour(colour))
                    S.keep(tv.DrawLatex(x_ - 0.015 * (EDGES[-1] - EDGES[0]), hi_ - 0.055 * (hi_ - lo_), f"{a_:+.0f}%"))
        else:
            ta = R_.TLatex(); ta.SetTextSize(0.5 * TS); ta.SetTextAlign(23); ta.SetTextColor(R_.kGray + 2)
            for c in range(ncol):                        # say so when a column's points are above the frame
                n_off = sum(vals[v][0][c][0] > hi_ for v, *_ in ROWS)
                if n_off: S.keep(ta.DrawLatex(0.5 * (EDGES[c] + EDGES[c + 1]), hi_ - 0.03 * (hi_ - lo_), f"#uparrow {n_off} of {len(ROWS)} off scale"))
        S.header(0.19, 0.955 if header_above else 0.88, [f"{ENERGY}", "S: VBF H#rightarrowinv., B: Z+jets",
                              ("Forward-forward pairs" if cat == 2 else "Central-forward pairs" if cat == 1 else "#geq 1 VBS jet in HGTD")
                              + (f", {FC.jvt_short(JVT_WP)}" if JVT_WP != "none" else ""),   # short: the legend sits beside this line
                              "No BDT: the veto alone" if EMJ is None else f"BDT at #varepsilon_{{S}} = {100 * EMJ:.0f}%, then the veto"], size=0.8 * TS)   # short: clear of column 3
        R_.gPad.RedrawAxis()
        mc_panel(cvm, stem, bS, bB, xs_, xb_, EDGES, ncol)
        R_.gPad.RedrawAxis()
        S.save(cvm, f"{args.out}_{stem}")

    draw_mjj("veto_mjj", "S/sqrt(B)", "S/#sqrt{B} change vs no veto [%]",
             lambda es_, eb_, var, var_s=0.0: (es_ / np.sqrt(eb_), es_ / np.sqrt(eb_) * 0.5 * np.sqrt(var) / eb_), 260.0)
    draw_mjj("veto_mjj_sb", "S/B", "S/B change vs no veto [%]",
             lambda es_, eb_, var, var_s=0.0: (es_ / eb_, es_ / eb_ * np.sqrt(var) / eb_), 420.0,
             legpos=(0.565, 0.37, 0.93, 0.62), legsize=0.55)            # upper right: column 3's points are low or off scale
    # the compact pair (the advisor's request via the user, 2026-11-05): no BDT and the 80% cut side by side; WAVeS, the
    # truth t0, jet vs jet and the perfect veto only; m_jj < 1 TeV and > 1 TeV (the 15-event > 2 TeV column merged in)
    # One y range per quantity, shared by the no-BDT and BDT figures so they read on the same scale (the user,
    # 2026-11-05), with the legend in the top-right corner: the frame is tall enough that every point sits below
    # the header and the legend (the top ~35% of the frame).
    ROWS4 = [r for r in ROWS if r[0] in (("waves", "wfix_time", "wfix_sel", "wfix_both", "truth", "perfect") if args.fixes
                                         else ("waves", "truth", "jj", "perfect"))]
    E2 = [0.0, 1000.0, 4000.0]
    for stem_, what_, ytitle_, fac_ in (
            ("veto_mjj2", "S/sqrt(B)", "S/#sqrt{B} change vs no veto [%]",
             lambda es_, eb_, var, var_s=0.0: (es_ / np.sqrt(eb_), es_ / np.sqrt(eb_) * 0.5 * np.sqrt(var) / eb_)),
            ("veto_mjj2_sb", "S/B", "S/B change vs no veto [%]",
             lambda es_, eb_, var, var_s=0.0: (es_ / eb_, es_ / eb_ * np.sqrt(var) / eb_))):
        ext = [mjj_extent(fac_, ROWS4, E2, e_) for e_ in (None, EMJ)]
        lo_s = min(-20.0, 10 * np.floor(min(a for a, _ in ext) / 10) - 10)
        if what_ == "S/sqrt(B)" and args.mjj2_inframe:   # label and legend in the frame's top ~40%: every point below
            top_s = max(b for _, b in ext)
            hi_s = 10.0 * np.ceil((lo_s + (top_s - lo_s) / 0.58) / 10.0)
        elif what_ == "S/sqrt(B)":   # a fixed top (the user, 2026-11-05: "keep the maximum to 150%, noting when a point is
            hi_s = args.mjj2_ymax  # off scale"): off-scale points become arrows with their value, in a strip above the legend
        else:
            top_s = max(b for _, b in ext)
            step = 50.0 if top_s < 300 else 100.0
            hi_s = step * np.ceil((lo_s + (top_s - lo_s) / 0.62) / step)
        LEG_ABOVE = (0.61, 0.715, 0.94, 0.975)               # the band above the frame (header_above), right of the header
        above = what_ == "S/sqrt(B)" and not args.mjj2_inframe
        inframe = what_ == "S/sqrt(B)" and args.mjj2_inframe
        for emj_, tag_ in ((None, "nobdt"), (EMJ, f"bdt{round(100 * EMJ)}")):
            draw_mjj(f"{stem_}_{tag_}", what_, ytitle_, fac_, 1e9,
                     legpos=LEG_ABOVE if above else (0.62, 0.625, 0.94, 0.895) if inframe else args.mjj2_leg,
                     legsize=0.5 if (above or inframe) else 0.46, rows=ROWS4, edges=E2, emj=emj_, ylim=(lo_s, hi_s), offscale_arrows=True,
                     header_above=above)
            if what_ == "S/sqrt(B)":           # the efficiencies behind it (the user, 2026-11-06), one figure per sample
                for cat_, ctag in ((None, ""), (2, "_ff")) if ONE else ((None, ""),):
                    draw_mjj(f"veto_mjj2_effs{ctag}_{tag_}", "VBF kept", "VBF H#rightarrowinv. kept [%]",
                             lambda es_, eb_, var, var_s=0.0: (es_, np.sqrt(var_s)), 1e9, legpos=LEG_ABOVE, legsize=0.5,
                             rows=ROWS4, edges=E2, emj=emj_, ylim=(60.0, 102.0), cat=cat_, header_above=True, as_kept=True)
                    draw_mjj(f"veto_mjj2_effb{ctag}_{tag_}", "Z+jets kept", "Z+jets kept [%]",
                             lambda es_, eb_, var, var_s=0.0: (eb_, np.sqrt(var)), 1e9, legpos=LEG_ABOVE, legsize=0.5,
                             rows=ROWS4, edges=E2, emj=emj_, ylim=(0.0, 105.0), cat=cat_, header_above=True, as_kept=True)
            if what_ == "S/sqrt(B)" and ONE:   # the forward-forward pairs alone (the user, 2026-11-06), where jet vs jet can act
                draw_mjj(f"{stem_}_ff_{tag_}", what_, ytitle_, fac_, 1e9, legpos=LEG_ABOVE, legsize=0.5,
                         rows=ROWS4, edges=E2, emj=emj_, ylim=(lo_s, hi_s), offscale_arrows=True, cat=2, header_above=True)

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

