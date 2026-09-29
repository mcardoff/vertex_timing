#!/usr/bin/env python3
"""HGTD timing tests inside the VBF H->inv signal region, weighted: the full
Run 2 analysis' SR (ANA-EXOT-2020-11-INT1, Sec. 4.5.2, and its Z->ll control
region, Sec. 4.5.6) emulated on util/vbs_time_veto --sr outputs.

The selection, applied here per timing variant (a variant that removes jets
re-runs all of it on the jets that are left):
  jets          the stored jets: EMTopo, >= 20 GeV, |eta| < 4.5, not lepton-
                overlap removed. JVT is the in-code R_pT proxy at the loose
                point (src/jet_tagging.h; the Run-2 Medium equivalent) inside
                its window, |eta| < 2.5 and pT < 60 GeV. Counted jets: > 25 GeV
                and passing JVT.
  N_jet         2, 3 or 4 counted jets
  j1, j2        the two leading counted jets: pT > 80 / 50 GeV. Each inside the
                fJVT window (2.5 <= |eta| < 4.5, pT < 120 GeV) needs fJVT < 0.5
                when E_T^miss > 200 GeV, else < 0.2 (the note's Very Tight)
  j3, j4        max centrality < 0.6 and max m_rel < 0.05
  E_T^jet       |vector sum of EVERY stored jet| > 140 GeV: JVT-rejected and
                timing-removed jets count, as in the note's E_T^{jet,no-jvt}
  pair          |dphi_jj| < 2, eta1 * eta2 < 0, |deta_jj| > 3.8, m_jj > 800 GeV
  E_T^miss      the union of the 16 SR bins: 2 jets and E_T^miss > 200 GeV;
                2 jets, 160-200 GeV and m_jj > 1.5 TeV; 3-4 jets, > 200 GeV and
                m_jj > 1.5 TeV
  Z->ll         (background) |m_ll - 91.2| < 25 GeV
E_T^miss is truth |sum pT(nu)| for VBF H->inv and pT(ll) for Z(->ll)+jets, the
Z->ll CR's "E_T^miss without leptons": the ntuples have no reco E_T^miss. Not
emulated: the soft-term cut, E_T^miss (with leptons) < 70 GeV, the lepton,
photon and b-jet vetoes, and PFlow jets (this production's PFlow is
miscalibrated forward: HS-jet response 1.14 at 2.4-4.0, 1.31 beyond).

Timing. A jet is timed when its core time has >= --min-trk tracks; its sigma is
f(n) x sigma, f measured on the signal file's HS legs (vbs_time_veto_plot.py's
default), and it fails against a t0 when |t - t0| / sqrt(sigma^2 + (infl
sigma_t0)^2) >= --nsigma:
  none          the SR above
  jj            event veto: j1 and j2 both timed and mutually incompatible
  t0_<src>      event veto: j1 or j2 timed and incompatible with the t0
  rm_<src>      every timed jet incompatible with the t0 is removed, then the
                whole selection runs on the rest (j1, j2, N_jet, j3/j4 re-formed)
Every number is relative to `none`. Removal can ADD events (a removed pileup
jet no longer fails N_jet or the extra-jet cuts), so they are yield ratios, not
efficiencies; their errors come from the events the two selections do not
share. Events are binned in the m_jj of the pair they end up with.

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_time_veto_sr.py \\
        --sig condor/vbf/vbf_sr_vbs_time_veto.root \\
        --bkg condor/zjets/zjets_sr_vbs_time_veto.root --out figs/time_veto_sr/sr

Writes <out>_{sig,bkg}_ratio_mjj, <out>_s_over_sqrtb (% change per m_jj column),
<out>_s_over_sqrtb_yields (at --lumi), <out>_tradeoff (threshold scan), each
.pdf and .png, and <out>_summary.md (cut flow, method table, SR bins, scan).
"""
import argparse, os
import numpy as np
import uproot
import awkward as ak
import vbs_time_veto_common as C
from vbf_hinv_sr import METHODS, load, take, sr_select, jet_sigma, outcomes, ratio

MJJ_EDGES = [800, 1000, 1500, 2000, 3500, 5000]     # the SR's m_jj bins; the last one includes the overflow
SCAN = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0])

STYLE = dict(C.STYLE)
STYLE.update({f"rm_{s}": C.STYLE[f"rp_{s}"] for s in C.T0_SRC})

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True, help="signal (VBF) vbs_time_veto --sr output, from the RAW ntuples")
ap.add_argument("--bkg", help="background (Z+jets) vbs_time_veto --sr output")
ap.add_argument("--sig-norm", default="vbf", choices=list(C.NORM) + ["none"])
ap.add_argument("--bkg-norm", default="zjets", choices=list(C.NORM) + ["none"])
ap.add_argument("--sig-label", default="#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.")
ap.add_argument("--bkg-label", default="#sqrt{s} = 14 TeV, HL-LHC, Z(#rightarrowll)+jets")
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--min-trk", type=int, default=1, help="tracks a jet-time estimate needs for the jet to count as timed")
ap.add_argument("--methods", default="t0_trkptz,t0_waves,t0_hgtd,t0_truth,rm_trkptz,rm_waves,rm_hgtd,rm_truth",
                help="drawn methods; the summary has all of " + ",".join(METHODS))
ap.add_argument("--lumi", type=float, default=3000.0, help="fb^-1")
ap.add_argument("--bf-hinv", type=float, default=1.0, help="B(H->inv) for the signal yields; cancels in every ratio")
ap.add_argument("--out", required=True, help="output stem")
args = ap.parse_args()
DRAWN = args.methods.split(",")
for m in DRAWN:
    if m not in METHODS: raise SystemExit(f"unknown method {m}; choose from {METHODS}")
EDGES = np.array(MJJ_EDGES, float)


def k_factor(s):
    if s["norm"] == "none": return None
    n = C.NORM[s["norm"]]
    # the local sample is not the grid production: its own sum of weights, same cross-section
    sumw = s["info"]["sumw_read"] if s["info"]["sample"] == "local" else n["sumw"]
    return args.lumi * 1e3 * n["xs_pb"] / sumw * (args.bf_hinv if s["norm"] == "vbf" else 1.0)


samples = {"sig": load(args.sig, "sig", args.sig_norm, "met_truth", False)}
if args.bkg: samples["bkg"] = load(args.bkg, "bkg", args.bkg_norm, "z_pt", True)
calib = C.measure(samples["sig"]["e"], True, "the sig file's HS legs")   # one procedure for S and B
for s in samples.values():
    s["calib"] = calib
    s["sig"] = jet_sigma(s)
    s["none"] = sr_select(s, ak.ones_like(s["J"].jet_pt, dtype=bool))
    s["out"] = outcomes(s, args.nsigma, args.min_trk)
    print(f"[{s['role']}] {s['path']}: {s['info']['njobs']} job(s), {s['info']['n_read']:,} read, "
          f"{len(s['w']):,} SR candidates stored, {int(s['none']['ok'].sum()):,} in the SR")

# ── Summary ──────────────────────────────────────────────────────────────────
both = "bkg" in samples
S_, B_ = samples["sig"], samples.get("bkg")
md = [f"# VBS timing tests in the VBF H->inv SR: {os.path.basename(args.out)}", "",
      "ANA-EXOT-2020-11-INT1 Sec. 4.5.2 / 4.5.6 emulated (see the script's docstring for what is and is not). "
      f"Timing test at {args.nsigma:g} sigma; a jet is timed with >= {args.min_trk} track(s); jet-time sigma "
      f"{calib.md()}; yields at {args.lumi:g} fb^-1.", ""]
for s in samples.values():
    i = s["info"]
    md.append(f"- **{s['role']}**: `{s['path']}` -- {i['njobs']} job(s), {i['n_read']:,} events read (sum w "
              f"{i['sumw_read']:.6g}), {len(s['w']):,} SR candidates ({s['rcol']} >= {i['sr_recoil_min']:g}); "
              f"normalisation: {C.NORM[s['norm']]['note'] if s['norm'] != 'none' else 'none'}")
md += ["", "**Cut flow of the SR without timing** (cumulative; MC events and the weighted fraction of the stored "
       "candidates)", "", "| cut | " + " | ".join(f"{r} MC | {r} frac." for r in samples) + " |",
       "|---|" + "---|---|" * len(samples)]
masks = {r: np.ones(len(s["w"]), bool) for r, s in samples.items()}
names = [n for n, _ in (B_ or S_)["none"]["steps"]]
for k, name in enumerate(names):
    cells = []
    for r, s in samples.items():
        st = s["none"]["steps"]
        if k < len(st): masks[r] &= st[k][1]
        cells.append(f"{int(masks[r].sum()):,} | {100 * s['w'][masks[r]].sum() / s['w'].sum():.2f}%")
    md.append(f"| {name} | " + " | ".join(cells) + " |")
md.append("")

hdr = "| method | S / S(none) | S MC |" + (" B / B(none) | B MC | S/sqrt(B) change |" if both else "")
md += ["**Timing relative to the SR without it** (errors: the events the two selections do not share)", "",
       hdr, "|" + "---|" * (hdr.count("|") - 1)]
summ = {}
for m in ["none"] + METHODS:
    rs, es = ratio(S_["w"], S_["out"][m]["ok"], S_["none"]["ok"])
    line = f"| {m} | {rs:.3f} +- {es:.3f} | {int(S_['out'][m]['ok'].sum()):,} |"
    if both:
        rb, eb = ratio(B_["w"], B_["out"][m]["ok"], B_["none"]["ok"])
        z = rs / np.sqrt(rb) if rb > 0 else np.nan
        ez = z * np.hypot(es / rs, 0.5 * eb / rb) if rb > 0 else np.nan
        summ[m] = (rs, es, rb, eb, z, ez)
        line += f" {rb:.3f} +- {eb:.3f} | {int(B_['out'][m]['ok'].sum()):,} | {100 * (z - 1):+.1f} +- {100 * ez:.1f}% |"
    md.append(line)
md.append("")
if both and k_factor(S_) and k_factor(B_):
    ks, kb = k_factor(S_), k_factor(B_)
    S0, B0 = ks * S_["w"][S_["none"]["ok"]].sum(), kb * B_["w"][B_["none"]["ok"]].sum()
    md += [f"Without timing: S = {S0:,.0f}, B = {B0:,.0f} (Z->ll+jets standing in for Z->nunu), "
           f"S/sqrt(B) = {S0 / np.sqrt(B0):.1f} at {args.lumi:g} fb^-1, B(H->inv) = {100 * args.bf_hinv:g}%.", ""]
    md += ["**SR bins without timing** (S, B yields; MC events in brackets)", "",
           "| bin | S | B |", "|---|---|---|"]
    for b in range(1, 17):
        ms, mb = S_["none"]["bin"] == b, B_["none"]["bin"] == b
        md.append(f"| SR{b} | {ks * S_['w'][ms].sum():,.0f} ({int(ms.sum()):,}) | {kb * B_['w'][mb].sum():,.0f} ({int(mb.sum())}) |")
    md.append("")
md += ["**Leading pair of the SR without timing** (weighted share; paper labels)", ""]
for r, s in samples.items():
    J, ok, i1, i2 = s["J"], s["none"]["ok"], s["none"]["i1"], s["none"]["i2"]
    h1, h2 = take(J.jet_hs, i1) > 0, take(J.jet_hs, i2) > 0
    p1, p2 = take(J.jet_pu, i1) > 0, take(J.jet_pu, i2) > 0
    w = s["w"][ok]
    comp = {"HS+HS": h1 & h2, "HS+PU": (h1 & p2) | (p1 & h2), "PU+PU": p1 & p2}
    md.append(f"- {r}: " + ", ".join(f"{k} {100 * w[v[ok]].sum() / w.sum():.1f}%" for k, v in comp.items())
              + f", other {100 * (1 - sum(w[v[ok]].sum() for v in comp.values()) / w.sum()):.1f}%")
md.append("")

# ── Threshold scan ───────────────────────────────────────────────────────────
scan = {m: [] for m in METHODS}
if both:
    for ns in SCAN:
        res = {r: outcomes(s, ns, args.min_trk) for r, s in samples.items()}
        for m in METHODS:
            rs = ratio(S_["w"], res["sig"][m]["ok"], S_["none"]["ok"])[0]
            rb = ratio(B_["w"], res["bkg"][m]["ok"], B_["none"]["ok"])[0]
            scan[m].append((rs, rb))
    md += ["**Threshold scan: S / S(none) and background rejection B(none) / B**", "",
           "| method | " + " | ".join(f"{t:g} sigma" for t in SCAN) + " |", "|" + "---|" * (len(SCAN) + 1)]
    md += [f"| {m} | " + " | ".join(f"{rs:.3f} / {1 / rb:.2f}" if rb > 0 else f"{rs:.3f} / inf" for rs, rb in scan[m])
           + " |" for m in METHODS]
    md.append("")
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
open(args.out + "_summary.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
print("wrote", args.out + "_summary.md")

# ── Figures ──────────────────────────────────────────────────────────────────
ROOT = C.setup_root(args.out)
JDR, JCUT = S_["info"]["jet_time_dr"], S_["info"]["jet_time_dist_cut"]
SEL1 = ("SR: j_{1,2} > 80, 50 GeV, #eta_{1}#eta_{2} < 0, |#Delta#eta_{jj}| > 3.8, m_{jj} > 0.8 TeV, "
        "|#Delta#phi_{jj}| < 2")
SEL2 = "2-4 jets > 25 GeV (JVT), j_{3,4}: C < 0.6, m_{rel} < 0.05; fJVT on j_{1,2}; E_{T}^{jet} > 140 GeV"
TIMING = (f"timing |#Deltat| / #sigma #geq {args.nsigma:g}: filled = event veto on j_{{1,2}}, "
          "open = failing jets removed first")
CALIB = (f"#sigma_{{jet}} #times {calib.fac[0]:.2f} (1 track) #rightarrow {calib.fac[-1]:.2f} "
         f"(#geq {C.NBINS[-2]} tracks), from signal HS jets, for S and B")


def met_line(role):
    if role == "sig": return "E_{T}^{miss,truth} > 200 GeV, or 160-200 GeV with 2 jets and m_{jj} > 1.5 TeV"
    return "p_{T}^{ll} > 200 GeV, or 160-200 GeV with 2 jets and m_{jj} > 1.5 TeV; |m_{ll} - m_{Z}| < 25 GeV"


def per_column(s, m):
    """per final-m_jj column: (ratio to none, err, err, MC events of none)"""
    res = []
    ci_m, ci_0 = C.mjj_index(s["out"][m]["mjj"], EDGES), C.mjj_index(s["none"]["mjj"], EDGES)
    for b in range(len(EDGES) - 1):
        r, e = ratio(s["w"], ci_m == b, ci_0 == b)
        res.append((r, e, e, int((ci_0 == b).sum())))
    return res


for r, ylab, label in [("sig", "Signal yield / no timing", args.sig_label),
                       ("bkg", "Background yield / no timing", args.bkg_label)]:
    if r not in samples: continue
    s = samples[r]
    series = {m: per_column(s, m) for m in DRAWN}
    extra = {m: f"  ({ratio(s['w'], s['out'][m]['ok'], s['none']['ok'])[0]:.3f})" for m in DRAWN}
    counts = [n for _, _, _, n in series[DRAWN[0]]]
    C.column_figure(f"{args.out}_{r}_ratio_mjj", ylab, series, extra, label,
                    [SEL1, SEL2, met_line(r), TIMING, CALIB + ";  all bins in brackets"],
                    counts, "MC events (no timing)", EDGES, STYLE, DRAWN, xnote="final pair; last bin includes overflow")

if both:
    change = {}
    for m in DRAWN:
        pts = []
        for (rs, es, _, _), (rb, eb, _, _) in zip(per_column(S_, m), per_column(B_, m)):
            if not (np.isfinite(rs) and np.isfinite(rb)) or rb <= 0 or rs <= 0:
                pts.append((np.nan, np.nan, np.nan, 0)); continue
            z = rs / np.sqrt(rb)
            ez = z * np.hypot(es / rs, 0.5 * eb / rb)
            pts.append((100 * (z - 1), 100 * ez, 100 * ez, 0))
        change[m] = pts
    extra = {m: f"  ({100 * (summ[m][4] - 1):+.1f}%)" for m in DRAWN}
    counts = [n for _, _, _, n in per_column(B_, DRAWN[0])]
    C.column_figure(f"{args.out}_s_over_sqrtb", "S/#sqrt{B} change vs no timing [%]", change, extra,
                    "S: VBF H#rightarrowinv.,  B: Z(#rightarrowll)+jets",
                    [SEL1, SEL2, met_line("bkg").replace("p_{T}^{ll}", "E_{T}^{miss} (S: truth, B: p_{T}^{ll})"),
                     TIMING, "all m_{jj} bins in brackets; > 0 = the timing test raises S/#sqrt{B}"],
                    counts, "B MC events", EDGES, STYLE, DRAWN, ref=0.0, pct=True,
                    xnote="final pair; last bin includes overflow")
    ks, kb = k_factor(S_), k_factor(B_)
    if ks and kb:
        absol, extra = {}, {}
        for m in ["none"] + DRAWN:
            ci_s, ci_b = C.mjj_index(S_["out"][m]["mjj"], EDGES), C.mjj_index(B_["out"][m]["mjj"], EDGES)
            pts = []
            for b in range(len(EDGES) - 1):
                ws, wb = S_["w"][ci_s == b], B_["w"][ci_b == b]
                if ws.sum() <= 0 or wb.sum() <= 0:
                    pts.append((np.nan, np.nan, np.nan, 0)); continue
                Sy, By = ks * ws.sum(), kb * wb.sum()
                z = Sy / np.sqrt(By)
                ez = z * np.hypot(np.sqrt((ws ** 2).sum()) / ws.sum(), 0.5 * np.sqrt((wb ** 2).sum()) / wb.sum())
                pts.append((z, ez, ez, 0))
            absol[m] = pts
            extra[m] = f"  ({ks * S_['w'][S_['out'][m]['ok']].sum() / np.sqrt(kb * B_['w'][B_['out'][m]['ok']].sum()):.0f})"
        C.column_figure(f"{args.out}_s_over_sqrtb_yields", "S / #sqrt{B}  per m_{jj} bin", absol, extra,
                        f"S: VBF H#rightarrowinv.,  B: Z(#rightarrowll)+jets,  {args.lumi:g} fb^{{-1}}",
                        [SEL1, SEL2, TIMING,
                         f"S at B(H#rightarrowinv) = {100 * args.bf_hinv:g}%;  B = Z(#rightarrowll)+jets alone "
                         "(stand-in for Z(#rightarrow#nu#nu));  all bins in brackets"],
                        counts, "B MC events", EDGES, STYLE, ["none"] + DRAWN, ref=None,
                        xnote="final pair; last bin includes overflow")
    i3 = int(np.argmin(np.abs(SCAN - args.nsigma)))
    xy = {m: (np.array([p[0] for p in scan[m]]), np.array([1 / p[1] if p[1] > 0 else np.nan for p in scan[m]]))
          for m in DRAWN}
    C.roc_figure(f"{args.out}_tradeoff", xy, DRAWN, STYLE, i3,
                 ["VBF H#rightarrowinv. SR (ANA-EXOT-2020-11 selection), S: VBF H#rightarrowinv., B: Z(#rightarrowll)+jets",
                  f"threshold scanned {SCAN[0]:g}-{SCAN[-1]:g}#sigma; large markers at {args.nsigma:g}#sigma; "
                  "dashed: failing jets removed first"],
                 "Signal yield / no timing", "Background rejection vs no timing",
                 dashed=[m for m in DRAWN if m.startswith("rm_")])
