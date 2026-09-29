#!/usr/bin/env python3
"""Pileup-tagged background in the VBF H->inv signal region: how often a pileup
jet is one of the two tagging jets (HS+PU), where it sits, whether it carries an
HGTD time, and whether the jet-vs-t0 test catches it. For the signal and the two
backgrounds this production has, from util/vbs_time_veto --sr outputs:
  VBF H->inv        E_T^miss = truth |sum pT(nu)|
  Z(->ll)+jets      E_T^miss = pT(ll), the Z->ll CR's "E_T^miss without leptons"
  multijet (dijet)  E_T^miss = met_jets (--sr-met=jets): its truth E_T^miss is
                    ~0, so it reaches the SR only through FAKE E_T^miss. The JZ
                    slices are shown separately: they are not normalised to
                    each other (the ntuple weight is the generator weight).
The selection is vbf_hinv_sr.py's (ANA-EXOT-2020-11-INT1), in stages from the
tagging pair towards the full SR (each cumulative):
  pair100 j1 > 80, j2 > 50 GeV (the two leading JVT-passing jets > 25 GeV),
          eta1 eta2 < 0, |deta_jj| > 3.8, m_jj > 0.8 TeV, E_T^miss > 100 GeV
          (the --sr floor; below the SR's 160, for statistics), and the Z
          window for Z+jets
  pair    E_T^miss > 160 GeV
  +fJVT   fJVT on j1, j2
  +dphi   |dphi_jj| < 2
  +jets   2-4 jets, j3/j4 centrality < 0.6 and m_rel < 0.05
  SR      E_T^jet,no-jvt > 140 GeV and the 16 SR bins
Labels are the paper's (isJetPaperHS / isJetPaperPU): a pair is HS+HS, HS+PU,
PU+PU, or "other" (a leg with neither label). A pileup leg is central (|eta| <
2.4, no HGTD), beyond |eta| = 4.0, or in HGTD (2.4-4.0) untimed, timed and
compatible with t0, or timed and failing; an event is classed by its most
catchable pileup leg. Timing is vbs_time_veto_sr.py's: the event veto on j1, j2
at --nsigma, with the signal's f(n) jet-time calibration for every sample.

The break-even: with the event veto keeping eps_S of the signal, a background
whose PU-tagged fraction is f is kept at eps_B = f (1 - c) + (1 - f)(1 - h), c
the catch rate on its PU-tagged events and h the veto rate on the rest.
S/sqrt(B) rises when eps_B < eps_S^2, i.e. f > f* = (1 - h - eps_S^2) / (c - h).

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_hspu_study.py \\
        --sig condor/vbf/vbf_sr_vbs_time_veto.root --zjets condor/zjets/zjets_sr_vbs_time_veto.root \\
        --dijet condor/dijet/dijet_sr_vbs_time_veto.root --out figs/time_veto/hspu

Writes <out>_summary.md and <out>_{composition_dphi,composition_sr,puleg,breakeven}.(pdf|png).
"""
import argparse, os, re
import numpy as np
import awkward as ak
import vbs_time_veto_common as C
from vbf_hinv_sr import load, take, sr_select, jet_sigma, fails_t0, outcomes

STAGES = [("pair100", ["E_T^miss > 100 GeV", "j1 > 80, j2 > 50 GeV", "eta1 eta2 < 0, |deta_jj| > 3.8", "m_jj > 800 GeV",
                       "|m_ll - 91.2| < 25 GeV"]),
          ("pair",  ["E_T^miss > 160 GeV"]),
          ("+fJVT", ["fJVT on j1, j2"]),
          ("+dphi", ["|dphi_jj| < 2"]),
          ("+jets", ["2-4 jets > 25 GeV (JVT)", "j3, j4: C < 0.6, m_rel < 0.05"]),
          ("SR",    ["E_T^jet,no-jvt > 140 GeV", "in an SR bin (E_T^miss / N_jet / m_jj)"])]
COMPOSITION = [("HS+HS", "#3F7FBF"), ("HS+PU", "#E8A33D"), ("PU+PU", "#C0392B"), ("other", "#B8B8B8")]
PULEG = [("central, |#eta| < 2.4", "#7F7F7F"), ("|#eta| > 4.0", "#B8B8B8"), ("HGTD, untimed", "#9DC3E6"),
         ("HGTD, passes t_{0} test", "#E8A33D"), ("HGTD, fails t_{0} test", "#2E8B57")]
PULEG_MD = ["central", "beyond 4.0", "HGTD untimed", "HGTD timed, passes", "HGTD timed, fails"]
JZ_RE = re.compile(r"\.(8011\d\d)\.Py8EG_A14NNPDF23LO_jj_(JZ\d)\.")
SIG = "VBF H#rightarrowinv."

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True, help="VBF H->inv vbs_time_veto --sr output (RAW ntuples)")
ap.add_argument("--zjets", help="Z+jets vbs_time_veto --sr output")
ap.add_argument("--dijet", help="dijet vbs_time_veto --sr --sr-met=jets output (RAW ntuples)")
ap.add_argument("--slices", default="JZ2,JZ3,JZ4", help="dijet JZ slices to show")
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--min-trk", type=int, default=1)
ap.add_argument("--t0", default="waves", choices=C.T0_SRC, help="the t0 of the pileup-leg figure")
ap.add_argument("--out", required=True)
args = ap.parse_args()


def subset(s, mask):
    """the sample restricted to the events in mask"""
    return dict(s, e={k: v[mask] for k, v in s["e"].items()}, J=s["J"][mask], w=s["w"][mask], met=s["met"][mask],
                jvt_ok=s["jvt_ok"][mask], etjet=s["etjet"][mask], zwin=s["zwin"][mask])


samples = {SIG: load(args.sig, "sig", "vbf", "met_truth", False)}
if args.zjets: samples["Z(#rightarrowll)+jets"] = load(args.zjets, "zjets", "zjets", "z_pt", True)
if args.dijet:
    dj = load(args.dijet, "dijet", "none", "met_jets", False)
    tags = np.array([(lambda m: m.group(2) if m else "?")(JZ_RE.search(str(p))) for p in dj["e"]["file_path"]])
    for jz in args.slices.split(","):
        if (tags == jz).any(): samples[f"multijet {jz}"] = subset(dj, tags == jz)
    if not any(k.startswith("multijet") for k in samples):   # not the JZ production: one multijet sample
        samples["multijet"] = dj
BKG = [k for k in samples if k != SIG]
calib = C.measure(samples[SIG]["e"], True, "the VBF file's HS legs")

for name, s in samples.items():
    s["calib"] = calib
    s["sig"] = jet_sigma(s)
    s["none"] = sr_select(s, ak.ones_like(s["J"].jet_pt, dtype=bool))
    s["out"] = outcomes(s, args.nsigma, args.min_trk)

    def stages(steps):
        steps, m, out = dict(steps, **{"E_T^miss > 100 GeV": s["met"] > 100.0}), np.ones(len(s["w"]), bool), {}
        for st, names in STAGES:
            for n in names:
                if n in steps: m = m & steps[n]
            out[st] = m.copy()
        return out
    s["stage"] = stages(s["none"]["steps"])
    s["stage_rm"] = {src: stages(s["out"][f"rm_{src}"]["steps"]) for src in C.T0_SRC}
    J, i1, i2 = s["J"], s["none"]["i1"], s["none"]["i2"]
    # the event veto's decision on the no-timing pair, for every candidate (so any stage can use it)
    fail = {src: fails_t0(s, src, args.nsigma, args.min_trk) for src in C.T0_SRC}
    s["vetoed"] = {src: take(f, i1) | take(f, i2) for src, f in fail.items()}
    h1, h2 = take(J.jet_hs, i1) > 0, take(J.jet_hs, i2) > 0
    p1, p2 = take(J.jet_pu, i1) > 0, take(J.jet_pu, i2) > 0
    s["comp"] = {"HS+HS": h1 & h2, "HS+PU": (h1 & p2) | (p1 & h2), "PU+PU": p1 & p2}
    s["comp"]["other"] = ~(s["comp"]["HS+HS"] | s["comp"]["HS+PU"] | s["comp"]["PU+PU"])
    s["pu_tag"] = s["comp"]["HS+PU"] | s["comp"]["PU+PU"]
    cls = np.full(len(s["w"]), -1)
    for i, p in [(i1, p1), (i2, p2)]:
        aeta, timed, f = np.abs(take(J.jet_eta, i)), take(J.jet_n, i) >= args.min_trk, take(fail[args.t0], i)
        c = np.where(aeta < 2.4, 0, np.where(aeta >= 4.0, 1, np.where(~timed, 2, np.where(f, 4, 3))))
        cls = np.where(p, np.maximum(cls, c), cls)
    s["puleg"] = cls
    print(f"[{name}] {len(s['w']):,} candidates; " + ", ".join(f"{st} {int(m.sum()):,}" for st, m in s["stage"].items()))


def neff(s, m):
    w = s["w"][m]
    return w.sum() ** 2 / (w ** 2).sum() if (w ** 2).sum() > 0 else 0.0


def frac(s, m, sub):
    d = s["w"][m].sum()
    return s["w"][m & sub].sum() / d if d > 0 else np.nan


def kept_veto(s, m, src, sub=None):
    """weighted fraction of the events in m (and sub) the event veto keeps"""
    sel = m if sub is None else m & sub
    return frac(s, sel, ~s["vetoed"][src])


# ── Summary ──────────────────────────────────────────────────────────────────
md = [f"# Pileup-tagged background in the VBF H->inv SR: {os.path.basename(args.out)}", "",
      f"Stages cumulative from the tagging pair to the full SR (see the docstring). Event veto on j1, j2 at "
      f"{args.nsigma:g} sigma; jet-time sigma {calib.md()}. Multijet: the JZ slices are not normalised to "
      "each other.", ""]
for name, s in samples.items():
    i = s["info"]
    md.append(f"- **{name}**: `{s['path']}` ({i['njobs']} job(s), {i['n_read']:,} read), E_T^miss = {s['rcol']}, "
              f"{len(s['w']):,} SR candidates")
md += ["", "**MC events per stage** (effective count (sum w)^2 / sum w^2 in brackets: the dijet slices are "
       "generated with a pT bias, their weights span orders of magnitude)", "",
       "| sample | " + " | ".join(st for st, _ in STAGES) + " |", "|---|" + "---|" * len(STAGES)]
md += [f"| {n} | " + " | ".join(f"{int(s['stage'][st].sum()):,} ({neff(s, s['stage'][st]):,.0f})" for st, _ in STAGES) + " |"
       for n, s in samples.items()]
md += ["", "**Tagging-pair composition** (weighted %)", "",
       "| sample | stage | " + " | ".join(c for c, _ in COMPOSITION) + " |", "|---|---|" + "---|" * len(COMPOSITION)]
for n, s in samples.items():
    for st, _ in STAGES:
        m = s["stage"][st]
        if m.any():
            md.append(f"| {n} | {st} | " + " | ".join(f"{100 * frac(s, m, s['comp'][c]):.1f}" for c, _ in COMPOSITION) + " |")
md += ["", f"**PU-tagged events by their most catchable PU leg, t0 = {args.t0}** (weighted %)", "",
       "| sample | stage | PU-tagged MC | " + " | ".join(PULEG_MD) + " |", "|---|---|---|" + "---|" * len(PULEG_MD)]
for n, s in samples.items():
    for st, _ in STAGES:
        m = s["stage"][st] & s["pu_tag"]
        if m.any():
            md.append(f"| {n} | {st} | {int(m.sum()):,} | "
                      + " | ".join(f"{100 * frac(s, m, s['puleg'] == k):.1f}" for k in range(len(PULEG))) + " |")
md += ["", "**Event veto on j1, j2: kept fraction** of all events / of the PU-tagged ones (1 - c) / of the rest "
       "(1 - h), per t0", "", "| sample | stage | " + " | ".join(C.T0_SRC) + " |", "|---|---|" + "---|" * len(C.T0_SRC)]
for n, s in samples.items():
    for st, _ in STAGES:
        m = s["stage"][st]
        if m.any():
            md.append(f"| {n} | {st} | " + " | ".join(
                f"{kept_veto(s, m, src):.3f} / {kept_veto(s, m, src, s['pu_tag']):.3f} / {kept_veto(s, m, src, ~s['pu_tag']):.3f}"
                for src in C.T0_SRC) + " |")
md += ["", "**Jet removal first: yield relative to no timing at each stage** (the stage re-run on the kept jets)", "",
       "| sample | stage | " + " | ".join(C.T0_SRC) + " |", "|---|---|" + "---|" * len(C.T0_SRC)]
for n, s in samples.items():
    for st, _ in STAGES:
        m = s["stage"][st]
        if m.any():
            md.append(f"| {n} | {st} | " + " | ".join(f"{s['w'][s['stage_rm'][src][st]].sum() / s['w'][m].sum():.3f}"
                                                      for src in C.T0_SRC) + " |")

# break-even: eps_S of the signal in the SR; c and h, the test's response to a PU-tagged
# and to an HS+HS pair, from the signal's own pairs at pair100 (thousands of each; the
# backgrounds have a few), assuming that response does not depend on the process
S = samples[SIG]
be_stage = "pair100"
pool = {src: (1 - kept_veto(S, S["stage"][be_stage], src, S["pu_tag"]), 1 - kept_veto(S, S["stage"][be_stage], src, ~S["pu_tag"]))
        for src in C.T0_SRC}
eps_s = {src: kept_veto(S, S["stage"]["SR"], src) for src in C.T0_SRC}
if BKG:
    md += ["", f"**Break-even PU-tagged background fraction f\\*** (eps_S: the signal in the SR; c, h: the test on the "
           f"signal's own PU-tagged and HS+HS pairs at {be_stage})", "", "| t0 | eps_S | c | h | f* |", "|---|---|---|---|---|"]
    for src in C.T0_SRC:
        c, h = pool[src]
        fstar = (1 - h - eps_s[src] ** 2) / (c - h) if c > h else np.nan
        md.append(f"| {src} | {eps_s[src]:.3f} | {c:.3f} | {h:.3f} | {100 * fstar:.0f}% |")
    md += ["", "Measured PU-tagged fraction f (weighted %; MC events in brackets):", "",
           "| sample | " + " | ".join(st for st, _ in STAGES) + " |", "|---|" + "---|" * len(STAGES)]
    for n in BKG:
        s = samples[n]
        md.append(f"| {n} | " + " | ".join(f"{100 * frac(s, s['stage'][st], s['pu_tag']):.0f} ({int((s['stage'][st] & s['pu_tag']).sum())})"
                                           if s["stage"][st].any() else "-" for st, _ in STAGES) + " |")
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
open(args.out + "_summary.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
print("wrote", args.out + "_summary.md")

# ── Figures ──────────────────────────────────────────────────────────────────
ROOT = C.setup_root(args.out)
KEEP = C.KEEP


def stacked(stem, names, cats, value, counts, ytitle, title_lines):
    """one bar per sample, stacked categories (value(sample, k) = fraction), MC counts above"""
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 900, 650); KEEP.append(c)
    c.SetLeftMargin(0.11); c.SetRightMargin(0.32); c.SetBottomMargin(0.16)
    stack = ROOT.THStack("st_" + os.path.basename(stem), ""); KEEP.append(stack)
    leg = ROOT.TLegend(0.70, 0.30, 0.99, 0.30 + 0.06 * len(cats)); ROOT.StyleLegend(leg, 0.032); KEEP.append(leg)
    hs = []
    for k, (lab, col) in enumerate(cats):
        h = ROOT.TH1D(f"h_{os.path.basename(stem)}_{k}", "", len(names), 0, len(names)); KEEP.append(h)
        for b, n in enumerate(names):
            v = value(n, k)
            h.SetBinContent(b + 1, 0.0 if not np.isfinite(v) else v)
            h.GetXaxis().SetBinLabel(b + 1, n)
        h.SetFillColor(C.colour(col)); h.SetLineColor(ROOT.kWhite); h.SetLineWidth(1)
        stack.Add(h); hs.append(h)
    stack.Draw("HIST")
    stack.SetMaximum(1.2); stack.SetMinimum(0.0)
    stack.GetYaxis().SetTitle(ytitle); stack.GetXaxis().SetLabelSize(0.045); stack.GetYaxis().SetTitleOffset(1.0)
    for h, (lab, _) in reversed(list(zip(hs, cats))): leg.AddEntry(h, lab, "f")
    leg.Draw()
    t = ROOT.TLatex(); t.SetTextFont(42); t.SetTextSize(0.035); t.SetTextAlign(21); KEEP.append(t)
    for b, n in enumerate(names):
        nmc, ne = counts[b]
        t.DrawLatex(b + 0.5, 1.07, f"{nmc / 1000:.0f}k MC" if nmc >= 10000 else f"{nmc:,} MC")
        if abs(ne - nmc) > 0.05 * max(nmc, 1): t.DrawLatex(b + 0.5, 1.025, f"(eff. {ne:,.0f})")
    ROOT.ATLASLabel(0.70, 0.90, "")
    a = ROOT.TLatex(); a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.04); KEEP.append(a)
    a.DrawLatex(0.70, 0.85, "Simulation Internal")
    g = ROOT.TLatex(); g.SetNDC(); g.SetTextFont(42); g.SetTextSize(0.028); g.SetTextColor(ROOT.kGray + 2); KEEP.append(g)
    for i, ln in enumerate(title_lines):
        g.DrawLatex(0.70, 0.79 - 0.035 * i, ln)
    c.Print(stem + ".pdf"); c.Print(stem + ".png")
    print("wrote", stem + ".pdf / .png")


names = list(samples)
for st, tag in [("pair100", "pair100"), ("+dphi", "dphi"), ("SR", "sr")]:
    stacked(f"{args.out}_composition_{tag}", names, COMPOSITION,
            lambda n, k: frac(samples[n], samples[n]["stage"][st], samples[n]["comp"][COMPOSITION[k][0]]),
            [(int(samples[n]["stage"][st].sum()), neff(samples[n], samples[n]["stage"][st])) for n in names],
            "Tagging-pair composition (weighted)",
            [{"SR": "VBF H#rightarrowinv. SR selection,", "+dphi": "VBF H#rightarrowinv. SR tagging pair,",
              "pair100": "VBF H#rightarrowinv. SR tagging pair,"}[st],
             {"SR": "full (ANA-EXOT-2020-11)", "+dphi": "fJVT, |#Delta#phi_{jj}| < 2; no jet veto",
              "pair100": "E_{T}^{miss} > 100 GeV only"}[st],
             "multijet: E_{T}^{miss} from jets", "JZ slices not normalised"])
PUL = [n for n in BKG if (samples[n]["stage"]["pair100"] & samples[n]["pu_tag"]).any()] + [SIG]
stacked(f"{args.out}_puleg", PUL, PULEG,
        lambda n, k: frac(samples[n], samples[n]["stage"]["pair100"] & samples[n]["pu_tag"], samples[n]["puleg"] == k),
        [(int((samples[n]["stage"]["pair100"] & samples[n]["pu_tag"]).sum()),
          neff(samples[n], samples[n]["stage"]["pair100"] & samples[n]["pu_tag"])) for n in PUL],
        "PU-tagged events (weighted)",
        ["events with a pileup tagging jet,", "SR tagging pair, E_{T}^{miss} > 100 GeV",
         f"t_{{0}}: {args.t0.upper() if args.t0 != 'waves' else 'WAVeS'}, test at {args.nsigma:g}#sigma",
         "most catchable pileup leg"])

if BKG:
    # break-even curves: S/sqrt(B) change vs f, per t0
    c = ROOT.TCanvas("c_be", "", 800, 650); KEEP.append(c)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.05)
    fr = ROOT.TH1D("fr_be", "", 10, 0.0, 100.0); KEEP.append(fr)
    fr.SetMinimum(-15); fr.SetMaximum(40)
    fr.GetXaxis().SetTitle("Background with a pileup tagging jet  f [%]")
    fr.GetYaxis().SetTitle("S/#sqrt{B} change, event veto [%]")
    fr.Draw("AXIS")
    zero = ROOT.TLine(0, 0, 100, 0); zero.SetLineStyle(2); zero.SetLineColor(ROOT.kGray + 1); zero.Draw(); KEEP.append(zero)
    leg = ROOT.TLegend(0.17, 0.60, 0.55, 0.80); ROOT.StyleLegend(leg, 0.03); KEEP.append(leg)
    fgrid = np.linspace(0.0, 1.0, 101)
    for src in ["trkptz", "waves", "hgtd", "truth"]:
        cc, hh = pool[src]
        z = 100 * (eps_s[src] / np.sqrt(fgrid * (1 - cc) + (1 - fgrid) * (1 - hh)) - 1)
        g = ROOT.TGraph(len(fgrid), 100 * fgrid, z); KEEP.append(g)
        lab, col, _ = C.STYLE[f"t0_{src}"]
        g.SetLineColor(C.colour(col)); g.SetLineWidth(3); g.Draw("L SAME")
        leg.AddEntry(g, lab.replace("Jet vs t_{0}: ", "t_{0}: ") + f"  (f* = {100 * (1 - hh - eps_s[src] ** 2) / (cc - hh):.0f}%)", "l")
    leg.Draw()
    # each background, measured directly: its own eps_B under the event veto (t0 = --t0);
    # open marker at pair100, filled at +dphi (the last stage before the jet veto); the
    # signal's eps_S taken at the same stage
    leg2 = ROOT.TLegend(0.56, 0.20, 0.93, 0.21 + 0.045 * (len(BKG) + 1)); ROOT.StyleLegend(leg2, 0.027); KEEP.append(leg2)
    t0name = {"waves": "WAVeS", "trkptz": "TRKPTZ", "hgtd": "HGTD", "tzp": "TZP", "truth": "truth"}[args.t0]
    leg2.SetHeader(f"measured, {t0name} t_{{0}}:", "")
    for k, n in enumerate(BKG):
        s = samples[n]
        filled, hollow = [20, 21, 22, 23, 33, 34][k % 6], [24, 25, 26, 32, 27, 28][k % 6]
        for st, mk in [("pair100", hollow), ("+dphi", filled)]:
            m = s["stage"][st]
            if not m.any(): continue
            f = frac(s, m, s["pu_tag"])
            eb, eerr, _, _ = C.eff(s["w"][m], ~s["vetoed"][args.t0][m])
            if not eb > 0: continue
            g = ROOT.TGraphErrors(1); KEEP.append(g)
            es = kept_veto(S, S["stage"][st], args.t0)
            g.SetPoint(0, 100 * f, 100 * (es / np.sqrt(eb) - 1))
            g.SetPointError(0, 100 * np.sqrt(max(f * (1 - f), 1e-12) / max(neff(s, m), 1)), 100 * es * 0.5 * eb ** -1.5 * eerr)
            g.SetMarkerStyle(mk); g.SetMarkerSize(1.5); g.SetMarkerColor(ROOT.kBlack); g.SetLineColor(ROOT.kGray + 2)
            g.Draw("P SAME")
            if st == "+dphi": leg2.AddEntry(g, n, "p")
    leg2.Draw()
    ROOT.ATLASLabel(0.17, 0.88, "")
    a = ROOT.TLatex(); a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.04); KEEP.append(a)
    a.DrawLatex(0.31, 0.88, "Simulation Internal")
    g2 = ROOT.TLatex(); g2.SetNDC(); g2.SetTextFont(42); g2.SetTextSize(0.028); g2.SetTextColor(ROOT.kGray + 2); KEEP.append(g2)
    g2.DrawLatex(0.17, 0.84, "lines: #varepsilon_{S} in the SR; the test's response to PU / HS tagging jets from the signal")
    g2.DrawLatex(0.17, 0.81, "points: each background directly (open: E_{T}^{miss} > 100 pair stage, filled: + fJVT, |#Delta#phi_{jj}| < 2)")
    stem = f"{args.out}_breakeven"
    c.Print(stem + ".pdf"); c.Print(stem + ".png")
    print("wrote", stem + ".pdf / .png")
