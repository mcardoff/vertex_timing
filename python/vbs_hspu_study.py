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

The break-even PU-tagged share f*, and S/sqrt(B), are in python/vbs_sr_minimal.py:
there S and B are always measured in the same selection (the minimal set).

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_hspu_study.py \\
        --sig condor/vbf/vbf_sr_vbs_time_veto.root --zjets condor/zjets/zjets_sr_vbs_time_veto.root \\
        --dijet condor/dijet/dijet_sr_vbs_time_veto.root --out figs/time_veto/hspu

Writes <out>_summary.md and <out>_{composition_pair100,composition_dphi,composition_sr,puleg}.(pdf|png).
"""
import argparse, os
import numpy as np
import vbs_time_veto_common as C
from vbf_hinv_sr import SIG, take, load_samples, prepare

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


samples = load_samples(args.sig, args.zjets, args.dijet, args.slices)
BKG = [k for k in samples if k != SIG]
calib = C.measure(samples[SIG]["e"], True, "the VBF file's HS legs")

for name, s in samples.items():
    prepare(s, calib, args.nsigma, args.min_trk)

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
    cls = np.full(len(s["w"]), -1)
    for i, p in zip([i1, i2], s["pu_leg"]):
        aeta, timed, f = np.abs(take(J.jet_eta, i)), take(J.jet_n, i) >= args.min_trk, take(s["fail"][args.t0], i)
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

S = samples[SIG]
if BKG:
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
