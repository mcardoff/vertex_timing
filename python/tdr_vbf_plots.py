#!/usr/bin/env python3
"""TDR Sec. 3.4.1 figures from util/tdr_vbf's event tree, in the ATLAS style (python/tdr_root_style.py).

The TDR preselection: >= 2 jets, leading pT > 75 GeV, second > 50 GeV, |Delta eta(j1, j2)| > 3, after
"all jets are required to pass an ITk-only R_pT pileup jet tagger operating at 85% hard-scatter
efficiency". The tagger is applied here as a jet cleaning: a jet failing it is dropped, and the
selection runs on the jets that remain (jets beyond |eta| = 4.0 have no tracks and are kept). A timing
tagger replaces the ITk one at its own 85% working point, from tdr_pujet_plots.py's json.

  <out>_cat_mjj   TDR Fig 3.24: the fraction of preselected events with m_jj above a threshold (dashed)
                  and its split into central-central / central-forward / forward-forward (forward:
                  |eta| > 2.4), ITk-only tagger; signal solid, background (--bkg) dotted
  <out>_pupair    per category and tagger: the share of preselected events whose VBF pair (the two
                  leading surviving jets) contains a pileup jet (TDR labels) -- hatched where that
                  jet is at |eta| >= 4.0, beyond tracking and HGTD, where no R_pT tagger can act --
                  (left) and the preselected events relative to the ITk-only tagger (right)

Events are weighted by `weight`; errors are binomial on the effective entries.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_vbf_plots.py --sig ../figs/hists/tdr_vbf.root \\
        --wp ../figs/tdr/local_tdr_pujet_wp85.json --note "local sample" --out ../figs/tdr/local_tdr_vbf
"""
import argparse, json, os
import numpy as np
import uproot
import tdr_root_style as S

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True, help="tdr_vbf output, VBF H->inv")
ap.add_argument("--bkg", help="tdr_vbf output, Z+jets (optional)")
ap.add_argument("--wp", required=True, help="tdr_pujet_plots.py's <out>_wp85.json")
ap.add_argument("--sig-label", default="VBF H#rightarrowinv.", help="ROOT TLatex")
ap.add_argument("--bkg-label", default="Z+jets")
ap.add_argument("--note", default="", help="appended to the sample line, e.g. 'local sample'")
ap.add_argument("--out", required=True)
args = ap.parse_args()
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
R_ = S.setup()

J1, J2, DETA = 75.0, 50.0, 3.0
# (key, legend, colour, ROOT marker) -- tdr_pujet_plots.py's colours; None = no pileup tagging at all
TAGGERS = [(None, "No pileup jet tagging", S.COL["none"], 0),
           ("rpt_itk", "ITk only (TDR)", S.COL["itk"], 24),
           ("rpt_t0_hgtd", "HGTD (Athena) t_{0}", S.COL["hgtd"], 23),
           ("rpt_t0_trkptz", "TRKPTZ t_{0}", S.COL["trkptz"], 21),
           ("rpt_t0_waves", "WAVeS t_{0}", S.COL["waves"], 20),
           ("rpt_t0_truth", "Truth t_{0}", S.COL["truth"], 25)]
plain = lambda s: s.replace("_{0}", "0").replace("#", "")
ENERGY = "#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200"
CATS = ["CC", "CF", "FF"]
WP = json.load(open(args.wp))
gate_wp = WP["gate_sigma"]
clean_wp = WP.get("time_cleaning", "off")
CLEAN_TXT = {"all": None, "t0": "Cleaned t_{0} only", "off": "Track times as in the ntuple"}   # the default (all) unlabelled
clean_line = [CLEAN_TXT[clean_wp]] if CLEAN_TXT[clean_wp] else []



def load(path):
    t = uproot.open(path)
    E = t["events"].arrays(library="np")
    meta = t["meta"].arrays(library="np")
    g = float(np.unique(meta["gate_sigma"])[0])
    cl = str(meta["time_cleaning"][0]) if "time_cleaning" in meta else "off"
    if abs(g - gate_wp) > 1e-9:
        raise SystemExit(f"{path}: gate {g} sigma, but the working points were made at {gate_wp}")
    if cl != clean_wp:
        raise SystemExit(f"{path}: track-time cleaning '{cl}', but the working points were made with '{clean_wp}'")
    return E


def thresholds(col):
    w = WP["working_points"][col]
    return {(r, s): w[f"{r}_{s}"]["threshold"] for r in ("central", "forward") for s in ("30_50", "50_inf")}


def select(E, col):
    """per event: passes, m_jj, category, pileup jet in the pair, one beyond |eta| = 4.0 --
    after dropping the jets failing `col`"""
    n = len(E["weight"])
    ok = np.zeros(n, bool); mjj = np.zeros(n); cat = np.full(n, -1); pupair = np.zeros(n, bool); pu40 = np.zeros(n, bool)
    th = thresholds(col) if col else None
    for i in range(n):
        pt, eta, phi = E["jet_pt"][i], E["jet_eta"][i], E["jet_phi"][i]
        keep = np.ones(len(pt), bool)
        if col:
            r = E["jet_" + col][i]
            ae = np.abs(eta)
            thr = np.where(ae < 2.4, np.where(pt < 50, th["central", "30_50"], th["central", "50_inf"]),
                           np.where(pt < 50, th["forward", "30_50"], th["forward", "50_inf"]))
            keep = (r < 0) | (r >= thr)                       # r < 0: beyond |eta| = 4.0, no tracks
        s = np.flatnonzero(keep)                              # the tree's jets are pT-ordered
        if len(s) < 2: continue
        a, b = s[0], s[1]
        if pt[a] <= J1 or pt[b] <= J2 or abs(eta[a] - eta[b]) <= DETA: continue
        ok[i] = True
        v = [pt[a] * np.cosh(eta[a]), pt[b] * np.cosh(eta[b])]
        px = pt[a] * np.cos(phi[a]) + pt[b] * np.cos(phi[b]); py = pt[a] * np.sin(phi[a]) + pt[b] * np.sin(phi[b])
        pz = pt[a] * np.sinh(eta[a]) + pt[b] * np.sinh(eta[b])
        mjj[i] = np.sqrt(max((v[0] + v[1]) ** 2 - px * px - py * py - pz * pz, 0.0))
        cat[i] = int(abs(eta[a]) > 2.4) + int(abs(eta[b]) > 2.4)
        pupair[i] = bool(E["jet_pu"][i][a] or E["jet_pu"][i][b])
        pu40[i] = any(E["jet_pu"][i][j] and abs(eta[j]) >= 4.0 for j in (a, b))
    return ok, mjj, cat, pupair, pu40


samples = [("sig", args.sig_label, load(args.sig), 1)]
if args.bkg: samples.append(("bkg", args.bkg_label, load(args.bkg), 3))
R = {}
for key, lab, E, _ in samples:
    for col, *_ in TAGGERS:
        R[key, col] = select(E, col)
    # the all-jets selection must reproduce tdr_vbf's own flag
    assert np.array_equal(R[key, None][0], E["presel_all"]), "offline preselection disagrees with presel_all"
sample_line = lambda lab: lab + (f" ({args.note})" if args.note else "")
gate_txt = f"t_{{0}} gate {gate_wp:g}#sigma"

# ── Fig 3.24 analogue ─────────────────────────────────────────────────────────
CC = {"CC": S.COL["cc"], "CF": S.COL["cf"], "FF": S.COL["ff"]}  # the TDR's black / red / blue, in the usual palette
CNAME = {"CC": "Central-central", "CF": "Central-forward", "FF": "Forward-forward"}
X = np.linspace(0, 3000, 301)
cv = S.keep(R_.TCanvas("catmjj", "", 800, 600))
fr = S.keep(cv.DrawFrame(0.0, 0.0, 3000.0, 1.8))
S.axes(fr, "m_{jj} threshold [GeV]", "Fraction of preselected events")
fr.GetYaxis().ChangeLabel(1, -1, 0)
leg = S.legend(0.58, 0.58 - 0.05 * 4 * len(samples) / (2 if len(samples) > 1 else 1), 0.93, 0.58,
               size=0.035, ncol=2 if len(samples) > 1 else 1)
print("Fig 3.24 analogue (ITk-only tagger): fraction of preselected events with m_jj above the threshold")
for key, lab, E, ls in samples:
    ok, mjj, cat, *_ = R[key, "rpt_itk"]
    w = E["weight"].astype(float)
    tot = w[ok].sum()
    frac = lambda m: np.array([w[m & (mjj > x)].sum() for x in X]) / tot
    pre = (plain(lab) if False else lab) + ", " if len(samples) > 1 else ""
    g = S.style(S.graph(X, frac(ok)), "#777777", 0, 2 if key == "sig" else 4, 3); g.Draw("L SAME")
    leg.AddEntry(g, f"{pre}All", "L")
    for c, name in enumerate(CATS):
        g = S.style(S.graph(X, frac(ok & (cat == c))), CC[name], 0, ls, 3); g.Draw("L SAME")
        leg.AddEntry(g, f"{pre}{CNAME[name]}", "L")
    print(f"  {plain(lab)}: {ok.sum()} events;  " + "  ".join(f"{name} {100 * w[ok & (cat == c)].sum() / tot:.1f}%" for c, name in enumerate(CATS))
          + "  | m_jj > 1 TeV: " + "  ".join(f"{name} {100 * w[ok & (cat == c) & (mjj > 1000)].sum() / tot:.1f}%" for c, name in enumerate(CATS)))
leg.Draw()
S.header(0.20, 0.87, [ENERGY, sample_line(", ".join(s[1] for s in samples)) + ", VBF preselection",
                      "ITk-only R_{p_{T}} tagger at 85%, forward: |#eta| > 2.4"])
cv.RedrawAxis()
S.save(cv, f"{args.out}_cat_mjj")

# ── pileup jets in the VBF pair, per category and tagger ──────────────────────
LABELS = ["Central-central", "Central-forward", "Forward-forward"]
for key, lab, E, _ in samples:
    w = E["weight"].astype(float)
    ref = R[key, "rpt_itk"]
    rows = []
    print(f"\n{plain(lab)}: per category -- events (weighted share with a pileup jet in the pair; of which beyond |eta| 4.0) [relative to ITk only]")
    print("  " + f"{'tagger':40s}" + "".join(f"{c:>30s}" for c in CATS))
    for col, name, colour, mk in TAGGERS:
        ok, _, cat, pp, p40 = R[key, col]
        frac, f40, err, rel, rerr, line = [], [], [], [], [], ""
        for c in range(3):
            m = ok & (cat == c)
            sw, sw2 = w[m].sum(), (w[m] ** 2).sum()
            f = w[m & pp].sum() / sw if sw > 0 else np.nan
            neff = sw * sw / sw2 if sw2 > 0 else 0
            frac.append(100 * f); err.append(100 * np.sqrt(f * (1 - f) / neff) if neff > 0 else np.nan)
            f40.append(100 * w[m & p40].sum() / sw if sw > 0 else np.nan)
            sref = w[ref[0] & (ref[2] == c)].sum()
            rel.append(sw / sref if sref > 0 else np.nan)
            rerr.append(np.sqrt(sw2) / sref if sref > 0 else np.nan)
            line += f"{int(m.sum()):>7d} ({100 * f:5.2f}%; {f40[-1]:5.2f}%) [{rel[-1]:5.3f}]"
        rows.append((frac, f40, err, rel, rerr))
        print("  " + f"{plain(name):40s}" + line)
    cv = S.keep(R_.TCanvas(f"pupair_{key}", "", 1600, 600))

    # left: the share of events with a pileup jet in the pair, bars per tagger, hatched beyond |eta| 4.0
    cv.cd(); S.pad(f"bars_{key}", 0.0, 0.0, 0.5, 1.0)
    ymax = 1.9 * np.nanmax([np.nanmax(r[0]) for r in rows])
    hf = S.keep(R_.TH1D(f"frame_{key}", "", 3, 0, 3))
    for b, t in enumerate(LABELS): hf.GetXaxis().SetBinLabel(b + 1, t)
    hf.SetMinimum(0); hf.SetMaximum(ymax)
    S.axes(hf, "", "Pileup jet in the VBF pair [%]", size=0.05, ytoff=1.35)
    hf.GetXaxis().SetLabelSize(0.055); hf.GetYaxis().ChangeLabel(1, -1, 0)
    hf.Draw("AXIS")
    nT = len(TAGGERS); slot = 0.84 / nT; bw = 0.86 * slot
    leg = S.legend(0.63, 0.53, 0.93, 0.92, size=0.033, margin=0.25)
    for k, ((col, name, colour, mk), (frac, f40, err, _, _)) in enumerate(zip(TAGGERS, rows)):
        off = 0.08 + k * slot + (slot - bw) / 2
        h = S.keep(R_.TH1D(f"bar_{key}_{k}", "", 3, 0, 3))
        h4 = S.keep(R_.TH1D(f"bar40_{key}_{k}", "", 3, 0, 3))
        for b in range(3): h.SetBinContent(b + 1, frac[b]); h4.SetBinContent(b + 1, f40[b])
        for hh in (h, h4): hh.SetBarWidth(bw); hh.SetBarOffset(off)
        h.SetFillColor(S.colour(colour)); h.SetLineColor(S.colour(colour)); h.Draw("BAR SAME")
        h4.SetFillColor(R_.kWhite); h4.SetFillStyle(3354); h4.SetLineWidth(0); h4.Draw("BAR SAME")
        ge = S.graph(np.arange(3) + off + bw / 2, frac, err)
        ge.SetMarkerSize(0); ge.SetLineColor(R_.kGray + 3); ge.SetLineWidth(1); ge.Draw("P SAME")
        leg.AddEntry(h, name, "F")
    hatch = S.keep(R_.TH1D(f"hatch_{key}", "", 1, 0, 1)); hatch.SetFillColor(R_.kGray + 1); hatch.SetFillStyle(3354); hatch.SetLineColor(R_.kGray + 1)
    leg.AddEntry(hatch, "Pileup jet at |#eta| #geq 4.0", "F")
    leg.Draw()
    S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200", sample_line(lab), "R_{p_{T}} taggers at 85%, " + gate_txt] + clean_line)
    R_.gPad.RedrawAxis()

    # right: preselected events relative to the ITk-only tagger
    cv.cd(); S.pad(f"rel_{key}", 0.5, 0.0, 1.0, 1.0)
    hr = S.keep(R_.TH1D(f"frame_rel_{key}", "", 3, 0, 3))
    for b, t in enumerate(LABELS): hr.GetXaxis().SetBinLabel(b + 1, t)
    hr.SetMinimum(0.94); hr.SetMaximum(1.12)
    S.axes(hr, "", "Preselected events / ITk only", size=0.05, ytoff=1.45)
    hr.GetXaxis().SetLabelSize(0.055)
    hr.Draw("AXIS")
    one = S.keep(R_.TLine(0, 1, 3, 1)); one.SetLineStyle(2); one.SetLineColor(R_.kGray + 2); one.Draw()
    leg2 = S.legend(0.63, 0.64, 0.93, 0.92, size=0.033, margin=0.25)
    for k, ((col, name, colour, mk), (_, _, _, rel, rerr)) in enumerate(zip(TAGGERS, rows)):
        if col is None: continue
        g = S.style(S.graph(np.arange(3) + 0.5 + (k - 3) * 0.13, rel, rerr), colour, mk, 1, 2, 1.4)
        g.Draw("P SAME"); leg2.AddEntry(g, name, "P")
    leg2.Draw()
    S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200", sample_line(lab)] + clean_line)
    S.text(0.20, 0.20, "No tagging: " + " / ".join(f"{v:.2f}" for v in rows[0][3]) + " (CC / CF / FF)", 0.042)
    S.save(cv, f"{args.out}_pupair_{key}")
