#!/usr/bin/env python3
"""TDR Sec. 3.4.1 figures from util/tdr_vbf's event tree, in the ATLAS style (python/tdr_root_style.py).

The TDR preselection: >= 2 jets, leading pT > 75 GeV, second > 50 GeV, |Delta eta(j1, j2)| > 3, after
"all jets are required to pass an ITk-only R_pT pileup jet tagger operating at 85% hard-scatter
efficiency". The tagger is applied here as a jet cleaning: a jet failing it is dropped, and the
selection runs on the jets that remain (jets beyond |eta| = 4.0 have no tracks and are kept). A timing
tagger replaces the ITk one at its own 85% working point, from tdr_pujet_plots.py's json.

  <out>_cat_mjj   TDR Fig 3.24: the central-central / central-forward / forward-forward fractions
                  (forward: |eta| > 2.4) among the events with m_jj above each threshold -- normalised per
                  point, as the TDR's -- ITk-only tagger; signal solid, background (--bkg) dotted; lower
                  panel: the events above each threshold
  <out>_cat_mjj_stack_<sig|bkg>  the same categories per m_jj bin, each column normalised to 1, the
                  event count underneath (python/vbs_region_stack.py's logic)
  <out>_truth_mjj_stack_<sig|bkg>  the same stack by the truth identity of the VBF pair: HS + HS,
                  HS + PU, PU + PU, or a leg with neither TDR label
  <out>_reach_mjj_stack_<sig|bkg>  the same stack by WHERE the pair's pileup leg is: central (ITk only),
                  inside 2.4 < |eta| < 4.0 (the tracker and HGTD: timing can act), or only beyond |eta| = 4.0
                  (no tracks: no R_pT tagger can act)
  <out>_sb_gain   TDR Fig 3.25 (needs --bkg): the S/B gain over the ITk-only preselection against the additional
                  pileup-jet rejection x = 1 - eps_PU(HGTD) / eps_PU(ITk), the TDR's parametrised way -- every
                  pileup leg HGTD acts on scales its event's weight by (1 - x), hard-scatter efficiency fixed
                  at 85%; for CF events only, FF only, and both, with the rejection applied only to legs inside
                  2.4 < |eta| < 4.0 (what HGTD can reach), and -- dashed, the TDR's assumption -- to every
                  forward leg; the measured x of each real t0 (> 50 GeV forward jets at 85%) marked
  <backup-out>_pupair  per category and tagger: the share of preselected events whose VBF pair (the two
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
ap.add_argument("--stack-edges", default="0,1000,2000,4000",
                help="m_jj bin edges (GeV) of the stacks, overflow folded into the last bin (user's choice 2026-11-04: 3 bins)")
ap.add_argument("--backup-out", default="",
                help="prefix for the backup figures (the pileup-in-pair ones); default: --out")
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
match_wp = WP.get("rpt_match", "dr02")
MATCH_TXT = {"ghost": None, "dr02": "R_{p_{T}} tracks within #DeltaR < 0.2"}   # the default (ghost) unlabelled
CLEAN_TXT = {"all": None, "t0": "Cleaned t_{0} only", "off": "Track times as in the ntuple"}   # the default (all) unlabelled
clean_line = [t for t in (CLEAN_TXT[clean_wp], MATCH_TXT[match_wp]) if t]   # the non-default modes, one line each



def load(path):
    t = uproot.open(path)
    E = t["events"].arrays(library="np")
    meta = t["meta"].arrays(library="np")
    g = float(np.unique(meta["gate_sigma"])[0])
    cl = str(meta["time_cleaning"][0]) if "time_cleaning" in meta else "off"
    mt = str(meta["rpt_match"][0]) if "rpt_match" in meta else "dr02"
    if mt != match_wp:
        raise SystemExit(f"{path}: R_pT matching '{mt}', but the working points were made with '{match_wp}'")
    if abs(g - gate_wp) > 1e-9:
        raise SystemExit(f"{path}: gate {g} sigma, but the working points were made at {gate_wp}")
    if cl != clean_wp:
        raise SystemExit(f"{path}: track-time cleaning '{cl}', but the working points were made with '{clean_wp}'")
    return E


def thresholds(col):
    w = WP["working_points"][col]
    return {(r, s): w[f"{r}_{s}"]["threshold"] for r in ("central", "forward") for s in ("30_50", "50_inf")}


def select(E, col):
    """per event: passes, m_jj, category, pileup jet in the pair, one beyond |eta| = 4.0, and how many
    of the pair are hard-scatter / pileup by the TDR labels -- after dropping the jets failing `col`"""
    n = len(E["weight"])
    ok = np.zeros(n, bool); mjj = np.zeros(n); cat = np.full(n, -1); pupair = np.zeros(n, bool); pu40 = np.zeros(n, bool)
    nhs = np.zeros(n, int); npu = np.zeros(n, int)
    pu_cen = np.zeros(n, int); pu_reach = np.zeros(n, int); pu_far = np.zeros(n, int)
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
        nhs[i] = int(E["jet_hs"][i][a]) + int(E["jet_hs"][i][b])
        npu[i] = int(E["jet_pu"][i][a]) + int(E["jet_pu"][i][b])
        for j in (a, b):                                      # the pair's pileup legs, by where they are
            if not E["jet_pu"][i][j]: continue
            aej = abs(eta[j])
            if aej < 2.4: pu_cen[i] += 1                      # central: ITk tracks, no HGTD
            elif aej < 4.0: pu_reach[i] += 1                  # forward, inside the tracker and HGTD
            else: pu_far[i] += 1                              # beyond |eta| = 4.0: no tracks at all
    return ok, mjj, cat, pupair, pu40, nhs, npu, pu_cen, pu_reach, pu_far


samples = [("sig", args.sig_label, load(args.sig), 1)]
if args.bkg: samples.append(("bkg", args.bkg_label, load(args.bkg), 3))
R = {}
for key, lab, E, _ in samples:
    for col, *_ in TAGGERS:
        R[key, col] = select(E, col)
    # the all-jets selection must reproduce tdr_vbf's own flag
    assert np.array_equal(R[key, None][0], E["presel_all"]), "offline preselection disagrees with presel_all"
sample_line = lambda lab: lab + (f" ({args.note})" if args.note else "")
gate_txt = "" if gate_wp == 2.0 else f", t_{{0}} gate {gate_wp:g}#sigma"   # the TDR's 2 sigma is the standard: unlabelled

# ── Fig 3.24 analogue: category fractions among the events above each m_jj threshold ──
# Normalised PER POINT, as the TDR's: at each threshold the CC / CF / FF fractions of the events with
# m_jj above it sum to 1 (the TDR's forward-forward fraction rises with m_jj). Lower panel: the events
# passing each threshold, i.e. the statistics behind every point.
CC = {"CC": S.COL["cc"], "CF": S.COL["cf"], "FF": S.COL["ff"]}  # the TDR's black / red / blue, in the usual palette
CNAME = {"CC": "Central-central", "CF": "Central-forward", "FF": "Forward-forward"}
X = np.linspace(0, 2000, 201)
MAIN, RAT = 0.68, 0.32
TS_M, TS_R = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)
cv = S.keep(R_.TCanvas("catmjj", "", 800, 800))
cv.cd(); S.pad("catmjj_m", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
fr = S.keep(R_.gPad.DrawFrame(0.0, 0.0, 2000.0, 1.4 if len(samples) == 1 else 1.6))
S.axes(fr, "", "Fraction of events", size=TS_M, ytoff=1.35)
fr.GetXaxis().SetLabelSize(0); fr.GetYaxis().ChangeLabel(1, -1, 0)
leg = S.legend(0.20, 0.555, 0.93, 0.615, size=0.75 * TS_M, ncol=3, margin=0.3)   # one row under the header
counts = {}
print("Fig 3.24 analogue (ITk-only tagger): category fractions among the events above each m_jj threshold")
for key, lab, E, ls in samples:
    ok, mjj, cat, *_ = R[key, "rpt_itk"]
    w = E["weight"].astype(float)
    above = np.array([w[ok & (mjj > x)].sum() for x in X])
    counts[key] = np.array([np.sum(ok & (mjj > x)) for x in X])
    for c, name in enumerate(CATS):
        f = np.array([w[ok & (cat == c) & (mjj > x)].sum() for x in X]) / np.where(above > 0, above, np.nan)
        g = S.style(S.graph(X, f), CC[name], 0, ls, 3); g.Draw("L SAME")
        if key == "sig": leg.AddEntry(g, CNAME[name], "L")
    for x in (0, 500, 1000, 1500):
        i = int(np.argmin(np.abs(X - x))); m = ok & (mjj > X[i])
        print(f"  {plain(lab)} m_jj > {X[i]:.0f}: {counts[key][i]} events;  " + "  ".join(f"{name} {100 * w[m & (cat == c)].sum() / w[m].sum():.1f}%" for c, name in enumerate(CATS)))
leg.Draw()
if len(samples) > 1:   # the samples: a second row under the categories
    leg2 = S.legend(0.20, 0.49, 0.93, 0.55, size=0.75 * TS_M, ncol=2, margin=0.2)
    for key, lab, E, ls in samples:
        d = S.keep(R_.TGraph(1)); d.SetLineStyle(ls); d.SetLineWidth(3); d.SetLineColor(R_.kGray + 3)
        leg2.AddEntry(d, lab, "L")
    leg2.Draw()
S.header(0.20, 0.87, [ENERGY + (", " + sample_line(samples[0][1]) if len(samples) == 1 else ""),
                      "VBF preselection, ITk-only R_{p_{T}} at 85%", "Forward: |#eta| > 2.4"] + clean_line, size=TS_M)
R_.gPad.RedrawAxis()
cv.cd(); S.pad("catmjj_r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
R_.gPad.SetLogy()
cmax = max(c.max() for c in counts.values())
fr2 = S.keep(R_.gPad.DrawFrame(0.0, 3.0, 2000.0, 30.0 * cmax))
S.axes(fr2, "m_{jj} threshold [GeV]", "Events above", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2)
for key, lab, E, ls in samples:
    S.style(S.graph(X, counts[key].astype(float)), "#555555", 0, ls, 2).Draw("L SAME")
S.save(cv, f"{args.out}_cat_mjj")

# ── companions: per m_jj bin, each column normalised to 1 (vbs_region_stack.py's logic) ──
# A column per m_jj bin, normalised to 1, the event count underneath. Two band sets: the TDR
# preselection's CC / CF / FF categories, and the truth identity of the VBF pair (the two leading
# jets left by the ITk-only tagger) by the TDR labels -- HS within dR < 0.3 of a truth jet > 10 GeV,
# PU beyond dR > 0.6 of every truth jet > 4 GeV, and a leg that is neither.
EDGES = np.array([float(x) for x in args.stack_edges.split(",")])


def mjj_stack(key, lab, E, ok, mjj, bands, name, title, sel="VBF preselection, ITk-only R_{p_{T}} at 85%"):
    """bands: (legend, mask, colour), bottom to top; the masks must partition the selected events"""
    w = E["weight"].astype(float)
    xm = np.clip(mjj, EDGES[0], EDGES[-1] - 1e-6)            # overflow folded into the top bin
    tot = np.histogram(xm[ok], EDGES, weights=w[ok])[0]
    nraw = np.histogram(xm[ok], EDGES)[0]
    assert sum(int(np.sum(ok & m)) for _, m, _ in bands) == int(ok.sum()), f"{name}: bands do not partition the events"
    cv = S.keep(R_.TCanvas(f"{name}_{key}", "", 800, 800))
    cv.cd(); S.pad(f"{name}_m_{key}", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
    st = S.keep(R_.THStack(f"st_{name}_{key}", ""))
    hs = []
    for k, (leg_txt, m, colour) in enumerate(bands):
        h = S.keep(R_.TH1D(f"{name}_{key}_{k}", "", len(EDGES) - 1, EDGES))
        part = np.histogram(xm[ok & m], EDGES, weights=w[ok & m])[0]
        for b in range(len(EDGES) - 1):
            h.SetBinContent(b + 1, part[b] / tot[b] if tot[b] > 0 else 0.0)
        h.SetFillColor(S.colour(colour)); h.SetLineColor(R_.kBlack); h.SetLineWidth(1)
        st.Add(h); hs.append((leg_txt, h))
    fr = S.keep(R_.gPad.DrawFrame(EDGES[0], 0.0, EDGES[-1], 1.9))
    S.axes(fr, "", "Fraction of events", size=TS_M, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0); fr.GetYaxis().ChangeLabel(1, -1, 0)
    st.Draw("HIST SAME")
    nrow = (len(hs) + 1) // 2
    leg = S.legend(0.20, 0.525, 0.93, 0.525 + 0.055 * nrow, size=0.7 * TS_M, ncol=2, margin=0.15)
    for leg_txt, h in reversed(hs): leg.AddEntry(h, leg_txt, "F")
    leg.Draw()
    S.header(0.20, 0.87, [ENERGY + ", " + sample_line(lab), sel] + ([title] if title else []) + clean_line, size=TS_M)
    R_.gPad.RedrawAxis()
    cv.cd(); S.pad(f"{name}_r_{key}", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
    R_.gPad.SetLogy()
    hn = S.keep(R_.TH1D(f"{name}_n_{key}", "", len(EDGES) - 1, EDGES))
    for b in range(len(EDGES) - 1): hn.SetBinContent(b + 1, nraw[b])
    hn.SetFillColor(R_.TColor.GetColor("#DCE3EE")); hn.SetLineColor(R_.kBlack)
    hn.SetMinimum(max(1.0, 0.3 * nraw[nraw > 0].min())); hn.SetMaximum(30.0 * nraw.max())
    S.axes(hn, "m_{jj} [GeV]   (last bin includes overflow)", "Events", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2)
    hn.Draw("HIST")
    t = R_.TLatex(); t.SetTextSize(0.55 * TS_R); t.SetTextAlign(21)
    for b in range(len(EDGES) - 1):
        S.keep(t.DrawLatex(0.5 * (EDGES[b] + EDGES[b + 1]), nraw[b] * 1.8, f"{nraw[b]}"))
    print(f"  {plain(lab)} {name}, per m_jj bin (" + " / ".join(b[0].replace('#', '') for b in bands) + "): " +
          "  ".join(f"[{EDGES[b]:.0f},{EDGES[b+1]:.0f}) " + "/".join(f"{100 * w[ok & m & (xm >= EDGES[b]) & (xm < EDGES[b + 1])].sum() / tot[b]:.0f}" for _, m, _ in bands)
                    for b in range(len(EDGES) - 1) if tot[b] > 0))
    S.save(cv, f"{args.out}_{name}_{key}")


# truth-identity colours: python/vbs_region_stack.py's
IDENT = [("HS + HS (genuine pair)", "#7BCB6B"), ("HS + PU", "#F5A02A"), ("PU + PU", "#636363"), ("A leg carries neither label", "#C79BF5")]
# pileup legs by reach (an event with a leg inside 2.4 < |eta| < 4.0 counts as reachable first)
# |eta| < 2.4: ITk tracks only; 2.4-4.0: ITk and HGTD; >= 4.0: no tracks, no R_pT tagger can act
REACH = [("HS + HS", "#7BCB6B"), ("PU leg at |#eta| < 2.4", "#B3231C"),
         ("PU leg at 2.4#minus4.0 (HGTD)", "#3F90DA"), ("PU leg only at |#eta| #geq 4.0", "#636363"),
         ("A leg has neither label", "#C79BF5")]
for key, lab, E, ls in samples:
    ok, mjj, cat, _, _, nhs, npu, pu_cen, pu_reach, pu_far = R[key, "rpt_itk"]
    mjj_stack(key, lab, E, ok, mjj, [(CNAME[n], cat == c, CC[n]) for c, n in enumerate(CATS)],
              "cat_mjj_stack", "Forward: |#eta| > 2.4")
    ident = [nhs == 2, (nhs == 1) & (npu == 1), npu == 2, (nhs + npu) < 2]
    mjj_stack(key, lab, E, ok, mjj, [(t, m, c) for (t, c), m in zip(IDENT, ident)],
              "truth_mjj_stack", "VBF pair = the two leading jets")
    # where the pair's pileup leg is (what timing can and cannot reach): after the ITk-only tagger, and
    # before any tagging -- the pair of figures shows what the ITk tagger removed and what is left
    for col, name, sel in (("rpt_itk", "reach_mjj_stack", None), (None, "reach_mjj_stack_notag", "VBF preselection, no pileup-jet tagging")):
        ok_, mjj_, _, _, _, nhs_, npu_, pc_, pr_, pf_ = R[key, col]
        reach = pr_ > 0
        far = ~reach & (pf_ > 0)
        cen = ~reach & ~far & (pc_ > 0)
        bands = [(npu_ == 0) & (nhs_ == 2), cen, reach, far, (npu_ == 0) & (nhs_ < 2)]
        mjj_stack(key, lab, E, ok_, mjj_, [(t, m, c) for (t, c), m in zip(REACH, bands)], name, None,
                  **({"sel": sel} if sel else {}))

# ── pileup jets in the VBF pair, per category and tagger ──────────────────────
LABELS = ["Central-central", "Central-forward", "Forward-forward"]
for key, lab, E, _ in samples:
    w = E["weight"].astype(float)
    ref = R[key, "rpt_itk"]
    rows = []
    print(f"\n{plain(lab)}: per category -- events (weighted share with a pileup jet in the pair; of which beyond |eta| 4.0) [relative to ITk only]")
    print("  " + f"{'tagger':40s}" + "".join(f"{c:>30s}" for c in CATS))
    for col, name, colour, mk in TAGGERS:
        ok, _, cat, pp, p40, *_ = R[key, col]
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
    S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200", sample_line(lab), "R_{p_{T}} taggers at 85%" + gate_txt] + clean_line)
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
    S.save(cv, f"{args.backup_out or args.out}_pupair_{key}")

# ── TDR Fig 3.25 analogue: the S/B gain against the additional pileup-jet rejection ──
# The TDR's parametrised method (Sec. 3.4.1, PDF p. 78): with the ITk-only preselection at 85% hard-scatter
# efficiency, events are reweighted by the pileup-jet efficiency gain relative to ITk only: every pileup leg
# of the VBF pair that HGTD acts on scales its event's weight by (1 - x), x = 1 - eps_PU(HGTD) / eps_PU(ITk);
# the hard-scatter efficiency is held, so the signal changes only through its own pileup legs.
# Gain = (S(x) / S(0)) / (B(x) / B(0)). A removed leg drops the event (a different pair might pass instead:
# ignored, as the TDR's reweighting ignores it). The x axis follows the TDR's FIGURE (S/B = 1 at x = 0); its
# text describes the complement (S/B = 1 "when the x-axis is unity"), which the figure does not show.
# Two scopes: legs inside 2.4 < |eta| < 4.0 (HGTD's acceptance), and every forward leg (|eta| > 2.4, up to
# the jets' 4.5). The TDR does not say which legs it reweighted; its curves imply HGTD acting on ~half the
# Z+jets events after the ITk tagger, which our sample reaches only with the legs beyond |eta| = 4.0.
# ATLAS-TDR-031 Fig 3.25 (PDF p. 79), digitised from the page rendered at 200 dpi (+- 0.01):
TDR_FIG325 = {
    "both": ([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95],
             [1.0, 1.060, 1.124, 1.196, 1.274, 1.364, 1.464, 1.583, 1.711, 1.868, 1.955]),   # leaves the frame at 2
    "cf": ([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
           [1.0, 1.062, 1.097, 1.135, 1.176, 1.221, 1.271, 1.327, 1.390, 1.458]),
    "ff": ([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
           [1.0, 1.046, 1.065, 1.084, 1.107, 1.125, 1.142, 1.160, 1.175, 1.187])}
if len(samples) > 1:
    XS = np.linspace(0.0, 1.0, 101)
    CURV = [("cf", [1], "CF only", S.COL["cf"]), ("ff", [2], "FF only", S.COL["ff"]), ("both", [1, 2], "CF and FF", "#000000")]
    SCOPES = [("reach", "Pileup legs in 2.4 < |#eta| < 4.0 (HGTD acceptance)"),
              ("all", "Every forward pileup leg, incl. |#eta| #geq 4.0")]
    rng = np.random.default_rng(20261104)
    NBOOT = 400

    def leg_sums(key, scope, cats):
        """weight sums of the selected events with 0, 1, 2 affected pileup legs"""
        ok, mjj, cat, _, _, nhs, npu, pu_cen, pu_reach, pu_far = R[key, "rpt_itk"]
        legs = pu_reach + (pu_far if scope == "all" else 0)
        legs = np.where(np.isin(cat, cats), legs, 0)[ok]
        w = dict((k_, E_) for k_, _, E_, _ in samples)[key]["weight"].astype(float)[ok]
        return legs, w

    def gain_curve(scope, cats, boot=False):
        out = []
        for key in ("sig", "bkg"):
            legs, w = leg_sums(key, scope, cats)
            if boot:
                idx = rng.integers(0, len(w), (NBOOT, len(w)))
                W = np.stack([np.sum(np.where(legs[idx] == k, w[idx], 0.0), axis=1) for k in (0, 1, 2)], axis=1)
            else:
                W = np.array([[w[legs == k].sum() for k in (0, 1, 2)]])
            out.append((W[:, :1] + W[:, 1:2] * (1 - XS) + W[:, 2:3] * (1 - XS) ** 2) / W.sum(1, keepdims=True))
        return out[0] / out[1]                       # (1 or NBOOT, len(XS))

    pe = lambda col: WP["working_points"][col]["forward_50_inf"]["pu_eff"]
    XT = [(col, name, colour, mk, 1.0 - pe(col) / pe("rpt_itk")) for col, name, colour, mk in TAGGERS if col not in (None, "rpt_itk")]
    YMAX = 2.1
    cv = S.keep(R_.TCanvas("sbgain", "", 1600, 700))
    print("\nFig 3.25 analogue: S/B gain at x = 0.2 / 0.5 / 1.0 (+- bootstrap 68%); events with an affected leg (sig / bkg)")
    for c, (scope, stxt) in enumerate(SCOPES):
        cv.cd(); S.pad(f"sbg{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0, top=0.05, bottom=0.14, left=0.14, right=0.04)
        fr = S.keep(R_.gPad.DrawFrame(0.0, 0.95, 1.0, YMAX))
        S.axes(fr, "Additional pileup-jet rejection  1 #minus #varepsilon_{PU}^{HGTD} / #varepsilon_{PU}^{ITk}",
               "Normalised #varepsilon_{S} / #varepsilon_{B}", size=0.045, ytoff=1.35, xtoff=1.3)
        boot = gain_curve(scope, [1, 2], boot=True)
        lo, hi = np.percentile(boot, 16, axis=0), np.percentile(boot, 84, axis=0)
        band = S.keep(R_.TGraph(2 * len(XS)))
        for k, (xv, yv) in enumerate(list(zip(XS, hi)) + list(zip(XS[::-1], lo[::-1]))): band.SetPoint(k, xv, yv)
        band.SetFillColor(S.colour("#C8C8C8")); band.SetLineWidth(0); band.Draw("F SAME")
        for name, cats, txt, colour in CURV:
            g = S.style(S.graph(*TDR_FIG325[name]), colour, 0, 2, 2); g.Draw("L SAME")
        ours = {}
        for name, cats, txt, colour in CURV:
            ours[name] = gain_curve(scope, cats)[0]
            g = S.style(S.graph(XS, ours[name]), colour, 0, 1, 3); g.Draw("L SAME")
            legs_s, _ = leg_sums("sig", scope, cats); legs_b, _ = leg_sums("bkg", scope, cats)
            if name == "both":
                err = lambda xv: 0.5 * (np.interp(xv, XS, hi) - np.interp(xv, XS, lo))
                vals = " / ".join(f"{np.interp(xv, XS, ours[name]):.3f} +- {err(xv):.3f}" for xv in (0.2, 0.5, 1.0))
            else:
                vals = " / ".join(f"{np.interp(xv, XS, ours[name]):.3f}" for xv in (0.2, 0.5, 1.0))
            print(f"  {plain(stxt):52s} {txt:10s} {vals}   ({int((legs_s > 0).sum())} / {int((legs_b > 0).sum())})")
        for col, name, colour, mk, xv in XT:
            m = S.style(S.graph([xv], [float(np.interp(xv, XS, ours["both"]))]), colour, mk, 1, 1, 1.7); m.Draw("P SAME")
        S.header(0.19, 0.88, [ENERGY, "S: VBF H#rightarrowinv., B: Z+jets, ITk-only R_{p_{T}} at 85%", stxt] + clean_line, size=0.038)
        if c == 0:
            leg = S.legend(0.17, 0.48, 0.58, 0.71, size=0.034, margin=0.22)
            for name, cats, txt, colour in CURV:
                e = S.keep(R_.TLine()); e.SetLineColor(S.colour(colour)); e.SetLineWidth(3); leg.AddEntry(e, txt, "L")
            e1 = S.keep(R_.TGraph(1)); e1.SetLineColor(R_.kGray + 2); e1.SetLineWidth(3); e1.SetFillColor(S.colour("#C8C8C8"))
            leg.AddEntry(e1, "This study #pm 1#sigma stat.", "LF")
            e2 = S.keep(R_.TLine()); e2.SetLineColor(R_.kGray + 2); e2.SetLineWidth(2); e2.SetLineStyle(2); leg.AddEntry(e2, "ATLAS-TDR-031", "L")
            leg.Draw()
        else:
            leg = S.legend(0.17, 0.52, 0.58, 0.71, size=0.034, margin=0.22)
            leg.SetHeader("Measured x: forward jets > 50 GeV")
            for col, name, colour, mk, xv in XT:
                m = S.style(S.graph([0], [0]), colour, mk, 1, 1, 1.7); leg.AddEntry(m, f"{name}  ({xv:.2f})", "P")
            leg.Draw()
        R_.gPad.RedrawAxis()
    S.save(cv, f"{args.out}_sb_gain")
    print("  measured x, forward > 50 GeV at 85%: " + ", ".join(f"{plain(n)} {x:.2f}" for _, n, _, _, x in XT))
