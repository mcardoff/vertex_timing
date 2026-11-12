#!/usr/bin/env python3
"""The clustering as a process, step by step, for the Global Performance talk's algorithm slide
(slides/slides_gp_timing.tex): a visual aid to "1. Clustering in time" (and "2. Choosing a cluster").

Replays src/clustering_functions.h's doIterativeClustering on the tracks of one event dumped by
python/clustering_dump.cxx (talk inputs: TDR-cleaned times, 2.5 sigma z association), and checks that the replay
reproduces the C++ clusters track for track before drawing. Four strips on one time axis, track pT on the y axis:
  1  the seed: the highest-pT track not yet in a cluster, and its 3 sigma window (drawn for the median track sigma;
     the real test is |t_c - t_i| / sqrt(sigma_c^2 + sigma_i^2) < 3 per track)
  2  the first absorption: the seed plus the closest track within 3 sigma, the mean updated from the two, its new
     3 sigma window, and the next track to be absorbed ringed (tracks are absorbed one at a time)
  3  the next seed, until every track is in a cluster: all clusters, each with its t +- sigma
truth t_HS is the dashed line on every strip.
  4  the choice: each cluster's WAVeS score, the chosen one outlined

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_clustering_steps.py --json ../figs/tdr/gp/gp_clustering_event.json \\
      --out ../figs/tdr/gp/gp_clustering_steps
"""
import argparse
import json

import numpy as np

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--json", required=True, help="clustering_dump.cxx output (gp_clustering_figure.py writes it)")
ap.add_argument("--out", required=True)
ap.add_argument("--label", default="VBF H#rightarrowinv.")
args = ap.parse_args()
E = json.load(open(args.json))
C = E["clusters"]
T = {tr["idx"]: tr for c in C for tr in c["tracks"]}          # the clustering input: every timed track

# ---- replay doIterativeClustering ----------------------------------------------------------------------------
left = sorted(T)                                               # track indices not yet consumed
clusters, history = [], []                                     # history[k]: centroid (t, sigma) after each absorption
while left:
    seed = max(left, key=lambda i: T[i]["pt"])
    left.remove(seed)
    t, s, members, hist = T[seed]["t"], T[seed]["res"], [seed], [(T[seed]["t"], T[seed]["res"])]
    while True:
        best, dmin = None, 3.0
        for i in left:
            d = abs(t - T[i]["t"]) / np.hypot(s, T[i]["res"])
            if d < dmin: dmin, best = d, i
        if best is None: break
        left.remove(best)
        w1, w2 = 1 / s**2, 1 / T[best]["res"]**2
        t, s = (t * w1 + T[best]["t"] * w2) / (w1 + w2), 1 / np.sqrt(w1 + w2)
        members.append(best); hist.append((t, s))
    clusters.append(dict(t=t, s=s, members=members, seed=seed)); history.append(hist)
cpp = sorted(sorted(tr["idx"] for tr in c["tracks"]) for c in C)
ok = sorted(sorted(c["members"]) for c in clusters) == cpp
print(f"replay: {len(clusters)} clusters; reproduces the C++ clusters track for track: {ok}")
if not ok: raise SystemExit("the replay does not reproduce the C++ clustering -- not drawing")
# map each replayed cluster to its C++ cluster (for the WAVeS scores and the colours)
cpp_of = {k: next(j for j, c in enumerate(C) if sorted(tr["idx"] for tr in c["tracks"]) == sorted(cl["members"]))
          for k, cl in enumerate(clusters)}
for k, cl in enumerate(clusters):
    print(f"  cluster {k + 1}: seed pT {T[cl['seed']]['pt']:.2f} GeV, {len(cl['members'])} tracks, t {cl['t']:.1f} +- {cl['s']:.1f} ps "
          f"(C++ {C[cpp_of[k]]['t']:.1f})")

# ---- draw: three strips on one time axis --------------------------------------------------------------------
R = S.setup()
PAL = ["#BD1F01", "#3F90DA", "#94A4A2", "#FFA90E", "#832DB6", "#A96B59", "#E76300", "#B9AC70"]   # in seed order
col = [S.colour(PAL[k % len(PAL)]) for k in range(len(clusters))]
def pale(k, a):
    h = PAL[k % len(PAL)].lstrip("#"); rgb = [int(h[j:j + 2], 16) for j in (0, 2, 4)]
    return R.TColor.GetColor(*[int(round(255 - a * (255 - v))) for v in rgb])
GREY = R.TColor.GetColor("#9A9A9A")
ts = [tr["t"] for tr in T.values()]; pts = [tr["pt"] for tr in T.values()]
x0, x1 = min(ts) - 60, max(ts) + 60
y0, y1 = 0.7, max(pts) * 3.2
med_res = float(np.median([tr["res"] for tr in T.values()]))
cv = S.keep(R.TCanvas("steps", "", 1500, 1300))
EDGES = [1.0, 0.655, 0.355, 0.0]
c0 = clusters[0]; sd = T[c0["seed"]]
K = 1                                         # the snapshot: after the FIRST absorption (one track at a time), the next one ringed
nxt = c0["members"][K + 1]; tK, sK = history[0][K]
TAGS = ["1  Seed: the highest-p_{T} track; its 3#sigma window",
        "2  Absorb the closest track within 3#sigma, one at a time; update the mean",
        "3  Next seed, until every track is in a cluster"]

def strip(i):
    cv.cd()
    S.pad(f"s{i}", 0, EDGES[i + 1], 1, EDGES[i], top=0.30 if i == 0 else 0.15, bottom=0.26 if i == 2 else 0.03,
          left=0.10, right=0.03)
    R.gPad.SetLogy()
    h = R.gPad.GetHNDC()
    sz = 0.076 * 0.30 / h                       # the same text size on every strip
    fr = S.keep(R.gPad.DrawFrame(x0, y0, x1, y1))
    S.axes(fr, "Track time [ps]" if i == 2 else "", "p_{T} [GeV]", size=sz, ytoff=0.75 * h / 0.30, xtoff=1.05)
    fr.GetYaxis().SetMoreLogLabels(); fr.GetYaxis().SetNoExponent(); fr.GetXaxis().SetNdivisions(508)
    if i != 2: fr.GetXaxis().SetLabelSize(0)
    # the step's label in the strip's top margin, above the frame: nothing drawn later can cover it
    tag = S.keep(R.TLatex(0.10, 1 - R.gPad.GetTopMargin() + 0.25 * sz, TAGS[i]))
    tag.SetNDC(); tag.SetTextFont(62); tag.SetTextSize(sz); tag.SetTextAlign(11); tag.Draw()
    return sz

def tracks(colour_of):
    for i, tr in T.items():
        c, filled = colour_of(i)
        g = S.keep(R.TGraphErrors(1)); g.SetPoint(0, tr["t"], tr["pt"]); g.SetPointError(0, tr["res"], 0)
        g.SetMarkerStyle(20 if filled else 24); g.SetMarkerSize(1.6); g.SetMarkerColor(c); g.SetLineColor(c); g.SetLineWidth(2); g.Draw("P SAME")

def truth_line(label=False):   # truth t_HS: above the bands, under the tracks
    l = S.keep(R.TLine(E["t_hs"], y0, E["t_hs"], y1)); l.SetLineStyle(2); l.SetLineWidth(3); l.SetLineColor(R.kBlack); l.Draw()
    if label:
        t = S.keep(R.TLatex(E["t_hs"] + 6, y1 / 1.5, "truth t_{HS}")); t.SetTextAlign(13); t.SetTextSize(0.85 * sz); t.Draw()

def band(t, s, colour):
    b = S.keep(R.TBox(max(x0, t - s), y0, min(x1, t + s), y1)); b.SetFillColor(colour); b.SetLineWidth(0); b.Draw()

# 1: the seed and its 3 sigma window
sz = strip(0)
band(sd["t"], 3 * np.hypot(sd["res"], med_res), pale(0, 0.20))
truth_line(label=True)
tracks(lambda i: (col[0], True) if i == c0["seed"] else (GREY, False))
st = S.keep(R.TMarker(sd["t"], sd["pt"], 29)); st.SetMarkerSize(3.6); st.SetMarkerColor(col[0]); st.Draw()
S.header(0.10, 0.90, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}: forward tracks with a time"], size=0.85 * sz, step=1.15 * sz)
R.gPad.RedrawAxis()
# 2: a snapshot of the growth: K tracks absorbed, the mean, the next one within 3 sigma ringed
sz = strip(1)
band(tK, 3 * np.hypot(sK, T[nxt]["res"]), pale(0, 0.14))
band(tK, sK, pale(0, 0.55))
truth_line()
absorbed = set(c0["members"][:K + 1])
tracks(lambda i: (col[0], True) if i in absorbed else (GREY, False))
st2 = S.keep(R.TMarker(sd["t"], sd["pt"], 29)); st2.SetMarkerSize(3.6); st2.SetMarkerColor(col[0]); st2.Draw()
ring = S.keep(R.TMarker(T[nxt]["t"], T[nxt]["pt"], 24)); ring.SetMarkerSize(3.8); ring.SetMarkerColor(R.kBlack); ring.SetMarkerStyle(24); ring.Draw()
ar = S.keep(R.TLatex(T[nxt]["t"] - T[nxt]["res"] - 10, T[nxt]["pt"] / 1.22, "next"))   # below left of the ring
ar.SetTextSize(0.85 * sz); ar.SetTextAlign(33); ar.Draw()
note = S.keep(R.TLatex(0.95, 1 - R.gPad.GetTopMargin() - 0.06, f"seed + 1 track: mean {tK:.0f} #pm {sK:.0f} ps"))
note.SetNDC(); note.SetTextAlign(33); note.SetTextSize(0.8 * sz); note.SetTextColor(R.kGray + 2); note.Draw()
note2 = S.keep(R.TLatex(0.95, 1 - R.gPad.GetTopMargin() - 0.06 - 1.15 * 0.8 * sz, "dark: mean #pm #sigma;  light: its 3#sigma window"))
note2.SetNDC(); note2.SetTextAlign(33); note2.SetTextSize(0.8 * sz); note2.SetTextColor(R.kGray + 2); note2.Draw()
R.gPad.RedrawAxis()
# 3: all clusters
sz = strip(2)
member_of = {i: k for k, cl in enumerate(clusters) for i in cl["members"]}
for k, cl in enumerate(clusters): band(cl["t"], cl["s"], pale(k, 0.35))
truth_line()
tracks(lambda i: (col[member_of[i]], True))
R.gPad.RedrawAxis()
S.save(cv, args.out)
