"""ATLAS-style drawing for the TDR-structured figures (tdr_pujet_plots.py, tdr_vbf_plots.py).

ROOT with the repo's src/AtlasStyle.h and src/AtlasLabels.h, as the other ROOT figures here
(vbs_time_veto_common.setup_root). A figure with several panels is one canvas cut into pads of the
size a single 800 x 600 (or 800 x 800 with a ratio pad) ATLAS plot would have, so the text sizes,
which ROOT takes relative to the pad height, are the standard ones in every panel.
"""
import ctypes
import os
import subprocess

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = None
_keep = []   # PyROOT objects drawn on pads, kept alive until the canvas is saved


def setup():
    global ROOT
    import ROOT as R
    ROOT = R
    R.gROOT.SetBatch(True)
    R.gInterpreter.Declare(f'#include "{REPO}/src/AtlasStyle.h"\n#include "{REPO}/src/AtlasLabels.h"')
    R.SetAtlasStyle()
    R.gStyle.SetOptStat(0)
    R.gStyle.SetEndErrorSize(3)
    R.gStyle.SetHatchesLineWidth(2)
    R.gErrorIgnoreLevel = R.kWarning
    return R


def keep(obj):
    _keep.append(obj)
    return obj


def colour(hexcode):
    return ROOT.TColor.GetColor(hexcode)


def pad(name, x1, y1, x2, y2, top=0.05, bottom=0.16, left=0.16, right=0.05):
    """a pad of the current canvas, drawn and made current"""
    p = keep(ROOT.TPad(name, "", x1, y1, x2, y2))
    p.SetTopMargin(top); p.SetBottomMargin(bottom); p.SetLeftMargin(left); p.SetRightMargin(right)
    p.SetFillStyle(4000)
    p.Draw(); p.cd()
    return p


def text(x, y, s, size, font=42, align=11):
    t = ROOT.TLatex()
    t.SetNDC(); t.SetTextFont(font); t.SetTextSize(size); t.SetTextAlign(align)
    return keep(t.DrawLatex(x, y, s))


def header(x, y, lines=(), size=0.05, step=None):
    """'ATLAS Simulation Internal' at NDC (x, y) of the current pad, one line below it per entry of
    `lines`. The gap after 'ATLAS' is measured from the drawn word, so it is right on any pad."""
    step = step or 1.2 * size
    a = text(x, y, "ATLAS", size, font=72)
    w, h = ctypes.c_uint(0), ctypes.c_uint(0)
    a.GetBoundingBox(w, h)
    pw = ROOT.gPad.GetWw() * ROOT.gPad.GetAbsWNDC()
    dx = w.value / pw if w.value > 0 else 0.115 * 696 * ROOT.gPad.GetWh() / (472 * ROOT.gPad.GetWw())
    text(x + dx + 0.015, y, "Simulation Internal", size)
    for i, s in enumerate(lines):
        text(x, y - (i + 1) * step, s, 0.84 * size)


def legend(x1, y1, x2, y2, size=0.04, ncol=1, margin=0.2):
    l = keep(ROOT.TLegend(x1, y1, x2, y2))
    ROOT.StyleLegend(l, size)
    l.SetNColumns(ncol)
    l.SetMargin(margin)
    return l


def graph(x, y, ey=None):
    """a TGraph (or TGraphErrors with y errors) of the finite points"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(y) & np.isfinite(x)
    if ey is not None:
        ey = np.asarray(ey, float); m &= np.isfinite(ey)
    xs, ys = np.ascontiguousarray(x[m]), np.ascontiguousarray(y[m])
    if ey is None:
        return keep(ROOT.TGraph(int(m.sum()), xs, ys))
    return keep(ROOT.TGraphErrors(int(m.sum()), xs, ys, np.zeros(int(m.sum())), np.ascontiguousarray(ey[m])))


def style(obj, hexcode, marker=0, line=1, width=2, msize=1.3):
    c = colour(hexcode)
    obj.SetLineColor(c); obj.SetLineStyle(line); obj.SetLineWidth(width)
    obj.SetMarkerColor(c); obj.SetMarkerStyle(marker if marker else 1); obj.SetMarkerSize(msize if marker else 0)
    return obj


def axes(frame, xtitle, ytitle, size=0.05, ytoff=1.4, xtoff=1.4, ndiv_y=None):
    for ax, t, off in ((frame.GetXaxis(), xtitle, xtoff), (frame.GetYaxis(), ytitle, ytoff)):
        ax.SetTitle(t); ax.SetTitleSize(size); ax.SetLabelSize(size); ax.SetTitleOffset(off)
    if ndiv_y: frame.GetYaxis().SetNdivisions(ndiv_y)
    return frame


def save(canvas, stem, scale=1.5):
    """the PDF on a page of the canvas's own shape, via ROOT's EPS (whose bounding box is the canvas)
    and epstopdf -- ROOT 6.40's PDF writer puts every canvas on a fixed ~A4 page whatever the style's
    paper size, clipping a wide one -- and the PNG rendered from that PDF, because ROOT's own PNG
    writer restarts the dash pattern on every segment and breaks dashed curves of many segments"""
    os.makedirs(os.path.dirname(os.path.abspath(stem)), exist_ok=True)
    W = canvas.GetWw()
    try:
        canvas.SaveAs(f"{stem}.eps")
        subprocess.run(["epstopdf", f"{stem}.eps", f"--outfile={stem}.pdf"], check=True)
        os.remove(f"{stem}.eps")
        subprocess.run(["pdftoppm", "-png", "-singlefile", "-scale-to-x", str(int(scale * W)), "-scale-to-y", "-1",
                        f"{stem}.pdf", stem], check=True)
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"[tdr_root_style] {e}: falling back to ROOT's own PDF and PNG")
        canvas.SaveAs(f"{stem}.pdf"); canvas.SaveAs(f"{stem}.png")
    print("wrote", stem + ".png")
