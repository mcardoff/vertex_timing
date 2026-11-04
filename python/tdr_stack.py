"""The m_jj stack: the selected events' fractions per band in m_jj bins (overflow folded into the last bin),
with the events per bin in a lower panel. Shared by tdr_vbf_plots.py and tdr_fwdpair_kin.py so the stacks look
the same (moved out of tdr_vbf_plots.py unchanged, 2026-11-04)."""
import numpy as np

import tdr_root_style as S

MAIN, RAT = 0.68, 0.32
TS_M, TS_R = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)
# truth-identity colours: python/vbs_region_stack.py's
IDENT = [("HS + HS (genuine pair)", "#7BCB6B"), ("HS + PU", "#F5A02A"), ("PU + PU", "#636363"), ("A leg carries neither label", "#C79BF5")]


def mjj_stack(R_, stem, cname, weight, ok, mjj, bands, EDGES, header, label):
    """bands: (legend, mask, colour), bottom to top; the masks must partition the selected events.
    stem: the output path without extension; cname: a unique ROOT name; header: the lines under 'ATLAS'"""
    w = np.asarray(weight).astype(float)
    name, key = cname, ""
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
    S.header(0.20, 0.87, header, size=TS_M)
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
    print(f"  {label}, per m_jj bin (" + " / ".join(b[0].replace('#', '') for b in bands) + "): " +
          "  ".join(f"[{EDGES[b]:.0f},{EDGES[b+1]:.0f}) " + "/".join(f"{100 * w[ok & m & (xm >= EDGES[b]) & (xm < EDGES[b + 1])].sum() / tot[b]:.0f}" for _, m, _ in bands)
                    for b in range(len(EDGES) - 1) if tot[b] > 0))
    S.save(cv, stem)
