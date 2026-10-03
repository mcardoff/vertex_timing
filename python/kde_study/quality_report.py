"""Before/after the time-quality cut, from the C++ outputs.
  quality_report.py <label> <clustering_hist.root> <rpt_v5_hist.root>
Core fraction + mean cluster purity per row, and pileup rejection at MATCHED hard-scatter efficiency
(interpolated on the traced (eff, rej) curve, not the nearest-bin cells of rpt_v5_plot's console table)."""
import sys, numpy as np, uproot

lab, fc, fr = sys.argv[1:4]
c = uproot.open(fc); r = uproot.open(fr)
ROWS = [("HGTD (Athena)", "HGTD Algorithm", "hgtd"), ("TZP", "d0 precision]", "trkptztzp"),
        ("TZP_KDE_TZ (no cut)", "[(t,z) clusters]", "tzpkdetz"), ("TZP_KDE_TZ, Q >= 2", "[(t,z), Q #geq 2]", "tzpkdetzq")]
names = [k.split(';')[0] for k in c.keys()]
key = "n_Forward_Jets"
print(f"== {lab}: clustering ==")
print(f"{'row':28s} {'time provided':>14s} {'core of provided':>17s} {'core over all':>14s} {'mean cluster purity':>20s}")
ntot = None
for name, pat, pk in ROWS:
    n = [x for x in names if x.startswith("flat_" + key + "_") and pat in x][0]
    tot = c[n].values(flow=True).sum(); good = c["good_" + n[5:]].values(flow=True).sum()
    if ntot is None: ntot = tot
    h = c["purity_hgtdtimes_" + pk]; v = h.values(flow=True)[1:-1]; x = h.axis().centers()
    pur = (v * x).sum() / v.sum() if v.sum() else float('nan')
    if False:
        pass
    else:
        # every row shares the denominator (all selected events); the purity histogram is filled
        # only where a time is provided, so its entries give the acceptance
        nprov = h.values(flow=True).sum()
        print(f"{name:28s} {100*nprov/ntot:13.1f}% {100*good/nprov:16.2f}% {100*good/tot:13.2f}% {100*pur:19.1f}%")


def curve(hs, pu):
    """(eff, rej) for cuts R_pT > edge_k, all edges"""
    a = hs.values(flow=True); b = pu.values(flow=True)
    eff = (a[::-1].cumsum()[::-1] / a.sum())[1:]       # cut at lower edge of bin k (k>=1): excludes underflow
    mis = (b[::-1].cumsum()[::-1] / b.sum())[1:]
    return eff, mis


def rej_at(eff, mis, target):
    # eff decreases with the cut; find bracketing points
    i = np.where(eff >= target)[0]
    if len(i) == 0 or i[-1] == len(eff) - 1: return float('nan')
    k = i[-1]
    e0, e1 = eff[k], eff[k + 1]; m0, m1 = mis[k], mis[k + 1]
    if e0 == e1: return 1 / m0 if m0 > 0 else float('inf')
    f = (e0 - target) / (e0 - e1)
    m = m0 + f * (m1 - m0)
    return 1 / m if m > 0 else float('inf')


SC = ["zonly", "hgtd", "trkptz", "waves", "tzp", "kde", "kde_q1", "kde_q"]
SL = [("fwd >40 GeV", "_hi", "_hi"), ("fwd 30-40 GeV", "_lo", "_lo"), ("R1 (both tags fwd)", "_r1", "_r1"), ("R2 (fwd PU, HS from R1)", "_r1", "_r2")]
nk = r["meta_n_kde_time"].member("fVal"); nq = r["meta_n_kde_quality"].member("fVal")
print(f"\n== {lab}: R_pT, pileup rejection at matched HS efficiency ==   (quality flag passes in {100*nq/nk:.1f}% of events with a time, RpT selection)")
for title, shs, spu in SL:
    print(f"  {title}:   HS jets {r['HS_zonly'+shs].values(flow=True).sum():.0f}, PU jets {r['PU_zonly'+spu].values(flow=True).sum():.0f}")
    print(f"    {'scenario':10s} {'eff ceiling':>11s} {'rej@0.80':>9s} {'rej@0.85':>9s} {'rej@0.90':>9s}  {'<RpT> HS':>9s} {'<RpT> PU':>9s} {'PU at RpT=0':>12s}")
    for s in SC:
        hs = r["HS_" + s + shs]; pu = r["PU_" + s + spu]
        eff, mis = curve(hs, pu)
        a = hs.values(flow=True); b = pu.values(flow=True); x = hs.axis().centers()
        mh = (a[1:-1] * x).sum() / a[1:-1].sum(); mp = (b[1:-1] * x).sum() / b[1:-1].sum()
        print(f"    {s:10s} {100*eff[1]:10.1f}% " + " ".join(f"{rej_at(eff, mis, t):9.1f}" for t in (.80, .85, .90)) + f"  {mh:9.3f} {mp:9.4f} {100*b[1]/b.sum():11.1f}%")
