"""What else predicts that a track's HGTD time is RIGHT (truth-HS tracks, within 3 sigma of truth)?"""
from sel3 import *

for n in ('zjets_novbs', 'vbf_novbs'):
    d, c = D[n]; t = d['tracks']; hs = t['truth_is_hs'] > 0.5
    good = np.abs(t['time_d'] - d['ttruth'][t['ev']]) < 3 * t['timeRes_d']
    one = t['nhgtd_hits'] == 1
    ae = np.abs(t['eta'].astype(float)); pt = t['pt_d']
    print(f"\n{n}: P(time right) for truth-HS tracks, 1 hit / >=2 hits  [share of HS tracks]")
    def tab(name, x, edges):
        row = []
        for i in range(len(edges) - 1):
            s = hs & (x >= edges[i]) & (x < edges[i + 1])
            a = s & one; b = s & ~one
            row.append(f"{edges[i]:g}-{edges[i+1]:g}: {good[a].mean()*100:4.1f}/{good[b].mean()*100:4.1f} [{s.sum()/hs.sum()*100:4.1f}%, 1-hit share {a.sum()/max(s.sum(),1)*100:3.0f}%]")
        print(f"  {name:8s} " + "  ".join(row))
    tab("|eta|", ae, [2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0])
    tab("pT", pt, [1, 1.5, 2, 3, 5, 10, 30])
    tab("sigma_z0", t['sigma_z0_d'], [0, 0.3, 0.6, 1.2, 2.4, 99])
    tab("in jet", t['injet'].astype(float), [0, 0.5, 1.5])
    ntr = np.bincount(t['ev'], minlength=d['nev'])[t['ev']].astype(float)
    tab("n trk ev", ntr, [0, 15, 25, 35, 50, 999])
