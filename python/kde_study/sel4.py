"""Time reliability by number of HGTD hits: mechanism, then a scan of the per-hit-count factor in the selection density."""
from sel3 import *

print("truth-HS tracks: fraction whose HGTD time is within 3 sigma (and within 60 ps) of the truth HS time, by n HGTD hits")
for n, _ in SAMPLES:
    d, c = D[n]; t = d['tracks']; hs = t['truth_is_hs'] > 0.5
    dtt = np.abs(t['time_d'] - d['ttruth'][t['ev']])
    row = []
    for k in (1, 2, 3, 4):
        s = hs & (t['nhgtd_hits'] == k)
        a = (t['nhgtd_hits'] == k)
        row.append(f"n={k}: {np.mean(dtt[s] < 3*t['timeRes_d'][s])*100:4.1f}% / {np.mean(dtt[s] < 60)*100:4.1f}% (share of tracks {a.mean()*100:4.1f}%, sigma_t {np.median(t['timeRes_d'][a]):.0f} ps)")
    print(f"  {n:12s} " + "  ".join(row))

print("\n% of TZP fails removed:                          zjets  dijet    vbf  ttbar")
run("reference", lambda d, c: c['S'])


def fac(d, f):
    return np.asarray(f)[np.clip(d['tracks']['nhgtd_hits'].astype(int), 1, 4) - 1]


for f in ((1, 2, 2, 2), (1, 3, 3, 3), (1, 4, 4, 4), (1, 2, 3, 3), (1, 2, 3, 4), (1, 2.5, 4, 4), (1, 3, 5, 5), (1, 2, 2.5, 2.5), (1, 1.5, 2.5, 2.5), (1, 3, 3, 2), (1, 2, 4, 2)):
    run(f"factors by n hits {f}", lambda d, c, f=f: dens_with(d, c, fac(d, f)) * c['env'])
for k in (3.0, 4.0):
    run(f"(25/sigma_t)^{k}", lambda d, c, k=k: dens_with(d, c, (25. / d['tracks']['timeRes_d']) ** k) * c['env'])
