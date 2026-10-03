"""Ceiling test: how much could ANY clustering change buy?  Truth-assisted partitions, same scores."""
from reclus import *
sw = lambda d: d['tracks']['w_tzp']
print("% of production-TZP fails removed; columns zjets dijet vbf ttbar")
rows = {}
for n, m in SAMPLES:
    d = get(n, m); t = d['tracks']; hs = t['truth_is_hs'] > 0.5
    f0 = 1 - d['ref']['tzp'].mean(); b0 = d['ref']['tzp'].mean()
    tt = d['ttruth'][t['ev']]
    intime = hs & (np.abs(t['time_d'] - tt) < 3 * t['timeRes_d'])
    base = cluster(d, cut=3.0, useZ=True, seedw=sw(d))
    for lab_name, sel in (("A: all truth-HS tracks in one cluster", hs), ("B: in-time (3 sigma) truth-HS tracks in one cluster", intime)):
        lab = np.where(sel, 4000, base).astype(np.int32)
        r = evaluate(d, lab)
        rows.setdefault(lab_name, []).append(((r['tzp'].mean() - b0) / f0 * 100, (r['kde'].mean() - b0) / f0 * 100))
    # C: perfect partition AND perfect selection of that cluster (time = its plain mean) -> pure timing ceiling
    c, crow = features(d, np.where(intime, 4000, base).astype(np.int32))
    lab = np.where(intime, 4000, base).astype(np.int32)
    key = t['ev'].astype(np.int64) * 4096 + lab; uk = np.unique(key)
    is_hs = (uk % 4096) == 4000
    okc = np.zeros(d['nev'], bool); okc[c['ev'][is_hs]] = np.abs(c['time'][is_hs] - d['ttruth'][c['ev'][is_hs]]) < 60
    rows.setdefault("C: B + oracle selection of the HS cluster", []).append(((okc.mean() - b0) / f0 * 100, np.nan))
for k, v in rows.items():
    print(f"{k:52s} TZP: " + " ".join(f"{a:+6.2f}" for a, _ in v) + "  | KDE: " + " ".join(f"{b:+6.2f}" for _, b in v))
run("(t,z) 3.0 seed w (reco, for reference)", cut=3.0, useZ=True, seedw=sw)
