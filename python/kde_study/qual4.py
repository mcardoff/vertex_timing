"""Purity-targeted quality: which reco variable best selects events whose PICKED CLUSTER is pure?
Mean cluster purity (and core) of the best 70% / 50% of events by each candidate."""
from qual import *

Q = {}
for n, m in SAMPLES:
    d = get(n, m); q = qvars(d)
    S2 = q['S1'] - q['S1mS2']
    q['Q'] = q['S1mS2'] / np.sqrt(q['S1'] + S2)
    Q[n] = q
C = {
    'Q': lambda q: q['Q'],
    'Q / krms': lambda q: q['Q'] / np.maximum(q['krms'], 10),
    'Q / krms^2': lambda q: q['Q'] / np.maximum(q['krms'], 10) ** 2,
    'Q / (1+tchi2)': lambda q: q['Q'] / (1 + q['tchi2']),
    'Q / (1+tchi2)^0.5': lambda q: q['Q'] / (1 + q['tchi2']) ** 0.5,
    'Q * frac': lambda q: q['Q'] * q['frac'],
    'Q * frac^0.5': lambda q: q['Q'] * q['frac'] ** 0.5,
    'Q * nfrac^0.5': lambda q: q['Q'] * q['nfrac'] ** 0.5,
    'Q * (wz/neff)': lambda q: q['Q'] * q['wz'] / q['neff'],
    'Q * (wz/neff)^2': lambda q: q['Q'] * (q['wz'] / q['neff']) ** 2,
    'Q * e^-absdz': lambda q: q['Q'] * np.exp(-q['absdz']),
    'Q * (wz/neff) / krms': lambda q: q['Q'] * q['wz'] / q['neff'] / np.maximum(q['krms'], 10),
    'wz/neff alone': lambda q: q['wz'] / q['neff'],
    'krms alone (low)': lambda q: -q['krms'],
    'TRUTH: cluster purity itself': lambda q: np.nan_to_num(q['purity_cluster']),
}
print(f"{'candidate':28s} | " + " | ".join(f"{n.split('_')[0]:^31s}" for n, _ in SAMPLES))
print(f"{'':28s} | " + " | ".join(" purity@70/50     core@70/50   " for _ in SAMPLES))
for name, fn in C.items():
    row = []
    for n, _ in SAMPLES:
        q = Q[n]; x = np.nan_to_num(fn(q), nan=-1e9); pur = np.nan_to_num(q['purity_cluster'])
        r = []
        for f in (.7, .5):
            s = x >= np.quantile(x, 1 - f); r.append((pur[s].mean() * 100, q['ok'][s].mean() * 100))
        row.append(f"{r[0][0]:5.1f}/{r[1][0]:5.1f}    {r[0][1]:5.1f}/{r[1][1]:5.1f}  ")
    print(f"{name:28s} | " + " | ".join(row))
