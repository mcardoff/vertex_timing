"""Quality score candidates: combinations. Reports, for the best X% of events by the candidate,
core fraction (time purity) and mean cluster purity (pT-weighted HS fraction of the picked cluster)."""
from qual import *

Q = {}
for n, m in SAMPLES:
    d = get(n, m); Q[n] = qvars(d)


def S2(q): return q['S1'] - q['S1mS2']


CANDS = {
    'S1': lambda q: q['S1'],
    'S1 - S2': lambda q: q['S1mS2'],
    'S1 - 0.5 S2': lambda q: q['S1'] - 0.5 * S2(q),
    'S1 - 1.5 S2': lambda q: q['S1'] - 1.5 * S2(q),
    'S1 - 2 S2': lambda q: q['S1'] - 2 * S2(q),
    '(S1-S2)/sqrt(S1)': lambda q: q['S1mS2'] / np.sqrt(q['S1']),
    '(S1-S2)/sqrt(S1+S2)': lambda q: q['S1mS2'] / np.sqrt(q['S1'] + S2(q)),
    'S1 * margin^2': lambda q: q['S1'] * q['margin'] ** 2,
    '(S1-S2) / krms': lambda q: q['S1mS2'] / np.maximum(q['krms'], 10),
    '(S1-S2) * e^-absdz': lambda q: q['S1mS2'] * np.exp(-q['absdz']),
    '(S1-S2), x0.3 if in/out > 60ps': lambda q: q['S1mS2'] * np.where(np.nan_to_num(q['inout'], nan=0) > 60, 0.3, 1.0),
    '(S1-S2), x0.5 if in/out > 60ps': lambda q: q['S1mS2'] * np.where(np.nan_to_num(q['inout'], nan=0) > 60, 0.5, 1.0),
    '(S1-S2) * (1 + n_injet_k)^0.5': lambda q: q['S1mS2'] * (1 + q['n_injet_k']) ** 0.5,
    'dens - dens2-like: frac*S1': lambda q: q['frac'] * q['S1'],
    '(S1-S2) * frac': lambda q: q['S1mS2'] * q['frac'],
    '(S1-S2) * frac^0.5': lambda q: q['S1mS2'] * q['frac'] ** 0.5,
}
print(f"{'candidate':32s} | " + " | ".join(f"{n.split('_')[0]:^35s}" for n, _ in SAMPLES))
print(f"{'':32s} | " + " | ".join("AUC   core@90/70/50   pur@90/70/50  " for _ in SAMPLES))
for name, fn in CANDS.items():
    row = []
    for n, _ in SAMPLES:
        q = Q[n]; x = np.nan_to_num(fn(q), nan=-1e9); ok = q['ok']; pur = np.nan_to_num(q['purity_cluster'])
        a = auc(x, ok)
        cs = []; ps = []
        for f in (.9, .7, .5):
            s = x >= np.quantile(x, 1 - f); cs.append(ok[s].mean() * 100); ps.append(pur[s].mean() * 100)
        row.append(f"{a:.3f} {cs[0]:4.1f}/{cs[1]:4.1f}/{cs[2]:4.1f}  {ps[0]:4.1f}/{ps[1]:4.1f}/{ps[2]:4.1f}")
    print(f"{name:32s} | " + " | ".join(row))
print("\nno cut: " + " | ".join(f"{n.split('_')[0]} core {Q[n]['ok'].mean()*100:.1f} pur {np.nanmean(Q[n]['purity_cluster'])*100:.1f}" for n, _ in SAMPLES))
