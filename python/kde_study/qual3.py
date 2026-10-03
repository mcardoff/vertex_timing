"""Is the quality score calibrated across samples?  P(time correct) in bins of Q, and working points at fixed thresholds."""
from qual import *

SETS = SAMPLES + [('vbf_mu0_novbs', 20)]
Q = {}
for n, m in SETS:
    d = get(n, m); q = qvars(d)
    S2 = q['S1'] - q['S1mS2']
    q['Q'] = q['S1mS2'] / np.sqrt(q['S1'] + S2)
    q['tzp_ok'], q['tzp_pur'] = baseline_tzp(d)
    q['ath_valid'], q['ath_ok'] = athena(d)
    Q[n] = q
edges = [-1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 1e9]
print("P(|dt| < 60 ps) in bins of Q = (S1 - S2)/sqrt(S1 + S2)   [share of events]")
print(f"{'Q bin':12s} " + " ".join(f"{n.split('_nov')[0]:>18s}" for n, _ in SETS))
for i in range(len(edges) - 1):
    row = []
    for n, _ in SETS:
        q = Q[n]; s = (q['Q'] >= edges[i]) & (q['Q'] < edges[i + 1])
        row.append(f"{q['ok'][s].mean()*100:6.1f} [{s.mean()*100:5.1f}%]" if s.sum() > 30 else "        -        ")
    print(f"{edges[i]:5.2f}-{min(edges[i+1],99):5.2f} " + " ".join(f"{r:>18s}" for r in row))
print("\nWorking points: events kept / core fraction of kept / mean cluster purity of kept / core counted over ALL events (rejected = fail)")
for thr in (0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
    row = []
    for n, _ in SETS:
        q = Q[n]; s = q['Q'] >= thr
        row.append(f"{s.mean()*100:5.1f}/{q['ok'][s].mean()*100:5.1f}/{np.nanmean(q['purity_cluster'][s])*100:5.1f}/{(q['ok'] & s).mean()*100:5.1f}")
    print(f"Q >= {thr:4.2f}: " + "   ".join(row))
print("\nReferences (kept / core of kept / cluster purity):")
for n, _ in SETS:
    q = Q[n]; v = q['ath_valid']
    print(f"  {n:18s} TZP all {q['tzp_ok'].mean()*100:5.1f} pur {np.nanmean(q['tzp_pur'])*100:5.1f} | new all {q['ok'].mean()*100:5.1f} pur {np.nanmean(q['purity_cluster'])*100:5.1f} | "
          f"Athena provides {v.mean()*100:5.1f}%, core of those {q['ath_ok'][v].mean()*100:5.1f} | new at Athena's acceptance: core "
          f"{q['ok'][q['Q'] >= np.quantile(q['Q'], 1 - v.mean())].mean()*100:5.1f}")
