"""Quality-flag working points on the reliability-weighted configuration."""
from sel import *
from sel8 import zf, rel
from qual import athena
print("provided / right among provided   (zjets | dijet | vbf | ttbar)")
R = {}
for n, m in SAMPLES:
    d = get(n, m); c = cands(d, sfac=rel(d), tfac=zf(d)); p1, p2, Q = top2(d, c); R[n] = (Q, c['ok'][p1], athena(d))
for thr in (1.0, 1.25, 1.5, 1.75, 2.0):
    print(f"Q >= {thr:4.2f}: " + " | ".join(f"{np.mean(R[n][0]>=thr)*100:5.1f} / {R[n][1][R[n][0]>=thr].mean()*100:5.1f}" for n, _ in SAMPLES))
print("Athena   : " + " | ".join(f"{R[n][2][0].mean()*100:5.1f} / {R[n][2][1][R[n][2][0]].mean()*100:5.1f}" for n, _ in SAMPLES))
print("ours at Athena's acceptance: " + " | ".join(f"{R[n][1][R[n][0] >= np.quantile(R[n][0], 1 - R[n][2][0].mean())].mean()*100:5.1f}" for n, _ in SAMPLES))
edges = [-9, 0.35, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 1e9]
for i in range(len(edges) - 1):
    print(f"  Q {edges[i]:5.2f}-{min(edges[i+1],99):5.2f}: " + " | ".join(f"{R[n][1][(R[n][0]>=edges[i])&(R[n][0]<edges[i+1])].mean()*100:5.1f}" for n, _ in SAMPLES))
