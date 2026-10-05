"""Tie-break scan repeated on the reliability-weighted baseline."""
from sel2 import *
from sel8 import zf, rel
res = {}
for n, m in SAMPLES:
    d = get(n, m); c = cands(d, sfac=rel(d), tfac=zf(d)); p1, p2, Q = top2(d, c); v = cvars(d, c)
    has = p2 >= 0; a = p1[has]; b = p2[has]; q = Q[has]
    dec = c['ok'][a] ^ c['ok'][b]; first = c['ok'][a]; lowq = q < 0.7
    print(f"{n}: pick {c['ok'][p1].mean()*100:.2f}  pick|2nd {np.mean(c['ok'][p1] | np.where(p2>=0, c['ok'][np.maximum(p2,0)], False))*100:.2f}  score accuracy on decisive: all {first[dec].mean()*100:.1f}, Q<0.7 {first[dec & lowq].mean()*100:.1f} ({np.mean(dec&lowq)*100:.1f}% of ev)")
    for k, x in v.items():
        d1 = x[a] > x[b]; tie = x[a] == x[b]
        acc = lambda s: (np.mean(np.where(tie[s], 0.5, d1[s] == first[s])) * 100)
        res.setdefault(k, []).append((acc(dec), acc(dec & lowq)))
print("\naccuracy of 'take the larger X' on decisive pairs:  all / low Q      (zjets | dijet | vbf | ttbar)")
for k in sorted(res, key=lambda k: -min(max(r[1], 100 - r[1]) for r in res[k]))[:14]:
    print(f"{k:14s} " + " | ".join(f"{a:5.1f} / {b:5.1f}" for a, b in res[k]))
