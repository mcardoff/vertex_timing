"""Track-list scan 4: how to admit 0.5-1 GeV tracks, one rule for every sample.
   wide4.py <trackdump.root> [max_file_idx]      (all rows: no pT ceiling, >= 1 GeV tracks at z0 significance < 4)"""
import sys
from wide import *

path = sys.argv[1]; mf = int(sys.argv[2]) if len(sys.argv) > 2 else None
ev, tr = load(path, mf); nev = len(ev['vz'])
base = evaluate(ev, tr, nominal(tr))
hi = nominal(tr, 4.0, 1.0, 1e9)
ref = evaluate(ev, tr, hi)
print(f"== {path}: {nev} events; nominal {base.mean()*100:.2f}; no ceiling + signif<4: {ref.mean()*100:.2f} ({100*(ref.mean()-base.mean())/(1-base.mean()):+.2f}% of fails)")
low = lambda f: (lambda t: np.where(t['pt_d'] < 1.0, f, 1.0))
two = tr['nhgtd_hits'] >= 2
out = []
for ns in (3.0, 4.0):
    for hits, hl in ((np.ones(len(two), bool), "all"), (two, ">=2 hits")):
        for f in (1.0, 0.7, 0.5):
            m = hi | (nominal(tr, ns, 0.5, 1.0) & hits)
            ok = evaluate(ev, tr, m, extra=low(f))
            print(f"  sub-GeV: signif < {ns}, {hl:8s}, x{f}:  core {ok.mean()*100:6.2f}  fails removed vs nominal {100*(ok.mean()-base.mean())/(1-base.mean()):+6.2f}%  tracks/ev {m.sum()/nev:5.1f}")
for lo in (0.6, 0.7, 0.8):
    m = hi | nominal(tr, 3.0, lo, 1.0)
    ok = evaluate(ev, tr, m)
    print(f"  sub-GeV: pT > {lo}, signif < 3, all, x1:  core {ok.mean()*100:6.2f}  fails removed vs nominal {100*(ok.mean()-base.mean())/(1-base.mean()):+6.2f}%  tracks/ev {m.sum()/nev:5.1f}")
