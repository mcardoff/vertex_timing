"""Track-list scan 2: pT ceiling, z0 significance, sub-GeV tracks with a down-weight."""
import sys
from wide import *
path = sys.argv[1]; mf = int(sys.argv[2]) if len(sys.argv) > 2 else None
ev, tr = load(path, mf); nev = len(ev['vz'])
q = tr['quality'] > 0.5
base = evaluate(ev, tr, nominal(tr))
print(f"{path}: {nev} events; nominal list core {base.mean()*100:.2f}%")
def run(label, mask, **kw):
    ok = evaluate(ev, tr, mask, **kw)
    print(f"{label:58s} core {ok.mean()*100:6.2f}  fails removed vs nominal list {100*(ok.mean()-base.mean())/(1-base.mean()):+6.2f}%   tracks/ev {mask.sum()/nev:5.1f}")
    return ok
run("no pT ceiling", nominal(tr, 3.0, 1.0, 1e9))
for cap in (50., 100.):
    run(f"pT ceiling {cap:.0f}", nominal(tr, 3.0, 1.0, cap))
for ns in (3.5, 4.0, 5.0, 8.0):
    run(f"no ceiling, signif < {ns}", nominal(tr, ns, 1.0, 1e9))
low = lambda f, thr=1.0: (lambda t: np.where(t['pt_d'] < thr, f, 1.0))
for lo in (0.75, 0.5):
    for f in (1.0, 0.5, 0.25):
        run(f"no ceiling, signif < 4, pT > {lo}, sub-GeV x{f}", nominal(tr, 4.0, lo, 1e9), extra=low(f))
for pw in (1.0, 2.0):
    run(f"no ceiling, signif < 4, pT > 0.5, weight x min(1,pT)^{pw}", nominal(tr, 4.0, 0.5, 1e9), extra=lambda t, pw=pw: np.minimum(1.0, t['pt_d']) ** pw)
run("no ceiling, signif < 4, pT > 0.5, sub-GeV need >= 2 hits", nominal(tr, 4.0, 0.5, 1e9) & ((tr['pt'] >= 1) | (tr['nhgtd_hits'] >= 2)))
run("no ceiling, signif < 4, sub-GeV: >= 2 hits and x0.5", nominal(tr, 4.0, 0.5, 1e9) & ((tr['pt'] >= 1) | (tr['nhgtd_hits'] >= 2)), extra=low(0.5))
