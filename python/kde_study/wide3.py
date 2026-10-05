"""Track-list scan for one sample (compact): accounting of right-timed HS tracks, then the candidate lists.
   wide3.py <trackdump.root> [max_file_idx]"""
import sys
from wide import *

path = sys.argv[1]; mf = int(sys.argv[2]) if len(sys.argv) > 2 else None
ev, tr = load(path, mf); nev = len(ev['vz'])
hs = tr['truth_is_hs'] > 0.5
right = np.abs(tr['time'] - ev['ttruth'][tr['ev']]) < 3 * tr['timeRes']
nom = nominal(tr); q = tr['quality'] > 0.5
zero = np.bincount(tr['ev'][nom & hs & right], minlength=nev) == 0
print(f"== {path}: {nev} events; {nom.sum()/nev:.1f} nominal tracks/event; no right-timed HS track in the nominal list: {zero.mean()*100:.1f}% of events")
cats = [("nominal list", nom), ("pT > 30 GeV", nominal(tr, 3.0, 30.0, 1e9)), ("z0 signif 3-4", nominal(tr, 4.0, 1.0, 1e9) & ~nominal(tr, 3.0, 1.0, 1e9)),
        ("z0 signif 4-8", nominal(tr, 8.0, 1.0, 1e9) & ~nominal(tr, 4.0, 1.0, 1e9)), ("pT 0.5-1 GeV, signif < 3", nominal(tr, 3.0, 0.5, 1.0)),
        ("pT 0.5-1 GeV, signif 3-4.5", nominal(tr, 4.5, 0.5, 1.0) & ~nominal(tr, 3.0, 0.5, 1.0))]
print(f"  {'category':28s} {'tracks/ev':>9s} {'HS purity':>9s} {'HS time right':>13s}   zero-HS events with a right-timed HS track here")
for name, m in cats:
    nh = np.bincount(tr['ev'][m & hs & right], minlength=nev)
    print(f"  {name:28s} {m.sum()/nev:9.2f} {100*(m&hs).sum()/max(m.sum(),1):8.1f}% {100*(m&hs&right).sum()/max((m&hs).sum(),1):12.1f}%   {100*np.mean(nh[zero] > 0):5.1f}%")
base = evaluate(ev, tr, nom)
print(f"  nominal list: core {base.mean()*100:.2f}%  ({base.sum()} / {nev});  on the zero-HS events {base[zero].mean()*100:.1f}%")
low = lambda f: (lambda t: np.where(t['pt_d'] < 1.0, f, 1.0))


def run(label, mask, **kw):
    ok = evaluate(ev, tr, mask, **kw)
    print(f"  {label:54s} core {ok.mean()*100:6.2f}  fails removed {100*(ok.mean()-base.mean())/(1-base.mean()):+6.2f}%  tracks/ev {mask.sum()/nev:5.1f}  zero-HS ev {ok[zero].mean()*100:5.1f}")


run("no pT ceiling", nominal(tr, 3.0, 1.0, 1e9))
run("pT ceiling 100", nominal(tr, 3.0, 1.0, 100.))
for ns in (3.5, 4.0, 5.0, 8.0):
    run(f"no ceiling, signif < {ns}", nominal(tr, ns, 1.0, 1e9))
run("signif < 4 (ceiling kept)", nominal(tr, 4.0))
twohit = (tr['pt'] >= 1) | (tr['nhgtd_hits'] >= 2)
run("no ceiling, signif < 4, pT > 0.5", nominal(tr, 4.0, 0.5, 1e9))
run("no ceiling, signif < 4, pT > 0.5, sub-GeV x0.5", nominal(tr, 4.0, 0.5, 1e9), extra=low(0.5))
run("no ceiling, signif < 4, pT > 0.75, sub-GeV x0.5", nominal(tr, 4.0, 0.75, 1e9), extra=low(0.5))
run("no ceiling, signif < 4, sub-GeV: >= 2 hits", nominal(tr, 4.0, 0.5, 1e9) & twohit)
run("no ceiling, signif < 4, sub-GeV: >= 2 hits, x0.5", nominal(tr, 4.0, 0.5, 1e9) & twohit, extra=low(0.5))
run("no ceiling, signif < 4, sub-GeV: >= 2 hits, x0.25", nominal(tr, 4.0, 0.5, 1e9) & twohit, extra=low(0.25))
run("no ceiling, >=1 GeV signif < 4, sub-GeV signif < 3, 2 hits, x0.5", (nominal(tr, 4.0, 1.0, 1e9) | (nominal(tr, 3.0, 0.5, 1.0) & (tr['nhgtd_hits'] >= 2))), extra=low(0.5))
