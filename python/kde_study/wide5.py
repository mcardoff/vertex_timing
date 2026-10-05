"""Chosen wide list on one sample, with its pieces: wide5.py <trackdump.root> [max_file_idx]"""
import sys
from wide import *
path = sys.argv[1]; mf = int(sys.argv[2]) if len(sys.argv) > 2 else None
ev, tr = load(path, mf); nev = len(ev['vz'])
hs = tr['truth_is_hs'] > 0.5; right = np.abs(tr['time'] - ev['ttruth'][tr['ev']]) < 3 * tr['timeRes']
nom = nominal(tr); base = evaluate(ev, tr, nom)
zero = np.bincount(tr['ev'][nom & hs & right], minlength=nev) == 0
nhs = np.bincount(tr['ev'][nom & hs], minlength=nev)
print(f"== {path}: {nev} events; nominal list core {base.mean()*100:.2f} ({base.sum()}); events with no right-timed HS track {zero.mean()*100:.1f}%")
def run(label, m):
    ok = evaluate(ev, tr, m)
    z2 = np.bincount(tr['ev'][m & hs & right], minlength=nev) == 0
    print(f"  {label:50s} core {ok.mean()*100:6.2f} ({ok.sum()})  fails removed {100*(ok.mean()-base.mean())/(1-base.mean()):+6.2f}%  tracks/ev {m.sum()/nev:5.1f}  no-right-HS events {z2.mean()*100:4.1f}%  gained/lost {np.mean(ok&~base)*100:.2f}/{np.mean(~ok&base)*100:.2f}")
    return ok
run("no pT ceiling", nominal(tr, 3.0, 1.0, 1e9))
run("+ signif < 4 above 1 GeV", nominal(tr, 4.0, 1.0, 1e9))
wide = nominal(tr, 4.0, 1.0, 1e9) | (nominal(tr, 3.0, 0.5, 1.0) & (tr['nhgtd_hits'] >= 2))
ok = run("+ 0.5-1 GeV tracks with >= 2 hits, signif < 3", wide)
run("  (same, all hit counts)", nominal(tr, 4.0, 1.0, 1e9) | nominal(tr, 3.0, 0.5, 1.0))
print("  by nominal-list HS track count 0/1/2-3/4-7/8+: " + " ".join(f"{base[s].mean()*100:.1f}->{ok[s].mean()*100:.1f}" for s in (nhs == 0, nhs == 1, (nhs >= 2) & (nhs <= 3), (nhs >= 4) & (nhs <= 7), nhs >= 8)))
