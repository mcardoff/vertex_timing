"""Reliability in the TIME weight and the clustering seeds, with the selection fixed at hit x0.5 and 1/(1+sigma_z0/2)."""
from sel import *
def rel(d):
    t = d['tracks']
    return np.where(t['nhgtd_hits'] >= 2, 1.0, 0.5) / (1 + t['sigma_z0_d'] / 2.0)
def zf(d): return 1 / (1 + d['tracks']['sigma_z0_d'] / 2.0)
def hf(d): return np.where(d['tracks']['nhgtd_hits'] >= 2, 1.0, 0.5)
def run(label, **kw):
    out = []
    for n, m in SAMPLES:
        d = get(n, m)
        c = cands(d, **{k: v(d) for k, v in kw.items()})
        ok = c['ok'][evmax_pick(c['ev'], c['S'], d['nev'])]
        b = d['ref']['tzp'].mean(); out.append((ok.mean() - b) / (1 - b) * 100)
    print(f"{label:52s} " + " ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f} mean {np.mean(out):+.2f}")
print("% of TZP fails removed:                                zjets  dijet    vbf  ttbar")
run("selection reliability only", sfac=rel)
run("+ time weight x 1/(1+sigma_z0/2)", sfac=rel, tfac=zf)
run("+ time weight x hit factor", sfac=rel, tfac=hf)
run("+ time weight x both", sfac=rel, tfac=rel)
run("+ seeds x reliability", sfac=rel, seedfac=rel)
run("+ seeds x hit factor", sfac=rel, seedfac=hf)
