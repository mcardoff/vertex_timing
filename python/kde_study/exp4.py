"""Round 3: re-time the SELECTED cluster (selection fixed = TZP)."""
import sys, numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core

SAMPLES = [('zjets_novbs', None), ('dijet_novbs', None), ('vbf_novbs', 120), ('ttbar_novbs', 60)]
DS = {n: prep(n, m) for n, m in SAMPLES}
for n, _ in SAMPLES:
    d = DS[n]; c = d['clusters']; t = d['tracks']
    d['S0'] = tzp_score(d); d['T0'] = tzp_time(c)
    d['pick'] = evmax_pick(c['ev'], d['S0'], d['nev'])
    d['base_ok'] = core(d, d['pick'], d['T0'])


def wmean(t, w, nc, sel=None):
    if sel is None: sel = np.ones(len(w), bool)
    num = np.bincount(t['crow'][sel], (w * t['time'])[sel], nc)
    den = np.bincount(t['crow'][sel], w[sel], nc)
    with np.errstate(invalid='ignore', divide='ignore'):
        return num / den, den


def run(label, tfn):
    out = []
    for n, _ in SAMPLES:
        d = DS[n]
        T = tfn(d)
        T = np.where(np.isfinite(T), T, d['T0'])
        ok = core(d, d['pick'], T)
        base = d['base_ok']
        out.append((ok.mean() - base.mean()) * 100 / (1 - base.mean()))
    print(f"{label:58s} " + "  ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f}  mean {np.mean(out):+.2f}")


print("% of fails removed:                                              zjets   dijet   vbf    ttbar")
run("plain inv-var full cluster", lambda d: d['clusters']['cluster_time'])
def truth_hs(d):
    t = d['tracks']; nc = len(d['clusters']['ev'])
    return wmean(t, 1 / t['timeRes'] ** 2, nc, t['truth_is_hs'] > 0.5)[0]
run("ORACLE: inv-var over truth-HS tracks in the cluster", truth_hs)
for a in (0.5, 1.0, 2.0):
    for p in (0.5, 1.0):
        def f(d, a=a, p=p):
            t = d['tracks']; nc = len(d['clusters']['ev'])
            return wmean(t, t['pt'] ** p * np.exp(-a * t['dz']) / t['timeRes'] ** 2, nc)[0]
        run(f"w = pT^{p} e^-{a}dz / sig_t^2", f)
def f(d):
    t = d['tracks']; nc = len(d['clusters']['ev'])
    return wmean(t, np.exp(-1.0 * t['dz']) / t['timeRes'] ** 2, nc)[0]
run("w = e^-dz / sig_t^2", f)
# trimmed: iteratively drop tracks > k sigma from the weighted mean (2 passes)
for k in (1.5, 2.0):
    for wt in ('iv', 'tzp'):
        def f(d, k=k, wt=wt):
            t = d['tracks']; nc = len(d['clusters']['ev'])
            w = 1 / t['timeRes'] ** 2
            if wt == 'tzp': w = w * t['pt'] * np.exp(-t['dz'])
            m, den = wmean(t, w, nc)
            for _ in range(2):
                sel = np.abs(t['time'] - m[t['crow']]) < k * t['timeRes']
                m2, den2 = wmean(t, w, nc, sel)
                m = np.where(den2 > 0, m2, m)
            return m
        run(f"trim {k} sigma x2, w={wt}", f)
# guarded combos: in-jet (>=3) else weighted
def f(d):
    t = d['tracks']; c = d['clusters']; nc = len(c['ev'])
    m = wmean(t, t['pt'] * np.exp(-t['dz']) / t['timeRes'] ** 2, nc)[0]
    use = (c['n_in_jet_timed'] >= 3) & np.isfinite(c['t_in_jet'])
    return np.where(use, c['t_in_jet'], m)
run("in-jet(>=3) else TZP-weighted", f)
