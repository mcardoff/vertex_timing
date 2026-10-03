"""Round 2: per-track weight forms. Score = [SUM_t w_t] * exp(-b|dz_C|) * info^(g/2)."""
import sys, numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core

SAMPLES = [('zjets_novbs', None), ('dijet_novbs', None), ('vbf_novbs', 120), ('ttbar_novbs', 60)]
DS = {n: prep(n, m) for n, m in SAMPLES}
for n, _ in SAMPLES:
    d = DS[n]; c = d['clusters']; t = d['tracks']
    d['S0'] = tzp_score(d); d['T0'] = tzp_time(c)
    d['base_ok'] = core(d, evmax_pick(c['ev'], d['S0'], d['nev']), d['T0'])
    d['info'] = np.bincount(t['crow'], 1 / t['sigma_d0'] ** 2, len(c['ev']))


def run(label, wfn, b=0.6, g=0.45, quiet=False):
    out = []
    for n, _ in SAMPLES:
        d = DS[n]; c = d['clusters']; t = d['tracks']
        S = np.bincount(t['crow'], wfn(t), len(c['ev'])) * np.exp(-b * np.abs(c['delta_z'])) * d['info'] ** (0.5 * g)
        ok = core(d, evmax_pick(c['ev'], S, d['nev']), d['T0'])
        base = d['base_ok']
        out.append((ok.mean() - base.mean()) * 100 / (1 - base.mean()))
    if not quiet:
        print(f"{label:52s} " + "  ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f}  mean {np.mean(out):+.2f}")
    return out


if __name__ == "__main__":
    t = DS['zjets_novbs']['tracks']
    print("zjets tracks: HS frac", t['truth_is_hs'].mean(), " sigma_z0 quantiles", np.quantile(t['sigma_z0'], [.1, .5, .9]),
          " dz q", np.quantile(t['dz'], [.1, .5, .9]))
    # purity table in (dz, sigma_z0)
    for name in ('zjets_novbs', 'vbf_novbs'):
        t = DS[name]['tracks']
        print(name, "HS purity by sigma_z0 (rows) x dz (cols)")
        se = [0, .1, .2, .4, .8, 1.6, 99]; de = [0, .1, .25, .5, 1, 2, 99]
        for i in range(6):
            row = []
            for j in range(6):
                m = (t['sigma_z0'] >= se[i]) & (t['sigma_z0'] < se[i + 1]) & (t['dz'] >= de[j]) & (t['dz'] < de[j + 1])
                row.append(f"{t['truth_is_hs'][m].mean()*100:5.1f}({m.mean()*100:4.1f})" if m.sum() > 50 else "    -      ")
            print(f"  sig {se[i]:.1f}-{se[i+1]:.1f}: " + " ".join(row))
        pe = [1, 1.5, 2, 3, 5, 10, 30]
        print("  purity by pT:", [f"{t['truth_is_hs'][(t['pt']>=pe[i])&(t['pt']<pe[i+1])].mean()*100:.1f}" for i in range(6)])
        print("  purity by nhgtd:", [f"{t['truth_is_hs'][t['nhgtd_hits']==k].mean()*100:.1f}" for k in (1, 2, 3, 4)])
        print("  purity by quality:", [(q, f"{t['truth_is_hs'][t['quality']==q].mean()*100:.1f}", int((t['quality']==q).sum())) for q in np.unique(t['quality'])][:6])
    print("\n% of fails removed:                                          zjets   dijet   vbf    ttbar")
    run("baseline", lambda t: t['pt'] * np.exp(-t['dz']))
    for cap in (4, 6, 8):
        run(f"cap {cap}", lambda t, cap=cap: np.minimum(t['pt'], cap) * np.exp(-t['dz']))
    for bb in (0.25, 0.5, 1.0):
        run(f"pT e^-dz sigma_z0^-{bb}, g=0", lambda t, bb=bb: t['pt'] * np.exp(-t['dz']) * t['sigma_z0'] ** -bb, g=0)
        run(f"pT e^-dz sigma_z0^-{bb}, g=.45", lambda t, bb=bb: t['pt'] * np.exp(-t['dz']) * t['sigma_z0'] ** -bb)
    for k in (1.0, 1.5, 2.0):
        for cc in (0.3, 1.0, 3.0):
            def w(t, k=k, cc=cc):
                s = k * t['sigma_z0']
                L = np.exp(-0.5 * (t['dz'] / s) ** 2) / (np.sqrt(2 * np.pi) * s)
                return t['pt'] * L / (L + cc)
            run(f"LR prob k {k} c {cc}, g=.45", w)
            run(f"LR prob k {k} c {cc}, g=0, b=0", w, b=0, g=0)
    for a in (0.5, 1.0):
        def w(t, a=a):
            return t['pt'] * np.exp(-a * t['dz']) * np.exp(-0.5 * np.minimum(t['z0_pull_pv'] ** 2, 50) / 2.0 ** 2)
        run(f"pT e^-{a}dz * gauss(pull/2)", w)
