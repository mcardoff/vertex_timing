"""Selection terms suggested by the tie-break scan: time-reliability of the tracks in the density, |eta|, and a
beam-time prior on the candidate. Reports % of production-TZP fails removed (reference: 9.2 / 14.3 / 15.1 / 16.0)."""
from sel import *

D = {}
for n, m in SAMPLES:
    d = get(n, m); c = cands(d); D[n] = (d, c)
    d['tzp_ok'] = d['ref']['tzp']


def run(label, fn):
    """fn(d, c) -> candidate score"""
    out = []
    for n, _ in SAMPLES:
        d, c = D[n]
        ok = c['ok'][evmax_pick(c['ev'], fn(d, c), d['nev'])]
        b = d['tzp_ok'].mean()
        out.append((ok.mean() - b) / (1 - b) * 100)
    print(f"{label:46s} " + " ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f} mean {np.mean(out):+.2f}")
    return out


def dens_with(d, c, extra):
    pc, pt_, dt, K, base, wt = d['_pairs']
    return np.bincount(pc, (base * extra)[pt_] * K, len(c['ev']))


if __name__ == "__main__":
    print("% of TZP fails removed:                          zjets  dijet    vbf  ttbar")
    run("reference (TZP_KDE_TZ)", lambda d, c: c['S'])
    for k in (0.5, 1.0, 2.0):
        run(f"density weight x (25ps/sigma_t)^{k}", lambda d, c, k=k: dens_with(d, c, (25. / d['tracks']['timeRes_d']) ** k) * c['env'])
    for f2 in (1.25, 1.5, 2.0):
        run(f"density weight x {f2} if nHGTD >= 2", lambda d, c, f2=f2: dens_with(d, c, np.where(d['tracks']['nhgtd_hits'] >= 2, f2, 1.0)) * c['env'])
    for k in (0.5, 1.0, 2.0):
        run(f"density weight x e^(-{k}(|eta|-2.4))", lambda d, c, k=k: dens_with(d, c, np.exp(-k * (np.abs(d['tracks']['eta'].astype(float)) - 2.4))) * c['env'])
    for k in (0.1, 0.25, 0.5, 1.0):
        for sg in (175., 250.):
            run(f"S x exp(-{k} t^2 / 2 {sg:.0f}^2)", lambda d, c, k=k, sg=sg: c['S'] * np.exp(-k * 0.5 * (c['T'] / sg) ** 2))
    # event-centred time: candidate time relative to the weighted mean time of all event tracks
    for k in (0.25, 0.5):
        def f(d, c, k=k):
            t = d['tracks']; w = t['w_tzp']
            mu = np.bincount(t['ev'], w * t['time_d'], d['nev']) / np.bincount(t['ev'], w, d['nev'])
            return c['S'] * np.exp(-k * 0.5 * ((c['T'] - mu[c['ev']]) / 175.) ** 2)
        run(f"S x exp(-{k} (t - <t>_ev)^2 / 2 175^2)", f)
