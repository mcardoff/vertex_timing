"""Reliability factors on top of the single-hit down-weight (x0.5): sigma_z0 (extrapolation quality), pT power, combined."""
from sel4 import *

H = lambda d: np.where(d['tracks']['nhgtd_hits'] >= 2, 1.0, 0.5)
print("\n-- on top of single-hit x0.5 --")
run("single-hit x0.5 (new reference)", lambda d, c: dens_with(d, c, H(d)) * c['env'])
for s in (1.0, 2.0, 4.0):
    run(f"x 1/(1 + sigma_z0/{s})", lambda d, c, s=s: dens_with(d, c, H(d) / (1 + d['tracks']['sigma_z0_d'] / s)) * c['env'])
for cap in (1.2, 2.4):
    run(f"x 0.7 if sigma_z0 > {cap}", lambda d, c, cap=cap: dens_with(d, c, H(d) * np.where(d['tracks']['sigma_z0_d'] > cap, 0.7, 1.0)) * c['env'])
for p in (-0.15, 0.15, 0.25):
    run(f"x pT^{p} (total power {0.5+p})", lambda d, c, p=p: dens_with(d, c, H(d) * d['tracks']['pt_d'] ** p) * c['env'])
for a in (0.35, 0.65):
    run(f"single-hit x{a}", lambda d, c, a=a: dens_with(d, c, np.where(d['tracks']['nhgtd_hits'] >= 2, 1.0, a)) * c['env'])
# reliability also in the cluster ENVELOPE's d0 information sum and in the z centroid? (hard cluster) -- test info with the hit factor
def f(d, c):
    t = d['tracks']; nc = len(c['ev'])
    info = np.bincount(c['crow'], H(d) / t['sigma_d0_d'] ** 2, nc)
    return dens_with(d, c, H(d)) * np.exp(-0.6 * np.abs(c['dz'])) * info ** 0.225
run("+ hit factor inside the d0 information sum", f)
# single-hit tracks whose time has no multi-hit neighbour within 2 sigma: orphan single hits count less
def g(d, c, lone=0.25):
    t = d['tracks']; pc, pt_, dt, K, base, wt = d['_pairs']
    multi = (t['nhgtd_hits'] >= 2).astype(float)
    return None
