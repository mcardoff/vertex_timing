"""Frozen reliability-weighted configuration: in-sample, out-of-sample, mu=0, by n HS tracks, and the quality flag."""
from sel import *


def zf(d): return 1 / (1 + d['tracks']['sigma_z0_d'] / 2.0)
def rel(d): return np.where(d['tracks']['nhgtd_hits'] >= 2, 1.0, 0.5) * zf(d)


SETS = [('zjets_novbs', None, 0), ('dijet_novbs', None, 0), ('vbf_novbs', 120, 0), ('ttbar_novbs', 60, 0),
        ('vbf_novbs', 300, 120), ('ttbar_novbs', 180, 60), ('zjets_mjj500p0', None, 0), ('dijet_mjj500p0', None, 0),
        ('vbf_mjj500p0', 500, 380), ('ttbar_mjj500p0', 400, 300),
        ('vbf_mu0_novbs', 20, 0), ('ttbar_mu0_novbs', 20, 0), ('zeejets_mu0_novbs', None, 0)]
print(f"{'set':28s} {'events':>7s} {'TZP':>6s} {'KDE_TZ':>7s} {'+rel':>6s} | fails removed KDE_TZ / +rel | Q>=2: provided / right (old) -> (new)")
for name, mx, mn in SETS:
    d = setup(name, mx, mn); t = d['tracks']
    b = evaluate(d, cluster(d), kde=False)['tzp']
    c0 = cands(d); p0, _, Q0 = top2(d, c0); ok0 = c0['ok'][p0]
    c1 = cands(d, sfac=rel(d), tfac=zf(d)); p1, _, Q1 = top2(d, c1); ok1 = c1['ok'][p1]
    f = 1 - b.mean()
    print(f"{name+f'[{mn}:{mx}]':28s} {d['nev']:7d} {b.mean()*100:6.2f} {ok0.mean()*100:7.2f} {ok1.mean()*100:6.2f} | {(ok0.mean()-b.mean())/f*100:+6.1f}% / {(ok1.mean()-b.mean())/f*100:+6.1f}% | "
          f"{np.mean(Q0>=2)*100:5.1f} / {ok0[Q0>=2].mean()*100:5.1f} -> {np.mean(Q1>=2)*100:5.1f} / {ok1[Q1>=2].mean()*100:5.1f}   pass counts {b.sum()} {ok0.sum()} {ok1.sum()}")
    if mn == 0 and 'mu0' not in name and 'mjj' not in name:
        nhs = np.bincount(t['ev'], t['truth_is_hs'], d['nev'])
        print("     by nHS 0/1/2-3/4-7/8+ (TZP -> new): " + " ".join(f"{b[s].mean()*100:.1f}->{ok1[s].mean()*100:.1f}" for s in (nhs == 0, nhs == 1, (nhs >= 2) & (nhs <= 3), (nhs >= 4) & (nhs <= 7), nhs >= 8))
              + "   new, right in Q bins <1/1-2/2-3/3-4/>4: " + " ".join(f"{ok1[(Q1>=a)&(Q1<bb)].mean()*100:.0f}" for a, bb in ((-9, 1), (1, 2), (2, 3), (3, 4), (4, 1e9))))
    del d
