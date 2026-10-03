"""Frozen configuration (no Athena term): KDE on time-only clusters vs KDE on (t,z) clusters. In-sample, out-of-sample, mu=0."""
from reclus import *
SETS = [('zjets_novbs', None, 0), ('dijet_novbs', None, 0), ('vbf_novbs', 120, 0), ('ttbar_novbs', 60, 0),
        ('vbf_novbs', 300, 120), ('ttbar_novbs', 180, 60), ('zjets_mjj500p0', None, 0), ('dijet_mjj500p0', None, 0),
        ('vbf_mjj500p0', 500, 380), ('ttbar_mjj500p0', 400, 300),
        ('vbf_mu0_novbs', 20, 0), ('ttbar_mu0_novbs', 20, 0), ('zeejets_mu0_novbs', None, 0)]
print(f"{'set':30s} {'events':>7s} {'TZP':>6s} | {'TZP on (t,z)':>12s} {'KDE':>6s} {'KDE on (t,z)':>12s} | fails removed: TZP(t,z) / KDE / KDE(t,z)")
for name, mx, mn in SETS:
    d = setup(name, mx, mn); t = d['tracks']
    r0 = evaluate(d, cluster(d))
    r1 = evaluate(d, cluster(d, cut=3.0, useZ=True, seedw=t['w_tzp']))
    b = r0['tzp'].mean(); f = 1 - b
    print(f"{name+f'[{mn}:{mx}]':30s} {d['nev']:7d} {b*100:6.2f} | {r1['tzp'].mean()*100:12.2f} {r0['kde'].mean()*100:6.2f} {r1['kde'].mean()*100:12.2f} | "
          f"{(r1['tzp'].mean()-b)/f*100:+6.1f}% / {(r0['kde'].mean()-b)/f*100:+6.1f}% / {(r1['kde'].mean()-b)/f*100:+6.1f}%   "
          f"pass counts {r0['tzp'].sum()} -> {r1['kde'].sum()}")
    if mx in (None, 120, 60) and mn == 0 and 'mu0' not in name and 'mjj' not in name:
        c, crow = features(d, cluster(d, cut=3.0, useZ=True, seedw=t['w_tzp'])); S, T = kde_eval(d, c, crow)
        dt1 = T[evmax_pick(c['ev'], S, d['nev'])] - d['ttruth']
        c0, crow0 = features(d, cluster(d)); dt0 = c0['t_tzp'][evmax_pick(c0['ev'], c0['S_tzp'], d['nev'])] - d['ttruth']
        nhs = np.bincount(t['ev'], t['truth_is_hs'], d['nev'])
        print("     windows 30/60/90/150: " + " ".join(f"{(np.mean(np.abs(dt0)>=w)-np.mean(np.abs(dt1)>=w))/np.mean(np.abs(dt0)>=w)*100:+.1f}%" for w in (30, 60, 90, 150))
              + "   by nHS 0/1/2-3/4-7/8+: " + " ".join(f"{np.mean(np.abs(dt0[s])<60)*100:.1f}->{np.mean(np.abs(dt1[s])<60)*100:.1f}" for s in (nhs == 0, nhs == 1, (nhs >= 2) & (nhs <= 3), (nhs >= 4) & (nhs <= 7), nhs >= 8))
              + f"   core RMS {dt0[np.abs(dt0)<60].std():.1f}->{dt1[np.abs(dt1)<60].std():.1f}")
    del d
