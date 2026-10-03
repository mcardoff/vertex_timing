"""Clustering scan 4: envelope re-tune on (t,z) clusters."""
from reclus import *
sw = lambda d: d['tracks']['w_tzp']
TZ = dict(cut=3.0, useZ=True, seedw=sw)
print("% of production-TZP fails removed; columns zjets dijet vbf ttbar")
for b in (0.3, 0.6, 1.0, 1.5, 2.5):
    for g in (0.45, 0.9):
        run(f"(t,z): beta {b} gamma {g}", kde_kw=dict(feat=dict(b=b, g=g)), **TZ)
for zw in ('pt', 'ptprec'):
    run(f"(t,z): z centroid weight {zw}", kde_kw=dict(feat=dict(zw=zw)), **TZ)
    run(f"(t,z): z centroid weight {zw}, beta 1.5", kde_kw=dict(feat=dict(zw=zw, b=1.5)), **TZ)
