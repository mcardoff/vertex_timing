"""Clustering scan 3: z-aware kernel on top of time-only and (t,z) clusterings."""
from reclus import *
sw = lambda d: d['tracks']['w_tzp']
TZ = dict(cut=3.0, useZ=True, seedw=sw)
print("% of production-TZP fails removed; columns zjets dijet vbf ttbar")
run("(t,z) 3.0 seed w, time-only kernel", **TZ)
for zk in (1.0, 2.0, 3.0, 5.0):
    run(f"(t,z), z-kernel {zk} sigma (density)", kde_kw=dict(zk=zk), **TZ)
    run(f"(t,z), z-kernel {zk} sigma (density+time)", kde_kw=dict(zk=zk, zk_time=True), **TZ)
for zk in (2.0, 3.0):
    run(f"(t,z), z-kernel {zk} sigma + 0.5mm floor (both)", kde_kw=dict(zk=zk, zk_time=True, zfloor=0.5), **TZ)
    run(f"production clusters, z-kernel {zk} (both)", kde_kw=dict(zk=zk, zk_time=True))
for W in (45., 75.):
    run(f"(t,z), W {W}", kde_kw=dict(W=W), **TZ)
