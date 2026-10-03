"""Clustering scan 2: joint (t,z) clustering details."""
from reclus import *
print("% of production-TZP fails removed; columns zjets dijet vbf ttbar")
run("iterative (t,z) 3.0")
for cut in (2.5, 3.5):
    run(f"(t,z) {cut}", cut=cut, useZ=True)
for sc in (0.5, 0.75, 1.15, 1.5, 2.5):
    for cut in (3.0, 3.5):
        run(f"(t,z) {cut}, sigma_z x {sc}", cut=cut, useZ=True, szscale=sc)
for fl in (0.1, 0.3, 1.0):
    run(f"(t,z) 3.0, sigma_z floor {fl} mm", cut=3.0, useZ=True, szfloor=fl)
run("(t,z) 3.0, seed pT e^-dz", cut=3.0, useZ=True, seedw=lambda d: d['tracks']['w_tzp'])
run("simultaneous (t,z) 3.0", method=SIMUL, cut=3.0, useZ=True)
run("cone (t,z) 3.0", method=CONE, cut=3.0, useZ=True)
