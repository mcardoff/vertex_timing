"""Clustering scan 1: distance cut, method, seed ordering, (t,z), resolution floor, fixed windows, mean-shift modes."""
from reclus import *

print("% of production-TZP fails removed; columns zjets dijet vbf ttbar")
run("iterative 3.0 (production)")
for cut in (1.5, 2.0, 2.5, 3.5, 4.0):
    run(f"iterative {cut}", cut=cut)
for cut in (2.0, 3.0, 4.0):
    run(f"simultaneous {cut}", method=SIMUL, cut=cut)
    run(f"cone {cut}", method=CONE, cut=cut)
run("iterative 3.0, seed = pT e^-dz", seedw=lambda d: d['tracks']['w_tzp'])
run("iterative 3.0, seed = sqrt(pT) e^-dz", seedw=lambda d: d['tracks']['pt_d'] ** 0.5 * np.exp(-d['tracks']['dz_d']))
run("iterative 3.0, seed = e^-dz", seedw=lambda d: np.exp(-d['tracks']['dz_d']))
for cut in (3.0, 4.0):
    run(f"iterative (t,z) {cut}", cut=cut, useZ=True)
for fl in (15., 30.):
    for cut in (2.0, 3.0):
        run(f"iterative {cut}, floor {fl} ps", cut=cut, floor=fl)
for w in (40., 60., 90., 120.):
    run(f"fixed window {w} ps, weighted centroid", method=WINDOW, cut=w)
for W in (20., 30., 45.):
    for eps in (5., 15.):
        run(f"mean-shift modes W {W} merge {eps}", method=MSHIFT, par=W, cut=eps)
