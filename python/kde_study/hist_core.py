"""Inclusive core fraction per score from a clustering_hist file (good_/flat_ histograms, all bins incl. flow)."""
import sys, uproot
f = uproot.open(sys.argv[1])
key = sys.argv[2] if len(sys.argv) > 2 else "n_Forward_Jets"
keys = [k.split(';')[0] for k in f.keys()]
ref = None
for k in keys:
    if not k.startswith("flat_" + key + "_"): continue
    name = k[len("flat_" + key + "_"):]
    tot = f[k].values(flow=True).sum()
    good = f["good_" + key + "_" + name].values(flow=True).sum()
    print(f"{name:70s} pass {good:8.0f} total {tot:8.0f} fail {tot-good:7.0f}  core {100*good/tot:6.2f}%")
