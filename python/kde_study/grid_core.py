"""Inclusive core-fraction table for several clustering_hist files: grid_core.py <label>=<file> ..."""
import sys, uproot
ROWS = [("HGTD (Athena)", "HGTD Algorithm"), ("TRKPTZ", "[Baseline Algorithm]"), ("WAVeS", "WAVeS Score_"),
        ("TZP", "d0 precision]"), ("TZP_KDE", "[mean-shift t]"), ("TZP_KDE_TZ", "[(t,z) clusters]")]
key = "n_Forward_Jets"
for arg in sys.argv[1:]:
    lab, path = arg.split("=")
    f = uproot.open(path)
    names = [k.split(';')[0] for k in f.keys() if k.startswith("flat_" + key + "_")]
    res = {}
    for short, pat in ROWS:
        n = [x for x in names if pat in x and (short != "WAVeS" or x.startswith("flat_" + key + "_WAVeS Score_"))][0]
        tot = f[n].values(flow=True).sum(); good = f["good_" + n[5:]].values(flow=True).sum()
        res[short] = (good, tot)
    tz = res["TZP"]
    print(f"{lab}: events {tz[1]:.0f}")
    for short, _ in ROWS:
        g, t = res[short]
        extra = f"   fails removed vs TZP {((t - tz[0]) - (t - g)) / (tz[1] - tz[0]) * 100:+.1f}%" if short.startswith("TZP_") and t == tz[1] else ""
        print(f"   {short:14s} pass {g:8.0f}  fail {t-g:7.0f}  core {100*g/t:6.2f}%{extra}")
