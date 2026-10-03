"""Round 7: combinations."""
import itertools, numpy as np
from exp7 import run, kde

print("% of fails removed:                                                      zjets   dijet   vbf    ttbar")
res = []
for cap, p, hg, ptt in itertools.product((None, 5, 8), (1.0, 0.5), (0.0, 0.5), (1.0, 0.5)):
    if cap is not None and p != 1.0: continue
    lab = f"cap {cap} p {p} hgtd {hg} time-pT^{ptt}"
    res.append((lab, run(lab, lambda d: kde(d, cap=cap, p=p, hg=hg, pt_t=ptt))))
print("-- around cap 5 hg .5 ptt .5 --")
for kw in (dict(W=50.), dict(W=75.), dict(Wt=75.), dict(a=0.7), dict(a=1.3), dict(b=0.4), dict(b=0.8), dict(g=0.3), dict(g=0.7), dict(tau=60.),
           dict(pd=0.25), dict(pd=-0.25), dict(soft_env=True)):
    base = dict(cap=5, hg=0.5, pt_t=0.5); base.update(kw)
    run(str(kw), lambda d: kde(d, **base))
