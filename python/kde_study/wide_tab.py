from wide import *

if __name__ == "__main__":
    path = sys.argv[1]; mf = int(sys.argv[2]) if len(sys.argv) > 2 else None
    ev, tr = load(path, mf); nev = len(ev['vz'])
    hs = tr['truth_is_hs'] > 0.5
    right = np.abs(tr['time'] - ev['ttruth'][tr['ev']]) < 3 * tr['timeRes']
    nom = nominal(tr)
    print(f"{path}: {nev} events in the nominal population, {len(tr['ev'])/nev:.0f} dumped tracks/event, {nom.sum()/nev:.1f} nominal")
    q = tr['quality'] > 0.5
    cats = [("nominal list", nom),
            ("z0 signif 3-4 (else nominal)", nominal(tr, 4.0) & ~nom),
            ("z0 signif 4-5", nominal(tr, 5.0) & ~nominal(tr, 4.0)),
            ("z0 signif 5-8", nominal(tr, 8.0) & ~nominal(tr, 5.0)),
            ("z0 signif > 8", nominal(tr, 1e9) & ~nominal(tr, 8.0)),
            ("pT 0.5-1 GeV (signif < 3, quality)", nominal(tr, 3.0, 0.5, 1.0)),
            ("pT < 0.5 GeV", nominal(tr, 3.0, 0.0, 0.5)),
            ("pT > 30 GeV", nominal(tr, 3.0, 30.0, 1e9)),
            ("fails quality (signif < 3, 1-30 GeV)", (tr['nsig'] < 3) & (tr['pt'] > 1) & (tr['pt'] < 30) & ~q)]
    nh_nom = np.bincount(tr['ev'][nom & hs & right], minlength=nev)
    zero = nh_nom == 0
    print(f"events with NO right-timed HS track in the nominal list: {zero.mean()*100:.1f}%")
    print(f"{'category':38s} {'tracks/ev':>9s} {'HS':>7s} {'HS+right':>9s} {'purity':>7s} {'HS time right':>13s}   of the zero-HS events, has a right-timed HS track here")
    for name, m in cats:
        nh = np.bincount(tr['ev'][m & hs & right], minlength=nev)
        print(f"{name:38s} {m.sum()/nev:9.2f} {(m&hs).sum()/nev:7.2f} {(m&hs&right).sum()/nev:9.2f} {100*(m&hs).sum()/max(m.sum(),1):6.1f}% {100*(m&hs&right).sum()/max((m&hs).sum(),1):12.1f}%   {100*np.mean(nh[zero] > 0):5.1f}%")
    anyh = np.bincount(tr['ev'][hs & right & ~nom], minlength=nev)
    print(f"of the zero-HS events, any right-timed HS track anywhere in the dump: {100*np.mean(anyh[zero]>0):.1f}%")

    base = evaluate(ev, tr, nom)
    print(f"\nnominal list, TZP_KDE_TZ + reliability: core {base.mean()*100:.2f}%  ({base.sum()} / {nev})")
    def run(label, mask, **kw):
        ok = evaluate(ev, tr, mask, **kw)
        print(f"{label:46s} core {ok.mean()*100:6.2f}  fails removed vs nominal list {100*(ok.mean()-base.mean())/(1-base.mean()):+6.2f}%   tracks/ev {mask.sum()/nev:5.1f}   zero-HS events: {ok[zero].mean()*100:.1f} (was {base[zero].mean()*100:.1f})")
    for ns in (2.0, 2.5, 3.5, 4.0, 5.0, 8.0):
        run(f"z0 signif < {ns}", nominal(tr, ns))
    for lo in (0.5, 0.75):
        run(f"pT > {lo}", nominal(tr, 3.0, lo))
    run("no pT ceiling", nominal(tr, 3.0, 1.0, 1e9))
    run("no quality cut", nominal(tr, quality=False))
    run("signif < 5, pT > 0.5, no ceiling", nominal(tr, 5.0, 0.5, 1e9))
    run("|dz| < 2 mm (no significance cut)", (tr['dz'] < 2.0) & (tr['pt'] > 1) & (tr['pt'] < 30) & q)
    run("|dz| < 3 mm", (tr['dz'] < 3.0) & (tr['pt'] > 1) & (tr['pt'] < 30) & q)
    run("signif < 3 OR |dz| < 1 mm", (nominal(tr, 1e9) & ((tr['nsig'] < 3) | (tr['dz'] < 1.0))))
    run("signif < 3 AND |dz| < 5 mm", nom & (tr['dz'] < 5.0))
