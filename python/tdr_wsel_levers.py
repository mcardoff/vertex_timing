"""The selection levers against the > 50 GeV capture (todo vt-5t0, results/wsel_levers.md): each t0 source's
core fraction, and its pileup-jet rejection paired against WAVeS and against perfect selection.

Two inputs, both written by the commit that added the waves_* lever t0 sources to util/tdr_common.h:
  --t0     a tdr_t0 tree: core fraction (|t0 - t_HS| < 60 ps over ALL events, a missing t0 a failure) of each cleaned
           t0, inclusively and by the event's hardest forward truth-pileup jet (pu_fwd_ptmax: none / 30-50 / > 50 GeV,
           the capture's population). Paired against WAVeS: the events each source rescues and breaks, and the net
           change's binomial-free error sqrt(rescued + broken).
  --pujet  a tdr_pujet tree: rejection (1 / pileup efficiency) at matched hard-scatter efficiency, forward jets
           (2.4 < |eta| < 4.0), 30-50 and > 50 GeV, unweighted and interpolated as python/tdr_pujet_plots.py does.
           Ratios between t0s carry an event bootstrap (Poisson weights per event, --folds), since every t0 gates the
           same jets: lever / WAVeS and lever / wfix_sel (perfect selection, the ceiling of a selector).

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_wsel_levers.py \\
        --t0 <tdr_t0.root> --pujet <tdr_pujet.root> [--sources waves_a0,waves_pv0,...]
"""
import argparse

import numpy as np
import uproot

LEVERS = ("waves_a0", "waves_a025", "waves_a05", "waves_pv0", "waves_npu")
LABEL = {"itk": "ITk only", "waves": "WAVeS", "wfix_sel": "perfect selection", "truth": "truth t0",
         "waves_a0": "WAVeS, pT_jet^0", "waves_a025": "WAVeS, pT_jet^0.25", "waves_a05": "WAVeS, pT_jet^0.5",
         "waves_pv0": "WAVeS x PV0 share", "waves_npu": "WAVeS x (1 - other-vertex share)",
         "wfix_time": "perfect timing", "wfix_both": "perfect timing + selection (C3)", "wfix_c4": "C3 + misclustering (C4)",
         "wfix_c5": "C4 + HS track efficiency (C5)", "wfix_c6": "C5 + HS time resolution (C6)",
         "truth_wsig": "truth t0, WAVeS sigma_t0 + inflation"}
PASS_PS = 60.0
EFFS = np.array([0.80, 0.85, 0.90, 0.95])
SLICES = [("30-50 GeV", 30.0, 50.0), ("> 50 GeV", 50.0, 1e9)]
PU_CATS = [("all events", None), ("no forward PU jet", (-1.0, 0.0)), ("forward PU jet 30-50 GeV", (0.0, 50.0)),
           ("forward PU jet > 50 GeV", (50.0, 1e9))]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--t0", help="tdr_t0 tree (core fractions)")
ap.add_argument("--pujet", help="tdr_pujet tree (rejection)")
ap.add_argument("--sources", default=",".join(LEVERS), help="comma list of lever t0s (default all five)")
ap.add_argument("--refs", default="waves,wfix_sel", help="comma list of the t0s each source is paired against (rejection)")
ap.add_argument("--folds", type=int, default=150, help="bootstrap folds (default 150, as results/lojo_t0.md)")
ap.add_argument("--seed", type=int, default=20261114)
args = ap.parse_args()
SRC = args.sources.split(",")
REFS = args.refs.split(",")


def core_fractions(path):
    E = uproot.open(path)["events"].arrays(library="np")
    dt = {"waves": E["dt_waves_clean"], "trkptz": E["dt_trkptz_clean"], "hgtd": E["dt_hgtd"]}
    for s in dict.fromkeys(["wfix_sel"] + SRC):
        dt[s] = E[f"dt_{s}_clean"]
    ok = {s: np.abs(np.nan_to_num(v, nan=1e9)) < PASS_PS for s, v in dt.items()}
    pu = E["pu_fwd_ptmax"]
    print(f"core fraction (|t0 - t_HS| < {PASS_PS:g} ps; no t0 = fail), {path}")
    order = list(dict.fromkeys(["hgtd", "waves", "wfix_sel"] + SRC))
    for name, rng in PU_CATS:
        m = np.ones(len(pu), bool) if rng is None else (pu > rng[0]) & (pu <= rng[1]) if rng[0] >= 0 else (pu <= 0)
        n = int(m.sum())
        print(f"\n  {name}: {n} events")
        print(f"  {'t0':36s} {'core %':>8s} {'vs WAVeS':>10s} {'rescued':>8s} {'broken':>7s} {'net / sigma':>12s}")
        for s in order:
            c = 100.0 * ok[s][m].mean()
            if s == "waves":
                print(f"  {LABEL.get(s, s):36s} {c:8.2f}")
                continue
            r = int(np.sum(ok[s][m] & ~ok["waves"][m]))
            b = int(np.sum(~ok[s][m] & ok["waves"][m]))
            sig = np.sqrt(r + b)
            d = c - 100.0 * ok["waves"][m].mean()
            print(f"  {LABEL.get(s, s):36s} {c:8.2f} {d:+10.2f} {r:8d} {b:7d} {((r - b) / sig if sig else 0):+12.1f}")


def roc_pe(hs, pu, w_hs, w_pu, effs):
    """pileup efficiency at the given hard-scatter efficiencies (pass = R_pT >= threshold): tdr_pujet_plots.at_eff
    with per-jet weights (all 1 reproduces it exactly)"""
    th = np.unique(np.concatenate([hs, pu]))
    def eff(x, w):
        o = np.argsort(x, kind="stable")
        xs, cw = x[o], np.concatenate([[0.0], np.cumsum(w[o])])
        below = cw[np.searchsorted(xs, th, side="left")]
        return 1.0 - below / cw[-1]
    e_hs, e_pu = eff(hs, w_hs), eff(pu, w_pu)
    u, inv = np.unique(e_hs, return_inverse=True)
    pmin = np.full(len(u), np.inf)
    np.minimum.at(pmin, inv, e_pu)
    return np.interp(effs, u, pmin)


def rejection(path):
    J = uproot.open(path)["jets"].arrays(library="np")
    cols = {"itk": "rpt_itk", "waves": "rpt_t0_waves", "wfix_sel": "rpt_t0_wfix_sel", "truth": "rpt_t0_truth"}
    cols.update({s: f"rpt_t0_{s}" for s in SRC + REFS if s != "itk"})
    ev = np.unique(np.stack([J["file_idx"].astype(np.int64), J["entry"].astype(np.int64)], 1), axis=0, return_inverse=True)[1].ravel()
    nev = ev.max() + 1
    rng = np.random.default_rng(args.seed)
    W = rng.poisson(1.0, size=(args.folds, nev)).astype(float)
    ae = np.abs(J["eta"])
    print(f"\npileup-jet rejection at matched HS efficiency, 2.4 < |eta| < 4.0, {path}")
    for lab, lo, hi in SLICES:
        m = (J["pt"] > lo) & (J["pt"] < hi) & (ae >= 2.4) & (ae < 4.0)
        H, P = m & (J["hs"] == 1), m & (J["pu"] == 1)
        print(f"\n  {lab}: {int(H.sum())} HS / {int(P.sum())} PU jets;  HS eff " + "  ".join(f"{e:.2f}" for e in EFFS))
        rej, boot = {}, {}
        for s, c in cols.items():
            x = J[c]
            inv = lambda pe: np.where(pe > 0, 1.0 / np.where(pe > 0, pe, 1.0), np.nan)   # a fold without a surviving PU jet: NaN
            rej[s] = inv(roc_pe(x[H], x[P], np.ones(H.sum()), np.ones(P.sum()), EFFS))
            boot[s] = np.array([inv(roc_pe(x[H], x[P], W[f][ev[H]], W[f][ev[P]], EFFS)) for f in range(args.folds)])
        for s in dict.fromkeys(["itk", "waves", "wfix_sel", "truth"] + SRC):
            print(f"  {LABEL.get(s, s):36s}" + "".join(f"{r:8.1f} (x{r / q:4.2f})" for r, q in zip(rej[s], rej["itk"])))
        for ref in REFS:
            print(f"  paired / {LABEL[ref]}:")
            for s in [x for x in dict.fromkeys(SRC + ["wfix_sel", "truth"]) if x != ref]:
                ratio = rej[s] / rej[ref]
                err = np.nanstd(boot[s] / boot[ref], axis=0)
                print(f"    {LABEL.get(s, s):34s}" + "".join(f"  x{q:5.3f} +- {e:5.3f}" for q, e in zip(ratio, err)))


if args.t0:
    core_fractions(args.t0)
if args.pujet:
    rejection(args.pujet)
