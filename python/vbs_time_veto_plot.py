#!/usr/bin/env python3
"""Signal / background efficiency of event-level VBS timing tests, weighted,
vs m_jj -- the region plots' columns, with efficiency in place of composition.

Reads util/vbs_time_veto outputs (tree `events`, one row per event passing the
region-plot selection) and applies, offline:
  * the recoil cut: truth E_T^miss > --met-min for the signal file (needs the
    raw-ntuple VBF run -- the skims have no truth record), pT(ll) > --zpt-min
    for the background file;
  * the jet-time calibration: sigma_jet x f(n), n = tracks in the jet-time
    estimate, f = 1.4826 x median |t_jet - t_HS| / sigma_jet over paper-HS
    legs in bins of n. By default (--jet-infl sig-ntrk) f(n) is measured ONCE,
    on the signal file, and applied to both files: the test must be one
    procedure for signal and background, and a per-file calibration would give
    the background a looser test than the signal (see results/vbs_time_veto.md);
  * a timing test at --nsigma (default 3), in two families:
    EVENT VETO -- the pair stays, the event is rejected:
      jj        both legs timed and |t_A - t_B| / sqrt(sA^2 + sB^2) >= n
      t0_<src>  the event has a t0 and ANY timed leg has
                |t_leg - t0| / sqrt(s_leg^2 + (infl sigma_t0)^2) >= n
    JET REMOVAL + RE-PAIRING (needs the jet_* arrays, vbs_time_veto >= the
    version that writes them) -- the jets failing the test go, the selection
    is re-applied to the rest, and the pair is re-formed:
      rp_<src>  every timed jet incompatible with t0 is removed; then >= 2
                jets, the max-m_jj opposite-hemisphere pair, and the m_jj and
                |Deta| cuts, as before. The region plots' ">= 1 jet in
                2.38 < |eta| < 4.0" preselection defines the starting sample
                and is NOT re-applied after removal (--repair-reapply-fwd
                does): it is not an analysis cut, and re-applying it drops
                ~4% of the signal whose valid pair survives (local VBF)
      rp_jj     no jet is removed; the pair is the max-m_jj pair whose two
                jets are time-compatible (untimed jets are compatible)
    A jet is "timed" when its jet-time estimate has >= --min-trk tracks.
Every efficiency is sum(w, kept) / sum(w, all) with the generator weight.
Re-paired events move in m_jj: efficiencies are binned in each event's NO-VETO
m_jj (the fraction of that column that survives anywhere), yields and S/sqrt(B)
in the m_jj of the pair the event ends up with.

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_time_veto_plot.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root \\
        --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root \\
        --tag "JVT + fJVT loose before pairing" --out figs/time_veto/jvtLoose

Writes <out>_{sig,bkg}_eff_mjj.(pdf|png), and with both files
<out>_s_over_sqrtb.(pdf|png) (S/sqrt(B) change vs no veto in %, per column),
<out>_s_over_sqrtb_yields.(pdf|png) (S/sqrt(B) from the yields at --lumi), the same two for S/B
(<out>_s_over_b, <out>_s_over_b_yields), <out>_jet_time_{pulls,calib} (the HS-leg jet-time pulls and f(n)),
<out>_t0_pulls(_bkg).(pdf|png) (the event t0s against the true vertex time, signal / background),
<out>_tradeoff.(pdf|png) (signal efficiency against background rejection 1/eps_B,
threshold scanned) and, with --compare, <out>_veto_vs_repair_tradeoff, plus
<out>_summary.md (with the scan's points appended).
"""
import argparse, os
import numpy as np
import uproot
import awkward as ak
import vbs_time_veto_common as C
from vbs_time_veto_common import NORM, EXPECTED_COMPLETE, T0_SRC, STYLE, NBINS, Calib, measure, eff, mjj_index

# 500 GeV columns above 1 TeV: 250 GeV ones there hold only 50-130 Z+jets MC events
# after the recoil cut. The last column includes the overflow.
MJJ_EDGES = [200, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000]
SPLIT = 2.4   # forward = |eta| >= this, as vbs_region_stack.py
# the event preselection re-applied after jet removal: >= MIN_PASSPT_JETS jets,
# >= 1 of them in MIN_ABS_ETA_JET < |eta| < MAX_ABS_ETA_JET (clustering_constants.h)
MIN_PASSPT_JETS, MIN_ABS_ETA_JET, MAX_ABS_ETA_JET = 2, 2.38, 4.00

VETO_METHODS = ["jj", "jj_all"] + ["t0_" + s for s in T0_SRC]
REPAIR_METHODS = ["rp_jj"] + ["rp_" + s for s in T0_SRC]
BANDS_MD = ["R1: fwd HS + fwd PU", "R2: fwd PU + central HS", "HS + PU, other eta",
            "HS + HS", "PU + PU, both fwd", "PU + PU, other eta", "neither label"]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", help="signal (VBF) vbs_time_veto output")
ap.add_argument("--bkg", help="background (Z+jets) vbs_time_veto output")
ap.add_argument("--sig-norm", default="vbf", choices=list(NORM) + ["none"])
ap.add_argument("--bkg-norm", default="zjets", choices=list(NORM) + ["none"])
ap.add_argument("--sig-label", default="#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.")
ap.add_argument("--bkg-label", default="#sqrt{s} = 14 TeV, HL-LHC, Z(#rightarrowll)+jets")
ap.add_argument("--met-min", type=float, default=200.0, help="signal: truth E_T^miss cut (GeV); 0 = none")
ap.add_argument("--zpt-min", type=float, default=200.0, help="background: pT(ll) cut (GeV); 0 = none")
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--jet-infl", default="sig-ntrk",
                help="sig-ntrk: f(n tracks) from the signal file's HS legs, for both files (default); "
                     "sig: one global factor from the signal file; per-file / per-file-ntrk: each file its own; "
                     "or a number")
ap.add_argument("--min-trk", type=int, default=1, help="tracks a jet-time estimate needs for the jet to count as timed")
ap.add_argument("--methods", default="jj,t0_trkptz,t0_waves,t0_hgtd,t0_truth",
                help="drawn methods, e.g. t0_trkptz,t0_waves (event veto) or rp_trkptz,rp_waves (jet "
                     "removal + re-pairing); the summary always has every method the files support")
ap.add_argument("--compare", default="",
                help="comma-separated t0 sources (e.g. trkptz,waves,hgtd,truth): also draw event veto "
                     "against jet removal + re-pairing for them (<out>_veto_vs_repair_tradeoff)")
ap.add_argument("--repair-reapply-fwd", action="store_true",
                help="after jet removal, re-require >= 1 jet in 2.38 < |eta| < 4.0 (the region plots' preselection)")
ap.add_argument("--mjj-edges", default=",".join(str(e) for e in MJJ_EDGES))
ap.add_argument("--lumi", type=float, default=3000.0, help="fb^-1")
ap.add_argument("--bf-hinv", type=float, default=1.0,
                help="B(H->inv) the signal yields assume (the VBF sample is H->ZZ->4nu kinematics, "
                     "normalised to sigma_VBF alone); scales S linearly, cancels in every relative number")
ap.add_argument("--tag", default="", help="selection note for the third label line (e.g. the JVT working point)")
ap.add_argument("--out", required=True, help="output stem")
ap.add_argument("--compact", action="store_true",
                help="the m_jj-column figures carry only the category line (e.g. 'both tagging jets in HGTD'), "
                     "not the selection / test / jet-time / calibration lines -- for slides that state those once")
ap.add_argument("--jet-pulls", default="both", choices=["both", "quoted"],
                help="<out>_jet_time_pulls: the quoted sigma and sigma x f(n) (both), or the quoted sigma alone")
C.add_acceptance_args(ap)
args = ap.parse_args()
if not args.sig and not args.bkg:
    raise SystemExit("give --sig and/or --bkg")
EDGES = np.array([float(x) for x in args.mjj_edges.split(",")])
DRAWN = args.methods.split(",")
for m in DRAWN:
    if m not in STYLE or m == "none": raise SystemExit(f"unknown method {m}; choose from {[k for k in STYLE if k != 'none']}")
COMPARE = [s for s in args.compare.split(",") if s]
for s in COMPARE:
    if s not in T0_SRC: raise SystemExit(f"--compare: unknown t0 source {s}; choose from {T0_SRC}")


def load(path, recoil_col, recoil_min, what):
    t = uproot.open(path)
    tree = t["events"]
    keys = tree.keys()
    e = tree.arrays([k for k in keys if not k.startswith("jet_")], library="np")
    jets = tree.arrays([k for k in keys if k.startswith("jet_")], library="ak") if "jet_pt" in keys else None
    meta = t["meta"].arrays(library="np")
    info = {k: meta[k].sum() for k in ["n_read", "n_selected", "sumw_read"]}
    for k in ["jvt_wp", "vbs_mjj_cut", "vbs_deta_cut", "jet_time_dr", "jet_time_dist_cut"]:
        u = np.unique(meta[k])
        if len(u) != 1: raise SystemExit(f"{path}: meta {k} differs between merged jobs: {u}")
        info[k] = u[0]
    info["sample"] = str(meta["sample"][0])
    info["njobs"] = len(meta["n_read"])
    exp = EXPECTED_COMPLETE.get(info["sample"])
    if exp:
        okn = info["n_read"] == exp["n_read"]
        oks = "sumw" not in exp or abs(info["sumw_read"] / exp["sumw"] - 1) < 1e-6
        print(f"[{what}] completeness: {info['n_read']:,} events read (complete sample {exp['n_read']:,})"
              + (f", sum w {info['sumw_read']:.10g} (raw {exp['sumw']:.10g})" if "sumw" in exp else "")
              + ("  OK" if okn and oks else "  WARNING: INCOMPLETE -- a shard is missing or doubled"))
    if recoil_min > 0:
        if np.all(e[recoil_col] < 0):
            raise SystemExit(f"{path}: {recoil_col} is unavailable (all -1) but the {what} cut asks for it"
                             + (" -- run VBF on the RAW ntuples (--ntuple-dir=.../highstats_vbf/)"
                                if recoil_col == "met_truth" else ""))
        keep = e[recoil_col] > recoil_min
    else:
        keep = np.ones(len(e["weight"]), bool)
    if C.ACCEPTANCE[args.acceptance] is not None:
        col = C.acc_column(args.acc_count)
        if "a_" + col not in e:
            raise SystemExit(f"{path}: no {col} columns -- rerun util/vbs_time_veto for --acceptance")
        acc = C.acceptance_mask(e["a_" + col], e["b_" + col], args.acceptance, args.acc_min_trk)
        info["acc_frac"] = (e["weight"][keep & acc].sum() / e["weight"][keep].sum(), int((keep & acc).sum()),
                            int(keep.sum()))
        keep = keep & acc
    ev = {k: v[keep] for k, v in e.items()}
    ev["jets"] = jets[keep] if jets is not None else None
    # the calibration uses EVERY selected event's HS legs, not only the recoil-passing ones
    return ev, e, info


def band_of(fA, fB, hA, hB, pA, pB):
    """vbs_region_stack.py's seven-band ladder; -1 where no band applies"""
    r1 = fA & fB & ((hA & pB) | (hB & pA))
    r2 = (fA & pA & ~fB & hB) | (fB & pB & ~fA & hA)
    hspu = ((hA & pB) | (hB & pA)) & ~r1 & ~r2
    hshs = hA & hB
    ppff = pA & pB & fA & fB
    ppo = pA & pB & ~ppff
    nei = (~hA & ~pA) | (~hB & ~pB)
    b = np.full(len(fA), -1)
    for i, m in enumerate([r1, r2, hspu, hshs, ppff, ppo, nei]):
        b[(b < 0) & m] = i
    return b


def bands(e):
    b = band_of(np.abs(e["a_eta"]) >= SPLIT, np.abs(e["b_eta"]) >= SPLIT,
                e["a_hs"], e["b_hs"], e["a_pu"], e["b_pu"])
    assert (b >= 0).all(), "band ladder not exhaustive"
    return b


def decisions(e, calib, nsigma, min_trk):
    """event vetoes: kept[method] -> bool per event, and which legs were timed (core).
    The literal all-ghost estimator (jj_all) gets the global factor: its n is
    ~23 tracks of mostly pileup, where the core's f(n) means nothing."""
    out = {}
    for est, key in [("core", "jj"), ("all", "jj_all")]:
        fA = calib(e[f"a_{est}_n"]) if est == "core" else calib.glob
        fB = calib(e[f"b_{est}_n"]) if est == "core" else calib.glob
        tA, sA, hA = e[f"a_{est}_t"], fA * e[f"a_{est}_sig"], e[f"a_{est}_n"] >= min_trk
        tB, sB, hB = e[f"b_{est}_t"], fB * e[f"b_{est}_sig"], e[f"b_{est}_n"] >= min_trk
        both = hA & hB
        den = np.sqrt(np.where(both, sA**2 + sB**2, 1.0))
        out[key] = ~(both & (np.abs(tA - tB) / den >= nsigma))
    tA, sA, hA = e["a_core_t"], calib(e["a_core_n"]) * e["a_core_sig"], e["a_core_n"] >= min_trk
    tB, sB, hB = e["b_core_t"], calib(e["b_core_n"]) * e["b_core_sig"], e["b_core_n"] >= min_trk
    for s in T0_SRC:
        ok, t0, s0 = e[f"ok_{s}"], e[f"t0_{s}"], e[f"infl_{s}"] * e[f"sig0_{s}"]
        rej = np.zeros(len(t0), bool)
        for t, sg, h in [(tA, sA, hA), (tB, sB, hB)]:
            rej |= ok & h & (np.abs(t - t0) / np.sqrt(sg**2 + s0**2) >= nsigma)
        out["t0_" + s] = ~rej
    return out, hA, hB


def take(arr, idx):
    """arr[event][idx[event]] as numpy; idx < 0 gives element 0 (callers mask those)"""
    j = ak.from_regular(ak.Array(np.where(idx >= 0, idx, 0).reshape(-1, 1)))
    return ak.to_numpy(ak.firsts(arr[j]))


def best_pair(J, jet_ok=None, pair_ok=None):
    """calcBestVbsPair over the jets with jet_ok: the max-m_jj opposite-hemisphere
    pair, first one on ties, among pairs pair_ok(i, k) allows. Returns numpy
    (mjj, i, k) with i, k local indices into the ORIGINAL jet list; -1 if none."""
    loc = ak.local_index(J.jet_pt)
    if jet_ok is not None: loc = loc[jet_ok]
    comb = ak.combinations(loc, 2)
    i, k = comb["0"], comb["1"]
    pt = ak.values_astype(J.jet_pt, np.float64)
    eta = ak.values_astype(J.jet_eta, np.float64)
    phi = ak.values_astype(J.jet_phi, np.float64)
    m2 = 2 * pt[i] * pt[k] * (np.cosh(eta[i] - eta[k]) - np.cos(phi[i] - phi[k]))
    ok = eta[i] * eta[k] < 0
    if pair_ok is not None: ok = ok & pair_ok(i, k)
    mjj = ak.where(ok, np.sqrt(ak.where(m2 > 0, m2, 0.0)), -1.0)
    b = ak.argmax(mjj, axis=1, keepdims=True)
    return (ak.to_numpy(ak.fill_none(ak.firsts(mjj[b]), -1.0)),
            ak.to_numpy(ak.fill_none(ak.firsts(i[b]), -1)),
            ak.to_numpy(ak.fill_none(ak.firsts(k[b]), -1)))


def repair(e, calib, nsigma, min_trk, source, cuts):
    """Jet removal + re-pairing (see the docstring). Returns (kept, final m_jj,
    final band) per event; final m_jj NaN and band -1 where not kept."""
    J = e["jets"]
    sig = calib.jagged(J.jet_n) * J.jet_sig
    timed = J.jet_n >= min_trk
    t = J.jet_t
    if source == "jj":
        jet_ok = None
        def pair_ok(i, k):
            both = timed[i] & timed[k]
            pull = np.abs(t[i] - t[k]) / np.sqrt(sig[i] ** 2 + sig[k] ** 2)
            return ~(both & (pull >= nsigma))
    else:
        s0 = e[f"infl_{source}"] * e[f"sig0_{source}"]
        pull = np.abs(t - e[f"t0_{source}"]) / np.sqrt(sig ** 2 + s0 ** 2)
        jet_ok = ~(timed & e[f"ok_{source}"] & (pull >= nsigma))
        pair_ok = None
    bm, bi, bk = best_pair(J, jet_ok, pair_ok)
    keepj = jet_ok if jet_ok is not None else ak.ones_like(timed)
    aeta = np.abs(J.jet_eta)
    nk = ak.to_numpy(ak.sum(keepj, axis=1))
    nf = ak.to_numpy(ak.sum(keepj & (aeta > MIN_ABS_ETA_JET) & (aeta < MAX_ABS_ETA_JET), axis=1))
    # the original pair, found again: its stored double-precision m_jj is kept, so
    # an event nothing was removed from is selected exactly as before
    same = (bi == e["la"]) & (bk == e["lb"])
    mjj = np.where(same, e["mjj"], bm)
    eta_i, eta_k = take(J.jet_eta, bi), take(J.jet_eta, bk)
    deta = np.abs(eta_i - eta_k)
    kept = (nk >= MIN_PASSPT_JETS) & (bi >= 0) & (mjj >= cuts[0]) & (deta >= cuts[1])
    if args.repair_reapply_fwd: kept &= nf >= 1
    fi, fk = np.abs(eta_i) >= SPLIT, np.abs(eta_k) >= SPLIT
    hi, hk = take(J.jet_hs, bi) > 0, take(J.jet_hs, bk) > 0
    pi, pk = take(J.jet_pu, bi) > 0, take(J.jet_pu, bk) > 0
    band = np.where(kept, band_of(fi, fk, hi, hk, pi, pk), -1)
    return kept, np.where(kept, mjj, np.nan), band


def binned(e, kept, edges):
    """per NO-VETO m_jj column: (e, lo, hi, n MC events)"""
    idx, w = mjj_index(e["mjj"], edges), e["weight"].astype(float)
    res = []
    for b in range(len(edges) - 1):
        m = idx == b
        v, _, lo, hi = eff(w[m], kept[m])
        res.append((v, lo, hi, int(m.sum())))
    return res


loaded = {}
for role, path, col, cut, norm in [("sig", args.sig, "met_truth", args.met_min, args.sig_norm),
                                   ("bkg", args.bkg, "z_pt", args.zpt_min, args.bkg_norm)]:
    if path: loaded[role] = (path, col, cut, norm) + load(path, col, cut, role)
HAS_JETS = all(v[4]["jets"] is not None for v in loaded.values())
ALL_METHODS = VETO_METHODS + (REPAIR_METHODS if HAS_JETS else [])
for m in DRAWN + [f"rp_{s}" for s in COMPARE]:
    if m.startswith("rp_") and not HAS_JETS:
        raise SystemExit(f"{m} needs the jet_* arrays, which these files predate -- rerun util/vbs_time_veto")
# one calibration for every file unless asked otherwise (see the docstring)
own = {r: measure(v[5], True, f"{r} HS legs") for r, v in loaded.items()}   # each file's own, for the record
mode = args.jet_infl
if mode in ("sig-ntrk", "sig"):
    src = "sig" if "sig" in loaded else "bkg"
    c = measure(loaded[src][5], mode == "sig-ntrk", f"the {src} file's HS legs")
    calibs = {r: c for r in loaded}
elif mode in ("per-file", "per-file-ntrk"):
    calibs = {r: measure(loaded[r][5], mode == "per-file-ntrk", f"the {r} file's own HS legs") for r in loaded}
else:
    calibs = {r: Calib(float(mode)) for r in loaded}


def outcomes(s, nsigma, methods):
    """kept / final m_jj / final band for each method at a threshold"""
    e, cal = s["e"], s["calib"]
    kept, hA, hB = decisions(e, cal, nsigma, args.min_trk)
    out = {}
    for m in methods:
        if m in kept:
            out[m] = (kept[m], np.where(kept[m], e["mjj"], np.nan), np.where(kept[m], s["band"], -1))
        else:
            out[m] = repair(e, cal, nsigma, args.min_trk, m[3:], s["cuts"])
    return out, hA, hB


samples = {}
for role, (path, col, cut, norm, ev, e_all, info) in loaded.items():
    cal = calibs[role]
    s = dict(e=ev, info=info, calib=cal, own=own[role], band=bands(ev), path=path, norm=norm, cutcol=col, cut=cut,
             cuts=(info["vbs_mjj_cut"], info["vbs_deta_cut"]))
    if HAS_JETS:
        J = ev["jets"]
        ev["la"] = ak.to_numpy(ak.argmax(J.jet_idx == ev["a_idx"], axis=1))
        ev["lb"] = ak.to_numpy(ak.argmax(J.jet_idx == ev["b_idx"], axis=1))
        # the stored pair must be the one calcBestVbsPair finds in the arrays (the
        # arrays are float, the C++ pairing double: only an exact m_jj tie could differ)
        bm, bi, bk = best_pair(J)
        bad = int(np.sum((bi != ev["la"]) | (bk != ev["lb"])))
        print(f"[{role}] pair re-formed from the jet arrays == stored legs in {len(bi) - bad:,} of {len(bi):,} events"
              + ("" if bad == 0 else "  WARNING: check the jet_* arrays"))
    res, hA, hB = outcomes(s, args.nsigma, ALL_METHODS)
    s["kept"] = {m: r[0] for m, r in res.items()}
    s["mjjf"] = {m: r[1] for m, r in res.items()}
    s["bandf"] = {m: r[2] for m, r in res.items()}
    s["mjjf"]["none"], s["kept"]["none"], s["bandf"]["none"] = ev["mjj"], np.ones(len(ev["mjj"]), bool), s["band"]
    s["hA"], s["hB"] = hA, hB
    samples[role] = s
    w = ev["weight"].astype(float)
    print(f"[{role}] {path}: {info['njobs']} job(s), {info['n_read']:,} read, {info['n_selected']:,} selected, "
          f"{len(w):,} after {col} > {cut:g}; jet-time sigma {cal.md()}; "
          f"this file's own HS-leg width {own[role].glob:.3f}"
          + ("; jet arrays present: re-pairing available" if HAS_JETS else ""))

# ── ROOT drawing ─────────────────────────────────────────────────────────────
ROOT = C.setup_root(args.out)


def sel_line(s):
    parts = [f"m_{{jj}} > {s['info']['vbs_mjj_cut']:g} GeV"]
    parts.append("no |#Delta#eta| cut" if s["info"]["vbs_deta_cut"] <= 0 else f"|#Delta#eta| > {s['info']['vbs_deta_cut']:g}")
    if s["cut"] > 0:
        parts.append(("E_{T}^{miss,truth}" if s["cutcol"] == "met_truth" else "p_{T}^{ll}") + f" > {s['cut']:g} GeV")
    if args.tag: parts.append(args.tag)
    return ",  ".join(parts)


def veto_line(methods=None):
    methods = methods or DRAWN
    rp = [m for m in methods if m.startswith("rp_")]
    what = ("jets failing it are removed and the pair re-formed" if len(rp) == len(methods)
            else "the event is vetoed if a pair jet fails it" if not rp
            else "filled: event veto;  open: jets removed, pair re-formed")
    return f"timing test |#Deltat| / #sigma #geq {args.nsigma:g}: {what}"


def jet_line():
    return (f"jet time: timed ghost tracks, #DeltaR < {JDR:g}, {JCUT:g}#sigma time clustering, "
            f"highest-#Sigmap_{{T}} cluster")


def calib_line(cal):
    if cal.fac is None:
        return f"#sigma_{{jet}} #times {cal.glob:.2f}" + ("" if cal.source == "fixed" else f"  (from {cal.source})")
    src = "signal" if "sig" in cal.source else "background"
    return (f"#sigma_{{jet}} #times {cal.fac[0]:.2f} (1 track) #rightarrow {cal.fac[-1]:.2f} (#geq {NBINS[-2]} tracks), "
            f"measured on {src} HS jets, applied to S and B alike")


# the bottom panel's MC count, split by the talk's cases (>= 1 timed HGTD track per jet)
CASE_PARTS = [("Forward-Forward", "#5B8DD6", "both"), ("Central-Forward", "#E8A33D", "one"), ("Other", "#BDBDBD", "neither")]


def case_counts(e):
    """per no-veto m_jj column: MC events of each case, in CASE_PARTS order"""
    if "a_hgtd_ntimed" not in e: return None
    idx = mjj_index(e["mjj"], EDGES)
    out = []
    for lab, col, mode in CASE_PARTS:
        m = C.acceptance_mask(e["a_hgtd_ntimed"], e["b_hgtd_ntimed"], mode, 1)
        out.append((lab, col, [int(np.sum(m & (idx == b))) for b in range(len(EDGES) - 1)]))
    return out


def column_figure(stem, ylab, series, legend_extra, label, lines, counts, count_title, ref=1.0, methods=None,
                  xnote="last bin includes overflow", pct=False, yrange=None, parts=None):
    """vbs_time_veto_common.column_figure on this script's columns; `methods` defaults to DRAWN"""
    acc = C.acceptance_line(args.acceptance, args.acc_min_trk, args.acc_count)
    if args.compact:
        lines = ["**" + C.CASE_TEXT[args.acceptance]]
    elif acc:
        lines = [lines[0], acc] + lines[1:]
    C.column_figure(stem, ylab, series, legend_extra, label, lines, counts, count_title, EDGES, STYLE,
                    methods or DRAWN, ref=ref, xnote=xnote, pct=pct, yrange=yrange, compact=args.compact,
                    count_parts=parts)


def fmt_pct(v): return f"  ({100 * v:.1f}%)"


any_s = next(iter(samples.values()))
JDR, JCUT = any_s["info"]["jet_time_dr"], any_s["info"]["jet_time_dist_cut"]
ANY_RP = any(m.startswith("rp_") for m in DRAWN)
EFF_SERIES = {role: {m: binned(samples[role]["e"], samples[role]["kept"][m], EDGES) for m in DRAWN} for role in samples}
EFF_RANGE = C.y_range(list(EFF_SERIES.values()), DRAWN)   # signal and background shown side by side: one axis
for role, ylab, label in [("sig", "Signal efficiency", args.sig_label), ("bkg", "Background efficiency", args.bkg_label)]:
    if role not in samples: continue
    s = samples[role]
    w = s["e"]["weight"].astype(float)
    series = EFF_SERIES[role]
    extra = {m: fmt_pct(eff(w, s["kept"][m])[0]) for m in DRAWN}
    counts = [n for _, _, _, n in series[DRAWN[0]]]
    column_figure(f"{args.out}_{role}_eff_mjj", ylab, series, extra, label,
                  [sel_line(s), veto_line(), jet_line(), calib_line(s["calib"]),
                   "all m_{jj} columns in brackets; generator-weighted"
                   + ("; events binned in their no-veto m_{jj}" if ANY_RP else "")],
                  counts, "MC events", yrange=EFF_RANGE, parts=case_counts(s["e"]))

# ── Yields and the summary ───────────────────────────────────────────────────
def k_factor(s):
    if s["norm"] == "none": return None
    n = NORM[s["norm"]]
    # the local sample is not the grid production: normalise it with its OWN sum of
    # weights (every event of it is read), assuming the same process cross-section
    sumw = s["info"]["sumw_read"] if s["info"]["sample"] == "local" else n["sumw"]
    return args.lumi * 1e3 * n["xs_pb"] / sumw * (args.bf_hinv if s["norm"] == "vbf" else 1.0)


md = [f"# VBS timing tests: {os.path.basename(args.out)}", "",
      f"timing test at {args.nsigma:g} sigma; a jet is timed with >= {args.min_trk} track(s); "
      f"yields at {args.lumi:g} fb^-1", "",
      "Methods: `jj` / `t0_<src>` veto the EVENT (pair unchanged); `rp_<src>` REMOVE the jets that fail "
      "the t0 test and re-form the pair from the rest (event lost only if no pair passes the selection); "
      "`rp_jj` keeps every jet and takes the max-m_jj time-compatible pair. After removal the "
      "forward-jet preselection is " + ("re-applied" if args.repair_reapply_fwd else "NOT re-applied (--repair-reapply-fwd)")
      + ".", ""]
for role, s in samples.items():
    i = s["info"]
    md.append(f"- **{role}**: `{s['path']}` -- {i['njobs']} job(s), {i['n_read']:,} events read "
              f"(sum w {i['sumw_read']:.6g}), {i['n_selected']:,} selected, {len(s['e']['weight']):,} after "
              f"{s['cutcol']} > {s['cut']:g}; jet-time sigma {s['calib'].md()} (this file's own HS-leg width "
              f"{s['own'].glob:.3f}); "
              f"normalisation: {NORM[s['norm']]['note'] if s['norm'] != 'none' else 'none'}"
              + (f"; scaled to B(H->inv) = {100 * args.bf_hinv:g}%" if s["norm"] == "vbf" and args.bf_hinv != 1.0 else "")
              + (" (LOCAL sample: its own sum w, same cross-section assumed)" if i["sample"] == "local" and s["norm"] != "none" else ""))
    if "acc_frac" in i:
        f, n, n0 = i["acc_frac"]
        md.append(f"  - **acceptance `{args.acceptance}`** (>= {args.acc_min_trk} HGTD track(s) per leg, {args.acc_count}): "
                  f"{n:,} of {n0:,} MC events after the recoil cut, {100 * f:.1f}% of their weight")
md.append("")
both = len(samples) == 2
hdr = "| method | " + " | ".join(("eff. " + r) for r in samples) + (" | S | B | S/sqrt(B) | vs no veto |" if both else " |")
md += [hdr, "|" + "---|" * (hdr.count("|") - 1)]
rows = {}
for m in ["none"] + ALL_METHODS:
    cells, Y = [], {}
    for role, s in samples.items():
        w = s["e"]["weight"].astype(float)
        kept = s["kept"][m]
        v, err, _, _ = eff(w, kept)
        cells.append(f"{100 * v:.1f} +- {100 * err:.1f}%")
        k = k_factor(s)
        Y[role] = (k * w[kept].sum(), k * np.sqrt((w[kept] ** 2).sum())) if k else (np.nan, np.nan)
    rows[m] = Y
    line = f"| {m} | " + " | ".join(cells)
    if both:
        S, B = Y["sig"][0], Y["bkg"][0]
        z0 = rows["none"]["sig"][0] / np.sqrt(rows["none"]["bkg"][0])
        line += f" | {S:,.0f} | {B:,.0f} +- {Y['bkg'][1]:,.0f} | {S / np.sqrt(B):.1f} | {S / np.sqrt(B) / z0:.3f} |"
    else:
        line += " |"
    md.append(line)
md.append("")
for role, s in samples.items():
    w, B = s["e"]["weight"].astype(float), s["band"]
    md += [f"**{role}: efficiency (%) by NO-VETO pair composition** (share = weighted fraction of the events)", "",
           "| composition | share | MC events | " + " | ".join(ALL_METHODS) + " |",
           "|" + "---|" * (3 + len(ALL_METHODS))]
    for i, name in enumerate(BANDS_MD):
        mm = B == i
        if w[mm].sum() <= 0: continue
        md.append(f"| {name} | {100 * w[mm].sum() / w.sum():.1f}% | {int(mm.sum()):,} | "
                  + " | ".join(f"{100 * eff(w[mm], s['kept'][m][mm])[0]:.1f}" for m in ALL_METHODS) + " |")
    tA, tB = s["hA"], s["hB"]
    md += ["", f"timed legs (core, >= {args.min_trk} track): both {100 * w[tA & tB].sum() / w.sum():.1f}%, "
           f"one {100 * w[tA ^ tB].sum() / w.sum():.1f}%, neither {100 * w[~tA & ~tB].sum() / w.sum():.1f}% (weighted)", ""]
    if HAS_JETS:
        rp = [m for m in ALL_METHODS if m.startswith("rp_")]
        md += [f"**{role}: composition (%) of the pair each surviving event ends up with** "
               "(weighted share of the no-veto yield; columns sum to that method's efficiency)", "",
               "| final composition | none | " + " | ".join(rp) + " |", "|" + "---|" * (2 + len(rp))]
        for i, name in enumerate(BANDS_MD):
            md.append(f"| {name} | " + " | ".join(f"{100 * w[s['bandf'][m] == i].sum() / w.sum():.1f}"
                                                   for m in ["none"] + rp) + " |")
        md.append("")
if both:
    # the analysis cuts at m_jj > X rather than using columns: cumulative thresholds, on
    # each event's FINAL m_jj (re-paired events move), against the no-veto yield above X
    sg_, bk_ = samples["sig"], samples["bkg"]
    ws_, wb_ = sg_["e"]["weight"].astype(float), bk_["e"]["weight"].astype(float)
    md += ["**S/sqrt(B) relative to no veto, by m_jj threshold** (final m_jj; no-veto B MC events in brackets)", "",
           "| m_jj > | B MC | " + " | ".join(ALL_METHODS) + " |", "|" + "---|" * (2 + len(ALL_METHODS))]
    for x in EDGES[:-1]:
        s0, b0 = ws_[sg_["e"]["mjj"] > x].sum(), wb_[bk_["e"]["mjj"] > x].sum()
        if b0 <= 0: continue
        cells = []
        for m in ALL_METHODS:
            sm = ws_[sg_["kept"][m] & (sg_["mjjf"][m] > x)].sum()
            bm_ = wb_[bk_["kept"][m] & (bk_["mjjf"][m] > x)].sum()
            cells.append(f"{(sm / s0) / np.sqrt(bm_ / b0):.3f}" if bm_ > 0 else "inf")
        md.append(f"| {x:g} GeV | {int((bk_['e']['mjj'] > x).sum())} | " + " | ".join(cells) + " |")
    md.append("")
open(args.out + "_summary.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
print("wrote", args.out + "_summary.md")

if not both:
    raise SystemExit(0)
sg, bk = samples["sig"], samples["bkg"]
ws, wb = sg["e"]["weight"].astype(float), bk["e"]["weight"].astype(float)
ks, kb = k_factor(sg), k_factor(bk)


def column_yields(s, w, k, m):
    """per FINAL m_jj column: (yield, MC-stat error)"""
    idx = mjj_index(s["mjjf"][m], EDGES)
    out = []
    for b in range(len(EDGES) - 1):
        mm = (idx == b) & s["kept"][m]
        out.append((k * w[mm].sum(), k * np.sqrt((w[mm] ** 2).sum())))
    return out


def label_lines(methods):
    return [f"S: {sel_line(sg)}",
            f"B: {sel_line(bk)}",
            veto_line(methods),
            jet_line(),
            calib_line(sg["calib"]) if sg["calib"] is bk["calib"]
            else f"#sigma_{{jet}}: {sg['calib'].short()} (S) / {bk['calib'].short()} (B)"]


XNOTE = "final pair; last bin includes overflow" if ANY_RP else "last bin includes overflow"
B_COUNTS = [int(np.sum(mjj_index(bk["e"]["mjj"], EDGES) == b)) for b in range(len(EDGES) - 1)]

# ── S / sqrt(B) and S / B per m_jj column relative to no veto, as a change in % ──
# (S/S0) / f(B/B0) per column, f = sqrt for S/sqrt(B) and the identity for S/B. For
# event vetoes the column content does not move and this is exactly eps_S / f(eps_B),
# with the efficiencies' binomial errors; re-paired events move between columns, so
# for those the ratio comes from the yields, with their MC-stat errors (conservative:
# numerator and denominator share events). The figures show 100 x (ratio - 1).
def relative(bfun, bpow):
    out = {}
    for m in DRAWN:
        out[m] = []
        if not m.startswith("rp_"):
            for (es, sl, sh, _), (eb, bl, bh, _) in zip(binned(sg["e"], sg["kept"][m], EDGES),
                                                       binned(bk["e"], bk["kept"][m], EDGES)):
                if not (np.isfinite(es) and np.isfinite(eb)) or eb <= 0:
                    out[m].append((np.nan, np.nan, np.nan, 0)); continue
                r = es / bfun(eb)
                err = r * np.hypot(0.5 * (sl + sh) / es, bpow * 0.5 * (bl + bh) / eb)
                out[m].append((r, err, err, 0))
        else:
            for (S, dS), (B, dB), (S0, _), (B0, _) in zip(column_yields(sg, ws, 1.0, m), column_yields(bk, wb, 1.0, m),
                                                        column_yields(sg, ws, 1.0, "none"), column_yields(bk, wb, 1.0, "none")):
                if S <= 0 or B <= 0 or S0 <= 0 or B0 <= 0:
                    out[m].append((np.nan, np.nan, np.nan, 0)); continue
                r = (S / S0) / bfun(B / B0)
                err = r * np.hypot(dS / S, bpow * dB / B)
                out[m].append((r, err, err, 0))
    return out


ratio = relative(np.sqrt, 0.5)
change = {m: [(100 * (r - 1), 100 * lo, 100 * hi, n) for r, lo, hi, n in ratio[m]] for m in DRAWN}
extra = {m: f"  ({100 * (eff(ws, sg['kept'][m])[0] / np.sqrt(eff(wb, bk['kept'][m])[0]) - 1):+.1f}%)" for m in DRAWN}
column_figure(f"{args.out}_s_over_sqrtb", "S/#sqrt{B} change vs no veto [%]", change, extra,
              f"S: {args.sig_label.split(', ')[-1]},  B: {args.bkg_label.split(', ')[-1]}",
              label_lines(DRAWN) + ["all m_{jj} columns in brackets; > 0 = the timing test raises S/#sqrt{B}"],
              B_COUNTS, "B MC events", ref=0.0, pct=True, xnote=XNOTE, parts=case_counts(bk["e"]))

# ── S / sqrt(B) from the YIELDS per m_jj column (absolute, at --lumi) ───────
# The same information as the change above -- which is exactly this divided by
# its no-veto value, minus 1, in % -- but on the absolute scale, which carries the two
# normalisation assumptions (B(H->inv), Z->ll standing in for Z->nunu) that
# cancel in the ratio. Errors: MC statistics of both yields.
if ks and kb:
    SB_METHODS = ["none"] + DRAWN
    absol, extra = {}, {}
    for m in SB_METHODS:
        absol[m] = []
        for (S, dS), (B, dB) in zip(column_yields(sg, ws, ks, m), column_yields(bk, wb, kb, m)):
            if S <= 0 or B <= 0:
                absol[m].append((np.nan, np.nan, np.nan, 0)); continue
            z = S / np.sqrt(B)
            err = z * np.hypot(dS / S, 0.5 * dB / B)
            absol[m].append((z, err, err, 0))
        extra[m] = f"  ({ks * ws[sg['kept'][m]].sum() / np.sqrt(kb * wb[bk['kept'][m]].sum()):.0f})"
    column_figure(f"{args.out}_s_over_sqrtb_yields", "S / #sqrt{B}  per m_{jj} column", absol, extra,
                  f"S: {args.sig_label.split(', ')[-1]},  B: {args.bkg_label.split(', ')[-1]},  {args.lumi:g} fb^{{-1}}",
                  label_lines(DRAWN) +
                  [f"S at B(H#rightarrowinv) = {100 * args.bf_hinv:g}%;  B = Z(#rightarrowll)+jets alone (stand-in for "
                   f"Z(#rightarrow#nu#nu));  all columns in brackets"],
                  B_COUNTS, "B MC events", ref=None, methods=SB_METHODS, xnote=XNOTE, parts=case_counts(bk["e"]))


# ── S / B per m_jj column: the change vs no veto, and from the yields ───────
ratio_b = relative(lambda x: x, 1.0)
change_b = {m: [(100 * (r - 1), 100 * lo, 100 * hi, n) for r, lo, hi, n in ratio_b[m]] for m in DRAWN}
extra = {m: f"  ({100 * (eff(ws, sg['kept'][m])[0] / eff(wb, bk['kept'][m])[0] - 1):+.1f}%)" for m in DRAWN}
column_figure(f"{args.out}_s_over_b", "S/B change vs no veto [%]", change_b, extra,
              f"S: {args.sig_label.split(', ')[-1]},  B: {args.bkg_label.split(', ')[-1]}",
              label_lines(DRAWN) + ["all m_{jj} columns in brackets; > 0 = the timing test raises S/B"],
              B_COUNTS, "B MC events", ref=0.0, pct=True, xnote=XNOTE, parts=case_counts(bk["e"]))
if ks and kb:
    absol, extra = {}, {}
    for m in SB_METHODS:
        absol[m] = []
        for (S, dS), (B, dB) in zip(column_yields(sg, ws, ks, m), column_yields(bk, wb, kb, m)):
            if S <= 0 or B <= 0:
                absol[m].append((np.nan, np.nan, np.nan, 0)); continue
            z = S / B
            absol[m].append((z, z * np.hypot(dS / S, dB / B), z * np.hypot(dS / S, dB / B), 0))
        extra[m] = f"  ({ws[sg['kept'][m]].sum() * ks / (kb * wb[bk['kept'][m]].sum()):.3f})"
    column_figure(f"{args.out}_s_over_b_yields", "S / B  per m_{jj} column", absol, extra,
                  f"S: {args.sig_label.split(', ')[-1]},  B: {args.bkg_label.split(', ')[-1]}",
                  label_lines(DRAWN) +
                  [f"S at B(H#rightarrowinv) = {100 * args.bf_hinv:g}%;  B = Z(#rightarrowll)+jets alone (stand-in for "
                   f"Z(#rightarrow#nu#nu));  all columns in brackets"],
                  B_COUNTS, "B MC events", ref=None, methods=SB_METHODS, xnote=XNOTE, parts=case_counts(bk["e"]))


# ── How a jet gets a time: the HS-leg pulls against the true vertex time ─────
def jet_time_figures(stem, e_all, cal):
    """<stem>_jet_time_pulls: (t_jet - t_HS^truth) / sigma for the signal's paper-HS
    legs (reco tagging jets truth-matched to a hard-scatter jet: dR < 0.3 to a truth
    HS jet above 10 GeV), with the quoted sigma and with sigma x f(n);
    <stem>_jet_time_calib: f(n)"""
    ROOT_ = C.ROOT
    p, n = C.hs_pulls(e_all)
    pc = p / cal(n)
    c = ROOT_.TCanvas("c_jtp", "", 800, 650); C.KEEP.append(c)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.05); c.SetLogy(True)
    hs = []
    drawn = [("raw", p, "#9E9E9E", True)] + ([("cal", pc, "kP10Blue", False)] if args.jet_pulls == "both" else [])
    for name, v, col, fill in drawn:
        h = ROOT_.TH1D(f"h_jtp_{name}", "", 64, -8, 8); C.KEEP.append(h)
        for x in np.clip(v, -7.99, 7.99): h.Fill(x)
        h.Scale(1.0 / h.Integral())
        h.SetLineColor(C.colour(col)); h.SetLineWidth(2)
        if fill: h.SetFillColor(C.colour(col)); h.SetFillStyle(3004)
        hs.append(h)
    hs[0].SetMinimum(2e-5); hs[0].SetMaximum(200.0)   # two decades of headroom: the label and legend sit above the peak
    hs[0].GetXaxis().SetTitle("(t_{jet} - t_{HS}^{truth}) / #sigma_{jet}   (overflow in the edge bins)")
    hs[0].GetYaxis().SetTitle("Fraction of truth-matched HS tagging jets")
    hs[0].Draw("HIST")
    for h in hs[1:]: h.Draw("HIST SAME")
    g = ROOT_.TF1("f_unit", "[0]*exp(-0.5*x*x)", -8, 8); C.KEEP.append(g)
    g.SetParameter(0, 0.25 / np.sqrt(2 * np.pi)); g.SetLineColor(ROOT_.kBlack); g.SetLineStyle(2); g.Draw("SAME")
    leg = ROOT_.TLegend(0.16, 0.735 - 0.045 * (len(hs) + 1), 0.60, 0.735); ROOT_.StyleLegend(leg, 0.03); C.KEEP.append(leg)
    leg.AddEntry(hs[0], f"quoted #sigma: {100 * np.mean(np.abs(p) >= args.nsigma):.0f}% beyond {args.nsigma:g}", "f")
    if len(hs) > 1:
        leg.AddEntry(hs[1], f"#sigma #times f(n): {100 * np.mean(np.abs(pc) >= args.nsigma):.1f}% beyond {args.nsigma:g}", "l")
    leg.AddEntry(g, "unit Gaussian", "l")
    leg.Draw()
    ROOT_.ATLASLabel(0.17, 0.88, "")
    t = ROOT_.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04); C.KEEP.append(t)
    t.DrawLatex(0.31, 0.88, "Simulation Internal")
    t2 = ROOT_.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.028); t2.SetTextColor(ROOT_.kGray + 2); C.KEEP.append(t2)
    t2.DrawLatex(0.17, 0.835, "VBF H#rightarrowinv. reco tagging jets truth-matched to a hard-scatter jet")
    t2.DrawLatex(0.17, 0.80, f"(#DeltaR < 0.3 to a truth HS jet with p_{{T}} > 10 GeV), with a time: {len(p):,} jets")
    t2.DrawLatex(0.17, 0.765, jet_line())
    c.Print(stem + "_jet_time_pulls.pdf"); c.Print(stem + "_jet_time_pulls.png")
    print("wrote", stem + "_jet_time_pulls.pdf / .png")
    # f(n), with each bin's share of the jets and its tail after the correction
    c2 = ROOT_.TCanvas("c_jtc", "", 800, 650); C.KEEP.append(c2)
    c2.SetLeftMargin(0.14); c2.SetRightMargin(0.05); c2.SetBottomMargin(0.14)
    labels = [f"{lo}" if hi == lo + 1 else (f"{lo}-{hi - 1}" if hi < 10**6 else f"#geq{lo}") for lo, hi in zip(C.NBINS[:-1], C.NBINS[1:])]
    fr = ROOT_.TH1D("fr_jtc", "", len(labels), 0, len(labels)); C.KEEP.append(fr)
    for i, lab in enumerate(labels): fr.GetXaxis().SetBinLabel(i + 1, lab)
    fr.SetMinimum(0.0); fr.SetMaximum(3.4)
    fr.GetXaxis().SetTitle("tracks in the jet time, n"); fr.GetXaxis().SetLabelSize(0.06); fr.GetXaxis().SetTitleOffset(1.1)
    fr.GetYaxis().SetTitle("f(n) = width of the pull / quoted #sigma")
    fr.Draw("AXIS")
    gr = ROOT_.TGraph(); C.KEEP.append(gr)
    tx = ROOT_.TLatex(); tx.SetTextFont(42); tx.SetTextSize(0.03); tx.SetTextAlign(21); C.KEEP.append(tx)
    for i, (lo, hi) in enumerate(zip(C.NBINS[:-1], C.NBINS[1:])):
        m = (n >= lo) & (n < hi)
        f = cal.fac[i] if cal.fac is not None else cal.glob
        gr.SetPoint(i, i + 0.5, f)
        tx.DrawLatex(i + 0.5, f + 0.12, f"{f:.2f}")
        tx.DrawLatex(i + 0.5, 0.25, f"{100 * m.mean():.0f}% of jets")
        tx.DrawLatex(i + 0.5, 0.10, f"{100 * np.mean(np.abs(pc[m]) >= args.nsigma):.1f}% > {args.nsigma:g}#sigma")
    gr.SetMarkerStyle(20); gr.SetMarkerSize(1.6); gr.SetMarkerColor(C.colour("kP10Blue")); gr.Draw("P SAME")
    one = ROOT_.TLine(0, 1, len(labels), 1); one.SetLineStyle(2); one.SetLineColor(ROOT_.kGray + 1); one.Draw(); C.KEEP.append(one)
    ROOT_.ATLASLabel(0.17, 0.88, "")
    t.DrawLatex(0.31, 0.88, "Simulation Internal")
    t2.DrawLatex(0.17, 0.83, "f(n) = 1.4826 #times median |pull|, hard-scatter tagging jets of VBF H#rightarrowinv.")
    t2.DrawLatex(0.17, 0.795, "applied to signal and background alike")
    c2.Print(stem + "_jet_time_calib.pdf"); c2.Print(stem + "_jet_time_calib.png")
    print("wrote", stem + "_jet_time_calib.pdf / .png")


jet_time_figures(args.out, loaded["sig"][5], calibs["sig"])


# ── The other side of the test: the event t0 against the true vertex time ────
def t0_figure(stem, e_all, what):
    """<stem>_t0_pulls: (t0 - t_HS^truth) / (infl x sigma_t0) for every selected event
    of the file with that t0, one line per reco t0 -- the t0 half of the jet-vs-t0
    pull, with the sigma the test uses (inflation included)"""
    ROOT_ = C.ROOT
    c = ROOT_.TCanvas("c_t0p_" + what, "", 800, 650); C.KEEP.append(c)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.05); c.SetLogy(True)
    leg = ROOT_.TLegend(0.16, 0.565, 0.93, 0.735); ROOT_.StyleLegend(leg, 0.027); C.KEEP.append(leg)
    first = True
    nall = len(e_all["t0_truth"])
    for k, src in enumerate(["trkptz", "waves", "hgtd"]):
        ok = e_all[f"ok_{src}"].astype(bool)
        sg = e_all[f"infl_{src}"][ok] * e_all[f"sig0_{src}"][ok]
        good = sg > 0
        pull = (e_all[f"t0_{src}"][ok][good] - e_all["t0_truth"][ok][good]) / sg[good]
        h = ROOT_.TH1D(f"h_t0p_{what}_{src}", "", 64, -8, 8); C.KEEP.append(h)
        for x in np.clip(pull, -7.99, 7.99): h.Fill(x)
        h.Scale(1.0 / h.Integral())
        lab, col, _ = STYLE["t0_" + src]
        h.SetLineColor(C.colour(col)); h.SetLineWidth(2); h.SetLineStyle(1 if src != "hgtd" else 7)
        if first:
            h.SetMinimum(2e-5); h.SetMaximum(200.0)
            h.GetXaxis().SetTitle("(t_{0} - t_{HS}^{truth}) / #sigma_{t_{0}}   (overflow in the edge bins)")
            h.GetYaxis().SetTitle("Fraction of events with that t_{0}")
            h.Draw("HIST"); first = False
        else:
            h.Draw("HIST SAME")
        name = lab.replace("Jet vs t_{0}: ", "")
        leg.AddEntry(h, f"{name}: {100 * np.mean(np.abs(pull) >= args.nsigma):.1f}% beyond {args.nsigma:g};  "
                     f"#sigma {np.median(sg[good]):.0f} ps (median);  valid in {100 * ok.sum() / nall:.0f}%", "l")
    g = ROOT_.TF1("f_unit_t0_" + what, "[0]*exp(-0.5*x*x)", -8, 8); C.KEEP.append(g)
    g.SetParameter(0, 0.25 / np.sqrt(2 * np.pi)); g.SetLineColor(ROOT_.kBlack); g.SetLineStyle(2); g.Draw("SAME")
    leg.AddEntry(g, "unit Gaussian", "l")
    leg.Draw()
    ROOT_.ATLASLabel(0.17, 0.88, "")
    t = ROOT_.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04); C.KEEP.append(t)
    t.DrawLatex(0.31, 0.88, "Simulation Internal")
    t2 = ROOT_.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.028); t2.SetTextColor(ROOT_.kGray + 2); C.KEEP.append(t2)
    t2.DrawLatex(0.17, 0.835, f"{'VBF H#rightarrowinv.' if what == 'sig' else 'Z(#rightarrowll)+jets'}: every selected event, {nall:,}")
    t2.DrawLatex(0.17, 0.80, "#sigma: the vertex time's quoted #sigma #times its per-sample inflation, as in the test")
    t2.DrawLatex(0.17, 0.765, "HGTD: Athena's RecoVtx_time, where valid")
    c.Print(f"{stem}_t0_pulls{'' if what == 'sig' else '_bkg'}.pdf"); c.Print(f"{stem}_t0_pulls{'' if what == 'sig' else '_bkg'}.png")
    print("wrote", f"{stem}_t0_pulls{'' if what == 'sig' else '_bkg'}.pdf / .png")


for r in loaded: t0_figure(args.out, loaded[r][5], r)

# ── The trade-off: eps_S against background rejection as the threshold moves ─
# (all m_jj columns). Background rejection = 1 / eps_B, the ROC convention.
SCAN = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0])


def tradeoff_curves(methods):
    """per method: (eps_S, eps_B) at each SCAN threshold"""
    curves = {m: [] for m in methods}
    for ns in SCAN:
        rs, rb = outcomes(sg, ns, methods)[0], outcomes(bk, ns, methods)[0]
        for m in methods:
            curves[m].append((eff(ws, rs[m][0])[0], eff(wb, rb[m][0])[0]))
    return {m: np.array(v) for m, v in curves.items()}


def tradeoff_figure(stem, curves, methods, dashed=()):
    """Signal efficiency against background rejection 1/eps_B, one point per
    scanned threshold, the --nsigma working point as a large marker."""
    xy = {m: (curves[m][:, 0], 1.0 / curves[m][:, 1]) for m in methods}
    C.roc_figure(stem, xy, methods, STYLE, int(np.argmin(np.abs(SCAN - args.nsigma))),
                 [f"S: {sel_line(sg)}", f"B: {sel_line(bk)}",
                  f"threshold scanned {SCAN[0]:g}-{SCAN[-1]:g}#sigma; large markers at {args.nsigma:g}#sigma"
                  + (";  dashed: jets removed, pair re-formed" if dashed else "")]
                 + ([C.acceptance_line(args.acceptance, args.acc_min_trk, args.acc_count)] if C.ACCEPTANCE[args.acceptance] is not None else []),
                 "Signal efficiency #varepsilon_{S}", "Background rejection 1/#varepsilon_{B}", dashed=dashed)


TRADEOFF = [(args.out + "_tradeoff", DRAWN, [m for m in DRAWN if m.startswith("rp_")])]
# ── Event veto against jet removal + re-pairing, same t0 sources ────────────
if COMPARE:
    cmp_methods = [f"t0_{s}" for s in COMPARE] + [f"rp_{s}" for s in COMPARE]
    TRADEOFF.append((args.out + "_veto_vs_repair_tradeoff", cmp_methods, [f"rp_{s}" for s in COMPARE]))
SCANNED = list(dict.fromkeys(m for _, ms, _ in TRADEOFF for m in ms))
CURVES = tradeoff_curves(SCANNED)
for stem, ms, dashed in TRADEOFF:
    tradeoff_figure(stem, CURVES, ms, dashed)
# the points behind the trade-off figures, appended to the summary
scan = ["**Threshold scan: signal efficiency (%) / background rejection 1/eps_B** (all m_jj columns)", "",
        "| method | " + " | ".join(f"{t:g} sigma" for t in SCAN) + " |", "|" + "---|" * (len(SCAN) + 1)]
scan += [f"| {m} | " + " | ".join(f"{100 * es:.1f} / {1 / eb:.3f}" for es, eb in CURVES[m]) + " |" for m in SCANNED]
open(args.out + "_summary.md", "a").write("\n".join(scan) + "\n")
print("\n".join(scan))
print("appended the threshold scan to", args.out + "_summary.md")
