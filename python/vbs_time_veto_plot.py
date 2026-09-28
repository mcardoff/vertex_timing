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
<out>_s_over_sqrtb.(pdf|png) (S/sqrt(B) relative to no veto, per column),
<out>_s_over_sqrtb_yields.(pdf|png) (S/sqrt(B) from the yields at --lumi),
<out>_tradeoff.(pdf|png) and, with --compare, <out>_veto_vs_repair_tradeoff,
plus <out>_summary.md.
"""
import argparse, os
import numpy as np
import uproot
import awkward as ak

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MJJ_EDGES = [200, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000]
SPLIT = 2.4   # forward = |eta| >= this, as vbs_region_stack.py
# the event preselection re-applied after jet removal: >= MIN_PASSPT_JETS jets,
# >= 1 of them in MIN_ABS_ETA_JET < |eta| < MAX_ABS_ETA_JET (clustering_constants.h)
MIN_PASSPT_JETS, MIN_ABS_ETA_JET, MAX_ABS_ETA_JET = 2, 2.38, 4.00

# Normalisation to --lumi (3000 fb^-1 by default): yield = w * L * sigma * eff / sum(w_raw).
# sigma x filter eff from AMI's 14 TeV evgen metadata (e8481), sum(w) over EVERY event of
# the RAW (pre-skim) sample: util/scratch/region_weights/sumw.txt, and the table in
# results/rpt_region_distributions.md. The raw sums are what the normalisation needs --
# a skimmed run's own meta sum is of the skim, not the sample.
NORM = {
    "vbf":   dict(xs_pb=4.4316 * 0.5148835, sumw=6497726.944462,
                  note="600026 VBF H->ZZ->4nu MET75, sigma_VBF only: B(H->inv) = 100%"),
    "zjets": dict(xs_pb=2067.2 + 2062.1, sumw=2067163337.673096 + 2061509864.110840,
                  note="601189 + 601190 Z->ee + Z->mumu"),
}

# What a COMPLETE run reads, so a missing shard is caught before it shows up as a
# yield: VBF must be read from the RAW ntuples (the truth record is needed for
# E_T^miss), so its meta sum of weights must equal the raw sum above; Z+jets reads
# the skims, whose event count after the vertex + Z->ll cuts is fixed.
EXPECTED_COMPLETE = {"vbf": dict(n_read=1464400, sumw=6497726.944462), "zjets": dict(n_read=530962)}

T0_SRC = ["trkptz", "waves", "hgtd", "tzp", "truth"]
# (legend, colour, marker): the region talk's case colours; jet-jet in black.
# Event vetoes filled, jet removal + re-pairing open in the same colour.
STYLE = {
    "jj":        ("Jet vs jet (no t_{0})",                     "#000000",    20),
    "jj_all":    ("Jet vs jet, all timed ghost tracks",        "#8C8C8C",    24),
    "t0_trkptz": ("Jet vs t_{0}: TRKPTZ",                      "kP10Red",    21),
    "t0_waves":  ("Jet vs t_{0}: WAVeS",                       "kP10Yellow", 22),
    "t0_hgtd":   ("Jet vs t_{0}: HGTD (Athena)",               "kP10Blue",   23),
    "t0_tzp":    ("Jet vs t_{0}: TZP",                         "kP10Green",  33),
    "t0_truth":  ("Jet vs t_{0}: truth (perfect t_{0})",       "kP10Brown",  34),
    "rp_jj":     ("Re-pair: time-compatible jets",             "#000000",    24),
    "rp_trkptz": ("Jet removal: TRKPTZ t_{0}",                 "kP10Red",    25),
    "rp_waves":  ("Jet removal: WAVeS t_{0}",                  "kP10Yellow", 26),
    "rp_hgtd":   ("Jet removal: HGTD t_{0} (Athena)",          "kP10Blue",   32),
    "rp_tzp":    ("Jet removal: TZP t_{0}",                    "kP10Green",  27),
    "rp_truth":  ("Jet removal: truth t_{0}",                  "kP10Brown",  28),
    "none":      ("No timing veto",                            "#7F7F7F",    25),
}
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
    ev = {k: v[keep] for k, v in e.items()}
    ev["jets"] = jets[keep] if jets is not None else None
    # the calibration uses EVERY selected event's HS legs, not only the recoil-passing ones
    return ev, e, info


# bins of n = tracks in the core jet-time estimate: the pull width falls from
# 2.3 (one track) to 1.3 (nine or more) on VBF, which a single factor averages over
NBINS = [1, 2, 3, 4, 6, 9, 10**6]


def mad(x):
    return 1.4826 * float(np.median(np.abs(x)))


def hs_pulls(e_all):
    """(t_jet - t_HS) / sigma_jet and n for every paper-HS leg with a core time"""
    p, n = [], []
    for L in "ab":
        m = e_all[f"{L}_hs"] & (e_all[f"{L}_core_n"] > 0)
        p.append((e_all[f"{L}_core_t"][m] - e_all["t0_truth"][m]) / e_all[f"{L}_core_sig"][m])
        n.append(e_all[f"{L}_core_n"][m])
    return np.concatenate(p), np.concatenate(n)


class Calib:
    """sigma_jet -> f(n) x sigma_jet; `fac` None means one global factor."""
    def __init__(self, glob, fac=None, source="fixed", nlegs=0):
        self.glob, self.fac, self.source, self.nlegs = glob, fac, source, nlegs

    def __call__(self, n):
        n = np.asarray(n)
        if self.fac is None: return np.full(len(n), self.glob)
        out = np.full(len(n), self.glob)
        for (lo, hi), f in zip(zip(NBINS[:-1], NBINS[1:]), self.fac):
            out[(n >= lo) & (n < hi)] = f
        return out

    def jagged(self, n):
        """the same, for a jagged (per-jet) array of n"""
        return ak.unflatten(self(ak.to_numpy(ak.flatten(n))), ak.num(n))

    def short(self):
        if self.fac is None: return f"{self.glob:.2f}"
        return f"f(n_{{trk}}): {self.fac[0]:.2f} (1 trk) - {self.fac[-1]:.2f} (#geq {NBINS[-2]})"

    def md(self):
        if self.fac is None: return f"x {self.glob:.3f} ({self.source})"
        bins = ", ".join(f"n={lo}" + ("" if hi == lo + 1 else (f"-{hi - 1}" if hi < 10**6 else "+")) + f": {f:.2f}"
                         for (lo, hi), f in zip(zip(NBINS[:-1], NBINS[1:]), self.fac))
        return f"x f(n) from {self.source} ({self.nlegs:,} legs; {bins}; global {self.glob:.3f})"


def measure(e_all, binned, source):
    p, n = hs_pulls(e_all)
    g = mad(p)
    if not binned: return Calib(g, None, source, len(p))
    fac = []
    for lo, hi in zip(NBINS[:-1], NBINS[1:]):
        m = (n >= lo) & (n < hi)
        fac.append(mad(p[m]) if m.sum() >= 20 else g)
    return Calib(g, fac, source, len(p))


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


def eff(w, kept):
    """Weighted efficiency, its weighted-binomial error, and a Wilson interval.

    The error is sqrt(sum_kept w^2 (1-e)^2 + sum_rejected w^2 e^2) / sum w. It
    vanishes at e = 0 or 1, which in a column of a few Z+jets events is a lie,
    so the drawn interval is Wilson's (1 sigma) at the effective event count
    that reproduces that variance -- or, at the boundary, at (sum w)^2 / sum w^2.
    Returns (e, err, lo, hi) with lo/hi the distances below/above e."""
    W = w.sum()
    if W <= 0: return np.nan, np.nan, np.nan, np.nan
    e = w[kept].sum() / W
    var = ((w[kept] ** 2).sum() * (1 - e) ** 2 + (w[~kept] ** 2).sum() * e ** 2) / W ** 2
    n = e * (1 - e) / var if (0 < e < 1 and var > 0) else W ** 2 / (w ** 2).sum()
    ec = min(max(e, 0.0), 1.0)
    cen = (ec + 0.5 / n) / (1 + 1 / n)
    half = np.sqrt(ec * (1 - ec) / n + 0.25 / n ** 2) / (1 + 1 / n)
    return e, np.sqrt(var), max(e - (cen - half), 0.0), max(cen + half - e, 0.0)


def mjj_index(mjj, edges):
    """column of each event; overflow into the top bin, NaN (not kept) -> -1"""
    idx = np.digitize(np.minimum(mjj, 0.5 * (edges[-2] + edges[-1])), edges) - 1
    return np.where(np.isfinite(mjj), idx, -1)


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
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare(f'#include "{REPO}/src/AtlasStyle.h"\n#include "{REPO}/src/AtlasLabels.h"')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
KEEP = []   # PyROOT objects must outlive the canvas print


def colour(c):
    return getattr(ROOT, c) if c.startswith("k") else ROOT.TColor.GetColor(c)


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


def graph(series, k, n_meth, key):
    """markers dodged inside each m_jj column so the methods do not overlap"""
    g = ROOT.TGraphAsymmErrors()
    for b, (v, lo, hi, _) in enumerate(series):
        if not np.isfinite(v): continue
        x0, x1 = EDGES[b], EDGES[b + 1]
        x = x0 + (x1 - x0) * (0.5 + 0.55 * ((k + 0.5) / n_meth - 0.5))
        i = g.GetN()
        g.SetPoint(i, x, v)
        g.SetPointError(i, 0, 0, lo, hi)
    lab, c, mk = STYLE[key]
    g.SetMarkerColor(colour(c)); g.SetLineColor(colour(c)); g.SetMarkerStyle(mk)
    g.SetMarkerSize(1.15 if key == "jj" else 1.0); g.SetLineWidth(2)
    KEEP.append(g)
    return g


def atlas_block(label, lines):
    ROOT.ATLASLabel(0.18, 0.895, "")
    al = ROOT.TLatex(); al.SetNDC(); al.SetTextFont(42); al.SetTextSize(0.045)
    al.DrawLatex(0.305, 0.895, "Simulation Internal")
    ROOT.ATLASEnergyLabel(0.18, 0.85, label)
    p = ROOT.TLatex(); p.SetNDC(); p.SetTextFont(42); p.SetTextSize(0.028); p.SetTextColor(ROOT.kGray + 2)
    for i, ln in enumerate(lines):
        p.DrawLatex(0.18, 0.81 - 0.036 * i, ln)
    KEEP.extend([al, p])
    return 0.81 - 0.036 * len(lines)


def column_figure(stem, ylab, series, legend_extra, label, lines, counts, count_title, ref=1.0, methods=None,
                  xnote="last bin includes overflow"):
    """The region plots' layout: per-m_jj-column values on top (one dodged marker
    per method), the MC event count of each column underneath on a log scale.
    `ref` draws a dashed reference line (None: none); `methods` defaults to DRAWN."""
    methods = methods or DRAWN
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 800, 780)
    pTop = ROOT.TPad("pTop", "", 0.0, 0.30, 1.0, 1.0)
    pBot = ROOT.TPad("pBot", "", 0.0, 0.00, 1.0, 0.30)
    pTop.SetBottomMargin(0.025); pTop.SetTopMargin(0.06)
    pBot.SetTopMargin(0.04); pBot.SetBottomMargin(0.40)
    pTop.Draw(); pBot.Draw(); KEEP.extend([c, pTop, pBot])
    pTop.cd()
    frame = ROOT.TH1D("frame_" + os.path.basename(stem), "", len(EDGES) - 1, EDGES); KEEP.append(frame)
    lows = [v - lo for m in methods for v, lo, _, _ in series[m] if np.isfinite(v)]
    highs = [v + hi for m in methods for v, _, hi, _ in series[m] if np.isfinite(v)]
    ytop = max([max(highs) if highs else 1.0] + ([ref] if ref is not None else []))
    if ref is not None:   # efficiency-like: round down to the next 0.1
        ymin = max(0.0, np.floor((min(lows) - 0.03) * 10) / 10) if lows else 0.0
    else:                 # absolute: a margin below the lowest point
        ymin = max(0.0, min(lows) - 0.06 * (ytop - min(lows))) if lows else 0.0
    frame.SetMinimum(ymin); frame.SetMaximum(ytop + (ytop - ymin) * (0.95 + 0.13 * max(0, len(lines) - 3)))
    frame.GetXaxis().SetLabelSize(0); frame.GetYaxis().SetTitle(ylab); frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw("AXIS")
    if ref is not None:
        one = ROOT.TLine(EDGES[0], ref, EDGES[-1], ref); one.SetLineStyle(2); one.SetLineColor(ROOT.kGray + 1); one.Draw()
        KEEP.append(one)
    ytxt = atlas_block(label, lines)
    nrow = (len(methods) + 1) // 2
    leg = ROOT.TLegend(0.17, ytxt - 0.055 * nrow, 0.93, ytxt - 0.005); ROOT.StyleLegend(leg); leg.SetNColumns(2); KEEP.append(leg)
    for k, m in enumerate(methods):
        g = graph(series[m], k, len(methods), m)
        g.Draw("P SAME")
        leg.AddEntry(g, STYLE[m][0] + legend_extra.get(m, ""), "p")
    leg.Draw()
    pTop.RedrawAxis()
    pBot.cd(); pBot.SetLogy(True)
    tot = ROOT.TH1D("tot_" + os.path.basename(stem), "", len(EDGES) - 1, EDGES); KEEP.append(tot)
    for b, n in enumerate(counts): tot.SetBinContent(b + 1, n)
    tot.SetFillColor(ROOT.TColor.GetColor("#DCE6F2")); tot.SetLineColor(ROOT.kBlack); tot.SetLineWidth(1)
    tot.GetXaxis().SetTitle(f"m_{{jj}} [GeV]   ({xnote})")
    tot.GetYaxis().SetTitle(count_title)
    tot.GetXaxis().SetTitleSize(0.13); tot.GetXaxis().SetLabelSize(0.11); tot.GetXaxis().SetTitleOffset(1.25)
    tot.GetYaxis().SetTitleSize(0.11); tot.GetYaxis().SetLabelSize(0.10); tot.GetYaxis().SetTitleOffset(0.5)
    tot.GetYaxis().SetNdivisions(503)
    top = max(max(counts), 1)
    tot.SetMinimum(max(0.5, min([n for n in counts if n > 0] or [1]) / 3.0)); tot.SetMaximum(top * 60.0)
    tot.Draw("HIST")
    q = ROOT.TLatex(); q.SetTextFont(42); q.SetTextSize(0.075); KEEP.append(q)
    for b, n in enumerate(counts):
        if n <= 0: continue
        q.SetTextAlign(11 if b == 0 else 21)
        x = EDGES[0] + 0.015 * (EDGES[-1] - EDGES[0]) if b == 0 else 0.5 * (EDGES[b] + EDGES[b + 1])
        q.DrawLatex(x, top * (12.0 if b % 2 == 0 else 3.0), f"{int(n):,}")
    pBot.RedrawAxis()
    c.Print(stem + ".pdf"); c.Print(stem + ".png")
    print("wrote", stem + ".pdf / .png")


def fmt_pct(v): return f"  ({100 * v:.1f}%)"


any_s = next(iter(samples.values()))
JDR, JCUT = any_s["info"]["jet_time_dr"], any_s["info"]["jet_time_dist_cut"]
ANY_RP = any(m.startswith("rp_") for m in DRAWN)
for role, ylab, label in [("sig", "Signal efficiency", args.sig_label), ("bkg", "Background efficiency", args.bkg_label)]:
    if role not in samples: continue
    s = samples[role]
    w = s["e"]["weight"].astype(float)
    series = {m: binned(s["e"], s["kept"][m], EDGES) for m in DRAWN}
    extra = {m: fmt_pct(eff(w, s["kept"][m])[0]) for m in DRAWN}
    counts = [n for _, _, _, n in series[DRAWN[0]]]
    column_figure(f"{args.out}_{role}_eff_mjj", ylab, series, extra, label,
                  [sel_line(s), veto_line(), jet_line(), calib_line(s["calib"]),
                   "all m_{jj} columns in brackets; generator-weighted"
                   + ("; events binned in their no-veto m_{jj}" if ANY_RP else "")],
                  counts, "MC events")

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

# ── S / sqrt(B) per m_jj column relative to no veto ──────────────────────────
# (S/sqrt(B))_method / (S/sqrt(B))_no veto per column. For event vetoes the
# column content does not move and this is exactly eps_S / sqrt(eps_B), with the
# efficiencies' binomial errors; re-paired events move between columns, so for
# those the ratio comes from the yields, with their MC-stat errors (conservative:
# numerator and denominator share events).
ratio = {}
for m in DRAWN:
    ratio[m] = []
    if not m.startswith("rp_"):
        for (es, sl, sh, _), (eb, bl, bh, _) in zip(binned(sg["e"], sg["kept"][m], EDGES),
                                                   binned(bk["e"], bk["kept"][m], EDGES)):
            if not (np.isfinite(es) and np.isfinite(eb)) or eb <= 0:
                ratio[m].append((np.nan, np.nan, np.nan, 0)); continue
            r = es / np.sqrt(eb)
            err = r * np.hypot(0.5 * (sl + sh) / es, 0.25 * (bl + bh) / eb)
            ratio[m].append((r, err, err, 0))
    else:
        for (S, dS), (B, dB), (S0, _), (B0, _) in zip(column_yields(sg, ws, 1.0, m), column_yields(bk, wb, 1.0, m),
                                                    column_yields(sg, ws, 1.0, "none"), column_yields(bk, wb, 1.0, "none")):
            if S <= 0 or B <= 0 or S0 <= 0 or B0 <= 0:
                ratio[m].append((np.nan, np.nan, np.nan, 0)); continue
            r = (S / S0) / np.sqrt(B / B0)
            err = r * np.hypot(dS / S, 0.5 * dB / B)
            ratio[m].append((r, err, err, 0))
extra = {m: f"  ({eff(ws, sg['kept'][m])[0] / np.sqrt(eff(wb, bk['kept'][m])[0]):.2f})" for m in DRAWN}
column_figure(f"{args.out}_s_over_sqrtb", "(S/#sqrt{B}) / (S/#sqrt{B})_{no veto}", ratio, extra,
              f"S: {args.sig_label.split(', ')[-1]},  B: {args.bkg_label.split(', ')[-1]}",
              label_lines(DRAWN) + ["all m_{jj} columns in brackets; > 1 = the timing test raises S/#sqrt{B}"],
              B_COUNTS, "B MC events", xnote=XNOTE)

# ── S / sqrt(B) from the YIELDS per m_jj column (absolute, at --lumi) ───────
# The same information as the ratio above -- that ratio is exactly this divided
# by its no-veto value -- but on the absolute scale, which carries the two
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
                  B_COUNTS, "B MC events", ref=None, methods=SB_METHODS, xnote=XNOTE)


# ── The trade-off: eps_B vs eps_S as the threshold moves (all columns) ───────
SCAN = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0])


def tradeoff_curves(methods):
    curves = {m: [] for m in methods}
    for ns in SCAN:
        rs, rb = outcomes(sg, ns, methods)[0], outcomes(bk, ns, methods)[0]
        for m in methods:
            curves[m].append((eff(ws, rs[m][0])[0], eff(wb, rb[m][0])[0]))
    return {m: np.array(v) for m, v in curves.items()}


def tradeoff_figure(stem, methods, dashed=()):
    curves = tradeoff_curves(methods)
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 800, 700); KEEP.append(c)
    c.SetLeftMargin(0.15); c.SetRightMargin(0.05)
    xs = np.concatenate([p[:, 0] for p in curves.values()]); ys = np.concatenate([p[:, 1] for p in curves.values()])
    ylo = max(0.0, ys.min() - 0.05)
    fr = ROOT.TH1D("fr_" + os.path.basename(stem), "", 10, max(0.0, xs.min() - 0.03), 1.005); KEEP.append(fr)
    fr.SetMinimum(ylo); fr.SetMaximum(1.0 + (1.0 - ylo) * (0.7 + 0.08 * max(0, len(methods) - 5)))
    fr.GetXaxis().SetTitle("Signal efficiency #varepsilon_{S}")
    fr.GetYaxis().SetTitle("Background efficiency #varepsilon_{B}")
    fr.Draw("AXIS")
    leg = ROOT.TLegend(0.18, 0.77 - 0.04 * len(methods), 0.70, 0.77); ROOT.StyleLegend(leg, 0.026); KEEP.append(leg)
    i3 = int(np.argmin(np.abs(SCAN - args.nsigma)))
    for m in methods:
        lab, col, mk = STYLE[m]
        g = ROOT.TGraph(len(SCAN), curves[m][:, 0].copy(), curves[m][:, 1].copy()); KEEP.append(g)
        g.SetLineColor(colour(col)); g.SetLineWidth(2); g.SetMarkerColor(colour(col)); g.SetMarkerStyle(mk); g.SetMarkerSize(0.8)
        if m in dashed: g.SetLineStyle(2)
        g.Draw("LP SAME")
        g3 = ROOT.TGraph(1, curves[m][i3:i3 + 1, 0].copy(), curves[m][i3:i3 + 1, 1].copy()); KEEP.append(g3)
        g3.SetMarkerColor(colour(col)); g3.SetMarkerStyle(mk); g3.SetMarkerSize(2.0); g3.Draw("P SAME")
        leg.AddEntry(g, lab, "lp")
    leg.Draw()
    ROOT.ATLASLabel(0.18, 0.895, "")
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04); KEEP.append(t)
    t.DrawLatex(0.315, 0.895, "Simulation Internal")
    t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.026); t2.SetTextColor(ROOT.kGray + 2); KEEP.append(t2)
    t2.DrawLatex(0.18, 0.855, f"S: {sel_line(sg)}")
    t2.DrawLatex(0.18, 0.825, f"B: {sel_line(bk)}")
    t2.DrawLatex(0.18, 0.795, f"threshold scanned {SCAN[0]:g}-{SCAN[-1]:g}#sigma; large markers at {args.nsigma:g}#sigma"
                 + (";  dashed: jets removed, pair re-formed" if dashed else ""))
    c.Print(stem + ".pdf"); c.Print(stem + ".png")
    print("wrote", stem + ".pdf / .png")
    return curves


tradeoff_figure(args.out + "_tradeoff", DRAWN, dashed=[m for m in DRAWN if m.startswith("rp_")])

# ── Event veto against jet removal + re-pairing, same t0 sources ────────────
if COMPARE:
    cmp_methods = [f"t0_{s}" for s in COMPARE] + [f"rp_{s}" for s in COMPARE]
    tradeoff_figure(args.out + "_veto_vs_repair_tradeoff", cmp_methods, dashed=[f"rp_{s}" for s in COMPARE])
