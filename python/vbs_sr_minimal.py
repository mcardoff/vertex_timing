#!/usr/bin/env python3
"""A minimal cut set for measuring the background's response to timing in the
VBF H->inv signal region (ANA-EXOT-2020-11-INT1 Sec. 4.5.2, python/vbf_hinv_sr.py).

The full SR leaves 4 Z+jets MC events, too few to measure B(timing) / B(no
timing). That ratio, not B itself, is what S/sqrt(B) needs, and a cut that is
uncorrelated with the timing decision cancels in it exactly. So the background
can be measured with only the cuts the timing decision depends on, and the
dropped ones' share of the yield never enters.

The SR's cuts, in groups (the tagging pair: the two leading JVT-passing jets
> 25 GeV):
  pair     j1 > 80, j2 > 50 GeV; eta1 eta2 < 0, |deta_jj| > 3.8; m_jj > 800 GeV;
           for Z+jets |m_ll - 91.2| < 25 GeV; E_T^miss above --met-floor (100 GeV).
           Always applied: vbs_time_veto --sr keeps only events with such a pair
           among all jets, and E_T^miss above its --sr-recoil-min (100 GeV unless
           the file is tagged sr_rec<v>), so neither can be loosened below that
  fjvt     fJVT on j1, j2 (< 0.5 above E_T^miss = 200 GeV, < 0.2 below)
  dphi     |dphi_jj| < 2
  met      E_T^miss > 200 GeV, or > 160 GeV with m_jj > 1.5 TeV (the 2-jet SR bins)
  etjet    E_T^jet,no-jvt > 140 GeV
  jetveto  2-4 jets; j3, j4: C < 0.6 and m_rel < 0.05; with 3-4 jets, E_T^miss >
           200 GeV and m_jj > 1.5 TeV (the 3-4 jet SR bins)
All six together are exactly the SR (at the region's thresholds): checked against
sr_select's cut flags for every sample and every timing variant.

Tests, with the event veto (j1, j2 against t0 at --nsigma) and jet removal
(failing jets removed, the selection re-run):
  * factorisation: the event veto's kept fraction and the pileup-tagged share f
    for the events passing each group against those failing it, on the pair
    alone and on the minimal set (--minimal) for the groups it drops. A group
    that selects a different response must stay;
  * the pair's own variables and E_T^miss split at a value above their cut
    (m_jj 1.5 TeV, |deta_jj| 5, j2 80 GeV, E_T^miss 200 GeV): whether loosening
    the ntuple preselection could add statistics without changing the response;
  * the jet veto's removal transfer, measured on every sample;
  * the ladder: cumulative sets from the pair to the SR, each group alone, N-1.
The signal (18k MC events in the SR) shows each cut's effect precisely; Z+jets
and the multijet JZ slices within their statistics.

S/sqrt(B) is only ever formed with S and B in the SAME set M:
  event veto   eps_S(M) / sqrt(eps_B(M));   jet removal   r_S(M) / sqrt(r_B(M)).
Nothing is transferred between sets. The jet veto's rescue of jet removal in the
full SR (a removed out-of-time jet stops failing the extra-jet cuts) is therefore
not in the minimal set's numbers; it is measured per sample, for the record.

The break-even, also in the minimal set: the event veto keeps a background at
eps_B = f (1 - c) + (1 - f)(1 - h), c the catch rate on pileup-tagged pairs and h
the false-veto rate on the rest, both from the signal's own pairs there, so
S/sqrt(B) rises when f > f* = (1 - h - eps_S^2) / (c - h).

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_sr_minimal.py \\
        --sig condor/vbf/vbf_sr_vbs_time_veto.root --zjets condor/zjets/zjets_sr_vbs_time_veto.root \\
        --dijet condor/dijet/dijet_sr_vbs_time_veto.root --out figs/time_veto/srcuts
(the minimal set's own figures, S and B both in it, are vbs_time_veto_sr.py --region minimal's)

--region relaxed runs every test on the relaxed set instead (vbf_hinv_sr.py's REGIONS):
j1 > 40, j2 > 30 GeV and |m_ll - 91.2| < 40 GeV, on the vbs_time_veto --sr-j1-pt=40
--sr-j2-pt=30 outputs. Its splits then include the SR's own thresholds (j1 > 80, j2 > 50
GeV, the 25 GeV Z window): whether the events the relaxation adds respond like the rest.

Writes <out>_summary.md and <out>_breakeven.(pdf|png).
"""
import argparse, os
import numpy as np
import vbs_time_veto_common as C
from vbf_hinv_sr import (SIG, MET_LO, MET_HI, MJJ_LOWBINS, MZ, MZ_WIN, SR_J1_PT, SR_J2_PT, REGIONS,
                         load_samples, prepare, take, ratio)

GROUPS = [("fjvt", "fJVT on j1, j2"),
          ("dphi", "|dphi_jj| < 2"),
          ("met", "E_T^miss > 200 GeV, or > 160 GeV with m_jj > 1.5 TeV"),
          ("etjet", "E_T^jet,no-jvt > 140 GeV"),
          ("jetveto", "2-4 jets; j3, j4 C < 0.6, m_rel < 0.05 (3-4 jets: E_T^miss > 200, m_jj > 1.5 TeV)")]
KEYS = [k for k, _ in GROUPS]
LADDER = [KEYS[:i] for i in range(len(KEYS) + 1)]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True, help="VBF H->inv vbs_time_veto --sr output (RAW ntuples)")
ap.add_argument("--zjets", required=True, help="Z+jets vbs_time_veto --sr output")
ap.add_argument("--dijet", help="dijet vbs_time_veto --sr --sr-met=jets output (RAW ntuples)")
ap.add_argument("--slices", default="JZ2,JZ3,JZ4", help="dijet JZ slices to show")
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--min-trk", type=int, default=1)
ap.add_argument("--t0", default="waves", choices=C.T0_SRC, help="the t0 of the per-cut tables")
ap.add_argument("--minimal", default="fjvt,dphi", help="the groups the minimal set keeps besides the pair")
ap.add_argument("--region", default="minimal", choices=["minimal", "relaxed"],
                help="the pair's thresholds (vbf_hinv_sr.py's REGIONS): the SR's, or the relaxed set's")
ap.add_argument("--met-floor", type=float, default=None,
                help="E_T^miss floor of the pair (default: the region's; a sample's own ntuple floor applies if higher)")
ap.add_argument("--calib-file", help="measure the jet-time correction f(n) on this signal file's HS legs instead of --sig's (to compare two productions under one calibration)")
ap.add_argument("--out", required=True)
args = ap.parse_args()
REG = REGIONS[args.region]
if args.met_floor is None: args.met_floor = REG["met_floor"]
J1, J2, ZW = REG.get("j1_pt", SR_J1_PT), REG.get("j2_pt", SR_J2_PT), REG.get("z_win", MZ_WIN)
MINIMAL = [k for k in KEYS if k in args.minimal.split(",")]
DROPPED = [k for k in KEYS if k not in MINIMAL]


def groups(sel, s):
    """per-event pass flags of each cut group, for one selection (none or a removal)"""
    fl = sel["flags"]   # at the region's thresholds
    n, met = sel["n"], s["met"]
    with np.errstate(invalid="ignore"):
        wide = sel["mjj_raw"] >= MJJ_LOWBINS
    return {"pair": fl["pt"] & fl["deta"] & fl["mjj"] & fl.get("zwin", True) & (met >= args.met_floor),
            "fjvt": fl["fjvt"],
            "dphi": fl["dphi"],
            "met": (met > MET_HI) | ((met > MET_LO) & wide),
            "etjet": fl["etjet"],
            "jetveto": fl["njet"] & fl["extra"] & ((n == 2) | ((met > MET_HI) & wide))}


def passes(g, keys):
    m = g["pair"].copy()
    for k in keys: m &= g[k]
    return m


samples = load_samples(args.sig, args.zjets, args.dijet, args.slices)
calib = C.measure(samples[SIG]["e"] if not args.calib_file else C.legs(args.calib_file), True,
                  "the VBF file's HS legs" if not args.calib_file else f"{args.calib_file}'s HS legs")
for name, s in samples.items():
    s["floor"] = max(args.met_floor, s["info"]["sr_recoil_min"])
    prepare(s, calib, args.nsigma, args.min_trk, args.region)
    s["g"] = groups(s["none"], s)
    s["g_rm"] = {src: groups(s["out"][f"rm_{src}"], s) for src in C.T0_SRC}
    # the decomposition must be every SR cut (at the region's thresholds) exactly, with and without timing
    every = lambda sel: np.logical_and.reduce(list(sel["flags"].values())) & (s["met"] >= args.met_floor)
    assert np.array_equal(passes(s["g"], KEYS), every(s["none"])), name
    for src in C.T0_SRC:
        assert np.array_equal(passes(s["g_rm"][src], KEYS), every(s["out"][f"rm_{src}"])), (name, src)
    # and the default minimal set is the region's selection, which vbs_time_veto_sr.py --region uses
    if MINIMAL == ["fjvt", "dphi"] and args.met_floor == REG["met_floor"]:
        assert np.array_equal(passes(s["g"], MINIMAL), s["none"]["ok"]), name
    print(f"[{name}] {len(s['w']):,} candidates, E_T^miss floor {s['floor']:g} GeV, pair {int(s['g']['pair'].sum()):,}, "
          f"{args.region} set {int(s['none']['ok'].sum()):,}")
S = samples[SIG]
BKG = [k for k in samples if k != SIG]


def neff(s, m):
    w = s["w"][m]
    return w.sum() ** 2 / (w ** 2).sum() if (w ** 2).sum() > 0 else 0.0


def eff(s, m, kept):
    """weighted fraction of m kept, and its 1-sigma error (Wilson's, so a column of a few
    events whose error the binomial formula sets to zero still gets one)"""
    if not m.any(): return np.nan, np.nan
    e, _, lo, hi = C.eff(s["w"][m], kept[m])
    return e, 0.5 * (lo + hi)


def response(s, keys, src):
    """the set's MC count, f, the event veto's kept fraction and jet removal's ratio"""
    m = passes(s["g"], keys)
    f, ferr = eff(s, m, s["pu_tag"])
    v, verr = eff(s, m, ~s["vetoed"][src])
    mr = passes(s["g_rm"][src], keys)
    r, rerr = ratio(s["w"], mr, m) if m.any() else (np.nan, np.nan)
    # the ratio's error (from the events removal changes) vanishes when it changes none of
    # a few events: floor it at the Wilson error of the fraction it leaves unchanged
    if m.any(): rerr = max(rerr, eff(s, m, mr)[1])
    return dict(n=int(m.sum()), neff=neff(s, m), f=f, ferr=ferr, v=v, verr=verr, r=r, rerr=rerr)


def label(keys):
    if not keys: return "pair"
    if keys == KEYS: return "SR"
    return "pair + " + " + ".join(keys)


def fmt(x, e, d=2):
    return "-" if not np.isfinite(x) else f"{x:.{d}f} +- {e:.{d}f}"


def count(s, m):
    n = int(m.sum())
    return f"{n:,}" if s["role"] != "dijet" else f"{n:,} ({neff(s, m):,.0f})"


md = [f"# A minimal SR cut set for the background's timing response: {os.path.basename(args.out)}", "",
      f"Event veto on j1, j2 at {args.nsigma:g} sigma; jet removal re-runs the selection on the kept jets. "
      f"Jet-time sigma {calib.md()}. Per-cut tables use the {args.t0} t0. Multijet: MC events with the "
      "effective count (sum w)^2 / sum w^2 in brackets; the JZ slices are not normalised to each other.", ""]
for name, s in samples.items():
    i = s["info"]
    md.append(f"- **{name}**: `{s['path']}` ({i['njobs']} job(s), {i['n_read']:,} read), E_T^miss = {s['rcol']} "
              f"> {s['floor']:g} GeV, {len(s['w']):,} SR candidates")
md += ["", "**The SR's cuts.** `pair` is always applied (the ntuple keeps only events with such a pair among all "
       "jets, above its E_T^miss floor); the other groups together with it are exactly the SR.", "",
       "| group | cuts |", "|---|---|",
       f"| pair | j1 > {J1:g}, j2 > {J2:g} GeV (the two leading JVT-passing jets > 25 GeV); eta1 eta2 < 0, "
       f"|deta_jj| > 3.8; m_jj > 800 GeV; Z+jets: |m_ll - 91.2| < {ZW:g} GeV; E_T^miss above the floor "
       "(per sample, above) |"]
md += [f"| {k} | {t} |" for k, t in GROUPS]

# ── MC events per set ────────────────────────────────────────────────────────
SETS = ([("ladder", k) for k in LADDER] + [("alone", [k]) for k in KEYS[1:]]
        + [("N-1", [x for x in KEYS if x != k]) for k in KEYS])
md += ["", "**MC events per cut set** (ladder: cumulative; alone: the pair and one group; N-1: the SR without one "
       "group)", "", "| | set | " + " | ".join(samples) + " |", "|---|---|" + "---|" * len(samples)]
for kind, keys in SETS:
    md.append(f"| {kind} | {label(keys)} | " + " | ".join(count(s, passes(s["g"], keys)) for s in samples.values()) + " |")

# ── Factorisation: does a cut select a different timing response? ───────────
def split(s, base, flag):
    """f and the event veto's kept fraction for base & flag against base & ~flag"""
    p, q = base & flag, base & ~flag
    fp, fq = eff(s, p, s["pu_tag"]), eff(s, q, s["pu_tag"])
    vp, vq = eff(s, p, ~s["vetoed"][args.t0]), eff(s, q, ~s["vetoed"][args.t0])
    return (f"{count(s, p)} / {count(s, q)} | {fmt(*fp)} | {fmt(*fq)} | {delta(s, p, q, s['pu_tag']):+.1f} | "
            f"{fmt(*vp)} | {fmt(*vq)} | {delta(s, p, q, ~s['vetoed'][args.t0]):+.1f} |")


def delta(s, p, q, flag):
    """(fraction of flag in p) - (in q), in standard deviations of the pooled two-proportion
    test at the effective counts: unlike the per-side errors, it does not collapse when one
    side is a few events all on one side"""
    if not (p.any() and q.any()): return np.nan
    w = s["w"]
    fp, fq = w[p & flag].sum() / w[p].sum(), w[q & flag].sum() / w[q].sum()
    pool = w[(p | q) & flag].sum() / w[p | q].sum()
    var = pool * (1 - pool) * (1 / neff(s, p) + 1 / neff(s, q))
    return (fp - fq) / np.sqrt(var) if var > 0 else 0.0


SPLIT_HEAD = ["| sample | split | MC pass / fail | f pass | f fail | Delta f | veto kept pass | veto kept fail "
              "| Delta kept |", "|---|---|---|---|---|---|---|---|---|"]
md += ["", f"**Factorisation on the pair** ({args.t0} t0): for the events passing / failing each group, the "
       "PU-tagged share f and the event veto's kept fraction. Delta = pass - fail, in standard deviations of the "
       "pooled two-proportion test. A group with no significant Delta in either can be dropped from the "
       "background measurement.", ""] + SPLIT_HEAD
for name, s in samples.items():
    md += [f"| {name} | {k} | " + split(s, s["g"]["pair"], s["g"][k]) for k in KEYS]
md += ["", f"**Factorisation on the {args.region} set** ({label(MINIMAL)}) for the groups it drops", ""] + SPLIT_HEAD
for name, s in samples.items():
    md += [f"| {name} | {k} | " + split(s, passes(s["g"], MINIMAL), s["g"][k]) for k in DROPPED]

# the pair's own variables (cut at the ntuple, so only their range above the cut is
# here) and E_T^miss, each split at a value inside that range, on the minimal set
md += ["", f"**The preselected variables, split above their cut**, on the {args.region} set ({label(MINIMAL)}). No "
       "Delta in a split is a hint, not a proof, that loosening that cut in the ntuple production would add "
       "statistics with the same response.", ""] + SPLIT_HEAD
for name, s in samples.items():
    J, i1, i2 = s["J"], s["none"]["i1"], s["none"]["i2"]
    with np.errstate(invalid="ignore"):
        pv = [("m_jj > 1.5 TeV", s["none"]["mjj_raw"] > 1500.0),
              ("|deta_jj| > 5", np.abs(take(J.jet_eta, i1) - take(J.jet_eta, i2)) > 5.0),
              ("E_T^miss > 200 GeV", s["met"] > MET_HI)]
        if J2 < SR_J2_PT:   # relaxed: do the events below the SR's own thresholds respond like the rest?
            pv += [(f"j1 > {SR_J1_PT:g} GeV", take(J.jet_pt, i1) > SR_J1_PT),
                   (f"j2 > {SR_J2_PT:g} GeV", take(J.jet_pt, i2) > SR_J2_PT),
                   (f"j1 > {SR_J1_PT:g} and j2 > {SR_J2_PT:g} GeV",
                    (take(J.jet_pt, i1) > SR_J1_PT) & (take(J.jet_pt, i2) > SR_J2_PT))]
        else:
            pv += [("j2 > 80 GeV", take(J.jet_pt, i2) > 80.0)]
        if s["zll"] and ZW > MZ_WIN:
            pv += [(f"|m_ll - 91.2| < {MZ_WIN:g} GeV", np.abs(s["e"]["z_mass"] - MZ) < MZ_WIN)]
    md += [f"| {name} | {lab} | " + split(s, passes(s["g"], MINIMAL), flag) for lab, flag in pv]

# E_T^miss binned on the minimal set: whether the floor could go lower (a file with a
# lower ntuple floor, sr_rec<v>, fills the bins below 100 GeV)
MET_BINS = [0.0, 50.0, 100.0, 160.0, 200.0, np.inf]
md += ["", f"**E_T^miss, binned**, on the {args.region} set ({label(MINIMAL)}, {args.t0} t0; Z+jets: pT(ll))", "",
       "| sample | E_T^miss [GeV] | MC | f | veto kept | removal ratio |", "|---|---|---|---|---|---|"]
for name, s in samples.items():
    for lo, hi in zip(MET_BINS[:-1], MET_BINS[1:]):
        if hi <= s["floor"]: continue
        inbin = (s["met"] >= lo) & (s["met"] < hi)
        m = passes(s["g"], MINIMAL) & inbin
        f, v = eff(s, m, s["pu_tag"]), eff(s, m, ~s["vetoed"][args.t0])
        r = ratio(s["w"], passes(s["g_rm"][args.t0], MINIMAL) & inbin, m) if m.any() else (np.nan, np.nan)
        md.append(f"| {name} | {lo:g}-{hi:g} | {count(s, m)} | {fmt(*f)} | {fmt(*v, 3)} | {fmt(*r, 3)} |")

# the jet veto's rescue in jet removal, measured on each sample: r(M + jetveto) / r(M)
md += ["", f"**The jet veto's removal transfer** r(M + jetveto) / r(M), measured on every sample ({args.t0} t0): "
       f"the rescue jet removal gets from the jet veto in the full SR, which the {args.region} set does not include. "
       "For the record only: no S/sqrt(B) below uses it.", "",
       "| sample | M | r(M) | MC in M + jetveto | r(M + jetveto) | transfer |", "|---|---|---|---|---|---|"]
TRANSFER = {}
for name, s in samples.items():
    for keys in [[], MINIMAL]:
        a, b = response(s, keys, args.t0), response(s, keys + ["jetveto"], args.t0)
        t = b["r"] / a["r"]
        t_err = t * np.hypot(b["rerr"] / b["r"], a["rerr"] / a["r"]) if np.isfinite(t) else np.nan
        TRANSFER[name, label(keys)] = (t, t_err)
        md.append(f"| {name} | {label(keys)} | {fmt(a['r'], a['rerr'], 3)} | {count(s, passes(s['g'], keys + ['jetveto']))} "
                  f"| {fmt(b['r'], b['rerr'], 3)} | {fmt(t, t_err, 3)} |")

# ── The response as composition x per-type response ──────────────────────────
# The veto keeps a background at eps_B = f (1 - c) + (1 - f)(1 - h): c its catch
# rate on PU-tagged pairs, h its false-veto rate on the rest. If c and h are the
# signal's (thousands of pairs of each type), the background MC is needed for f
# alone. Checked on the minimal set: each background's c against the signal's,
# and its measured eps_B against the model at its own f.
M = passes(S["g"], MINIMAL)
md += ["", f"**Composition x per-type response**, on the {args.region} set ({label(MINIMAL)}): c = 1 - kept on PU-tagged "
       "pairs, h = 1 - kept on the rest; the model eps_B = f (1 - c_S) + (1 - f)(1 - h_S) with the signal's c_S, h_S "
       "and the background's own f, against its measured eps_B.", "",
       "| sample | t0 | f | c | h | eps_B measured | eps_B model (signal c, h) | Delta / sigma |",
       "|---|---|---|---|---|---|---|---|"]
for src in C.T0_SRC:
    kc = eff(S, M & S["pu_tag"], ~S["vetoed"][src]); kh = eff(S, M & ~S["pu_tag"], ~S["vetoed"][src])
    cS, hS = 1 - kc[0], 1 - kh[0]
    md.append(f"| {SIG} | {src} | {fmt(*eff(S, M, S['pu_tag']), 3)} | {fmt(cS, kc[1], 3)} | {fmt(hS, kh[1], 3)} "
              f"| {fmt(*eff(S, M, ~S['vetoed'][src]), 3)} | - | - |")
    for name in BKG:
        s = samples[name]
        m = passes(s["g"], MINIMAL)
        if not m.any(): continue
        f, ferr = eff(s, m, s["pu_tag"])
        kc_b, kh_b = eff(s, m & s["pu_tag"], ~s["vetoed"][src]), eff(s, m & ~s["pu_tag"], ~s["vetoed"][src])
        v, verr = eff(s, m, ~s["vetoed"][src])
        model = f * (1 - cS) + (1 - f) * (1 - hS)
        model_err = ferr * abs(hS - cS)
        d = (v - model) / np.hypot(verr, model_err)
        md.append(f"| {name} | {src} | {fmt(f, ferr, 3)} | {fmt(1 - kc_b[0], kc_b[1], 3)} | {fmt(1 - kh_b[0], kh_b[1], 3)} "
                  f"| {fmt(v, verr, 3)} | {fmt(model, model_err, 3)} | {d:+.1f} |")

# ── The response along the ladder, each group alone, and N-1 ──────────────────
md += ["", f"**The timing response per cut set** ({args.t0} t0): f, the event veto's kept fraction and jet "
       "removal's yield ratio (removal / no timing, the selection re-run on the kept jets)", ""]
for name, s in samples.items():
    md += [f"*{name}*", "", "| | set | MC | f | veto kept | removal ratio |", "|---|---|---|---|---|---|"]
    for kind, keys in SETS:
        r = response(s, keys, args.t0)
        md.append(f"| {kind} | {label(keys)} | {count(s, passes(s['g'], keys))} | {fmt(r['f'], r['ferr'])} "
                  f"| {fmt(r['v'], r['verr'], 3)} | {fmt(r['r'], r['rerr'], 3)} |")
    md.append("")

# ── S/sqrt(B): S and B in the same set ───────────────────────────────────────
def zchange(rs, rb, key):
    """S/sqrt(B) relative to no timing for S and B measured in the same set, and its error"""
    a, ea, b, eb = rs[key], rs[key + "err"], rb[key], rb[key + "err"]
    if not (a > 0 and b > 0): return np.nan, np.nan
    z = a / np.sqrt(b)
    return z, z * np.hypot(ea / a, 0.5 * eb / b)


md += ["", "**S/sqrt(B) relative to no timing, S and B in the same set** (errors: both samples' MC statistics; "
       "the multijet slices each on their own, as their normalisations cancel within a slice)", "",
       "| set | B | B MC | t0 | eps_S veto | eps_B veto | S/sqrt(B) veto | r_S removal | r_B removal | S/sqrt(B) removal |",
       "|---|---|---|---|---|---|---|---|---|---|"]
ZTAB = {}
for keys in LADDER:
    for name in BKG:
        if not passes(samples[name]["g"], keys).any(): continue
        for src in C.T0_SRC:
            rs, rb = response(S, keys, src), response(samples[name], keys, src)
            zv, zv_err = zchange(rs, rb, "v")
            zr, zr_err = zchange(rs, rb, "r")
            ZTAB[label(keys), name, src] = dict(zv=zv, zv_err=zv_err, zr=zr, zr_err=zr_err, n=rb["n"])
            md.append(f"| {label(keys)}{f' ({args.region})' if keys == MINIMAL else ''} | {name} | {count(samples[name], passes(samples[name]['g'], keys))} "
                      f"| {src} | {fmt(rs['v'], rs['verr'], 3)} | {fmt(rb['v'], rb['verr'])} | {100 * (zv - 1):+.1f} +- "
                      f"{100 * zv_err:.1f}% | {fmt(rs['r'], rs['rerr'], 3)} | {fmt(rb['r'], rb['rerr'])} | "
                      f"{100 * (zr - 1):+.1f} +- {100 * zr_err:.1f}% |")

# ── Break-even in the minimal set ────────────────────────────────────────────
M = passes(S["g"], MINIMAL)
BE = {}
for src in C.T0_SRC:
    eS = eff(S, M, ~S["vetoed"][src])[0]
    c = 1 - eff(S, M & S["pu_tag"], ~S["vetoed"][src])[0]
    h = 1 - eff(S, M & ~S["pu_tag"], ~S["vetoed"][src])[0]
    BE[src] = dict(eS=eS, c=c, h=h, fstar=(1 - h - eS ** 2) / (c - h) if c > h else np.nan)
md += ["", f"**Break-even PU-tagged background share f\\*, in the {args.region} set** ({label(MINIMAL)}): eps_S, c and h "
       "from the signal there; the event veto raises S/sqrt(B) for a background with f > f*.", "",
       "| t0 | eps_S | c | h | f* |", "|---|---|---|---|---|"]
md += [f"| {src} | {v['eS']:.3f} | {v['c']:.3f} | {v['h']:.3f} | {100 * v['fstar']:.0f}% |" for src, v in BE.items()]
md += ["", f"Each background's f in the {args.region} set: " + "; ".join(
    f"{n} {fmt(*eff(samples[n], passes(samples[n]['g'], MINIMAL), samples[n]['pu_tag']))}"
    for n in BKG if passes(samples[n]["g"], MINIMAL).any())]

os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
open(args.out + "_summary.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
print("wrote", args.out + "_summary.md")

# ── Figure: the break-even curves and each background, all in the minimal set ─
ROOT = C.setup_root(args.out)
KEEP = C.KEEP
c = ROOT.TCanvas("c_be", "", 800, 650); KEEP.append(c)
c.SetLeftMargin(0.14); c.SetRightMargin(0.05)
fr = ROOT.TH1D("fr_be", "", 10, 0.0, 100.0); KEEP.append(fr)
fr.SetMinimum(-25); fr.SetMaximum(25)
fr.GetXaxis().SetTitle("Background with a pileup tagging jet  f [%]")
fr.GetYaxis().SetTitle("S/#sqrt{B} change, event veto [%]")
fr.Draw("AXIS")
zero = ROOT.TLine(0, 0, 100, 0); zero.SetLineStyle(2); zero.SetLineColor(ROOT.kGray + 1); zero.Draw(); KEEP.append(zero)
leg = ROOT.TLegend(0.17, 0.56, 0.55, 0.74); ROOT.StyleLegend(leg, 0.03); KEEP.append(leg)
fgrid = np.linspace(0.0, 1.0, 101)
for src in ["trkptz", "waves", "hgtd", "truth"]:
    v = BE[src]
    z = 100 * (v["eS"] / np.sqrt(fgrid * (1 - v["c"]) + (1 - fgrid) * (1 - v["h"])) - 1)
    g = ROOT.TGraph(len(fgrid), 100 * fgrid, z); KEEP.append(g)
    lab, col, _ = C.STYLE[f"t0_{src}"]
    g.SetLineColor(C.colour(col)); g.SetLineWidth(3); g.Draw("L SAME")
    leg.AddEntry(g, lab.replace("Jet vs t_{0}: ", "t_{0}: ") + f"  (f* = {100 * v['fstar']:.0f}%)", "l")
leg.Draw()
# each background measured directly in the minimal set, with the signal's eps_S there (t0 = --t0)
t0name = {"waves": "WAVeS", "trkptz": "TRKPTZ", "hgtd": "HGTD", "tzp": "TZP", "truth": "truth"}[args.t0]
# lower left: below the curves (S/sqrt(B) at small f is the signal loss alone) and clear of every point
leg2 = ROOT.TLegend(0.20, 0.19, 0.56, 0.20 + 0.04 * (len(BKG) + 1)); ROOT.StyleLegend(leg2, 0.027); KEEP.append(leg2)
leg2.SetHeader(f"measured, {t0name} t_{{0}}:", "")
for k, n in enumerate(BKG):
    s = samples[n]
    m = passes(s["g"], MINIMAL)
    if not m.any(): continue
    f, _ = eff(s, m, s["pu_tag"])
    eb, eerr = eff(s, m, ~s["vetoed"][args.t0])
    if not eb > 0: continue
    es = BE[args.t0]["eS"]
    g = ROOT.TGraphErrors(1); KEEP.append(g)
    g.SetPoint(0, 100 * f, 100 * (es / np.sqrt(eb) - 1))
    g.SetPointError(0, 100 * np.sqrt(max(f * (1 - f), 1e-12) / max(neff(s, m), 1)), 100 * es * 0.5 * eb ** -1.5 * eerr)
    g.SetMarkerStyle([20, 21, 22, 23, 33, 34][k % 6]); g.SetMarkerSize(1.5)
    g.SetMarkerColor(ROOT.kBlack); g.SetLineColor(ROOT.kGray + 2)
    g.Draw("P SAME")
    leg2.AddEntry(g, f"{n}  ({count(s, m)} MC)", "p")
leg2.Draw()
ROOT.ATLASLabel(0.17, 0.88, "")
a = ROOT.TLatex(); a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.04); KEEP.append(a)
a.DrawLatex(0.31, 0.88, "Simulation Internal")
g2 = ROOT.TLatex(); g2.SetNDC(); g2.SetTextFont(42); g2.SetTextSize(0.028); g2.SetTextColor(ROOT.kGray + 2); KEEP.append(g2)
g2.DrawLatex(0.17, 0.84, f"S and B both in the {args.region} set: SR pair with j_{{1,2}} > {J1:g}, {J2:g} GeV,")
g2.DrawLatex(0.17, 0.81, f"fJVT, |#Delta#phi_{{jj}}| < 2, E_{{T}}^{{miss}} > {args.met_floor:g} GeV.  Lines: #varepsilon_{{S}} and "
             "the test's response,")
g2.DrawLatex(0.17, 0.78, "pileup / hard-scatter tagging pairs, from the signal")
stem = f"{args.out}_breakeven"
c.Print(stem + ".pdf"); c.Print(stem + ".png")
print("wrote", stem + ".pdf / .png")
