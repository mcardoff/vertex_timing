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
           for Z+jets |m_ll - 91.2| < 25 GeV. Always applied: vbs_time_veto --sr
           keeps only events with such a pair among all jets, and E_T^miss above
           a 100 GeV floor, so neither can be loosened here
  fjvt     fJVT on j1, j2 (< 0.5 above E_T^miss = 200 GeV, < 0.2 below)
  dphi     |dphi_jj| < 2
  met      E_T^miss > 200 GeV, or > 160 GeV with m_jj > 1.5 TeV (the 2-jet SR bins)
  etjet    E_T^jet,no-jvt > 140 GeV
  jetveto  2-4 jets; j3, j4: C < 0.6 and m_rel < 0.05; with 3-4 jets, E_T^miss >
           200 GeV and m_jj > 1.5 TeV (the 3-4 jet SR bins)
All six together are exactly the SR: checked against sr_select's decision for
every sample and every timing variant.

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

S/sqrt(B) with the background measured in a set M (S always the signal in the SR):
  event veto   eps_S(SR) / sqrt(eps_B(M))
  jet removal  r_S(SR) / sqrt(r_B(M) T(M)), T(M) = r_S(SR) / r_S(M),
the signal's removal transfer from M to the SR: the jet veto's rescue of events
whose failing extra jet is removed, when M lacks the jet veto.

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_sr_minimal.py \\
        --sig condor/vbf/vbf_sr_vbs_time_veto.root --zjets condor/zjets/zjets_sr_vbs_time_veto.root \\
        --dijet condor/dijet/dijet_sr_vbs_time_veto.root --out figs/time_veto/srmin

Writes <out>_summary.md.
"""
import argparse, os
import numpy as np
import vbs_time_veto_common as C
from vbf_hinv_sr import SIG, ZLL, MET_LO, MET_HI, MJJ_LOWBINS, load_samples, prepare, take, ratio

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
ap.add_argument("--out", required=True)
args = ap.parse_args()
MINIMAL = [k for k in KEYS if k in args.minimal.split(",")]
DROPPED = [k for k in KEYS if k not in MINIMAL]


def groups(sel, s):
    """per-event pass flags of each cut group, for one selection (none or a removal)"""
    st = dict(sel["steps"])
    n, met = sel["n"], s["met"]
    with np.errstate(invalid="ignore"):
        wide = sel["mjj_raw"] >= MJJ_LOWBINS
    return {"pair": st["j1 > 80, j2 > 50 GeV"] & st["eta1 eta2 < 0, |deta_jj| > 3.8"] & st["m_jj > 800 GeV"] & s["zwin"],
            "fjvt": st["fJVT on j1, j2"],
            "dphi": st["|dphi_jj| < 2"],
            "met": (met > MET_HI) | ((met > MET_LO) & wide),
            "etjet": st["E_T^jet,no-jvt > 140 GeV"],
            "jetveto": st["2-4 jets > 25 GeV (JVT)"] & st["j3, j4: C < 0.6, m_rel < 0.05"]
                       & ((n == 2) | ((met > MET_HI) & wide))}


def passes(g, keys):
    m = g["pair"].copy()
    for k in keys: m &= g[k]
    return m


samples = load_samples(args.sig, args.zjets, args.dijet, args.slices)
calib = C.measure(samples[SIG]["e"], True, "the VBF file's HS legs")
for name, s in samples.items():
    prepare(s, calib, args.nsigma, args.min_trk)
    s["g"] = groups(s["none"], s)
    s["g_rm"] = {src: groups(s["out"][f"rm_{src}"], s) for src in C.T0_SRC}
    # the decomposition must be the SR exactly, with and without timing
    assert np.array_equal(passes(s["g"], KEYS), s["none"]["ok"]), name
    for src in C.T0_SRC:
        assert np.array_equal(passes(s["g_rm"][src], KEYS), s["out"][f"rm_{src}"]["ok"]), (name, src)
    print(f"[{name}] {len(s['w']):,} candidates, pair {int(s['g']['pair'].sum()):,}, SR {int(s['none']['ok'].sum()):,}")
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
    r, rerr = ratio(s["w"], passes(s["g_rm"][src], keys), m) if m.any() else (np.nan, np.nan)
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
    md.append(f"- **{name}**: `{s['path']}` ({i['njobs']} job(s), {i['n_read']:,} read), E_T^miss = {s['rcol']}, "
              f"{len(s['w']):,} SR candidates")
md += ["", "**The SR's cuts.** `pair` is always applied (the ntuple keeps only events with such a pair among all "
       "jets, above an E_T^miss floor of 100 GeV); the other groups together with it are exactly the SR.", "",
       "| group | cuts |", "|---|---|",
       "| pair | j1 > 80, j2 > 50 GeV (the two leading JVT-passing jets > 25 GeV); eta1 eta2 < 0, |deta_jj| > 3.8; "
       "m_jj > 800 GeV; Z+jets: |m_ll - 91.2| < 25 GeV; E_T^miss > 100 GeV |"]
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
    sig = lambda a, b: (a[0] - b[0]) / np.hypot(a[1], b[1]) if np.hypot(a[1], b[1]) > 0 else np.nan
    return (f"{count(s, p)} / {count(s, q)} | {fmt(*fp)} | {fmt(*fq)} | {sig(fp, fq):+.1f} | {fmt(*vp)} | "
            f"{fmt(*vq)} | {sig(vp, vq):+.1f} |")


SPLIT_HEAD = ["| sample | split | MC pass / fail | f pass | f fail | Delta f | veto kept pass | veto kept fail "
              "| Delta kept |", "|---|---|---|---|---|---|---|---|---|"]
md += ["", f"**Factorisation on the pair** ({args.t0} t0): for the events passing / failing each group, the "
       "PU-tagged share f and the event veto's kept fraction. Delta = pass - fail, in standard deviations. A "
       "group with no significant Delta in either can be dropped from the background measurement.", ""] + SPLIT_HEAD
for name, s in samples.items():
    md += [f"| {name} | {k} | " + split(s, s["g"]["pair"], s["g"][k]) for k in KEYS]
md += ["", f"**Factorisation on the minimal set** ({label(MINIMAL)}) for the groups it drops", ""] + SPLIT_HEAD
for name, s in samples.items():
    md += [f"| {name} | {k} | " + split(s, passes(s["g"], MINIMAL), s["g"][k]) for k in DROPPED]

# the pair's own variables (cut at the ntuple, so only their range above the cut is
# here) and E_T^miss, each split at a value inside that range, on the minimal set
md += ["", f"**The preselected variables, split above their cut**, on the minimal set ({label(MINIMAL)}). No "
       "Delta in a split is a hint, not a proof, that loosening that cut in the ntuple production would add "
       "statistics with the same response.", ""] + SPLIT_HEAD
for name, s in samples.items():
    J, i1, i2 = s["J"], s["none"]["i1"], s["none"]["i2"]
    with np.errstate(invalid="ignore"):
        pv = [("m_jj > 1.5 TeV", s["none"]["mjj_raw"] > 1500.0),
              ("|deta_jj| > 5", np.abs(take(J.jet_eta, i1) - take(J.jet_eta, i2)) > 5.0),
              ("j2 > 80 GeV", take(J.jet_pt, i2) > 80.0),
              ("E_T^miss > 200 GeV", s["met"] > MET_HI)]
    md += [f"| {name} | {lab} | " + split(s, passes(s["g"], MINIMAL), flag) for lab, flag in pv]

# the jet veto's rescue in jet removal, measured on each sample: r(M + jetveto) / r(M)
md += ["", f"**The jet veto's removal transfer** r(M + jetveto) / r(M), measured on every sample ({args.t0} t0): "
       "jet removal's S/sqrt(B) takes the signal's for the background.", "",
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

# ── S/sqrt(B) with Z+jets measured in each set ───────────────────────────────
Z = samples[ZLL]
sr = KEYS
md += ["**S/sqrt(B) relative to no timing, Z+jets measured in each set**: S is the signal in the full SR; "
       "jet removal's B is r_B(M) x T(M), T(M) = r_S(SR) / r_S(M) from the signal. Errors: Z+jets MC statistics.", "",
       "| B measured in | Z+jets MC | t0 | eps_B veto | S/sqrt(B) veto | r_B(M) | T(M) | S/sqrt(B) removal |",
       "|---|---|---|---|---|---|---|---|"]
ZTAB = {}
for keys in LADDER:
    for src in C.T0_SRC:
        rs, rb = response(S, sr, src), response(Z, keys, src)
        T = rs["r"] / response(S, keys, src)["r"]
        zv = rs["v"] / np.sqrt(rb["v"])
        zv_err = 0.5 * zv * rb["verr"] / rb["v"]
        zr = rs["r"] / np.sqrt(rb["r"] * T)
        zr_err = 0.5 * zr * rb["rerr"] / rb["r"]
        ZTAB[label(keys), src] = dict(zv=zv, zv_err=zv_err, zr=zr, zr_err=zr_err, n=rb["n"])
        md.append(f"| {label(keys)}{' (minimal)' if keys == MINIMAL else ''} | {rb['n']} | {src} | {fmt(rb['v'], rb['verr'])} | {100 * (zv - 1):+.1f} +- "
                  f"{100 * zv_err:.1f}% | {fmt(rb['r'], rb['rerr'])} | {T:.3f} | {100 * (zr - 1):+.1f} +- "
                  f"{100 * zr_err:.1f}% |")

os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
open(args.out + "_summary.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
print("wrote", args.out + "_summary.md")
