"""The timing study's event vetoes on util/tdr_fwdpair.cxx's trees (numpy only, so the BDT environment and the
ROOT one can both import it): the jet-time calibration, the pulls, and the 3-sigma decisions, exactly as
python/vbs_time_veto_plot.py decisions() makes them, on tdr_fwdpair's core jet time (core_t / core_sig /
core_n, commit 53d16c3 on). Used by tdr_fwdpair_kinplots.py (stacks after the veto), tdr_fwdpair_bdt.py
(the pulls as BDT inputs) and tdr_fwdpair_plots.py (the rejection curves).

  jet vs jet  both legs timed: |t_1 - t_2| / sqrt(s_1^2 + s_2^2)
  jet vs t0   per leg and t0 source: |t_jet - t0| / sqrt(s_jet^2 + (infl sigma_t0)^2)
s_jet = f(n) x sigma_jet; f(n) is the MAD of (t_core - t_truth) / sigma over the VBF file's hard-scatter legs, per
core-track bin (vbs_time_veto_common's NBINS; the global MAD below 20 legs). A leg is timed with >= 1 core track.
A pull is NaN where it cannot be formed (an untimed leg, a missing t0).
"""
import numpy as np

def jvt_label(wp):
    """the pileup-tagging baseline in a figure's header, from the trees' meta jvt_wp (tdr_fwdpair --jvt)"""
    return {"loose": "JVT + fJVT loose", "tight": "JVT + fJVT tight",
            "jvt4pt60": "JVT to |#eta| < 4.0, no fJVT",            # JVT's own p_T < 60 GeV range: on the slides, not in
            "jvt4pt120": "JVT to |#eta| < 4.0 (p_{T} < 120), no fJVT"}.get(wp, f"JVT + fJVT {wp}")   # the figures (too long)


def jvt_short(wp):
    """a legend-length form of jvt_label"""
    return {"loose": "JVT + fJVT loose", "tight": "JVT + fJVT tight", "jvt4pt60": "JVT to |#eta| < 4.0",
            "jvt4pt120": "JVT to |#eta| < 4.0 (p_{T} < 120)"}.get(wp, f"JVT + fJVT {wp}")


NBINS = [1, 2, 3, 4, 6, 9, 10**6]                    # python/vbs_time_veto_common.py's
T0_SOURCES = ("hgtd", "trkptz", "waves", "truth")
# util/tdr_common.h's truth-fixed WAVeS t0s (2026-11-13, todo vt-23j): perfect HS track times, perfect selection, both.
# Vetoed on whenever a tree carries them (sources(E)); T0_SOURCES itself is unchanged, so the BDT features are too.
T0_FIXED = ("wfix_time", "wfix_sel", "wfix_both")


def sources(E):
    """the t0 sources the trees E carry: T0_SOURCES, plus the truth-fixed ones when present"""
    return T0_SOURCES + tuple(s for s in T0_FIXED if f"t0_{s}" in E)


def calibration(E, verbose=True):
    """f(n) from the hard-scatter legs of the trees E (the VBF file)"""
    p, n = [], []
    for k in (1, 2):
        m = (E[f"j{k}_hs"] == 1) & (E[f"j{k}_core_n"] > 0)
        p.append((E[f"j{k}_core_t"][m] - E["t0_truth"][m]) / E[f"j{k}_core_sig"][m]); n.append(E[f"j{k}_core_n"][m])
    p, n = np.concatenate(p), np.concatenate(n)
    mad = lambda x: 1.4826 * float(np.median(np.abs(x)))
    g = mad(p)
    fac = [mad(p[(n >= lo) & (n < hi)]) if np.sum((n >= lo) & (n < hi)) >= 20 else g for lo, hi in zip(NBINS[:-1], NBINS[1:])]
    if verbose:
        print(f"jet-time calibration from {len(p):,} VBF hard-scatter legs: global {g:.3f}; per n bin " +
              ", ".join(f"{lo}{'' if hi == lo + 1 else '+' if hi > 10**5 else f'-{hi - 1}'}: {x:.2f}"
                        for (lo, hi), x in zip(zip(NBINS[:-1], NBINS[1:]), fac)))

    def f(nn):
        out = np.full(len(nn), g)
        for (lo, hi), x in zip(zip(NBINS[:-1], NBINS[1:]), fac): out[(nn >= lo) & (nn < hi)] = x
        return out
    return f


def pulls(E, f):
    """{"jj": ..., "<src>_1": ..., "<src>_2": ...}: the veto pulls per event, NaN where undefined"""
    t = {k: E[f"j{k}_core_t"] for k in (1, 2)}
    sg = {k: f(E[f"j{k}_core_n"]) * E[f"j{k}_core_sig"] for k in (1, 2)}
    h = {k: E[f"j{k}_core_n"] >= 1 for k in (1, 2)}
    both = h[1] & h[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        out = {"jj": np.where(both, np.abs(t[1] - t[2]) / np.sqrt(np.where(both, sg[1] ** 2 + sg[2] ** 2, 1.0)), np.nan)}
        for s in sources(E):
            ok0, t0, s0 = E[f"ok_{s}"].astype(bool), E[f"t0_{s}"], E[f"infl_{s}"] * E[f"sig0_{s}"]
            for k in (1, 2):
                out[f"{s}_{k}"] = np.where(ok0 & h[k], np.abs(t[k] - t0) / np.sqrt(sg[k] ** 2 + s0 ** 2), np.nan)
    return out


def vetoes(E, f, nsigma=3.0):
    """kept[method] per event, method in jj and the t0 sources: vbs_time_veto_plot.decisions()"""
    P = pulls(E, f)
    kept = {"jj": ~(np.nan_to_num(P["jj"], nan=0.0) >= nsigma)}
    for s in sources(E):
        kept[s] = ~((np.nan_to_num(P[f"{s}_1"], nan=0.0) >= nsigma) | (np.nan_to_num(P[f"{s}_2"], nan=0.0) >= nsigma))
    return kept
