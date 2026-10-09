#!/usr/bin/env python3
"""Event cards for the jet-vs-t0 timing veto: one figure per event, drawn from the
util/vbs_time_veto trees alone (no ntuple access).

  left    the event's jets in eta-phi (radius 0.4), hard-scatter green, pileup grey,
          the tagging pair outlined; HGTD's 2.38 < |eta| < 4.0 shaded
  right   the time axis: the truth hard-scatter time, each reco t0 (TRKPTZ, WAVeS,
          HGTD) and each tagging jet's time, with the +-1 sigma the test uses (t0:
          quoted x inflation; jet: quoted x f(n), f measured on the signal's HS legs
          as in vbs_time_veto_plot.py), and each jet-vs-t0 pull with its verdict

Events are chosen by fixed rules, one per finding of the study (see CASES), after
the same recoil cuts as the slides; the rule and the event are printed.

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_veto_event_cards.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root \\
        --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root --out figs/time_veto/cards/card
"""
import argparse, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 13})   # cards are shown one per slide at ~half size
from matplotlib.patches import Circle, Rectangle
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--met-min", type=float, default=100.0)
ap.add_argument("--zpt-min", type=float, default=100.0)
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--out", required=True, help="output stem: <out>_<case>.(pdf|png)")
ap.add_argument("--set", default="findings", choices=["findings", "cases"],
                help="findings: the five one-per-finding rules (A-E); cases: the talk's case-ordered set -- "
                     "Forward-Forward (jet time compatibility vetoes where t0 keeps, and the reverse) and "
                     "Central-Forward (a pileup jet in HGTD caught, and missed), each on --bkg AND --local")
ap.add_argument("--local", help="apply every rule to THIS VBF file instead (e.g. the local sample, which keeps "
                "every truth vertex for the displays); the jet-time calibration still comes from --sig")
args = ap.parse_args()

SRC = [("trkptz", "TRKPTZ $t_0$", "#C0392B"), ("waves", "WAVeS $t_0$", "#E69F00"), ("hgtd", "HGTD $t_0$ (Athena)", "#3B8ED0")]


def load(path, col, cut):
    t = uproot.open(path)["events"]
    keys = [k for k in t.keys() if not k.startswith("jet_")]
    e_all = t.arrays(keys, library="np")
    J = t.arrays(["jet_idx", "jet_pt", "jet_eta", "jet_phi", "jet_hs", "jet_pu", "jet_t", "jet_sig", "jet_n"], library="np")
    k = e_all[col] > cut
    e = {c: v[k] for c, v in e_all.items()}
    Jk = {c: v[k] for c, v in J.items()}
    return e, Jk, e_all


S, SJ, S_all = load(args.sig, "met_truth", args.met_min)
B, BJ, _ = load(args.bkg, "z_pt", args.zpt_min)
cal = C.measure(S_all, True, "signal HS legs")   # one jet-time calibration for S and B, as the slides


def derived(e):
    d = {}
    for L in "ab":
        d[L + "_timed"] = e[L + "_core_n"] >= 1
        d[L + "_s"] = cal(e[L + "_core_n"]) * e[L + "_core_sig"]
        d[L + "_in"] = e[L + "_hgtd_ntimed"] >= 1
    for s, _, _ in SRC + [("truth", "", "")]:
        s0 = e[f"infl_{s}"] * e[f"sig0_{s}"]
        veto = np.zeros(len(s0), bool)
        for L in "ab":
            pull = np.abs(e[L + "_core_t"] - e[f"t0_{s}"]) / np.sqrt(d[L + "_s"] ** 2 + s0 ** 2)
            d[f"{L}_pull_{s}"] = np.where(d[L + "_timed"] & e[f"ok_{s}"], pull, np.nan)
            veto |= e[f"ok_{s}"] & d[L + "_timed"] & (pull >= args.nsigma)
        d["veto_" + s] = veto
    d["nin"] = d["a_in"].astype(int) + d["b_in"].astype(int)
    bothT = d["a_timed"] & d["b_timed"]
    with np.errstate(divide="ignore", invalid="ignore"):
        pj = np.abs(e["a_core_t"] - e["b_core_t"]) / np.sqrt(d["a_s"] ** 2 + d["b_s"] ** 2)
    d["pull_jj"] = np.where(bothT, pj, np.nan)
    d["veto_jj"] = bothT & (pj >= args.nsigma)
    return d


SD, BD = derived(S), derived(B)
if args.local:
    L, LJ, _ = load(args.local, "met_truth", args.met_min)
    LD = derived(L)
SAMPLES = {"sig": (S, SD, SJ), "bkg": (B, BD, BJ)}
if args.local:
    SAMPLES["local"] = (L, LD, LJ)


def pick(e, d, mask, score, rule):
    """the event in mask with the largest score (deterministic)"""
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        raise SystemExit(f"no event for rule: {rule}")
    i = idx[np.argmax(score[idx])]
    print(f"  {len(idx):,} candidates; picked entry {e['entry'][i]} of {os.path.basename(str(e['file_path'][i]))}")
    return i


def short_file(path):
    """zjets_skim_shard3of20_0.root stays; the raw VBF user.mcardiff.<task>.Output._002625.SuperNtuple.root -> raw _002625"""
    b = os.path.basename(str(path))
    if b.startswith("user.") and "._" in b:
        return "raw ntuple _" + b.split("._")[-1].split(".")[0]
    return b


def leg_role(e, L, i):
    return "HS" if e[L + "_hs"][i] else ("PU" if e[L + "_pu"][i] else "unlabelled")


# ── the cases: (key, sample, title, rule, mask, score) ───────────────────────
def build_cases(srole, brole):
    """the five rules. C (a hard-scatter pair lost) is evaluated on srole; A, B, D, E
    (pileup jets caught, a wrong t0) on brole. On the grid that is VBF and Z+jets; with
    --local, both are the local VBF file, and A drops its 'no valid Athena t0' clause
    (that clause described Z+jets: every such event a valid Athena t0 keeps has none)"""
    S, SD, _ = SAMPLES[srole]
    B, BD, _ = SAMPLES[brole]
    bname = "VBF $H\\to$inv." if brole != "bkg" else "$Z$+jets"
    bmd = "VBF" if brole != "bkg" else "Z+jets"
    out = []
    both, one = BD["nin"] == 2, BD["nin"] == 1
    hspu = (B["a_hs"] & B["b_pu"]) | (B["a_pu"] & B["b_hs"])
    a_mask = both & hspu & BD["a_timed"] & BD["b_timed"] & BD["veto_trkptz"] & BD["veto_waves"]
    if brole == "bkg":
        a_mask = a_mask & ~B["ok_hgtd"]
    out.append(("A_both_pu_caught", brole,
                f"{bname}, both tagging jets in HGTD: the pileup jet is caught by the reco $t_0$s",
                f"{bmd}, both legs in HGTD, HS+PU pair, both legs timed; vetoed by TRKPTZ and WAVeS"
                + ("; no valid HGTD t0 (every such event a valid HGTD t0 keeps has none: 5 of 5)" if brole == "bkg" else "")
                + "; highest PU-leg pT",
                a_mask, np.where(B["a_pu"], B["a_pt"], B["b_pt"])))
    pu_in = np.where(BD["a_in"], B["a_pu"], B["b_pu"])
    oth_hs_cen = np.where(BD["a_in"], B["b_hs"] & (np.abs(B["b_eta"]) < 2.4), B["a_hs"] & (np.abs(B["a_eta"]) < 2.4))
    out.append(("B_one_r2_caught", brole,
                f"{bname}, one tagging jet in HGTD: a forward pileup jet beside a central hard-scatter jet, caught",
                f"{bmd}, exactly one leg in HGTD and it is PU, the other HS with |eta| < 2.4; vetoed by TRKPTZ, WAVeS "
                "and HGTD; highest PU-leg pT",
                one & pu_in & oth_hs_cen & BD["veto_trkptz"] & BD["veto_waves"] & BD["veto_hgtd"],
                np.where(BD["a_in"], B["a_pt"], B["b_pt"])))
    sboth = (SD["nin"] == 2) & S["a_hs"] & S["b_hs"] & SD["a_timed"] & SD["b_timed"]
    fail1 = np.where(SD["a_pull_truth"] >= args.nsigma, S["a_core_n"],
                     np.where(SD["b_pull_truth"] >= args.nsigma, S["b_core_n"], 0))
    worst = np.fmax(SD["a_pull_truth"], SD["b_pull_truth"])
    out.append(("C_sig_jet_tail", srole,
                "VBF $H\\to$inv., both jets hard-scatter and in HGTD: lost to the jet-time tail, even with the true $t_0$",
                "VBF, both legs HS, in HGTD and timed; vetoed by WAVeS AND by the truth t0, the failing leg single-track; "
                "pull vs truth closest to 4",
                sboth & SD["veto_waves"] & SD["veto_truth"] & (fail1 == 1), -np.abs(worst - 4.0)))
    hshs = B["a_hs"] & B["b_hs"] & (BD["a_timed"] | BD["b_timed"])
    with np.errstate(divide="ignore", invalid="ignore"):
        t0off = np.abs(B["t0_waves"] - B["t0_truth"]) / (B["infl_waves"] * B["sig0_waves"])
    out.append(("D_bkg_wrong_t0", brole,
                f"{bname}, both tagging jets hard-scatter: vetoed only because the reco $t_0$ is wrong",
                f"{bmd}, both legs HS, >= 1 timed; vetoed by WAVeS and TRKPTZ, kept by the truth t0, WAVeS t0 > 3 sigma "
                "from truth; highest m_jj",
                hshs & BD["veto_waves"] & BD["veto_trkptz"] & ~BD["veto_truth"] & (t0off > 3.0), B["mjj"]))
    edgeA = B["a_pu"] & (np.abs(B["a_eta"]) > 4.0) & BD["a_in"] & BD["a_timed"]
    edgeB = B["b_pu"] & (np.abs(B["b_eta"]) > 4.0) & BD["b_in"] & BD["b_timed"]
    edge_pull = np.where(edgeA, BD["a_pull_waves"], BD["b_pull_waves"])
    out.append(("E_edge_pu_caught", brole,
                f"{bname}: a pileup tagging jet with its axis past $|\\eta| = 4.0$, timed by HGTD and caught",
                f"{bmd}, a PU leg with |eta_axis| > 4.0 and a timed HGTD track; vetoed by TRKPTZ and WAVeS on THAT leg; "
                "highest edge-leg |eta|",
                (edgeA | edgeB) & BD["veto_waves"] & BD["veto_trkptz"] & (edge_pull >= args.nsigma),
                np.where(edgeA, np.abs(B["a_eta"]), np.abs(B["b_eta"]))))
    return out


def build_case_set(role):
    """the talk's case-ordered set on one sample (Z+jets or the local VBF file). Within a rule
    the pair with a hard-scatter AND a pileup leg is preferred, so the decision is about a
    pileup jet; the score picks the clearest such event"""
    E, D, _ = SAMPLES[role]
    nm = {"bkg": "Z+jets", "local": "VBF (local)", "sig": "VBF"}[role]
    hspu = (E["a_hs"] & E["b_pu"]) | (E["a_pu"] & E["b_hs"])
    bothT = (D["nin"] == 2) & D["a_timed"] & D["b_timed"]
    reco_veto = D["veto_waves"] & D["veto_trkptz"]
    reco_keep = ~D["veto_waves"] & ~D["veto_trkptz"]
    t0pull = np.fmax(np.nan_to_num(D["a_pull_waves"]), np.nan_to_num(D["b_pull_waves"]))
    one = D["nin"] == 1
    pu_in = np.where(D["a_in"], E["a_pu"], E["b_pu"])
    in_timed = np.where(D["a_in"], D["a_timed"], D["b_timed"])
    pu_pt = np.where(D["a_in"], E["a_pt"], E["b_pt"])
    out = []
    def add(key, title, rule, mask, score, prefer):
        m = mask & prefer if (mask & prefer).any() else mask
        out.append((f"{key}_{role}", role, title, f"{nm}: {rule}" + ("" if (mask & prefer).any() else " (no HS+PU pair passes)"),
                    m, score))
    add("FF1_jj_vetoes", "Forward-Forward: jet time compatibility rejects where $t_0$ does not",
        "both legs in HGTD and timed; jet time compatibility vetoes, TRKPTZ and WAVeS keep; HS+PU preferred; largest jet-jet pull",
        bothT & D["veto_jj"] & reco_keep, np.nan_to_num(D["pull_jj"]), hspu)
    add("FF2_t0_vetoes", "Forward-Forward: jet time compatibility misses where $t_0$ rejects",
        "both legs in HGTD and timed; jet time compatibility keeps, TRKPTZ and WAVeS veto; HS+PU preferred; largest WAVeS pull",
        bothT & ~D["veto_jj"] & reco_veto, t0pull, hspu)
    add("CF1_caught", "Central-Forward: a pileup jet in HGTD, rejected",
        "one leg in HGTD, timed and pileup; TRKPTZ and WAVeS veto; highest pileup-leg pT",
        one & pu_in & in_timed & reco_veto, pu_pt, np.ones(len(pu_pt), bool))
    add("CF2_missed", "Central-Forward: a pileup jet in HGTD, not rejected",
        "one leg in HGTD, timed and pileup; TRKPTZ and WAVeS keep; highest pileup-leg pT",
        one & pu_in & in_timed & reco_keep, pu_pt, np.ones(len(pu_pt), bool))
    return out


if args.set == "cases":
    if not args.local:
        raise SystemExit("--set cases needs --local (one event per case from Z+jets and one from the local VBF file)")
    cases = [c for pair in zip(build_case_set("bkg"), build_case_set("local")) for c in pair]
else:
    cases = build_cases("local", "local") if args.local else build_cases("sig", "bkg")


def draw_times(at, e, d, i, names):
    """the time axis: one row per t0 (truth, TRKPTZ, WAVeS, HGTD), then one per
    tagging jet, each with the +-1 sigma the test uses, and each t0's verdict"""
    tt = e["t0_truth"][i]
    rows = [("truth $t_{\\rm HS}$", tt, 0.0, "#5B3A29", "truth")]
    for s, lab, col in SRC:
        ok = e[f"ok_{s}"][i]
        rows.append((lab + ("" if ok else "\n(no time)"), e[f"t0_{s}"][i] if ok else np.nan,
                     e[f"infl_{s}"][i] * e[f"sig0_{s}"][i] if ok else np.nan, col, s))
    jets = []
    for L, name in zip("ab", names):
        role_ = leg_role(e, L, i)
        # in HGTD: every timed ghost track (the acceptance count); time: the size of the
        # highest-sum-pT time cluster of the timed tracks within dR < 0.2 (the jet-time estimate)
        txt = f"{name}: {role_}, $\\eta$ = {e[L + '_eta'][i]:+.2f}\n{e[L + '_hgtd_ntimed'][i]} timed tracks in HGTD"
        if d[L + "_timed"][i]:
            n = e[L + '_core_n'][i]
            jets.append((txt + f"\njet time: {n}-track cluster" if n > 1 else txt + "\njet time: 1 track", e[L + "_core_t"][i], d[L + "_s"][i],
                         "#2E9E4F" if role_ == "HS" else "#555555", L))
        else:
            jets.append((txt + "\nno jet time", np.nan, np.nan, "#999999", L))
    allrows = rows + jets
    # t0 rows one unit apart; the three-line jet rows 1.6 apart
    ys = np.array([1.6 * len(jets) - 0.3 + (len(rows) - 1 - k) for k in range(len(rows))]
                  + [1.6 * (len(jets) - 1 - k) for k in range(len(jets))])
    vals = [v for _, v, _, _, _ in allrows if np.isfinite(v)]
    sigs = [x for _, v, x, _, _ in allrows if np.isfinite(v) and np.isfinite(x)]
    lo = min([v - 3 * x for _, v, x, _, _ in allrows if np.isfinite(v) and np.isfinite(x)] + [min(vals) - 40])
    hi = max([v + 3 * x for _, v, x, _, _ in allrows if np.isfinite(v) and np.isfinite(x)] + [max(vals) + 40])
    pad = 0.08 * (hi - lo)
    for y, (lab, v, sg, col, tag) in zip(ys, allrows):
        if np.isfinite(v):
            if sg > 0:   # thick: +-2 sigma; thin, capped: out to +-3 sigma (the test's threshold)
                at.plot([v - 3 * sg, v + 3 * sg], [y, y], color=col, lw=1.2, alpha=0.8, solid_capstyle="butt")
                at.plot([v - 3 * sg] * 2, [y - 0.13, y + 0.13], color=col, lw=1.2, alpha=0.8)
                at.plot([v + 3 * sg] * 2, [y - 0.13, y + 0.13], color=col, lw=1.2, alpha=0.8)
                at.plot([v - 2 * sg, v + 2 * sg], [y, y], color=col, lw=5, solid_capstyle="butt")
            at.plot([v], [y], "s" if tag in [x[0] for x in SRC] else "o", color=col, ms=8,
                    markeredgecolor="white", markeredgewidth=0.8, zorder=5)
        if tag in [x[0] for x in SRC] + ["truth"] and tag != "truth" or tag == "truth":
            if not e[f"ok_{tag}"][i]:
                verdict, vcol = "no $t_0$: kept", "#777777"
            else:
                p = [d[f"{L}_pull_{tag}"][i] for L in "ab"]
                ps = ", ".join("–" if not np.isfinite(x) else f"{x:.1f}$\\sigma$" for x in p)
                verdict = ("VETOED" if d["veto_" + tag][i] else "kept") + f"  ({ps})"
                vcol = "#B03A2E" if d["veto_" + tag][i] else "#2E7D32"
            at.text(1.015, y, verdict, transform=at.get_yaxis_transform(), ha="left", va="center", fontsize=12,
                    color=vcol, fontweight="bold")
    at.set_yticks(ys)
    at.set_yticklabels([r[0] for r in allrows], fontsize=12)
    for tl, r in zip(at.get_yticklabels(), allrows): tl.set_color(r[3])
    at.axhline(1.6 * (len(jets) - 1) + 0.85, color="#BBBBBB", lw=0.8)
    # jet time compatibility: the two jets against each other, no t0
    pj = d["pull_jj"][i]
    if np.isfinite(pj):
        jtxt, jcol = ("VETOED" if d["veto_jj"][i] else "kept") + f"  ({pj:.1f}$\\sigma$)", ("#B03A2E" if d["veto_jj"][i] else "#2E7D32")
    else:
        jtxt, jcol = "n/a: one jet timed", "#777777"
    at.text(1.015, 0.8, "jet time compatibility:\n" + jtxt, transform=at.get_yaxis_transform(), ha="left", va="center",
            fontsize=12, color=jcol, fontweight="bold")
    at.axvline(tt, color="#5B3A29", lw=0.8, ls=":")
    at.set_xlim(lo - pad, hi + pad); at.set_ylim(-0.9, ys.max() + 0.6)
    at.set_xlabel("time [ps]   (thick: $\\pm 2\\sigma$, thin: $\\pm 3\\sigma$; $\\sigma$ as used by the test)")
    at.text(1.015, ys.max() + 0.75, f"jet vs $t_0$ at {args.nsigma:g}$\\sigma$\n({names[0]}, {names[1]} pulls)",
            transform=at.get_yaxis_transform(), ha="left", va="bottom", fontsize=11, color="#555555")


def card(key, role, title, rule, mask, score):
    e, d, J = SAMPLES[role]
    print(f"[{key}] {rule}")
    i = pick(e, d, mask, score, rule)
    fig = plt.figure(figsize=(13.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.75, left=0.05, right=0.83, top=0.82, bottom=0.13)
    # ── eta-phi ──
    ax = fig.add_subplot(gs[0])
    for sgn in (-1, 1):
        ax.add_patch(Rectangle((min(sgn * 2.38, sgn * 4.0), -np.pi), 4.0 - 2.38, 2 * np.pi, color="#DCE6F2", zorder=0))
        ax.text(sgn * 3.19, -np.pi * 0.95, "HGTD", ha="center", va="bottom", fontsize=12, color="#4A6FA5")
    legs = {int(e["a_idx"][i]): "a", int(e["b_idx"][i]): "b"}
    for k in range(len(J["jet_pt"][i])):
        eta, phi, pt = J["jet_eta"][i][k], J["jet_phi"][i][k], J["jet_pt"][i][k]
        hs, pu = J["jet_hs"][i][k] > 0, J["jet_pu"][i][k] > 0
        col = "#2E9E4F" if hs else ("#8C8C8C" if pu else "#C9B458")
        leg = legs.get(int(J["jet_idx"][i][k]))
        for dphi in (0.0, 2 * np.pi, -2 * np.pi):   # wrap round phi = +-pi
            ax.add_patch(Circle((eta, phi + dphi), 0.4, facecolor=col, alpha=0.55 if leg else 0.25,
                                edgecolor="black" if leg else col, lw=2.2 if leg else 0.8, zorder=3))
        lab = (("jet 1" if leg == "a" else "jet 2") + f"\n{pt:.0f} GeV") if leg else f"{pt:.0f}"
        ty = min(max(phi, -np.pi + 0.45), np.pi - 0.45)
        ax.text(eta, ty, lab, ha="center", va="center", fontsize=11 if leg else 9, zorder=4,
                fontweight="bold" if leg else "normal")
    ax.set_xlim(-4.6, 4.6); ax.set_ylim(-np.pi, np.pi); ax.set_aspect("equal")
    ax.set_xlabel("jet $\\eta$"); ax.set_ylabel("jet $\\phi$")
    ax.set_title("green: hard scatter, grey: pileup\noutlined: tagging pair", fontsize=12)
    # ── time axis: one row per t0, then one per tagging jet ──
    draw_times(fig.add_subplot(gs[1]), e, d, i, ("jet 1", "jet 2"))
    rec = (f"truth $E_T^{{\\rm miss}}$ = {e['met_truth'][i]:.0f} GeV" if role != "bkg" else f"$p_T^{{\\ell\\ell}}$ = {e['z_pt'][i]:.0f} GeV")
    # the slide carries the title; the card keeps only the event's identity
    fig.text(0.02, 0.955, f"{'VBF H→inv.' if role == 'sig' else ('VBF H→inv. (local sample)' if role == 'local' else 'Z(→ℓℓ)+jets')}  ·  {short_file(e['file_path'][i])}  "
             f"entry {e['entry'][i]}  ·  $m_{{jj}}$ = {e['mjj'][i]:.0f} GeV  ·  {rec}  ·  ATLAS Simulation Internal",
             fontsize=12, color="#555555")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}_{key}.{ext}", dpi=150)
    plt.close(fig)
    print(f"  wrote {args.out}_{key}.pdf / .png")
    # the time panel alone, beside an event display: its jets are named as event_display.py
    # numbers them (ntuple jet index + 1)
    f2 = plt.figure(figsize=(8.6, 5.0))
    at = f2.add_axes([0.35, 0.15, 0.40, 0.72])
    draw_times(at, e, d, i, (f"Jet {int(e['a_idx'][i]) + 1}", f"Jet {int(e['b_idx'][i]) + 1}"))
    for ext in ("pdf", "png"):
        f2.savefig(f"{args.out}_{key}_times.{ext}", dpi=170, bbox_inches="tight", pad_inches=0.08)
    plt.close(f2)
    print(f"  wrote {args.out}_{key}_times.pdf / .png")
    return dict(key=key, rule=rule, file=os.path.basename(str(e["file_path"][i])), entry=int(e["entry"][i]),
                role=role, i=int(i))


picked = [card(*c) for c in cases]
with open(args.out + "_events.md", "w") as f:
    f.write("# Timing-veto event cards\n\n| card | rule | file | entry |\n|---|---|---|---|\n")
    for p in picked:
        f.write(f"| {p['key']} | {p['rule']} | {p['file']} | {p['entry']} |\n")
print("wrote", args.out + "_events.md")
# the candidates CSV vbs_veto_displays.py renders (the display's t0 is WAVeS')
import csv
COLS = ["file_path", "entry", "a_idx", "b_idx", "a_hs", "a_pu", "b_hs", "b_pu", "t0_waves", "sig0_waves", "infl_waves",
        "t0_trkptz", "t0_hgtd", "ok_hgtd", "t0_truth", "mjj"]
with open(args.out + "_candidates.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["card", "sample"] + COLS)
    w.writeheader()
    for p in picked:
        e = SAMPLES[p["role"]][0]
        row = {c: (e[c][p["i"]].item() if hasattr(e[c][p["i"]], "item") else str(e[c][p["i"]])) for c in COLS}
        w.writerow(dict(card=p["key"], sample={"sig": "vbf", "bkg": "zjets", "local": "vbf_local"}[p["role"]], **row))
print("wrote", args.out + "_candidates.csv")
