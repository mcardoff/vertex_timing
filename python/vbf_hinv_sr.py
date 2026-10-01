"""The VBF H->inv signal region of the full Run 2 analysis (ANA-EXOT-2020-11-INT1,
Sec. 4.5.2, and its Z->ll control region, Sec. 4.5.6), applied to
util/vbs_time_veto --sr outputs, and the HGTD timing variants on top of it:
python/vbs_time_veto_sr.py's docstring describes both. Shared by that script
and python/vbs_hspu_study.py.

A sample is a dict made by load(): its events `e`, jets `J` (awkward), weights
`w`, the E_T^miss proxy `met` (truth |sum pT(nu)| for VBF, pT(ll) for Z+jets,
the jet-based met_jets for multijet) and whether the Z->ll window applies."""
import re
import numpy as np
import uproot
import awkward as ak
import vbs_time_veto_common as C

# ANA-EXOT-2020-11-INT1 Sec. 4.5.2 (and 4.5.6 for the Z->ll window)
SR_JET_PT, SR_J1_PT, SR_J2_PT = 25.0, 80.0, 50.0
SR_DPHI, SR_DETA, SR_MJJ = 2.0, 3.8, 800.0
SR_ETJET = 140.0
SR_CEN, SR_MREL = 0.6, 0.05
FJVT_LOOSE, FJVT_VTIGHT, FJVT_MET = 0.5, 0.2, 200.0
MET_LO, MET_HI, MJJ_LOWBINS = 160.0, 200.0, 1500.0
MZ, MZ_WIN = 91.2, 25.0
SR_MJJ5 = [800.0, 1000.0, 1500.0, 2000.0, 3500.0]   # 2-jet, E_T^miss > 200 GeV bins (x2 in dphi)
SR_MJJ3 = [1500.0, 2000.0, 3500.0]                  # 3-4 jets, and 2 jets at 160-200 GeV

METHODS = ["jj"] + [f"t0_{s}" for s in C.T0_SRC] + [f"rm_{s}" for s in C.T0_SRC]
# The selection a timing variant is measured in; S and B are always measured in the SAME one.
#   sr       the full SR: every cut sr_select lists
#   minimal  only the cuts the background's timing response depends on (python/vbs_sr_minimal.py):
#            the tagging pair, fJVT on it, |dphi_jj| < 2, E_T^miss > 100 GeV, the Z window
REGIONS = {"sr": dict(steps=None, met_floor=None),
           "minimal": dict(steps=["j1 > 80, j2 > 50 GeV", "fJVT on j1, j2", "|dphi_jj| < 2",
                                  "eta1 eta2 < 0, |deta_jj| > 3.8", "m_jj > 800 GeV", "|m_ll - 91.2| < 25 GeV"],
                           met_floor=100.0)}
SIG, ZLL = "VBF H#rightarrowinv.", "Z(#rightarrowll)+jets"
JZ_RE = re.compile(r"\.(8011\d\d)\.Py8EG_A14NNPDF23LO_jj_(JZ\d)\.")


# ── Loading ─────────────────────────────────────────────────────────────────
def load(path, role, norm, rcol, zll):
    f = uproot.open(path)
    meta = f["meta"].arrays(library="np")
    if "sr_mode" not in meta or not meta["sr_mode"].all():
        raise SystemExit(f"{path}: not a vbs_time_veto --sr output")
    info = {k: meta[k].sum() for k in ["n_read", "n_selected", "sumw_read"]}
    for k in ["leg_jvt_rpt_min", "sr_recoil_min", "jet_time_dr", "jet_time_dist_cut"]:
        u = np.unique(meta[k])
        if len(u) != 1: raise SystemExit(f"{path}: meta {k} differs between merged jobs: {u}")
        info[k] = u[0]
    info["sample"], info["njobs"] = str(meta["sample"][0]), len(meta["n_read"])
    exp = C.EXPECTED_COMPLETE.get(info["sample"])
    if exp:
        okn = info["n_read"] == exp["n_read"]
        oks = "sumw" not in exp or abs(info["sumw_read"] / exp["sumw"] - 1) < 1e-6
        print(f"[{role}] completeness: {info['n_read']:,} events read (complete sample {exp['n_read']:,})"
              + (f", sum w {info['sumw_read']:.10g} (raw {exp['sumw']:.10g})" if "sumw" in exp else "")
              + ("  OK" if okn and oks else "  WARNING: INCOMPLETE -- a shard is missing or doubled"))
    t = f["events"]
    keys = t.keys()
    e = t.arrays([k for k in keys if not k.startswith("jet_")], library="np")
    J = t.arrays([k for k in keys if k.startswith("jet_")], library="ak")
    if np.any(e[rcol] < 0):
        raise SystemExit(f"{path}: {rcol} missing for some events (VBF: read the RAW ntuples)")
    w = e["weight"].astype(float)
    s = dict(path=path, role=role, norm=norm, info=info, e=e, J=J, w=w, met=e[rcol].astype(float), rcol=rcol,
             zll=zll)
    # fixed per event, whatever is removed: JVT decisions, E_T^jet, the Z window
    s["jvt_ok"] = ~((J.jet_jvtwin > 0) & (J.jet_rpt < info["leg_jvt_rpt_min"]))
    s["etjet"] = np.hypot(ak.to_numpy(ak.sum(J.jet_pt * np.cos(J.jet_phi), axis=1)),
                          ak.to_numpy(ak.sum(J.jet_pt * np.sin(J.jet_phi), axis=1)))
    s["zwin"] = (np.abs(e["z_mass"] - MZ) < MZ_WIN) if zll else np.ones(len(w), bool)
    return s


def subset(s, mask):
    """the sample restricted to the events in mask"""
    return dict(s, e={k: v[mask] for k, v in s["e"].items()}, J=s["J"][mask], w=s["w"][mask], met=s["met"][mask],
                jvt_ok=s["jvt_ok"][mask], etjet=s["etjet"][mask], zwin=s["zwin"][mask])


def load_samples(sig, zjets=None, dijet=None, slices="JZ2,JZ3,JZ4"):
    """{name: sample}, the signal first. The dijet JZ slices are separate samples:
    they are not normalised to each other (the ntuple weight is the generator weight)."""
    samples = {SIG: load(sig, "sig", "vbf", "met_truth", False)}
    if zjets: samples[ZLL] = load(zjets, "zjets", "zjets", "z_pt", True)
    if dijet:
        dj = load(dijet, "dijet", "none", "met_jets", False)
        tags = np.array([(lambda m: m.group(2) if m else "?")(JZ_RE.search(str(p))) for p in dj["e"]["file_path"]])
        for jz in slices.split(","):
            if (tags == jz).any(): samples[f"multijet {jz}"] = subset(dj, tags == jz)
        if not any(k.startswith("multijet") for k in samples):   # not the JZ production: one multijet sample
            samples["multijet"] = dj
    return samples


def take(arr, idx):
    """arr[event][idx[event]] as numpy; idx < 0 gives element 0 (callers mask those)"""
    j = ak.from_regular(ak.Array(np.where(idx >= 0, idx, 0).reshape(-1, 1)))
    return ak.to_numpy(ak.firsts(arr[j]))


def m2(pt1, eta1, phi1, pt2, eta2, phi2):
    return 2 * pt1 * pt2 * (np.cosh(eta1 - eta2) - np.cos(phi1 - phi2))


def sr_bin(n, met, mjj, dphi):
    """1-10: 2 jets, E_T^miss > 200 (m_jj bins x dphi < 1 | 1-2); 11-13: 3-4 jets,
    > 200, m_jj > 1.5 TeV; 14-16: 2 jets, 160-200 GeV, m_jj > 1.5 TeV; 0: none"""
    b5 = np.digitize(mjj, SR_MJJ5) - 1          # 0-4 above 0.8 TeV
    b3 = np.digitize(mjj, SR_MJJ3) - 1          # 0-2 above 1.5 TeV
    two, hi, lo = n == 2, met > MET_HI, (met > MET_LO) & (met <= MET_HI)
    out = np.zeros(len(n), int)
    m = two & hi & (b5 >= 0)
    out[m] = 1 + b5[m] + 5 * (dphi[m] >= 1.0)
    m = (n >= 3) & (n <= 4) & hi & (b3 >= 0)
    out[m] = 11 + b3[m]
    m = two & lo & (b3 >= 0)
    out[m] = 14 + b3[m]
    return out


def sr_select(s, keep, region="sr"):
    """The SR on the jets with `keep`. Returns every cut's flag (in the cut-flow
    order), the pass flag of `region` (REGIONS), the SR bin, m_jj and the j1 / j2
    indices."""
    J, met = s["J"], s["met"]
    cnt = keep & s["jvt_ok"] & (J.jet_pt > SR_JET_PT)
    loc = ak.local_index(J.jet_pt)[cnt]                # stored jets are pT-ordered
    n = ak.to_numpy(ak.num(loc))
    pad = ak.fill_none(ak.pad_none(loc, 2, axis=1), -1)
    i1, i2 = ak.to_numpy(pad[:, 0]), ak.to_numpy(pad[:, 1])
    has2 = i2 >= 0
    pt1, pt2 = take(J.jet_pt, i1), take(J.jet_pt, i2)
    e1, e2 = take(J.jet_eta, i1), take(J.jet_eta, i2)
    f1, f2 = take(J.jet_phi, i1), take(J.jet_phi, i2)
    with np.errstate(divide="ignore", invalid="ignore"):
        mjj = np.sqrt(np.maximum(m2(pt1, e1, f1, pt2, e2, f2), 0.0))
        dphi = np.abs((f1 - f2 + np.pi) % (2 * np.pi) - np.pi)
        thr = np.where(met > FJVT_MET, FJVT_LOOSE, FJVT_VTIGHT)
        fj_ok = ~((take(J.jet_fjvtwin, i1) > 0) & (take(J.jet_fjvt, i1) >= thr)) \
                & ~((take(J.jet_fjvtwin, i2) > 0) & (take(J.jet_fjvt, i2) >= thr))
        # j3, j4: per-event quantities as 1-D awkward arrays so they broadcast along the jets
        x = loc[:, 2:]
        xpt, xeta, xphi = J.jet_pt[x], J.jet_eta[x], J.jet_phi[x]
        A = lambda v: ak.Array(v)
        cen = np.exp(-4.0 / A((e1 - e2) ** 2) * (xeta - A(0.5 * (e1 + e2))) ** 2)
        mrel = np.minimum(np.sqrt(np.maximum(m2(A(pt1), A(e1), A(f1), xpt, xeta, xphi), 0.0)),
                          np.sqrt(np.maximum(m2(A(pt2), A(e2), A(f2), xpt, xeta, xphi), 0.0))) / A(mjj)
        extra_ok = ak.to_numpy(ak.all((cen < SR_CEN) & (mrel < SR_MREL), axis=1))
    bins = sr_bin(n, met, mjj, dphi)
    steps = [("E_T^miss > 160 GeV", met > MET_LO),
             ("2-4 jets > 25 GeV (JVT)", (n >= 2) & (n <= 4)),
             ("j1 > 80, j2 > 50 GeV", has2 & (pt1 > SR_J1_PT) & (pt2 > SR_J2_PT)),
             ("fJVT on j1, j2", fj_ok),
             ("j3, j4: C < 0.6, m_rel < 0.05", extra_ok),
             ("E_T^jet,no-jvt > 140 GeV", s["etjet"] > SR_ETJET),
             ("|dphi_jj| < 2", dphi < SR_DPHI),
             ("eta1 eta2 < 0, |deta_jj| > 3.8", (e1 * e2 < 0) & (np.abs(e1 - e2) > SR_DETA)),
             ("m_jj > 800 GeV", mjj > SR_MJJ),
             ("in an SR bin (E_T^miss / N_jet / m_jj)", bins > 0)]
    if s["zll"]: steps.append(("|m_ll - 91.2| < 25 GeV", s["zwin"]))
    R = REGIONS[region]
    ok = np.ones(len(n), bool)
    for name, f in steps:
        if R["steps"] is None or name in R["steps"]: ok = ok & f
    if R["met_floor"] is not None: ok = ok & (met >= R["met_floor"])
    return dict(steps=steps, ok=ok, bin=np.where(ok, bins, 0), mjj=np.where(ok, mjj, np.nan), i1=i1, i2=i2,
                n=n, mjj_raw=mjj)


# ── Timing ───────────────────────────────────────────────────────────────────
def jet_sigma(s):
    return s["calib"].jagged(s["J"].jet_n) * s["J"].jet_sig


def fails_t0(s, src, nsigma, min_trk):
    """per jet: timed, the event has this t0, and incompatible with it"""
    J, e = s["J"], s["e"]
    s0 = e[f"infl_{src}"] * e[f"sig0_{src}"]
    pull = abs(J.jet_t - e[f"t0_{src}"]) / np.sqrt(s["sig"] ** 2 + s0 ** 2)
    return (J.jet_n >= min_trk) & e[f"ok_{src}"] & (pull >= nsigma)


def outcomes(s, nsigma, min_trk, region="sr"):
    """per method: dict(ok, bin, mjj) in `region`; `none` (s["none"], made in the same region) included"""
    base = s["none"]
    out = {"none": base}
    J, i1, i2 = s["J"], base["i1"], base["i2"]
    t1, t2 = take(J.jet_t, i1), take(J.jet_t, i2)
    s1, s2 = take(s["sig"], i1), take(s["sig"], i2)
    both = (take(J.jet_n, i1) >= min_trk) & (take(J.jet_n, i2) >= min_trk)
    with np.errstate(divide="ignore", invalid="ignore"):
        jj = both & (np.abs(t1 - t2) / np.sqrt(s1 ** 2 + s2 ** 2) >= nsigma)
    out["jj"] = veto(base, jj)
    for src in C.T0_SRC:
        bad = fails_t0(s, src, nsigma, min_trk)
        out[f"t0_{src}"] = veto(base, take(bad, i1) | take(bad, i2))
        out[f"rm_{src}"] = sr_select(s, ~bad, region)
    return out


def prepare(s, calib, nsigma, min_trk, region="sr"):
    """The timing on one sample: the no-timing selection, every method's
    outcome, each t0's per-jet failures and the event veto's decision on the
    no-timing pair (for every candidate, so any cut set can use it), and the
    pair's truth composition."""
    s["calib"] = calib
    s["sig"] = jet_sigma(s)
    s["none"] = sr_select(s, ak.ones_like(s["J"].jet_pt, dtype=bool), region)
    s["out"] = outcomes(s, nsigma, min_trk, region)
    J, i1, i2 = s["J"], s["none"]["i1"], s["none"]["i2"]
    s["fail"] = {src: fails_t0(s, src, nsigma, min_trk) for src in C.T0_SRC}
    s["vetoed"] = {src: take(f, i1) | take(f, i2) for src, f in s["fail"].items()}
    h1, h2 = take(J.jet_hs, i1) > 0, take(J.jet_hs, i2) > 0
    p1, p2 = take(J.jet_pu, i1) > 0, take(J.jet_pu, i2) > 0
    s["pu_leg"] = (p1, p2)
    s["comp"] = {"HS+HS": h1 & h2, "HS+PU": (h1 & p2) | (p1 & h2), "PU+PU": p1 & p2}
    s["comp"]["other"] = ~(s["comp"]["HS+HS"] | s["comp"]["HS+PU"] | s["comp"]["PU+PU"])
    s["pu_tag"] = s["comp"]["HS+PU"] | s["comp"]["PU+PU"]
    return s


def veto(base, rejected):
    ok = base["ok"] & ~rejected
    return dict(ok=ok, bin=np.where(ok, base["bin"], 0), mjj=np.where(ok, base["mjj"], np.nan))


def ratio(w, a, b):
    """sum w over mask a / over mask b, and its error from the events the two do
    not share (those in both cancel exactly)"""
    yb = w[b].sum()
    if yb <= 0: return np.nan, np.nan
    return w[a].sum() / yb, np.sqrt((w[a ^ b] ** 2).sum()) / yb
