"""The offline TDR VBF preselection on util/tdr_vbf.cxx's trees, shared by tdr_vbf_plots.py (ROOT) and
tdr_vbf_bdt.py (xgboost): numpy / uproot only, so either environment can import it.

The preselection (TDR p. 52): drop the jets failing the pileup tagger `col` (None = no tagging), then the
two leading surviving jets need pT > 75 / > 50 GeV and |Delta eta| > 3. One routine (_pair) finds that pair,
and both select() (the stacks and figures) and pair_features() (the BDT inputs) use it, so the two can never
disagree on which events pass or which jets form the pair.
"""
import json

import numpy as np
import uproot

J1, J2, DETA = 75.0, 50.0, 3.0


def load_wp(path):
    return json.load(open(path))


def load(path, WP):
    """the events tree, after checking its gate / cleaning / matching against the working points'"""
    t = uproot.open(path)
    E = t["events"].arrays(library="np")
    meta = t["meta"].arrays(library="np")
    g = float(np.unique(meta["gate_sigma"])[0])
    cl = str(meta["time_cleaning"][0]) if "time_cleaning" in meta else "off"
    mt = str(meta["rpt_match"][0]) if "rpt_match" in meta else "dr02"
    if mt != WP.get("rpt_match", "dr02"):
        raise SystemExit(f"{path}: R_pT matching '{mt}', but the working points were made with '{WP.get('rpt_match', 'dr02')}'")
    if abs(g - WP["gate_sigma"]) > 1e-9:
        raise SystemExit(f"{path}: gate {g} sigma, but the working points were made at {WP['gate_sigma']}")
    if cl != WP.get("time_cleaning", "off"):
        raise SystemExit(f"{path}: track-time cleaning '{cl}', but the working points were made with '{WP.get('time_cleaning', 'off')}'")
    return E


def thresholds(WP, col):
    w = WP["working_points"][col]
    return {(r, s): w[f"{r}_{s}"]["threshold"] for r in ("central", "forward") for s in ("30_50", "50_inf")}


def _surviving(E, i, th, col):
    """indices of event i's jets passing the tagger `col` (all of them for col = None), pT-ordered"""
    pt, eta = E["jet_pt"][i], E["jet_eta"][i]
    if not col:
        return np.arange(len(pt))
    r = E["jet_" + col][i]
    ae = np.abs(eta)
    thr = np.where(ae < 2.4, np.where(pt < 50, th["central", "30_50"], th["central", "50_inf"]),
                   np.where(pt < 50, th["forward", "30_50"], th["forward", "50_inf"]))
    return np.flatnonzero((r < 0) | (r >= thr))               # r < 0: beyond |eta| = 4.0, no tracks


def _pair(E, i, th, col):
    """(leading, second) jet indices of the preselected pair, or None; plus the surviving jets"""
    s = _surviving(E, i, th, col)                             # the tree's jets are pT-ordered
    if len(s) < 2:
        return None, s
    pt, eta = E["jet_pt"][i], E["jet_eta"][i]
    a, b = s[0], s[1]
    if pt[a] <= J1 or pt[b] <= J2 or abs(eta[a] - eta[b]) <= DETA:
        return None, s
    return (a, b), s


def _mjj(pt, eta, phi, a, b):
    e = pt[a] * np.cosh(eta[a]) + pt[b] * np.cosh(eta[b])
    px = pt[a] * np.cos(phi[a]) + pt[b] * np.cos(phi[b]); py = pt[a] * np.sin(phi[a]) + pt[b] * np.sin(phi[b])
    pz = pt[a] * np.sinh(eta[a]) + pt[b] * np.sinh(eta[b])
    return np.sqrt(max(e * e - px * px - py * py - pz * pz, 0.0))


def select(E, WP, col):
    """per event: passes, m_jj, category, pileup jet in the pair, one beyond |eta| = 4.0, how many of the pair
    are hard-scatter / pileup by the TDR labels, and the pair's pileup legs by where they are -- after
    dropping the jets failing `col`"""
    n = len(E["weight"])
    ok = np.zeros(n, bool); mjj = np.zeros(n); cat = np.full(n, -1); pupair = np.zeros(n, bool); pu40 = np.zeros(n, bool)
    nhs = np.zeros(n, int); npu = np.zeros(n, int)
    pu_cen = np.zeros(n, int); pu_reach = np.zeros(n, int); pu_far = np.zeros(n, int)
    th = thresholds(WP, col) if col else None
    for i in range(n):
        p, _ = _pair(E, i, th, col)
        if p is None: continue
        a, b = p
        pt, eta, phi = E["jet_pt"][i], E["jet_eta"][i], E["jet_phi"][i]
        ok[i] = True
        mjj[i] = _mjj(pt, eta, phi, a, b)
        cat[i] = int(abs(eta[a]) > 2.4) + int(abs(eta[b]) > 2.4)
        pupair[i] = bool(E["jet_pu"][i][a] or E["jet_pu"][i][b])
        pu40[i] = any(E["jet_pu"][i][j] and abs(eta[j]) >= 4.0 for j in (a, b))
        nhs[i] = int(E["jet_hs"][i][a]) + int(E["jet_hs"][i][b])
        npu[i] = int(E["jet_pu"][i][a]) + int(E["jet_pu"][i][b])
        for j in (a, b):                                      # the pair's pileup legs, by where they are
            if not E["jet_pu"][i][j]: continue
            aej = abs(eta[j])
            if aej < 2.4: pu_cen[i] += 1                      # central: ITk tracks, no HGTD
            elif aej < 4.0: pu_reach[i] += 1                  # forward, inside the tracker and HGTD
            else: pu_far[i] += 1                              # beyond |eta| = 4.0: no tracks at all
    return ok, mjj, cat, pupair, pu40, nhs, npu, pu_cen, pu_reach, pu_far


FEATURES = ["mjj", "deta", "dphi", "pt1", "pt2", "aeta1", "aeta2", "eta_prod", "n_jets", "recoil"]


def pair_features(E, WP, col):
    """the BDT inputs of the preselected pair (NaN where the event fails): m_jj, |Delta eta|, |Delta phi|, both
    jets' pT and |eta|, eta1 * eta2 (< 0: opposite hemispheres), the surviving jets with pT > 30 GeV, and the
    recoil -- the Z pT for Z+jets (its leptons stand in for Z->nunu) and the truth MET otherwise"""
    n = len(E["weight"])
    F = {k: np.full(n, np.nan) for k in FEATURES}
    th = thresholds(WP, col) if col else None
    recoil = np.where(E["z_pt"] >= 0, E["z_pt"], E["met_truth"]).astype(float)
    recoil = np.where(recoil >= 0, recoil, np.nan)
    for i in range(n):
        p, s = _pair(E, i, th, col)
        if p is None: continue
        a, b = p
        pt, eta, phi = E["jet_pt"][i], E["jet_eta"][i], E["jet_phi"][i]
        F["mjj"][i] = _mjj(pt, eta, phi, a, b)
        F["deta"][i] = abs(eta[a] - eta[b])
        F["dphi"][i] = abs((phi[a] - phi[b] + np.pi) % (2 * np.pi) - np.pi)
        F["pt1"][i], F["pt2"][i] = pt[a], pt[b]
        F["aeta1"][i], F["aeta2"][i] = abs(eta[a]), abs(eta[b])
        F["eta_prod"][i] = eta[a] * eta[b]
        F["n_jets"][i] = len(s)
        F["recoil"][i] = recoil[i]
    return F
