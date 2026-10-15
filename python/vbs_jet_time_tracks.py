"""The tracks a tagging jet's time is built from, recovered exactly from a util/vbs_time_veto
output: the subset of the jet's time candidates (timed ghost tracks passing the track selection
in HGTD within dR < jet_time_dr) of size core_n whose inverse-variance mean IS the tree's core_t.
Shared by vbs_veto_displays.py (drawn on the event displays) and vbs_case_schematic.py."""
import sys
import numpy as np


def jet_time_tracks(veto_tree, fpath, entry, src, ent):
    """{jet index: [track indices]} -- the tracks each tagging leg's time is built from: the
    subset of the leg's time candidates (timed ghost tracks passing the track selection in
    HGTD within dR < jet_time_dr) of size core_n whose inverse-variance mean IS core_t"""
    import itertools, uproot
    T = uproot.open(veto_tree)
    e = T["events"].arrays(["file_path", "entry", "a_idx", "b_idx", "a_core_t", "b_core_t", "a_core_n", "b_core_n"],
                           library="np")
    drmax = float(T["meta"].arrays(["jet_time_dr"], library="np")["jet_time_dr"][0])
    i = np.flatnonzero((e["file_path"] == fpath) & (e["entry"] == entry))[0]
    br = ["AntiKt4EMTopoJets_eta", "AntiKt4EMTopoJets_phi", "AntiKt4EMTopoJets_ghostTrack_idx", "Track_pt", "Track_eta",
          "Track_phi", "Track_quality", "Track_time", "Track_timeRes", "Track_hasValidTime"]
    g = {k: v[0] for k, v in uproot.open(src)["ntuple"].arrays(br, library="np", entry_start=ent, entry_stop=ent + 1).items()}
    out = {}
    for L in "ab":
        j, n, t = int(e[L + "_idx"][i]), int(e[L + "_core_n"][i]), float(e[L + "_core_t"][i])
        if j < 0 or n == 0: continue
        je, jp = g["AntiKt4EMTopoJets_eta"][j], g["AntiKt4EMTopoJets_phi"][j]
        cand = [int(k) for k in g["AntiKt4EMTopoJets_ghostTrack_idx"][j]
                if g["Track_hasValidTime"][k] == 1 and g["Track_timeRes"][k] > 0 and g["Track_quality"][k]
                and 2.38 < abs(g["Track_eta"][k]) < 4.0 and 1.0 < g["Track_pt"][k] < 30.0
                and np.hypot(g["Track_eta"][k] - je, (g["Track_phi"][k] - jp + np.pi) % (2 * np.pi) - np.pi) < drmax]
        best, bd = None, 1e9
        for sub in itertools.combinations(cand, n):
            w = np.array([1 / g["Track_timeRes"][k] ** 2 for k in sub])
            tt = float(np.sum(w * np.array([g["Track_time"][k] for k in sub])) / w.sum())
            if abs(tt - t) < bd: best, bd = sub, abs(tt - t)
        if bd > 0.01:
            sys.exit(f"{fpath}#{entry}: jet {j} time cluster not recovered (closest subset off by {bd:.3f} ps)")
        out[j] = list(best)
    return out

