#!/usr/bin/env python3
"""WHY each displayed timing-veto event's times are where they are.

For every event of the candidates CSV(s) (python/vbs_veto_event_cards.py) this
  * recovers the tracks behind each tagging jet's time: the subset of the jet's
    time candidates (timed ghost tracks passing the track selection in HGTD within
    dR < JET_TIME_DR) of size core_n whose inverse-variance mean IS the tree's
    core_t -- so the cluster is identified exactly, not re-clustered;
  * reports per track: pT, time and resolution, HGTD hits (primary hits), its truth
    vertex (0 = hard scatter; grid samples link only the hard scatter, so -1 there
    means pileup or unlinked) and the truth time of its particle / vertex, and the
    pull of its HGTD time against that truth time (a large pull with a matched track
    = a mis-assigned HGTD hit);
  * parses the event display's log for the t0 clusters and says which one TRKPTZ
    and WAVeS selected and how hard-scatter it is.

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_veto_event_why.py \\
        figs/time_veto/case_cards/zjets_case_candidates.csv:condor/zjets/zjets_case_display_events:condor/zjets/zjets_jvtLoose_vbs_time_veto.root \\
        figs/time_veto/case_cards/local_case_candidates.csv:-:figs/hists/jvtLoose_vbs_time_veto.root \\
        --logs figs/time_veto/case_displays --out figs/time_veto/case_cards/why.md
"""
import argparse, csv, itertools, os, re
import numpy as np
import uproot

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("inputs", nargs="+", help="<candidates.csv>:<picked stem or ->:<vbs_time_veto tree>")
ap.add_argument("--logs", required=True, help="directory of the event displays' .log files")
ap.add_argument("--out", required=True)
args = ap.parse_args()

BR = ["AntiKt4EMTopoJets_pt", "AntiKt4EMTopoJets_eta", "AntiKt4EMTopoJets_phi", "AntiKt4EMTopoJets_ghostTrack_idx",
      "Track_pt", "Track_eta", "Track_phi", "Track_quality", "Track_time", "Track_timeRes", "Track_hasValidTime",
      "Track_nHGTDHits", "Track_nHGTDPrimaryHits", "Track_truthVtx_idx", "Track_truthPart_idx", "Track_z0",
      "TruthVtx_time", "TruthVtx_z", "TruthPart_prodVtx_time"]


def truth_time(g, k):
    """(label, truth time) of track k: its particle's production time, else its vertex's"""
    tv, tp = g["Track_truthVtx_idx"][k], g["Track_truthPart_idx"][k]
    tpt = g["TruthPart_prodVtx_time"]
    t = float(tpt[tp]) if 0 <= tp < len(tpt) else (float(g["TruthVtx_time"][tv]) if tv >= 0 else np.nan)
    lab = "HS" if tv == 0 else (f"PU vtx {tv}" if tv > 0 else "no truth link")
    return lab, t


def jet_cluster(g, j, core_t, core_n, dr_max):
    je, jp = g["AntiKt4EMTopoJets_eta"][j], g["AntiKt4EMTopoJets_phi"][j]
    cand = []
    for k in g["AntiKt4EMTopoJets_ghostTrack_idx"][j]:
        if g["Track_hasValidTime"][k] != 1 or not g["Track_timeRes"][k] > 0 or not g["Track_quality"][k]:
            continue
        if not (2.38 < abs(g["Track_eta"][k]) < 4.0 and 1.0 < g["Track_pt"][k] < 30.0):
            continue
        if np.hypot(g["Track_eta"][k] - je, (g["Track_phi"][k] - jp + np.pi) % (2 * np.pi) - np.pi) >= dr_max:
            continue
        cand.append(int(k))
    best, bd = None, 1e9
    for sub in itertools.combinations(cand, core_n):
        w = np.array([1 / g["Track_timeRes"][k] ** 2 for k in sub])
        t = float(np.sum(w * np.array([g["Track_time"][k] for k in sub])) / w.sum())
        if abs(t - core_t) < bd:
            best, bd = sub, abs(t - core_t)
    return best, bd, cand


def t0_clusters(log):
    """[(t, score_trkptz, score_waves, [track idx])] from an event display log"""
    out, cur = [], None
    for ln in open(log):
        ln = ln.strip()
        if ln.startswith("t: "):
            cur = [float(ln[3:]), None, None, []]
        elif cur is not None and ln.startswith("score_trkptz:"):
            cur[1] = float(ln.split(":")[1])
        elif cur is not None and ln.startswith("score_waves:"):
            cur[2] = float(ln.split(":")[1])
        elif cur is not None and re.match(r"^\d+,", ln):
            cur[3].append(int(ln.split(",")[0]))
        elif cur is not None and ln.startswith("passes?"):
            out.append(tuple(cur)); cur = None
    return out


md = ["# Why each displayed event's times are where they are", ""]
for spec in args.inputs:
    cand_csv, stem, tree = spec.split(":")
    T = uproot.open(tree)["events"].arrays(library="np", filter_name=lambda k: not k.startswith("jet_"))
    dr_max = float(uproot.open(tree)["meta"].arrays(["jet_time_dr"], library="np")["jet_time_dr"][0])
    remap = ({(m["file_path"], int(m["entry"])): int(m["picked_entry"]) for m in csv.DictReader(open(stem + "_map.csv"))}
             if stem != "-" else None)
    for r in csv.DictReader(open(cand_csv)):
        key = (r["file_path"], int(r["entry"]))
        src, ent = (stem + ".root", remap[key]) if remap else key
        g = {k: v[0] for k, v in uproot.open(src)["ntuple"].arrays(BR, library="np", entry_start=ent, entry_stop=ent + 1).items()}
        i = np.flatnonzero((T["file_path"] == r["file_path"]) & (T["entry"] == int(r["entry"])))[0]
        tHS = float(T["t0_truth"][i])
        md += [f"## {r['card']}  ({os.path.basename(r['file_path'])} #{r['entry']}; t_HS = {tHS:.1f} ps)", ""]
        for L in "ab":
            j, n, t = int(T[L + "_idx"][i]), int(T[L + "_core_n"][i]), float(T[L + "_core_t"][i])
            role = "HS" if T[L + "_hs"][i] else ("PU" if T[L + "_pu"][i] else "unlabelled")
            head = (f"**Jet {j + 1}** ({role}, eta {T[L + '_eta'][i]:+.2f}, pT {T[L + '_pt'][i]:.0f} GeV, "
                    f"{T[L + '_hgtd_ntimed'][i]} timed HGTD tracks)")
            if n == 0:
                md.append(head + ": no jet time"); continue
            sub, bd, cand = jet_cluster(g, j, t, n, dr_max)
            md.append(head + f": time {t:.1f} ps from {n} of {len(cand)} core candidates"
                      + ("" if bd < 0.01 else f"  (closest subset off by {bd:.2f} ps -- CHECK)"))
            md.append("")
            md.append("| track | pT | time | res | HGTD hits (primary) | truth | truth time | pull vs own truth |")
            md.append("|---|---|---|---|---|---|---|---|")
            for k in sorted(sub, key=lambda k: -g["Track_pt"][k]):
                lab, tt = truth_time(g, k)
                pull = (g["Track_time"][k] - tt) / g["Track_timeRes"][k] if np.isfinite(tt) else np.nan
                md.append(f"| {k} | {g['Track_pt'][k]:.1f} | {g['Track_time'][k]:.1f} | {g['Track_timeRes'][k]:.1f} | "
                          f"{g['Track_nHGTDHits'][k]} ({g['Track_nHGTDPrimaryHits'][k]}) | {lab} | "
                          f"{tt:.1f} | {pull:+.1f} |")
            others = [k for k in cand if k not in sub]
            if others:
                md.append("")
                md.append("other core candidates (not in the time cluster): " + ", ".join(
                    f"{k}: {g['Track_time'][k]:.0f} ps {truth_time(g, k)[0]}" for k in others))
            md.append("")
        log = os.path.join(args.logs, r["card"] + ".log")
        if os.path.exists(log):
            cl = t0_clusters(log)
            for name, si in (("TRKPTZ", 1), ("WAVeS", 2)):
                c = max(cl, key=lambda x: x[si])
                hs = [k for k in c[3] if g["Track_truthVtx_idx"][k] == 0]
                pt = np.array([g["Track_pt"][k] for k in c[3]])
                hsf = sum(g["Track_pt"][k] for k in hs) / pt.sum() if pt.sum() > 0 else 0
                md.append(f"- {name} cluster: t = {c[0]:.1f} ps, {len(c[3])} tracks, {len(hs)} hard-scatter "
                          f"({100 * hsf:.0f}% of its pT); tree t0_{name.lower()} = {float(T['t0_' + name.lower()][i]):.1f} ps")
        md.append("")
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
open(args.out, "w").write("\n".join(md) + "\n")
print("\n".join(md))
