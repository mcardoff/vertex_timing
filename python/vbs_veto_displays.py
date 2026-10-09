#!/usr/bin/env python3
"""Standard event displays (event_display.py) of the timing-veto example events.

Reads the candidates CSV(s) written for the veto events (one row per event:
card, sample, file_path, entry, the tagging pair a_idx / b_idx with their truth
labels, and each t0 from the util/vbs_time_veto tree) and the files
`rpt_region_displays.py pick` copied the events into (with their _map.csv), and
renders each with the WAVeS t0 as the display's t0 -- its value, quoted sigma
and inflation exactly as the veto used them -- and the tagging pair as the legs,
coloured by their truth labels.

    # AF (where the ntuples are), once per sample:
    python3 python/rpt_region_displays.py pick condor/<s>_veto_display_candidates.csv \\
        --out condor/<s>_veto_display_events
    # laptop, from python/:
    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python vbs_veto_displays.py \\
        ../figs/time_veto/displays/vbf_veto_display_candidates.csv:../condor/vbf/vbf_veto_display_events \\
        ../figs/time_veto/displays/zjets_veto_display_candidates.csv:../condor/zjets/zjets_veto_display_events \\
        --out-dir ../figs/time_veto/displays
"""
import argparse, csv, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TITLE = {"vbf": "VBF H->inv.", "zjets": "Z(->ll)+jets", "vbf_local": "VBF H->inv. (local sample)"}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("inputs", nargs="+", help="<candidates.csv>:<picked stem> (stem.root and stem_map.csv), or "
                "<candidates.csv> alone to read each event from its own file_path (files that are on this machine)")
ap.add_argument("--veto-tree", help="the vbs_time_veto output the events came from: work out each event's "
                "analysis-removed jets (tagged_out) when the CSV has no such column")
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()

def tagged_out(veto_tree, fpath, entry, src, ent):
    """jets >= 30 GeV of the event that the vbs_time_veto tree's jet_* arrays (the jets left
    after the lepton overlap removal and JVT/fJVT) do not hold: within dR < 0.5 of a selected
    lepton (Track_leptonID, pT > 20 GeV) it is the overlap removal, otherwise JVT/fJVT"""
    import uproot
    e = uproot.open(veto_tree)["events"].arrays(["file_path", "entry", "jet_idx", "n_rm_jvt", "n_rm_fjvt"], library="np")
    i = np.flatnonzero((e["file_path"] == fpath) & (e["entry"] == entry))[0]
    n = uproot.open(src)["ntuple"]
    want = ["AntiKt4EMTopoJets_pt", "AntiKt4EMTopoJets_eta", "AntiKt4EMTopoJets_phi", "Track_pt", "Track_eta", "Track_phi"]
    lep = "Track_leptonID" in n.keys()
    a = n.arrays(want + (["Track_leptonID"] if lep else []), library="np", entry_start=ent, entry_stop=ent + 1)
    g = lambda k: a[k][0]
    kept = set(int(x) for x in e["jet_idx"][i])
    tags = []
    for k in range(len(g("AntiKt4EMTopoJets_pt"))):
        if g("AntiKt4EMTopoJets_pt")[k] < 30 or k in kept:
            continue
        near = lep and any(g("Track_leptonID")[t] != 0 and g("Track_pt")[t] > 20 and
                           np.hypot(g("AntiKt4EMTopoJets_eta")[k] - g("Track_eta")[t],
                                    (g("AntiKt4EMTopoJets_phi")[k] - g("Track_phi")[t] + np.pi) % (2 * np.pi) - np.pi) < 0.5
                           for t in range(len(g("Track_leptonID"))))
        tags.append(f"{k}:{'lepton' if near else 'JVT/fJVT'}")
    njvt = sum(t.endswith("JVT/fJVT") for t in tags)
    if njvt != e["n_rm_jvt"][i] + e["n_rm_fjvt"][i]:
        sys.exit(f"{fpath}#{entry}: {njvt} jets attributed to JVT/fJVT, the tree counts {e['n_rm_jvt'][i] + e['n_rm_fjvt'][i]}")
    return ",".join(tags)


bad = 0
for spec in args.inputs:
    cand, stem = spec.split(":") if ":" in spec else (spec, None)
    remap = ({(m["file_path"], int(m["entry"])): int(m["picked_entry"]) for m in csv.DictReader(open(stem + "_map.csv"))}
             if stem else None)
    for r in csv.DictReader(open(cand)):
        key = (r["file_path"], int(r["entry"]))
        src, ent = (os.path.abspath(stem + ".root"), remap[key]) if stem else (r["file_path"], int(r["entry"]))
        if "tagged_out" not in r and args.veto_tree:
            r["tagged_out"] = tagged_out(args.veto_tree, r["file_path"], int(r["entry"]), src, ent)
        lab = lambda L: "HS" if r[f"{L}_hs"] == "True" else ("PU" if r[f"{L}_pu"] == "True" else "other")
        legs, labels = f"{r['a_idx']},{r['b_idx']}", f"{lab('a')},{lab('b')}"
        t0 = float(r["t0_waves"])
        cmd = [sys.executable, "event_display.py",
               "--file_path", src, "--event_num", str(ent),
               "--extra_time", repr(t0), "--t0_sigma", r["sig0_waves"], "--infl", r["infl_waves"],
               "--gate_sigma", "3.0", "--assoc_nsigma", "2.5",
               "--legs", legs, "--leg_labels", labels,
               "--case_label", f"WAVeS t0 = {t0:.1f} ps   (truth t_HS = {float(r['t0_truth']):.1f} ps)",
               "--t0_label", "WAVeS t0",
               "--source_label", f"{TITLE[r['sample']]}   {os.path.basename(r['file_path'])} #{r['entry']}   |   JVT + fJVT loose",
               "--output_dir", os.path.abspath(args.out_dir), "--output_name", r["card"]]
        if r.get("tagged_out"):   # jets the analysis removed (JVT/fJVT, lepton overlap): not drawn
            cmd += ["--tagged_out", r["tagged_out"], "--drop_tagged"]
        os.makedirs(args.out_dir, exist_ok=True)
        p = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
        with open(os.path.join(args.out_dir, r["card"] + ".log"), "w") as f:
            f.write(" ".join(cmd) + "\n\n" + p.stdout + "\n" + p.stderr)
        pdf = os.path.join(args.out_dir, r["card"] + ".pdf")
        ok = p.returncode == 0 and os.path.exists(pdf)
        if ok:
            for page in (1, 2, 3):
                subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "110", "-png", "-singlefile",
                                pdf, os.path.join(args.out_dir, f"{r['card']}_p{page}")], check=False)
        notes = [ln for ln in p.stdout.splitlines() if ln.startswith("LEG_LABEL_NOTE")]
        print(f"{r['card']}: {'ok' if ok else 'FAILED rc=%d' % p.returncode}" + (f"  {notes}" if notes else ""))
        bad += not ok
sys.exit(1 if bad else 0)
