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

HERE = os.path.dirname(os.path.abspath(__file__))
TITLE = {"vbf": "VBF H->inv.", "zjets": "Z(->ll)+jets"}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("inputs", nargs="+", help="<candidates.csv>:<picked stem> (stem.root and stem_map.csv)")
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()

bad = 0
for spec in args.inputs:
    cand, stem = spec.split(":")
    remap = {(m["file_path"], int(m["entry"])): int(m["picked_entry"]) for m in csv.DictReader(open(stem + "_map.csv"))}
    for r in csv.DictReader(open(cand)):
        key = (r["file_path"], int(r["entry"]))
        lab = lambda L: "HS" if r[f"{L}_hs"] == "True" else ("PU" if r[f"{L}_pu"] == "True" else "other")
        legs, labels = f"{r['a_idx']},{r['b_idx']}", f"{lab('a')},{lab('b')}"
        t0 = float(r["t0_waves"])
        cmd = [sys.executable, "event_display.py",
               "--file_path", os.path.abspath(stem + ".root"), "--event_num", str(remap[key]),
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
