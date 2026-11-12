#!/usr/bin/env python3
"""The event displays behind the Global Performance talk's failure-mode backup slides (slides/slides_gp_timing.tex,
linked from "Where the t0 fails"): one local-VBF event per WAVeS failure category, drawn by event_display.py (page 2:
the R-Z view and the time histogram), with the tracks that make the event belong to its category highlighted.

Inputs: python/clustering_dump.cxx JSONs, figs/tdr/gp/failures/gp_fail_<mode>.json (the talk's processing: TDR-cleaned
times, 2.5 sigma z association). The events were chosen from util/failure_decomposition's local tree
(figs/hists/failure_decomposition_local_talk.root), classified with the talk ladder (python/gp_failure_ladder.py
--no-ideal), each the clearest of its category; notes/gp_talk.md has the selection.

Highlighted, per mode:
  retime     the chosen cluster's tracks within dR < 0.4 of a forward jet: the WAVeS t0 is built from these alone
  badtime    the hard-scatter tracks >= 3 sigma from the truth vertex time
  selection  the chosen (pileup) cluster's tracks, and the hard scatter's own cluster's (most HS pT, > 50% HS)
  miscl      nothing: the display's solid / dashed tracks already say which of the chosen cluster's are pileup
  nousable   nothing: the display shows it (no hard-scatter track has a time)

  cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_failure_event_displays.py
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "figs", "tdr", "gp", "failures")
TITLE = {"retime": "re-timing", "badtime": "misassigned HS track times", "selection": "selection",
         "miscl": "misclustering", "nousable": "no usable HS cluster"}


def highlights(m, E):
    C = E["clusters"]; c = C[E["i_waves"]]
    if m == "retime":
        return [("chosen cluster, in a jet: the WAVeS t0", [t["idx"] for t in c["tracks"] if t["injet"]])]
    if m == "badtime":
        return [("HS track, time ≥ 3σ from $t_{\\rm HS}$", [t["idx"] for x in C for t in x["tracks"] if t["mistimed"]])]
    if m == "selection":
        k = int(np.argmax([x["purity"] * x["sumpt"] for x in C]))
        return [("chosen cluster (pileup)", [t["idx"] for t in c["tracks"]]),
                ("the hard scatter's cluster", [t["idx"] for t in C[k]["tracks"]])]
    return []   # miscl, nousable (user, 2026-11-12: the display already shows it)


for m in TITLE:
    E = json.load(open(os.path.join(D, f"gp_fail_{m}.json")))
    hl = highlights(m, E)
    cmd = [sys.executable, "event_display.py", "--file_path", E["file"], "--event_num", str(E["entry"]),
           "--extra_time", repr(E["t0_waves"]), "--assoc_nsigma", "2.5", "--clean_clusters", "--clean_times",
           "--t0_label", "WAVeS t0",
           "--case_label", f"WAVeS t0 = {E['t0_waves']:.1f} ps   (truth t_HS = {E['t_hs']:.1f} ps)   |   failure: {TITLE[m]}",
           "--source_label", f"VBF H->inv. (local)   {os.path.basename(E['file'])} #{E['entry']}",
           "--output_dir", os.path.abspath(D), "--output_name", f"gp_evd_{m}"]
    if any(idx for _, idx in hl):
        cmd += ["--highlight_tracks", ";".join(f"{lab}={','.join(map(str, idx))}" for lab, idx in hl if idx)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE)
    with open(os.path.join(D, f"gp_evd_{m}.log"), "w") as f:
        f.write(" ".join(cmd) + "\n\n" + r.stdout[-4000:] + r.stderr[-3000:])
    if r.returncode:
        sys.exit(f"{m}: event_display.py failed (rc {r.returncode}), see gp_evd_{m}.log")
    subprocess.run([os.path.join(HERE, "pdfpage2png.sh"), os.path.join(D, f"gp_evd_{m}.pdf"), "2",
                    os.path.join(D, f"gp_evd_{m}.png")], check=True, capture_output=True)
    print(f"{m:9s} {os.path.basename(E['file'])} #{E['entry']}: t0 - t_HS {E['t0_waves'] - E['t_hs']:7.1f} ps; "
          + "; ".join(f"{lab}: {len(idx)}" for lab, idx in hl))
