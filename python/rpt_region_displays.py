#!/usr/bin/env python3
"""Event displays of the VBS-region events whose R_pT changes most under each
timing case, each display reproducing that case's gate exactly.

    select  rank rpt_v5_hist's region tree per case -- the R1 events whose
            HS-minus-PU R_pT margin moves most (either sign), and the R2
            events whose forward PU leg drops most -- and write a CSV holding
            every argument a display needs (t0, sigma_t0, inflation, gate,
            legs, labels, the analysis' own R_pT values).
    pick    run where the ntuples are (the skims live only on the AF): copy
            the selected events into one small file plus a mapping CSV, so
            the displays can be rendered anywhere.
    render  run event_display.py for every candidate. Each display
            recomputes its two legs' R_pT and must match the tree
            (RPT_CHECK). The idealised cases must also replay the same t0
            from src/idealised_timing.h (REPLAY_CHECK). Any mismatch fails
            the run.

PyROOT only -- no uproot -- so `pick` runs under the AF's `lsetup root`
python. Laptop:

    export PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root
    PY=~/.venv-hgtd/bin/python
    $PY python/rpt_region_displays.py select condor/vbf/vbf_rpt_v5_regions.root --sample vbf
    # AF:  python3 python/rpt_region_displays.py pick <candidates.csv>
    $PY python/rpt_region_displays.py render <candidates.csv> \\
        --picked vbf_region_display_events.root --pick-map vbf_region_display_events_map.csv
"""
import argparse, csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))

SCEN = ["zonly", "hgtd", "trkptz", "waves", "waves_ideal", "truth", "tzp"]
# (scenario row, short name for paths, display title, ideal_t0 for REPLAY_CHECK)
CASES = [
    ("trkptz",      "trkptz",      "TRKPTZ t0",                   None),
    ("waves",       "waves",       "WAVeS t0",                    None),
    ("hgtd",        "hgtd",        "HGTD t0 (Athena)",            None),
    ("truth",       "truth",       "Truth vertex t0",             "truth"),
    ("waves_ideal", "ideal_trk",   "Ideal track time assignment", "cluster"),
]
FIELDS = ["sample", "case", "case_title", "region", "rank", "tag", "delta",
          "file_path", "entry", "idx_hs", "idx_pu", "hs_pt", "hs_eta", "pu_pt", "pu_eta",
          "rpt_hs_z", "rpt_hs_t", "rpt_pu_z", "rpt_pu_t",
          "t0", "sig0", "infl", "gate_sigma", "ideal_times", "ideal_t0"]


def cmd_select(args):
    cols = (["file_path", "entry", "region", "idx_hs", "idx_pu", "hs_pt", "hs_eta",
             "pu_pt", "pu_eta", "n_jets_fwd_acc", "gate_sigma", "rpt_dzpara"]
            + [f"rpt_{leg}_{s}" for leg in ("hs", "pu") for s in SCEN]
            + [f"{q}_{s}" for q in ("t0", "sig0", "infl", "ok") for s in SCEN])
    a = ROOT.RDataFrame("regions", args.file).AsNumpy(cols)
    fp = [str(x) for x in a["file_path"]]
    n = len(fp)
    if n == 0:
        sys.exit(f"{args.file}: the regions tree is empty")
    if not all(bool(x) for x in a["rpt_dzpara"]):
        # event_display.py only implements the getNewDzpara association
        sys.exit("region tree was made with --rpt-signif; the display cannot reproduce it")
    keep = [True] * n
    if args.require_fwd_acc_jet:
        keep = [int(x) >= 1 for x in a["n_jets_fwd_acc"]]

    rows = []
    for s, short, title, ideal_t0 in CASES:
        for region in (1, 2):
            cand = []
            for i in range(n):
                if not keep[i] or int(a["region"][i]) != region:
                    continue
                zpu, tpu = float(a["rpt_pu_zonly"][i]), float(a[f"rpt_pu_{s}"][i])
                if region == 1:
                    zhs, ths = float(a["rpt_hs_zonly"][i]), float(a[f"rpt_hs_{s}"][i])
                    delta = (ths - tpu) - (zhs - zpu)   # margin change: >0 timing helped
                    if abs(delta) <= 1e-12:
                        continue
                    key = -abs(delta)
                else:
                    zhs = ths = float("nan")
                    delta = tpu - zpu                   # <= 0 by construction
                    if delta >= -1e-12:
                        continue
                    key = delta
                cand.append((key, fp[i], int(a["entry"][i]), i, delta, zhs, ths, zpu, tpu))
            cand.sort()
            picked = cand[:args.n_per_region]
            last = abs(picked[-1][4]) if picked else float("nan")
            print(f"  {title:30s} R{region}: {len(cand):6d} events change, "
                  f"keeping {len(picked):2d} (smallest kept |delta| = {last:.3f})")
            for rank, (_, f, e, i, delta, zhs, ths, zpu, tpu) in enumerate(picked, start=1):
                tag = ("helped" if delta > 0 else "hurt") if region == 1 else "suppressed"
                rows.append(dict(
                    sample=args.sample, case=short, case_title=title, region=region,
                    rank=rank, tag=tag, delta=repr(delta), file_path=f, entry=e,
                    idx_hs=int(a["idx_hs"][i]), idx_pu=int(a["idx_pu"][i]),
                    hs_pt=f"{float(a['hs_pt'][i]):.2f}", hs_eta=f"{float(a['hs_eta'][i]):.3f}",
                    pu_pt=f"{float(a['pu_pt'][i]):.2f}", pu_eta=f"{float(a['pu_eta'][i]):.3f}",
                    rpt_hs_z=repr(zhs), rpt_hs_t=repr(ths), rpt_pu_z=repr(zpu), rpt_pu_t=repr(tpu),
                    t0=repr(float(a[f"t0_{s}"][i])), sig0=repr(float(a[f"sig0_{s}"][i])),
                    infl=repr(float(a[f"infl_{s}"][i])), gate_sigma=repr(float(a["gate_sigma"][i])),
                    ideal_times=int(ideal_t0 is not None), ideal_t0=ideal_t0 or ""))
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.file)),
                                   f"{args.sample}_region_display_candidates.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    events = {(r["file_path"], r["entry"]) for r in rows}
    print(f"wrote {out}: {len(rows)} displays over {len(events)} distinct events")


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def cmd_pick(args):
    rows = read_csv(args.candidates)
    pairs = sorted({(r["file_path"], int(r["entry"])) for r in rows})
    files = sorted({p[0] for p in pairs})
    stem = args.out or os.path.splitext(os.path.abspath(args.candidates))[0].replace(
        "_region_display_candidates", "_region_display_events")
    chain = ROOT.TChain("ntuple")
    for f in files:
        chain.Add(f)
    fout = ROOT.TFile(stem + ".root", "RECREATE")
    # A chain's clone follows it across file switches (TChain re-points every
    # clone's branch addresses on LoadTree), so one clone serves every file.
    clone = chain.CloneTree(0)
    chain.GetEntries()          # loads each file's entry offset
    offs = chain.GetTreeOffset()
    fidx = {f: i for i, f in enumerate(files)}
    mapping = []
    for new, (f, e) in enumerate(pairs):
        if chain.GetEntry(offs[fidx[f]] + e) <= 0:
            sys.exit(f"could not read {f} entry {e}")
        clone.Fill()
        mapping.append((f, e, new))
    clone.Write()
    fout.Close()
    with open(stem + "_map.csv", "w", newline="") as fm:
        w = csv.writer(fm)
        w.writerow(["file_path", "entry", "picked_entry"])
        w.writerows(mapping)
    print(f"wrote {stem}.root ({len(pairs)} events from {len(files)} files) and {stem}_map.csv")


def render_one(job):
    r, path, entry, out_dir, png = job
    region, legs = int(r["region"]), f"{r['idx_hs']},{r['idx_pu']}"
    expect = ",".join([r["rpt_hs_z"], r["rpt_hs_t"], r["rpt_pu_z"], r["rpt_pu_t"]])
    t0, sig0 = float(r["t0"]), float(r["sig0"])
    region_txt = (f"R1 #{r['rank']} (timing {r['tag']}, margin {float(r['delta']):+.2f})" if region == 1
                  else f"R2 #{r['rank']} (PU leg {float(r['rpt_pu_z']):.2f} -> {float(r['rpt_pu_t']):.2f})")
    stem = os.path.splitext(os.path.basename(r["file_path"]))[0]
    name = f"{int(r['rank']):02d}_{r['tag']}_{stem}_ev{r['entry']}"
    cmd = [sys.executable, "event_display.py",
           "--file_path", path, "--event_num", str(entry),
           "--extra_time", repr(t0), "--t0_sigma", repr(sig0),
           "--infl", r["infl"], "--gate_sigma", r["gate_sigma"], "--assoc_nsigma", "2.5",
           "--legs", legs, "--leg_labels", "HS,PU",
           "--case_label", f"{r['case_title']} = {t0:.1f} ps   {region_txt}",
           "--t0_label", r["case_title"],
           "--source_label", f"{os.path.basename(r['file_path'])} #{r['entry']}",
           "--expect_rpt", expect,
           "--output_dir", out_dir, "--output_name", name]
    if int(r["ideal_times"]):
        cmd += ["--ideal_times", "--ideal_t0", r["ideal_t0"]]
    p = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    log = os.path.join(out_dir, name + ".log")
    with open(log, "w") as f:
        f.write(" ".join(cmd) + "\n\n" + p.stdout + "\n" + p.stderr)
    checks = [ln for ln in p.stdout.splitlines()
              if ln.startswith(("RPT_CHECK", "REPLAY_CHECK", "LEG_LABEL_NOTE"))]
    pdf = os.path.join(out_dir, name + ".pdf")
    ok = (p.returncode == 0 and os.path.exists(pdf)
          and any(c.startswith("RPT_CHECK OK") for c in checks)
          and not any("MISMATCH" in c for c in checks)
          and (not int(r["ideal_times"]) or any(c.startswith("REPLAY_CHECK OK") for c in checks)))
    if ok and png:
        subprocess.run(["pdftoppm", "-f", "2", "-l", "2", "-r", "110", "-png", "-singlefile",
                        pdf, os.path.join(out_dir, name + "_p2")], check=False)
    return r, name, ok, p.returncode, checks


def cmd_render(args):
    rows = read_csv(args.candidates)
    if args.case:
        rows = [r for r in rows if r["case"] in args.case]
    if args.region:
        rows = [r for r in rows if int(r["region"]) == args.region]
    if args.max_rank:
        rows = [r for r in rows if int(r["rank"]) <= args.max_rank]
    remap = {}
    if args.pick_map:
        if not args.picked:
            sys.exit("--pick-map needs --picked")
        for m in read_csv(args.pick_map):
            remap[(m["file_path"], int(m["entry"]))] = int(m["picked_entry"])
    jobs = []
    for r in rows:
        key = (r["file_path"], int(r["entry"]))
        if remap:
            if key not in remap:
                sys.exit(f"{key} is not in {args.pick_map}")
            path, entry = os.path.abspath(args.picked), remap[key]
        else:
            path, entry = r["file_path"], int(r["entry"])
        out_dir = os.path.abspath(os.path.join(args.out_dir, r["sample"], r["case"], f"r{r['region']}"))
        os.makedirs(out_dir, exist_ok=True)
        jobs.append((r, path, entry, out_dir, not args.no_png))
    print(f"rendering {len(jobs)} displays with {args.jobs} workers ...", flush=True)
    bad, notes = [], 0
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for k, (r, name, ok, rc, checks) in enumerate(ex.map(render_one, jobs), start=1):
            notes += sum(c.startswith("LEG_LABEL_NOTE") for c in checks)
            if not ok:
                bad.append((r["case"], r["region"], name, rc, checks))
            if k % 20 == 0 or not ok:
                print(f"  [{k}/{len(jobs)}] {r['case']} R{r['region']} {name}: "
                      f"{'ok' if ok else 'FAILED rc=%d %s' % (rc, checks)}", flush=True)
    summary = os.path.join(os.path.abspath(args.out_dir), "render_summary.txt")
    with open(summary, "a") as f:
        f.write(f"{args.candidates}: {len(jobs)} rendered, {len(bad)} failed, "
                f"{notes} leg-label notes\n")
        for b in bad:
            f.write(f"  FAILED {b}\n")
    print(f"{len(jobs) - len(bad)}/{len(jobs)} displays OK "
          f"(RPT_CHECK{' + REPLAY_CHECK' if any(int(j[0]['ideal_times']) for j in jobs) else ''}); "
          f"{notes} leg-label notes; summary appended to {summary}")
    if bad:
        sys.exit(1)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = ap.add_subparsers(dest="cmd", required=True)
s = sub.add_parser("select", help="rank the region tree, write the candidates CSV")
s.add_argument("file", help="merged <prefix>rpt_v5_regions.root")
s.add_argument("--sample", required=True)
s.add_argument("--n-per-region", type=int, default=20)
s.add_argument("--require-fwd-acc-jet", action="store_true",
               help="only events with >= 1 jet > 30 GeV in 2.38 < |eta| < 4.0 (composition-plot set)")
s.add_argument("--out", help="CSV path (default beside the tree)")
p = sub.add_parser("pick", help="copy the candidates' events into one small file (run where the ntuples are)")
p.add_argument("candidates")
p.add_argument("--out", help="output stem (default <sample>_region_display_events beside the CSV)")
r = sub.add_parser("render", help="run event_display.py for every candidate, with cross-checks")
r.add_argument("candidates")
r.add_argument("--picked", help="picked-events ROOT file (from `pick`)")
r.add_argument("--pick-map", help="its mapping CSV")
r.add_argument("--out-dir", default="figs/rpt_regions/event_displays")
r.add_argument("--jobs", type=int, default=4)
r.add_argument("--case", nargs="*", help="only these cases (trkptz waves hgtd truth ideal_trk)")
r.add_argument("--region", type=int, choices=(1, 2))
r.add_argument("--max-rank", type=int, help="only the top N per case/region (quick checks)")
r.add_argument("--no-png", action="store_true", help="skip the page-2 PNG")
args = ap.parse_args()
ROOT.gROOT.SetBatch(True)
{"select": cmd_select, "pick": cmd_pick, "render": cmd_render}[args.cmd](args)
