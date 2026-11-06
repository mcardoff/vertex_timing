#!/usr/bin/env python3
"""Event displays for slides_tdr_timing, chosen from the deck's own trees by fixed rules (the user's request,
2026-11-05: "sprinkle in some event displays throughout").

Trees: util/tdr_fwdpair --pair=one --jvt=loose (section 3's: >= 1 VBS jet in HGTD, JVT + fJVT loose before pairing,
cleaned track times), with the kinematic BDT's held-out scores (tdr_fwdpair_kinbdt.py). Cases, one per point of the
talk; each rule is printed with the events it found:
  t0         section 1  Z+jets: the WAVeS t0 within 15 ps of the truth t_HS, Athena's HGTD t0 150-350 ps off
  pujet      section 2  VBF: a forward pileup VBS jet (30-50 GeV) that passes the ITk-only R_pT tagger at 85% and
                        that the WAVeS t0 gate (2 sigma, the TDR's) empties; the other VBS jet hard-scatter
  pujet_cf   section 2  VBF: a forward pileup VBS jet (30-50 GeV) the WAVeS gate flips from passing the ITk-only tagger
                        at 85% to failing it, the partner a CENTRAL hard-scatter jet (the central-forward pairs)
  pujet_z    section 2  the same, Z+jets, any partner
  jj_miss    section 3  Z+jets forward-forward, no BDT: two pileup jets that agree in time (jet-vs-jet pull < 2) while both
                        are >= 3 sigma from the WAVeS t0 -- the t0 veto fires, jet vs jet does not
  jj_catch   section 3  Z+jets forward-forward, no BDT: HS + PU, the jets >= 3 sigma apart, the HS jet on the t0
  veto_bkg   section 3  Z+jets after the BDT at 80%: a pileup VBS jet in HGTD that the 3-sigma WAVeS veto catches
  irreducible section 3 Z+jets after the BDT at 80%: two hard-scatter VBS jets, both timed, the veto keeps it
  sig_cost   section 3  VBF after the BDT at 80%: two hard-scatter VBS jets, the veto rejects it
Each case also requires a correct WAVeS t0 (within 20 ps of truth) except where stated, so the display shows the
case's own point and not a t0 failure.

  select  (laptop)  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_display_events.py select \\
                        --sig ../condor/vbf/vbf_onefwd_jvtLoose_tdr_fwdpair.root \\
                        --bkg ../condor/zjets/zjets_onefwd_jvtLoose_tdr_fwdpair.root \\
                        --bdt ../condor/zjets/vbf_vs_zjets_onefwd_jvtLoose_kinbdt.npz --out ../figs/tdr/displays/candidates.csv
  pick    (AF, PyROOT) python3 tdr_display_events.py pick --csv candidates.csv --out tdr_display_events
                        copies each candidate out of the RAW ntuples (the skims drop Track_hitInLastLayer), found by
                        the grid file number the trees store (file_idx) and the entry, checked by the VBS jets' pT/eta
  render  (laptop)  ... tdr_display_events.py render --csv ../figs/tdr/displays/candidates.csv \\
                        --picked ../condor/tdr_display_events --out-dir ../figs/tdr/displays
                        python/event_display.py with the row's WAVeS t0, the TDR's 2 sigma gate, ghost-matched R_pT
                        and cleaned times (--rpt_match ghost --clean_times), the core-time tracks of both VBS jets
                        marked, and --expect_rpt: every display must print RPT_CHECK OK (the tree's R_pT, exactly)
"""
import argparse
import csv
import glob
import itertools
import os
import subprocess
import sys

import numpy as np

RAW = {"vbf": "/data/mcardiff/exotic_superntuples/highstats_vbf", "zjets": "/data/mcardiff/exotic_superntuples/zjets"}
WP85_FWD = {"30_50": 0.1871491444419373, "50_inf": 0.16511387904857375}   # figs/tdr/grid/grid_tdr_pujet_wp85.json
GATE = 2.0                                                                # tdr::GATE_SIGMA_TDR
FIELDS = ["case", "rank", "sample", "file_idx", "entry", "mjj", "bdt", "t0_waves", "sig0_waves", "infl_waves",
          "ok_hgtd", "t0_hgtd", "t0_truth"] + [f"j{k}_{v}" for k in (1, 2) for v in
          ("pt", "eta", "phi", "hs", "pu", "core_t", "core_sig", "core_n", "rpt_itk", "rpt_t0_waves", "pull_waves")] + ["pull_jj"]


def cmd_select(a):
    import uproot
    import tdr_fwdpair_common as FC
    E = {"vbf": uproot.open(a.sig)["events"].arrays(library="np"), "zjets": uproot.open(a.bkg)["events"].arrays(library="np")}
    D = np.load(a.bdt)
    sc = {"vbf": D["sig_score"], "zjets": D["bkg_score"]}
    ws = E["vbf"]["weight"].astype(float)
    o = np.argsort(-sc["vbf"]); cs = np.cumsum(ws[o]) / ws.sum()
    thr80 = sc["vbf"][o][np.searchsorted(cs, 0.80)]
    cal = FC.calibration(E["vbf"], verbose=False)
    P = {s: FC.pulls(E[s], cal) for s in E}
    K = {s: FC.vetoes(E[s], cal)["waves"] for s in E}
    print(f"BDT cut at 80% signal efficiency: score > {thr80:.3f}")

    def common(s):
        e = E[s]
        return e, np.abs(e["t0_waves"] - e["t0_truth"]), np.abs(e["j1_eta"]), np.abs(e["j2_eta"])

    rows = []
    # ---- t0: Z+jets, WAVeS right, Athena's HGTD t0 wrong (but on the drawn time axis)
    e, dw, a1, a2 = common("zjets")
    dh = np.abs(e["t0_hgtd"] - e["t0_truth"])
    m = (e["ok_waves"] == 1) & (dw < 15) & (e["ok_hgtd"] == 1) & (dh > 150) & (dh < 350)
    rows += [("t0", "zjets", i) for i in np.flatnonzero(m)[np.argsort(dw[m])][: a.n]]
    print(f"t0: {m.sum()} Z+jets events (WAVeS within 15 ps, HGTD 150-350 ps off)")
    # ---- pujet: VBF, a forward PU VBS jet passing the ITk tagger, emptied by the WAVeS gate
    e, dw, a1, a2 = common("vbf")
    cand = []
    for k, other in ((1, 2), (2, 1)):
        ae, pt = np.abs(e[f"j{k}_eta"]), e[f"j{k}_pt"]
        thr = np.where(pt < 50, WP85_FWD["30_50"], WP85_FWD["50_inf"])
        m = ((e[f"j{k}_pu"] == 1) & (ae > 2.4) & (pt > 30) & (pt < 50) & (e[f"j{k}_rpt_itk"] >= thr)
             & (e[f"j{k}_rpt_t0_waves"] <= 1e-9) & (e[f"j{other}_hs"] == 1) & (e["ok_waves"] == 1) & (dw < 20))
        cand += [(i, e[f"j{k}_rpt_itk"][i]) for i in np.flatnonzero(m)]
    cand.sort(key=lambda x: -x[1])
    rows += [("pujet", "vbf", i) for i, _ in cand[: a.n]]
    print(f"pujet: {len(cand)} VBF events (forward PU VBS jet 30-50 GeV, ITk R_pT >= the 85% cut, gated to 0)")
    # ---- pujet_cf / pujet_z (the user's request, 2026-11-05: "a couple more exemplary event displays from VBF and
    # Z+jets" after the first): the same jet, in a VBF central-forward pair (79% of the events; the first is
    # forward-forward) and in Z+jets with any partner. Exemplary rather than extreme: the gate FLIPS the tagger (ITk
    # R_pT passes the 85% cut, the WAVeS-gated one fails it; no Z+jets candidate is gated all the way to 0), ranked
    # by the R_pT drop
    for case, s_, cf in (("pujet_cf", "vbf", True), ("pujet_z", "zjets", False)):
        e, dw, a1, a2 = common(s_)
        cand = []
        for k, other in ((1, 2), (2, 1)):
            ae, pt = np.abs(e[f"j{k}_eta"]), e[f"j{k}_pt"]
            thr = np.where(pt < 50, WP85_FWD["30_50"], WP85_FWD["50_inf"])
            m = ((e[f"j{k}_pu"] == 1) & (ae > 2.4) & (pt > 30) & (pt < 50) & (e[f"j{k}_rpt_itk"] >= thr)
                 & (e[f"j{k}_rpt_t0_waves"] < thr) & (e["ok_waves"] == 1) & (dw < 20))
            if cf: m &= (e[f"j{other}_hs"] == 1) & (np.abs(e[f"j{other}_eta"]) < 2.4)
            cand += [(i, e[f"j{k}_rpt_itk"][i] - e[f"j{k}_rpt_t0_waves"][i]) for i in np.flatnonzero(m)]
        cand.sort(key=lambda x: -x[1])
        rows += [(case, s_, i) for i, _ in cand[: a.n]]
        print(f"{case}: {len(cand)} {s_} events (forward PU VBS jet 30-50 GeV: ITk R_pT passes the 85% cut, gated fails it"
              + (", partner a central HS jet)" if cf else ")"))
    # ---- veto_bkg: Z+jets after the BDT, PU jet in HGTD caught by the WAVeS veto
    e, dw, a1, a2 = common("zjets")
    bdt = sc["zjets"] >= thr80
    cand = []
    for k in (1, 2):
        pk = P["zjets"][f"waves_{k}"]
        m = (bdt & (e[f"j{k}_pu"] == 1) & (np.abs(e[f"j{k}_eta"]) > 2.4) & ~K["zjets"] & (np.nan_to_num(pk) >= 3)
             & (e[f"j{k}_core_n"] >= 2) & (dw < 20))
        cand += [(i, e[f"j{k}_core_n"][i] + 0.01 * min(pk[i], 20)) for i in np.flatnonzero(m)]
    cand.sort(key=lambda x: -x[1])
    rows += [("veto_bkg", "zjets", i) for i, _ in cand[: a.n]]
    print(f"veto_bkg: {len(cand)} Z+jets events after the BDT at 80% (PU jet in HGTD, >= 2 core tracks, caught)")
    # ---- irreducible: Z+jets after the BDT, HS + HS, both timed, kept
    m = bdt & (e["j1_hs"] == 1) & (e["j2_hs"] == 1) & K["zjets"] & (e["j1_core_n"] >= 2) & (e["j2_core_n"] >= 2) & (dw < 20)
    rows += [("irreducible", "zjets", i) for i in np.flatnonzero(m)[np.argsort(-e["mjj"][m])][: a.n]]
    print(f"irreducible: {m.sum()} Z+jets events after the BDT at 80% (HS + HS, both timed, kept)")
    # ---- sig_cost: VBF after the BDT, HS + HS, rejected by the veto with a correct t0
    e, dw, a1, a2 = common("vbf")
    bdt = sc["vbf"] >= thr80
    m = bdt & (e["j1_hs"] == 1) & (e["j2_hs"] == 1) & ~K["vbf"] & (dw < 20)
    pmax = np.fmax(np.nan_to_num(P["vbf"]["waves_1"]), np.nan_to_num(P["vbf"]["waves_2"]))
    order = np.flatnonzero(m)[np.argsort((pmax[m] >= 12) * 1e4 + np.abs(e["mjj"][m] - 1500.0))]   # typical m_jj; off-scale pulls last
    rows += [("sig_cost", "vbf", i) for i in order[: a.n]]
    print(f"sig_cost: {m.sum()} VBF events after the BDT at 80% (HS + HS, correct t0, vetoed)")

    # ---- jj_miss / jj_catch (the user, 2026-11-06: "a couple examples of the jet-v-jet cases"): Z+jets forward-forward
    # pairs, both jets timed by >= 2 core tracks, a correct WAVeS t0, no BDT (as the pull figure)
    e, dw, a1, a2 = common("zjets")
    Pz = P["zjets"]; jj = np.nan_to_num(Pz["jj"], nan=-1.0)
    p1, p2 = np.nan_to_num(Pz["waves_1"], nan=0.0), np.nan_to_num(Pz["waves_2"], nan=0.0)
    base = (e["cat"] == 2) & (e["j1_core_n"] >= 2) & (e["j2_core_n"] >= 2) & (e["ok_waves"] == 1) & (dw < 20)
    # jet vs jet misses: two pileup jets that agree with each other (pull < 2), both >= 3 sigma from t0 (the t0 veto fires)
    m = base & (e["j1_pu"] == 1) & (e["j2_pu"] == 1) & (jj >= 0) & (jj < 2) & (p1 >= 3) & (p2 >= 3)
    sc_ = np.fmin(p1, p2) - jj                           # together, and far from t0; off-scale times last
    sc_ = np.where(np.fmax(p1, p2) > 25, sc_ - 100, sc_)
    rows += [("jj_miss", "zjets", i) for i in np.flatnonzero(m)[np.argsort(-sc_[m])][: a.n]]
    print(f"jj_miss: {m.sum()} Z+jets FF events (PU + PU, jet-vs-jet pull < 2, both >= 3 sigma from the WAVeS t0)")
    # jet vs jet catches: HS + PU, the jets >= 3 sigma apart, the HS jet on t0 (< 2) and the PU jet off it (>= 3)
    hs1 = (e["j1_hs"] == 1) & (e["j2_pu"] == 1); hs2 = (e["j2_hs"] == 1) & (e["j1_pu"] == 1)
    m = base & (jj >= 3) & ((hs1 & (p1 < 2) & (p2 >= 3)) | (hs2 & (p2 < 2) & (p1 >= 3)))
    sc_ = e["j1_core_n"] + e["j2_core_n"] - 100 * (jj > 20)
    rows += [("jj_catch", "zjets", i) for i in np.flatnonzero(m)[np.argsort(-sc_[m])][: a.n]]
    print(f"jj_catch: {m.sum()} Z+jets FF events (HS + PU, jet-vs-jet pull >= 3, HS jet on t0, PU jet off it)")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    rank = {}
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader()
        for case, s, i in rows:
            rank[case] = rank.get(case, 0) + 1
            e = E[s]
            r = dict(case=case, rank=rank[case], sample=s, file_idx=int(e["file_idx"][i]), entry=int(e["entry"][i]),
                     mjj=float(e["mjj"][i]), bdt=float(sc[s][i]))
            for v in ("t0_waves", "sig0_waves", "infl_waves", "ok_hgtd", "t0_hgtd", "t0_truth"): r[v] = repr(float(e[v][i]))
            for k in (1, 2):
                for v in ("pt", "eta", "phi", "hs", "pu", "core_t", "core_sig", "core_n", "rpt_itk", "rpt_t0_waves"):
                    r[f"j{k}_{v}"] = repr(float(e[f"j{k}_{v}"][i]))
                r[f"j{k}_pull_waves"] = repr(float(P[s][f"waves_{k}"][i]))
            r["pull_jj"] = repr(float(P[s]["jj"][i]))
            w.writerow(r)
    print(f"wrote {a.out}: {len(rows)} candidates")


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def cmd_pick(a):
    import ROOT
    rows = read_csv(a.csv)
    out = []
    for s in sorted({r["sample"] for r in rows}):
        want = sorted({(int(r["file_idx"]), int(r["entry"]), float(r["j1_pt"]), float(r["j1_eta"]))
                       for r in rows if r["sample"] == s})
        chain = ROOT.TChain("ntuple")
        found = []
        for fi, ent, pt, eta in want:
            hits = []
            for path in glob.glob(f"{RAW[s]}/**/*._{fi:06d}.SuperNtuple.root", recursive=True):
                tf = ROOT.TFile.Open(path); t = tf.Get("ntuple"); t.GetEntry(ent)
                if any(abs(p - pt) < 1e-3 and abs(h - eta) < 1e-4 for p, h in zip(t.AntiKt4EMTopoJets_pt, t.AntiKt4EMTopoJets_eta)):
                    hits.append(path)
                tf.Close()
            if len(hits) != 1:
                sys.exit(f"{s} file {fi} entry {ent}: {len(hits)} raw files hold the event ({hits})")
            found.append((hits[0], ent, fi))
        files = sorted({h for h, _, _ in found})
        offs, n0 = [], 0
        for f in files:                                  # each file's chain offset, counted here: TChain::GetTreeOffset
            chain.Add(f)                                 # can come back unfilled (kMaxEntries) for later files
            tf = ROOT.TFile.Open(f); offs.append(n0); n0 += tf.Get("ntuple").GetEntries(); tf.Close()
        fidx = {f: i for i, f in enumerate(files)}
        fout = ROOT.TFile(f"{a.out}_{s}.root", "RECREATE")
        clone = chain.CloneTree(0)
        for new, (path, ent, fi) in enumerate(found):
            if chain.GetEntry(offs[fidx[path]] + ent) <= 0: sys.exit(f"could not read {path} entry {ent}")
            clone.Fill()
            out.append((s, fi, ent, path, new))
        clone.Write(); fout.Close()
        print(f"wrote {a.out}_{s}.root ({len(found)} events)")
    with open(f"{a.out}_map.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sample", "file_idx", "entry", "file_path", "picked_entry"]); w.writerows(out)


def jet_index(g, pt, eta):
    m = [j for j, (p, h) in enumerate(zip(g["AntiKt4EMTopoJets_pt"], g["AntiKt4EMTopoJets_eta"]))
         if abs(p - pt) < 1e-3 and abs(h - eta) < 1e-4]
    if len(m) != 1: sys.exit(f"jet (pt {pt}, eta {eta}): {len(m)} matches")
    return m[0]


def core_tracks(g, j, core_t, core_n):
    """the core-time tracks of jet j: util/tdr_fwdpair.cxx's candidates (cleaned timed ghost tracks, 1 < pT < 30 GeV,
    quality, 2.38 < |eta| < 4.0, dR < 0.2), the subset of size core_n whose inverse-variance mean IS core_t"""
    if core_n < 1: return []
    je, jp = g["AntiKt4EMTopoJets_eta"][j], g["AntiKt4EMTopoJets_phi"][j]
    cand = []
    for k in g["AntiKt4EMTopoJets_ghostTrack_idx"][j]:
        ae = abs(g["Track_eta"][k])
        if g["Track_hasValidTime"][k] != 1 or not g["Track_hitInLastLayer"][k] or (3.5 < ae < 3.9 and g["Track_nHGTDHits"][k] < 2):
            continue
        if not g["Track_timeRes"][k] > 0 or not g["Track_quality"][k] or not (1.0 <= g["Track_pt"][k] <= 30.0) or not (2.38 <= ae <= 4.0):
            continue
        if np.hypot(g["Track_eta"][k] - je, (g["Track_phi"][k] - jp + np.pi) % (2 * np.pi) - np.pi) >= 0.2:
            continue
        cand.append(int(k))
    best, bd = None, np.inf
    for sub in itertools.combinations(cand, core_n):
        w = np.array([1 / float(g["Track_timeRes"][k]) ** 2 for k in sub])
        t = np.sum(w * np.array([float(g["Track_time"][k]) for k in sub])) / w.sum()
        if abs(t - core_t) < bd: best, bd = sub, abs(t - core_t)
    if bd > 1e-3: sys.exit(f"jet {j}: no {core_n}-track subset of {cand} has the core time {core_t} (best off by {bd})")
    return list(best)


def cmd_render(a):
    import uproot
    rows, mp = read_csv(a.csv), read_csv(a.picked + "_map.csv")
    where = {(m["sample"], int(m["file_idx"]), int(m["entry"])): (m["file_path"], int(m["picked_entry"])) for m in mp}
    os.makedirs(a.out_dir, exist_ok=True)
    sel = set(a.only.split(",")) if a.only else None
    BR = ["AntiKt4EMTopoJets_pt", "AntiKt4EMTopoJets_eta", "AntiKt4EMTopoJets_phi", "AntiKt4EMTopoJets_ghostTrack_idx",
          "Track_pt", "Track_eta", "Track_phi", "Track_quality", "Track_time", "Track_timeRes", "Track_hasValidTime",
          "Track_hitInLastLayer", "Track_nHGTDHits"]
    TITLE = {"t0": "the vertex time", "pujet": "a forward pileup jet", "veto_bkg": "Z+jets caught by the veto",
             "irreducible": "Z+jets the veto keeps", "sig_cost": "VBF lost to the veto"}
    summary = []
    for r in rows:
        name = f"{r['case']}_{int(r['rank'])}"
        if sel and name not in sel and r["case"] not in sel: continue
        orig, pe = where[(r["sample"], int(r["file_idx"]), int(r["entry"]))]
        pf = os.path.abspath(f"{a.picked}_{r['sample']}.root")
        g = uproot.open(pf)["ntuple"].arrays(BR, entry_start=pe, entry_stop=pe + 1)
        g = {k: g[k][0].to_list() for k in BR}
        J = [jet_index(g, float(r[f"j{k}_pt"]), float(r[f"j{k}_eta"])) for k in (1, 2)]
        lab = ["HS" if float(r[f"j{k}_hs"]) == 1 else "PU" if float(r[f"j{k}_pu"]) == 1 else "other" for k in (1, 2)]
        jt = [f"{J[k - 1]}:" + ",".join(str(i) for i in core_tracks(g, J[k - 1], float(r[f"j{k}_core_t"]), int(float(r[f"j{k}_core_n"]))))
              for k in (1, 2)]
        t0, tt = float(r["t0_waves"]), float(r["t0_truth"])
        pulls = ", ".join(f"jet {J[k - 1] + 1} pull {float(r[f'j{k}_pull_waves']):.1f}" for k in (1, 2)
                          if np.isfinite(float(r[f"j{k}_pull_waves"])) and r["case"] not in ("t0", "pujet", "pujet_cf", "pujet_z"))
        cmd = [sys.executable, "event_display.py", "--file_path", pf, "--event_num", str(pe),
               "--extra_time", repr(t0), "--t0_sigma", r["sig0_waves"], "--infl", r["infl_waves"],
               "--gate_sigma", repr(GATE), "--assoc_nsigma", "2.5", "--rpt_match", "ghost", "--clean_times",
               "--legs", f"{J[0]},{J[1]}", "--leg_labels", ",".join(lab),
               "--expect_rpt", ",".join(r[f"j{k}_{v}"] for k in (1, 2) for v in ("rpt_itk", "rpt_t0_waves")),
               "--jet_time_tracks", ";".join(jt),
               "--case_label", f"WAVeS t0 = {t0:.1f} ps   (truth t_HS = {tt:.1f} ps)" + (f"   |   {pulls}" if pulls else "")
                               + (f"   |   jet vs jet pull {float(r['pull_jj']):.1f}" if r["case"].startswith("jj_") else ""),
               "--t0_label", "WAVeS t0",
               "--source_label", f"{'VBF H->inv.' if r['sample'] == 'vbf' else 'Z(->ll)+jets'}   "
                                 f"{os.path.basename(orig)} #{r['entry']}   |   JVT + fJVT loose",
               "--output_dir", os.path.abspath(a.out_dir), "--output_name", name]
        res = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(a.out_dir, name + ".log"), "w") as f:
            f.write(" ".join(cmd) + "\n\n" + res.stdout + res.stderr)
        chk = [l for l in res.stdout.splitlines() if l.startswith("RPT_CHECK")]
        summary.append(f"{name:14s} {r['sample']:5s} {os.path.basename(orig)} #{r['entry']}: rc {res.returncode}, "
                       f"{chk[0] if chk else 'no RPT_CHECK'}; legs {lab}, core tracks {jt}")
        print(summary[-1], flush=True)
    with open(os.path.join(a.out_dir, "render_summary.txt"), "a") as f:
        f.write("\n".join(summary) + "\n")


def cmd_times(a):
    """the time panel beside a display, as slides_time_veto_money's (python/vbs_veto_event_cards.py draw_times): one row
    per t0 (truth t_HS, WAVeS, Athena's HGTD) and one per VBS jet, each with the +-2 / +-3 sigma the test uses
    (t0: quoted x inflation; jet: core sigma x f(n), the veto's calibration), each t0's verdict with the two jets' pulls,
    and jet vs jet. Read from the trees, the event found by (file_idx, entry) and checked by the leading jet's pT;
    the jets named as event_display.py numbers them (the render log's --legs)."""
    import re
    import uproot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 13})
    import tdr_fwdpair_common as FC
    E = {"vbf": uproot.open(a.sig)["events"].arrays(library="np"), "zjets": uproot.open(a.bkg)["events"].arrays(library="np")}
    cal = FC.calibration(E["vbf"], verbose=False)
    P = {s_: FC.pulls(e, cal) for s_, e in E.items()}
    K = {s_: FC.vetoes(e, cal, a.nsigma) for s_, e in E.items()}
    rows = {f"{r['case']}_{r['rank']}": r for r in read_csv(a.csv)}
    SRC = [("waves", "WAVeS $t_0$", "#E69F00"), ("hgtd", "HGTD $t_0$ (Athena)", "#3F90DA")]
    for name in a.only.split(","):
        r = rows[name]; s_ = r["sample"]; e = E[s_]
        i = np.flatnonzero((e["file_idx"] == int(r["file_idx"])) & (e["entry"] == int(r["entry"])) & (np.abs(e["j1_pt"] - float(r["j1_pt"])) < 1e-3))
        assert len(i) == 1, (name, i)
        i = int(i[0])
        legs = re.search(r"--legs (\d+),(\d+)", open(os.path.join(a.logs, name + ".log")).read())
        names = [f"Jet {int(legs.group(k)) + 1}" for k in (1, 2)]
        tt = e["t0_truth"][i]
        t0rows = [("truth $t_{\\rm HS}$", tt, 0.0, "#5B3A29", "truth")]
        for src, lab, col in SRC:
            ok = e[f"ok_{src}"][i] == 1
            t0rows.append((lab + ("" if ok else "\n(no time)"), e[f"t0_{src}"][i] if ok else np.nan,
                           e[f"infl_{src}"][i] * e[f"sig0_{src}"][i] if ok else np.nan, col, src))
        jets = []
        for k, nm in zip((1, 2), names):
            role = "HS" if e[f"j{k}_hs"][i] == 1 else "PU" if e[f"j{k}_pu"][i] == 1 else "neither label"
            n = int(e[f"j{k}_core_n"][i])
            txt = f"{nm}: {role}, $\\eta$ = {e[f'j{k}_eta'][i]:+.2f}, {e[f'j{k}_pt'][i]:.0f} GeV"
            if n >= 1:
                jets.append((txt + f"\njet time: {n} track{'s' if n > 1 else ''}", e[f"j{k}_core_t"][i],
                             cal(np.array([n]))[0] * e[f"j{k}_core_sig"][i], "#2E9E4F" if role == "HS" else "#555555", f"j{k}"))
            else:
                jets.append((txt + ("\nno jet time (central)" if abs(e[f"j{k}_eta"][i]) < 2.4 else "\nno jet time"),
                             np.nan, np.nan, "#999999", f"j{k}"))
        allrows = t0rows + jets
        ys = np.array([1.6 * len(jets) - 0.3 + (len(t0rows) - 1 - k) for k in range(len(t0rows))]
                      + [1.6 * (len(jets) - 1 - k) for k in range(len(jets))])
        fin = [(v, x) for _, v, x, _, _ in allrows if np.isfinite(v)]
        lo = min([v - 3 * x for v, x in fin if np.isfinite(x)] + [min(v for v, _ in fin) - 40])
        hi = max([v + 3 * x for v, x in fin if np.isfinite(x)] + [max(v for v, _ in fin) + 40])
        pad = 0.08 * (hi - lo)
        fig = plt.figure(figsize=(8.6, 4.6))
        at = fig.add_axes([0.35, 0.16, 0.40, 0.70])
        for y, (lab, v, sg, col, tag) in zip(ys, allrows):
            if np.isfinite(v):
                if sg > 0:   # thick: +-2 sigma; thin, capped: out to +-3 sigma (the veto's threshold)
                    at.plot([v - 3 * sg, v + 3 * sg], [y, y], color=col, lw=1.2, alpha=0.8, solid_capstyle="butt")
                    for xx in (v - 3 * sg, v + 3 * sg): at.plot([xx] * 2, [y - 0.13, y + 0.13], color=col, lw=1.2, alpha=0.8)
                    at.plot([v - 2 * sg, v + 2 * sg], [y, y], color=col, lw=5, solid_capstyle="butt")
                at.plot([v], [y], "s" if tag in ("waves", "hgtd") else "o", color=col, ms=8, markeredgecolor="white",
                        markeredgewidth=0.8, zorder=5)
            if tag in ("truth", "waves", "hgtd"):
                if e[f"ok_{tag}"][i] != 1:
                    verdict, vcol = "no $t_0$: kept", "#777777"
                else:
                    pp = [P[s_][f"{tag}_{k}"][i] for k in (1, 2)]
                    ps = ", ".join("\u2013" if not np.isfinite(x) else f"{x:.1f}$\\sigma$" for x in pp)
                    vet = not K[s_][tag][i]
                    verdict, vcol = ("VETOED" if vet else "kept") + f"  ({ps})", ("#B03A2E" if vet else "#2E7D32")
                at.text(1.015, y, verdict, transform=at.get_yaxis_transform(), ha="left", va="center", fontsize=12,
                        color=vcol, fontweight="bold")
        at.set_yticks(ys); at.set_yticklabels([x[0] for x in allrows], fontsize=12)
        for tl, x in zip(at.get_yticklabels(), allrows): tl.set_color(x[3])
        at.axhline(1.6 * (len(jets) - 1) + 0.85, color="#BBBBBB", lw=0.8)
        pj = P[s_]["jj"][i]
        if np.isfinite(pj):
            vet = not K[s_]["jj"][i]
            jtxt, jcol = ("VETOED" if vet else "kept") + f"  ({pj:.1f}$\\sigma$)", ("#B03A2E" if vet else "#2E7D32")
        else:
            jtxt, jcol = "n/a: one jet timed", "#777777"
        at.text(1.015, 0.8, "jet vs jet:\n" + jtxt, transform=at.get_yaxis_transform(), ha="left", va="center",
                fontsize=12, color=jcol, fontweight="bold")
        at.axvline(tt, color="#5B3A29", lw=0.8, ls=":")
        at.set_xlim(lo - pad, hi + pad); at.set_ylim(-0.9, ys.max() + 0.6)
        at.set_xlabel("time [ps]   (thick: $\\pm 2\\sigma$, thin: $\\pm 3\\sigma$; $\\sigma$ as used by the veto)")
        at.text(1.015, ys.max() + 0.75, f"jet vs $t_0$ at {a.nsigma:g}$\\sigma$\n({names[0]}, {names[1]} pulls)",
                transform=at.get_yaxis_transform(), ha="left", va="bottom", fontsize=11, color="#555555")
        out = os.path.join(a.out_dir, f"times_{name}")
        fig.savefig(out + ".png", dpi=170, bbox_inches="tight", pad_inches=0.08)
        fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
        print(f"wrote {out}.png / .pdf")


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = ap.add_subparsers(dest="cmd", required=True)
p = sub.add_parser("select"); p.add_argument("--sig", required=True); p.add_argument("--bkg", required=True)
p.add_argument("--bdt", required=True); p.add_argument("--out", required=True); p.add_argument("--n", type=int, default=4)
p = sub.add_parser("pick"); p.add_argument("--csv", required=True); p.add_argument("--out", required=True)
p = sub.add_parser("render"); p.add_argument("--csv", required=True); p.add_argument("--picked", required=True)
p.add_argument("--out-dir", required=True); p.add_argument("--only", help="comma list of names (case_rank) or cases")
p = sub.add_parser("times"); p.add_argument("--csv", required=True); p.add_argument("--sig", required=True)
p.add_argument("--bkg", required=True); p.add_argument("--logs", required=True, help="the render logs (for the jet numbering)")
p.add_argument("--out-dir", required=True); p.add_argument("--only", required=True, help="comma list of names (case_rank)")
p.add_argument("--nsigma", type=float, default=3.0)
a = ap.parse_args()
{"select": cmd_select, "pick": cmd_pick, "render": cmd_render, "times": cmd_times}[a.cmd](a)
