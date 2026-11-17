---
name: condor-grid-submission
description: Submit *_hist / skim jobs to the UChicago AF HTCondor pool — sharding, transfer_output_remaps, hist_merge vs hadd, selection tags, and measured job-cost data.
---

## Grid Submission (condor)

Templates in `condor/` run `clustering_hist` / `rpt_v5_hist` on the UChicago AF HTCondor pool, one job per `--sample` value (`vbf`, `zjets`, `dijet`). Condor only runs the **histogramming** stage — the plotting stage is fast/local and runs afterward on whatever machine the `.root` file lands on, with no condor submission needed:

- **`condor/run_analysis.sh`** — generic executable, invoked as `run_analysis.sh <executable> <sample> [<threads>]`. Sources `atlasLocalSetup.sh` + `lsetup "root <version>"` (relaxing `set -e`/`-u` around just that call — those scripts reference unset internal variables like `ALRB_frontlineSite` and aren't strict-mode safe), then runs the transferred `<executable>` **directly from the job's scratch directory** with `--sample=<sample>` and, if given, `--threads=<N>`. No longer moves the binary into a fake `build/` subdir first (see `src/sample_config.h`'s output-path note below) — that indirection existed only to make the old `../<sample>` output convention resolve to the scratch root; the grid-sample output paths are cwd-relative now, so it's unnecessary and was a latent risk (if that step were ever skipped, `../<sample>` would climb out of the scratch directory into its parent, which may not be writable). Unchanged by the histogramming/plotting split otherwise — generic enough for both the legacy and `*_hist` executables.
- **`condor/clustering_hist.sub`** / **`condor/rpt_v5_hist.sub`** — one `.sub` file per histogramming executable, each queuing 3 jobs via `queue sample in (vbf, zjets, dijet)`. Submit each separately (own `condor_submit` call) to get independent job clusters. `arguments` passes `4` as the `<threads>` arg, matching `request_cpus = 4` — both `*_hist` executables parallelize their event loop via `TTreeProcessorMT` (see Multicore in Development Notes), and condor's cgroup throttles any threads spawned beyond `request_cpus`, so the two must track each other.
- **`condor/skim_ntuples.sub`** — regenerates the skimmed+slimmed sample copies (see "Skimmed + slimmed copies" under Data Format). One sample per submit, `$(nshards)` jobs; writes straight to `/data/mcardiff/skimmed_ntuples/<sample>/` and sets `transfer_output_files = ""`, since the output goes to a ceph path the execute hosts already mount rather than back to the submit node. Belongs on condor rather than the login node: the four samples are ~2.1 TB read twice, and 24 concurrent readers on a login node pushed cold reads from 0.3 s/file to 8 s/file and the node to load 94.
- **`condor/clustering_dt.sub`** / **`condor/rpt_v5.sub`** — legacy templates for the monolithic executables; still work, kept for reference/fallback, not the primary workflow.

Execute hosts do **not** share a filesystem with the submit host on this pool — confirmed via a `TARGET.FileSystemDomain == MY.FileSystemDomain` match failure when `should_transfer_files = NO` was tried first (0 of 157 slots matched). The `*_hist.sub` files use `transfer_input_files` (stages the prebuilt binary in) and `transfer_output_files = $(sample)/$(sample)_<name>.root` (stages back a single small histogram file — much smaller than the full PDF tree the legacy `.sub` files' `transfer_output_files = $(sample)` transfers). Ntuple inputs are read via the absolute `/data/mcardiff/exotic_superntuples/...` path baked into `--sample=` (see `src/sample_config.h`), so they don't need transferring.

**`transfer_output_remaps` is used** on both `*_hist.sub` files, though its role has narrowed since output file names became sample-prefixed (see "Sample-prefixed output paths" above): condor's default output-transfer behavior flattens any nested `transfer_output_files` path to its basename when staging back, it does not recreate `$(sample)/` in the submit directory on its own. Historically (when every sample wrote a bare `clustering_hist.root`/`rpt_v5_hist.root`), skipping the remap meant all 3 queued jobs silently overwrote each other's output — confirmed by running the real `rpt_v5_hist` jobs once without one: only the last-to-finish sample's histogram survived, and no `vbf/`, `zjets/`, `dijet/` subdirectories were created at all. Now that each job's basename is already unique (`vbf_clustering_hist.root`, `zjets_clustering_hist.root`, `dijet_clustering_hist.root`, ...), that particular collision can no longer happen — the remap is kept purely so the file lands nested under `$(sample)/` on the submit side, matching the local `--sample=` convention (so `clustering_plot --sample=<sample>` finds it without needing `--hist-file=`). (The legacy `.sub` files' `transfer_output_files = $(sample)` isn't affected either way — that argument has no nested subdirectory, so its basename already equals the whole thing they want back.) The obvious alternative, `preserve_relative_paths = True`, was tried first and rejected: it's a submit-wide knob that also governs `transfer_input_files`, and turned `../build/<exe>` (which deliberately climbs *out* of the sandbox to reach the built binary) from silently-flattened into a hard `Attempt to write to illegal sandbox path` hold — confirmed by running the real jobs with it set. `transfer_output_remaps` only touches the output side, so it can't have that side effect.

**`transfer_output_remaps`'s LHS (`name`) must be a bare filename** — per `condor_submit(1)`: "You cannot specify directories to be remapped." A first attempt using a full nested path (`$(sample)/hists/clustering_hist.root`, back when the histogram still lived under a `hists/` subdirectory) on the LHS silently failed to match anything and fell back to the same basename-flattening bug it was meant to fix — confirmed against the real jobs. Directories *are* allowed on the RHS (`newname`) and are created as needed, so the working form is `transfer_output_remaps = "$(sample)_clustering_hist.root=$(sample)/$(sample)_clustering_hist.root"` (bare sample-prefixed basename on the left, full destination path on the right).

```bash
cd build && cmake .. && make clustering_hist rpt_v5_hist   # build before submitting; condor doesn't rebuild
cd condor
condor_submit -dump dryrun.ad clustering_hist.sub          # dry-run: validate without queuing
condor_submit clustering_hist.sub
condor_submit rpt_v5_hist.sub                               # separate job cluster
condor_q
# After jobs complete, the .root files land in
# ../<sample>/<sample>_{clustering,rpt_v5}_hist.root. Plot locally — fast, no
# further condor submission:
cd ../build
./clustering_plot --sample=<sample>
./rpt_v5_plot --sample=<sample>
# A tagged run (selargs below) is plotted with the SAME flags, e.g.
# ./rpt_v5_plot --sample=<sample> --jvt=loose: the plot stage builds the tag
# from them to find the file, then checks the file's record against them.
```

Logs land in `condor/logs/clustering_hist_<sample>.<ClusterId>.<ProcId>.{out,err,log}` (prefixed `rpt_v5_hist_` for the rpt_v5_hist cluster). `stream_output`/`stream_error` are **on**, so `.out`/`.err` are written back to the submit node as the job runs rather than only on exit — without that, `condor_tail` mid-job shows only whatever happened to be flushed.

### Job cost: where the time actually goes

Measured from the three-sample `rpt_v5_hist` clusters of 2026-07-03 (`condor/logs/`, three repeats each). **The event counts here are stale — see the 2026-08-26 table below for current values.** Kept because the *wall/CPU* analysis that motivated the sharding work is still the record of how the problem was found:

| sample | files | events | wall | CPU | wall/event |
|---|---:|---:|---:|---:|---:|
| vbf   | 1002 |   499,100 |  24.9 min | 12.2 min | 2.99 ms |
| dijet | 1414 |   699,500 |  25.2 min |  8.4 min | 2.16 ms |
| zjets | 2005 | 8,084,000 | 183.3 min | 108.0 min | 1.36 ms |

**Z+jets is not inefficient — it is 16.2x larger than vbf** (both DSIDs, Zee+Zmumu). Per event it is the *fastest* of the three. File count is not a driver either: dijet has 1.4x vbf's files and 1.4x its events at the same wall time.

**Startup dead time is the real inefficiency, and it is not proportional to anything.** From the condor `.log` image-size timeline (memory sits at ~60 MB with ROOT loaded and nothing allocated, then jumps to ~500 MB when the worker threads book their histograms):

| sample | run 1 | run 2 | run 3 |
|---|---|---|---|
| vbf   | 20.2 min |  5.1 min | 15.2 min |
| dijet | 10.2 min | 20.2 min | 10.1 min |
| zjets | 30.2 min | 25.2 min | 20.2 min |

That is 5–30 min, varying **4x run-to-run on the same sample** — filesystem contention, not a per-sample property. It is 15–20% of a zjets job but **40–80% of a vbf or dijet job**. Sampled at condor's 5-minute granularity, so ±5 min.

Mechanism: `TChain::Add` is lazy (opens nothing), but `TChain::GetEntries()` opens **every file** in the chain to sum entry counts — 0.3–1.2 s per file on `/data` against 4.3 ms on local SSD. That call has been removed from `setupChain` (the "did we find files?" question is answered by `GetListOfFiles()->GetEntries()` with no I/O) and from the `*_hist` progress denominators, leaving only `TTreeProcessorMT`'s own pass. **Do not reintroduce it** — see the note on `setupChain` in `src/event_processing.h`.

How much that actually saved, measured on the local 33-file sample (time from start to the first event reaching a worker):

| | serial startup |
|---|---:|
| with `GetEntries()` | 0.205 s (6.1 ms/file) |
| deferred | 0.103 s (3.1 ms/file) |

**~50%**, not more: `TTreeProcessorMT` still opens every file serially for cluster boundaries before dispatching any task. In the old ordering that second pass was nearly free (0.004 s) only because `GetEntries` had already warmed the files.

**What is left, and what does NOT fix it.** The remaining pass is latency-bound, so overlapping it looks attractive — opening every file concurrently first does collapse the serial scan to 0.004 s. But the prewarm then costs as much as it saves (0.132 s at 8 threads vs a 0.117 s scan): on local SSD the opens are bandwidth-bound and concurrency buys nothing. Whether the AF's shared storage behaves differently is untested and should not be assumed — an earlier plausible-sounding fix here (hoisting `passBasicCuts` above the Z+jets lepton work, so rejected events never touch per-track branches) measured at **1 µs/event**, i.e. 0.1 min over 8M events, and was dropped.

### Sharding a sample across jobs

`--file-shard=<i>/<N>` restricts a `*_hist` process to files i, i+N, i+2N, … of the **sorted** file list, and `hist_merge` combines the per-shard outputs. N jobs each pay 1/N of the per-file startup *and* run 1/N of the event loop concurrently, so it attacks both costs at once.

Both `*_hist.sub` files now take `sample` and `nshards` as submit macros and queue `$(nshards)` jobs for **one** sample, with `$(ProcId)` as the shard index — so each sample gets its own ClusterId and its own shard count:

```bash
cd build && cmake .. && make clustering_hist rpt_v5_hist hist_merge
cd condor
condor_submit -a sample=vbf   -a nshards=4  clustering_hist.sub
condor_submit -a sample=dijet -a nshards=4  clustering_hist.sub
condor_submit -a sample=zjets -a nshards=12 clustering_hist.sub
# after they land -- note condor/<sample>/, NOT <repo-root>/<sample>/:
cd ..
./build/hist_merge condor/zjets/zjets_clustering_hist.root \
                   condor/zjets/zjets_clustering_hist.shard*of12.root
./build/clustering_plot --sample=zjets --hist-file=condor/zjets/zjets_clustering_hist.root
```

**Outputs land in `condor/<sample>/`.** `transfer_output_remaps`' RHS resolves
against the **submit** directory, and submission happens from `condor/`. An
earlier version of this section said `<repo-root>/<sample>/`; that is wrong,
and since `clustering_plot --sample=` looks under `OUTPUT_DIR` it needs
`--hist-file=` unless the merged file is moved.

A non-default VBS selection needs **two** matching macros, because
`resolveSelection` stamps a `SELECTION_TAG` onto the output filename (so a
loosened run cannot overwrite the standard one) and `transfer_output_files`
must name that same file:

| `selargs` | `seltag` (note trailing `_`) |
|---|---|
| *(empty)* — no \|Deta\|, m_jj >= 200 | *(empty)* |
| `--vbs-mjj=500` | `mjj500p0_` |
| `--vbs-deta=3 --vbs-mjj=0` (pre-2026-08-26 default) | `deta3p0_mjj0p0_` |
| `--vbs-deta=-1` — VBS **pair** requirement dropped | `novbs_` |

The other result-changing flags tag the name the same way (`--jvt` through
`resolveRptJvt`, the rest through `recordRunSetting` in `src/sample_config.h`;
since 2026-11-17 for `--dzpara`, `--rpt-signif` and `--fix-select-share`):

| `selargs` | `seltag` | executables |
|---|---|---|
| `--jvt=loose` (any `JVT_WPS` point but the t0-gated / TDR ones) | `jvtLoose_` (the point's tag) | `rpt_v5_hist` |
| `--dzpara` | `dzpara_` | `clustering_hist`, `rpt_v5_hist` |
| `--rpt-signif` | `rptsignif_` | `rpt_v5_hist` |
| `--clean-times` (raw ntuples: add `--ntuple-dir=`, which does not tag) | `cleant0_` | `clustering_hist` |
| `--assoc-nsigma=2.5` | `assoc2p5_` | `clustering_hist` |
| `--fix-select-share=0.3` | `fixshare0p30_` | `clustering_hist` |

Combined, the tags join in a fixed order whatever the order of the flags:
`<sel>_<jvt>_dzpara_rptsignif_` for `rpt_v5_hist` and
`<sel>_cleant0_assoc<X>_dzpara_fixshare<X>_` for `clustering_hist`, e.g.
`--vbs-mjj=500 --jvt=loose --dzpara` gives `seltag=mjj500p0_jvtLoose_dzpara_`.

**Plot a tagged run with the same flags.** `rpt_v5_plot` takes `--jvt=` like
`rpt_v5_hist`, and both plot stages take the run-setting flags: from them they
build the tag that finds the merged file, and then they check the file's own
record (`meta_vbs_jvt_wp`, `meta_vbs_<setting>`). A mismatch aborts, which also
covers a `--hist-file=` whose name says nothing. Files made before 2026-11-17
carry no setting record: a missing key reads as the default, and a non-default
flag against such a file gets a warning instead of the check.

**The pair is the leading-jet pair** (since 2026-11-17, `calcLeadingVbsPair`):
the two highest-pT jets with pT > 30 GeV and |eta| < 4.0, required in opposite
hemispheres, with no fallback to another pair. Files from that date carry
`meta_vbs_pair_lead = 1` (a `meta_vbs_` run constant, so `hist_merge` keeps it
once and checks the shards agree; `hadd` would sum it and the plot stage then
aborts), and the plot stages warn on a file without it: that file used the old
max-m_jj pair, and its untagged name cannot say so.

**A negative `--vbs-deta` is a third mode, not a smaller threshold.** It drops
the VBS candidate-*pair* requirement outright, keeping events with no valid
pair. No combination of the two thresholds expresses
this: `calcLeadingVbsPair` leaves both `dEta` and `mjj` at their `-1` sentinel when
no pair exists, so each fails any `>= 0` test — lowering a threshold relaxes only
that quantity's *magnitude* cut while the pair requirement survives. Setting
`--vbs-deta=0` was measured worth **+4% on Z+jets** (67,074 → 69,745) where
dropping the pair requirement was expected to be worth several-fold; the loose
mode is +13% on local VBF (48,314 → 54,659). Both numbers are from the old
max-m_jj pair.

`passJetPtCut` bypasses the m_jj test as well as the |Δη| one in this mode —
bypassing only |Δη| would let the default `m_jj >= 200` reject every no-pair
event on its own, reproducing the original bug exactly. Tagged `novbs` alone,
with no m_jj component, since naming a cut that is not applied would describe a
selection that does not exist. `--vbs-mjj` still rejects negatives; there is
deliberately one flag for disabling the pair, not two.

```bash
condor_submit -a sample=zjets -a nshards=12 \
              -a 'selargs=--vbs-mjj=500' \
              -a seltag=mjj500p0_ clustering_hist.sub
```

A mismatch is not silent — the job completes and then goes on hold failing to
transfer a file that does not exist — but `condor_submit -dump` shows the
expanded name, and running the executable locally with the same `selargs`
prints `[selection] ... (tag: ...)`.

**Shard counts should track EVENT count, not file count.** Beware that the
counts in the "Job cost" table above are from the 2026-07-03 runs, and the
zjets dataset has since been **replaced** (`...20260712_Output`, produced
2026-07-12): the current zjets is **2.0M events over 2000 files**, not 8.08M.
Re-measure before sizing shards for a sample you have not run recently -- the
`[chain] N files added` line plus the run's own event total give it for free.
At 2.0M, zjets wants ~12 shards, not 32; the 32-shard run below spent 49 s of
cvmfs setup against 257 s of actual work per job.

`--file-shard` tags the output filename whenever it is passed, **including at N = 1** (`Shard::given`, not `count > 1`). That looks redundant but is what lets one `.sub` template the output name on `$(ProcId)/$(nshards)` for any shard count: tying the tag to `count > 1` would make `nshards=1` write an untagged file the submit file then fails to transfer back — a foot-gun that only fires at the least-tested value.

The sort in `setupChain` is load-bearing, not tidiness: `directory_iterator` order is filesystem-defined and not guaranteed stable, so without it two shards could disagree about which file is "3rd" and silently double-process or skip data. It also makes an unsharded run deterministic, which it previously was not — harmless for full runs, but it changes *which* events a `--max-events` run sees.

**`hadd` will not do this job.** Every scalar `HistWriter` emits is a `TParameter`, and `TParameter<T>::Merge` exists: `hadd` **sums** them (ROOT 6.40.04: `Long64_t` 10 + 11 → 21, `Double_t` 1.5 + 1.5 → 3.0). That is right for the accumulators — event counts, rejection tallies, the pT accumulators `rpt_v5_plot` turns into untimed-fraction percentages — and wrong for the **run constants**, which are `TParameter`s too: `hadd` sums `meta_vbs_deta`, `meta_vbs_mjj` and `meta_vbs_jvt_wp` as well, silently. It also keeps every copy of `meta_energy_label` without checking that they agree. `util/hist_merge.cxx` sums the accumulators, keeps each run constant once after checking the inputs agree (next paragraph but one), and refuses (rather than passing through) any key type it does not recognise.

**`rpt_v5_hist` shards come in pairs.** Each shard also writes `<sample>_rpt_v5_regions.shard<i>of<N>.root`, the wide-window VBS-region tree, and `condor/rpt_v5_hist.sub` transfers and remaps both files. That file is TTree-only, so it is the export rule, not the histogram rule: merge it with `hadd -f`. `hist_merge` would refuse its key type.

**Not every scalar is additive, though.** `meta_vbs_*` describe the RUN rather than accumulate during it, so `hist_merge` treats them like `meta_energy_label`: must agree across inputs, kept once. Summing them — which is what `hadd` does — is not merely meaningless: m_jj = 200 in every shard summed to 2400 over 12 shards, which then tripped `HistReader::CheckSelection` and made every merged file unreadable by the plot stage. (Found exactly that way on the first three-sample run, when `hist_merge` still summed every scalar, by `CheckSelection` itself.) The `meta_vbs_` prefix IS the rule, so a future selection constant gets the right behaviour by being named that way.

**Full three-sample run at the new default selection** (clusters 2928810–15,
2026-08-26 09:20 CST, both `clustering_hist` and `rpt_v5_hist`): **40/40 jobs
completed, none held.** These are the current numbers — the 2026-07-03 table
below them is superseded for vbf and zjets.

| sample | files | events | accepted | shards | wall | loop/shard | cvmfs |
|---|---:|---:|---:|---:|---:|---:|---:|
| vbf   | 1001 | **1,464,400** | 630,720 (43.1%) | 4  | 21.5 min | 1092 s | 52 s |
| dijet | 1400 |   699,500     | 105,137 (15.0%) | 4  | 14.0 min |  781 s | 54 s |
| zjets | 2000 | **2,000,000** |  64,950 ( 3.2%) | 12 | 13.7 min |  759 s | 50 s |

**Every event count in the older table was wrong except dijet's.** vbf is
~2.9x larger than recorded (1.46M vs 499k) and zjets ~4x smaller (2.0M vs
8.08M) — the datasets were re-made after the July baseline. Re-measure rather
than trusting any recorded count; one `--max-events=1` run prints
`[chain] N files` and the sample's total.

**Rebalance shards accordingly**: aim for **shards ≈ events / 150,000**, which
puts every job near 10 min. On the counts above that is vbf 8–10 (not 4),
dijet 5–6, zjets 13–16. There is a ~50 s cvmfs floor per job, so splitting
below ~5 min of work per shard buys little.

A useful accident: zjets accepted **64,950** in this 12-shard default run and
**64,950** in the earlier 32-shard run that requested `--vbs-deta=0
--vbs-mjj=200` explicitly. Exact agreement confirms both that sharding is
invariant on real grid data at a different shard count, and that the new
default is precisely the selection that run asked for.

**Merging: glob on the shard count, not on `shard*`.** A sample directory can
accumulate several shardings of the same data, and `*clustering_hist.shard*.root`
matches all of them — merging that glob double-counts every event. Always pin
the count: `*clustering_hist.shard*of12.root`.

As of the 2026-09-08 prune **every stem holds exactly one sharding**, so the
glob is currently safe; it will not stay that way the next time a sample is
re-run at a different shard count. Surviving sets: dijet clustering `of4` +
rpt_v5 `of6`, ttbar rpt_v5 `of16`, vbf clustering `of4` + rpt_v5 `of4`, zjets
clustering `of12` + rpt_v5 `of14`.

**Identify a shard generation by CONTENT, not by mtime and not by event count.**
Two traps, both hit on 2026-09-08:
- `meta_n_total` counts events *read*, so **any complete sharding sums to the
  full sample** — dijet's `of4` and `of6` sets both summed to exactly 699,500,
  matching the merged file, and the scalar discriminated nothing.
- What does discriminate is the histogram set the code version wrote. The
  pre-TZP `rpt_v5_hist` shards carry 91 keys and zero `tzp`-named histograms;
  post-TZP shards and the merged files carry 103 keys including 12 `tzp` ones.
  That is what proved dijet `of4` and zjets `of12` were superseded and could
  not have been the merge source.

**`vbf/vbf_rpt_v5_hist.root` predates TZP** (2026-08-27, merged from `of4`).
dijet/zjets/ttbar were re-run 2026-09-01 with the tzp scenario; VBF was not, so
plotting it yields six scenarios, not seven. It needs a fresh submission before
VBF's RpT can be read beside the other three.

### Cleaning up after a sharded campaign

Shards and archives are intermediates. Once a merge is verified they are pure
duplication of something already on `/data`, and on `/home` they are actively
expensive. Precedent: on 2026-09-08 `condor/` held **50 G** — half the `$HOME`
quota — and cleaning it left **73 M** with no loss.

**Verify before deleting, and verify with entry counts, not sizes.** `hadd`
header overhead makes a merged file's *size* match the shard sum approximately
even when a shard was silently dropped, so size agreement proves nothing. The
check that does: open the merged file and its shards, compare
`clusters->GetEntries()` against the sum over the shard set, and require exact
equality. All 14 exports in `training_new` passed (e.g. vbf_novbs 3,994,782 =
3,994,782) before anything was removed. For an archived directory holding its
own merged copy, compare `stat -c %s` against the `/data` copy — 14 were
byte-identical.

That pass also surfaces what is NOT duplicated: three
`*_deta0p0_mjj500p0_training.root` files (vbf_mu0, ttbar_mu0, zeejets_mu0,
~2.4 G) existed only inside the archives with no `/data` copy. They were the
pre-vtx column set at a selection that `mjj500p0` now expresses identically, and
were deleted deliberately after being flagged — the point is that a wholesale
`rm -rf` of the archives would have destroyed them silently.

What must survive a cleanup: the merged `*_hist.root` files (the `*_plot` stage
inputs, ~9 M total), the PDF output subdirectories, the `.sub`/`.sh` templates,
and `logs/`.

**First sharded grid run** (cluster 2928803, zjets `clustering_hist`,
32 shards, `--vbs-deta=0 --vbs-mjj=200`): **32/32 shards completed, none held**,
end-to-end **5.3 min**. Per shard: cvmfs setup 49 s mean / 58 s max, `[chain]`
build **0.1 s for 63 files**, event loop 257 s mean / 265 s max.

Two things that confirms:

- **The startup diagnosis was right.** cvmfs is only ~50 s, so the 5–30 min of
  pre-loop dead time seen in July was almost entirely the `TChain::GetEntries()`
  file-open pass — 0.1 s for 63 files now, against 0.3–1.2 s *per file* before.
- **Sharding delivers.** The 32 shards' event loops total 137 min, which is
  what one unsharded job would have spent (same 4 threads, 32x the events),
  plus the 10–40 min serial open pass that no longer exists. So ~150–180 min
  unsharded against 5.3 min sharded: **~30x**.

**Validation** (local VBF, `util/scratch/compare_hists.py`, unsharded vs sharded-and-merged):

| executable | shards | histograms bit-identical | worst relative bin diff |
|---|---|---|---|
| `rpt_v5_hist`      | 3 | **72 / 72**       | **0** (exactly) |
| `clustering_hist`  | 4 | 461 / 477         | 5.6e-13 |
| `assoc_study_hist` | 5 | 1326 / 1378       | 1.1e-12 |

Every scalar matched exactly in all three, and shard event counts summed to the unsharded total exactly (12353+11785+11635+11762 = 47535). `assoc_study_plot`'s full numeric output was **byte-identical** between the unsharded and merged files.

The only objects that differ are `TProfile2D`s, which accumulate sums of floating-point values whose addition order changes with the partition — machine epsilon, ~12 orders of magnitude below any physics tolerance. Every other histogram is bit-identical because all fills in this codebase are unweighted, so bin contents are integer-valued sums. `rpt_v5_hist` matching exactly also independently confirms that its per-event RNG re-seeding really is independent of how work is distributed (see the Diagnostics section) — sharding is a stronger test of that than varying the thread count.

The remaining dead time is the ATLAS/cvmfs `lsetup` in `run_analysis.sh`, which now prints `[job t+NNNNs]` stamps around itself so setup and analysis are separately attributable; the executables print `[t+NNN.N s]` phase marks (`MyUtl::PhaseTimer`). Both flush immediately — a block-buffered stdout with no flushes is why a working job looked hung under `condor_tail`, which is what sent an earlier investigation after "file loading" in the first place.

**Caveat:** the build uses `-march=native` (see Compiler Flags above), tuned to whichever machine ran `make`. If execute hosts have a different CPU microarchitecture than the build host, jobs can crash with an "Illegal instruction" (SIGILL) rather than a normal error.

