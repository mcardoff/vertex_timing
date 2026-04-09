#!/usr/bin/env python3
"""
clustering_animation.py
Animates the iterative HGTD clustering algorithm step-by-step for one event.

Phases shown:
  1. All tracks unclustered (initial state)
  2. Seed selection — highest-pT track highlighted
  3. Distance probe — nearest unconsumed track tested against 3σ cut
  4. Merge — centroid shifts, tracks absorbed
  5. Cluster finalised — moves to next seed
  6. End — all final centroids labelled

Usage:
    python clustering_animation.py [--file_num 000001] [--event_num 2808]
                                   [--output figs/clustering_animation.mp4]

Dependencies: uproot, numpy, matplotlib (with ffmpeg for MP4 output)
"""

import argparse
import colorsys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import uproot

# ── Physics constants (mirror C++ DIST_CUT_CONE, MIN_TRACK_PT, etc.) ─────────
DIST_CUT     = 3.0   # σ  — iterative clustering distance cut
MIN_PT       = 1.0   # GeV
MAX_PT       = 30.0  # GeV
MIN_ABS_ETA  = 2.38
MAX_ABS_ETA  = 4.0
MAX_NSIGMA_Z = 3.0

# ── Colour palette ────────────────────────────────────────────────────────────
BG          = 'white'
GRID        = '#dedede'
TEXT        = '#1a1a1a'
DIM         = '#aaaaaa'
UNCLUST     = '#1f77b4'
SEED_C      = '#e6a817'   # amber — readable on white
PROBE_OK    = '#2ca02c'
PROBE_BAD   = '#d62728'
TRUTH_C     = '#1f77b4'


def golden_palette(n: int):
    phi, hue, out = 0.618033988749895, 0.07, []
    for _ in range(n):
        r, g, b = colorsys.hsv_to_rgb(hue, 0.78, 0.90)
        out.append((r, g, b))
        hue = (hue + phi) % 1.0
    return out


# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--file_num',  default='000001',
                    help='6-digit file number (default: 000001)')
parser.add_argument('--event_num', type=int, default=2808,
                    help='Event entry number within the file (default: 2808)')
parser.add_argument('--output',    default=None,
                    help='Output path (.gif or .mp4); default: figs/clustering_animation_<file>_<event>.gif')
args   = parser.parse_args()
F_NUM  = args.file_num
EVT    = args.event_num
OUT    = Path(args.output) if args.output else \
         Path(f'figs/clustering_animation_{F_NUM}_{EVT}.gif')
OUT.parent.mkdir(parents=True, exist_ok=True)


# ── Load event data ───────────────────────────────────────────────────────────
ntuple = f'../ntuple-hgtd/user.mcardiff.45809429.Output._{F_NUM}.SuperNtuple.root'
print(f'Loading {ntuple}  event {EVT}')

tree = uproot.open(ntuple)['ntuple']
b    = tree.arrays([
    'RecoVtx_z',
    'TruthVtx_z', 'TruthVtx_time', 'TruthVtx_isHS',
    'Track_z0', 'Track_var_z0', 'Track_pt', 'Track_eta',
    'Track_time', 'Track_timeRes',
    'Track_hasValidTime', 'Track_quality', 'Track_truthVtx_idx',
])

reco_z     = float(b.RecoVtx_z[EVT][0])
truth_hs_t = float(b.TruthVtx_time[EVT][0])


# ── Track selection (mirrors C++ getAssociatedTracks + makeSimpleClusters) ───
# hgtd_tracks : participate in clustering (have valid HGTD time)
# bg_tracks   : z-selected only, shown faint in background
hgtd_tracks, bg_tracks = [], []

for i in range(len(b.Track_eta[EVT])):
    eta    = float(b.Track_eta[EVT][i])
    pt     = float(b.Track_pt[EVT][i])
    z0     = float(b.Track_z0[EVT][i])
    vz0    = float(b.Track_var_z0[EVT][i])
    t      = float(b.Track_time[EVT][i])
    tres   = float(b.Track_timeRes[EVT][i])
    valid  = int(b.Track_hasValidTime[EVT][i]) == 1
    qual   = int(b.Track_quality[EVT][i]) == 1
    tvtx   = int(b.Track_truthVtx_idx[EVT][i])
    is_hs  = bool(b.TruthVtx_isHS[EVT][tvtx])

    in_hgtd  = MIN_ABS_ETA < abs(eta) < MAX_ABS_ETA
    nsigma_z = abs(z0 - reco_z) / np.sqrt(max(vz0, 1e-12))
    z_ok     = in_hgtd and MIN_PT <= pt < MAX_PT and nsigma_z < MAX_NSIGMA_Z and qual

    if z_ok and valid:
        hgtd_tracks.append(dict(
            idx=i, t=t, tres=tres,
            z0=z0, z0res=np.sqrt(vz0),
            pt=pt, is_hs=is_hs,
        ))
    elif z_ok:
        bg_tracks.append(dict(idx=i, z0=z0, pt=pt, is_hs=is_hs))

N = len(hgtd_tracks)
print(f'  {N} HGTD tracks  |  {len(bg_tracks)} z-only background tracks')


# ── Reproduce doIterativeClustering, recording every decision ─────────────────
def dist(a, b):
    return abs(a['t'] - b['t']) / np.sqrt(a['tres']**2 + b['tres']**2)

def merge(a, b):
    wa, wb   = 1/a['tres']**2, 1/b['tres']**2
    t_new    = (a['t']*wa + b['t']*wb) / (wa + wb)
    tres_new = np.sqrt(1/(wa + wb))
    return dict(t=t_new, tres=tres_new,
                pt=a['pt'] + b['pt'],
                z0_num=a['z0_num'] + b['z0_num'],
                z0_den=a['z0_den'] + b['z0_den'],
                indices=a['indices'] + b['indices'])

def trkptz(c):
    """TRKPTZ score for a cluster: TRKPT * exp(-1.5 * |mean_z0 - reco_z|)."""
    mean_z0 = c['z0_num'] / c['z0_den']
    return c['pt'] * np.exp(-1.5 * abs(mean_z0 - reco_z))

# Initialise: one cluster per track
C    = [dict(t=tr['t'], tres=tr['tres'],
             pt=tr['pt'],
             z0_num=tr['z0'] / tr['z0res']**2,
             z0_den=1.0 / tr['z0res']**2,
             indices=[tr['idx']])
        for tr in hgtd_tracks]
done   = [False] * N
finals = []   # completed clusters, in order of finalisation

steps = []  # list of step-dicts

def snap(kind, seed=-1, probe=-1, d=None, ok=None):
    steps.append(dict(
        kind=kind,
        C=[{**c} for c in C],
        done=done[:],
        finals=[{**f} for f in finals],
        seed=seed, probe=probe, d=d, ok=ok,
    ))

snap('init')

while True:
    available = [i for i in range(N) if not done[i]]
    if not available:
        break

    seed = max(available, key=lambda i: C[i]['pt'])
    snap('seed', seed=seed)

    while True:
        neighbors = [j for j in range(N) if j != seed and not done[j]]
        if not neighbors:
            break
        nn = min(neighbors, key=lambda j: dist(C[seed], C[j]))
        d  = dist(C[seed], C[nn])
        snap('probe', seed=seed, probe=nn, d=d, ok=(d < DIST_CUT))

        if d < DIST_CUT:
            C[seed] = merge(C[seed], C[nn])
            done[nn] = True
            snap('merge', seed=seed, d=d)
        else:
            break

    finals.append(C[seed].copy())
    done[seed] = True
    snap('finalize', seed=seed)

snap('end')

print(f'  {len(finals)} clusters  |  {len(steps)} animation steps')

# Palette: one colour per final cluster
pal = golden_palette(len(finals))


def get_final_color(indices_list):
    """Return the palette colour for the cluster that owns this track-index list."""
    for ci, fc in enumerate(finals):
        if fc['indices'] == indices_list:
            return pal[ci]
    return '#888888'

def track_dot_color(tr_idx, step):
    """Colour for one track dot given the current step state."""
    # In a finalised cluster?
    for fc in step['finals']:
        if tr_idx in fc['indices']:
            return get_final_color(fc['indices'])
    # In the growing seed cluster?
    seed = step['seed']
    C_now = step['C']
    if seed >= 0 and not step['done'][seed]:
        if tr_idx in C_now[seed]['indices']:
            return SEED_C
    return UNCLUST


# ── Figure setup ──────────────────────────────────────────────────────────────
# Extra bottom margin gives space for the legend + status text row below the axes
fig, ax = plt.subplots(figsize=(10, 10))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)
fig.subplots_adjust(bottom=0.18)   # room for status line + distance text below axes
ax.tick_params(labelsize=11)
ax.grid(True, color=GRID, linewidth=0.5, alpha=0.9)

all_t  = [tr['t']  for tr in hgtd_tracks]
all_pt = [tr['pt'] for tr in hgtd_tracks]
pt_pad  = max(max(all_pt) * 0.12, 0.5)
_t_buf  = (max(all_t) - min(all_t)) * 0.10
T_MIN, T_MAX = min(all_t) - _t_buf, max(all_t) + _t_buf
Z_MIN, Z_MAX = 0.5, max(all_pt) + pt_pad
T_RANGE = T_MAX - T_MIN
Z_RANGE = Z_MAX - Z_MIN

ax.set_xlim(T_MIN, T_MAX)
ax.set_ylim(Z_MIN, Z_MAX)
ax.set_yscale('log')
ax.set_xlabel('Track Time  (ps)', fontsize=13)
ax.set_ylabel('Track $p_T$  (GeV)', fontsize=13)
ax.set_title(
    f'HGTD Iterative Clustering  ·  Event {EVT}  ·  File {F_NUM}',
    fontsize=14, pad=10)

# Static: background z-only tracks — show at their pT on the y-axis (faint ticks)
for tr in bg_tracks:
    ax.plot([T_MIN, T_MIN + 0.01*T_RANGE], [tr['pt'], tr['pt']],
            color='#cccccc', linewidth=1.5, alpha=0.8, zorder=1)

# Static: truth HS time line
ax.axvline(truth_hs_t, color=TRUTH_C, linewidth=1.5,
           linestyle='--', alpha=0.6, zorder=2)

# ── Per-track artists ─────────────────────────────────────────────────────────
dots, errs = [], []
for tr in hgtd_tracks:
    dot, = ax.plot(tr['t'], tr['pt'], 'o',
                   color=UNCLUST, markersize=9, zorder=5,
                   markeredgecolor='white', markeredgewidth=0.7)
    _, caps, bars = ax.errorbar(tr['t'], tr['pt'], xerr=tr['tres'],
                                fmt='none', ecolor=DIM,
                                elinewidth=1.2, capsize=3, zorder=4)
    dots.append(dot)
    errs.append((caps, bars))

# ── Dynamic overlay artists ───────────────────────────────────────────────────
seed_ring, = ax.plot([], [], 'o', markersize=26, markerfacecolor='none',
                     markeredgecolor=SEED_C, markeredgewidth=2.5,
                     zorder=9, alpha=0)

probe_line, = ax.plot([], [], '-', linewidth=2.2, zorder=8,
                      alpha=0, color=PROBE_OK)

centroid_mk, = ax.plot([], [], '*', color=SEED_C, markersize=22,
                       markeredgecolor='#555555', markeredgewidth=0.7,
                       zorder=11, alpha=0)

# Distance result text — centered just below the axes
dist_txt = fig.text(0.5, 0.115, '', ha='center', va='top',
                    color=PROBE_OK, fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                              edgecolor=PROBE_OK, alpha=0.9))

# Status text — below the distance text, full width
status_txt = fig.text(0.5, 0.065, '', ha='center', va='top',
                      color=TEXT, fontsize=12,
                      bbox=dict(boxstyle='round,pad=0.4', facecolor='#f5f5f5',
                                edgecolor='#cccccc', alpha=0.95))

# Step counter — top-right corner of axes (no overlap with legend)
step_txt = ax.text(0.985, 0.97, '', transform=ax.transAxes,
                   color=DIM, fontsize=10, ha='right', va='top')

# Pool for artists added/removed each frame (final cluster stars, band, etc.)
dynamic_pool = []


# ── Legend — horizontal, centred below the axes, above the status text ────────
leg_els = [
    Line2D([0],[0], marker='o', color='none', markerfacecolor=UNCLUST,
           markeredgecolor='white', markersize=9, label='HGTD track (unclustered)'),
    Line2D([0],[0], marker='o', color='none', markerfacecolor=SEED_C,
           markeredgecolor='white', markersize=9, label='Track in growing cluster'),
    Line2D([0],[0], marker='*', color='none', markerfacecolor=SEED_C,
           markersize=14, label='Growing cluster centroid'),
    Line2D([0],[0], marker='*', color='none', markerfacecolor=pal[0] if pal else 'grey',
           markersize=14, label='Finalised cluster centroid'),
    Line2D([0],[0], color=TRUTH_C, linewidth=1.5, linestyle='--',
           label='Truth HS vertex time'),
    Line2D([0],[0], color=PROBE_OK,  linewidth=2, label='Probe: within 3σ (merge)'),
    Line2D([0],[0], color=PROBE_BAD, linewidth=2, label='Probe: beyond 3σ (reject)'),
]
ax.legend(handles=leg_els, loc='upper left',
          ncol=1, fontsize=9,
          facecolor='white', edgecolor=DIM, framealpha=0.95)


def centroid_y(c):
    """Mean pT of all constituent tracks — used as the y position for a centroid."""
    return c['pt'] / len(c['indices'])


# ── Frame update ──────────────────────────────────────────────────────────────
STATUS = {
    'init':     '① Initial state — all tracks unclustered',
    'seed':     '② Seed selected — highest pT track',
    'probe':    '③ Probing nearest unconsumed track',
    'merge':    '④ Within 3σ — merging, centroid shifts',
    'finalize': '⑤ No more neighbours — cluster finalised',
    'end':      '⑥ Done — all tracks clustered',
}

def update(frame_idx):
    step   = steps[frame_idx]
    kind   = step['kind']
    C_now  = step['C']
    done_n = step['done']
    seed   = step['seed']
    probe  = step['probe']

    # Remove last frame's dynamic artists
    for a in dynamic_pool:
        try: a.remove()
        except Exception: pass
    dynamic_pool.clear()

    # ── Track dot colours ─────────────────────────────────────────────────────
    for i_dot, tr in enumerate(hgtd_tracks):
        col = track_dot_color(tr['idx'], step)
        dots[i_dot].set_color(col)
        for cap in errs[i_dot][0]:
            cap.set_color(col if col != UNCLUST else DIM)
        for bar in errs[i_dot][1]:
            bar.set_color(col if col != UNCLUST else DIM)

    # ── Seed highlight ring ───────────────────────────────────────────────────
    if seed >= 0 and not done_n[seed]:
        sc = C_now[seed]
        seed_ring.set_data([sc['t']], [centroid_y(sc)])
        seed_ring.set_alpha(0.88)
    else:
        seed_ring.set_alpha(0)

    # ── Distance probe line ───────────────────────────────────────────────────
    if probe >= 0 and seed >= 0 and not done_n[seed]:
        sc, pc = C_now[seed], C_now[probe]
        probe_line.set_data([sc['t'], pc['t']],
                            [centroid_y(sc), centroid_y(pc)])
        ok = step['ok']
        c  = PROBE_OK if ok else PROBE_BAD
        probe_line.set_color(c); probe_line.set_alpha(0.95)
        dist_txt.set_text(f"d = {step['d']:.2f} σ   {'✓  MERGE' if ok else '✗  too far'}")
        dist_txt.set_color(c)
        dist_txt.get_bbox_patch().set_edgecolor(c)
    else:
        probe_line.set_alpha(0)
        dist_txt.set_text('')

    # ── Growing centroid star + 3σ acceptance band ───────────────────────────
    if seed >= 0 and not done_n[seed]:
        sc  = C_now[seed]
        cy  = centroid_y(sc)
        centroid_mk.set_data([sc['t']], [cy])
        centroid_mk.set_alpha(1.0)
        # Shade the ±3σ acceptance band (in time only)
        band_hw = DIST_CUT * sc['tres']
        rect = mpatches.Rectangle(
            (sc['t'] - band_hw, Z_MIN), 2*band_hw, Z_RANGE,
            facecolor=SEED_C, alpha=0.06, zorder=3)
        ax.add_patch(rect)
        dynamic_pool.append(rect)
        # label next to centroid
        lbl = ax.text(sc['t'] + 0.015*T_RANGE, cy,
                      f" ΣpT={sc['pt']:.1f} GeV\n σ_t={sc['tres']:.0f} ps",
                      color=SEED_C, fontsize=8.5, va='center', zorder=12,
                      bbox=dict(boxstyle='round,pad=0.2', facecolor=BG,
                                alpha=0.75, edgecolor='none'))
        dynamic_pool.append(lbl)
    else:
        centroid_mk.set_alpha(0)

    # ── Finalised cluster centroids (permanent) ───────────────────────────────
    fc_scores = [(trkptz(fc), ci) for ci, fc in enumerate(step['finals'])]
    winner_ci = max(fc_scores, key=lambda x: x[0])[1] if fc_scores else -1

    for ci, fc in enumerate(step['finals']):
        col      = pal[ci]
        cy       = centroid_y(fc)
        score    = trkptz(fc)
        is_winner = (ci == winner_ci) and (kind == 'end')
        ms    = 28 if is_winner else 22
        edge  = '#111111' if is_winner else '#333333'
        ew    = 1.8 if is_winner else 0.8
        star, = ax.plot(fc['t'], cy, '*', color=col,
                        markersize=ms, markeredgecolor=edge,
                        markeredgewidth=ew, zorder=11)
        n     = len(fc['indices'])
        label_lines = [f" n={n}  t={fc['t']:.0f} ps",
                       f" TRKPTZ={score:.2f}"]
        if is_winner:
            label_lines.append(" ◀ SELECTED")
        lbl   = ax.text(fc['t'] + 0.012*T_RANGE, cy,
                        '\n'.join(label_lines),
                        color=col, fontsize=9, va='center', zorder=12,
                        fontweight='bold' if is_winner else 'normal',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor=BG,
                                  alpha=0.80, edgecolor=col if is_winner else 'none',
                                  linewidth=1.2))
        dynamic_pool.extend([star, lbl])

    # ── Status bar ───────────────────────────────────────────────────────────
    extra = ''
    if kind == 'seed' and seed >= 0:
        extra = f'  (pT = {C_now[seed]["pt"]:.2f} GeV)'
    if kind == 'finalize' and seed >= 0:
        extra = f'  ({len(C_now[seed]["indices"])} tracks,  {C_now[seed]["pt"]:.1f} GeV)'
    if kind == 'end':
        extra = f'  ({len(finals)} clusters)'
    status_txt.set_text(STATUS.get(kind, kind) + extra)
    step_txt.set_text(f'step {frame_idx+1} / {len(steps)}')


# ── Build frame sequence (hold each step for a duration appropriate to its type)
# GIFs run at lower fps to keep file size manageable; MP4 can afford higher fps.
is_gif = OUT.suffix.lower() == '.gif'
FPS    = 12 if is_gif else 20
HOLD   = {'init': 18, 'seed': 12, 'probe': 10, 'merge': 9, 'finalize': 15, 'end': 36} \
         if is_gif else \
         {'init': 30, 'seed': 20, 'probe': 16, 'merge': 14, 'finalize': 25, 'end': 50}

frame_seq = []
for i, step in enumerate(steps):
    frame_seq.extend([i] * HOLD.get(step['kind'], 10 if is_gif else 14))

print(f'  Rendering {len(frame_seq)} frames at {FPS} fps  '
      f'(≈ {len(frame_seq)/FPS:.1f} s)  →  {OUT}')

ani = animation.FuncAnimation(
    fig, update, frames=frame_seq,
    interval=1000//FPS, blit=False, repeat=False)

if is_gif:
    writer = animation.PillowWriter(fps=FPS)
else:
    writer = animation.FFMpegWriter(
        fps=FPS, bitrate=3000,
        extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
ani.save(str(OUT), writer=writer)
plt.close(fig)
print(f'Done → {OUT}')
