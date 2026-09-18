#!/usr/bin/env python3
"""The two-sample 21-cell VBS pair composition figure (ranked bars + 6x6
half-matrices), for any set of vbs_region_diag outputs -- e.g. the JVT/fJVT
tagged ones. Generalises figs/diagnostics/vbs_region_diag/plot_decomposition.py,
which was hard-wired to the untagged grid files.

    python/vbs_pair_composition_plot.py \
        --file zjets=condor/zjets/zjets_jvtLoose_vbs_region_diag.root \
        --file vbf=figs/local_jvtLoose_vbs_region_diag.root \
        --label zjets="Z+jets" --label vbf="VBF (local)" \
        --tag "JVT + fJVT loose applied before pairing" \
        --out figs/vbs_pair_composition_jvtLoose [--mjj 500 --deta 2.5]
"""
import argparse
import numpy as np
import uproot
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

SPLIT = 2.4
ORDER = ['F-HS', 'F-PU', 'F-X', 'C-HS', 'C-PU', 'C-X']
RANK  = {k: i for i, k in enumerate(ORDER)}
NAMED = {('F-HS', 'F-PU'): 'R1', ('F-PU', 'C-HS'): 'R2'}
C_HS, C_PU, C_X, C_STEEL = '#2E7D5B', '#78838F', '#B3B9BF', '#2E5F8C'

ap = argparse.ArgumentParser()
ap.add_argument('--file', action='append', required=True, help='key=path, in drawing order')
ap.add_argument('--label', action='append', default=[], help='key=display label')
ap.add_argument('--tag', default='', help='subtitle fragment describing the tagging')
ap.add_argument('--mjj', type=float, default=500.0)
ap.add_argument('--deta', type=float, default=2.5)
ap.add_argument('--out', required=True, help='output stem')
args = ap.parse_args()
files  = [f.split('=', 1) for f in args.file]
labels = dict(l.split('=', 1) for l in args.label)
SAMPLES = [(k, labels.get(k, k)) for k, _ in files]
PATHS   = dict(files)


def load(path):
    t = uproot.open(path)['events']
    a = t.arrays(['wide_pair_mjj', 'wide_pair_deta',
                  'wide_legA_hs', 'wide_legA_pu', 'wide_legA_abseta',
                  'wide_legB_hs', 'wide_legB_pu', 'wide_legB_abseta'], library='np')
    def cat(p):
        ae  = a[f'wide_leg{p}_abseta']
        reg = np.where(ae >= SPLIT, 'F', 'C')
        idn = np.where(a[f'wide_leg{p}_hs'] > 0.5, 'HS',
              np.where(a[f'wide_leg{p}_pu'] > 0.5, 'PU', 'X'))
        return np.char.add(np.char.add(reg, '-'), idn)
    ra = np.array([RANK[x] for x in cat('A')])
    rb = np.array([RANK[x] for x in cat('B')])
    sw = rb < ra
    lo, hi = np.where(sw, rb, ra), np.where(sw, ra, rb)
    m = (a['wide_pair_mjj'] >= args.mjj) & (a['wide_pair_deta'] >= args.deta)
    n = int(m.sum())
    M = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            M[i, j] = 100.0 * ((lo == i) & (hi == j) & m).sum() / max(n, 1)
    return M, n, int(len(m))


data = {s: load(PATHS[s]) for s, _ in SAMPLES}

cells = []
for i in range(6):
    for j in range(i, 6):
        f = {s: data[s][0][i, j] for s, _ in SAMPLES}
        if max(f.values()) < 0.05:
            continue
        cells.append((ORDER[i], ORDER[j], f))
cells.sort(key=lambda c: -sum(c[2].values()))


def identity_colour(a, b):
    ia, ib = a.split('-')[1], b.split('-')[1]
    if 'X' in (ia, ib):           return C_X
    if ia == 'HS' and ib == 'HS': return C_HS
    if ia == 'PU' and ib == 'PU': return C_PU
    return C_STEEL


plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.linewidth': 0.8, 'xtick.direction': 'in',
    'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'axes.labelsize': 10, 'figure.dpi': 160,
})
fig = plt.figure(figsize=(12.4, 7.4))
gs  = fig.add_gridspec(2, 2, height_ratios=[1.42, 1.0], width_ratios=[1, 1],
                       hspace=0.30, wspace=0.13, left=0.115, right=0.985, top=0.885, bottom=0.075)

axA = fig.add_subplot(gs[0, :])
y, h = np.arange(len(cells))[::-1], 0.38
xmax = max(max(c[2].values()) for c in cells) * 1.15 + 3
for k, (s, lbl) in enumerate(SAMPLES):
    vals = [c[2][s] for c in cells]
    cols = [identity_colour(c[0], c[1]) for c in cells]
    off  = h / 2 if k == 0 else -h / 2
    axA.barh(y + off, vals, height=h, color=cols, edgecolor='#1B2027', linewidth=0.5,
             alpha=1.0 if k == 0 else 0.42, zorder=3)
    for yy, v in zip(y + off, vals):
        if v >= 1.0:
            axA.text(v + 0.6, yy, f'{v:.1f}', va='center', ha='left', fontsize=7.4, color='#3D4753', zorder=4)
lab = []
for a, b, _ in cells:
    tag = NAMED.get((a, b)) or NAMED.get((b, a))
    lab.append(f'{a} + {b}' + (f'   [{tag}]' if tag else ''))
axA.set_yticks(y); axA.set_yticklabels(lab, fontsize=8.6, family='monospace')
for t, (a, b, _) in zip(axA.get_yticklabels(), cells):
    if NAMED.get((a, b)) or NAMED.get((b, a)):
        t.set_color(C_STEEL); t.set_fontweight('bold')
axA.set_xlabel('fraction of VBS pairs  [%]        (cells below 0.05% in both samples omitted)')
axA.set_xlim(0, xmax); axA.set_ylim(-0.9, len(cells) - 0.1)
axA.grid(axis='x', color='#E4E8EC', lw=0.7, zorder=0); axA.set_axisbelow(True)
for sp in ('top', 'right'): axA.spines[sp].set_visible(False)
axA.tick_params(top=False, right=False)
for yy, (a, b, _) in zip(y, cells):
    if NAMED.get((a, b)) or NAMED.get((b, a)):
        axA.add_patch(Rectangle((0, yy - 0.48), xmax, 0.96, facecolor=C_STEEL, alpha=0.07, zorder=1, lw=0))
leg1 = axA.legend(handles=[
    Line2D([], [], marker='s', ls='', ms=8, mfc='#55606B', mec='#1B2027', label=f'{SAMPLES[0][1]}  (solid)'),
    Line2D([], [], marker='s', ls='', ms=8, mfc='#55606B', mec='#1B2027', alpha=0.42, label=f'{SAMPLES[1][1]}  (faded)')],
    loc='lower right', frameon=False, fontsize=8.6, handletextpad=0.5, bbox_to_anchor=(1.0, 0.02))
axA.add_artist(leg1)
axA.legend(handles=[
    Line2D([], [], marker='s', ls='', ms=8, mfc=C_HS,    mec='none', label='HS + HS'),
    Line2D([], [], marker='s', ls='', ms=8, mfc=C_STEEL, mec='none', label='HS + PU'),
    Line2D([], [], marker='s', ls='', ms=8, mfc=C_PU,    mec='none', label='PU + PU'),
    Line2D([], [], marker='s', ls='', ms=8, mfc=C_X,     mec='none', label='either leg unlabelled')],
    loc='lower right', frameon=False, fontsize=8.6, handletextpad=0.5, ncol=1, bbox_to_anchor=(1.0, 0.20))

for k, (s, lbl) in enumerate(SAMPLES):
    ax = fig.add_subplot(gs[1, k])
    M, n, nsel = data[s]
    Mm = np.ma.masked_where(np.tril(np.ones((6, 6)), -1) > 0, M)
    ax.imshow(Mm ** 0.5, cmap='Blues', vmin=0, vmax=np.sqrt(M.max()), interpolation='nearest', aspect='auto')
    for i in range(6):
        for j in range(i, 6):
            v = M[i, j]
            ax.text(j, i, f'{v:.1f}' if v >= 0.05 else '·', ha='center', va='center', fontsize=7.8,
                    family='monospace', color='white' if np.sqrt(v) > 0.58 * np.sqrt(M.max()) else '#20262D')
            nm = NAMED.get((ORDER[i], ORDER[j])) or NAMED.get((ORDER[j], ORDER[i]))
            if nm:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor='#B4522E', lw=1.9, zorder=5))
                ax.text(j - 0.42, i - 0.34, nm, fontsize=6.4, color='#B4522E', fontweight='bold', zorder=6)
    ax.set_xticks(range(6)); ax.set_xticklabels(ORDER, fontsize=8, family='monospace')
    ax.set_yticks(range(6)); ax.set_yticklabels(ORDER, fontsize=8, family='monospace')
    ax.set_title(f'{lbl}   —   {n:,} pairs  ({100.0 * n / nsel:.0f}% of selected)', fontsize=9.6, pad=7, loc='left')
    ax.tick_params(length=0, top=False, right=False)
    for sp in ax.spines.values(): sp.set_visible(False)

fig.suptitle('VBS pair composition — truth identity × detector region, both legs'
             + (f'   ·   {args.tag}' if args.tag else ''),
             x=0.115, ha='left', y=0.968, fontsize=13, fontweight='bold', color='#161A20')
fig.text(0.115, 0.925,
         f'forward = |$\\eta$| $\\geq$ {SPLIT} (no upper edge)   ·   $m_{{jj}} \\geq$ {args.mjj:.0f} GeV'
         + (f', $|\\Delta\\eta| \\geq$ {args.deta}' if args.deta > 0 else '')
         + '   ·   paper HS / PU truth labels   ·   lepton–jet overlap removal ON',
         ha='left', fontsize=9, color='#69747F')
for ext in ('pdf', 'png'):
    fig.savefig(f'{args.out}.{ext}', bbox_inches='tight', facecolor='white')
print('wrote', args.out + '.png')
