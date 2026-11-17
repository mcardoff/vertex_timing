#!/usr/bin/env python3
"""VBS pair composition under JVT/fJVT working points.

Reads the three vbs_region_diag outputs of one sample (--jvt=none / loose /
tight), re-bins the wide_ block's forward/central split at |eta| = 2.4 with no
upper edge (as results/vbs_pair_composition.md), and shows how the 21-cell
truth-identity x region composition of the chosen tagging pair moves when the
standard pileup-jet taggers are applied BEFORE pairing. Also prints the
marginals per working point and what each tagger removed, by truth identity.

    python/vbs_jvt_plot.py --dir figs --sample local --label "VBF (local)" [--mjj 200]
    python/vbs_jvt_plot.py --dir condor/zjets --sample zjets --label "Z+jets"

Expects <dir>/<sample>_vbs_region_diag.root, <sample>_jvtLoose_..., <sample>_jvtTight_...
"""
import argparse, os, sys
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
WPS   = [('none', '', 'no tagger'), ('loose', '_jvtLoose', 'JVT + fJVT loose'),
         ('tight', '_jvtTight', 'JVT + fJVT tight')]
C_HS, C_PU, C_X, C_STEEL = '#2E7D5B', '#78838F', '#B3B9BF', '#2E5F8C'
WP_ALPHA = {'none': 1.0, 'loose': 0.62, 'tight': 0.32}

ap = argparse.ArgumentParser()
ap.add_argument('--dir', default='figs')
ap.add_argument('--sample', default='local')
ap.add_argument('--label', default=None)
ap.add_argument('--mjj', type=float, default=200.0)
ap.add_argument('--deta', type=float, default=0.0)
ap.add_argument('--out', default=None, help='output stem (default <dir>/<sample>_vbs_jvt)')
args = ap.parse_args()
label = args.label or args.sample
stem  = args.out or os.path.join(args.dir, f'{args.sample}_vbs_jvt')


def load(tag):
    fn = os.path.join(args.dir, f'{args.sample}{tag}_vbs_region_diag.root')
    t = uproot.open(fn)['events']
    a = t.arrays(['wide_pair_mjj', 'wide_pair_deta',
                  'wide_legA_hs', 'wide_legA_pu', 'wide_legA_abseta',
                  'wide_legB_hs', 'wide_legB_pu', 'wide_legB_abseta',
                  'n_rm_jvt', 'n_rm_jvt_hs', 'n_rm_jvt_pu',
                  'n_rm_fjvt', 'n_rm_fjvt_hs', 'n_rm_fjvt_pu'], library='np')
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
    nhs = np.isin(lo, [0, 3]).astype(int) + np.isin(hi, [0, 3]).astype(int)
    nf  = (lo < 3).astype(int) + (hi < 3).astype(int)
    marg = {
        'R1': M[RANK['F-HS'], RANK['F-PU']], 'R2': M[RANK['F-PU'], RANK['C-HS']],
        'no HS leg': 100.0 * (nhs[m] == 0).mean(), 'one HS': 100.0 * (nhs[m] == 1).mean(),
        'two HS': 100.0 * (nhs[m] == 2).mean(),
        'PU+PU (any region)': sum(M[RANK[x], RANK[y]] for x in ('F-PU', 'C-PU') for y in ('F-PU', 'C-PU') if RANK[x] <= RANK[y]),
        'both fwd': 100.0 * (nf[m] == 2).mean(),
    }
    rm = {k: int(a[k].sum()) for k in ('n_rm_jvt', 'n_rm_jvt_hs', 'n_rm_jvt_pu',
                                        'n_rm_fjvt', 'n_rm_fjvt_hs', 'n_rm_fjvt_pu')}
    return dict(M=M, n=n, nsel=int(len(m)), marg=marg, rm=rm, fn=fn)


data = {}
for wp, tag, _ in WPS:
    try:
        data[wp] = load(tag)
    except FileNotFoundError as e:
        print(f'[skip] {wp}: {e}', file=sys.stderr)
if 'none' not in data:
    sys.exit('need at least the --jvt=none file')
wps = [w for w in WPS if w[0] in data]

# ── console: marginals + removals ────────────────────────────────────────────
print(f'\n{label}: wide_ block, m_jj >= {args.mjj:.0f} GeV, |d-eta| >= {args.deta}')
print(f'{"":22s}' + ''.join(f'{d:>18s}' for _, _, d in wps))
print(f'{"selected events":22s}' + ''.join(f'{data[w]["nsel"]:>18,d}' for w, _, _ in wps))
print(f'{"pairs in window":22s}' + ''.join(f'{data[w]["n"]:>18,d}' for w, _, _ in wps))
for k in data['none']['marg']:
    print(f'{k:22s}' + ''.join(f'{data[w]["marg"][k]:>17.2f}%' for w, _, _ in wps))
print('\nremoved by the tagger (jets, over selected events; HS/PU by paper label):')
for k in ('n_rm_jvt', 'n_rm_jvt_hs', 'n_rm_jvt_pu', 'n_rm_fjvt', 'n_rm_fjvt_hs', 'n_rm_fjvt_pu'):
    print(f'{k:22s}' + ''.join(f'{data[w]["rm"][k]:>18,d}' for w, _, _ in wps))

# ── figure: ranked cells, one bar per working point ─────────────────────────
cells = []
for i in range(6):
    for j in range(i, 6):
        f = {w: data[w]['M'][i, j] for w, _, _ in wps}
        if max(f.values()) < 0.05:
            continue
        cells.append((ORDER[i], ORDER[j], f))
cells.sort(key=lambda c: -c[2]['none'])

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
fig, ax = plt.subplots(figsize=(10.5, 0.36 * len(cells) + 2.2))
fig.subplots_adjust(left=0.20, right=0.98, top=0.86, bottom=0.10)
y = np.arange(len(cells))[::-1]
nw = len(wps)
h = 0.8 / nw
xmax = max(max(c[2].values()) for c in cells) * 1.18
for k, (w, _, desc) in enumerate(wps):
    vals = [c[2][w] for c in cells]
    cols = [identity_colour(c[0], c[1]) for c in cells]
    off = (nw - 1) / 2 * h - k * h
    ax.barh(y + off, vals, height=h, color=cols, edgecolor='#1B2027', linewidth=0.5,
            alpha=WP_ALPHA[w], zorder=3)
    for yy, v in zip(y + off, vals):
        if v >= 1.0:
            ax.text(v + 0.4, yy, f'{v:.1f}', va='center', ha='left', fontsize=6.8,
                    color='#3D4753', zorder=4)
labels = []
for a, b, _ in cells:
    tag = NAMED.get((a, b)) or NAMED.get((b, a))
    labels.append(f'{a} + {b}' + (f'   [{tag}]' if tag else ''))
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.4, family='monospace')
for t, (a, b, _) in zip(ax.get_yticklabels(), cells):
    if NAMED.get((a, b)) or NAMED.get((b, a)):
        t.set_color(C_STEEL); t.set_fontweight('bold')
for yy, (a, b, _) in zip(y, cells):
    if NAMED.get((a, b)) or NAMED.get((b, a)):
        ax.add_patch(Rectangle((0, yy - 0.48), xmax, 0.96, facecolor=C_STEEL, alpha=0.07, zorder=1, lw=0))
ax.set_xlim(0, xmax); ax.set_ylim(-0.9, len(cells) - 0.1)
ax.set_xlabel('fraction of VBS pairs  [%]        (cells below 0.05% under every working point omitted)')
ax.grid(axis='x', color='#E4E8EC', lw=0.7, zorder=0); ax.set_axisbelow(True)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.tick_params(top=False, right=False)
leg1 = ax.legend(handles=[Line2D([], [], marker='s', ls='', ms=8, mfc='#55606B', mec='#1B2027',
                                 alpha=WP_ALPHA[w], label=f'{desc}   ({data[w]["n"]:,} pairs)')
                          for w, _, desc in wps],
                 loc='lower right', frameon=False, fontsize=8.4, handletextpad=0.5)
ax.add_artist(leg1)
ax.legend(handles=[Line2D([], [], marker='s', ls='', ms=8, mfc=C_HS, mec='none', label='HS + HS'),
                   Line2D([], [], marker='s', ls='', ms=8, mfc=C_STEEL, mec='none', label='HS + PU'),
                   Line2D([], [], marker='s', ls='', ms=8, mfc=C_PU, mec='none', label='PU + PU'),
                   Line2D([], [], marker='s', ls='', ms=8, mfc=C_X, mec='none', label='either leg unlabelled')],
          loc='lower right', frameon=False, fontsize=8.4, handletextpad=0.5,
          bbox_to_anchor=(1.0, 0.16 + 0.02 * nw))
fig.suptitle(f'{label}: VBS pair composition with pileup-jet tagging applied before pairing',
             x=0.20, ha='left', y=0.975, fontsize=12.5, fontweight='bold', color='#161A20')
fig.text(0.20, 0.925,
         f'forward = |$\\eta$| $\\geq$ {SPLIT} (no upper edge)   ·   $m_{{jj}} \\geq$ {args.mjj:.0f} GeV'
         + (f', $|\\Delta\\eta| \\geq$ {args.deta}' if args.deta > 0 else '')
         + '   ·   JVT: $R_{pT}$ proxy (|$\\eta$|<2.5, $p_T$<60)   ·   fJVT (2.5$\\leq$|$\\eta$|<4.5, $p_T$<120)',
         ha='left', fontsize=8.6, color='#69747F')
for ext in ('pdf', 'png'):
    fig.savefig(f'{stem}.{ext}', bbox_inches='tight', facecolor='white')
print(f'\nwrote {stem}.pdf / .png')
