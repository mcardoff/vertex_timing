#!/usr/bin/env python3
"""Event-by-event migration of the VBS pair's region under JVT/fJVT.

Joins the --jvt=none and --jvt=<wp> outputs of vbs_region_diag on
(file_idx, entry) and, for every event that is R1 / R2 / other under no
tagger, reports where it lands once the tagger runs: same region, another
region, or dropped (the event no longer has the jets to form a pair, or the
pair falls out of the m_jj window). Answers whether tagging REMOVES the
R1/R2 events -- by tagging their pileup leg -- or leaves them alone.

    python/vbs_jvt_migration.py --dir figs --sample local [--mjj 200]
"""
import argparse, os
import numpy as np
import uproot

SPLIT = 2.4
ap = argparse.ArgumentParser()
ap.add_argument('--dir', default='figs'); ap.add_argument('--sample', default='local')
ap.add_argument('--mjj', type=float, default=200.0)
args = ap.parse_args()

def region(a, m):
    """0 = not in window, 1 = R1, 2 = R2, 3 = both HS, 4 = other."""
    fA, fB = a['wide_legA_abseta'] >= SPLIT, a['wide_legB_abseta'] >= SPLIT
    hA, hB = a['wide_legA_hs'] > .5, a['wide_legB_hs'] > .5
    pA, pB = a['wide_legA_pu'] > .5, a['wide_legB_pu'] > .5
    r1 = fA & fB & ((hA & pB) | (hB & pA))
    r2 = (fA & pA & ~fB & hB) | (fB & pB & ~fA & hA)
    hh = hA & hB
    r = np.where(r1, 1, np.where(r2, 2, np.where(hh, 3, 4)))
    return np.where(m, r, 0)

def load(tag):
    t = uproot.open(os.path.join(args.dir, f'{args.sample}{tag}_vbs_region_diag.root'))['events']
    a = t.arrays(['file_idx', 'entry', 'wide_pair_mjj', 'wide_legA_abseta', 'wide_legB_abseta',
                  'wide_legA_hs', 'wide_legA_pu', 'wide_legB_hs', 'wide_legB_pu'], library='np')
    key = a['file_idx'].astype(np.int64) * 10_000_000 + a['entry'].astype(np.int64)
    return dict(zip(key, region(a, a['wide_pair_mjj'] >= args.mjj)))

NAMES = {1: 'R1', 2: 'R2', 3: 'both HS', 4: 'other', 0: 'out of window', -1: 'event dropped'}
base = load('')
for wp, tag in (('loose', '_jvtLoose'), ('tight', '_jvtTight')):
    try:
        after = load(tag)
    except FileNotFoundError:
        continue
    print(f'\n{args.sample}, none -> {wp}   (m_jj >= {args.mjj:.0f}; rows = region under no tagger, cols = after)')
    cols = [1, 2, 3, 4, 0, -1]
    print(f'{"":14s}' + ''.join(f'{NAMES[c]:>15s}' for c in cols) + f'{"n":>9s}')
    for r in (1, 2, 3, 4):
        keys = [k for k, v in base.items() if v == r]
        dest = np.array([after.get(k, -1) for k in keys])
        n = len(dest)
        print(f'{NAMES[r]:14s}' + ''.join(f'{100.0*(dest==c).mean():>14.1f}%' for c in cols) + f'{n:>9,d}')
    # and the reverse question: where do the AFTER R1/R2 come from?
    print(f'  composition of the post-tagger R1/R2 by their pre-tagger region:')
    for r in (1, 2):
        keys = [k for k, v in after.items() if v == r]
        src = np.array([base.get(k, -1) for k in keys])
        parts = ', '.join(f'{NAMES[c]} {100.0*(src==c).mean():.1f}%' for c in (1, 2, 3, 4, 0) if (src == c).any())
        print(f'    {NAMES[r]} ({len(keys):,}): {parts}')
