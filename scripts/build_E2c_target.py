#!/usr/bin/env python3
"""Rebuilt from data: how much of the convergence is carried by the gene SET and how much by the GROUPING.

The previous PNG had no builder in the repository, carried a suptitle asserting that the convergence "is
driven by the biology of the selected pathway genes, not by any gene set", and titled its first panel
"Convergence TARGET is specific to the curated gene set" — while the verdict written by the analysis that
produced its own data (results/E2c_target_out.txt) is NOT-specific on both metrics. This version plots the
same two nulls and states what they show.

Reads results/E2c_target_null.csv (150 draws of each null) and the real values from E2c_target_out.txt.
Panel letters only, nothing above the axes (house rule A5).
"""
import os, re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import figstyle as F  # noqa: E402

d = pd.read_csv(ROOT / 'results/E2c_target_null.csv')
txt = (ROOT / 'results/E2c_target_out.txt').read_text()
REAL_HIT = float(re.search(r'trophoblast-hit=([\d.]+)%', txt).group(1))
REAL_STAY = float(re.search(r'Stayton-ratio=([\d.]+)', txt).group(1))
NCELL = len(re.search(r'target set \((\d+) cells\)', txt).group(1)) and int(
    re.search(r'target set \((\d+) cells\)', txt).group(1))

BASE = F.use()
fig, (axA, axB) = plt.subplots(1, 2, figsize=F.panel_size(1.0, 3.05))

GREY, BLUE, RED = '#b6b6b6', '#6ba3d6', '#c0392b'


def panel(ax, rand, shuf, real, xlab, real_lab):
    lo = min(rand.min(), shuf.min(), real); hi = max(rand.max(), shuf.max(), real)
    bins = np.linspace(lo, hi, 26)
    ax.hist(rand, bins=bins, color=GREY, alpha=0.85, label='random-genome panels')
    ax.hist(shuf, bins=bins, color=BLUE, alpha=0.70, label='same genes, shuffled into random modules')
    ax.axvline(real, color=RED, lw=BASE * 0.16, zorder=6, label=real_lab)
    ax.set_xlabel(xlab, fontsize=BASE * 0.95)
    ax.set_ylabel('null catalogs', fontsize=BASE * 0.95)
    ax.legend(fontsize=BASE * 0.76, frameon=False, loc='upper right',
              handlelength=1.1, labelspacing=0.3)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


panel(axA, 100 * d.rand_hit.values, 100 * d.shuf_hit.values, REAL_HIT,
      'malignant cells reaching the target set (%)', f'real 44 modules ({REAL_HIT:.0f}%)')
axA.annotate(f'vs random genome  z = +3.1,  P = 0.013\n'
             f'vs shuffled grouping  z = +1.3,  P = 0.11  (n.s.)',
             xy=(0.03, 0.40), xycoords='axes fraction', fontsize=BASE * 0.78, color='0.25',
             bbox=dict(boxstyle='round,pad=0.3', fc='#f7f7f5', ec='#ccc', lw=0.6))

panel(axB, d.rand_stay.values, d.shuf_stay.values, REAL_STAY,
      'convergence-concentration ratio', f'real ({REAL_STAY:.1f})')
axB.annotate('the real catalog concentrates\nLESS than either null\n(P = 0.92, 0.93)',
             xy=(0.52, 0.55), xycoords='axes fraction', fontsize=BASE * 0.78, color='0.25',
             bbox=dict(boxstyle='round,pad=0.3', fc='#f7f7f5', ec='#ccc', lw=0.6))

for ax, L in ((axA, 'A'), (axB, 'B')):
    ax.text(-0.12, 1.03, L, transform=ax.transAxes, fontsize=BASE * 1.2, fontweight='bold', va='bottom')

F.report('E2c target-specificity', F.panel_size(1.0, 3.05), 1.0)
plt.tight_layout()
out = ROOT / 'manuscript/figures/v7_fig_E2c_target.png'
plt.savefig(out, dpi=200, bbox_inches='tight')
print('wrote', out, f'| target set = {NCELL} cells, 17 cancers')
