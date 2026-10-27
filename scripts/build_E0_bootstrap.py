#!/usr/bin/env python3
"""Rebuilt from the reference itself: how stable is each principal axis under resampling of the cell types?

The previous PNG had no builder AND no data file anywhere in the repository — its numbers lived only in the
caption and in the image. This recomputes them: resample the 154 reference cell types with replacement,
rebuild the eigenspace by the canonical recipe, and record |cos| between each bootstrap eigenvector and the
corresponding full-data one. Writes results/E0_bootstrap_stability.csv and the figure.

Panel letters only; nothing above the axes (house rule A5).
"""
import os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'perceptome_validation' / 'scripts'))
import figstyle as F            # noqa: E402
import prereg_lib as L          # noqa: E402

NBOOT = 500
RNG = np.random.default_rng(0)
OUT = ROOT / 'results' / 'E0_bootstrap_stability.csv'


def main():
    panel = L.perceptome_panel()
    hpa = L.hpa_gene_level()
    R = L.score_ref(hpa['H_pb'], hpa['genes'], panel)          # 154 x 44, the eigenspace calibration
    full = L.eigenspace(R, n_pc=L.NPC)
    V0 = full['L']                                             # 44 x 9
    n = R.shape[0]
    print(f'reference {R.shape}, {NBOOT} bootstraps of the {n} cell types')

    cos = np.zeros((NBOOT, L.NPC))
    for b in range(NBOOT):
        idx = RNG.integers(0, n, n)
        try:
            e = L.eigenspace(R[idx], n_pc=L.NPC)
        except Exception:
            cos[b] = np.nan; continue
        # match each full-data axis to its best bootstrap partner, then record |cos| for that pairing
        M = np.abs(V0.T @ e['L'])                              # 9 x 9
        used = set()
        for i in range(L.NPC):
            j = int(np.argmax([M[i, k] if k not in used else -1 for k in range(L.NPC)]))
            used.add(j); cos[b, i] = M[i, j]
        if (b + 1) % 100 == 0:
            print(f'  {b+1}/{NBOOT}', flush=True)

    df = pd.DataFrame(cos, columns=[f'PC{i+1}' for i in range(L.NPC)])
    df.to_csv(OUT, index=False)
    m = df.mean(); lo = df.quantile(0.05); hi = df.quantile(0.95)
    print('\nmean |cos| to the full-data axis (5th-95th percentile):')
    for c in df.columns:
        print(f'  {c}: {m[c]:.2f}  [{lo[c]:.2f}, {hi[c]:.2f}]')

    BASE = F.use()
    fig, ax = plt.subplots(figsize=F.panel_size(1.0, 2.60))
    x = np.arange(L.NPC)
    colours = ['#2166ac' if v >= 0.80 else '#f0a202' if v >= 0.60 else '#b2182b' for v in m]
    ax.bar(x, m.values, color=colours, width=0.68,
           yerr=[m.values - lo.values, hi.values - m.values], capsize=BASE * 0.18,
           error_kw=dict(lw=BASE * 0.07, ecolor='0.35'))
    for i, v in enumerate(m.values):
        ax.text(i, v + 0.035, f'{v:.2f}', ha='center', fontsize=BASE * 0.80, color='0.25')
    ax.axhline(0.80, color='0.45', ls=':', lw=BASE * 0.07)
    ax.text(L.NPC - 0.45, 0.805, 'stable', fontsize=BASE * 0.74, color='0.45', va='bottom', ha='right')
    ax.set_xticks(x); ax.set_xticklabels(df.columns)
    ax.set_ylabel('|cos| to the full-data axis\n(mean of 500 bootstraps, 5–95%)', fontsize=BASE * 0.95)
    ax.set_xlabel('principal axis of the reference', fontsize=BASE * 0.95)
    ax.set_ylim(0, 1.12)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    F.report('E0 bootstrap', F.panel_size(1.0, 2.60), 1.0)
    plt.tight_layout()
    p = ROOT / 'manuscript/figures/v7_figS1_bootstrap.png'
    plt.savefig(p, dpi=200, bbox_inches='tight')
    print('wrote', p, '\nwrote', OUT)


if __name__ == '__main__':
    main()
