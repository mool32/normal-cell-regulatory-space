"""Fig S26 — the coordinate moves under therapy, with the transcriptome-complexity control shown.

Two panels (house rule A5: panel letters only, no titles above the axes):
  A  the course, as shipped (QC floor 500) against three common absolute detected-gene windows, so the
     reader sees directly how much of the level is the threshold and how much is the trajectory;
  B  core fraction by detected-gene decile within each timepoint — the control E26 never ran.
"""
import sys, os
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import figstyle as F  # noqa: E402

R = ROOT / 'results/referee2_newdata'
d = pd.read_csv(R / 'E26_dynamics.csv').sort_values('day')
win = pd.read_csv(R / 'E26b_window.csv')
dec = pd.read_csv(R / 'E26b_decile.csv')

base = F.use()
fig, (axA, axB) = plt.subplots(1, 2, figsize=F.panel_size(1.0, 2.30),
                               gridspec_kw=dict(width_ratios=[1.18, 1.0]))

# ---------------------------------------------------------------- A: the course
axA.plot(d.day, d.SHAPE_core, '--o', color='0.55', lw=base * 0.10, ms=base * 0.34,
         label='as shipped (QC floor 500)', zorder=3)
WINS = [(1000, 3000, '#92c5de'), (1500, 6000, '#2166ac'), (2000, 8000, '#053061')]
for lo, hi, c in WINS:
    w = win[(win.lo == lo) & (win.hi == hi)].sort_values('day')
    axA.plot(w.day, w.SHAPE, '-o', color=c, lw=base * 0.13, ms=base * 0.36,
             label=f'{lo}–{hi} detected genes', zorder=5)
w = win[(win.lo == 1500) & (win.hi == 6000)].sort_values('day')
axA.axhline(w.SHAPE.iloc[0], color='#2166ac', ls=':', lw=base * 0.07, zorder=1)
for _, r in w.iterrows():                      # n for the primary window only
    dx, ha = (6, 'left') if r.day == 0 else (0, 'center')
    axA.annotate(f'n={int(r.n)}', (r.day, r.SHAPE), xytext=(dx, base * 0.62),
                 textcoords='offset points', ha=ha, fontsize=base * 0.72, color='0.45', zorder=8,
                 bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
axA.set_xlabel('days of targeted therapy (erlotinib)')
axA.set_ylabel('shape-only core fraction (%)')
axA.legend(fontsize=base * 0.80, frameon=False, loc='lower center', ncol=2,
           handlelength=1.3, columnspacing=1.1, bbox_to_anchor=(0.52, -0.02))
axA.set_ylim(8, 74); axA.set_xlim(-0.7, 12.4); axA.grid(alpha=0.22)

# ---------------------------------------------------------------- B: the complexity control
days = sorted(dec.day.unique())
cmap = plt.get_cmap('viridis')
for i, day in enumerate(days):
    s = dec[dec.day == day].sort_values('decile')
    rho = spearmanr(s.decile, s.SHAPE).statistic
    axB.plot(s.decile, s.SHAPE, '-o', color=cmap(i / max(1, len(days) - 1)),
             lw=base * 0.11, ms=base * 0.30, label=f'D{day}  ρ={rho:+.2f}')
axB.set_xlabel('detected-gene decile within timepoint')
axB.set_ylabel('shape-only core fraction (%)')
axB.set_xticks(range(1, 11))
axB.legend(fontsize=base * 0.76, frameon=False, loc='upper left', ncol=2,
           handlelength=1.3, columnspacing=1.0)
axB.set_ylim(-3, 84); axB.grid(alpha=0.22)

for ax, L in ((axA, 'A'), (axB, 'B')):
    ax.text(-0.10, 1.02, L, transform=ax.transAxes, fontsize=base * 1.25, fontweight='bold', va='bottom')

F.report('S26 A (course)', F.panel_size(1.0, 2.30), 1.0)
plt.tight_layout()
out = ROOT / 'manuscript/figures/v7_figS26_dynamics.png'
plt.savefig(out, dpi=200, bbox_inches='tight')
print('wrote', out)
