"""Fig. S29 — in-vitro transformation from a known precursor: a two-stage displacement, and what the
cross-origin convergence does once a non-transforming control is available.

A  GSE240057, ten donors, paired: cosine to the malignant attractor by stage, component-abundance scoring.
B  the same under regulatory-output (activity) scoring, where the transduced-not-yet-tumour arm contributes
   nothing.
C  region membership by stage — the discrete readout moves late where the continuous one moves early.
D  cross-origin convergence in the four two-origin series, each beside its non-transforming control twin.

House rule A5: panel letters only, no titles above the axes.
"""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import figstyle as F  # noqa: E402

R = ROOT / 'perceptome_validation/results'
STAGES = ['prec_basal', 'ctrl_transduced', 'tr_early', 'tr_transition', 'tr_endpoint']
LAB = ['primary\nbasal', 'all five\nfactors,\nno tumour', 'early\ntumour', 'transition', 'endpoint\ntumour']
R2 = __import__('json').load(open(R / 'T9_attractor.json'))['regions']['R2']

base = F.use()
fig = plt.figure(figsize=F.panel_size(1.0, 1.30))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.42,
                      left=0.11, right=0.985, top=0.95, bottom=0.085)
axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])


def per_donor(path):
    d = pd.read_csv(path)
    d = d[d.arm.isin(STAGES)]
    w = d.groupby(['donor', 'arm']).cos_A.mean().unstack().reindex(columns=STAGES)
    return w


def spaghetti(ax, w, colour, note):
    x = np.arange(len(STAGES))
    for _, row in w.iterrows():
        ax.plot(x, row.values, '-', color=colour, lw=base * 0.045, alpha=0.40, zorder=2)
    m = w.mean(0).values
    ax.plot(x, m, '-o', color=colour, lw=base * 0.15, ms=base * 0.40, zorder=5)
    ax.axvspan(-0.35, 1.0, color='#fde0c5', alpha=0.55, zorder=0, lw=0)
    ax.text(0.32, ax.get_ylim()[1], note, fontsize=base * 0.74, color='#9a5b1e',
            ha='center', va='top', zorder=6)
    ax.set_xticks(x); ax.set_xticklabels(LAB, fontsize=base * 0.72)
    ax.set_xlim(-0.35, len(STAGES) - 0.6)
    ax.grid(alpha=0.20, axis='y')
    return m


# ---------------------------------------------------------------- A: component-abundance scoring
wA = per_donor(R / 'T9c_GSE240057_per_sample.csv')
axA.set_ylim(-0.08, 1.02)
mA = spaghetti(axA, wA, '#2166ac', 'before any\ntumour exists')
dA = wA.tr_endpoint - wA.prec_basal
cA = wA.ctrl_transduced - wA.prec_basal
axA.set_ylabel('cosine to the malignant attractor')
axA.annotate(f'{100*cA.mean()/dA.mean():.0f}% of the total move,\n'
             f'{int((cA>0).sum())} of {len(cA)} donors',
             xy=(1, mA[1]), xytext=(1.55, 0.10), fontsize=base * 0.78, color='0.25',
             arrowprops=dict(arrowstyle='->', lw=base * 0.055, color='0.45'))

# ---------------------------------------------------------------- B: activity scoring
wB = per_donor(R / 'T9d_GSE240057_activity_viper.csv')
# A and B plot the same quantity under two scorings, so they must carry the same units per
# pixel or the reader compares sizes against two rulers (referee item 19). A spans 1.10; B is
# centred on its own data with the same span. Tightening B to its data instead would magnify a
# null, which is the opposite error.
axB.set_ylim(-1.055, 0.045)
mB = spaghetti(axB, wB, '#b2182b', 'before any\ntumour exists')
dB = wB.tr_endpoint - wB.prec_basal
cB = wB.ctrl_transduced - wB.prec_basal
axB.set_ylabel('cosine, regulatory-output scoring', fontsize=base * 0.96)
axB.annotate(f'{int((cB>0).sum())} of {len(cB)} donors;\nno net move',
             xy=(1, mB[1]), xytext=(1.6, -0.82), fontsize=base * 0.78, color='0.25',
             arrowprops=dict(arrowstyle='->', lw=base * 0.055, color='0.45'))

# ---------------------------------------------------------------- C: region membership ladder
d = pd.read_csv(R / 'T9c_GSE240057_per_sample.csv')
d = d[d.arm.isin(STAGES)]
frac, nn = [], []
for s in STAGES:
    k = d[d.arm == s]
    frac.append(100 * k.top1.isin(R2).mean()); nn.append(len(k))
axC.bar(np.arange(len(STAGES)), frac, color='#2166ac', width=0.62,
        edgecolor='white', lw=base * 0.05)
for i, (f_, n_) in enumerate(zip(frac, nn)):
    axC.text(i, f_ + 2.5, f'{int(round(f_*n_/100))}/{n_}', ha='center',
             fontsize=base * 0.74, color='0.35')
axC.set_xticks(np.arange(len(STAGES))); axC.set_xticklabels(LAB, fontsize=base * 0.72)
axC.set_ylabel('samples inside the region (%)')
axC.set_ylim(0, 118); axC.grid(alpha=0.20, axis='y')

# ---------------------------------------------------------------- D: cross-origin, with control twins
co = pd.read_csv(R / 'T9_cross_origin.csv')
sub = __import__('json').load(open(R / 'T9e_summary.json'))['substitution_all']
g = lambda s, k='cross': float(co[(co.series == s) & (co.kind == k)].delta.iloc[0])
ROWS = [
    ('GSE222990\n2 primary origins', sub['acinar_KPTS_organoid']['delta'],
     sub['acinar_organoid']['delta'], ''),
    ('GSE118206\n2 primary epithelia', g('GSE118206'), np.nan, 'no control possible'),
    ('GSE17941\n2 immortalised lines', g('GSE17941'), g('GSE17941', 'cross_ctrl'), ''),
    ('GSE6885\n2 mammary lines', g('GSE6885'), np.nan, 'precursors 0.02 apart'),
]
y = np.arange(len(ROWS))[::-1]
h = 0.33
for yi, (lab, dlt, ctl, note) in zip(y, ROWS):
    axD.barh(yi + h / 1.8, dlt, height=h, color='#2166ac', zorder=3)
    if not np.isnan(ctl):
        axD.barh(yi - h / 1.8, ctl, height=h, color='0.68', zorder=3)
        axD.text(ctl + 0.012, yi - h / 1.8, f'{ctl:+.2f}', va='center',
                 fontsize=base * 0.70, color='0.40')
    else:
        axD.text(0.015, yi - h / 1.8, note, va='center', fontsize=base * 0.68,
                 color='0.55', style='italic')
    axD.text(dlt + 0.012, yi + h / 1.8, f'{dlt:+.2f}', va='center',
             fontsize=base * 0.70, color='#2166ac')
axD.set_yticks(y); axD.set_yticklabels([r[0] for r in ROWS], fontsize=base * 0.72)
axD.axvline(0, color='0.3', lw=base * 0.06)
axD.set_xlabel('cross-origin contraction Δ', fontsize=base * 0.98)
axD.set_xlim(-0.05, 0.66); axD.set_ylim(-0.62, 3.62); axD.grid(alpha=0.20, axis='x')
from matplotlib.patches import Patch
axD.legend(handles=[Patch(facecolor='#2166ac', label='transformed'),
                    Patch(facecolor='0.68', label='non-transforming control')],
           fontsize=base * 0.74, frameon=False, loc='lower right',
           bbox_to_anchor=(1.0, 0.02))

for ax, L in ((axA, 'A'), (axB, 'B'), (axC, 'C'), (axD, 'D')):
    ax.text(-0.19, 1.06, L, transform=ax.transAxes, fontsize=base * 1.22,
            fontweight='bold', va='bottom')

F.report('S29', F.panel_size(1.0, 1.72), 1.0)
out = ROOT / 'manuscript/figures/v11_figS29_invitro.png'
plt.savefig(out, dpi=200, bbox_inches='tight', pad_inches=0.04)
print('wrote', out)
print(f"  A: delta_endpoint {dA.mean():+.3f} ({int((dA>0).sum())}/{len(dA)}), "
      f"delta_cocktail {cA.mean():+.3f} ({int((cA>0).sum())}/{len(cA)}), f = {cA.mean()/dA.mean():.3f}")
print(f"  B: delta_endpoint {dB.mean():+.3f} ({int((dB>0).sum())}/{len(dB)}), "
      f"delta_cocktail {cB.mean():+.3f} ({int((cB>0).sum())}/{len(cB)})")
print(f"  C: region ladder {[f'{int(round(f_*n_/100))}/{n_}' for f_, n_ in zip(frac, nn)]}")
