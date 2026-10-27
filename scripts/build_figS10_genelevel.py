"""Fig. S10 — the per-gene suppression signature, redrawn.

Two reasons this needed a new generator. (1) The producer of v7_fig_E4_genelevel.png is gone from the
repo; only the PNG survived, so the figure could not be corrected at all. (2) The surviving render put
one point at cancer_z = +62 against a y-axis otherwise spanning -2 to +3, compressing every other gene
into a sliver. That point is UGT1A1, a liver UDP-glucuronosyltransferase that is near-zero across most
of the 154 normal cell types, so its reference SD is tiny and the z-score explodes. The published
statistic is Spearman and is untouched by it (rho 0.558 with, 0.557 without), but a reader measuring a
LINEAR correlation off the old panel would get r = 0.10 instead of 0.41 — so the outlier is now drawn
off-scale with its value stated rather than silently setting the axis.
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
sys.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix()); import figstyle as F
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path
R=Path(__file__).resolve().parents[1]; FIG=R/'manuscript/figures'
d=pd.read_csv(R/'results/E4_genelevel.csv')
rho,pv=spearmanr(d.target_distinct,d.cancer_z)
CLIP=3.2
off=d[d.cancer_z>CLIP]
BASE=F.use(); FS=F.panel_size(1.0,2.30)
fig,ax=plt.subplots(figsize=FS)
for gold,col,lab in [(True,'#1e8449','genes of internally coherent modules'),
                     (False,'#e67e22','genes of flagged modules')]:
    g=d[(d.gold==gold)&(d.cancer_z<=CLIP)]
    r=spearmanr(d[d.gold==gold].target_distinct,d[d.gold==gold].cancer_z)[0]
    ax.scatter(g.target_distinct,g.cancer_z,s=46,color=col,edgecolor='white',lw=.4,
               label=f'{lab}  (ρ = {r:+.2f})',zorder=3)
b=np.polyfit(d.target_distinct,d.cancer_z.clip(upper=CLIP),1)
xs=np.linspace(d.target_distinct.min(),d.target_distinct.max(),50)
ax.plot(xs,np.polyval(b,xs),'--',color='0.45',lw=1.8,zorder=2)
ax.axhline(0,color='0.7',lw=.8); ax.axvline(0,color='0.7',lw=.8)
ax.set_ylim(d.cancer_z.min()-0.4, CLIP)
for _,r_ in off.iterrows():
    ax.annotate(f'{r_.gene} off scale at {r_.cancer_z:.0f}',(r_.target_distinct,CLIP),
                xytext=(0,-14),textcoords='offset points',ha='center',va='top',
                fontsize=BASE*0.78,color='#8e1b1b',
                arrowprops=dict(arrowstyle='->',color='#8e1b1b',lw=1.2))
ax.set_xlabel('cytotrophoblast gene suppression (per gene, ← more suppressed)',fontsize=BASE*0.92)
ax.set_ylabel('cancer gene suppression (mean z)',fontsize=BASE*0.92)
ax.legend(fontsize=BASE*0.82,loc='lower right',frameon=False)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(alpha=.2,zorder=0)
plt.savefig(FIG/'v7_fig_E4_genelevel.png',dpi=200,bbox_inches='tight'); plt.close()
F.report('Fig S10 gene-level',FS,1.0)
print(f'Fig. S10: n={len(d)}, Spearman rho={rho:+.3f} (p={pv:.1e}); {len(off)} point(s) off scale: '
      + ', '.join(f'{r.gene} {r.cancer_z:.1f}' for _,r in off.iterrows()))
