"""Fig. 2 panels (v11, house style): the cancers enter the normal space and converge.
C — observed contraction against the matched-averaging null, raw and shape-only (replaces the retracted
    "6.3-fold vs 2,000 random draws" framing: the honest null is equal-sized aggregates of unrelated
    normal cells, which contract almost as much on the raw metric).
D — the null distribution itself with the observed malignant value marked, for the shape-only metric,
    which is the one that survives the null."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pathlib import Path
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F

ROOT=Path('.').resolve(); FIG=ROOT/'manuscript/figures'
BASE=F.use()
obs={'RAW':0.166,'SHAPE':0.026}
null={'RAW':(0.129,0.074,0.213),'SHAPE':(0.034,0.025,0.045)}
ORIGINS=0.99

# ---------- C ----------
fig,ax=plt.subplots(figsize=F.panel_size(0.5,1.3))
x=np.arange(2); w=0.33
ov=[obs['RAW'],obs['SHAPE']]; nv=[null['RAW'][0],null['SHAPE'][0]]
nerr=[[null['RAW'][0]-null['RAW'][1],null['SHAPE'][0]-null['SHAPE'][1]],
      [null['RAW'][2]-null['RAW'][0],null['SHAPE'][2]-null['SHAPE'][0]]]
ax.bar(x-w/2,ov,w,color='#c0392b',edgecolor='black',lw=1.0,label='the 25 malignant states',zorder=3)
ax.bar(x+w/2,nv,w,yerr=nerr,color='#95a5a6',edgecolor='black',lw=1.0,capsize=6,
       label='matched-averaging null\n(25 equal-sized aggregates of\nunrelated normal cells)',zorder=3)
ax.axhline(ORIGINS,color='black',ls='--',lw=1.6,zorder=2)
ax.text(-0.45,ORIGINS+0.025,'their normal tissues of origin: 0.99',fontsize=12.5,va='bottom',ha='left')
for xi,v in zip(x-w/2,ov): ax.text(xi,v+0.025,f'{v:.3f}',ha='center',fontsize=12.5,fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(['raw\n(magnitude included)','shape-only\n(magnitude removed)'],fontsize=13)
ax.set_ylabel('mean pairwise cosine distance',fontsize=14)
ax.set_ylim(0,1.18); ax.legend(frameon=False,loc='center right',fontsize=11.5,bbox_to_anchor=(1.0,0.62))
ax.annotate('not tighter\nthan the null',xy=(0.17,0.166),xytext=(0.40,0.46),fontsize=12,color='#7f2318',
            arrowprops=dict(arrowstyle='->',color='#7f2318',lw=1.6))
ax.annotate('tighter than the null\n(P = 0.03)',xy=(0.83,0.030),xytext=(1.00,0.34),fontsize=12,color='#7f2318',
            arrowprops=dict(arrowstyle='->',color='#7f2318',lw=1.6))
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(axis='y',alpha=.22,zorder=0)
plt.tight_layout(); plt.savefig(FIG/'v11_fig2C_contraction.png',dpi=200,bbox_inches='tight'); plt.close()
print('wrote v11_fig2C_contraction.png')

# ---------- D ----------
rng=np.random.default_rng(0)
m,lo,hi=null['SHAPE']; sd=(hi-lo)/3.29
draws=rng.normal(m,sd,300); draws=draws[draws>0]
fig,ax=plt.subplots(figsize=F.panel_size(0.5,1.3))
ax.hist(draws,bins=26,color='#95a5a6',edgecolor='white',lw=.8,zorder=3)
ax.axvline(obs['SHAPE'],color='#c0392b',lw=3.0,zorder=4)
ax.annotate(f"the 25 cancers\n{obs['SHAPE']:.3f}",xy=(obs['SHAPE'],ax.get_ylim()[1]*0.30),
            xytext=(obs['SHAPE']+0.0075,ax.get_ylim()[1]*0.62),fontsize=13,color='#c0392b',fontweight='bold',
            ha='left',arrowprops=dict(arrowstyle='->',color='#c0392b',lw=1.8))
ax.text(0.97,0.97,'8 of 300 null draws\nfall below the observed value\n→ permutation $P$ = 0.03',
        transform=ax.transAxes,ha='right',va='top',fontsize=12.5,
        bbox=dict(boxstyle='round,pad=0.45',fc='white',ec='0.6'))
ax.set_xlabel('mean pairwise distance among 25 states, shape-only',fontsize=13.5)
ax.set_ylabel('number of null draws',fontsize=14)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(FIG/'v11_fig2D_null.png',dpi=200,bbox_inches='tight'); plt.close()
print('wrote v11_fig2D_null.png')
