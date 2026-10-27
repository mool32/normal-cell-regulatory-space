"""Fig S9b (full-fidelity) — robustness of the convergence to the exact 44-module choice (Gudkov G6b).
Full 25-cancer canonical pipeline (pct.score_modules + rebuild() eigenspace); reads results/E2e_fullfid_*.csv.
A: drop 5 random modules x200 -> core-fraction magnitude stable.
B: the cytotrophoblast ANCHOR (rank by cosine to the malignant attractor) stays at the top under the same drops.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
d5=pd.read_csv(ROOT+'/results/E2e_fullfid_leave5.csv')
d1=pd.read_csv(ROOT+'/results/E2e_fullfid_leave1.csv')
b=pd.read_csv(ROOT+'/results/E2e_fullfid_baseline.csv').iloc[0]
REAL=b.corefrac*100; REAL_CYTO=int(b.cyto_attr)
AC='#c0392b'; DK='#7b241c'; INK='#222'; MUT='#666'; GREY='#c8ccd1'

fig=plt.figure(figsize=F.panel_size(1.0, 2.739)); gs=fig.add_gridspec(1,2,width_ratios=[1.05,1.0],wspace=0.26)

# ---- A: leave-5-out core-fraction distribution ----
axA=fig.add_subplot(gs[0,0]); cf=d5.corefrac.values*100
axA.hist(cf,bins=22,color='#dfe6ec',edgecolor='white',lw=.6,zorder=2)
axA.axvline(REAL,color=AC,lw=2.0,zorder=4,label=f'44-module result ({REAL:.1f}%)')
axA.axvline(cf.min(),color=MUT,lw=1.1,ls='--',zorder=3,label=f'worst of 200 draws ({cf.min():.1f}%)')
axA.axvline(100*3/154,color='#888',lw=1.0,ls=':',zorder=3,label='chance (1.9%)')
axA.set_ylim(0,axA.get_ylim()[1]*1.28)
axA.set_xlabel('core fraction (% of malignant cells)\n→ cytotrophoblast / megakaryocyte',fontsize=BASE*0.95)
axA.set_ylabel('draws',fontsize=BASE*0.976)
axA.set_title('A',fontsize=BASE*1.002,loc='left',weight='bold')
axA.text(0.985,0.97,f'mean {cf.mean():.1f} ± {cf.std():.1f}%\nnever < {cf.min():.1f}%  (≥ {d5.enrich.min():.1f}× chance)\nleave-one-out {d1.corefrac.min()*100:.1f}–{d1.corefrac.max()*100:.1f}%',
         transform=axA.transAxes,fontsize=BASE*0.839,color=MUT,ha='right',va='top')
axA.legend(fontsize=BASE*0.825,loc='upper left',framealpha=.95)
for s in ['top','right']: axA.spines[s].set_visible(False)

# ---- B: cytotrophoblast anchor rank (cosine-to-attractor) under leave-5-out ----
axB=fig.add_subplot(gs[0,1])
r=d5.cyto_attr.values; N=len(r)
cats=['1','2','3','4','5','≥6']
counts=[ (r==1).sum(),(r==2).sum(),(r==3).sum(),(r==4).sum(),(r==5).sum(),(r>=6).sum() ]
pct=[100*c/N for c in counts]
x=np.arange(len(cats))
cols=[AC if i<3 else GREY for i in range(len(cats))]
axB.bar(x,pct,color=cols,edgecolor='white',lw=.5,width=0.78,zorder=3)
for xi,p in zip(x,pct):
    if p>1: axB.text(xi,p+0.6,f'{p:.0f}%',ha='center',va='bottom',fontsize=BASE*0.868,color=INK)
axB.set_xticks(x); axB.set_xticklabels(cats,fontsize=BASE*0.943)
axB.set_xlabel('cytotrophoblast rank among 154 cells\n(cosine to the malignant attractor)',fontsize=BASE*0.93)
axB.set_ylabel('% of 200 draws',fontsize=BASE*0.976)
axB.set_ylim(0,max(pct)*1.22)
axB.set_title('B',fontsize=BASE*1.002,loc='left',weight='bold')
top3=100*(r<=3).mean(); top5=100*(r<=5).mean()
axB.text(0.985,0.95,f'unperturbed: rank {REAL_CYTO}\nin top-3: {top3:.0f}% of draws\nin top-5: {top5:.0f}%\nleave-one-out: top-3 in 100% of 44',transform=axB.transAxes,fontsize=BASE*0.854,color=MUT,ha='right',va='top')
for s in ['top','right']: axB.spines[s].set_visible(False)
axB.axvspan(-0.5,2.5,color='#fbeee9',zorder=0)

# suptitle removed (house rule A5)
plt.tight_layout()
plt.savefig(FIG+'/v7_figS9b_module_robustness.png',dpi=190,bbox_inches='tight'); plt.close()
print(f'saved | core {cf.mean():.1f}±{cf.std():.1f}% (min {cf.min():.1f}, enrich min {d5.enrich.min():.1f}x) | cyto top-3 {top3:.0f}% top-5 {top5:.0f}% (real rank {REAL_CYTO})')
