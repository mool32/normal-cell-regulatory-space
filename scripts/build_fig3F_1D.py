"""Two new panels requested in the figure plan.
Fig. 3F — enrichment over chance at three successively broader definitions of the convergence region.
Fig. 1D — nearest-neighbour biological-family concordance, observed against the two relevant nulls.
House style for v11: no title above the axes, large fonts, panel letter drawn by the figure assembler."""
import sys, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix()); import figstyle as F
from pathlib import Path
ROOT=Path('.').resolve(); FIG=ROOT/'manuscript/figures'

# ---------------- Fig 3F ----------------
B=pd.read_csv(ROOT/'results/cancer_v6_6_extras/table1_boundaries.csv')
BP=F.use(); FS3F=F.panel_size(0.40,1.25)
fig,ax=plt.subplots(figsize=FS3F)
x=np.arange(len(B)); col=['#8e1b1b','#c0392b','#e8846b']
bars=ax.bar(x,B.fold,color=col,edgecolor='black',lw=1.0,width=.62,zorder=3)
for i,r in B.iterrows():
    ax.text(i,r.fold+0.55,f"{r.fold:.1f}\u00d7",ha='center',fontsize=BP*0.80,fontweight='bold')
    # occupancy inside the bar on two lines: three lines overflowed the narrowest bar
    ax.text(i,r.fold/2,f"{r.occupancy_pct:.1f}%\nof cells",ha='center',va='center',
            fontsize=BP*0.62,color='white',fontweight='bold',linespacing=1.15)
ax.axhline(1,color='0.35',ls='--',lw=1.3,zorder=2)
ax.text(len(B)-0.45,1.5,'chance',fontsize=BP*0.58,color='0.35',ha='right')
ax.set_xticks(x); ax.set_xticklabels([f"top-{int(k)}\n({c:.1f}% by chance)" for k,c in zip(B.k,B.chance_pct)],fontsize=BP*0.62)
ax.set_xlabel('definition of the convergence region\n(rank in the cross-cancer recurrence ranking)',fontsize=BP*0.66)
ax.set_ylabel('enrichment over chance (fold)',fontsize=BP*0.70)
ax.set_ylim(0,max(B.fold)*1.18); ax.grid(axis='y',alpha=.25,zorder=0)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(FIG/'v11_fig3F_boundaries.png',bbox_inches='tight'); plt.close()
F.report('Fig 3F boundaries',FS3F,0.40)
print('wrote v11_fig3F_boundaries.png', list(zip(B.k,B.occupancy_pct,B.fold)))

# ---------------- Fig 1D ----------------
obs, ci = 62.0, (55.0, 70.0)
perm_mean, perm_max = 14.0, 28.0
matched_med, matched_p95 = 61.2, 67.1
hvg = 68.4
FS1D=F.panel_size(0.37,1.35)
fig,ax=plt.subplots(figsize=FS1D)
ypos=[2.6,1.6,0.6]
ax.barh(ypos[0],obs,height=.46,color='#2166ac',edgecolor='black',lw=1.0,zorder=3)
ax.plot([ci[0],ci[1]],[ypos[0]]*2,color='black',lw=2.2,zorder=4)
for c in ci: ax.plot([c,c],[ypos[0]-.1,ypos[0]+.1],color='black',lw=2.2,zorder=4)
ax.barh(ypos[1],matched_med,height=.46,color='#95a5a6',edgecolor='black',lw=1.0,zorder=3)
ax.plot([matched_med,matched_p95],[ypos[1]]*2,color='0.25',lw=2.2,zorder=4)
ax.plot([matched_p95]*2,[ypos[1]-.1,ypos[1]+.1],color='0.25',lw=2.2,zorder=4)
ax.barh(ypos[2],perm_mean,height=.46,color='#d8d8d8',edgecolor='black',lw=1.0,zorder=3)
ax.plot([perm_mean,perm_max],[ypos[2]]*2,color='0.45',lw=2.2,zorder=4)
ax.plot([perm_max]*2,[ypos[2]-.1,ypos[2]+.1],color='0.45',lw=2.2,zorder=4)
ax.axvline(hvg,color='#7b241c',ls=':',lw=2.0,zorder=2)
ax.text(hvg+0.8,3.08,'generic HVG\ndecomposition',fontsize=BP*0.56,color='#7b241c',va='top')
for y,v in zip(ypos,[obs,matched_med,perm_mean]):
    ax.text(2.5,y,f"{v:.0f}%",ha='left',va='center',fontsize=15,fontweight='bold',
            color='white' if v>20 else 'black',zorder=5)
ax.set_yticks(ypos)
ax.set_yticklabels(['the 44 modules\n(95% CI)','size- and expression-\nmatched random panels\n(median, 95th pct)',
                    'family labels permuted\n(mean, max)'],fontsize=BP*0.62)
ax.set_xlabel("nearest neighbour shares the\ncell's biological family (%)",fontsize=BP*0.68)
ax.set_xlim(0,80); ax.set_ylim(0.1,3.3); ax.grid(axis='x',alpha=.25,zorder=0)
for sp in ('top','right','left'): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(FIG/'v11_fig1D_nnfamily.png',bbox_inches='tight'); plt.close()
F.report('Fig 1C nn-family',FS1D,0.37)
print('wrote v11_fig1D_nnfamily.png')
