"""SUPERSEDED. This script and scripts/build_E2e_robustness.py both wrote
v7_figS9b_module_robustness.png, so whichever ran last won. The Fig. S9b caption describes the
FULL-FIDELITY run (E2e_fullfid_*: 17.5 +/- 2.5%, never below 11.7%), but the shipped PNG carried
this script's earlier numbers (18.1 +/- 3.0%, 10.8%). E2e is now the sole writer; this one saves
under a _SUPERSEDED name and is kept only for provenance.
"""
"""Fig S9b — robustness of the convergence to the exact 44-module choice (Gudkov G6b).
Reads results/E2d_robustness_{leave5,leave1,adddecoy}.csv (from scripts/E2d_module_robustness.py).
A: drop 5 random modules x200 -> distribution of core fraction (magnitude stable).
B: drop each single module (44 rebuilds) -> none is load-bearing.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
d5=pd.read_csv(ROOT+'/results/E2d_robustness_leave5.csv')
d1=pd.read_csv(ROOT+'/results/E2d_robustness_leave1.csv').sort_values('corefrac').reset_index(drop=True)
d3=pd.read_csv(ROOT+'/results/E2d_robustness_adddecoy.csv')
REAL=18.5; REAL_COMPACT=5.63   # from the unperturbed 44-module run (E2d_module_robustness.py)
AC='#c0392b'; DK='#7b241c'; INK='#222'; MUT='#666'; GREY='#c8ccd1'

fig=plt.figure(figsize=F.panel_size(1.0, 2.739)); gs=fig.add_gridspec(1,2,width_ratios=[1.05,1.0],wspace=0.24)

# ---- A: leave-5-out distribution ----
axA=fig.add_subplot(gs[0,0])
cf=d5.corefrac.values*100
axA.hist(cf,bins=22,color='#dfe6ec',edgecolor='white',lw=.6,zorder=2)
axA.axvline(REAL,color=AC,lw=2.0,zorder=4,label=f'44-module result ({REAL:.1f}%)')
axA.axvline(cf.min(),color=MUT,lw=1.1,ls='--',zorder=3,label=f'worst of 200 draws ({cf.min():.1f}%)')
axA.axvline(100*3/154,color='#888',lw=1.0,ls=':',zorder=3,label='chance (1.9%)')
axA.set_xlabel('core fraction (%)',fontsize=BASE*0.967)
axA.set_ylabel('draws',fontsize=BASE*0.994)
axA.set_title('A',fontsize=BASE*1.02,loc='left',weight='bold')
axA.set_ylim(0, axA.get_ylim()[1]*1.28)
axA.text(0.985,0.60,f'mean {cf.mean():.1f} ± {cf.std():.1f}%\nnever < {cf.min():.1f}%  (≥ {d5.enrich.min():.1f}× chance)\ncompactness {d5.ratio.mean():.1f} (real {REAL_COMPACT:.1f})\nanchor a top-5 target in 69%',
         transform=axA.transAxes,fontsize=BASE*0.855,color=MUT,ha='right',va='top')
axA.legend(fontsize=BASE*0.84,loc='upper left',framealpha=.95)
for s in ['top','right']: axA.spines[s].set_visible(False)

# ---- B: leave-1-out (each single module dropped) ----
axB=fig.add_subplot(gs[0,1])
y=np.arange(len(d1)); cf1=d1.corefrac.values*100
col=[AC if v<REAL-2.5 else GREY for v in cf1]
axB.barh(y,cf1,color=col,edgecolor='white',lw=.4,height=0.82,zorder=3)
axB.axvline(REAL,color=AC,lw=1.6,ls='-',zorder=4)
axB.text(REAL,len(d1)+0.6,f'44-module result {REAL:.1f}%',fontsize=BASE*0.84,color=AC,ha='center',va='bottom')
# label only the 3 lowest and 2 highest
for i,dy in [(0,0.0),(len(d1)-2,-0.9),(len(d1)-1,0.9)]:   # the top two labels overlapped
    axB.text(cf1[i]+0.3,y[i]+dy,d1.dropped.iloc[i],va='center',ha='left',fontsize=BASE*0.825,
             color=(AC if i==0 else INK))
axB.set_yticks([]); axB.set_xlim(0,max(cf1)+4); axB.set_ylim(-1,len(d1)+2)
axB.set_xlabel('core fraction after the drop (%)',fontsize=BASE*0.967)
axB.set_title('B',fontsize=BASE*1.02,loc='left',weight='bold')
axB.text(0.985,0.42,f'range {cf1.min():.1f}–{cf1.max():.1f}%\n(real {REAL:.1f}%)\nanchor top-5 in\n95% of 44 drops',transform=axB.transAxes,fontsize=BASE*0.855,color=MUT,ha='right',va='center')
for s in ['top','right','left']: axB.spines[s].set_visible(False)

# suptitle removed: it duplicated the figure caption and no longer fits the printed width
plt.tight_layout(rect=(0,0,1,0.94))
plt.savefig(FIG+'/v7_figS9b_module_robustness_E2d_SUPERSEDED.png',dpi=190,bbox_inches='tight'); plt.close()
print(f'saved v7_figS9b_module_robustness.png | leave5 core {cf.mean():.1f}±{cf.std():.1f}% (min {cf.min():.1f}) | leave1 {cf1.min():.1f}-{cf1.max():.1f}% | add-decoy {d3.corefrac.mean()*100:.1f}%')
