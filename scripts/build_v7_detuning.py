"""Fig 9: within-patient difference-in-occupancy (Pelka colorectal).
Both malignant and matched-normal cells are IN the duet state; we ask how the MALIGNANT occupant
differs from the NORMAL occupant of the same state, n_genes-controlled (within-patient residualization
on log10 n_genes; reproduced by n_genes-matched downsampling).
A = controlled per-module delta (malignant - normal) for modules that survive the detection control,
    with the raw (uncontrolled) value shown as a ghost so the shrinkage from the n_genes control is visible.
B = the delta vs the duet's released-brake axis: mostly off-axis (rho=+0.36) -> a qualitative detuning,
    not 'deeper on the same axis'. IFN/JAK-STAT/BMP do not survive the control (named).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
from matplotlib.lines import Line2D
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
T=pd.read_csv(ROOT+'/results/pelka_diff_occupancy_controlled.csv')

CLASS={**{m:'stress' for m in ['HIF','NRF2','HSF1','UPR-PERK','p53','cGAS-STING','Autophagy']},
       **{m:'growth' for m in ['mTOR','Cell Cycle','Wnt','SREBP','ER','PI3K/PTEN']},
       **{m:'metab'  for m in ['AhR','PPARα','MR']}}
COL={'stress':'#c0392b','growth':'#8e44ad','metab':'#d68910','other':'#c8ccd1'}
LAB={'stress':'integrity / stress sensor','growth':'growth / anabolic','metab':'metabolic / hormonal'}
DISP={'p53':'p53*'}  # flag p53: rise is target-driven (stress footprint), not TP53 / retained function

surv=T[(T.resid_delta.abs()>0.25)&(T.resid_pW<0.05)].copy()
surv=surv.sort_values('resid_delta')
surv['cls']=surv.module.map(lambda m:CLASS.get(m,'other'))

fig=plt.figure(figsize=F.panel_size(1.0, 1.88)); gs=fig.add_gridspec(1,2,width_ratios=[1.0,1.12],wspace=0.26)

# ---- Panel A: controlled delta bars + raw ghost ----
axA=fig.add_subplot(gs[0,0])
y=np.arange(len(surv))
axA.barh(y,surv.resid_delta,color=[COL[c] for c in surv.cls],edgecolor='white',lw=.5,zorder=3,height=0.72)
# raw (uncontrolled) value as a hollow marker to show the n_genes-control shrinkage
axA.scatter(surv.raw_delta,y,facecolor='none',edgecolor='#555',s=26,lw=.9,zorder=4)
for yi,(_,r) in zip(y,surv.iterrows()):
    if abs(r.raw_delta-r.resid_delta)>0.12:
        axA.plot([r.resid_delta,r.raw_delta],[yi,yi],color='#999',lw=.7,ls=':',zorder=2)
axA.axvline(0,color='#888',lw=.8)
axA.set_yticks(y); axA.set_yticklabels([DISP.get(m,m) for m in surv.module],fontsize=BASE*0.966)
for t,c in zip(axA.get_yticklabels(),surv.cls): t.set_color(COL[c]); t.set_weight('bold')
axA.set_xlabel('malignant − normal occupant (module z)',fontsize=BASE*0.994)
axA.set_title('A',fontsize=BASE*1.00,loc='left',weight='bold')
axA.text(0.98,0.12,'○ = raw Δ, before the detection-depth control',transform=axA.transAxes,fontsize=BASE*0.80,color='#555',ha='right',va='bottom')
for s in ['top','right']: axA.spines[s].set_visible(False)
axA.set_ylim(-0.7,len(surv)-0.3)

# ---- Panel B: delta vs released-brake axis ----
axB=fig.add_subplot(gs[0,1])
ns=T[~T.module.isin(surv.module)]
axB.axvspan(T.distinct.min()-0.1,-0.5,color='#eef2f6',zorder=0)
axB.text(T.distinct.min()-0.05,axB.get_ylim()[1] if False else 1.18,'released-brake axis\n(core-defining) →',fontsize=BASE*0.864,color='#5d6d7e',va='top',ha='left',style='italic')
axB.axhline(0,color='#bbb',lw=.7); axB.axvline(0,color='#bbb',lw=.7,ls=':')
axB.scatter(ns.distinct,ns.resid_delta,s=16,color='#c8ccd1',zorder=2)
# the three crowded survivors get explicit data-space label anchors + short leader lines
# hand-placed label anchors did not survive the move to the printed width; the labels are now
# de-collided against their real extents (scripts/figstyle.py declutter)
_t18=[]
for _,r in surv.iterrows():
    c=COL[CLASS.get(r.module,'other')]
    axB.scatter(r.distinct,r.resid_delta,s=58,color=c,edgecolor='white',lw=.5,zorder=4)
    _t18.append(axB.text(r.distinct,r.resid_delta,DISP.get(r.module,r.module),
                         fontsize=BASE*0.871,color=c,weight='bold',zorder=6,ha='center',va='bottom'))
for m in ['Type I IFN','JAK-STAT','BMP']:
    r=T[T.module==m].iloc[0]
    _t18.append(axB.text(r.distinct,r.resid_delta,m,fontsize=BASE*0.825,color='#7f8c8d',
                         style='italic',zorder=6,ha='center',va='bottom'))
import numpy as _np
F.declutter(axB,_t18,avoid=_np.c_[T.distinct.values,T.resid_delta.values])
axB.set_xlabel('core distinctiveness\n(← reduced signature      retained/high →)',fontsize=BASE*0.994)
axB.set_ylabel('malignant − normal occupant  Δ  (controlled)',fontsize=BASE*0.994)
axB.set_title('B',fontsize=BASE*1.00,loc='left',weight='bold')
hd=[Line2D([0],[0],marker='o',color='w',markerfacecolor=COL[k],markersize=8,label=LAB[k]) for k in ['stress','growth','metab']]
axB.legend(handles=hd,fontsize=BASE*0.879,loc='upper left',bbox_to_anchor=(0.0,0.80),framealpha=.95)
for s in ['top','right']: axB.spines[s].set_visible(False)

# suptitle removed: it duplicated the figure caption and no longer fits the printed width
# the two footnotes moved to the Fig. S18 caption: at the printed width they overlapped
plt.tight_layout(rect=(0,0.03,1,0.96))
plt.savefig(FIG+'/v7_fig_detuning.png',dpi=195,bbox_inches='tight'); plt.close()
print('saved v7_fig_detuning.png; survivors:',len(surv),'| up:',list(surv[surv.resid_delta>0].module),'| down:',list(surv[surv.resid_delta<0].module))
