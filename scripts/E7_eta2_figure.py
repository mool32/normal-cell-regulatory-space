"""Standalone depiction of eta^2: how much the independent biological family explains cell position on each PC,
versus a permutation null (shuffled family labels). Real eta^2 (0.31-0.70) towers over the chance level (~0.09) —
the spread along every axis is structured by biology, not uniform. Also annotates each axis' variance %."""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
sys.path.insert(0,'tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
PCall=(Rm.values-mu)/sd@load; cells=list(Rm.index)
fam=json.load(open(ROOT+'/results/fam.json'))
POST={'erythrocytes','platelets','reticulocytes'}
keep=[i for i in range(len(cells)) if cells[i] not in POST]
PC=PCall[keep]; famarr=np.array([fam[cells[i]] for i in keep]); ufam=np.unique(famarr)
def eta2(x,g):
    gm=x.mean(); ssb=sum((g==u).sum()*(x[g==u].mean()-gm)**2 for u in ufam); return ssb/((x-gm)**2).sum()
real=[eta2(PC[:,k],famarr) for k in range(6)]
varpct=[PCall[:,k].var()/PCall.var(0).sum()*100 for k in range(6)]
rng=np.random.default_rng(0); NPERM=2000
null=np.array([[eta2(PC[:,k],rng.permutation(famarr)) for k in range(6)] for _ in range(NPERM)])  # 2000 x 6
pval=[(null[:,k]>=real[k]).mean() for k in range(6)]
nlo,nhi=np.percentile(null,2.5,0),np.percentile(null,97.5,0); nmean=null.mean(0)

fig,ax=plt.subplots(figsize=F.panel_size(1.0, 2.60)); xs=np.arange(6)
cols=plt.cm.viridis(np.linspace(0.15,0.85,6))
ax.bar(xs,real,width=0.62,color=cols,edgecolor='white',zorder=3)
# permutation-null band (chance level)
ax.fill_between([-0.6,5.6],[nlo.mean()]*2,[nhi.mean()]*2,color='#bbbbbb',alpha=0.45,zorder=1,label='permutation null (shuffled families, 95%)')
ax.axhline(nmean.mean(),color='#777',lw=1,ls='--',zorder=2)
for k in xs:
    ax.text(k,real[k]+0.015,f'η²={real[k]:.2f}',ha='center',fontsize=BASE*0.968,weight='bold')
    ax.text(k,real[k]/2,f'{varpct[k]:.0f}%\nvar',ha='center',va='center',fontsize=BASE*0.899,color='white',weight='bold')
    ax.text(k,-0.045,'p<0.001' if pval[k]<1e-3 else f'p={pval[k]:.3f}',ha='center',fontsize=BASE*0.825,color='#444')
ax.set_xticks(xs); ax.set_xticklabels([f'PC{k+1}' for k in xs],fontsize=BASE*1.067)
ax.set_ylabel('η²  — fraction of cell position\nexplained by biological family',fontsize=BASE*1.036)
ax.set_ylim(-0.07,0.80); ax.set_xlim(-0.6,5.6)
ax.text(5.45,nhi.mean()+0.01,f'chance ≈ {nmean.mean():.2f}',ha='right',fontsize=BASE*0.899,color='#555',style='italic')
ax.legend(fontsize=BASE*0.906,loc='upper right',framealpha=.95)
# the long two-line ax title duplicated the caption and set the tight bbox (and so the printed type size)
for s in ['top','right']: ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(FIG+'/v7_fig_eta2.png',dpi=190,bbox_inches='tight'); plt.close()
print('saved v7_fig_eta2.png'); print('real eta2:',[round(r,2) for r in real],'| chance~',round(nmean.mean(),3),'| max null',round(null.max(),3))
