"""Fig S20 — pre-invasive prediction test (GSE161277 colon, ask #4).
Per-patient within-patient progression normal -> adenoma -> carcinoma of the epithelial core fraction;
(A) full core (cyto+mega), (B) cytotrophoblast-specific fraction. Reads E9_colon_preinvasive.csv."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/'results/referee2_newdata/E9_colon_preinvasive.csv')
STAGES=['normal','adenoma','carcinoma']
def admean(g,col): return np.nanmean([g[g.stage==s][col].mean() for s in ['adenoma','adenoma_1','adenoma_2']])
def patseries(pat,col):
    g=df[df.patient==pat]; d={}
    for s in STAGES:
        d[s]= admean(g,col) if s=='adenoma' else (g[g.stage==s][col].mean() if s in set(g.stage) else np.nan)
    return [d[s] for s in STAGES]
pats=[p for p in sorted(df.patient.unique()) if not np.isnan(patseries(p,'core')[0]) and not np.isnan(patseries(p,'core')[1]) and not np.isnan(patseries(p,'core')[2])]
fig,(axA,axB)=plt.subplots(1,2,figsize=F.panel_size(1.0, 2.326))
x=np.arange(3); cols=plt.cm.viridis(np.linspace(0.15,0.75,len(pats)))
for ax,metric,title in [(axA,'core','A'),(axB,'cyto','B')]:
    M=[]
    for pat,c in zip(pats,cols):
        y=patseries(pat,metric); M.append(y)
        ax.plot(x,y,'-o',color=c,lw=1.6,ms=6,alpha=0.85,label=pat)
    M=np.array(M); ax.plot(x,np.nanmean(M,0),'-s',color='crimson',lw=3,ms=9,label='mean',zorder=5)
    ax.axhline(2.0,color='k',ls=':',lw=1,label='chance ~2%')
    ax.set_xticks(x); ax.set_xticklabels(['normal','adenoma\n(pre-invasive)','carcinoma'])
    ax.set_ylabel('epithelial core fraction (%)'); ax.set_title(title,fontsize=BASE*0.999,loc='left')
    ax.set_ylim(bottom=-1)
axA.legend(fontsize=BASE*0.825,frameon=False,loc='upper left')
axA.annotate('pre-invasive adenoma\nalready elevated\n(3/3 patients)',xy=(1,11),xytext=(0.15,17),fontsize=BASE*0.855,color='crimson',
             arrowprops=dict(arrowstyle='->',color='crimson',lw=1))
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS20_preinvasive.png'
plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
