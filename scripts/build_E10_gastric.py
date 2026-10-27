"""Fig S21 — held-out gastric within-patient deepening (GSE150290, ask #1).
Per-patient paired normal->cancer epithelial core fraction (slope plot); reads E10_gastric_deepening.csv."""
import sys, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts')); import figstyle as F
BASE=F.use()   # re-authored at HALF the printed width so it can sit in a merged figure at 8 pt
df=pd.read_csv(ROOT/'results/referee2_newdata/E10_gastric_deepening.csv')
ok=df[(df.n_epi_normal>=20)&(df.n_epi_cancer>=20)].copy(); ok['diff']=ok.cancer_core-ok.normal_core
fig,ax=plt.subplots(figsize=F.panel_size(0.5,1.55))
for _,r in ok.iterrows():
    up=r['diff']>0
    ax.plot([0,1],[r.normal_core,r.cancer_core],'-',color=('#c0392b' if up else '#7f8c8d'),lw=1.2,alpha=0.6,zorder=2)
    ax.plot([0,1],[r.normal_core,r.cancer_core],'o',color=('#c0392b' if up else '#7f8c8d'),ms=4,alpha=0.7,zorder=3)
ax.plot([0,1],[ok.normal_core.mean(),ok.cancer_core.mean()],'-',color='k',lw=3.2,zorder=5)
ax.plot([0,1],[ok.normal_core.mean(),ok.cancer_core.mean()],'s',color='k',ms=11,zorder=6,label='mean')
ax.annotate(f'mean {ok.normal_core.mean():.1f}%',xy=(0,ok.normal_core.mean()),xytext=(-0.34,ok.normal_core.mean()),fontsize=BASE*0.82,va='center',fontweight='bold')
ax.annotate(f'mean {ok.cancer_core.mean():.1f}%',xy=(1,ok.cancer_core.mean()),xytext=(1.06,ok.cancer_core.mean()),fontsize=BASE*0.82,va='center',fontweight='bold')
ax.axhline(2.0,color='k',ls=':',lw=1)
ax.text(1.02,2.0,'chance ~2%',fontsize=BASE*0.72,va='bottom')
ax.set_xticks([0,1]); ax.set_xticklabels(['matched normal\n(Pat##-A)','cancer\n(Pat##-B)'])
ax.set_ylabel('epithelial core fraction (%)'); ax.set_xlim(-0.45,1.4)
ax.text(0.5,ax.get_ylim()[1]*0.94,'red = deeper in cancer   grey = not',fontsize=BASE*0.78,ha='center',color='0.3')
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS21_gastric.png'
plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
