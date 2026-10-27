"""Fig S27 — the pre-invasive prediction across three organs: holds in colorectum and breast, fails in esophagus."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'
fig,(axA,axB)=plt.subplots(1,2,figsize=(11,4.3),gridspec_kw={'width_ratios':[1,1.1]})
# A: breast ladder
br=pd.read_csv(DD/'E22_dcis.csv')
order=['normal','DCIS','IDC','LN-met']; cols=['#7f8c8d','#e67e22','#c0392b','#8e44ad']
for i,(st,c) in enumerate(zip(order,cols)):
    g=br[br.stage==st]
    if not len(g): continue
    axA.scatter(np.full(len(g),i)+np.random.RandomState(0).uniform(-.09,.09,len(g)),g.core,color=c,s=42,alpha=.8,zorder=3)
    axA.plot([i-.22,i+.22],[g.core.mean()]*2,color=c,lw=3,zorder=4)
    axA.text(i,g.core.mean()+2.4,f'{g.core.mean():.1f}%',ha='center',fontsize=9,color=c,fontweight='bold',
             zorder=8,bbox=dict(boxstyle='round,pad=0.14',fc='white',ec='none',alpha=0.85))
axA.set_xticks(range(4)); axA.set_xticklabels(['normal\n(n=1)','DCIS\n(pre-invasive, n=7)','IDC\n(invasive, n=6)','LN-met\n(n=6)'],fontsize=8.5)
axA.set_ylabel('epithelial core fraction (%)'); axA.axhline(2,color='k',ls=':',lw=1)
axA.set_title('A',fontsize=10.5,loc='left')
# B: esophagus within-patient
be=pd.read_csv(DD/'E22_barrett.csv'); agg=be.groupby(['patient','tissue']).core.mean().unstack()
pairs=agg[['be','es']].dropna()
for _,r in pairs.iterrows(): axB.plot([0,1],[r['es'],r['be']],'-o',color='#7f8c8d',lw=1.4,ms=5,alpha=.75)
axB.plot([0,1],[pairs['es'].mean(),pairs['be'].mean()],'-s',color='#c0392b',lw=3,ms=10,label='mean',zorder=5)
axB.text(-0.28,pairs['es'].mean(),f"{pairs['es'].mean():.1f}%",fontsize=9,va='center',fontweight='bold')
axB.text(1.06,pairs['be'].mean(),f"{pairs['be'].mean():.1f}%",fontsize=9,va='center',fontweight='bold')
axB.set_xticks([0,1]); axB.set_xticklabels(["patient's own\nnormal squamous","Barrett's metaplasia\n(pre-invasive?)"],fontsize=8.5)
axB.set_ylabel('epithelial core fraction (%)'); axB.set_xlim(-0.45,1.4); axB.axhline(2,color='k',ls=':',lw=1)
axB.set_title('B',fontsize=10.5,loc='left')
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS27_preinvasive2.png'; plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
