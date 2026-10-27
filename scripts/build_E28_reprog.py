"""Fig S28 — the cancer-blind coordinate detects trophoblast-like cells in mouse reprogramming."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path('.').resolve(); d=pd.read_csv(ROOT/'results/referee2_newdata/E28_reprogramming.csv').sort_values('day')
fig,ax=plt.subplots(figsize=(7.4,4.3))
ax.plot(d.day,d.SHAPE_cyto,'-o',color='#2166ac',lw=2.2,ms=5,label='cytotrophoblast-assigned (shape-only)')
ax.plot(d.day,d.RAW_cyto,'-',color='#c0392b',lw=1.2,alpha=.6,label='cytotrophoblast-assigned (raw)')
ax.axvspan(-0.3,0.8,color='0.85',alpha=.6); ax.text(0.25,25.0,'starting\nMEFs*',fontsize=7.5,ha='center',color='0.4')
ax.axvspan(8,18.4,color='#2166ac',alpha=.07)
ax.text(13,25.0,'late phase: extraembryonic /\ntrophoblast-like fates emerge',fontsize=8.5,ha='center',color='#2166ac')
ax.annotate('trough 1.8%\n(day 4)',xy=(4,1.8),xytext=(2.4,7.5),fontsize=8,color='0.35',arrowprops=dict(arrowstyle='->',lw=1,color='0.45'))
ax.annotate('plateau 16–20%\n(~10× the trough)',xy=(16,18.4),xytext=(9.6,4.5),fontsize=8,color='#2166ac',arrowprops=dict(arrowstyle='->',lw=1,color='#2166ac'))
ax.set_xlabel('day of reprogramming (MEF → iPSC)'); ax.set_ylabel('cells assigned to cytotrophoblast (%)')
ax.set_title('',fontsize=10,loc='left')
ax.legend(fontsize=8,frameon=False,loc='upper left',bbox_to_anchor=(0.11,0.80)); ax.set_ylim(0,27); ax.grid(alpha=.22)
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS28_reprogramming.png'; plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
