"""Fig S23 — the deepening is a genuine shape shift (survives magnitude removal), and is ~constant across baseline.
(A) raw vs shape-only within-patient deepening across 3 cohorts (E17/E17b). (B) additive deepening vs baseline (3.3)."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
from pathlib import Path
ROOT=Path('.').resolve()
coh=['colorectal\n(Pelka)','gastric\n(held-out)','lung\n(corrected)']
raw=[15.5,5.5,15.7]; shp=[10.5,3.4,12.1]; retain=[68,63,77]
fig,(axA,axB)=plt.subplots(1,2,figsize=F.panel_size(1.0, 2.558),gridspec_kw={'width_ratios':[1.1,1]})
x=np.arange(3); w=0.36
axA.bar(x-w/2,raw,w,color='#c0392b',label='raw core')
axA.bar(x+w/2,shp,w,color='#2166ac',label='shape-only core\n(magnitude removed)')
for i in range(3): axA.text(i+w/2,shp[i]+0.3,f'{retain[i]}%',ha='center',fontsize=BASE*0.858,color='#2166ac')
axA.axhline(0,color='k',lw=0.8); axA.set_xticks(x); axA.set_xticklabels(coh,fontsize=BASE*0.921)
axA.set_ylabel('within-patient deepening (points)')
axA.set_title('A',fontsize=BASE*0.95,loc='left',weight='bold')
axA.legend(fontsize=BASE*0.858,frameon=False,loc='upper right')
axA.annotate('% = shape-only\nretained',xy=(1.18,3.4),xytext=(0.30,13.0),fontsize=BASE*0.825,color='#2166ac',
             arrowprops=dict(arrowstyle='->',color='#2166ac',lw=0.8))
# B: additive deepening vs baseline (3.3)
base=[3.2,5.3,8.8,10.9,23.0]; add=[15.7,6.2,5.4,8.5,16.0]; lab=['colon','kidney','gastric','liver','lung']
axB.scatter(base,add,s=70,color='#8e44ad',zorder=3)
for i,l in enumerate(lab): axB.annotate(l,(base[i],add[i]),textcoords='offset points',xytext=(6,4),fontsize=BASE*0.858)
z=np.polyfit(base,add,1); xs=np.linspace(2,25,20); axB.plot(xs,np.polyval(z,xs),'--',color='#8e44ad',lw=1.2)
axB.axhline(np.mean(add),color='0.6',ls=':',lw=1); axB.text(24,np.mean(add)+0.4,f'mean {np.mean(add):.0f}pp',fontsize=BASE*0.858,ha='right',color='0.4')
axB.set_xlabel('normal-tissue baseline (%)'); axB.set_ylabel('deepening (points)')
axB.set_ylim(0,20); axB.grid(alpha=0.25)
axB.set_title('B',fontsize=BASE*0.95,loc='left',weight='bold')
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS23_shapeonly.png'
plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
