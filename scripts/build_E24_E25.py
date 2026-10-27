"""Fig S24 — malignancy-specificity: fibrosis (non-malignant) does NOT enter the state, cancer does.
Fig S25 — the state is relinquished on colonisation: within-patient primary vs lymph-node metastasis."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
from pathlib import Path
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'

# ---------------- Fig S24: IPF vs lung cancer ----------------
ipf=pd.read_csv(DD/'E24_ipf_specificity.csv')
fig,(axA,axB)=plt.subplots(1,2,figsize=F.panel_size(1.0, 2.558),gridspec_kw={'width_ratios':[1,1.15]})
labels=['fibrosis (IPF)\nnon-malignant','lung cancer\nmalignant']
raw=[3.8, 15.7]; shp=[8.9-7.6, 30.3-18.2]
x=np.arange(2); w=0.36
axA.bar(x-w/2,raw,w,color='#c0392b',label='raw core'); axA.bar(x+w/2,shp,w,color='#2166ac',label='shape-only core')
for i,(r,s) in enumerate(zip(raw,shp)):
    # the value sits inside its bar and the count/P above it, so the two labels no longer overlap
    axA.text(i-w/2,r/2,f'+{r:.1f}',ha='center',va='center',color='white',fontweight='bold',fontsize=BASE*0.86)
    axA.text(i+w/2,s/2,f'+{s:.1f}',ha='center',va='center',color='white',fontweight='bold',fontsize=BASE*0.86)
axA.axhline(0,color='k',lw=0.8); axA.set_xticks(x); axA.set_xticklabels(labels,fontsize=BASE*0.921)
axA.set_ylabel('deepening vs the organ\'s\nown control (pts)')  # one line ran off the top of the canvas and was clipped to "(pt"
axA.set_title('A',fontsize=BASE*1.01,loc='left')
axA.legend(fontsize=BASE*0.858,frameon=False); axA.set_ylim(0,21)
axA.text(0,5.2,'n.s.\n(P=0.12 / 0.21)',ha='center',fontsize=BASE*0.858,color='0.35')
# raw and shape-only differ in BOTH the count and the P value; one bracket over the pair read as
# a single '10/10' and contradicted Fig. S23 (review finding 14)
axA.text(1-w/2,raw[1]+0.5,'10/10 patients\nP = 0.002',ha='center',va='bottom',fontsize=BASE*0.858,color='#c0392b')
axA.text(1+w/2,shp[1]+0.5,'9/10 patients\nP = 0.004',ha='center',va='bottom',fontsize=BASE*0.858,color='#2166ac')
# B: cell types within fibrotic lung
ct=[('Transitional AT2',18.4),('Proliferating epith.',16.2),('AT2',15.5),('Basal',13.4),
    ('SCGB3A2+',12.7),('KRT5-/KRT17+',8.3),('MUC5B+',7.9),('Ciliated',5.2)]
names=[c[0] for c in ct][::-1]; vals=[c[1] for c in ct][::-1]
cols=['#c0392b' if v>=13 else '#95a5a6' for v in vals]
axB.barh(range(len(vals)),vals,color=cols)
axB.set_yticks(range(len(vals))); axB.set_yticklabels(names,fontsize=BASE*0.89)
axB.set_xlabel('shape-only core fraction (%)')
axB.axvline(30.3,color='#c0392b',ls='--',lw=1.2); axB.text(30.0,0.2,'malignant lung 30.3%',fontsize=BASE*0.825,color='#c0392b',ha='right',rotation=90)
axB.set_title('B',fontsize=BASE*0.981,loc='left')
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS24_fibrosis.png'; plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)

# ---------------- Fig S25: primary vs metastasis ----------------
d=pd.read_csv(DD/'E25_primary_vs_met.csv')
b=d[d.cohort=='breast']; piv=b.pivot_table(index='patient',columns='site',values=['RAW','SHAPE'])
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 2.60))
for m,c,off in [('RAW','#c0392b',0),('SHAPE','#2166ac',0.06)]:
    ok=piv[m][['primary','metastasis']].dropna()
    for _,r in ok.iterrows():
        ax.plot([0+off,1+off],[r['primary'],r['metastasis']],'-o',color=c,lw=1.1,ms=4,alpha=0.55)
    ax.plot([0+off,1+off],[ok['primary'].mean(),ok['metastasis'].mean()],'-s',color=c,lw=3,ms=10,
            label=f"{'raw' if m=='RAW' else 'shape-only'} core (mean)")
ax.set_xticks([0.03,1.03]); ax.set_xticklabels(['primary tumour','lymph-node\nmetastasis'])
ax.set_ylabel('core fraction (%)'); ax.set_xlim(-0.35,1.45)
ax.set_title('',fontsize=BASE*0.981,loc='left')
ax.legend(fontsize=BASE*0.858,frameon=False); ax.axhline(2,color='k',ls=':',lw=1)
plt.tight_layout(); out2=ROOT/'manuscript/figures/v7_figS25_metastasis.png'; plt.savefig(out2,dpi=200,bbox_inches='tight'); print('wrote',out2)
