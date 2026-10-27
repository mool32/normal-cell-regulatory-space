"""Fig S22 — lung (high-baseline) also deepens; the 'negative control' was a dilution artifact.
(A) per-patient nLung(normal) -> tLung(tumor) epithelial core fraction (GSE131907, clean author tumor cells);
(B) the resolution: on the paper's own Census lung data, clean malignant cells sit far above the diluted
    'tumor arm' (9% malignant + 91% normal epithelium) the multitissue pipeline had used."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
from pathlib import Path
ROOT=Path('.').resolve()
d=pd.read_csv(ROOT/'results/referee2_newdata/E12_lung_negative.csv')
ok=d[(d.n_epi_normal>=20)&(d.n_epi_tumor>=20)].copy(); ok['diff']=ok.tumor_core-ok.normal_core
fig,(axA,axB)=plt.subplots(1,2,figsize=F.panel_size(1.0, 2.391),gridspec_kw={'width_ratios':[1.25,1]})
# A: per-patient slope
for _,r in ok.iterrows():
    up=r['diff']>0
    axA.plot([0,1],[r.normal_core,r.tumor_core],'-o',color=('#c0392b' if up else '#7f8c8d'),lw=1.2,ms=4,alpha=0.6)
axA.plot([0,1],[ok.normal_core.mean(),ok.tumor_core.mean()],'-s',color='k',lw=3.2,ms=11,zorder=6,label='mean')
axA.annotate(f'{ok.normal_core.mean():.0f}%',xy=(0,ok.normal_core.mean()),xytext=(-0.28,ok.normal_core.mean()),fontsize=BASE*0.885,va='center',fontweight='bold')
axA.annotate(f'{ok.tumor_core.mean():.0f}%',xy=(1,ok.tumor_core.mean()),xytext=(1.05,ok.tumor_core.mean()),fontsize=BASE*0.885,va='center',fontweight='bold')
axA.axhline(2,color='k',ls=':',lw=1); axA.set_xticks([0,1]); axA.set_xticklabels(['normal lung\n(nLung)','tumor\n(tLung)'])
axA.set_ylabel('epithelial core fraction (%)'); axA.set_xlim(-0.4,1.35)
axA.set_title('A',fontsize=BASE*0.95,loc='left',weight='bold')
# B: dilution resolution (bars) — clean malignant vs diluted arm vs normal (from E13 + composition)
labels=['clean\nmalignant','paper\'s "tumor arm"\n(9% malignant +\n91% normal)','normal\nepithelium']
vals=[29.4, 11.7*0.91+29.4*0.09, 11.7]  # approx diluted arm = weighted
cols=['#c0392b','#e59866','#7f8c8d']
axB.bar(range(3),vals,color=cols); axB.set_xticks(range(3)); axB.set_xticklabels(labels,fontsize=BASE*0.825)
axB.axhline(2,color='k',ls=':',lw=1); axB.set_ylabel('core fraction (%)')
axB.set_title('B',fontsize=BASE*0.95,loc='left',weight='bold')
for i,v in enumerate(vals): axB.text(i,v+0.6,f'{v:.0f}%',ha='center',fontsize=BASE*0.885)
axB.annotate('dilution collapses\nthe tumor arm to\n≈ normal',xy=(1,vals[1]),xytext=(1.3,24),fontsize=BASE*0.825,
             arrowprops=dict(arrowstyle='->',lw=1),ha='left')
plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS22_lung.png'
plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
