"""Fig S19 — specification curve of the core convergence (referee round-2 ask #6, multiverse).
Reads results/E7_spec_curve.csv (36 specs). Two panels:
 A) sorted specification curve: core fraction per spec, coloured by calibration, each spec's own
    chance baseline as a tick; the single-cell specs all clear chance, the (independently degenerate)
    pseudobulk specs collapse.
 B) median core fraction by analytic choice, split by calibration — shows which choices move the estimate."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; df=pd.read_csv(ROOT/'results/E7_spec_curve.csv')
SC='#2166ac'; PB='#b2182b'
fig,(axA,axB)=plt.subplots(1,2,figsize=F.panel_size(1.00, 1.75),gridspec_kw={'width_ratios':[1.7,1]})

# --- A: sorted specification curve ---
d=df.sort_values('core_pct',ascending=False).reset_index(drop=True)
x=np.arange(len(d)); col=[SC if c=='single_cell_scaled' else PB for c in d.calib]
axA.bar(x,d.core_pct,color=col,width=0.8,zorder=3)
axA.scatter(x,d.chance,marker='_',s=90,color='k',lw=1.4,zorder=4,label='per-spec chance baseline')
axA.set_ylabel('core convergence (%)'); axA.set_xlabel('analytic specification (sorted)')
axA.set_title('A',fontsize=BASE*1.038,loc='left')
axA.set_xlim(-0.7,len(d)-0.3)
from matplotlib.patches import Patch
axA.legend(loc='upper center',bbox_to_anchor=(0.5,-0.19),ncol=2,handles=[Patch(color=SC,label='single-cell calibration'),
                    Patch(color=PB,label='pseudobulk calibration (degenerate)'),
                    plt.Line2D([0],[0],marker='_',color='k',lw=1.4,ls='none',label='per-spec chance baseline')],
           fontsize=BASE*0.825,frameon=False)
sc=d[d.calib=='single_cell_scaled']; n_below=int((d.core_pct<=d.chance).sum())
axA.text(0.985,0.99,f'single-cell (blue): above chance in\nall 18 specs — median {sc.core_pct.median():.1f}%, ≥{sc.enrich.min():.1f}× chance',
         transform=axA.transAxes,fontsize=BASE*0.889,color=SC,va='top',ha='right')
axA.annotate(f'the {n_below} specs that collapse to chance\nare ALL pseudobulk (degenerate\ncalibration, Note S-B)',
             xy=(len(d)-3,0.8),xytext=(0.42,0.55),textcoords=axA.transAxes,fontsize=BASE*0.858,color=PB,va='top',
             arrowprops=dict(arrowstyle='->',color=PB,lw=1))

# --- B: median core% by analytic choice, split by calibration ---
choices=[('magnitude',['raw','center','rank']),('coredef',['3cell','2cell','cyto']),('eryth',['excl','incl'])]
labels=[]; scv=[]; pbv=[]
for col_,levels in choices:
    for lv in levels:
        labels.append(f'{col_}={lv}')
        scv.append(df[(df.calib=='single_cell_scaled')&(df[col_]==lv)].core_pct.median())
        pbv.append(df[(df.calib=='pseudobulk')&(df[col_]==lv)].core_pct.median())
y=np.arange(len(labels)); h=0.38
axB.barh(y+h/2,scv,height=h,color=SC,label='single-cell'); axB.barh(y-h/2,pbv,height=h,color=PB,label='pseudobulk')
axB.set_yticks(y); axB.set_yticklabels(labels,fontsize=BASE*0.858); axB.invert_yaxis()
axB.set_xlabel('median core convergence (%)'); axB.set_title('B',fontsize=BASE*1.038,loc='left')
axB.axvline(2.0,color='k',ls=':',lw=1); axB.legend(fontsize=BASE*0.858,loc='lower right',frameon=False)
for i,l in enumerate(labels):
    if l.split('=')[0]!=(labels[i-1].split('=')[0] if i else None): axB.axhline(i-0.5,color='0.85',lw=0.6)

plt.tight_layout(); out=ROOT/'manuscript/figures/v7_figS19_speccurve.png'
out.parent.mkdir(exist_ok=True); plt.savefig(out,dpi=200,bbox_inches='tight'); print('wrote',out)
