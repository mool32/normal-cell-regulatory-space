"""E1 control figure (Gudkov #1): is megakaryocyte a low-complexity artifact target?
A = assignment vs n_genes decile (mega tracks eryth DOWN with complexity; cytotrophoblast OPPOSITE).
B = assignment among cells above an n_genes threshold (cyto strengthens, mega decays toward eryth residual).
C = causal downsampling: thin the SAME high-complexity cells; cyto COLLAPSES, eryth SURGES (positive control),
    mega rises only modestly -> partly complexity-sensitive, far less than erythrocyte, but real.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 0.99 x 6.5 in
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
PC=pd.read_csv(f'{ROOT}/results/E1_percell_sub.csv.gz')
D=pd.read_csv(f'{ROOT}/results/E1b_downsample.csv')
CY='#1e8449'; MG='#e67e22'; ER='#7f8c8d'   # cyto green(real), mega orange(questioned), eryth grey(artifact)

fig=plt.figure(figsize=F.panel_size(0.99, 2.35)); gs=fig.add_gridspec(1,3,wspace=0.34,width_ratios=[1,0.92,1.05])

# A: gradient by decile
axA=fig.add_subplot(gs[0,0])
PC['dec']=pd.qcut(PC.ng,10,labels=False,duplicates='drop')
g=PC.groupby('dec').apply(lambda d:pd.Series({'ng':d.ng.median(),
    'mega':100*np.mean(d.code==1),'cyto':100*np.mean(d.code.isin([2,3])),'eryth':100*np.mean(d.code==4)}))
axA.plot(g.ng,g.cyto,'-o',color=CY,ms=5,label='cytotrophoblast')
axA.plot(g.ng,g.mega,'-o',color=MG,ms=5,label='megakaryocyte')
axA.plot(g.ng,g.eryth,'-o',color=ER,ms=5,label='erythrocyte (quarantined artifact)')
axA.set_xlabel('median detected genes, by decile',fontsize=BASE*0.88)
axA.set_ylabel('% of malignant cells assigned',fontsize=BASE*0.88)
axA.set_title('A',fontsize=BASE*0.92,loc='left',weight='bold')
for s_ in ['top','right']: axA.spines[s_].set_visible(False)
axA.grid(alpha=.2)

# B: threshold
axB=fig.add_subplot(gs[0,1])
thr=[1000,1500,2000,3000,4000]
mega=[100*np.mean(PC[PC.ng>=t].code==1) for t in thr]
cyto=[100*np.mean(PC[PC.ng>=t].code.isin([2,3])) for t in thr]
eryth=[100*np.mean(PC[PC.ng>=t].code==4) for t in thr]
axB.plot(thr,cyto,'-o',color=CY,ms=5); axB.plot(thr,mega,'-o',color=MG,ms=5); axB.plot(thr,eryth,'-o',color=ER,ms=5)
axB.set_xlabel('cells kept: \u2265 N genes',fontsize=BASE*0.88)
axB.set_ylabel('% assigned',fontsize=BASE*0.88)
axB.set_title('B',fontsize=BASE*0.92,loc='left',weight='bold')
for s_ in ['top','right']: axB.spines[s_].set_visible(False)
axB.grid(alpha=.2)

# C: causal downsampling (Melanoma + Glioblastoma)
axC=fig.add_subplot(gs[0,2])
for name,ls in [('Melanoma','-'),('Glioblastoma','--')]:
    d=D[D.cancer==name].sort_values('ng_med')
    axC.plot(d.ng_med,d.cyto,ls,color=CY,marker='o',ms=4)
    axC.plot(d.ng_med,d.mega,ls,color=MG,marker='o',ms=4)
    axC.plot(d.ng_med,d.eryth,ls,color=ER,marker='o',ms=4)
axC.invert_xaxis()
axC.set_xlabel('\u2190 thinned   (median genes)',fontsize=BASE*0.88)
axC.set_ylabel('% assigned after downsampling',fontsize=BASE*0.88)
axC.set_title('C',pad=4,fontsize=BASE*0.92,loc='left',weight='bold')
from matplotlib.lines import Line2D
axC.legend(handles=[Line2D([0],[0],color='#444',ls='-',label='Melanoma'),
                    Line2D([0],[0],color='#444',ls='--',label='Glioblastoma')],
           fontsize=BASE*0.82,loc='upper left',frameon=False)
for s_ in ['top','right']: axC.spines[s_].set_visible(False)
axC.grid(alpha=.2)

# one shared cell-type key under the row; three in-panel legends covered the data at print size
h,l=axA.get_legend_handles_labels()
fig.legend(h,l,fontsize=BASE*0.86,loc='lower center',bbox_to_anchor=(0.5,-0.09),ncol=3,frameon=False)
plt.savefig(FIG+'/v7_fig_E1_complexity.png',dpi=190,bbox_inches='tight'); plt.close()
F.report('Fig S7 complexity',F.panel_size(0.99,2.35),0.99)
print('saved v7_fig_E1_complexity.png')
print('threshold mega:',[round(x,1) for x in mega],'cyto:',[round(x,1) for x in cyto])
