"""Fig. 3 panels C and E (v11, figstyle) — rebuilt from scratch.

Why they are rebuilt rather than reused: both were still the legacy v7 PNGs
(v7_figS2_2axis.png, v7_fig1C_heatmap_z.png), produced by scripts/build_v7_fig1.py, which reads
results/cancer_v6_6_extras/panel_25_top1_matrix.csv and multiplies it by 100. Since
scripts/E41_regen_ranking.py regenerated that matrix it is already in PERCENT, so the legacy panels
were off by 100x: panel C drew 13 of its 14 points outside the axes (one grey dot survived), and the
heatmap's "single-cancer resonance" crosses fired at >0.25% instead of >25%, covering the panel.
That is review blocker #6 (Fig. 3B vs 3C vs Table S2). Here the matrix is read as percent, once.

C  breadth (how many of the 25 cancers rank a normal type in their top-5) x depth (its median
   per-cancer share of malignant cells).  E  per-cancer z-score heat of the top candidates.
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
sys.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix()); import figstyle as F
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; FIG=ROOT/'manuscript/figures'; TAB=ROOT/'results/cancer_v6_6_extras'

M=pd.read_csv(TAB/'panel_25_top1_matrix.csv',index_col=0)   # PERCENT of each cancer's malignant cells
assert M.max().max()>1.5, 'panel_25_top1_matrix.csv is expected in percent (see E41_regen_ranking.py)'
POSTCELL={'erythrocytes','reticulocytes','platelets'}
rec=pd.Series(0,index=M.index)
for c in M.columns: rec[M[c].sort_values(ascending=False).head(5).index]+=1
med=M.median(axis=1)                                         # already percent
DUET=['megakaryocytes','cytotrophoblasts']
NEIGH=['gastric progenitor cells','esophageal apical cells','pancreatic acinar cells',
       'gastric chief cells','epididymal clear cells','late spermatids']
def col(c): return '#8e1b1b' if c in DUET else '#e69138' if c in NEIGH else '#b9c0c6'
SHORT={'megakaryocytes':'megakaryocyte','cytotrophoblasts':'cytotrophoblast',
 'gastric progenitor cells':'gastric progenitor','esophageal apical cells':'esophageal apical',
 'pancreatic acinar cells':'pancreatic acinar','late spermatids':'late spermatid',
 'gastric chief cells':'gastric chief','epididymal clear cells':'epididymal clear',
 'hofbauer cells':'hofbauer','migrating cytotrophoblasts':'migrating cytotrophoblast',
 'differentiating spermatogonia':'diff. spermatogonia','early primary spermatocytes':'early spermatocyte'}

# ======================= C: breadth x depth =======================
BP=F.use(); FS=F.panel_size(0.47,1.30)
show=[c for c in rec.sort_values(ascending=False).index if c not in POSTCELL and rec[c]>=4]
# only the named candidates carry a label; the tail is drawn but left unlabelled (it is listed in
# panel B by rank and in full in Table S2) — labelling all twelve collided with the axes.
LABEL=['cytotrophoblasts','megakaryocytes','gastric progenitor cells','late spermatids',
       'gastric chief cells','esophageal apical cells','pancreatic acinar cells']
fig,ax=plt.subplots(figsize=FS)
xs=np.array([rec[c] for c in show],float); ys=np.array([med[c] for c in show],float)
ax.scatter(xs,ys,s=[170 if c in DUET else 100 for c in show],
           color=[col(c) for c in show],edgecolor='black',lw=1.0,zorder=4)
ax.set_xlim(0,24.5); ax.set_ylim(0,13.6)
ax.add_patch(Ellipse((18.0,(med['cytotrophoblasts']+med['megakaryocytes'])/2),
                     3.2,5.6,fill=False,ec='#8e1b1b',lw=2.2,ls='--',zorder=3))
ax.annotate('the two-cell core leads\non both axes',(19.6,10.2),xytext=(12.2,13.2),
            fontsize=BP*0.80,color='#8e1b1b',fontweight='bold',ha='center',va='top',
            arrowprops=dict(arrowstyle='->',color='#8e1b1b',lw=1.5))
# labelled points get a label to the LOWER-RIGHT, de-collided downward within the axes
rowh=0.80
placed=[]
for c in sorted(LABEL,key=lambda c:-med[c]):
    ty=med[c]-0.25
    for _,py in placed:
        if abs(ty-py)<rowh: ty=py-rowh
    placed.append((c,ty))
LX=float(max(rec[c] for c in LABEL if c not in DUET))+1.4   # one label column, clear of every point
for c,ty in placed:
    clr='#8e1b1b' if c in DUET else '#b06000'
    ax.annotate(SHORT.get(c,c.replace(' cells','')),(rec[c],med[c]),
                xytext=((rec[c]+0.8) if c in DUET else LX,ty),
                ha='left',va='center',fontsize=BP*0.74,color=clr,
                fontweight=('bold' if c in DUET else 'normal'),
                arrowprops=dict(arrowstyle='-',color=clr,lw=0.7,shrinkA=1,shrinkB=4))
ax.text(0.025,0.60,'unlabelled points:\nthe tail of the ranking\n(by rank in panel B,\nin full in Table S2)',
        transform=ax.transAxes,ha='left',va='top',fontsize=BP*0.68,color='#777',style='italic')
ax.set_xlabel('breadth — cross-cancer recurrence\n(a top-5 target in how many of the 25 cancers)',fontsize=BP*0.90)
ax.set_ylabel('depth — median per-cancer share\nof malignant cells (%)',fontsize=BP*0.90)
ax.grid(alpha=.22,zorder=0)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
plt.savefig(FIG/'v11_fig3C_breadthdepth.png',bbox_inches='tight'); plt.close()
F.report('Fig 3C breadth x depth',FS,0.47)
print('   points drawn: %d, labelled: %d (the legacy panel drew 1 of 14)'%(len(show),len(LABEL)))

# ======================= E: per-cancer z-score heat =======================
rows=['cytotrophoblasts','megakaryocytes','gastric progenitor cells','esophageal apical cells',
      'pancreatic acinar cells','gastric chief cells','late spermatids','epididymal clear cells']
H=M.loc[rows].copy()                                          # percent
Z=(H-H.mean(axis=0))/(H.std(axis=0)+1e-9)
LIN={'B-ALL':'blood','AML':'blood','DLBCL':'blood','MM':'blood',
 'HCC':'epithelial','BRCA':'epithelial','RCC':'epithelial','Prostate':'epithelial','HNSCC':'epithelial',
 'SqLung':'epithelial','TNBC':'epithelial','PDAC':'epithelial','NSCLC':'epithelial','Colon':'epithelial',
 'Gastric':'epithelial','OvarianSerous':'epithelial','BCC':'epithelial',
 'Osteosarcoma':'sarcoma','Ewing':'sarcoma','GermCell':'germ cell','Seminoma':'germ cell',
 'SCLC':'neuroendocrine','Neuroblastoma':'neural crest','Melanoma':'neural crest','Glioblastoma':'glial'}
megv=M.loc['megakaryocytes']
glist=sorted(set(LIN[c] for c in M.columns),
             key=lambda L:-megv[[c for c in M.columns if LIN[c]==L]].mean())
order_c=[]; bounds=[]
for L in glist:
    cs=sorted([c for c in M.columns if LIN[c]==L],key=lambda c:-megv[c])
    order_c+=cs; bounds.append((L,len(order_c)))
order_c=pd.Index(order_c); Z=Z[order_c]; Hc=H[order_c]

BP=F.use(); FS=F.panel_size(0.60,2.05)
fig,ax=plt.subplots(figsize=FS)
im=ax.imshow(Z.values,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
CROSS=25.0     # percent of that cancer's malignant cells — a single-cancer resonance
nx=0
for i,r in enumerate(rows):
    if r in DUET: continue
    for jx,c in enumerate(order_c):
        if Hc.loc[r,c]>CROSS:
            nx+=1
            ax.text(jx,i,'×',ha='center',va='center',color='white',fontsize=BP*1.05,fontweight='bold',
                    path_effects=[pe.withStroke(linewidth=2.0,foreground='black')],zorder=6)
# lineage shown as a colour key above the panel: bracket labels could not be fitted over the
# one-column blocks (glial, neuroendocrine) without stacking them four levels deep
from matplotlib.patches import Patch
LCOL={'blood':'#c0392b','epithelial':'#2e86c1','sarcoma':'#27ae60','germ cell':'#8e44ad',
      'neuroendocrine':'#e67e22','neural crest':'#16a085','glial':'#7f8c8d'}
x0=0
for L,end in bounds:
    if end<len(order_c): ax.axvline(end-0.5,color='#333',lw=1.3,zorder=5)
    ax.plot([x0-0.40,end-0.60],[-0.72]*2,color=LCOL[L],lw=4.0,clip_on=False,
            zorder=6,solid_capstyle='butt')
    x0=end
ax.legend(handles=[Patch(fc=LCOL[L],ec='none',label=L) for L,_ in bounds],
          loc='lower center',bbox_to_anchor=(0.5,1.015),ncol=len(bounds),frameon=False,
          fontsize=BP*0.74,columnspacing=1.0,handlelength=1.1,handletextpad=0.4,handleheight=0.9)
ax.set_ylim(len(rows)-0.5,-0.5)
ax.set_xticks(range(len(order_c))); ax.set_xticklabels(order_c,rotation=55,ha='right',fontsize=BP*0.78)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([SHORT.get(r,r.replace(' cells','')) for r in rows],fontsize=BP*0.84)
for i,r in enumerate(rows):
    if r in DUET:
        ax.get_yticklabels()[i].set_color('#8e1b1b'); ax.get_yticklabels()[i].set_weight('bold')
for t,c in zip(ax.get_xticklabels(),order_c): t.set_color(LCOL[LIN[c]])
cb=plt.colorbar(im,ax=ax,fraction=0.022,pad=0.012); cb.set_label('per-cancer z-score',fontsize=BP*0.84)
cb.ax.tick_params(labelsize=BP*0.78)
plt.savefig(FIG/'v11_fig3E_heatmap.png',bbox_inches='tight'); plt.close()
F.report('Fig 3E heatmap',FS,0.60)
print('   crosses drawn: %d at the >%.0f%% rule (the legacy panel fired at >0.25%%)'%(nx,CROSS))
