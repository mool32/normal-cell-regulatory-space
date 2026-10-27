"""Fig. 2A/2B (v11) — the same normal-cell map twice, side by side, never overlaid.
A: the 154 normal cell types, coloured by biological family (the space, before cancer).
B: the 25 malignant states projected into that same space, as diamonds coloured by origin class,
   with a dashed outline of the region they fall into. Per the figure plan: no open circles for the
   tissues of origin, no cytotrophoblast or megakaryocyte labels here, large fonts throughout."""
import sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
import matplotlib.pyplot as plt, matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import umap
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
sys.path.insert(0, ROOT+'/tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
BASE=F.use()
e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
PC=(Rm.values-mu)/sd@load; cells=Rm.index.tolist()
POST={'erythrocytes','platelets','reticulocytes'}
_k=[i for i in range(len(cells)) if cells[i] not in POST]
PC=PC[_k]; cells=[cells[i] for i in _k]; N=len(cells)
fam=json.load(open(ROOT+'/results/fam.json'))
FAMS=['germline / embryonic','stem / progenitor','erythroid / megakaryocyte','blood / immune',
 'epithelial — secretory/glandular','epithelial — barrier/squamous','respiratory','hepatic','renal',
 'endocrine','neural','stromal / mesenchymal','endothelial','muscle','melanocyte']
palette=dict(zip(FAMS,[cm.tab20(i/20) for i in range(len(FAMS))]))
CC=pd.read_csv(ROOT+'/results/cancer_v6_6_extras/cancer_centroids9d_25.csv',index_col=0)
if CC.shape[0]==9: CC=CC.T
OCLASS={'HCC':'epithelial','BRCA':'epithelial','RCC':'epithelial','Prostate':'epithelial','HNSCC':'epithelial',
 'SqLung':'epithelial','TNBC':'epithelial','PDAC':'epithelial','NSCLC':'epithelial','Colon':'epithelial',
 'Gastric':'epithelial','OvarianSerous':'epithelial','BCC':'epithelial','B-ALL':'hematopoietic','DLBCL':'hematopoietic',
 'AML':'hematopoietic','MM':'hematopoietic','Osteosarcoma':'mesenchymal','Ewing':'mesenchymal','GermCell':'germ-cell',
 'Seminoma':'germ-cell','SCLC':'neuroendocrine','Neuroblastoma':'neural-crest','Melanoma':'neural-crest','Glioblastoma':'glial'}
OC_COL={'epithelial':'#2e86c1','hematopoietic':'#c0392b','mesenchymal':'#27ae60','germ-cell':'#8e44ad',
        'neuroendocrine':'#e67e22','neural-crest':'#16a085','glial':'#7f8c8d'}
reducer=umap.UMAP(n_neighbors=15,min_dist=0.4,metric='cosine',random_state=0)
E=reducer.fit_transform(PC); EC=reducer.transform(CC.values)

def base(ax, grey=False):
    for f in FAMS:
        idx=[k for k in range(N) if fam[cells[k]]==f]
        if not idx: continue
        ax.scatter(E[idx,0],E[idx,1],s=105 if not grey else 72,
                   color='#d6d6d6' if grey else palette[f], edgecolor='white',lw=.6,
                   label=None if grey else f, zorder=3)
    ax.set_xlabel('UMAP-1'); ax.set_ylabel('UMAP-2')
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

# ---- A ----
FS2A=F.panel_size(0.5,1.18)
fig,ax=plt.subplots(figsize=FS2A); base(ax)
# the fifteen family labels were drawn at each family's median and overlapped in a block of six;
# they are now de-collided against their real rendered extents, with a leader back to the median.
txt=[]; anc=[]
for f in FAMS:
    idx=[k for k in range(N) if fam[cells[k]]==f]
    if not idx: continue
    mx,my=np.median(E[idx,0]),np.median(E[idx,1])
    lab=f.replace('epithelial — secretory/glandular','epithelial\n(secretory)').replace('epithelial — barrier/squamous','epithelial\n(barrier)').replace('erythroid / megakaryocyte','erythroid /\nmegakaryocyte').replace('stromal / mesenchymal','stromal /\nmesenchymal').replace('germline / embryonic','germline /\nembryonic')
    txt.append(ax.text(mx,my,lab,fontsize=BASE*0.82,weight='bold',color='#222',ha='center',va='center',
            bbox=dict(boxstyle='round,pad=0.2',fc='white',ec=palette[f],lw=1.3,alpha=.92),zorder=10))
    anc.append((mx,my))
xl,yl=ax.get_xlim(),ax.get_ylim()
ax.set_xlim(xl[0]-0.9*(xl[1]-xl[0])*0.06, xl[1]+0.9*(xl[1]-xl[0])*0.06)
ax.set_ylim(yl[0]-(yl[1]-yl[0])*0.06, yl[1]+(yl[1]-yl[0])*0.10)
F.declutter(ax,txt,avoid=E)
plt.savefig(FIG+'/v11_fig2A_normalmap.png',dpi=200,bbox_inches='tight'); plt.close()
F.report('Fig 2A normal map',FS2A,0.5)
print('wrote v11_fig2A_normalmap.png')

# ---- B ----
FS2B=F.panel_size(0.5,1.18)
fig,ax=plt.subplots(figsize=FS2B); base(ax,grey=True)
cx,cy=EC[:,0].mean(),EC[:,1].mean()
rx,ry=2.6*EC[:,0].std(),2.6*EC[:,1].std()
ax.add_patch(Ellipse((cx,cy),2*rx,2*ry,fill=True,fc='#f7d7d2',alpha=.45,ec='#c0392b',lw=2.2,ls='--',zorder=2))
for i,c in enumerate(CC.index):
    k=OCLASS.get(c,'epithelial')
    ax.scatter(EC[i,0],EC[i,1],marker='D',s=190,color=OC_COL[k],edgecolor='black',lw=1.0,zorder=6)
ax.text(cx,cy+ry+0.45,'the 25 malignant states',fontsize=BASE*1.00,color='#c0392b',
        fontweight='bold',ha='center',va='bottom')
h=[Line2D([0],[0],marker='D',color='w',markerfacecolor=v,markeredgecolor='black',markersize=13,label=k)
   for k,v in OC_COL.items()]
h.append(Line2D([0],[0],marker='o',color='w',markerfacecolor='#d6d6d6',markersize=12,label='normal cell types'))
# the legend sat on the malignant diamonds and cut the 'malignant states' label in half
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]+1.2)
ax.legend(handles=h,fontsize=BASE*0.84,loc='lower right',framealpha=.95,
          title='lineage of origin',title_fontsize=BASE*0.88,ncol=2,columnspacing=0.9,
          handletextpad=0.4,borderpad=0.5)
plt.savefig(FIG+'/v11_fig2B_malignant.png',dpi=200,bbox_inches='tight'); plt.close()
F.report('Fig 2B malignant',FS2B,0.5)
print('wrote v11_fig2B_malignant.png')
from scipy.spatial.distance import cdist
sp=lambda P: np.median(cdist(P,P)[np.triu_indices(len(P),1)])
print(f'  spread on the map: normal cells {sp(E):.1f}  vs malignant states {sp(EC):.1f}')
