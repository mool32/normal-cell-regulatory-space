"""Honest convergence MAP (replaces the self-referential polar figure).
A  UMAP of the 154 normal cell types in the full 9-D program space (cosine metric) — position =
   genuine state similarity, NOT distance to the duet. Colour = convergence heat. The 25 malignant
   centroids are PROJECTED onto the same map (.transform) as diamonds coloured by origin class;
   landmark normal cells (neural, hepatic, secretory, immune) are labelled so the map can be read.
B  3-D landscape: convergence heat as height over the same UMAP plane — every peak labelled.
No polar coordinates, no duet-centring, no trajectories, no cell-of-origin in the layout. Run from repo root.
"""
import sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
import umap
sys.path.insert(0, __import__('pathlib').Path(__file__).resolve().parents[1].as_posix()+'/tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); TAB=ROOT+'/results/cancer_v6_6_extras'; FIG=ROOT+'/manuscript/figures'
e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
PC=(Rm.values-mu)/sd@load; cells=Rm.index.tolist(); N=len(cells)
CC=pd.read_csv(TAB+'/cancer_centroids9d_25.csv',index_col=0)
if CC.shape[0]==9: CC=CC.T
_M=pd.read_csv(TAB+'/panel_25_top1_matrix.csv',index_col=0)
assert _M.max().max()>1.5, 'panel_25_top1_matrix.csv is expected in PERCENT (see scripts/E41_regen_ranking.py)'
votes=pd.Series(_M.mean(1)).reindex(cells).fillna(0).values   # percent of malignant cells, already
POST={'erythrocytes','reticulocytes','platelets'}; postm=np.array([c in POST for c in cells])
DUET=['megakaryocytes','cytotrophoblasts']; GERM=['late spermatids','early spermatids','oocytes']
EMB=[c for c in cells if any(k in c for k in ['trophoblast','spermat','germ','oogon','oocyte']) and c not in DUET]
vv=votes.copy(); vv[postm]=np.nan

# --- UMAP on the 154 normals (full 9-D, cosine); project the 25 cancers onto it ---
reducer=umap.UMAP(n_neighbors=15,min_dist=0.45,metric='cosine',random_state=0)
E=reducer.fit_transform(PC)               # 154 x 2
EC=reducer.transform(CC.values)           # 25 x 2 (cancers placed on the normal map)

OCLASS={'HCC':'epithelial','BRCA':'epithelial','RCC':'epithelial','Prostate':'epithelial','HNSCC':'epithelial',
 'SqLung':'epithelial','TNBC':'epithelial','PDAC':'epithelial','NSCLC':'epithelial','Colon':'epithelial',
 'Gastric':'epithelial','OvarianSerous':'epithelial','BCC':'epithelial',
 'B-ALL':'hematopoietic','DLBCL':'hematopoietic','AML':'hematopoietic','MM':'hematopoietic',
 'Osteosarcoma':'mesenchymal','Ewing':'mesenchymal','GermCell':'germ-cell','Seminoma':'germ-cell',
 'SCLC':'neuroendocrine','Neuroblastoma':'neural-crest','Melanoma':'neural-crest','Glioblastoma':'glial'}
CLASS_COL={'epithelial':'#2e86c1','hematopoietic':'#c0392b','mesenchymal':'#27ae60','germ-cell':'#8e44ad',
 'neuroendocrine':'#e67e22','neural-crest':'#16a085','glial':'#7f8c8d'}
ccol=[CLASS_COL.get(OCLASS.get(cx,'epithelial'),'#999') for cx in CC.index]

# landmark cells to orient the reader (hot + cold reference points)
LAND={'megakaryocytes':('megakaryocyte','#c0392b'),'cytotrophoblasts':('cytotrophoblast','#c0392b'),
 'gastric progenitor cells':('gastric progenitor','#7a4a00'),'pancreatic acinar cells':('pancreatic acinar','#7a4a00'),
 'gastric chief cells':('gastric chief','#7a4a00'),'late spermatids':('late spermatids','#2471a3'),
 'oocytes':('oocytes','#2471a3'),'hofbauer cells':('hofbauer','#7a4a00'),
 'brain inhibitory neurons':('neurons','#555'),'hepatocytes':('hepatocytes','#555'),
 'basal keratinocytes':('keratinocytes','#555'),'b-cells':('B cells','#555'),
 'adrenal medulla cells':('adrenal medulla','#555')}

from scipy.spatial import ConvexHull
# geometric medoid of the manifold (the true centre) — cosine
Cn=PC/(np.linalg.norm(PC,axis=1,keepdims=True)+1e-12); Dm=1-Cn@Cn.T
medoid_i=int(np.argmin(np.where(postm,np.inf,Dm.mean(1)))); medoid_name=cells[medoid_i]
VMAX=float(np.nanpercentile(vv[~postm],97))
m=~postm
hotthr=np.nanpercentile(vv[m],88)
anch=np.array([c in DUET for c in cells])
hot=m & (vv>=hotthr) & ~anch
cold=m & ~hot & ~anch

def draw(ax,inset=False):
    ax.scatter(E[cold,0],E[cold,1],s=(14 if inset else 22),c='#dfe3e6',edgecolor='none',zorder=2)
    hp=E[m & (vv>=np.nanpercentile(vv[m],90))]          # basin = top-decile hot cells
    if len(hp)>=4:                                        # drop spatial outliers (MAD) so the hull stays tight
        c0=np.median(hp,axis=0); dd=np.hypot(*(hp-c0).T)
        hp=hp[dd<=np.median(dd)+2.0*(np.median(np.abs(dd-np.median(dd)))+1e-9)]
    if len(hp)>=3:
        pass   # convex hull removed in v11: it suggested hand-picked vertices (figure plan, Fig. 3A)
    sc=ax.scatter(E[hot,0],E[hot,1],c=vv[hot],s=(70 if inset else 100),cmap='YlOrRd',vmin=0,vmax=VMAX,
                  edgecolor='#8a5a00',lw=.4,zorder=4)
    ax.scatter(EC[:,0],EC[:,1],marker='D',s=(34 if inset else 60),c='#2c3e50',edgecolor='white',lw=.6,zorder=6)
    for a in DUET:
        i=cells.index(a); ax.scatter(E[i,0],E[i,1],marker='*',s=(230 if inset else 430),c='#7b241c',edgecolor='black',lw=1.0,zorder=8)
    return sc

BASE=F.use(); FSA=F.panel_size(0.56,1.26)
figA,axA=plt.subplots(figsize=FSA)
sc=draw(axA)
# manifold centre (medoid) — receives no convergence
axA.scatter(E[medoid_i,0],E[medoid_i,1],marker='P',s=150,c='white',edgecolor='#2471a3',lw=1.8,zorder=8)
axA.annotate('manifold centre → 0 votes\n('+medoid_name.replace(' cells','')+')',(E[medoid_i,0],E[medoid_i,1]),
    xytext=(E[medoid_i,0]+0.6,E[medoid_i,1]+1.6),fontsize=BASE*0.86,color='#2471a3',weight='bold',ha='center',
    arrowprops=dict(arrowstyle='->',color='#2471a3',lw=1.0))
# anchor labels with leader lines so they never float
for a,(lab2) in [('megakaryocytes','megakaryocyte'),('cytotrophoblasts','cytotrophoblast')]:
    i=cells.index(a); axA.annotate(lab2,(E[i,0],E[i,1]),xytext=(E[i,0]-2.4,E[i,1]-1.8),
        fontsize=BASE*0.96,color='#7b241c',weight='bold',ha='center',
        arrowprops=dict(arrowstyle='-',color='#7b241c',lw=0.8))
# a few FAR-FIELD orientation landmarks only (the basin itself is named by the stars) — keep the panel uncluttered
ORIENT={'brain inhibitory neurons':'neurons','hepatocytes':'hepatocytes','b-cells':'B cells',
        'basal keratinocytes':'keratinocytes'}
for c,lab in ORIENT.items():
    if c in cells:
        i=cells.index(c); axA.annotate(lab,(E[i,0],E[i,1]),fontsize=BASE*0.86,color='#666',
            xytext=(4,4),textcoords='offset points',zorder=8)

axA.set_xlabel('UMAP-1',fontsize=BASE*1.00); axA.set_ylabel('UMAP-2',fontsize=BASE*1.00)
cb=plt.colorbar(sc,ax=axA,fraction=0.040,pad=0.015)
cb.set_label('convergence heat\n(% of malignant cells, mean over the 25 cancers)',fontsize=BASE*0.80)
cb.ax.tick_params(labelsize=BASE*0.80)
leg=[Line2D([0],[0],marker='D',color='w',markerfacecolor='#2c3e50',markeredgecolor='white',markersize=9,label='25 cancers (malignant states)'),
     Line2D([0],[0],marker='*',color='w',markerfacecolor='#7b241c',markeredgecolor='k',markersize=14,label='two leading targets'),
     Line2D([0],[0],marker='P',color='w',markerfacecolor='white',markeredgecolor='#2471a3',markersize=10,label='manifold centre (0 votes)')]
axA.legend(handles=leg,fontsize=BASE*0.84,loc='upper left',framealpha=.95,borderpad=0.5)
figA.savefig(FIG+'/v11_fig3A_map.png',dpi=200,bbox_inches='tight'); plt.close(figA)
F.report('Fig 3A map',FSA,0.56)
# diagnostics: are cancers near the core on the map?
dx=np.array([np.hypot(*(EC[k]-E[cells.index('cytotrophoblasts')])) for k in range(len(EC))])
print('saved v11_fig3A_map.png (Fig 4) | medoid=%s | median cancer→cytotrophoblast UMAP dist %.2f (span %.1f)'%(
      medoid_name,np.median(dx),np.hypot(E[:,0].ptp(),E[:,1].ptp())))
