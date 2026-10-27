"""Figure 1 panels for paper3 v7 (platelet-removed FINAL panel, Gudkov reorder).
A = painted normal manifold (ENTRY, Gudkov)   B = recurrence step   C = heatmap (white-cross resonances)
SI = breadth x depth.  Run from repo root.
"""
import sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import sys as _sf; _sf.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix()
sys.path.insert(0, ROOT+'/tool/perceptome2'); import perceptome as pct
FIG=ROOT+'/manuscript/figures'; TAB=ROOT+'/results/cancer_v6_6_extras'

M=pd.read_csv(TAB+'/panel_25_top1_matrix.csv',index_col=0)   # 153 x 25, in PERCENT since E41_regen_ranking.py
# v11: one convention everywhere (Note S-K) — the quarantined post-cellular bodies are NOT candidates,
# so they are dropped BEFORE the per-cancer top-5 is taken, exactly as in scripts/E41_regen_ranking.py.
# the matrix is now generated under the same convention by scripts/E41_regen_ranking.py
# post-cellular bodies (no nucleus): excluded from ranking + drawn as quarantine squares
POSTCELL={'erythrocytes','reticulocytes','platelets'}  # platelets already absent from panel_25
rec=pd.Series(0,index=M.index)
for c in M.columns:
    rec[M[c].sort_values(ascending=False).head(5).index]+=1
med=M.median(axis=1)
DUET=['megakaryocytes','cytotrophoblasts']
NEIGH=['gastric progenitor cells','esophageal apical cells','pancreatic acinar cells','gastric chief cells','epididymal clear cells','late spermatids']
def col(c):
    if c in DUET: return '#c0392b'
    if c in NEIGH: return '#e69138'
    return '#bdc3c7'

# ---------- Fig 1A: convergence-painted normal manifold (ENTRY) ----------
e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json'))
mods=e['module_order']; load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
hpa=pct.load_hpa_perceptivity(mode='single_cell_scaled'); Rm=hpa['R'][mods]
muM=Rm.mean(0).values; sdM=Rm.std(0).values+1e-12; PC=(Rm.values-muM)/sdM@load; idxM=Rm.index.tolist()
votes=M.mean(axis=1).reindex(idxM).fillna(0).values
postm=np.array([c in POSTCELL for c in idxM]); vv=votes.copy(); vv[postm]=0
EMB=[c for c in idxM if any(k in c for k in ['trophoblast','spermat','germ','oogon','oocyte']) and c not in DUET]
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 1.22))
sc1=ax.scatter(PC[:,0],PC[:,1],c=vv,s=74,cmap='YlOrRd',edgecolor='#999',lw=.3,vmin=0,vmax=np.percentile(vv,98),zorder=2)
# embryonic: numbered crosses + legend list
ei=[idxM.index(c) for c in EMB]
ax.scatter(PC[ei,0],PC[ei,1],s=26,c='#2471a3',marker='x',lw=1.2,zorder=5)
for n,c in enumerate(EMB,1):
    j=idxM.index(c); ax.annotate(str(n),(PC[j,0],PC[j,1]),fontsize=BASE*0.86,color='#1a5276',xytext=(2.5,2.5),textcoords='offset points',weight='bold')
di=[idxM.index('cytotrophoblasts')]; ax.scatter(PC[di,0],PC[di,1],s=240,facecolor='none',edgecolor='#c0392b',lw=2.6,zorder=6,label='cytotrophoblast (anchor)')
pm=np.where(postm)[0]; ax.scatter(PC[pm,0],PC[pm,1],s=58,facecolor='none',edgecolor='#777',lw=1.3,marker='s',zorder=4,label='post-cellular (erythrocyte, platelet) — excluded')
labM={'megakaryocytes':('megakaryocyte',(2.0,-1.9)),'cytotrophoblasts':('cytotrophoblast',(-1.2,2.4)),
 'gastric progenitor cells':('gastric prog.',(-4.0,1.7)),'esophageal apical cells':('esoph. apical',(-4.0,-0.8)),
 'pancreatic acinar cells':('panc. acinar',(1.6,1.9)),'gastric chief cells':('gastric chief',(-1.0,-2.6)),
 'late spermatids':('late spermatid (no.7, low-complexity)',(-3.0,-1.6))}
for c,(t,(dx,dy)) in labM.items():
    if c in idxM:
        j=idxM.index(c); clr='#c0392b' if c in DUET else '#7a5200'
        ax.annotate(t,(PC[j,0],PC[j,1]),fontsize=BASE*0.943,xytext=(PC[j,0]+dx,PC[j,1]+dy),color=clr,arrowprops=dict(arrowstyle='-',color=clr,lw=0.5,alpha=0.7))
ax.set_xlabel('PC1',fontsize=BASE*1.097); ax.set_ylabel('PC2',fontsize=BASE*1.097)
ax.set_title('A',fontsize=BASE*1.168,loc='left',weight='bold')
ax.legend(fontsize=BASE*0.943,loc='lower left',framealpha=.95)
cb=plt.colorbar(sc1,ax=ax,fraction=0.045,pad=0.02); cb.set_label('convergence received\n(mean top-1 votes from 25 cancers, %)',fontsize=BASE*0.983)
emb_list='  '.join(f'{n}.{c.replace(" cells","").replace("cytotrophoblasts","cyto")[:14]}' for n,c in enumerate(EMB,1))
ax.text(0.5,-0.105,'Colour = convergence (top-1 votes received).  Position = PC1–PC2 (context).  Blue ×: 10 embryonic/germline cells (cytotrophoblast, the 11th, is a red core ring) —\n'+emb_list,
        transform=ax.transAxes,ha='center',va='top',fontsize=BASE*0.825,style='italic',color='#444')
plt.tight_layout(); plt.savefig(FIG+'/v7_fig1A_space.png',dpi=200,bbox_inches='tight'); plt.close()
print('Fig1A (manifold, ENTRY) done')

# ---------- Fig 1B: recurrence ranking ----------
order=[c for c in rec.sort_values(ascending=False).index if c not in POSTCELL][:11]
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 1.95))
ax.bar(range(len(order)),[rec[c] for c in order],color=[col(c) for c in order],edgecolor='black',lw=0.6)
ax.axhline(15,color='#c0392b',ls='--',lw=1,alpha=0.55)
ax.annotate('',xy=(1,17),xytext=(2,13),arrowprops=dict(arrowstyle='<->',color='#c0392b',lw=1.2))
ax.text(3.4,15.5,'gap 17 → 13',ha='left',va='center',fontsize=BASE*1.022,color='#c0392b',weight='bold')
ax.set_xticks(range(len(order))); ax.set_xticklabels([c.replace(' cells','').replace(' ','\n',1) for c in order],rotation=40,ha='right',fontsize=BASE*0.983)
ax.set_ylabel('cross-cancer recurrence\n(# of 25 cancers, cell in top-5 targets)',fontsize=BASE*1.097)
ax.set_title('',fontsize=BASE*1.168,loc='left',weight='bold')
ax.legend(handles=[Patch(fc='#c0392b',label='cytotrophoblast:\nleads on both axes'),Patch(fc='#e69138',label='enriched neighborhood'),Patch(fc='#bdc3c7',label='other targets (tail)')],fontsize=BASE*0.983,loc='upper right',framealpha=0.95)
ax.text(0.985,0.60,'post-cellular bodies (erythrocyte, platelet)\nexcluded — no nucleus (Methods)',transform=ax.transAxes,ha='right',va='top',fontsize=BASE*0.91,style='italic',color='#555')
plt.tight_layout(); plt.savefig(FIG+'/v7_fig1B_recurrence.png',dpi=200); plt.close()
print('Fig1B (recurrence) done')

# ---------- Fig 1C: heatmap with WHITE-CROSS resonances (no double color-code) ----------
rows=['cytotrophoblasts','megakaryocytes','gastric progenitor cells','esophageal apical cells','pancreatic acinar cells','gastric chief cells','late spermatids','epididymal clear cells']
# the SI version carries the quarantined post-cellular body that the main-figure panel (Fig. 3E) omits
rows=rows+[c for c in ['erythrocytes'] if c in M.index]
H=M.loc[rows].copy()
Z=(H-H.mean(axis=0))/(H.std(axis=0)+1e-9)
# --- order columns by cancer LINEAGE so same-type cancers are adjacent; lineages ordered by
#     mean megakaryocyte convergence, so the high-megakaryocyte blood cancers form one block ---
LIN={'B-ALL':'blood','AML':'blood','DLBCL':'blood','MM':'blood',
 'HCC':'epithelial','BRCA':'epithelial','RCC':'epithelial','Prostate':'epithelial','HNSCC':'epithelial',
 'SqLung':'epithelial','TNBC':'epithelial','PDAC':'epithelial','NSCLC':'epithelial','Colon':'epithelial',
 'Gastric':'epithelial','OvarianSerous':'epithelial','BCC':'epithelial',
 'Osteosarcoma':'sarcoma','Ewing':'sarcoma','GermCell':'germ-cell','Seminoma':'germ-cell',
 'SCLC':'neuroendo','Neuroblastoma':'neural-crest','Melanoma':'neural-crest','Glioblastoma':'glial'}
megv=M.loc['megakaryocytes']
glist=sorted(set(LIN.get(c,'other') for c in M.columns),
             key=lambda L:-megv[[c for c in M.columns if LIN.get(c,'other')==L]].mean())
order_c=[]; bounds=[]
for L in glist:
    cs=sorted([c for c in M.columns if LIN.get(c,'other')==L], key=lambda c:-megv[c])
    order_c+=cs; bounds.append((L,len(order_c)))
order_c=pd.Index(order_c)
Z=Z[order_c]; Hc=H[order_c]
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 1.85))
im=ax.imshow(Z.values,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
# white crosses on single-cancer resonances (raw fraction > 25% in a non-core row)
import matplotlib.patheffects as pe
for i,r in enumerate(rows):
    if r in DUET: continue
    for jx,c in enumerate(order_c):
        if Hc.loc[r,c]>25:
            ax.text(jx,i,'×',ha='center',va='center',color='white',fontsize=BASE*1.606,weight='bold',
                    path_effects=[pe.withStroke(linewidth=2.2,foreground='black')],zorder=6)
# lineage block separators + labels on top
# the lineage headers overlapped for the one- and two-column blocks; stagger them over as many
# levels as needed (the main-figure twin, Fig. 3E, solves this with a colour key instead)
x0=0; occupied=[]
for L,end in bounds:
    if end<len(order_c): ax.axvline(end-0.5,color='#333',lw=1.4,zorder=5)
    mid=(x0+end-1)/2.0; half=0.46*len(L)+0.4
    lvl=0
    while lvl<len(occupied) and mid-half < occupied[lvl]: lvl+=1
    if lvl==len(occupied): occupied.append(-9.9)
    occupied[lvl]=mid+half
    y=-0.85-0.95*lvl
    ax.plot([x0-0.40,end-0.60],[y-0.12]*2,color='#333',lw=1.0,clip_on=False,zorder=6)
    ax.text(mid,y,L,ha='center',va='bottom',fontsize=BASE*1.237,weight='bold',color='#333',clip_on=False)
    x0=end
ax.set_ylim(len(rows)-0.5,-0.5)
ax.set_xticks(range(len(order_c))); ax.set_xticklabels(order_c,rotation=45,ha='right',fontsize=BASE*1.237)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r.replace(' cells','') for r in rows],fontsize=BASE*1.271)
for i,r in enumerate(rows):
    if r in DUET: ax.get_yticklabels()[i].set_color('#c0392b'); ax.get_yticklabels()[i].set_weight('bold')
cb=plt.colorbar(im,fraction=0.025,pad=0.01); cb.set_label('per-cancer z-score',fontsize=BASE*1.271)
plt.tight_layout(); plt.savefig(FIG+'/v7_fig1C_heatmap_z.png',dpi=200); plt.close()
print('Fig1C (heatmap, white crosses) done')

# ---------- SI: breadth x depth (Fig. S2) ----------
# This is the SI companion to Fig. 3C: the same two axes, but with EVERY one of the fourteen
# candidates labelled (the main-figure panel labels seven and leaves the tail unlabelled).
import sys as _s2; _s2.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as _F
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 1.75))
show=[c for c in rec.sort_values(ascending=False).index if c not in POSTCELL][:14]
for c in show: ax.scatter(rec[c],med[c],s=70,color=col(c),edgecolor='black',lw=0.5,zorder=3)
ax.set_xlim(0,23); ax.set_ylim(0,max(13.0,float(med.max())*1.16))
ax.add_patch(Ellipse((18.0,float(med['cytotrophoblasts'])),2.4,2.0,fill=False,ec='#c0392b',lw=2.4,zorder=2))
ax.text(0.03,0.97,'leads on both axes',transform=ax.transAxes,ha='left',va='top',
        fontsize=BASE*1.336,color='#c0392b',weight='bold')
lab={'megakaryocytes':'megakaryocyte','cytotrophoblasts':'cytotrophoblast',
 'gastric progenitor cells':'gastric progenitor','esophageal apical cells':'esophageal apical',
 'pancreatic acinar cells':'pancreatic acinar','late spermatids':'late spermatid',
 'gastric chief cells':'gastric chief','epididymal clear cells':'epididymal clear',
 'hofbauer cells':'hofbauer','migrating cytotrophoblasts':'migrating cytotrophoblast',
 'differentiating spermatogonia':'diff. spermatogonia','early primary spermatocytes':'early spermatocyte',
 'undifferentiated spermatogonia':'undiff. spermatogonia','erythrocyte progenitors':'erythrocyte progenitor',
 'enteric stem cells':'enteric stem','syncytiotrophoblasts':'syncytiotrophoblast','paneth cells':'paneth'}
txt=[]
for c in show:
    clr=('#c0392b' if c in DUET else '#b06000' if c in NEIGH else '#666')
    txt.append(ax.text(rec[c],med[c],lab.get(c,c.replace(' cells','')),fontsize=BASE*1.203,color=clr,
                       weight=('bold' if c in DUET else 'normal'),zorder=8))
ax.set_xlabel('breadth — cross-cancer recurrence (# of 25)',fontsize=BASE*1.399)
ax.set_ylabel('median per-cancer fraction (%)',fontsize=BASE*1.399)
_F.declutter(ax,txt,avoid=[(rec[c],med[c]) for c in show])
plt.savefig(FIG+'/v7_figS2_2axis.png',dpi=200,bbox_inches='tight'); plt.close()
print('SI (2axis) done')
