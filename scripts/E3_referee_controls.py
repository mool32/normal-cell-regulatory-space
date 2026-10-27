"""Referee controls #3 (hard-assignment fragility) + #4 (magnitude/complexity artifact).
Re-scores the 25 cancers (canonical pipeline, bundled eigenspace) and asks whether the convergence onto
the cytotrophoblast+megakaryocyte core survives:
  (3a) SOFT assignment (softmax over cosine) instead of winner-take-all top-1;
  (3b) TOP-k membership and the RANK of the best core centroid per cell (is the core a real winner or a hair-margin one?);
  (3c) is the core just the most EXTREME/high-norm centroid attracting everything? (cosine ignores norm — reported to preempt);
  (4a) does core assignment track cell library COMPLEXITY (n_genes)? (the erythrocyte control, extended to all cells);
  (4b) MAGNITUDE-REMOVED control: re-project after within-cell centering of the 44 module scores -> does the core survive?
Reuses the full 25-cancer loader (backed + mask + to_memory; canonical score_modules; correlation-PCA bundled eigenspace).
"""
import scanpy as sc, anndata as ad, json, numpy as np, pandas as pd, sys, re, gc, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def stage(m): print(f"[stage] {m}",flush=True)
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']; NPC=9; CAP=4000
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
refPC=(Rm.values-mu)/sd@load; refPn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)
refnorm=np.linalg.norm(refPC,axis=1)                                   # ref centroid norms (for #3c)
print(f"core centroid 9-D norm percentile among 154: "+
      ", ".join(f"{refcells[i]} {100*(refnorm<refnorm[i]).mean():.0f}%" for i in core_idx),flush=True)

C=ROOT/'results'; CACHE=C/'cancer_scaled/cache'; V63=C/'cancer_v6_3_extras'; V64=C/'cancer_v6_4_extras'; V65=C/'cancer_v6_5_extras'; OUT=C/'cancer_v6_6_extras'
BCC=lambda o:((o['disease'].astype(str)=='basal cell carcinoma')&(o['cell_type'].astype(str).isin(['basal cell','basal cell of epidermis','suprabasal keratinocyte']))).values
DATA=[('HCC',C/'cancer_sun2021_hcc/sun2021.h5ad',lambda o:(o['celltype']=='Hepatocyte').values),
 ('BRCA',C/'cancer_validation/wu2021_brca.h5ad',None),('Glioblastoma',CACHE/'Glioblastoma_malignant.h5ad',None),
 ('OvarianSerous',CACHE/'OvarianSerous_malignant.h5ad',None),('RCC',CACHE/'RCC_clearcell_malignant.h5ad',None),
 ('Melanoma',CACHE/'Melanoma_malignant.h5ad',None),('Prostate',CACHE/'Prostate_malignant.h5ad',None),
 ('HNSCC',V63/'HNSCC_malignant.h5ad',None),('Neuroblastoma',V64/'Neuroblastoma_malignant.h5ad',None),
 ('B-ALL',V64/'B-ALL_malignant.h5ad',None),('DLBCL',V64/'DLBCL_malignant.h5ad',None),
 ('SqLung',V64/'SqLung_malignant.h5ad',None),('TNBC',V64/'TNBC_malignant.h5ad',None),
 ('PDAC',V65/'PDAC_Peng2019_malignant.h5ad',None),
 ('NSCLC',ROOT/'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',lambda o:(o['cell_type'].astype(str)=='malignant cell').values),
 ('AML',OUT/'vanGalen_AML_malignant.h5ad',None),('MM',OUT/'Ledergor_MM_malignant.h5ad',None),
 ('BCC',ROOT/'paper6_immune_3layer/data/bcc_full.h5ad',BCC),
 ('Colon',C/'cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',lambda o:(o['ClusterMidway'].astype(str)=='EpiT').values),
 ('SCLC',OUT/'SCLC_malignant.h5ad',None),('Gastric',OUT/'Gastric_malignant.h5ad',None),
 ('GermCell',OUT/'GermCell_malignant.h5ad',None),('Osteosarcoma',OUT/'Osteosarcoma_malignant.h5ad',None),
 ('Ewing',OUT/'Ewing_malignant.h5ad',None),('Seminoma',OUT/'Seminoma_malignant.h5ad',None)]
def nv(a):
    v0=str(a.var.index[0])
    if 'feature_name' in a.var.columns and (v0.startswith('ENSG') or not any(c.isalpha() for c in v0)): a.var_names=a.var['feature_name'].astype(str).values
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    base=pd.Series([re.sub(r'-\d+$','',str(v)) for v in a.var_names]); k=~base.duplicated().values
    a=a[:,k].copy(); a.var_names=base[k].values; return a
rng=np.random.default_rng(0); Zblocks=[]; ng=[]; ids=[]
for name,p,pf in DATA:
    try:
        a=sc.read_h5ad(p,backed='r'); a=(a[pf(a.obs)].to_memory() if pf else a.to_memory()); a=nv(a)
        g=np.asarray((a.X>0).sum(1)).ravel(); a=a[g>=1000].copy(); g=g[g>=1000]
        if a.n_obs>CAP:
            sel=rng.choice(a.n_obs,CAP,replace=False); a=a[sel].copy(); g=g[sel]
        if a.n_obs<50: continue
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0)
        Zblocks.append(((S.values-mu)/sd).astype(np.float32)); ng+=list(g); ids+=[name]*a.n_obs; stage(f"{name}: {a.n_obs}")
        del a,S; gc.collect()
    except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
Z=np.vstack(Zblocks); ng=np.array(ng); cid=np.array(ids); del Zblocks; gc.collect()
print(f"\nscored {Z.shape} cells, {len(set(cid))} cancers\n",flush=True)

def project(Zc):
    pc=Zc@load; return pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
def report(pcn, refPn, tag, chunk=20000):
    n=pcn.shape[0]; top1=np.empty(n,int); soft=0.0; core_top3=0; core_top5=0; core_rank=np.empty(n,int)
    core_set=set(core_idx.tolist()); tau=0.05
    for i in range(0,n,chunk):
        sim=pcn[i:i+chunk]@refPn.T; sim[:,EXCLUDE]=-1e9
        t1=np.argmax(sim,axis=1); top1[i:i+chunk]=t1
        order=np.argsort(-sim,axis=1)
        for r,row in enumerate(order):
            rk=next((k+1 for k,c in enumerate(row) if c in core_set),154); core_rank[i+r]=rk
        core_top3+=int(np.isin(order[:,:3],core_idx).any(1).sum()); core_top5+=int(np.isin(order[:,:5],core_idx).any(1).sum())
        w=np.exp((sim-sim.max(1,keepdims=True))/tau); w/=w.sum(1,keepdims=True); soft+=float(w[:,core_idx].sum())
    hard=100*np.isin(top1,core_idx).mean()
    print(f"[{tag}] HARD top-1 core={hard:.1f}% | SOFT(τ={tau}) core mass={100*soft/n:.1f}% | core in top-3={100*core_top3/n:.1f}% top-5={100*core_top5/n:.1f}% | median core rank={int(np.median(core_rank))}",flush=True)
    return top1, hard

# --- baseline (hard/soft/top-k/rank) ---
pcn=project(Z); top1,hard=report(pcn,refPn,'REAL projection')
# #3c: are the winning centroids the highest-norm? correlation of assignment with centroid norm
picks=np.bincount(top1,minlength=154).astype(float)
r_norm=np.corrcoef(picks[[i for i in range(154) if i not in set(EXCLUDE.tolist())]],
                   refnorm[[i for i in range(154) if i not in set(EXCLUDE.tolist())]])[0,1]
print(f"[#3c] corr(cells assigned to a centroid, centroid 9-D norm) = {r_norm:+.2f}  (near 0 => high-norm centroids are NOT just attracting everything)",flush=True)
# #4a: does core assignment track cell complexity (n_genes)?
isc=np.isin(top1,core_idx).astype(float)
print(f"[#4a] n_genes: core-cells median {np.median(ng[isc==1]):.0f} vs non-core {np.median(ng[isc==0]):.0f} | point-biserial r(core, n_genes)={np.corrcoef(isc,ng)[0,1]:+.2f}  (>0 => core is HIGH-complexity, i.e. NOT the low-complexity erythrocyte artifact)",flush=True)
# #4b: magnitude-removed — apply the SAME transform to BOTH reference and cancer (fair test),
#      rebuild the reference centroids under the transform, then project+assign.
Zref=(Rm.values-mu)/sd
def transform_both(fn, tag):
    zc=fn(Z); zr=fn(Zref)
    rpn=zr@load; rpn=rpn/(np.linalg.norm(rpn,axis=1,keepdims=True)+1e-12)
    report(project(zc), rpn, tag)
transform_both(lambda x: x-x.mean(1,keepdims=True), 'MAGNITUDE-REMOVED (both ref+cancer within-cell centered)')
def rankt(x):
    r=np.argsort(np.argsort(x,axis=1),axis=1).astype(np.float32); return r-r.mean(1,keepdims=True)
transform_both(rankt, 'RANK-TRANSFORMED (both ref+cancer, shape-only)')
print("\nDONE",flush=True)
