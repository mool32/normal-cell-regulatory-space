"""Referee control #2 — is the central result (convergence + compactness) calibration-specific?
Score the 25 cancers ONCE (raw mean_raw module scores; calibration-independent), then under EACH
reference calibration (pseudobulk vs single_cell_scaled) build that calibration's own eigenspace
(correlation-PCA, Kaiser 9) and reference centroids, project the cancers, and measure:
  - core fraction (top-1 in {cyto, mega, migrating}) — exclude erythrocyte/platelet;
  - malignant COMPACTNESS = mean pairwise cosine distance among the 25 malignant centroids,
    and its percentile vs 2,000 random 25-cell draws of normal cells (low = tighter than random).
Reuses the 25-cancer loader (backed+mask; canonical score_modules)."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, re, gc, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def stage(m): print(f"[stage] {m}",flush=True)
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']; NPC=9; CAP=4000
mods=list(pct.list_modules())
# ---- score 25 cancers once -> raw module scores S (calibration-independent) ----
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
rng=np.random.default_rng(0); Sb=[]; ids=[]
for name,p,pf in DATA:
    try:
        a=sc.read_h5ad(p,backed='r'); a=(a[pf(a.obs)].to_memory() if pf else a.to_memory()); a=nv(a)
        g=np.asarray((a.X>0).sum(1)).ravel(); a=a[g>=1000].copy()
        if a.n_obs>CAP: a=a[rng.choice(a.n_obs,CAP,replace=False)].copy()
        if a.n_obs<50: continue
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0)
        Sb.append(S.values.astype(np.float32)); ids+=[name]*a.n_obs; stage(f"{name}: {a.n_obs}")
        del a,S; gc.collect()
    except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
S=np.vstack(Sb); cid=np.array(ids); uc=sorted(set(cid)); del Sb; gc.collect()
print(f"\nscored {S.shape}, {len(uc)} cancers\n",flush=True)

def eigen(Rm):
    X=Rm.values.astype(float); mean=X.mean(0); std=X.std(0); std[std<1e-10]=1.0
    Z=(X-mean)/std; corr=np.corrcoef(Z.T); ev,vec=np.linalg.eigh(corr); idx=np.argsort(ev)[::-1]
    L=vec[:,idx[:NPC]]; refP=Z@L; refPn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
    return mean,std,L,refPn,int((ev>1).sum())
def mpd(X): s=X@X.T; iu=np.triu_indices(len(X),1); return float((1-s[iu]).mean())

for mode in ['single_cell_scaled','pseudobulk']:
    Rm=load_hpa_perceptivity(mode=mode)['R'][mods]; refcells=list(Rm.index)
    core_idx=np.array([refcells.index(c) for c in CORE]); EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
    mean,std,L,refPn,nk=eigen(Rm)
    z=(S-mean)/std; pc=z@L; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; top1=np.argmax(sim,axis=1)
    corefrac=100*np.isin(top1,core_idx).mean()
    cents=np.array([pcn[cid==c].mean(0) for c in uc]); cents=cents/(np.linalg.norm(cents,axis=1,keepdims=True)+1e-12)
    mal=mpd(cents)
    keep=np.array([i for i in range(len(refcells)) if i not in set(EXCLUDE.tolist())])
    rg2=np.random.default_rng(1); rnd=[]
    for _ in range(2000):
        d=refPn[rg2.choice(keep,len(uc),replace=False)]; d=d/(np.linalg.norm(d,axis=1,keepdims=True)+1e-12); rnd.append(mpd(d))
    rnd=np.array(rnd); pct_tile=(rnd<=mal).mean()
    print(f"[{mode:18s}] Kaiser PCs={nk} | core fraction={corefrac:.1f}% (enrich {corefrac/(3/152):.1f}×) | malignant compactness(mpd)={mal:.3f} vs random mean {rnd.mean():.3f} | percentile={pct_tile:.2f} (0=tighter than ALL random, ~0.5=diffuse) | contraction random/malig={rnd.mean()/mal:.2f}×",flush=True)
print("\nDONE",flush=True)
