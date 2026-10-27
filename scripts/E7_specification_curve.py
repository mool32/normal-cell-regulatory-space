"""Referee round-2 ask #6 (multiverse): specification curve of the core convergence.
One honest picture assembling every analytic branch we have run separately, so the reader
sees the whole distribution of the core-fraction estimate rather than one number.

Axes (fully crossed):
  - calibration:        {single_cell_scaled, pseudobulk}     (eigenspace rebuilt per calibration)
  - magnitude transform:{raw, within-cell-center, per-cell-rank}   (applied to BOTH reference & cancer)
  - core definition:    {3-cell cyto+mega+migrating, 2-cell cyto+mega, cyto-only}
  - erythroid handling: {excluded (default), included}
Primary outcome = hard top-1 core fraction (comparable across all specs). Reported vs each spec's
own chance baseline (n_core / n_candidates). Sensitivity table: assignment rule {hard, soft-softmax,
top-3, top-5} at the default spec.
Loads the 25 cancers ONCE, caches raw module scores S to results/E7_S_raw.npz for instant re-sweeps."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, re, gc, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def stage(m): print(f"[stage] {m}",flush=True)
CORE3=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']; NPC=9; CAP=4000
mods=list(pct.list_modules())
C=ROOT/'results'; CACHE=C/'cancer_scaled/cache'; V63=C/'cancer_v6_3_extras'; V64=C/'cancer_v6_4_extras'; V65=C/'cancer_v6_5_extras'; OUT=C/'cancer_v6_6_extras'
SCACHE=C/'E7_S_raw.npz'

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

if SCACHE.exists():
    d=np.load(SCACHE,allow_pickle=True); S=d['S']; cid=d['cid']; stage(f"loaded cached S {S.shape}")
else:
    rng=np.random.default_rng(0); Sb=[]; ids=[]
    for name,p,pf in DATA:
        try:
            a=sc.read_h5ad(p,backed='r'); a=(a[pf(a.obs)].to_memory() if pf else a.to_memory()); a=nv(a)
            g=np.asarray((a.X>0).sum(1)).ravel(); a=a[g>=1000].copy()
            if a.n_obs>CAP: a=a[rng.choice(a.n_obs,CAP,replace=False)].copy()
            if a.n_obs<50: continue
            if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
            Sc=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0)
            Sb.append(Sc.values.astype(np.float32)); ids+=[name]*a.n_obs; stage(f"{name}: {a.n_obs}")
            del a,Sc; gc.collect()
        except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
    S=np.vstack(Sb); cid=np.array(ids); np.savez(SCACHE,S=S,cid=cid); del Sb; gc.collect()
uc=sorted(set(cid.tolist())); print(f"\nscored {S.shape}, {len(uc)} cancers\n",flush=True)

# ---- magnitude transforms (applied identically to reference & cancer raw module-score matrices) ----
def t_raw(X):    return X
def t_center(X): return X - X.mean(1,keepdims=True)
def t_rank(X):
    r=np.argsort(np.argsort(X,axis=1),axis=1).astype(np.float64); return r-r.mean(1,keepdims=True)
TRANSFORMS={'raw':t_raw,'center':t_center,'rank':t_rank}

def eigen(Xref):
    mean=Xref.mean(0); std=Xref.std(0); std[std<1e-10]=1.0
    Z=(Xref-mean)/std; corr=np.corrcoef(Z.T); ev,vec=np.linalg.eigh(corr); idx=np.argsort(ev)[::-1]
    L=vec[:,idx[:NPC]]; refP=Z@L; refPn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
    return mean,std,L,refPn

rows=[]; default_key=None; default_bundle=None
for calib in ['single_cell_scaled','pseudobulk']:
    Rm=load_hpa_perceptivity(mode=calib)['R'][mods]; refcells=list(Rm.index); Xref0=Rm.values.astype(float)
    idx_map={c:refcells.index(c) for c in refcells}
    for tname,tfn in TRANSFORMS.items():
        Xref=tfn(Xref0); mean,std,L,refPn=eigen(Xref)
        Sc=tfn(S.astype(float)); z=(Sc-mean)/std; pc=z@L; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        sim_full=pcn@refPn.T  # cells x refcells
        for eryth in ['excl','incl']:
            EXC=[idx_map[c] for c in ['erythrocytes','platelets'] if c in idx_map] if eryth=='excl' else []
            sim=sim_full.copy()
            if EXC: sim[:,EXC]=-1e9
            ncand=len(refcells)-len(EXC)
            top1=np.argmax(sim,axis=1)
            for cname,coreset in [('3cell',CORE3),('2cell',['megakaryocytes','cytotrophoblasts']),('cyto',['cytotrophoblasts'])]:
                cidx=[idx_map[c] for c in coreset]
                corefrac=100*np.isin(top1,cidx).mean()
                chance=100*len(cidx)/ncand
                rows.append(dict(calib=calib,magnitude=tname,eryth=eryth,coredef=cname,assign='hard',
                                 core_pct=corefrac,chance=chance,enrich=corefrac/chance))
        # keep default bundle (single_cell_scaled, raw) for the assignment-rule sensitivity table
        if calib=='single_cell_scaled' and tname=='raw':
            default_bundle=(sim_full.copy(),idx_map,refcells)

df=pd.DataFrame(rows); df.to_csv(C/'E7_spec_curve.csv',index=False)
hard=df[df['assign']=='hard']
print("="*72)
print(f"SPECIFICATION CURVE — {len(hard)} specs (calibration×magnitude×erythroid×coredef), hard top-1")
print(f"  core fraction: median {hard.core_pct.median():.1f}%  IQR [{hard.core_pct.quantile(.25):.1f}, {hard.core_pct.quantile(.75):.1f}]  range [{hard.core_pct.min():.1f}, {hard.core_pct.max():.1f}]")
print(f"  every spec above its own chance baseline: {(hard.core_pct>hard.chance).mean()*100:.0f}%  (min enrichment {hard.enrich.min():.1f}×, median {hard.enrich.median():.1f}×)")
print(f"  worst-case spec: {hard.loc[hard.core_pct.idxmin(),['calib','magnitude','eryth','coredef']].to_dict()} -> {hard.core_pct.min():.1f}% ({hard.enrich.min():.1f}× chance)")
print("\n  by magnitude transform (median core% over the other axes):")
for t,g in hard.groupby('magnitude'): print(f"     {t:8s} median {g.core_pct.median():.1f}%  min-enrich {g.enrich.min():.1f}×")
print("  by calibration:")
for t,g in hard.groupby('calib'): print(f"     {t:20s} median {g.core_pct.median():.1f}%  min-enrich {g.enrich.min():.1f}×")
print("  by core definition:")
for t,g in hard.groupby('coredef'): print(f"     {t:8s} median {g.core_pct.median():.1f}%  min-enrich {g.enrich.min():.1f}×")

# ---- assignment-rule sensitivity at the default spec (single-cell, raw, excl-erythroid, 3-cell) ----
sim_full,idx_map,refcells=default_bundle
EXC=[idx_map[c] for c in ['erythrocytes','platelets'] if c in idx_map]
sim=sim_full.copy(); sim[:,EXC]=-1e9; ncand=len(refcells)-len(EXC)
cidx=np.array([idx_map[c] for c in CORE3])
top1=np.argmax(sim,axis=1)
hardf=100*np.isin(top1,cidx).mean()
# soft: softmax over similarities, mass on core
sm=np.exp((sim-sim.max(1,keepdims=True))/0.05); sm=sm/sm.sum(1,keepdims=True); softf=100*sm[:,cidx].sum(1).mean()
order=np.argsort(-sim,axis=1)
top3=100*np.isin(order[:,:3],cidx).any(1).mean(); top5=100*np.isin(order[:,:5],cidx).any(1).mean()
print("\n  assignment-rule sensitivity @ default (single-cell, raw, 3-cell, excl-erythroid):")
print(f"     hard top-1 {hardf:.1f}% | soft-softmax mass {softf:.1f}% | top-3 membership {top3:.1f}% | top-5 membership {top5:.1f}%  (chance top-1 {100*3/ncand:.1f}%)")
print("\nDONE",flush=True)
