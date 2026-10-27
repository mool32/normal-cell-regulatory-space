"""E2(e) — FULL-FIDELITY 44-module robustness (Gudkov G6b, ranking-level).
Unlike E2d (reduced gene-pool harness, 23 cancers), this reproduces the PAPER'S pipeline:
 - the exact 25 cancers with their real malignant-cell selection,
 - canonical scoring: pct.score_modules(gene_set='core') per cell (as Rm itself is built),
 - the canonical eigenspace RECIPE: eigh(corr(z-scored modules)) — pct.rebuild() — 9 PCs,
then perturbs the module set and re-derives BOTH rankings the core rests on:
 (i) recurrence  — # of the 25 cancers in which a cell is a top-5 top-1 target (Fig S12),
 (ii) cosine-to-attractor — rank by cosine to the mean of the 25 malignant centroids (Fig 6C/7A),
plus the core fraction. GATE: baseline (full 44) must reproduce cyto/mega at the top + core ~21%.
Perturbations: leave-5-out (44->39, N=200) and leave-1-out (44). Cap 6000 cells/cancer for stable centroids.
"""
import scanpy as sc, anndata as ad, json, numpy as np, pandas as pd, sys, re, gc, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import binomtest
ROOT=Path(__file__).resolve().parents[1]; OUTD=ROOT/'results'
def stage(m): print(f"[stage] {m}",flush=True)
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
mods=list(pct.list_modules()); NPC=9; CAP=6000
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]           # 154 x 44 canonical
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
RM=Rm.values.astype(float)

def eigenspace(cols):                # canonical rebuild() recipe on the module subset
    X=RM[:,cols]; mean=X.mean(0); std=X.std(0); std[std<1e-10]=1.0
    Z=(X-mean)/std; corr=np.corrcoef(Z.T); ev,vec=np.linalg.eigh(corr)
    idx=np.argsort(ev)[::-1]; L=vec[:,idx[:NPC]]
    refP=Z@L; refPn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
    return mean,std,L,refPn,int((ev>1).sum())

# ---- load + canonically score the 25 cancers on the 44 modules (once) ----
C=OUTD; CACHE=C/'cancer_scaled/cache'; V63=C/'cancer_v6_3_extras'; V64=C/'cancer_v6_4_extras'; V65=C/'cancer_v6_5_extras'; OUT=C/'cancer_v6_6_extras'
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
rng=np.random.default_rng(0); blocks=[]; ids=[]
for name,p,pf in DATA:
    try:
        a=sc.read_h5ad(p,backed='r'); a=(a[pf(a.obs)].to_memory() if pf else a.to_memory()); a=nv(a)  # backed: materialize only selected malignant cells (memory-safe)
        ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=1000].copy()
        if a.n_obs>CAP: a=a[rng.choice(a.n_obs,CAP,replace=False)].copy()
        if a.n_obs<50: print(f"  !!{name}: only {a.n_obs} cells",flush=True); continue
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0)
        blocks.append(S.values.astype(np.float32)); ids+=[name]*a.n_obs; stage(f"{name}: {a.n_obs} cells")
        del a,S; gc.collect()
    except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
Sall=np.vstack(blocks); cancer_ids=np.array(ids); del blocks; gc.collect()
uc=sorted(set(cancer_ids)); print(f"\nscored matrix {Sall.shape} | {len(uc)} cancers: {uc}\n",flush=True)

def metrics(cols):
    mean,std,L,refPn,nk=eigenspace(cols)
    z=(Sall[:,cols]-mean)/std; pc=z@L; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; top1=np.argmax(sim,axis=1)
    corefrac=float(np.isin(top1,core_idx).mean())
    recur=np.zeros(154); cents=[]
    for cid in uc:
        m=cancer_ids==cid; t=top1[m]; recur[np.argsort(-np.bincount(t,minlength=154))[:5]]+=1
        cents.append(pcn[m].mean(0))
    recur[EXCLUDE]=-1; ro=np.argsort(-recur); rrank={refcells[ro[i]]:i+1 for i in range(154)}
    CC=np.array(cents); CC=CC/(np.linalg.norm(CC,axis=1,keepdims=True)+1e-12); attr=CC.mean(0); attr/=np.linalg.norm(attr)+1e-12
    cosF=refPn@attr; cosF[EXCLUDE]=-1e9; ao=np.argsort(-cosF); arank={refcells[ao[i]]:i+1 for i in range(154)}
    return dict(nk=nk, corefrac=corefrac, enrich=corefrac/(len(core_idx)/154),
                cyto_rec=rrank['cytotrophoblasts'], mega_rec=rrank['megakaryocytes'], rec_top2=[refcells[ro[0]],refcells[ro[1]]],
                cyto_attr=arank['cytotrophoblasts'], mega_attr=arank['megakaryocytes'], attr_top3=[refcells[ao[i]] for i in range(3)])

# ---- GATE: baseline full 44 ----
allcols=list(range(len(mods))); R0=metrics(allcols)
print("="*70,flush=True)
print(f"[GATE baseline 44 modules] Kaiser PCs={R0['nk']}",flush=True)
print(f"  RECURRENCE ranking : cyto rank {R0['cyto_rec']}, mega rank {R0['mega_rec']}  | top-2={R0['rec_top2']}",flush=True)
print(f"  COSINE-to-attractor: cyto rank {R0['cyto_attr']}, mega rank {R0['mega_attr']}  | top-3={R0['attr_top3']}",flush=True)
print(f"  core fraction={R0['corefrac']*100:.1f}%  (enrich {R0['enrich']:.1f}x)",flush=True)
gate_ok = (R0['cyto_attr']<=3 and R0['mega_attr']<=4) or (R0['cyto_rec']<=2 and R0['mega_rec']<=3)
print(f"  GATE {'PASSED — baseline reproduces the core ranking' if gate_ok else 'WEAK — baseline does not cleanly reproduce (interpret with care)'}",flush=True)
print("="*70+"\n",flush=True)

def summ(rows,label):
    d=pd.DataFrame(rows); print(f"=== {label} (N={len(d)}) ===",flush=True)
    for nm,c in [('cyto','cyto_attr'),('mega','mega_attr')]:
        print(f"  {nm} (cosine-to-attractor): median rank {int(d[c].median())} | top-2 {100*(d[c]<=2).mean():.0f}% | top-3 {100*(d[c]<=3).mean():.0f}% | top-5 {100*(d[c]<=5).mean():.0f}%  (real {R0[c]})",flush=True)
    for nm,c in [('cyto','cyto_rec'),('mega','mega_rec')]:
        print(f"  {nm} (recurrence)        : median rank {int(d[c].median())} | top-2 {100*(d[c]<=2).mean():.0f}% | top-5 {100*(d[c]<=5).mean():.0f}%  (real {R0[c]})",flush=True)
    print(f"  BOTH cyto&mega top-3 (attractor): {100*((d.cyto_attr<=3)&(d.mega_attr<=3)).mean():.0f}%",flush=True)
    print(f"  core_fraction: {100*d.corefrac.mean():.1f} ± {100*d.corefrac.std():.1f}%  min {100*d.corefrac.min():.1f}%  (real {100*R0['corefrac']:.1f}%)",flush=True)
    print(f"  Kaiser PCs across draws: {int(d.nk.min())}–{int(d.nk.max())}\n",flush=True)
    return d

# ---- leave-5-out ----
N5=200; rows5=[]
for b in range(N5):
    rg=np.random.default_rng(7000+b); drop=set(rg.choice(len(mods),5,replace=False))
    m=metrics([j for j in range(len(mods)) if j not in drop]); m['dropped']=';'.join(mods[j] for j in sorted(drop)); rows5.append(m)
    if (b+1)%50==0: stage(f"leave-5-out {b+1}/{N5}")
d5=summ(rows5,"LEAVE-5-OUT (drop 5 random real modules, 44->39)")
# ---- leave-1-out ----
rows1=[]
for j in range(len(mods)):
    m=metrics([k for k in range(len(mods)) if k!=j]); m['dropped']=mods[j]; rows1.append(m)
d1=summ(rows1,"LEAVE-1-OUT (drop each single module)")
worst=d1.reindex(d1.corefrac.sort_values().index).head(3)
print("  3 lowest-core single drops:",[(r.dropped,round(100*r.corefrac,1),r.cyto_attr,r.mega_attr) for _,r in worst.iterrows()],flush=True)

for d,nm in [(d5,'leave5'),(d1,'leave1')]: d.to_csv(OUTD/f'E2e_fullfid_{nm}.csv',index=False)
pd.DataFrame([R0]).to_csv(OUTD/'E2e_fullfid_baseline.csv',index=False)
print("\nsaved results/E2e_fullfid_{baseline,leave5,leave1}.csv",flush=True)
