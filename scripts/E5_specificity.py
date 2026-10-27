"""Referee round-2: turn the magnitude concession into a specificity STRENGTH.
(A) SHAPE-ONLY ranking: after removing magnitude (rank-transform BOTH ref+cancer), does the convergence
    still NAME the cytotrophoblast/megakaryocyte, or smear across high-output cells? (recurrence ranking)
(B) HIGH-OUTPUT specificity: convergence RECEIVED by professional high-output secretory normals
    (hepatocytes, plasma cells, acinar, breast-secretory, gastric chief, parietal). If it were pure magnitude
    they'd attract too; if they get ~0, the address is specific, not "any bright cell".
(C) TOTAL-COUNTS control: is core assignment explained by library size (total counts), not just n_genes?
(E) ORTHOGONAL trophoblast markers OUTSIDE the 206 module genes (HLA-G, MMP2/9, EVT/syncytin, PLAC1, PEG10...):
    do core-assigned cancer cells express them more than non-core? (independent 'trophoblast-ness' check)
Reuses the 25-cancer loader."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, re, gc, warnings
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
MODGENES=set(g.upper() for m in mods for g in pct.get_genes(m))
MARKERS=[g for g in ['HLA-G','MMP2','MMP9','PLAC1','PEG10','ERVW-1','ERVFRD-1','ITGA5','ITGA1','HTRA1',
    'PAPPA','KRT7','TFAP2C','CGA','CGB','CGB5','CSH1','ADAM12','PRG2','GATA3','CDKN1C'] if g.upper() not in MODGENES]
print(f"orthogonal trophoblast markers used (all OUTSIDE the {len(MODGENES)} module genes): {MARKERS}",flush=True)

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
rng=np.random.default_rng(0); Sb=[]; ng=[]; tc=[]; ids=[]; mk=[]
for name,p,pf in DATA:
    try:
        a=sc.read_h5ad(p,backed='r'); a=(a[pf(a.obs)].to_memory() if pf else a.to_memory()); a=nv(a)
        g=np.asarray((a.X>0).sum(1)).ravel(); a=a[g>=1000].copy(); g=g[g>=1000]
        if a.n_obs>CAP:
            sel=rng.choice(a.n_obs,CAP,replace=False); a=a[sel].copy(); g=g[sel]
        if a.n_obs<50: continue
        raw_tot=np.asarray(a.X.sum(1)).ravel()
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0)
        present=[m2 for m2 in MARKERS if m2 in set(a.var_names)]
        M=np.zeros((a.n_obs,len(MARKERS)),np.float32); pos={g2:i for i,g2 in enumerate(MARKERS)}
        if present:
            sub=a[:,present].X; sub=np.asarray(sub.todense()) if hasattr(sub,'todense') else np.asarray(sub)
            for k2,g2 in enumerate(present): M[:,pos[g2]]=sub[:,k2]
        Sb.append(S.values.astype(np.float32)); ng+=list(g); tc+=list(raw_tot); ids+=[name]*a.n_obs; mk.append(M)
        stage(f"{name}: {a.n_obs}"); del a,S; gc.collect()
    except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
S=np.vstack(Sb); ng=np.array(ng); tc=np.array(tc); cid=np.array(ids); MK=np.vstack(mk); uc=sorted(set(cid)); del Sb,mk; gc.collect()
print(f"\nscored {S.shape}, {len(uc)} cancers\n",flush=True)

def proj(Z): pc=Z@load; return pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
def assign(pcn, rpn):
    sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; return np.argmax(sim,axis=1)
def recurrence(top1):
    rec=np.zeros(154)
    for c in uc:
        t=top1[cid==c]; rec[np.argsort(-np.bincount(t,minlength=154))[:5]]+=1
    rec[EXCLUDE]=-1; o=np.argsort(-rec); return [(refcells[o[i]],int(rec[o[i]])) for i in range(6)]

Z=(S-mu)/sd
# ---- REAL ----
pcnR=proj(Z); top1R=assign(pcnR,refPn); isc=np.isin(top1R,core_idx)
print("[A/REAL] top recurrence targets:", recurrence(top1R),flush=True)
# ---- (A) SHAPE-ONLY: rank-transform BOTH ref & cancer ----
def rankt(x):
    r=np.argsort(np.argsort(x,axis=1),axis=1).astype(np.float32); return r-r.mean(1,keepdims=True)
Zc=rankt(Z); Zrefr=rankt((Rm.values-mu)/sd); rpnr=proj_ref=(Zrefr@load); rpnr=rpnr/(np.linalg.norm(rpnr,axis=1,keepdims=True)+1e-12)
top1S=assign(proj(Zc),rpnr)
print("[A/SHAPE-ONLY] top recurrence targets (magnitude removed):", recurrence(top1S),flush=True)
print(f"   -> core fraction real {100*isc.mean():.1f}% vs shape-only {100*np.isin(top1S,core_idx).mean():.1f}%",flush=True)
# ---- (B) HIGH-OUTPUT specificity: convergence received by professional secretory normals ----
frac=np.bincount(top1R,minlength=154)/len(top1R)*100
HI=['cytotrophoblasts','megakaryocytes','hepatocytes','plasma cells','pancreatic acinar cells','gastric chief cells',
    'salivary gland mucous cells','breast glandular cells','parietal cells','prostatic glandular cells','proximal tubule cells']
print("[B] convergence RECEIVED (% of all cancer cells whose nearest normal it is):",flush=True)
for c in HI:
    if c in refcells: print(f"     {c:34s} {frac[refcells.index(c)]:5.2f}%",flush=True)
# ---- (C) total-counts control ----
print(f"[C] core vs non-core: total_counts median {np.median(tc[isc]):.0f} vs {np.median(tc[~isc]):.0f} | r(core, log10 total_counts)={np.corrcoef(isc.astype(float),np.log10(tc+1))[0,1]:+.2f} | r(core, n_genes)={np.corrcoef(isc.astype(float),ng)[0,1]:+.2f}",flush=True)
# ---- (E) orthogonal trophoblast markers OUTSIDE module genes ----
print("[E] orthogonal trophoblast/EVT markers — mean log-expr in core-assigned vs non-core cancer cells:",flush=True)
for i,g2 in enumerate(MARKERS):
    a_=MK[isc,i]; b_=MK[~isc,i]
    if a_.mean()+b_.mean()>0.01: print(f"     {g2:10s} core {a_.mean():.3f} vs non-core {b_.mean():.3f}  ({'HIGHER in core' if a_.mean()>b_.mean()*1.3 else ('~same' if abs(a_.mean()-b_.mean())<0.3*max(a_.mean(),1e-9) else 'lower')})",flush=True)
print("\nDONE",flush=True)
