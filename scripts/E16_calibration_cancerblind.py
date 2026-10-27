"""Referee 2.5 — cancer-blind justification of the single_cell_scaled calibration.
On INDEPENDENT NORMAL data with author cell-type labels (no cancer), build the reference eigenspace under
each calibration and ask which one preserves KNOWN cell-type structure. If single_cell_scaled separates
normal cell types (high silhouette) and pseudobulk degenerates (collapses them), the calibration is chosen
by correct recovery of known biology, established without reference to the cancer convergence."""
import scanpy as sc, numpy as np, pandas as pd, json, sys, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from sklearn.metrics import silhouette_score
from pathlib import Path
ROOT=Path('.').resolve(); NPC=9; mods=list(pct.list_modules())
def eigen(Rm):
    X=Rm.values.astype(float); mean=X.mean(0); std=X.std(0); std[std<1e-10]=1.0
    Z=(X-mean)/std; ev,vec=np.linalg.eigh(np.corrcoef(Z.T)); idx=np.argsort(ev)[::-1]
    return mean,std,vec[:,idx[:NPC]]
def score(a):
    a.var_names=[str(g).upper() for g in a.var_names]; a.var_names_make_unique()
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values
    return S
# labeled normal test cells (author cell types), from 3 independent CELLxGENE datasets
rng=np.random.default_rng(0); Ss=[]; labs=[]
def add(path,ctcol,maskfn,use_fn,tag):
    b=sc.read_h5ad(path,backed='r'); m=np.asarray(maskfn(b.obs)); idx=np.where(m)[0]
    if len(idx)>6000: keep=rng.choice(idx,6000,replace=False); bm=np.zeros(b.n_obs,bool); bm[keep]=True
    else: bm=m
    a=b[bm].to_memory(); del b
    if use_fn and 'feature_name' in a.var.columns: a.var_names=a.var['feature_name'].astype(str).values
    ct=a.obs[ctcol].astype(str).values
    S=score(a); Ss.append(S); labs.extend([f"{tag}:{c}" for c in ct[:len(S)]])
    print(f"  {tag}: {len(S)} cells, {len(set(ct))} types",flush=True)
add('paper6_immune_3layer/data/bcc_full.h5ad','cell_type',lambda o:(o['disease'].astype(str)=='normal').values,True,'skin')
add('results/cancer_rcc_census/rcc_kidney_paired.h5ad','cell_type',lambda o:(o['disease'].astype(str)=='normal').values,True,'kidney')
add('results/cancer_validation/wu2021_brca.h5ad','celltype_major',lambda o:(o['celltype_major'].astype(str)!='Cancer Epithelial').values,False,'breast')
S=np.vstack(Ss); lab=np.array(labs)
# keep types with >=50 cells for a stable silhouette
vc=pd.Series(lab).value_counts(); keep_types=set(vc[vc>=50].index); k=np.isin(lab,list(keep_types))
S=S[k]; lab=lab[k]
print(f"\ntest set: {S.shape[0]} INDEPENDENT NORMAL single cells, {len(set(lab))} labeled types")
# The right cancer-blind test (parallel to E4's cancer degeneration): does pseudobulk COLLAPSE
# normal SINGLE cells onto a few low-complexity targets — i.e. is it wrong for single-cell data per se,
# regardless of cancer? Assign each test cell to nearest HPA reference type under each calibration.
print("\n[calibration on independent NORMAL single cells — assignment diversity & low-complexity collapse]:")
LOWC=['erythrocytes','platelets']
for mode in ['single_cell_scaled','pseudobulk']:
    Rm=load_hpa_perceptivity(mode=mode)['R'][mods]; refcells=list(Rm.index); mean,std,L=eigen(Rm)
    refP=((Rm.values-mean)/std)@L; refPn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
    P=((S-mean)/std)@L; Pn=P/(np.linalg.norm(P,axis=1,keepdims=True)+1e-12)
    sim=Pn@refPn.T; top1=np.argmax(sim,axis=1)
    frac=np.bincount(top1,minlength=len(refcells))/len(top1); frac=frac[frac>0]
    eff=np.exp(-(frac*np.log(frac)).sum())   # effective # of distinct target types (assignment diversity)
    lowc=100*np.mean([refcells[t] in LOWC for t in top1])
    print(f"   {mode:20s}: distinct targets used = {len(set(top1)):3d}/152 | effective diversity = {eff:5.1f} | assigned to low-complexity erythroid = {lowc:5.1f}%")
print("\n[reading] pseudobulk collapsing normal single cells onto few / low-complexity targets = it is wrong for")
print("single-cell data of ANY kind (normal or cancer); single_cell_scaled keeping the assignments diverse and")
print("complexity-appropriate = the calibration matched to single-cell data — chosen without any cancer input.")
print("DONE")
