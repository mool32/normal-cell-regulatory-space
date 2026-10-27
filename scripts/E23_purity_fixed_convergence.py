"""Round-4 (1) — the purity audit (Table S-P) FLAGGED three contaminated cohorts but never fixed them.
Contamination with normal-adjacent cells pulls a malignant arm toward normal, so it UNDERSTATES convergence:
cleaning should STRENGTHEN the (currently modest, P~0.04) shape-only convergence. Here we fix what the labels
allow and recompute everything:
   HCC   celltype=='Hepatocyte'                -> + site=='Tumor'          (site: Tumor/Normal/PVTT/Lymph)
   BRCA  pf=None (ALL 100k cells of the atlas) -> celltype_major=='Cancer Epithelial'
   BCC   (no malignant annotation available; per-sample disease label only) -> left as is, reported honestly
Recomputes: per-cohort core fraction, the 25-centroid contraction (raw + shape-only), and the shape-only
contraction against the same matched-averaging null used in E14."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, re, gc, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']; CAP=4000
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
core_idx=np.array([refcells.index(c) for c in CORE])
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
demag=lambda X:X-X.mean(1,keepdims=True); ident=lambda X:X
def refset(tf):
    R=tf(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; return mu,sd,rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
REF={'RAW':refset(ident),'SHAPE':refset(demag)}; TF={'RAW':ident,'SHAPE':demag}
def project(S,m):
    mu,sd,rpn=REF[m]; pc=((TF[m](S)-mu)/sd)@load; return pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)

def score_file(path,maskfn,tag):
    b=sc.read_h5ad(path,backed='r')
    m=np.asarray(maskfn(b.obs)) if maskfn else np.ones(b.n_obs,bool)
    rng=np.random.default_rng(0); idx=np.where(m)[0]
    if len(idx)>CAP:
        sel=rng.choice(idx,CAP,replace=False); bm=np.zeros(b.n_obs,bool); bm[sel]=True
    else: bm=m
    a=b[bm].to_memory(); del b
    if 'feature_name' in a.var.columns and str(a.var_names[0]).startswith('ENSG'): a.var_names=a.var['feature_name'].astype(str).values
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    g=np.asarray((a.X>0).sum(1)).ravel(); a=a[g>=1000].copy()
    if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(np.float32)
    print(f"  {tag}: {len(S)} cells (of {int(m.sum())} matching)",flush=True); del a; gc.collect(); return S

# ---- rescore the fixable cohorts, clean ----
CLEAN={
 'HCC': ('results/cancer_sun2021_hcc/sun2021.h5ad',
         lambda o:((o['celltype'].astype(str)=='Hepatocyte')&(o['site'].astype(str)=='Tumor')).values),
 'BRCA':('results/cancer_validation/wu2021_brca.h5ad',
         lambda o:(o['celltype_major'].astype(str)=='Cancer Epithelial').values),
}
print("rescoring cleaned cohorts:",flush=True)
newS={k:score_file(ROOT/p,f,k) for k,(p,f) in CLEAN.items()}

# ---- combine with the cached S for the other cohorts ----
d=np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True); S0=d['S'].astype(float); cid0=d['cid']
uc=sorted(set(cid0.tolist()))
def build(which):   # 'orig' or 'clean'
    Ss=[]; ids=[]
    for c in uc:
        if which=='clean' and c in newS: X=newS[c].astype(float)
        else: X=S0[cid0==c]
        Ss.append(X); ids+=[c]*len(X)
    return np.vstack(Ss),np.array(ids)
def summarize(S,cid,label):
    out={}
    for m in ['RAW','SHAPE']:
        P=project(S,m); rpn=REF[m][2]
        sim=P@rpn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
        cf=100*np.isin(t,core_idx).mean()
        C=np.array([P[cid==c].mean(0) for c in uc]); C=C/np.linalg.norm(C,axis=1,keepdims=True)
        s=C@C.T; iu=np.triu_indices(len(C),1); contr=float((1-s[iu]).mean())
        out[m]=(cf,contr)
        percore={c:100*np.isin(t[cid==c],core_idx).mean() for c in uc}
        out[m+'_per']=percore
    print(f"\n[{label}] cells={len(S)}")
    for m in ['RAW','SHAPE']: print(f"   {m:6s}: core fraction {out[m][0]:5.1f}% | 25-centroid contraction {out[m][1]:.3f}")
    return out
So,co=build('orig'); Sc,cc=build('clean')
o=summarize(So,co,'ORIGINAL (as published)'); n=summarize(Sc,cc,'PURITY-FIXED (HCC tumor-site, BRCA cancer-epithelial)')
print("\nper-cohort core fraction change (fixed cohorts):")
for c in CLEAN: print(f"   {c}: RAW {o['RAW_per'][c]:.1f}% -> {n['RAW_per'][c]:.1f}% | SHAPE {o['SHAPE_per'][c]:.1f}% -> {n['SHAPE_per'][c]:.1f}%")

# ---- shape-only contraction vs the E14 matched-averaging null ----
pl=np.load(DD/'normalpool_local.npz',allow_pickle=True); Sn=pl['S'].astype(float)
sizes=[int((cc==c).sum()) for c in uc]
Pn=project(Sn,'SHAPE'); G=Pn@Pn.T; rng=np.random.default_rng(0); nn=[]
for _ in range(300):
    seeds=rng.choice(len(Sn),len(uc),replace=False); C=[]
    for s,k in zip(seeds,sizes):
        nb=np.argpartition(-G[s],min(k,len(Sn)-1))[:min(k,len(Sn))]
        v=Pn[nb].mean(0); C.append(v/np.linalg.norm(v))
    C=np.array(C); s_=C@C.T; iu=np.triu_indices(len(C),1); nn.append(float((1-s_[iu]).mean()))
nn=np.array(nn)
print(f"\nSHAPE-ONLY contraction vs matched-averaging null (n=300):")
print(f"   original  {o['SHAPE'][1]:.3f}  -> percentile {100*(nn<o['SHAPE'][1]).mean():.1f}%")
print(f"   purity-fixed {n['SHAPE'][1]:.3f} -> percentile {100*(nn<n['SHAPE'][1]).mean():.1f}%   (null mean {nn.mean():.3f})")
print(f"\n[VERDICT] cleaning {'STRENGTHENS' if n['SHAPE'][1]<o['SHAPE'][1] else 'does not strengthen'} the shape-only convergence "
      f"({o['SHAPE'][1]:.3f} -> {n['SHAPE'][1]:.3f}); core fraction {o['SHAPE'][0]:.1f}% -> {n['SHAPE'][0]:.1f}%")
pd.DataFrame({'cohort':uc,'orig_RAW':[o['RAW_per'][c] for c in uc],'clean_RAW':[n['RAW_per'][c] for c in uc],
              'orig_SHAPE':[o['SHAPE_per'][c] for c in uc],'clean_SHAPE':[n['SHAPE_per'][c] for c in uc]}).to_csv(DD/'E23_purity_fixed.csv',index=False)
print("DONE")
