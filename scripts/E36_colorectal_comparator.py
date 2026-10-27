"""Gudkov's comparator question (2026-10-12): cancer arises from a specific epithelial subtype (crypt-base
stem/transit-amplifying cells), not from 'the tissue'. Is the colorectal within-patient result comparing the
right cells? Pelka et al. annotate both arms with the same epithelial subtypes (ClusterFull: cE01-cE03
Stem/TA-like; cE04/05 enterocytes; cE06-08 goblet lineage; cE09 Best4; cE10 tuft; cE11 enteroendocrine).
  A  core occupancy of each NORMAL subtype (does the crypt-base compartment already sit higher?)
  B  per-patient tests: (i) all malignant vs all normal [gate: reproduces the paper];
     (ii) all malignant vs normal Stem/TA only; (iii) malignant Stem/TA vs normal Stem/TA;
     (iv) subtype-to-subtype (Tumor cE0X vs normal cE0X)
  C  whole-transcriptome (HVG, not modules): which normal subtype is each patient's malignant epithelium closest to?
Uses the validation-suite cache of the paper's exact Pelka epithelium (n_genes >= 1000)."""
import anndata as ad, scanpy as sc, json, numpy as np, pandas as pd, sys, warnings, gc
from scipy.stats import wilcoxon, binomtest
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve(); OUT=ROOT/'results/referee2_newdata'
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
CORE3=[refcells.index(c) for c in ['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']]
TROPH=[refcells.index(c) for c in ['cytotrophoblasts','migrating cytotrophoblasts']]
EX=[refcells.index(c) for c in ['erythrocytes','platelets']]
demag=lambda X:X-X.mean(1,keepdims=True); ident=lambda X:X
def refset(tf):
    R=tf(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; return mu,sd,rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
REF={'RAW':(refset(ident),ident),'SHAPE':(refset(demag),demag)}

# subtype labels from the original object via src_row
orig=sc.read_h5ad(ROOT/'results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',backed='r')
full=orig.obs['ClusterFull'].astype(str).values; del orig
C=ad.read_h5ad(ROOT/'perceptome_validation/cache/pelka_epi.h5ad',backed='r')
obs=C.obs.copy(); obs['sub']=full[obs['src_row'].astype(int).values]
obs['subtype']=obs['sub'].str.replace('Tumor ','',regex=False)
obs['arm2']=np.where(obs['sub'].str.startswith('Tumor'),'malig','normal')
print(f"cache {C.shape}; arm(cache) vs arm(ClusterFull) agreement: {np.mean((obs['arm'].astype(str).str.lower().str[:3]=='mal')==(obs['arm2']=='malig')):.4f}",flush=True)
STEM={'cE01 (Stem/TA-like)','cE02 (Stem/TA-like/Immature Goblet)','cE03 (Stem/TA-like prolif)'}

# HVGs from a normal-only subsample
rng=np.random.default_rng(0); nidx=np.where(obs['arm2'].values=='normal')[0]
s=np.sort(rng.choice(nidx,min(15000,len(nidx)),replace=False)); m=np.zeros(C.n_obs,bool); m[s]=True
a=C[m].to_memory()
if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat'); HVG=a.var_names[a.var['highly_variable']].tolist(); del a; gc.collect()
hv_idx=[C.var_names.get_loc(g) for g in HVG]
print(f"HVGs: {len(HVG)}",flush=True)

assign={k:np.zeros(C.n_obs,int) for k in REF}
key=obs['patient'].astype(str)+'|'+obs['arm2']+'|'+obs['subtype']
sums={}; cnt={}
CH=12000
for st in range(0,C.n_obs,CH):
    m=np.zeros(C.n_obs,bool); m[st:st+CH]=True
    a=C[m].to_memory()
    if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
    for k,((mu,sd,rpn),tf) in REF.items():
        pc=((tf(S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        sim=pcn@rpn.T; sim[:,EX]=-1e9; assign[k][st:st+CH]=np.argmax(sim,axis=1)
    X=a.X[:,hv_idx]; X=X.toarray() if hasattr(X,'toarray') else np.asarray(X)
    kk=key.values[st:st+CH]
    for u in np.unique(kk):
        r=kk==u; sums[u]=sums.get(u,0)+X[r].sum(0); cnt[u]=cnt.get(u,0)+int(r.sum())
    del a,S,X; gc.collect(); print(f"  chunk {st//CH+1}/{int(np.ceil(C.n_obs/CH))}",flush=True)

obs['core_RAW']=np.isin(assign['RAW'],CORE3); obs['core_SHAPE']=np.isin(assign['SHAPE'],CORE3)
obs['troph_RAW']=np.isin(assign['RAW'],TROPH); obs['troph_SHAPE']=np.isin(assign['SHAPE'],TROPH)
obs['stem']=obs['subtype'].isin(STEM)
obs[['patient','arm2','subtype','stem','core_RAW','core_SHAPE','troph_RAW','troph_SHAPE']].to_csv(OUT/'E36_cells.csv.gz',index=False)

print("\n"+"="*86); print("A. NORMAL epithelium: core occupancy by subtype (pooled normal cells)")
nm=obs[obs.arm2=='normal']
A=nm.groupby('subtype').agg(n=('patient','size'),core_raw=('core_RAW','mean'),core_shape=('core_SHAPE','mean'),troph_raw=('troph_RAW','mean'))
print((A.assign(core_raw=100*A.core_raw,core_shape=100*A.core_shape,troph_raw=100*A.troph_raw).round(1)).sort_values('core_raw',ascending=False).to_string())
print(f"  normal Stem/TA share of normal arm: {100*nm.stem.mean():.1f}% | malignant Stem/TA share: {100*obs[obs.arm2=='malig'].stem.mean():.1f}%")

def test(mal_mask,nor_mask,col,label,minn=20):
    rows=[]
    for p,g in obs.groupby('patient'):
        mi=g[mal_mask(g)]; ni=g[nor_mask(g)]
        if len(mi)>=minn and len(ni)>=minn: rows.append((100*mi[col].mean(),100*ni[col].mean()))
    R=pd.DataFrame(rows,columns=['m','n'])
    if len(R)<3: print(f"  {label:58s} n={len(R)} (too few)"); return None
    d=R.m-R.n; up=int((d>0).sum())
    try: pw=wilcoxon(R.m,R.n).pvalue
    except Exception: pw=float('nan')
    ps=binomtest(up,len(R),0.5).pvalue
    print(f"  {label:58s} n={len(R):2d}  malig {R.m.median():5.1f}% vs normal {R.n.median():5.1f}% (median)  Δmean {d.mean():+5.1f}pp  higher {up}/{len(R)}  Wilcoxon p={pw:.1e}  sign p={ps:.1e}")
    return dict(test=label,n=len(R),malig_median=R.m.median(),normal_median=R.n.median(),delta_mean=d.mean(),higher=up,p_wilcoxon=pw,p_sign=ps)

res=[]
for col in ['core_RAW','core_SHAPE','troph_RAW']:
    print(f"\nB. within-patient tests — metric {col}")
    M=lambda g:g.arm2=='malig'; N=lambda g:g.arm2=='normal'
    res.append(test(M,N,col,'(i) all malignant vs all normal  [paper gate]'))
    res.append(test(M,lambda g:(g.arm2=='normal')&g.stem,col,'(ii) all malignant vs normal Stem/TA (crypt base)'))
    res.append(test(lambda g:(g.arm2=='malig')&g.stem,lambda g:(g.arm2=='normal')&g.stem,col,'(iii) malignant Stem/TA vs normal Stem/TA'))
    res.append(test(M,lambda g:(g.arm2=='normal')&~g.stem,col,'(ref) all malignant vs normal DIFFERENTIATED only'))
    if col=='core_RAW':
        for st_ in sorted(STEM)+['cE04 (Enterocyte 1)','cE05 (Enterocyte 2)','cE06 (Immature Goblet)']:
            res.append(test(lambda g,s_=st_:(g.arm2=='malig')&(g.subtype==s_),lambda g,s_=st_:(g.arm2=='normal')&(g.subtype==s_),col,f'(iv) {st_}: tumor vs normal same subtype',minn=15))
for r,c in zip([x for x in res],[None]*len(res)): pass
pd.DataFrame([r for r in res if r]).to_csv(OUT/'E36_tests.csv',index=False)

print("\n"+"="*86); print("C. WHOLE-TRANSCRIPTOME (2,000 HVGs, not modules): nearest normal subtype")
cent={}
for u in sums:
    p,arm,stp=u.split('|')
    if arm=='normal': cent.setdefault(stp,[0,0]); cent[stp][0]=cent[stp][0]+sums[u]; cent[stp][1]+=cnt[u]
cent={k:v[0]/v[1] for k,v in cent.items() if v[1]>=200}
names=list(cent); Cm=np.vstack([cent[k] for k in names]); gm=Cm.mean(0)
Cc=Cm-gm; Cc=Cc/np.linalg.norm(Cc,axis=1,keepdims=True)
def nearest(vec):
    v=vec-gm; v=v/np.linalg.norm(v); r=Cc@v; o=np.argsort(-r); return names[o[0]],r[o[0]],names[o[1]],r[o[1]],r
pt={}
for u in sums:
    p,arm,stp=u.split('|')
    if arm=='malig': pt.setdefault(p,[0,0]); pt[p][0]=pt[p][0]+sums[u]; pt[p][1]+=cnt[u]
rows=[]
for p,(s_,n_) in pt.items():
    if n_<50: continue
    a1,r1,a2,r2,_=nearest(s_/n_); rows.append(dict(patient=p,n=n_,nearest=a1,r=r1,second=a2,r2=r2))
P=pd.DataFrame(rows); P.to_csv(OUT/'E36_nearest_subtype_per_patient.csv',index=False)
print(f"  per-patient malignant pseudobulk (n={len(P)} patients) — nearest normal subtype:")
print(P.nearest.value_counts().to_string())
print(f"  share nearest to a Stem/TA subtype: {100*P.nearest.isin(STEM).mean():.0f}%")
ts={}
for u in sums:
    p,arm,stp=u.split('|')
    if arm=='malig': ts.setdefault(stp,[0,0]); ts[stp][0]=ts[stp][0]+sums[u]; ts[stp][1]+=cnt[u]
print("\n  per TUMOR subtype pseudobulk -> nearest normal subtype (r):")
for stp,(s_,n_) in sorted(ts.items(),key=lambda x:-x[1][1]):
    if n_<200: continue
    a1,r1,a2,r2,_=nearest(s_/n_); print(f"     Tumor {stp:40s} n={n_:6d} -> {a1} (r={r1:.2f}); 2nd {a2} (r={r2:.2f})")
print("DONE")
