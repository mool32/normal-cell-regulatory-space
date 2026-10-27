"""Round-4 (2) — SPECIFICITY TEST the paper lost when the lung negative control fell:
is the convergent state specific to MALIGNANCY, or does any invasive/remodelling epithelium enter it?
IPF (Habermann 2020, GSE135893): fibrotic lung, 22 donors, IPF vs Control, author cell types.
Fibrosis is the strongest non-malignant comparator available: chronically injured epithelium, ECM-remodelling,
and it contains 'Transitional AT2' (the aberrant-basaloid/KRT17+ regenerative intermediate repeatedly
described as 'cancer-like'). If ANY benign state should occupy the core, it is this one.
Outcomes: (a) IPF epithelium NOT deeper than control -> the state is malignancy-specific (restores the lost
negative control); (b) IPF deeper -> the state is 'invasive remodelling', not malignancy — an honest reframe.
Reported RAW and SHAPE-ONLY (the magnitude-independent metric that carries the paper)."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from scipy.stats import mannwhitneyu
from pathlib import Path
ROOT=Path('.').resolve()
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
demag=lambda X:X-X.mean(1,keepdims=True); ident=lambda X:X
def refset(tf):
    R=tf(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; return mu,sd,rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
REF={'RAW':refset(ident),'SHAPE':refset(demag)}; TF={'RAW':ident,'SHAPE':demag}

b=sc.read_h5ad(ROOT/"results/diseases/D3_IPF/ipf_habermann2020.h5ad")
m=(b.obs['population'].astype(str)=='Epithelial').values
print(f"IPF dataset {b.shape}; epithelial cells {int(m.sum())}",flush=True)
a=b[m].copy(); del b; gc.collect()
a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
status=a.obs['Status'].astype(str).values; donor=a.obs['Sample_Name'].astype(str).values; ct=a.obs['celltype'].astype(str).values
print(f"scored {len(S)} epithelial cells ({int((status=='Disease').sum())} IPF / {int((status=='Control').sum())} control)",flush=True)

res={}
for mm in ['RAW','SHAPE']:
    mu,sd,rpn=REF[mm]; pc=((TF[mm](S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
    res[mm]=dict(core=np.isin(t,core_idx),cyto=np.isin(t,cyto_ids),top1=t)

print("\n"+"="*72); print("IS THE CONVERGENT STATE MALIGNANCY-SPECIFIC? (fibrotic vs control lung epithelium)")
rows=[]
for mm in ['RAW','SHAPE']:
    isc=res[mm]['core']
    per=[]
    for d in pd.unique(donor):
        for st in ['Control','Disease']:
            k=(donor==d)&(status==st)
            if k.sum()>=20: per.append(dict(donor=d,status=st,n=int(k.sum()),core=100*isc[k].mean()))
    P=pd.DataFrame(per); rows.append(P.assign(metric=mm))
    c=P[P.status=='Control'].core; i=P[P.status=='Disease'].core
    try: p=mannwhitneyu(i,c,alternative='greater').pvalue
    except Exception: p=float('nan')
    print(f"  [{mm:5s}] control lung {c.mean():5.1f}% (n={len(c)} donors)  vs  IPF {i.mean():5.1f}% (n={len(i)})   Δ {i.mean()-c.mean():+5.1f}pp   one-sided p={p:.3f}")
pd.concat(rows).to_csv(ROOT/'results/referee2_newdata/E24_ipf_specificity.csv',index=False)

print("\n  per epithelial cell type (SHAPE-ONLY core %, IPF cells only; is the 'cancer-like' aberrant pool in the core?):")
isc=res['SHAPE']['core']; dis=status=='Disease'
tab=[]
for c in pd.unique(ct):
    k=(ct==c)&dis
    if k.sum()>=30: tab.append((c,int(k.sum()),100*isc[k].mean()))
for c,n,v in sorted(tab,key=lambda x:-x[2])[:12]: print(f"     {c:26s} n={n:5d}  core {v:5.1f}%")
tr=[x for x in tab if 'Transitional' in x[0] or 'Basal' in x[0] or 'KRT5' in x[0]]
print("\n  the regenerative/aberrant pool specifically:")
for c,n,v in tr: print(f"     {c:26s} n={n:5d}  core {v:5.1f}%")
print("\n  (compare: malignant epithelium in the paper's cohorts sits at ~13-30% SHAPE-ONLY core;")
print("   matched-normal at ~3-18%. Where does fibrotic epithelium fall?)")
print("DONE")
