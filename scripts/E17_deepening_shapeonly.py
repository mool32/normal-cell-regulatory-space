"""Referee round-2 3.5 (MAKE-OR-BREAK) — does the within-patient deepening survive the magnitude/shape-only
decomposition? The n_genes control (2.3) removes a sequencing-depth confound, NOT the brightness confound:
malignant epithelium is transcriptionally hyperactive, so it may be assigned to the high-output core simply
because it is brighter. Here we recompute the per-patient core fraction with the overall MAGNITUDE removed
(per-cell centering across the 44 modules) and ask whether patient-matched malignant epithelium still sits
deeper in the SHAPE-ONLY core than matched-normal. Flagship cohort: Pelka colorectal.
If shape-only deepening persists -> the deepening is a genuine module-pattern shift, not just high output."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import wilcoxon
ROOT=Path('.').resolve()
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])

def build_ref(transform):
    R=transform(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; rpn=rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
    return mu,sd,rpn
def assign_core(S,mu,sd,rpn):
    pc=((S-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
    return np.isin(t,core_idx)
ident=lambda X:X
demag=lambda X:X-X.mean(1,keepdims=True)   # remove per-cell overall magnitude

# ---- load Pelka epithelial, split malig/normal ----
P=ROOT/'results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad'
b=sc.read_h5ad(P,backed='r'); cm=b.obs['ClusterMidway'].astype(str).values
a=b[np.isin(cm,['EpiT','Epi'])].to_memory(); del b
if 'feature_name' in a.var.columns and str(a.var_names[0]).startswith('ENSG'): a.var_names=a.var['feature_name'].astype(str).values
a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=1000].copy()
arm=np.where(a.obs['ClusterMidway'].astype(str).values=='EpiT','malig','normal')
pcol='donor_id' if 'donor_id' in a.obs else 'PatientTypeID'; pat=a.obs[pcol].astype(str).values
if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
print(f"Pelka epithelial {a.n_obs}: malig {int((arm=='malig').sum())} normal {int((arm=='normal').sum())}",flush=True)

rows=[]
for tag,tf in [('RAW',ident),('SHAPE-ONLY',demag)]:
    mu,sd,rpn=build_ref(tf); isc=assign_core(tf(S),mu,sd,rpn)
    for p in pd.unique(pat):
        mi=(pat==p)&(arm=='malig'); ni=(pat==p)&(arm=='normal')
        if mi.sum()<20 or ni.sum()<20: continue
        rows.append(dict(metric=tag,patient=p,malig=100*isc[mi].mean(),normal=100*isc[ni].mean()))
df=pd.DataFrame(rows); df.to_csv(ROOT/'results/referee2_newdata/E17_shapeonly_deepening.csv',index=False)
print("\n"+"="*66); print("SHAPE-ONLY DEEPENING (referee 3.5, make-or-break) — Pelka colorectal")
for tag in ['RAW','SHAPE-ONLY']:
    d=df[df.metric==tag]; dif=d.malig.values-d.normal.values
    p=wilcoxon(d.malig,d.normal).pvalue
    print(f"  [{tag:10s}] malig {d.malig.mean():5.1f}% vs normal {d.normal.mean():5.1f}%  ->  Δ {dif.mean():+5.1f}pp, deeper {int((dif>0).sum())}/{len(d)}, paired p={p:.2e}")
raw=df[df.metric=='RAW']; sh=df[df.metric=='SHAPE-ONLY']
rd=(raw.malig.values-raw.normal.values).mean(); sd_=(sh.malig.values-sh.normal.values).mean()
print(f"\n  shape-only retains {100*sd_/rd:.0f}% of the raw deepening")
print(f"  [VERDICT] deepening {'SURVIVES magnitude removal — a genuine module-shape shift, not just high output' if sd_>0 and (sh.malig.values>sh.normal.values).mean()>0.7 else 'is MOSTLY MAGNITUDE — collapses shape-only (bad for the paper)'}")
print("DONE")
