"""Referee 2.3 — does the within-patient deepening survive WITHIN-PAIR complexity matching?
Malignant epithelium carries more detected genes than matched-normal; the panel-wide |r|=0.03 is
across-cell scatter, not the paired offset. Here, PER PATIENT, we resample malignant cells to match the
matched-normal n_genes distribution (nearest-n_genes, with replacement), then recompute the core fraction.
If the deepening persists after matching, it is not a complexity artifact. Flagship cohort: Pelka colorectal."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import wilcoxon
ROOT=Path('.').resolve()
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
refPC=(Rm.values-mu)/sd@load; refPn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)

P=ROOT/'results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad'
b=sc.read_h5ad(P,backed='r'); cm=b.obs['ClusterMidway'].astype(str).values
a=b[np.isin(cm,['EpiT','Epi'])].to_memory(); del b
if 'feature_name' in a.var.columns and str(a.var_names[0]).startswith('ENSG'):
    a.var_names=a.var['feature_name'].astype(str).values
a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
ng=np.asarray((a.X>0).sum(1)).ravel(); keep=ng>=1000; a=a[keep].copy(); ng=ng[keep]
arm=np.where(a.obs['ClusterMidway'].astype(str).values=='EpiT','malig','normal')
pcol='donor_id' if 'donor_id' in a.obs else 'PatientTypeID'
pat=a.obs[pcol].astype(str).values
if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values
Z=(S-mu)/sd; pc=Z@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; top1=np.argmax(sim,axis=1)
isc=np.isin(top1,core_idx)
print(f"Pelka epithelial cells {a.n_obs}: malig {int((arm=='malig').sum())} normal {int((arm=='normal').sum())}",flush=True)

rng=np.random.default_rng(0); rows=[]
for p in pd.unique(pat):
    mi=np.where((pat==p)&(arm=='malig'))[0]; ni=np.where((pat==p)&(arm=='normal'))[0]
    if len(mi)<20 or len(ni)<20: continue
    mg_ng=ng[mi]; nr_ng=ng[ni]; order=np.argsort(mg_ng); mg_sorted=mg_ng[order]; mi_sorted=mi[order]
    # nearest-n_genes match: for each normal cell, pick malignant cell with closest n_genes (avg over 20 draws)
    mfr=[]; mcy=[]
    for _ in range(20):
        picks=[]
        for v in nr_ng:
            j=np.searchsorted(mg_sorted,v)
            cands=[k for k in (j-1,j) if 0<=k<len(mg_sorted)]
            best=min(cands,key=lambda k:abs(mg_sorted[k]-v)); picks.append(mi_sorted[best])
        picks=np.array(picks); mfr.append(100*isc[picks].mean()); mcy.append(100*np.isin(top1[picks],cyto_ids).mean())
    rows.append(dict(patient=p,n_malig=len(mi),n_normal=len(ni),
        ng_malig=float(np.median(mg_ng)),ng_normal=float(np.median(nr_ng)),
        normal_core=100*isc[ni].mean(), malig_core_raw=100*isc[mi].mean(), malig_core_matched=float(np.mean(mfr)),
        normal_cyto=100*np.isin(top1[ni],cyto_ids).mean(), malig_cyto_matched=float(np.mean(mcy))))
df=pd.DataFrame(rows); df.to_csv(ROOT/'results/referee2_newdata/E15_pelka_complexity_matched.csv',index=False)
print(f"\npatients: {len(df)} | median n_genes malig {df.ng_malig.median():.0f} vs normal {df.ng_normal.median():.0f} (the confound)")
print("\n[deepening RAW vs COMPLEXITY-MATCHED] (malig core - normal core):")
df['raw_diff']=df.malig_core_raw-df.normal_core; df['matched_diff']=df.malig_core_matched-df.normal_core
print(f"  RAW      malig {df.malig_core_raw.mean():.1f}% vs normal {df.normal_core.mean():.1f}%  -> Δ {df.raw_diff.mean():+.1f}pp, deeper {int((df.raw_diff>0).sum())}/{len(df)}, p={wilcoxon(df.malig_core_raw,df.normal_core).pvalue:.2e}")
print(f"  MATCHED  malig {df.malig_core_matched.mean():.1f}% vs normal {df.normal_core.mean():.1f}%  -> Δ {df.matched_diff.mean():+.1f}pp, deeper {int((df.matched_diff>0).sum())}/{len(df)}, p={wilcoxon(df.malig_core_matched,df.normal_core).pvalue:.2e}")
retention=100*df.matched_diff.mean()/df.raw_diff.mean()
print(f"\n[VERDICT] complexity-matched deepening retains {retention:.0f}% of the raw effect; "
      f"{'SURVIVES — not a complexity artifact' if df.matched_diff.mean()>0 and (df.matched_diff>0).sum()>len(df)*0.7 else 'attenuated — see numbers'}")
print("DONE")
