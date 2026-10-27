"""Round-4 BONUS — does the convergent state persist in METASTASIS, or is it a primary-tumour property?
Noticed in E22: breast IDC 21.8% -> lymph-node metastasis 13.4%. Test it properly, within patient, twice:
  BREAST  GSE195861: P1..P6 (primary IDC) vs M1..M6 (matched LN metastasis of the SAME patient)
  LIVER   HCC (GSE149614, Lu et al. 2022; the local path results/cancer_sun2021_hcc/ is a legacy misnomer): site=='Tumor' vs site=='Lymph' (LN met) and 'PVTT' (portal-vein tumour thrombus),
          same patients — an independent organ with intravascular and nodal spread.
If metastases sit LOWER, the state is tied to invasion/primary growth rather than colonisation — a new,
testable claim. RAW and SHAPE-ONLY (magnitude-removed) both reported."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, glob, os, re, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from scipy.stats import wilcoxon
from pathlib import Path
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
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
def corefracs(S):
    out={}
    for m in ['RAW','SHAPE']:
        mu,sd,rpn=REF[m]; pc=((TF[m](S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; out[m]=100*np.isin(np.argmax(sim,axis=1),core_idx).mean()
    return out
def score_ad(a,epcam_gate=True):
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    if epcam_gate:
        def col(g):
            if g not in set(a.var_names): return np.zeros(a.n_obs)
            v=a[:,g].X; return np.asarray(v.todense()).ravel() if hasattr(v,'todense') else np.asarray(v).ravel()
        a=a[(col('EPCAM')>0)&(col('PTPRC')==0)].copy()
    if a.n_obs<20: return None,0
    if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
    n=a.n_obs; del a; gc.collect(); return S,n

rows=[]
# ---------- BREAST: P# primary vs M# matched LN metastasis ----------
for d in sorted(glob.glob(str(DD/'dcis195861/GSM*'))):
    nm=os.path.basename(d); m=re.search(r'_(P|M)(\d)_',nm)
    if not m: continue
    try: a=sc.read_10x_mtx(str(d),var_names='gene_symbols',make_unique=True)
    except Exception: continue
    S,n=score_ad(a)
    if S is None: continue
    cf=corefracs(S)
    rows.append(dict(cohort='breast',patient=m.group(2),site=('primary' if m.group(1)=='P' else 'metastasis'),
                     n=n,RAW=cf['RAW'],SHAPE=cf['SHAPE']))
    print(f"  breast {m.group(1)}{m.group(2)}: n={n} RAW {cf['RAW']:.1f}% SHAPE {cf['SHAPE']:.1f}%",flush=True)

# ---------- LIVER HCC: same patients, Tumor vs Lymph / PVTT ----------
b=sc.read_h5ad(ROOT/'results/cancer_sun2021_hcc/sun2021.h5ad')
hep=(b.obs['celltype'].astype(str)=='Hepatocyte').values
site=b.obs['site'].astype(str).values; pat=b.obs['patient'].astype(str).values
for p in pd.unique(pat):
    for s,lab in [('Tumor','primary'),('Lymph','metastasis'),('PVTT','vascular-thrombus')]:
        k=hep&(pat==p)&(site==s)
        if k.sum()<20: continue
        S,n=score_ad(b[k].copy(),epcam_gate=False)
        if S is None: continue
        cf=corefracs(S)
        rows.append(dict(cohort='HCC',patient=p,site=lab,n=n,RAW=cf['RAW'],SHAPE=cf['SHAPE']))
        print(f"  HCC {p} {lab}: n={n} RAW {cf['RAW']:.1f}% SHAPE {cf['SHAPE']:.1f}%",flush=True)
del b; gc.collect()

df=pd.DataFrame(rows); df.to_csv(DD/'E25_primary_vs_met.csv',index=False)
print("\n"+"="*70); print("DOES THE CONVERGENT STATE PERSIST IN METASTASIS? (within patient)")
for coh in df.cohort.unique():
    g=df[df.cohort==coh]
    piv=g.pivot_table(index='patient',columns='site',values=['RAW','SHAPE'])
    for m in ['RAW','SHAPE']:
        if ('primary' in piv[m].columns) and ('metastasis' in piv[m].columns):
            ok=piv[m][['primary','metastasis']].dropna()
            if len(ok)>=3:
                dif=ok['metastasis']-ok['primary']
                try: p=wilcoxon(ok['metastasis'],ok['primary']).pvalue
                except Exception: p=float('nan')
                print(f"  [{coh:6s} {m:5s}] {len(ok)} patients: primary {ok['primary'].mean():5.1f}% -> metastasis {ok['metastasis'].mean():5.1f}%  "
                      f"(Δ {dif.mean():+5.1f}pp); LOWER in metastasis {int((dif<0).sum())}/{len(ok)}; paired p={p:.3f}")
    if 'vascular-thrombus' in set(g.site):
        v=g[g.site=='vascular-thrombus']; pr=g[g.site=='primary']
        print(f"  [{coh} PVTT] vascular thrombus RAW {v.RAW.mean():.1f}% vs primary {pr.RAW.mean():.1f}% (n={len(v)} samples)")
print("\nDONE")
