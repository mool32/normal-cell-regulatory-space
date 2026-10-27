"""Round-4 next — does the convergent state actually behave as an IMMUNE-EVADING state?
The paper's biological rationale for the cytotrophoblast is that it physiologically invades, proliferates
AND evades immunity. The immune-evasion leg has never been tested. Pelka et al. colorectal carries ALL cell
compartments per patient, so we can ask, across tumours: do patients whose malignant epithelium sits DEEPER
in the state have LESS immune infiltration?
Critical confounder: mismatch-repair status — MMRd tumours are hyper-infiltrated — so we stratify by
MMRStatusTumor and report within-stratum as well as pooled.
Infiltration is measured as the fraction of a sample's cells that are T/NK/ILC (and myeloid), acknowledging
dissociation bias. Deepening is the malignant-epithelium core fraction, raw and shape-only."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from scipy.stats import spearmanr
from pathlib import Path
ROOT=Path('.').resolve()
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

P=ROOT/'results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad'
b=sc.read_h5ad(P,backed='r')
obs=b.obs
pcol='donor_id' if 'donor_id' in obs else 'PatientTypeID'
# ---- composition per tumour sample from obs alone (no expression needed) ----
tum=(obs['disease'].astype(str)=='colon adenocarcinoma').values
comp=[]
for p in pd.unique(obs[pcol].astype(str)[tum]):
    k=tum&(obs[pcol].astype(str).values==p)
    top=obs['ClusterTop'].astype(str).values[k]
    n=len(top)
    if n<200: continue
    comp.append(dict(patient=p,n_cells=n,
                     tnk=100*np.mean(top=='TNKILC'),myeloid=100*np.mean(top=='Myeloid'),
                     bcell=100*np.mean(np.isin(top,['B','Plasma'])),
                     mmr=str(obs['MMRStatusTumor'].astype(str).values[k][0])))
C=pd.DataFrame(comp); print(f"tumour samples with composition: {len(C)}",flush=True)

# ---- malignant epithelium core fraction per patient ----
mal=tum&(obs['ClusterMidway'].astype(str)=='EpiT').values
a=b[mal].to_memory(); del b; gc.collect()
if 'feature_name' in a.var.columns and str(a.var_names[0]).startswith('ENSG'): a.var_names=a.var['feature_name'].astype(str).values
a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=1000].copy()
if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
pat=a.obs[pcol].astype(str).values
core={}
for m in ['RAW','SHAPE']:
    mu,sd,rpn=REF[m]; pc=((TF[m](S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; core[m]=np.isin(np.argmax(sim,axis=1),core_idx)
rows=[]
for p in pd.unique(pat):
    k=pat==p
    if k.sum()<20: continue
    rows.append(dict(patient=p,n_malig=int(k.sum()),RAW=100*core['RAW'][k].mean(),SHAPE=100*core['SHAPE'][k].mean()))
D=pd.DataFrame(rows).merge(C,on='patient',how='inner')
D.to_csv(ROOT/'results/referee2_newdata/E27_immune.csv',index=False)
print(f"patients with both malignant core and composition: {len(D)}")
print("\n"+"="*72); print("IS THE CONVERGENT STATE IMMUNE-EVADING? (deeper core <-> less infiltration?)")
for m in ['RAW','SHAPE']:
    for imm in ['tnk','myeloid']:
        r,p_=spearmanr(D[m],D[imm])
        print(f"  [{m:5s}] core vs {imm:7s} infiltration: Spearman ρ = {r:+.2f} (p={p_:.3f}, n={len(D)})   [evasion predicts ρ<0]")
print("\n  stratified by mismatch-repair status (MMRd tumours are hyper-infiltrated — the key confounder):")
for mm,g in D.groupby('mmr'):
    if len(g)<6: print(f"     {mm}: n={len(g)} (too few)"); continue
    r1,p1=spearmanr(g['SHAPE'],g['tnk'])
    print(f"     {mm:6s} n={len(g):2d}: median core(SHAPE) {g.SHAPE.median():5.1f}%, median T/NK {g.tnk.median():5.1f}%  |  within-stratum ρ(core,T/NK) = {r1:+.2f} (p={p1:.2f})")
print("\n  (dissociation bias makes cell-fraction infiltration noisy; treat as a first-pass test)")
print("DONE")
