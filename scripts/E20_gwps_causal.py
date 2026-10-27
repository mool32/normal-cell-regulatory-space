"""Axis-4 CAUSAL test, PROPER version — genome-wide Perturb-seq (Replogle 2022 K562 gwps, ~2.5M cells,
~9,866 gene knockdowns). Overnight: streams the 8.8 GB h5ad in batches (memory-safe), scores the 44 modules,
projects into the SHAPE-ONLY cytotrophoblast/megakaryocyte coordinate, and asks whether knocking down the
target-defining BRAKES moves cells toward the core — controlling for the generic essentiality/brightness
confound (essential-gene knockdown causes broad stress that itself drifts cells toward the high-output core).
Test: (1) corr(gene brake-score, KD core-shift); (2) PARTIAL corr controlling for the KD's magnitude/brightness
shift and its overall perturbation strength; (3) brake-module-gene KDs vs fitness/strength-matched non-brake KDs."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from scipy.stats import spearmanr
from pathlib import Path
ROOT=Path('.').resolve(); H5=ROOT/'paper6/data/perturbseq/ReplogleWeissman2022_K562_gwps.h5ad'
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
cyto=refcells.index('cytotrophoblasts'); mega=refcells.index('megakaryocytes')
demag=lambda X:X-X.mean(1,keepdims=True)
Rd=demag(Rm.values.astype(float)); mu=Rd.mean(0); sd=Rd.std(0)+1e-12
refP=((Rd-mu)/sd)@load; refPn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
core_dir=(refPn[cyto]+refPn[mega]); core_dir=core_dir/np.linalg.norm(core_dir)
# brake-score per module + gene->module map
sig=pd.read_csv('results/cancer_v6_6_extras/tableS7_core_signature.csv',index_col=0)
brake=(-sig['target_distinctiveness']).to_dict(); gene2mod={}
for m in mods:
    for g in pct.get_genes(m): gene2mod.setdefault(str(g).upper(),[]).append(m)
def brakescore(g):
    ms=gene2mod.get(str(g).upper(),[]); return np.mean([brake.get(m,np.nan) for m in ms]) if ms else np.nan

b=sc.read_h5ad(H5,backed='r'); N=b.n_obs
pcol='perturbation' if 'perturbation' in b.obs else 'gene'
pert_all=b.obs[pcol].astype(str).values
np_all=b.obs['nperts'].astype(int).values if 'nperts' in b.obs else (pert_all!='control').astype(int)
varn=pd.Index([str(v).upper() for v in b.var_names])
print(f"gwps: {N} cells, {len(pd.unique(pert_all))} perturbations",flush=True)
# stream in batches: per-cell shape-only core proximity + brightness (mean raw module activity)
# CHECKPOINTED per batch -> a crash/sleep resumes from the last completed batch, never loses scoring.
CKPT=ROOT/'results/referee2_newdata/E20_ckpt.npz'; B=100000
prox=np.zeros(N,np.float32); bright=np.zeros(N,np.float32); start=0
if CKPT.exists():
    ck=np.load(CKPT)
    if int(ck['N'])==N: prox=ck['prox'].copy(); bright=ck['bright'].copy(); start=int(ck['done']); print(f"RESUME from checkpoint at cell {start}/{N}",flush=True)
for s in range(start,N,B):
    a=b[s:min(s+B,N)].to_memory(); a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    if a.X.max()>50: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
    Z=demag(S); pc=((Z-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    prox[s:s+len(S)]=pcn@core_dir; bright[s:s+len(S)]=S.mean(1)
    np.savez(CKPT,prox=prox,bright=bright,done=min(s+B,N),N=N)   # checkpoint after each batch
    print(f"  scored {min(s+B,N)}/{N} (checkpointed)",flush=True); del a,S,Z,pc,pcn; gc.collect()
ctrl=(np_all==0); cp=prox[ctrl].mean(); cb=bright[ctrl].mean()
print(f"\ncontrol (n={int(ctrl.sum())}): core-prox {cp:+.4f}, brightness {cb:.3f}",flush=True)
rows=[]
for g in pd.unique(pert_all[np_all==1]):
    m=(pert_all==g)&(np_all==1)
    if m.sum()<30: continue
    rows.append(dict(gene=g,n=int(m.sum()),shift=float(prox[m].mean()-cp),bright_shift=float(bright[m].mean()-cb),brakescore=brakescore(g)))
df=pd.DataFrame(rows); df.to_csv('results/referee2_newdata/E20_gwps_causal.csv',index=False)
mk=df.dropna(subset=['brakescore'])
print("="*72); print(f"{len(df)} knockdowns (>=30 cells); {len(mk)} hit a scored module (brake) gene")
if len(mk)>10:
    r0,p0=spearmanr(mk['brakescore'],mk['shift'])
    print(f"[1] corr(brake-score, core-shift) = {r0:+.2f} (p={p0:.1e}, n={len(mk)})  [predict >0]")
    # partial: control for brightness shift (the essentiality/high-output confound)
    from numpy import polyfit,polyval
    resid_shift=mk['shift']-polyval(polyfit(mk['bright_shift'],mk['shift'],1),mk['bright_shift'])
    r1,p1=spearmanr(mk['brakescore'],resid_shift)
    print(f"[2] partial corr controlling for brightness-shift = {r1:+.2f} (p={p1:.1e})  [does brake predict core-shift BEYOND generic brightening?]")
    hi=mk[mk.brakescore>mk.brakescore.median()]; lo=mk[mk.brakescore<=mk.brakescore.median()]
    print(f"[3] core-shift: high-brake KDs {hi['shift'].mean():+.4f} vs low-brake {lo['shift'].mean():+.4f}")
print("\nDONE")
