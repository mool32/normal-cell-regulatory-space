"""Referee 2.1 (THE central issue) — is the 6.3x contraction convergence or an averaging/magnitude artifact?
The observed statistic (0.158) is the mean pairwise cosine distance among 25 malignant CENTROIDS (means over
thousands of cells). The paper compared it to single cell-type centroids (means-of-one) — not matched for averaging.
Here we build MATCHED-AVERAGING nulls (25 random size-matched cell SETS -> centroids -> pairwise distance) and
repeat everything MATCHED-MAGNITUDE (magnitude removed), to see whether 0.158 beats the shrinkage floor.
  Null-M: repartition the pooled malignant cells into 25 random size-matched groups (pure averaging floor of THESE cells).
  Null-N: draw 25 size-matched sets from a diverse pool of NORMAL single cells (does heterogeneous averaging alone
          reproduce the tightness?).
Observed origins spread (0.99) is the reference for 'no convergence'."""
import numpy as np, pandas as pd, json, sys, warnings, anndata as ad
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12

d=np.load(DD/'../E7_S_raw.npz',allow_pickle=True) if (DD/'../E7_S_raw.npz').exists() else np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True)
S=d['S'].astype(float); cid=d['cid']; uc=sorted(set(cid.tolist()))
sizes=[int((cid==c).sum()) for c in uc]; print(f"malignant cells {S.shape}, {len(uc)} cancers, sizes {min(sizes)}-{max(sizes)}",flush=True)
pl=np.load(DD/'normalpool_local.npz',allow_pickle=True); Sn=pl['S'].astype(float); ptis=pl['tissue']
print(f"normal pool {Sn.shape} across tissues {dict(pd.Series(ptis).value_counts())}",flush=True)

def proj(X):   # module scores -> normalized 9-D
    pc=((X-mu)/sd)@load; return pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
def demag(X):  # remove per-cell magnitude (mean across modules)
    return X-X.mean(1,keepdims=True)
def centroids_dist(P, groups):
    C=[]
    for g in groups:
        c=P[g].mean(0); C.append(c/(np.linalg.norm(c)+1e-12))
    C=np.array(C); s=C@C.T; iu=np.triu_indices(len(C),1); return float((1-s[iu]).mean())

rng=np.random.default_rng(0)
for tag,fn in [('RAW',lambda X:X),('SHAPE-ONLY (magnitude removed)',demag)]:
    Pm=proj(fn(S)); Pn=proj(fn(Sn))
    real_groups=[np.where(cid==c)[0] for c in uc]
    obs=centroids_dist(Pm,real_groups)
    # Null-M: random size-matched partitions of the malignant pool
    nm=[]
    for _ in range(1000):
        perm=rng.permutation(len(S)); gs=[]; i=0
        for n in sizes: gs.append(perm[i:i+n]); i+=n
        nm.append(centroids_dist(Pm,gs))
    nm=np.array(nm)
    # Null-N: 25 COHERENT normal micro-populations (seed + N_c-1 nearest neighbours) = distinct
    #          normal populations formed by the SAME averaging (means-of-many), representing "no convergence".
    G=Pn@Pn.T  # cosine sims in the normal pool (unit vectors)
    nn=[]
    for _ in range(300):
        seeds=rng.choice(len(Sn),25,replace=False); gs=[]
        for s,n in zip(seeds,sizes):
            nbr=np.argpartition(-G[s],min(n,len(Sn)-1))[:min(n,len(Sn))]; gs.append(nbr)
        nn.append(centroids_dist(Pn,gs))
    nn=np.array(nn)
    print(f"\n===== {tag} =====")
    print(f"  OBSERVED 25 malignant-centroid mean pairwise cosine distance = {obs:.3f}")
    print(f"  Null-M (shuffle malignant labels, pure averaging floor):        mean {nm.mean():.3f} [{np.percentile(nm,2.5):.3f},{np.percentile(nm,97.5):.3f}]")
    print(f"  Null-N (25 DISTINCT normal micro-populations, matched averaging): mean {nn.mean():.3f} [{np.percentile(nn,2.5):.3f},{np.percentile(nn,97.5):.3f}]  -> observed percentile {100*(nn<obs).mean():.1f}%")
    print(f"  => observed {'TIGHTER than' if obs<nn.mean() else 'NOT tighter than'} 25 distinct normal populations of matched size (the meaningful convergence test)")
print("\n(reference: origin spread ~0.99 among the 25 distinct assigned normal cell types)")
print("DONE")
