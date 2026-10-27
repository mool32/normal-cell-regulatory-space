"""Round-4 next — INDEPENDENT VALIDATION of the named target, in a system with no cancer in it.
Referee 2.6 asked whether 'cytotrophoblast' is real or an artifact of a sparse, cross-platform gallery.
Reprogramming offers an orthogonal check: MEF -> iPSC reprogramming is known to throw off trophoblast-like
and extraembryonic intermediates (Schiebinger et al. 2019, 236,285 cells across 39 timepoints, day 0-18).
If the cancer-blind coordinate is really reading trophoblast-ness, it should flag a rising cytotrophoblast-
assigned subpopulation during reprogramming — in mouse, in another lab's data, with no cancer anywhere.
The modules are already scored in this object (obsm['module_scores'], 43 of the 44 via mouse orthologs)."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore')
from pathlib import Path
ROOT=Path('.').resolve()
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
sys.path.insert(0,'tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
core_idx=np.array([refcells.index(c) for c in CORE])
cyto=refcells.index('cytotrophoblasts'); mig=refcells.index('migrating cytotrophoblasts'); mega=refcells.index('megakaryocytes')
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
demag=lambda X:X-X.mean(1,keepdims=True); ident=lambda X:X
def refset(tf):
    R=tf(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; return mu,sd,rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
REF={'RAW':refset(ident),'SHAPE':refset(demag)}; TF={'RAW':ident,'SHAPE':demag}

b=sc.read_h5ad(ROOT/'paper6/results/schiebinger_scored.h5ad',backed='r')
mn=[str(x) for x in b.uns['module_names']]; M=np.asarray(b.obsm['module_scores'],dtype=float)
missing=[m for m in mods if m not in mn]
print(f"cells {M.shape[0]}, modules scored {len(mn)}/44; missing: {missing}",flush=True)
# reindex to the reference module order; a missing module gets the reference mean (neutral, no artificial zero)
idx={m:i for i,m in enumerate(mn)}
S=np.zeros((M.shape[0],len(mods)))
for j,m in enumerate(mods):
    S[:,j]=M[:,idx[m]] if m in idx else Rm[m].mean()
day=b.obs['time_numeric'].astype(float).values
res={}
for mm in ['RAW','SHAPE']:
    mu,sd,rpn=REF[mm]; pc=((TF[mm](S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; res[mm]=np.argmax(sim,axis=1)
rows=[]
for d in sorted(set(day)):
    k=day==d
    if k.sum()<100: continue
    r=dict(day=d,n=int(k.sum()))
    for mm in ['RAW','SHAPE']:
        t=res[mm][k]
        r[f'{mm}_cyto']=100*np.mean(np.isin(t,[cyto,mig]))
        r[f'{mm}_core']=100*np.mean(np.isin(t,core_idx))
        r[f'{mm}_mega']=100*np.mean(t==mega)
    rows.append(r)
D=pd.DataFrame(rows); D.to_csv(ROOT/'results/referee2_newdata/E28_reprogramming.csv',index=False)
print("\n"+"="*74); print("DOES THE CANCER-BLIND COORDINATE FLAG TROPHOBLAST-LIKE CELLS DURING REPROGRAMMING?")
print(D.round(2).to_string(index=False))
d0=D[D.day<=0.5]; peak=D.loc[D.SHAPE_cyto.idxmax()]
print(f"\n  start (day<=0.5): cytotrophoblast-assigned {d0.SHAPE_cyto.mean():.2f}% (shape-only), {d0.RAW_cyto.mean():.2f}% (raw)")
print(f"  peak: day {peak.day} -> {peak.SHAPE_cyto:.2f}% shape-only ({peak.RAW_cyto:.2f}% raw), n={int(peak.n)}")
print(f"  fold rise over the course: {peak.SHAPE_cyto/max(d0.SHAPE_cyto.mean(),1e-9):.1f}x (shape-only)")
from scipy.stats import spearmanr
print(f"  correlation of cytotrophoblast fraction with reprogramming day: rho = {spearmanr(D.day,D.SHAPE_cyto)[0]:+.2f} (p={spearmanr(D.day,D.SHAPE_cyto)[1]:.1e})")
print("\n[READING] a rising trophoblast-assigned population during reprogramming — a known off-target fate of")
print("          this system — would be independent evidence that the coordinate reads trophoblast-ness,")
print("          in mouse, in another lab's data, with no cancer anywhere in the analysis.")
print("DONE")
