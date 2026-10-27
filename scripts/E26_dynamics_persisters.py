"""Round-4 next — DYNAMICS: the paper's biggest self-declared limitation is that a convergent endpoint in
static snapshots is not a proven dynamical attractor. Aissa et al. 2021 (GSE134839) gives the missing time
axis in the relevant system: human PC9 lung adenocarcinoma cells profiled through an erlotinib time course
(D0, D1, D2, D4, D9, D11) as they pass through drug-tolerant persistence and begin to expand again.
Question: does the population MOVE in the cancer-blind coordinate over time, and in which direction — do
surviving cells enter and hold the convergent state (attractor-like), or leave it?
Reports raw and shape-only core fraction and the continuous core proximity per timepoint."""
import scanpy as sc, anndata as ad, json, numpy as np, pandas as pd, sys, glob, os, re, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve()
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
core_idx=np.array([refcells.index(c) for c in CORE])
cyto=refcells.index('cytotrophoblasts'); mega=refcells.index('megakaryocytes')
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
demag=lambda X:X-X.mean(1,keepdims=True); ident=lambda X:X
def refset(tf):
    R=tf(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; return mu,sd,rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
REF={'RAW':refset(ident),'SHAPE':refset(demag)}; TF={'RAW':ident,'SHAPE':demag}
cdir={}
for m in ['RAW','SHAPE']:
    rpn=REF[m][2]; v=rpn[cyto]+rpn[mega]; cdir[m]=v/np.linalg.norm(v)

rows=[]
for f in sorted(glob.glob(str(ROOT/'paper6/data/aissa2021/extracted/*.dge.txt.gz'))):
    day=int(re.search(r'_D(\d+)\.',os.path.basename(f)).group(1))
    df=pd.read_csv(f,sep='\t',index_col=0)
    df.index=[str(g).upper() for g in df.index]
    if pd.Index(df.index).has_duplicates: df=df.groupby(level=0).sum()
    a=ad.AnnData(df.T.values.astype(np.float32)); a.var_names=list(df.index); a.var_names_make_unique(); del df
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
    r=dict(day=day,n=len(S))
    for m in ['RAW','SHAPE']:
        mu,sd,rpn=REF[m]; pc=((TF[m](S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
        r[f'{m}_core']=100*np.isin(t,core_idx).mean()
        r[f'{m}_prox']=float((pcn@cdir[m]).mean())
        r[f'{m}_cyto']=100*(t==cyto).mean()
    rows.append(r); print(f"  D{day}: n={r['n']}  RAW core {r['RAW_core']:.1f}% prox {r['RAW_prox']:+.3f} | SHAPE core {r['SHAPE_core']:.1f}% prox {r['SHAPE_prox']:+.3f}",flush=True)
    del a,S; gc.collect()
d=pd.DataFrame(rows).sort_values('day'); d.to_csv(ROOT/'results/referee2_newdata/E26_dynamics.csv',index=False)
print("\n"+"="*70); print("DYNAMICS UNDER TARGETED THERAPY (PC9 + erlotinib, human)")
print(d.round(3).to_string(index=False))
b=d.iloc[0]; e9=d[d.day>=9]
print(f"\n  baseline D0: RAW core {b.RAW_core:.1f}% / SHAPE {b.SHAPE_core:.1f}%")
print(f"  persister/regrowth (D9-11): RAW {e9.RAW_core.mean():.1f}% / SHAPE {e9.SHAPE_core.mean():.1f}%")
print(f"  Δ over the course: RAW {e9.RAW_core.mean()-b.RAW_core:+.1f}pp | SHAPE {e9.SHAPE_core.mean()-b.SHAPE_core:+.1f}pp")
mono=all(d.SHAPE_prox.values[i]<=d.SHAPE_prox.values[i+1]+1e-9 for i in range(len(d)-1))
print(f"  shape-only proximity monotone increasing across the course: {mono}")
print("\n[READING] a population that moves INTO the state and holds it is attractor-like in this system;")
print("          a population that leaves it argues the state is tied to proliferation the drug removes.")
print("DONE")
