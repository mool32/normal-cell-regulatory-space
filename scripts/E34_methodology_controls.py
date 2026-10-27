"""QC-review methodology items M2, M3, M4, M7, M9 — computed once, for the SI."""
import json, numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve()
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
CORE3=[refcells.index(c) for c in ['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']]
EX=[refcells.index(c) for c in ['erythrocytes','platelets']]
refP=((Rm.values-mu)/sd)@load

# ---------- M7: nearest-neighbour family recovery, with a family-size-preserving null ----------
fam=json.load(open(ROOT/'results/fam.json'))
keep=[i for i,c in enumerate(refcells) if c not in {'erythrocytes','platelets','reticulocytes'} and c in fam]
P=refP[keep]; cells=[refcells[i] for i in keep]; f=np.array([fam[c] for c in cells])
U=P/np.linalg.norm(P,axis=1,keepdims=True); Sim=U@U.T; np.fill_diagonal(Sim,-9)
nn=np.argmax(Sim,axis=1); obs=float(np.mean(f[nn]==f))
sizes=pd.Series(f).value_counts().values; N=len(f)
analytic=float(sum(n*(n-1) for n in sizes)/(N*(N-1)))          # family-size-preserving pair rate
rng=np.random.default_rng(0); null=[]
for _ in range(10000):
    fp=rng.permutation(f); null.append(np.mean(fp[nn]==fp))
null=np.array(null); pval=(np.sum(null>=obs)+1)/(len(null)+1)
lo,hi=np.percentile(rng.binomial(N,obs,20000)/N,[2.5,97.5])
print(f"[M7] NN same-family: observed {100*obs:.0f}% (95% CI {100*lo:.0f}–{100*hi:.0f}%), N={N} cell types")
print(f"     naive 1/K chance (as previously quoted): {100/len(set(f)):.0f}%   |   family-size-preserving chance: {100*analytic:.0f}%")
print(f"     label-permutation null mean {100*null.mean():.0f}% (max {100*null.max():.0f}%), permutation P = {pval:.1e} (10,000 permutations)")

# ---------- M2: modules that share genes; correlations on disjoint gene sets ----------
g2m={}
for m in mods:
    for g in pct.get_genes(m): g2m.setdefault(str(g).upper(),[]).append(m)
shared={g:ms for g,ms in g2m.items() if len(ms)>1}
print(f"\n[M2] genes shared between modules: {len(shared)} genes covering {len(set(sum(shared.values(),[])))} modules")
for g,ms in sorted(shared.items()): print(f"     {g:9s} -> {', '.join(ms)}")
Zr=(Rm.values-mu)/sd; C=np.corrcoef(Zr.T)
iu=np.triu_indices(len(mods),1); pairs=[(mods[i],mods[j],C[i,j]) for i,j in zip(*iu)]
hi_pairs=[(a,b,r) for a,b,r in pairs if abs(r)>0.7]
def shares(a,b): return any(a in ms and b in ms for ms in g2m.values())
print(f"     module pairs with |r|>0.7: {len(hi_pairs)}; of these, sharing >=1 gene: {sum(1 for a,b,_ in hi_pairs if shares(a,b))}")
for a,b,r in sorted(hi_pairs,key=lambda x:-abs(x[2])): print(f"       {a:14s} ~ {b:14s} r={r:+.2f}  {'SHARES GENES' if shares(a,b) else 'disjoint'}")
# recompute those correlations on disjoint gene sets only
rows=[]
for a,b,r in hi_pairs:
    ga=set(str(x).upper() for x in pct.get_genes(a)); gb=set(str(x).upper() for x in pct.get_genes(b))
    da,db=ga-gb,gb-ga
    if len(da)<2 or len(db)<2: rows.append((a,b,r,np.nan)); continue
    Rd=load_hpa_perceptivity(mode='single_cell_scaled')['R']
    # score disjoint sets directly from the reference gene-level matrix is unavailable here; use module minus shared proxy
    rows.append((a,b,r,np.nan))
print("     (disjoint-set recomputation requires gene-level reference scores; flagged as a production item)")

# ---------- M3: 6-D vs 9-D sensitivity ----------
d=np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True); S=d['S'].astype(float); cid=d['cid']; uc=sorted(set(cid.tolist()))
print("\n[M3] sensitivity of the core fraction and contraction to the number of retained components:")
for k in [6,9]:
    L=load[:,:k]; rp=((Rm.values-mu)/sd)@L; rpn=rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
    pc=((S-mu)/sd)@L; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@rpn.T; sim[:,EX]=-1e9; t=np.argmax(sim,axis=1)
    cf=100*np.isin(t,CORE3).mean()
    Cc=np.array([pcn[cid==c].mean(0) for c in uc]); Cc=Cc/np.linalg.norm(Cc,axis=1,keepdims=True)
    s_=Cc@Cc.T; iu2=np.triu_indices(len(Cc),1); contr=float((1-s_[iu2]).mean())
    rec=np.zeros(len(refcells))
    for c in uc:
        tt=t[cid==c]; rec[np.argsort(-np.bincount(tt,minlength=len(refcells)))[:5]]+=1
    rec[EX]=-1; top=[refcells[i] for i in np.argsort(-rec)[:3]]
    print(f"     {k}-D: core fraction {cf:5.1f}% | contraction {contr:.3f} | top recurrence {top}")

# ---------- M4: leave-one-cancer-out over all 25 ----------
print("\n[M4] leave-one-cancer-out (all 25):")
rp=((Rm.values-mu)/sd)@load; rpn=rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
pc=((S-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
sim=pcn@rpn.T; sim[:,EX]=-1e9; t=np.argmax(sim,axis=1)
res=[]
for drop in uc:
    k=cid!=drop; sub=sorted(set(cid[k].tolist()))
    cf=100*np.isin(t[k],CORE3).mean()
    rec=np.zeros(len(refcells))
    for c in sub:
        tt=t[cid==c]; rec[np.argsort(-np.bincount(tt,minlength=len(refcells)))[:5]]+=1
    rec[EX]=-1; order=np.argsort(-rec)
    cyto_rank=1+int(np.where(order==refcells.index('cytotrophoblasts'))[0][0])
    res.append((drop,cf,cyto_rank))
R=pd.DataFrame(res,columns=['dropped','core_pct','cyto_recurrence_rank'])
R.to_csv(ROOT/'results/referee2_newdata/E34_leave_one_out.csv',index=False)
print(f"     core fraction across the 25 leave-one-out panels: {R.core_pct.min():.1f}–{R.core_pct.max():.1f}% (full panel {100*np.isin(t,CORE3).mean():.1f}%)")
print(f"     cytotrophoblast recurrence rank: median {R.cyto_recurrence_rank.median():.0f}, worst {R.cyto_recurrence_rank.max()} — top-3 in {int((R.cyto_recurrence_rank<=3).sum())}/25 panels")
print(f"     most influential cohort (largest drop in core): {R.loc[R.core_pct.idxmin(),'dropped']}")
