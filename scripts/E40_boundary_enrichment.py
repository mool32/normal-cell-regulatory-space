"""N6 fix + data for Fig. 3F — recompute the three convergence-region boundaries ONCE, on the panel the
paper actually reports (83,819 malignant cells, 25 cohorts, 152 candidates after the erythrocyte/platelet
quarantine), so that Table 1, Fig. 3F and the main text can all quote the same numbers.
Boundaries are defined by cross-cancer recurrence, not by proximity:
  core           = the 2 cells leading both recurrence and depth
  neighbourhood  = the cells recurrent in a MAJORITY of cancers (>=13/25)
  broader region = the cells recurrent in >=5 of 25 cancers"""
import json, numpy as np, pandas as pd, sys
sys.path.insert(0,'tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve()
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; cells=list(Rm.index)
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refP=((Rm.values-mu)/sd)@load; refn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
EX=[cells.index(c) for c in ('erythrocytes','platelets')]
d=np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True); S=d['S'].astype(float); cid=d['cid']
uc=sorted(set(cid.tolist()))
pc=((S-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
sim=pcn@refn.T; sim[:,EX]=-1e9; top=np.argmax(sim,axis=1)
rec=np.zeros(len(cells))
for c in uc:
    t=top[cid==c]; b=np.bincount(t,minlength=len(cells)); b[EX]=-1
    rec[np.argsort(-b)[:5]]+=1
rec[EX]=-1
order=np.argsort(-rec)
NCAND=len(cells)-len(EX)
print(f"panel: {len(S):,} malignant cells, {len(uc)} cohorts, {NCAND} candidate cell types\n")
print("recurrence ranking (top 16):")
for i in order[:16]: print(f"   {cells[i]:34s} {int(rec[i]):2d}/25   occupancy {100*np.mean(top==i):5.2f}%")

core=list(order[:2]); nbh=list(order[:4]); brd=list(order[:14])
rows=[]
for nm,idx,rule in [('core (top-2)',core,'the two cells leading the cross-cancer recurrence ranking'),
                    ('neighbourhood (top-4)',nbh,'the four most cross-cancer-recurrent cell types'),
                    ('broader region (top-14)',brd,'the fourteen most cross-cancer-recurrent cell types')]:
    occ=100*np.isin(top,idx).mean(); ch=100*len(idx)/NCAND
    rows.append(dict(boundary=nm,k=len(idx),rule=rule,occupancy_pct=round(occ,1),chance_pct=round(ch,2),
                     fold=round(occ/ch,1),members='; '.join(cells[i] for i in idx)))
    print(f"\n{nm:16s} k={len(idx):2d}  occupancy {occ:5.1f}%  chance {ch:4.2f}%  enrichment {occ/ch:5.1f}x")
    print(f"   members: {', '.join(cells[i] for i in idx)}")
T=pd.DataFrame(rows); T.to_csv(ROOT/'results/cancer_v6_6_extras/table1_boundaries.csv',index=False)
per=[]
for c in uc:
    t=top[cid==c]
    per.append(dict(cancer=c,n_cells=int(len(t)),
                    core=round(100*np.isin(t,core).mean(),1),
                    neighbourhood=round(100*np.isin(t,nbh).mean(),1),
                    broader=round(100*np.isin(t,brd).mean(),1)))
P=pd.DataFrame(per).sort_values('core',ascending=False)
P.to_csv(ROOT/'results/cancer_v6_6_extras/table1_per_cancer.csv',index=False)
print(f"\nper-cancer table written ({len(P)} cohorts); cohorts above the 2-cell chance level: "
      f"{int((P.core>100*2/NCAND).sum())}/{len(P)}")
print("DONE")
