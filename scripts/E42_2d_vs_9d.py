"""Recompute the 2-D-versus-9-D comparator quoted in Note S1 (Fig. S15) on the CURRENT panel.

The figure in the SI and Note S1 carried '~7% versus 56%' and the Fig. S5 legend carried '40%' —
three numbers for one quantity, inherited from an older 19-cohort run and from the panel as it stood
before scripts/E41_regen_ranking.py fixed the quarantine convention. This recomputes both arms from
the one cached matrix, under the one convention, so the SI can quote a matched pair.
"""
import json, numpy as np, pandas as pd, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,(ROOT/'tool/perceptome2').as_posix())
from perceptome.perceptivity.reference import load_hpa_perceptivity
e=json.load(open(ROOT/'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; cells=list(Rm.index)
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refP=((Rm.values-mu)/sd)@load
EX=[cells.index(c) for c in ('erythrocytes','platelets')]
CORE2=[cells.index(c) for c in ('cytotrophoblasts','megakaryocytes')]
d=np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True); S=d['S'].astype(float); cid=d['cid']
pc=((S-mu)/sd)@load
def corefrac(k, sel):
    r=refP[:,:k]; r=r/(np.linalg.norm(r,axis=1,keepdims=True)+1e-12)
    q=pc[sel][:,:k]; q=q/(np.linalg.norm(q,axis=1,keepdims=True)+1e-12)
    sim=q@r.T; sim[:,EX]=-1e9
    return 100*np.isin(np.argmax(sim,axis=1),CORE2).mean()
rows=[]
for c in ['B-ALL','Neuroblastoma']:
    sel=cid==c
    rows.append(dict(cancer=c,n=int(sel.sum()),core_2d=round(corefrac(2,sel),1),core_9d=round(corefrac(9,sel),1)))
sel=np.ones(len(cid),bool)
rows.append(dict(cancer='ALL 25 (pooled)',n=int(sel.sum()),core_2d=round(corefrac(2,sel),1),core_9d=round(corefrac(9,sel),1)))
D=pd.DataFrame(rows); D.to_csv(ROOT/'results/cancer_v6_6_extras/E42_2d_vs_9d.csv',index=False)
print(D.to_string(index=False))
m=pd.read_csv(ROOT/'results/cancer_v6_6_extras/panel_25_top1_matrix.csv',index_col=0)
print('\ncross-check against panel_25_top1_matrix.csv (two-cell core):')
for c in ['B-ALL','Neuroblastoma']:
    print(f"   {c:15s} {m.loc[['cytotrophoblasts','megakaryocytes'],c].sum():5.1f}%  "
          f"(cyto {m.loc['cytotrophoblasts',c]:.1f} + mega {m.loc['megakaryocytes',c]:.1f})")
