"""Single source of truth for the recurrence-and-depth ranking (review blocker #6).
Fig. 3B, Fig. 3C, Fig. 3F, Fig. S12, Table 1 and Table S2 were computed under three different conventions
for the quarantined post-cellular bodies, and disagreed with one another for seven cell types. This script
fixes ONE convention — the one the paper documents in Note S-K, "152 candidates after the erythrocyte /
platelet quarantine", i.e. the quarantined types are not candidates at all and are removed BEFORE the
nearest-normal assignment — and regenerates the ranking table from the same 83,819-cell panel everything
else uses. Table S2 is rewritten from this; the quarantined pair is reported separately, pre-quarantine,
so nothing is hidden."""
import json, numpy as np, pandas as pd, sys
sys.path.insert(0,'tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve(); T=ROOT/'results/cancer_v6_6_extras'
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; cells=list(Rm.index)
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refn=((Rm.values-mu)/sd)@load; refn=refn/np.linalg.norm(refn,axis=1,keepdims=True)
QUAR=['erythrocytes','platelets']; EX=[cells.index(c) for c in QUAR]
d=np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True); S=d['S'].astype(float); cid=d['cid']; uc=sorted(set(cid.tolist()))
pc=((S-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
sim=pcn@refn.T

def ranking(quarantine_before):
    s=sim.copy()
    if quarantine_before: s[:,EX]=-1e9
    top=np.argmax(s,axis=1)
    rec=np.zeros(len(cells)); depth=[[] for _ in cells]
    for c in uc:
        t=top[cid==c]; b=np.bincount(t,minlength=len(cells))
        if not quarantine_before: b[EX]=-1
        for i in np.argsort(-b)[:5]: rec[i]+=1
        for i in range(len(cells)): depth[i].append(100*np.mean(t==i))
    return rec,top,depth

rec,top,depth=ranking(True)
rec_pre,top_pre,_=ranking(False)
rows=[]
for i,c in enumerate(cells):
    q=c in QUAR
    rows.append(dict(cell=c,
        recurrence_top5=('' if q else int(rec[i])),
        median_depth_pct=('' if q else round(float(np.median(depth[i])),2)),
        mean_depth_pct=('' if q else round(float(np.mean(depth[i])),2)),
        occupancy_pct=('' if q else round(100*float(np.mean(top==i)),2)),
        qc_status='quarantined (post-cellular, low-complexity)' if q else 'candidate',
        prequarantine_occupancy_pct=round(100*float(np.mean(top_pre==i)),2)))
D=pd.DataFrame(rows)
D['_r']=D.recurrence_top5.replace('',-1).astype(float)
D=D.sort_values(['_r','median_depth_pct'],ascending=False).drop(columns='_r')
D.to_csv(T/'tableS2_ranking154.csv',index=False)
print(f"rewrote tableS2_ranking154.csv ({len(D)} rows, {int((D.qc_status=='candidate').sum())} candidates)")
print("\ntop of the single-convention ranking:")
print(D.head(9)[['cell','recurrence_top5','median_depth_pct','occupancy_pct']].to_string(index=False))
print("\nquarantined pair, pre-quarantine occupancy (reported, not hidden):")
print(D[D.qc_status!='candidate'][['cell','prequarantine_occupancy_pct']].to_string(index=False))
# boundaries from the SAME ranking
order=[i for i in np.argsort(-rec) if cells[i] not in QUAR]
NC=len(cells)-len(QUAR); out=[]
for nm,k in [('core (top-2)',2),('neighbourhood (top-4)',4),('broader region (top-14)',14)]:
    idx=order[:k]; occ=100*np.isin(top,idx).mean(); ch=100*k/NC
    out.append(dict(boundary=nm,k=k,occupancy_pct=round(occ,1),chance_pct=round(ch,2),fold=round(occ/ch,1),
                    members='; '.join(cells[i] for i in idx)))
    print(f"\n{nm:24s} occupancy {occ:5.1f}%  chance {ch:4.2f}%  {occ/ch:5.1f}x")
    print(f"   {', '.join(cells[i] for i in idx)}")
pd.DataFrame(out).to_csv(T/'table1_boundaries.csv',index=False)

# ---- regenerate the per-cancer top-1 matrix from the SAME assignment, so Fig. 3C / Fig. S12 / the
#      heatmap derive from one object rather than from a pre-computed matrix built under another rule ----
cand=[i for i,c in enumerate(cells) if c not in QUAR]
Mt=pd.DataFrame(0.0,index=[cells[i] for i in cand],columns=uc)
for c in uc:
    t_=top[cid==c]
    for i in cand: Mt.loc[cells[i],c]=100*float(np.mean(t_==i))
Mt.to_csv(T/'panel_25_top1_matrix.csv')
chk={}
for c in uc:
    for i in Mt[c].sort_values(ascending=False).head(5).index: chk[i]=chk.get(i,0)+1
print("\nregenerated panel_25_top1_matrix.csv; its recurrence now matches the table:")
for k,v in sorted(chk.items(),key=lambda x:-x[1])[:6]: print(f"   {k:28s} {v}/25")

print("\nDONE — ranking, boundaries and the per-cancer matrix now share one convention")
