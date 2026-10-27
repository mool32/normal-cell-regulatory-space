"""QC-review fixes N4, N11, N12 — recompute Table S4 on all 25 cohorts, complete Tables S2/S6 to 154 rows
with explicit QC flags, and report the four cell-count denominators so Methods can define them."""
import json, numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve(); T=ROOT/'results/cancer_v6_6_extras'
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
Z=(Rm.values-mu)/sd                                   # 154 x 44 standardized reference
core2=[refcells.index('cytotrophoblasts'),refcells.index('megakaryocytes')]
u=Z[core2].mean(0); u=u/np.linalg.norm(u)             # two-component convergent signature (unit)

# ---------------- N4: Table S4 recomputed on all 25 cohorts ----------------
d=np.load(ROOT/'results/E7_S_raw.npz',allow_pickle=True); S=d['S'].astype(float); cid=d['cid']
uc=sorted(set(cid.tolist()))
C=np.array([((S[cid==c].mean(0))-mu)/sd for c in uc])          # 25 x 44 standardized cancer centroids
C=C/np.linalg.norm(C,axis=1,keepdims=True)
contrib=C*u                                                     # per-module terms of the cosine
S4=pd.DataFrame({'mean_contribution':contrib.mean(0).round(4),
                 'sd':contrib.std(0).round(4),
                 'pct_positive':(contrib>0).mean(0).round(3),
                 'n_cohorts':len(uc)},index=mods).sort_values('mean_contribution',ascending=False)
old=pd.read_csv(T/'tableS4_modules.csv',index_col=0)
S4.to_csv(T/'tableS4_modules.csv')
print(f"[N4] Table S4 recomputed on {len(uc)} cohorts (was implicitly 16: all old pct_positive were k/16)")
print("     old pct_positive values:",sorted(set(old.pct_positive))[:6],"...")
print("     new pct_positive values:",sorted(set(S4.pct_positive))[:6],"...  (k/25 grid)")
print("     top 5 modules by contribution:",list(S4.index[:5]))

# ---------------- N11: Tables S2 and S6 completed to 154 rows with QC flags ----------------
QUAR={'erythrocytes','platelets'}
s2=pd.read_csv(T/'tableS2_ranking154.csv',index_col=0)
missing=[c for c in refcells if c not in s2.index]
for c in missing: s2.loc[c]=[0,0.0,0.0]
s2['qc_status']=['quarantined (post-cellular, low-complexity)' if c in QUAR else 'candidate' for c in s2.index]
s2=s2.reindex(refcells); s2.to_csv(T/'tableS2_ranking154.csv')
print(f"\n[N11] Table S2: {len(s2)} rows (added: {missing}); qc_status column added, {int((s2.qc_status!='candidate').sum())} quarantined")
s6=pd.read_csv(T/'tableS6_complexity.csv',index_col=0)
miss6=[c for c in refcells if c not in s6.index]
for c in miss6: s6.loc[c]=[0,np.nan,np.nan]
s6['qc_status']=['quarantined (post-cellular, low-complexity)' if c in QUAR else 'candidate' for c in s6.index]
s6['note']=['no cells assigned' if s6.loc[c,'n_cells_assigned']==0 else '' for c in s6.index]
s6=s6.reindex(refcells); s6.to_csv(T/'tableS6_complexity.csv')
print(f"[N11] Table S6: {len(s6)} rows (added {len(miss6)} zero-assignment types); qc_status + note columns added")

# ---------------- N12: the four cell-count denominators ----------------
s1=pd.read_csv(T/'tableS1_datasets.csv')
print("\n[N12] denominators:")
print(f"   (a) full QC pool across the 25 cohorts (Table S1 sum of n_cells_QC): {int(s1.n_cells_QC.sum()):,}")
print(f"   (b) projection subset actually scored (cap 4,000 cells/cohort, n_genes>=1000): {S.shape[0]:,} over {len(uc)} cohorts")
print(f"   (c) assignments tabulated in Table S6 (includes the normal-atlas projections, not only cancer): {int(s6.n_cells_assigned.sum()):,}")
print(f"   (d) subset with per-cell gene counts used in the complexity analyses (Fig. S7): ~312,000")
