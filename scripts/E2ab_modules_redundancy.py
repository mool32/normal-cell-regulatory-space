"""E2 (a) clear module definition + (b) redundancy — Gudkov #2, parts a/b.
(a) Emit the explicit 44-module -> core-gene table (what is in the catalog).
(b) Redundancy: pairwise correlation of the 44 modules across the 154 HPA cell types (the matrix the
    eigenspace is built from); list and visualize pairs with |r|>0.7. Few, biologically expected => no
    hidden dependence driving the geometry.
"""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'; OUT=ROOT+'/results/cancer_v6_6_extras'
mods=list(pct.list_modules())

# (a) module -> genes table
rows=[]
for m in mods:
    g=list(pct.get_genes(m))
    rows.append({'module':m,'n_genes':len(g),'core_genes':', '.join(g)})
T=pd.DataFrame(rows).sort_values('module')
T.to_csv(OUT+'/tableS9_module_genes.csv',index=False)
allg=set(); [allg.update(pct.get_genes(m)) for m in mods]
print(f"(a) {len(mods)} modules, {len(allg)} unique core genes, sizes {T.n_genes.min()}-{T.n_genes.max()} (median {int(T.n_genes.median())})")
print("    saved tableS9_module_genes.csv")

# (b) redundancy across the 154 HPA cell types (single_cell_scaled, the eigenspace matrix)
R=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]   # 154 x 44
C=np.corrcoef(R.values.T)
pairs=[(mods[i],mods[j],C[i,j]) for i in range(len(mods)) for j in range(i+1,len(mods)) if abs(C[i,j])>0.7]
pairs.sort(key=lambda x:-abs(x[2]))
print(f"\n(b) redundancy: {len(pairs)} of {len(mods)*(len(mods)-1)//2} module pairs with |r|>0.7:")
for a,b,r in pairs: print(f"    {a:12s} ~ {b:12s} r={r:+.2f}")
offdiag=C[~np.eye(len(mods),dtype=bool)]
print(f"    mean |r| off-diagonal = {np.abs(offdiag).mean():.3f}; 95th pct |r| = {np.percentile(np.abs(offdiag),95):.2f}")
pd.DataFrame(C,index=mods,columns=mods).to_csv(OUT+'/module_correlation_matrix.csv')

# heatmap (clustered by hand: keep module order)
order=np.argsort([R[m].values.tolist() for m in mods],axis=0)  # not used; keep alpha order
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 1.15))
im=ax.imshow(C,cmap='RdBu_r',vmin=-1,vmax=1)
ax.set_xticks(range(len(mods))); ax.set_xticklabels(mods,rotation=90,fontsize=BASE*0.825)
ax.set_yticks(range(len(mods))); ax.set_yticklabels(mods,fontsize=BASE*0.825)
# mark the |r|>0.7 off-diagonal pairs
for a,b,r in pairs:
    i,j=mods.index(a),mods.index(b)
    for (x,y) in [(i,j),(j,i)]:
        ax.add_patch(plt.Rectangle((x-0.5,y-0.5),1,1,fill=False,edgecolor='black',lw=1.2))
cb=fig.colorbar(im,fraction=0.046,pad=0.04); cb.set_label('Pearson r across 154 normal cell types',fontsize=BASE*1.057)
plt.tight_layout(); plt.savefig(FIG+'/v7_fig_E2_redundancy.png',dpi=170,bbox_inches='tight'); plt.close()
print("    saved v7_fig_E2_redundancy.png + module_correlation_matrix.csv")
