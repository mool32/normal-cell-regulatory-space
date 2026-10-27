import sys, json, gc, numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
import warnings; warnings.filterwarnings('ignore')
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refPC=(Rm.values-mu)/sd@load; refcells=np.array(Rm.index.tolist()); rn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)
DUET=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
# backed: read obs cheaply, load only epithelial rows
b=sc.read_h5ad('results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',backed='r')
cm=b.obs['ClusterMidway'].astype(str).values
mask=np.isin(cm,['EpiT','Epi'])
print(f"epithelial cells: {mask.sum()} of {b.n_obs}",flush=True)
a=b[mask].to_memory(); del b; gc.collect()
if 'feature_name' in a.var.columns: a.var_names=a.var['feature_name'].astype(str)
a.var_names=[str(g).upper() for g in a.var_names]; a.var_names_make_unique()
sc.pp.filter_cells(a,min_genes=1000)
arm=np.where(a.obs['ClusterMidway'].astype(str).values=='EpiT','malig','normal')
pcol='donor_id' if 'donor_id' in a.obs else 'PatientTypeID'
pat=a.obs[pcol].astype(str).values
if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0)
del a; gc.collect()
cp=((S.values-mu)/sd)@load; cn=cp/(np.linalg.norm(cp,axis=1,keepdims=True)+1e-12)
top1=refcells[(cn@rn.T).argmax(1)]; isd=np.array([t in DUET for t in top1])
nm=arm=='normal'
print("normal-epi maps to (top6):", {k:f'{100*v/nm.sum():.0f}%' for k,v in pd.Series(top1[nm]).value_counts().head(6).items()},flush=True)
print(f"normal-epi eryth-default: {100*np.mean(top1[nm]=='erythrocytes'):.1f}%",flush=True)
rows=[]
for pp in sorted(set(pat)):
    m=(pat==pp)&(arm=='malig'); n=(pat==pp)&(arm=='normal')
    if m.sum()>=20 and n.sum()>=20: rows.append((pp,m.sum(),100*isd[m].mean(),n.sum(),100*isd[n].mean()))
df=pd.DataFrame(rows,columns=['pt','mn','mp','nn','np'])
print(f"\nper-patient (n>=20 each), patients with both arms: {len(df)}",flush=True)
for _,r in df.iterrows(): print(f"  {str(r.pt)[:18]:18s} malig {int(r.mn):5d} {r.mp:5.1f}%  normal {int(r.nn):5d} {r['np']:5.1f}%  ratio {r.mp/max(r['np'],0.1):.1f}",flush=True)
print(f"\nnormal-epi duet%: median {df['np'].median():.1f} range [{df['np'].min():.1f},{df['np'].max():.1f}]",flush=True)
print(f"malig-epi  duet%: median {df['mp'].median():.1f} range [{df['mp'].min():.1f},{df['mp'].max():.1f}]",flush=True)
print(f"malignant>normal in {(df.mp>df['np']).sum()}/{len(df)} patients",flush=True)
from scipy.stats import wilcoxon
if len(df)>=6: w,pv=wilcoxon(df.mp,df['np']); print(f"paired Wilcoxon malignant>normal: p={pv:.5f}",flush=True)
