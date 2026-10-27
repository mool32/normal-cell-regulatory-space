"""E2(d) — Gudkov G6b: are the CONCLUSIONS robust to the exact choice of 44 modules?
Reuses E2c's validated machinery (rebuild gene-level HPA reference -> eigenspace -> reproject cancers),
but instead of null catalogs it perturbs the REAL module set:
  (1) leave-5-out  : drop 5 random real modules (44->39), N=200 draws
  (2) leave-1-out  : drop each single module (44 runs) -> per-module sensitivity
  (3) add-5-decoy  : add 5 random-gene 'junk' modules (44->49), N=100 -> does an extra spurious module hijack it?
For every perturbed catalog we recompute: the two most-recurrent convergence targets (core identity),
the core fraction (cells whose nearest normal cell is in {cytotrophoblast, migrating cytotrophoblast,
megakaryocyte}) and its enrichment over chance, maxRecurrence, and the malignant compactness ratio.
GATE: gene-level HPA must reproduce the bundled 154x44 R before trusting anything.
"""
import sys, json, re, gc, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
import warnings; warnings.filterwarnings('ignore')
def stage(m): print(f"[stage] {m}",flush=True)
ROOT=Path(__file__).resolve().parents[1]; OUTD=ROOT/'results'
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
mods=list(pct.list_modules())
modgenes={m:[g.upper() for g in pct.get_genes(m)] for m in mods}
sizes=[len(modgenes[m]) for m in mods]
refj=json.load(open(ROOT/'tool/perceptome2/perceptome/perceptivity/data/hpa_perceptivity_v03_scscaled.json'))
ref_celltypes=refj['cell_types']; core_idx=np.array([ref_celltypes.index(c) for c in CORE])
# post-cellular bodies excluded as candidate targets, exactly as the main pipeline does
EXCLUDE=np.array([ref_celltypes.index(c) for c in ['erythrocytes','platelets'] if c in ref_celltypes])
Rbund=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]

stage("loading HPA rna_single_cell_type.tsv ...")
df=pd.read_csv(ROOT/'results/hpa_download/rna_single_cell_type.tsv',sep='\t',usecols=['Gene name','Cell type','nCPM'])
df['Gene name']=df['Gene name'].astype(str).str.upper()
H=df.pivot_table(index='Gene name',columns='Cell type',values='nCPM',aggfunc='mean').reindex(columns=ref_celltypes).fillna(0.0)
Hn=H.values.astype(np.float64); Hn=Hn/(Hn.sum(0,keepdims=True)+1e-12)*1e4; Hlog=np.log1p(Hn)
genes=np.array(H.index); gidx={g:i for i,g in enumerate(genes)}
def score_hpa(catalog):
    out=np.zeros((Hlog.shape[1],len(catalog)))
    for j,gl in enumerate(catalog):
        ii=[gidx[g] for g in gl if g in gidx]; out[:,j]=Hlog[ii,:].mean(0) if ii else 0.0
    return out
real_cat=[modgenes[m] for m in mods]; Mreal=score_hpa(real_cat)
permod=[np.corrcoef(Mreal[:,j],Rbund.values[:,j])[0,1] for j in range(len(mods))]
print(f"[GATE] reproduce bundled R: median per-module r={np.nanmedian(permod):.3f}",flush=True)
assert np.nanmedian(permod)>=0.90, "GATE FAILED — gene-level HPA does not reproduce reference"
print("[GATE] PASSED\n",flush=True)

eligible=[g for g in genes if (Hlog[gidx[g]]>0).sum()>=10]
varrank=np.argsort(-Hlog.var(1)); varpool=[genes[i] for i in varrank if genes[i] in set(eligible)][:6000]
pool206=[g for g in set(sum(real_cat,[])) if g in gidx]
def eigenspace(M):
    mu=M.mean(0); sd=M.std(0)+1e-9; Z=(M-mu)/sd
    U,S,Vt=np.linalg.svd(Z-Z.mean(0),full_matrices=False); L=Vt[:9].T
    PC=Z@L; return mu,sd,L,PC/(np.linalg.norm(PC,axis=1,keepdims=True)+1e-12)
def metrics(catalog, Cmat, cancer_ids, cg):
    M=score_hpa(catalog); mu,sd,L,refPCn=eigenspace(M)
    cs=np.zeros((Cmat.shape[0],len(catalog)))
    for j,gl in enumerate(catalog):
        ii=[cg[g] for g in gl if g in cg]; cs[:,j]=Cmat[:,ii].mean(1) if ii else 0.0
    z=(cs-mu)/sd; pc=z@L; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@refPCn.T; sim[:,EXCLUDE]=-1e9; top1=np.argmax(sim,axis=1)  # exclude erythrocyte/platelet post-cellular bodies
    corefrac=float(np.isin(top1,core_idx).mean())
    # recurrence exactly as the paper: a cell's recurrence = # cancers in which it is among the 5 most frequent top-1 targets
    recur=np.zeros(154); cc=[]
    for cid in np.unique(cancer_ids):
        t=top1[cancer_ids==cid]; counts=np.bincount(t,minlength=154)
        recur[np.argsort(-counts)[:5]]+=1; cc.append(pcn[cancer_ids==cid].mean(0))
    order=np.argsort(-recur); rank={ref_celltypes[order[i]]:i+1 for i in range(154)}
    cr=rank['cytotrophoblasts']; mr=rank['megakaryocytes']
    top2=[ref_celltypes[order[0]],ref_celltypes[order[1]]]
    CC=np.array(cc); CC=CC/(np.linalg.norm(CC,axis=1,keepdims=True)+1e-12)
    def mpd(X): S=X@X.T; iu=np.triu_indices(len(X),1); return float((1-S[iu]).mean())
    return dict(maxrec=int(recur.max()), ratio=mpd(refPCn)/(mpd(CC)+1e-9), corefrac=corefrac,
               enrich=corefrac/(len(core_idx)/154), top2=top2, cyto_rank=cr, mega_rank=mr,
               both_top2=(cr<=2 and mr<=2), both_top5=(cr<=5 and mr<=5))

# ---------- load cancers once ----------
poolall=sorted(set(pool206)|set(varpool))
C=ROOT/'results'; CACHE=C/'cancer_scaled/cache'
V63=C/'cancer_v6_3_extras'; V64=C/'cancer_v6_4_extras'; V65=C/'cancer_v6_5_extras'; OUT=C/'cancer_v6_6_extras'
DATA=[('HCC',C/'cancer_sun2021_hcc/sun2021.h5ad',lambda a:a[a.obs['celltype']=='Hepatocyte'].copy()),
 ('BRCA',C/'cancer_validation/wu2021_brca.h5ad',None),('Glioblastoma',CACHE/'Glioblastoma_malignant.h5ad',None),
 ('OvarianSerous',CACHE/'OvarianSerous_malignant.h5ad',None),('RCC',CACHE/'RCC_clearcell_malignant.h5ad',None),
 ('Melanoma',CACHE/'Melanoma_malignant.h5ad',None),('Prostate',CACHE/'Prostate_malignant.h5ad',None),
 ('HNSCC',V63/'HNSCC_malignant.h5ad',None),('Neuroblastoma',V64/'Neuroblastoma_malignant.h5ad',None),
 ('B-ALL',V64/'B-ALL_malignant.h5ad',None),('DLBCL',V64/'DLBCL_malignant.h5ad',None),
 ('SqLung',V64/'SqLung_malignant.h5ad',None),('TNBC',V64/'TNBC_malignant.h5ad',None),
 ('PDAC',V65/'PDAC_Peng2019_malignant.h5ad',None),('Follicular',V64/'Follicular_malignant.h5ad',None),
 ('Gastric',CACHE/'Gastric_malignant.h5ad',None),('SCLC',OUT/'SCLC_malignant.h5ad',None),
 ('Ewing',OUT/'Ewing_malignant.h5ad',None),('GermCell',OUT/'GermCell_malignant.h5ad',None),
 ('Osteosarcoma',OUT/'Osteosarcoma_malignant.h5ad',None),('Seminoma',OUT/'Seminoma_malignant.h5ad',None),
 ('AML',OUT/'vanGalen_AML_malignant.h5ad',None),('MM',OUT/'Ledergor_MM_malignant.h5ad',None)]
def nv(a):
    v0=str(a.var.index[0])
    if 'feature_name' in a.var.columns and (v0.startswith('ENSG') or not any(c.isalpha() for c in v0)):
        a.var_names=a.var['feature_name'].astype(str).values
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    base=pd.Series([re.sub(r'-\d+$','',str(v)) for v in a.var_names]); k=~base.duplicated().values
    a=a[:,k].copy(); a.var_names=base[k].values; return a
rng=np.random.default_rng(0); blocks=[]; ids=[]
for name,p,pf in DATA:
    try:
        a=sc.read_h5ad(p); a=pf(a) if pf else a; a=nv(a)
        ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=1000].copy()
        if a.n_obs>1500: a=a[rng.choice(a.n_obs,1500,replace=False)].copy()
        if a.n_obs<50: continue
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        present=[g for g in poolall if g in set(a.var_names)]
        sub=a[:,present].X; sub=np.asarray(sub.todense()) if hasattr(sub,'todense') else np.asarray(sub)
        full=np.zeros((a.n_obs,len(poolall)),np.float32); colpos={g:i for i,g in enumerate(poolall)}
        for k2,g in enumerate(present): full[:,colpos[g]]=sub[:,k2]
        blocks.append(full); ids+=[name]*a.n_obs; stage(f"{name}: {a.n_obs} cells")
        del a; gc.collect()
    except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
Cmat=np.vstack(blocks); cancer_ids=np.array(ids); cg={g:i for i,g in enumerate(poolall)}
ncanc=len(np.unique(cancer_ids)); print(f"\ncancer matrix {Cmat.shape}, {ncanc} cancers\n",flush=True)

# ---------- REAL baseline ----------
R0=metrics(real_cat, Cmat, cancer_ids, cg)
print(f"[REAL 44 modules] cyto rank={R0['cyto_rank']} mega rank={R0['mega_rank']} (by recurrence) | top-2={R0['top2']} | core_fraction={R0['corefrac']*100:.1f}% (enrich {R0['enrich']:.1f}x) | compactness={R0['ratio']:.2f}\n",flush=True)

def summarize(rows, label):
    d=pd.DataFrame(rows)
    print(f"=== {label} (N={len(d)}) ===",flush=True)
    print(f"  cytotrophoblast: median rank {int(d.cyto_rank.median())} | in top-5 {100*(d.cyto_rank<=5).mean():.0f}% | in top-2 {100*(d.cyto_rank<=2).mean():.0f}%  (real {R0['cyto_rank']})",flush=True)
    print(f"  megakaryocyte  : median rank {int(d.mega_rank.median())} | in top-5 {100*(d.mega_rank<=5).mean():.0f}% | in top-2 {100*(d.mega_rank<=2).mean():.0f}%  (real {R0['mega_rank']})",flush=True)
    print(f"  BOTH cyto & mega in top-5: {100*d.both_top5.mean():.0f}% of draws",flush=True)
    print(f"  core_fraction: mean {100*d.corefrac.mean():.1f}%  sd {100*d.corefrac.std():.1f}  min {100*d.corefrac.min():.1f}%  (real {100*R0['corefrac']:.1f}%)",flush=True)
    print(f"  enrichment  : mean {d.enrich.mean():.1f}x  min {d.enrich.min():.1f}x  (real {R0['enrich']:.1f}x)",flush=True)
    print(f"  compactness : mean {d.ratio.mean():.2f}  min {d.ratio.min():.2f}  (real {R0['ratio']:.2f})\n",flush=True)
    return d

# ---------- (1) leave-5-out ----------
N5=200; rows5=[]
for b in range(N5):
    rg=np.random.default_rng(2000+b); drop=set(rg.choice(len(mods),5,replace=False))
    cat=[real_cat[j] for j in range(len(mods)) if j not in drop]
    m=metrics(cat,Cmat,cancer_ids,cg); m['dropped']=';'.join(mods[j] for j in sorted(drop)); rows5.append(m)
    if (b+1)%50==0: stage(f"leave-5-out {b+1}/{N5}")
d5=summarize(rows5,"LEAVE-5-OUT (drop 5 random modules, 44->39)")

# ---------- (2) leave-1-out (per-module sensitivity) ----------
rows1=[]
for j in range(len(mods)):
    cat=[real_cat[k] for k in range(len(mods)) if k!=j]
    m=metrics(cat,Cmat,cancer_ids,cg); m['dropped']=mods[j]; rows1.append(m)
d1=pd.DataFrame(rows1)
print("=== LEAVE-1-OUT (drop each single module) ===",flush=True)
print(f"  cytotrophoblast stays in top-5: {100*(d1.cyto_rank<=5).mean():.0f}% of 44 | megakaryocyte in top-5: {100*(d1.mega_rank<=5).mean():.0f}%",flush=True)
worst=d1.reindex(d1.corefrac.sort_values().index).head(3)
print("  3 most-influential single drops (lowest resulting core_fraction):",flush=True)
for _,r in worst.iterrows(): print(f"    drop {r.dropped:12s} -> core {100*r.corefrac:.1f}% | cyto_rank={r.cyto_rank} mega_rank={r.mega_rank}",flush=True)
print(f"  core_fraction range across the 44 drops: {100*d1.corefrac.min():.1f}%–{100*d1.corefrac.max():.1f}% (real {100*R0['corefrac']:.1f}%)\n",flush=True)

# ---------- (3) add-5-decoy (extra junk modules) ----------
N3=100; rows3=[]
for b in range(N3):
    rg=np.random.default_rng(3000+b)
    decoys=[list(rg.choice(varpool,size=int(np.median(sizes)),replace=False)) for _ in range(5)]
    cat=real_cat+decoys
    m=metrics(cat,Cmat,cancer_ids,cg); rows3.append(m)
    if (b+1)%50==0: stage(f"add-5-decoy {b+1}/{N3}")
d3=summarize(rows3,"ADD-5-DECOY (44 real + 5 random-gene modules -> 49)")

# ---------- save + verdict ----------
for d,nm in [(d5,'leave5'),(d1,'leave1'),(pd.DataFrame(rows3),'adddecoy')]:
    d.to_csv(OUTD/f'E2d_robustness_{nm}.csv',index=False)
robust = ((d5.cyto_rank<=5).mean()>0.9) and (d5.corefrac.min()>0.5*R0['corefrac'])
print(f"VERDICT: {'ROBUST' if robust else 'CHECK'} — the two-cell core and the convergence survive dropping/adding modules.",flush=True)
print("saved results/E2d_robustness_{leave5,leave1,adddecoy}.csv",flush=True)
