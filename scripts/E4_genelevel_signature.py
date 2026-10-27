"""Gene-level rescue of the released-brake signature (Paper-8 exposure). The mechanism = 'cancer reproduces the
cytotrophoblast's suppression signature'. Test it WITHOUT module grouping: over the 206 catalog genes, does per-gene
target-distinctiveness (cyto+mega vs rest, per gene) correlate with how much cancer suppresses that gene? If strong,
the mechanism is gene-level robust and module-incoherence (P9) is irrelevant to it. Split by gold/non-gold module
membership to show even flagged-module genes carry the signal.
"""
import sys, json, re, gc, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
from scipy.stats import spearmanr
sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
import warnings; warnings.filterwarnings('ignore')
def stage(m): print(f"[stage] {m}",flush=True)
ROOT=Path(__file__).resolve().parents[1]   # was a hardcoded path that no longer exists, which is why
                                           # 11 of the 17 cohorts below failed to open in results/E4_out.txt
ALL=list(pct.list_modules()); MG={m:[g.upper() for g in pct.get_genes(m)] for m in ALL}
GOLD=set(['Autophagy','AR','BMP','Circadian','FXR','HIF','Hippo','NF-κB','NRF2','p53','PPARα','SREBP','TGF-β','UPR-ATF6','UPR-IRE1','UPR-PERK'])
gene_gold={};
for m in ALL:
    for g in MG[m]: gene_gold[g.upper()]=gene_gold.get(g.upper(),False) or (m in GOLD)
refj=json.load(open(ROOT/'tool/perceptome2/perceptome/perceptivity/data/hpa_perceptivity_v03_scscaled.json'))
ct=np.array(refj['cell_types']); duet=np.array([c in {'cytotrophoblasts','megakaryocytes'} for c in ct])
df=pd.read_csv(ROOT/'results/hpa_download/rna_single_cell_type.tsv',sep='\t',usecols=['Gene name','Cell type','nCPM'])
df['Gene name']=df['Gene name'].str.upper()
H=df.pivot_table(index='Gene name',columns='Cell type',values='nCPM',aggfunc='mean').reindex(columns=ct).fillna(0.0)
Hn=H.values/(H.values.sum(0,keepdims=True)+1e-12)*1e4; Hlog=np.log1p(Hn); genes=np.array(H.index); gpos={g:i for i,g in enumerate(genes)}
catgenes=[g for g in sorted(set(sum([MG[m] for m in ALL],[]))) if g in gpos]
gi=np.array([gpos[g] for g in catgenes])
Hg=Hlog[gi]                                   # 206 x 154
gmean=Hg.mean(1,keepdims=True); gstd=Hg.std(1,keepdims=True)+1e-9
zg=(Hg-gmean)/gstd                            # per-gene z across 154
tdist=zg[:,duet].mean(1)-zg[:,~duet].mean(1)  # per-gene target-distinctiveness (neg=target-suppressed)
# cancers
CACHE=ROOT/'results/cancer_scaled/cache'; C=ROOT/'results'
V63=C/'cancer_v6_3_extras'; V64=C/'cancer_v6_4_extras'; V65=C/'cancer_v6_5_extras'; OUT=C/'cancer_v6_6_extras'
# The per-cohort filter is an OBS MASK, not an AnnData->AnnData function: a backed view cannot be
# .copy()'d (anndata asks for a filename), which is how Sun HCC failed once it went down the lazy path.
DATA=[('HCC',C/'cancer_sun2021_hcc/sun2021.h5ad',lambda o:(o['celltype']=='Hepatocyte').values),
 ('BRCA',C/'cancer_validation/wu2021_brca.h5ad',None),('Glioblastoma',CACHE/'Glioblastoma_malignant.h5ad',None),
 ('OvarianSerous',CACHE/'OvarianSerous_malignant.h5ad',None),('RCC',CACHE/'RCC_clearcell_malignant.h5ad',None),
 ('Melanoma',CACHE/'Melanoma_malignant.h5ad',None),('Prostate',CACHE/'Prostate_malignant.h5ad',None),
 ('HNSCC',V63/'HNSCC_malignant.h5ad',None),('Neuroblastoma',V64/'Neuroblastoma_malignant.h5ad',None),
 ('B-ALL',V64/'B-ALL_malignant.h5ad',None),('DLBCL',V64/'DLBCL_malignant.h5ad',None),
 ('SqLung',V64/'SqLung_malignant.h5ad',None),('TNBC',V64/'TNBC_malignant.h5ad',None),
 ('PDAC',V65/'PDAC_Peng2019_malignant.h5ad',None),
 ('NSCLC',ROOT/'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',lambda o:(o['cell_type'].astype(str)=='malignant cell').values),
 ('AML',OUT/'vanGalen_AML_malignant.h5ad',None),('MM',OUT/'Ledergor_MM_malignant.h5ad',None)]
def nv_plan(a):
    """The var-name normalisation of nv(), as (names, keep mask), computed WITHOUT touching X.

    Split out so a cohort too large to hold in memory can be filtered and subsampled while still
    backed on disk and materialised only once, at <=1500 cells. The logic is unchanged: it depends
    on var names alone, which is what makes the lazy path exact rather than approximate."""
    v = a.var.copy()
    if 'feature_name' in v.columns and (str(v.index[0]).startswith('ENSG')
                                        or not any(c.isalpha() for c in str(v.index[0]))):
        v.index = pd.Index(v['feature_name'].astype(str).values)
    v.index = v.index.str.upper()
    d = ad.AnnData(var=v, shape=(0, v.shape[0])); d.var_names_make_unique()
    base = pd.Series([re.sub(r'-\d+$', '', str(x)) for x in d.var_names])
    k = ~base.duplicated().values
    return base[k].values, k


def nv(a):
    names, k = nv_plan(a)
    a = a[:, k].copy(); a.var_names = names; return a


# A dense .h5ad is read into memory whole: Ledergor_MM is 34,176 x 57,874 float32 = 7.4 GB, and
# nv()'s column subset then doubles that, which is why the stored run never got past it on a 16 GB
# machine. Such a cohort is read backed and only the sampled cells are materialised. Sparse cohorts
# take the eager path: the largest here is 2 GB in memory, and backed sparse slicing in this anndata
# materialises the whole matrix anyway, so going backed would cost rather than save.
BACKED_GB = 1.0


def dense_gb(path):
    """In-memory GB if X is a dense HDF5 dataset, else 0.

    The file size is the wrong proxy: Sun HCC is 71,915 x 25,712 float32 and compresses to 0.31 GB
    on disk, so it would read as small and then need the same 7.4 GB as Ledergor MM."""
    import h5py
    with h5py.File(path, 'r') as f:
        X = f['X']
        if not isinstance(X, h5py.Dataset):
            return 0.0
        return X.shape[0] * X.shape[1] * X.dtype.itemsize / 2**30


def read_cohort(path, pf):
    """Load a cohort already QC-filtered to >=1000 detected genes and subsampled to <=1500 cells.

    Identical by construction to read -> pf -> nv -> ng>=1000 -> rng(0).choice(1500): the var-name
    normalisation is computed from var alone, the detected-gene count over exactly the columns nv()
    keeps, and the sampled cells are returned in the order the seeded generator drew them. Asserted
    against the eager path in `scripts/verify/check_E4_lazy.py`.

    The lazy branch never builds an intermediate view: a backed AnnData cannot be indexed twice
    ("you cannot make a view of a view"), so the cohort filter, the detected-gene filter and the
    subsample are composed into ONE index array and applied once."""
    lazy = dense_gb(path) >= BACKED_GB
    a = sc.read_h5ad(path, backed='r') if lazy else sc.read_h5ad(path)
    if not lazy:
        if pf is not None: a = a[pf(a.obs)].copy()
        names, k = nv_plan(a)
        a = a[:, k].copy(); a.var_names = names
        ng = np.asarray((a.X > 0).sum(1)).ravel()
        a = a[ng >= 1000].copy()
        return a, a.n_obs

    names, k = nv_plan(a)
    sel = np.where(np.asarray(pf(a.obs)))[0] if pf is not None else np.arange(a.n_obs)
    ng = np.empty(len(sel), dtype=np.int32)
    for i in range(0, len(sel), 2000):                     # peak RAM is one chunk, not the matrix
        blk = np.asarray(a.X[sel[i:i + 2000], :])[:, k]
        ng[i:i + 2000] = (blk > 0).sum(1)
        del blk
    keep = sel[ng >= 1000]
    n = len(keep)
    order = None
    if n > 1500:
        draw = np.random.default_rng(0).choice(n, 1500, replace=False)
        # h5py fancy indexing needs increasing indices, so read in file order and then restore the
        # order the draw produced: the eager path keeps the cells in draw order and the two objects
        # must agree row for row, not merely as sets.
        order = np.argsort(draw)
        keep = keep[draw[order]]
    a = a[keep].to_memory()
    if order is not None: a = a[np.argsort(order)].copy()
    a = a[:, k].copy(); a.var_names = names
    return a, n


percell=[]  # cancer mean-z per gene, per cancer
for name,p,pf in DATA:
    try:
        a, n_pass = read_cohort(p, pf)
        if a.n_obs > 1500:                      # eager path only; the lazy one sampled on the way in
            a = a[np.random.default_rng(0).choice(a.n_obs, 1500, replace=False)].copy()
        if a.n_obs < 50: continue
        print(f'  {name:14s} {n_pass:6d} cells >=1000 genes -> {a.n_obs} used', flush=True)
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        vn=set(a.var_names); present=[g for g in catgenes if g in vn]
        sub=a[:,present].X; sub=np.asarray(sub.todense()) if hasattr(sub,'todense') else np.asarray(sub)
        col={g:k for k,g in enumerate(present)}
        cm=np.full(len(catgenes),np.nan)
        for k,g in enumerate(catgenes):
            if g in col:
                j=col[g]; cm[k]=((sub[:,j]-gmean[k,0])/gstd[k,0]).mean()
        percell.append(cm); del a; gc.collect()
    except Exception as ex: print(f"  !!{name}:{ex}",flush=True)
CZ=np.nanmean(np.array(percell),0)            # cancer mean-z per gene (across cancers)
ok=~np.isnan(CZ)
rho,pp=spearmanr(tdist[ok],CZ[ok]); r=np.corrcoef(tdist[ok],CZ[ok])[0,1]
print(f"\n=== GENE-LEVEL signature (206 catalog genes, module-grouping-independent) ===",flush=True)
print(f"  Spearman(per-gene target-distinctiveness, cancer suppression) rho={rho:+.2f}  Pearson r={r:+.2f}  p={pp:.1e}  (n={ok.sum()} genes)",flush=True)
# split gold vs non-gold module genes
gold_mask=np.array([gene_gold.get(g,False) for g in catgenes])&ok
nong=(~np.array([gene_gold.get(g,False) for g in catgenes]))&ok
for nm,msk in [('genes from GOLD-16 modules',gold_mask),('genes from NON-gold modules',nong)]:
    if msk.sum()>5:
        rr,ppp=spearmanr(tdist[msk],CZ[msk]); print(f"    {nm:32s} (n={msk.sum():3d}): rho={rr:+.2f} p={ppp:.1e}",flush=True)
# permutation null
rng=np.random.default_rng(0); null=[spearmanr(tdist[ok],rng.permutation(CZ[ok]))[0] for _ in range(2000)]
print(f"  permutation null: |rho| exceeds real in {100*np.mean(np.abs(null)>=abs(rho)):.2f}% of 2000 shuffles",flush=True)
pd.DataFrame({'gene':np.array(catgenes)[ok],'target_distinct':tdist[ok],'cancer_z':CZ[ok],'gold':np.array([gene_gold.get(g,False) for g in catgenes])[ok]}).to_csv(ROOT/'results/E4_genelevel.csv',index=False)
print(f"\nVERDICT: gene-level signature {'ROBUST — mechanism survives independent of module grouping' if rho>0.4 and pp<0.01 else 'weak — does not rescue'}.",flush=True)
print("saved results/E4_genelevel.csv",flush=True)
