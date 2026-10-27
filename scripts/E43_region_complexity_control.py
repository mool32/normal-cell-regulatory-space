"""Referee round 4, §2.2 + §2.4: run the transcriptome-complexity control on EVERY member of the
convergence region, not only on the two core cells and the erythrocyte, and audit gene-symbol
recovery per module per dataset while the data are loaded.

Why: Table S6 reports the median detected genes of the cells assigned to each target. Five members of
the top-14 region sit at or below the platelet (1,541) and three at or below the erythrocyte (890) that
is quarantined precisely for low complexity. Fig. 3D's decile control was only ever run for the
cytotrophoblast, the megakaryocyte and the erythrocyte. This runs it for all of them, with the same QC,
cap, seed and assignment rule as scripts/E1_megakaryocyte_complexity.py, so the series are comparable.

Outputs (results/cancer_v6_6_extras/):
  E43_decile_by_target.csv    % of malignant cells assigned to each tracked target, per complexity decile
  E43_target_verdicts.csv     per target: rho with decile, first/last decile, verdict
  E43_boundaries_clean.csv    boundary enrichments after dropping targets that fall with complexity
  E43_gene_coverage.csv       per dataset x module: genes in the catalog, found, missing
"""
import sys, gc, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,(R/'tool/perceptome2').as_posix())
import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity

e=json.load(open(R/'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refP=((Rm.values-mu)/sd)@load; refn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
idx=np.array(Rm.index.tolist())

C=R/'results'; CACHE=C/'cancer_scaled/cache'
V63=C/'cancer_v6_3_extras'; V64=C/'cancer_v6_4_extras'; V65=C/'cancer_v6_5_extras'; OUT=C/'cancer_v6_6_extras'
DATA=[('HCC',C/'cancer_sun2021_hcc/sun2021.h5ad',lambda a:a[a.obs['celltype']=='Hepatocyte'].copy()),
 ('BRCA',C/'cancer_validation/wu2021_brca.h5ad',None),
 ('Glioblastoma',CACHE/'Glioblastoma_malignant.h5ad',None),('OvarianSerous',CACHE/'OvarianSerous_malignant.h5ad',None),
 ('RCC',CACHE/'RCC_clearcell_malignant.h5ad',None),('Melanoma',CACHE/'Melanoma_malignant.h5ad',None),
 ('Prostate',CACHE/'Prostate_malignant.h5ad',None),('HNSCC',V63/'HNSCC_malignant.h5ad',None),
 ('Neuroblastoma',V64/'Neuroblastoma_malignant.h5ad',None),('B-ALL',V64/'B-ALL_malignant.h5ad',None),
 ('DLBCL',V64/'DLBCL_malignant.h5ad',None),('SqLung',V64/'SqLung_malignant.h5ad',None),
 ('TNBC',V64/'TNBC_malignant.h5ad',None),('PDAC',V65/'PDAC_Peng2019_malignant.h5ad',None),
 ('NSCLC',R/'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',
  lambda a:a[a.obs['cell_type'].astype(str)=='malignant cell'].copy()),
 ('AML',OUT/'vanGalen_AML_malignant.h5ad',None),('MM',OUT/'Ledergor_MM_malignant.h5ad',None)]
CAP=40000; rng=np.random.default_rng(0); CH=20000

def nv(a):
    """Exactly the normalisation scripts/E1_megakaryocyte_complexity.py uses, including the
    ENSEMBL->feature_name remap. A plain uppercase (as first written here) leaves the Census-derived
    cohorts on ENSEMBL ids, where NOTHING matches the catalog and the whole run is silently garbage."""
    import re as _re
    v0=str(a.var.index[0])
    if 'feature_name' in a.var.columns and (v0.startswith('ENSG') or not any(c.isalpha() for c in v0)):
        a.var_names=a.var['feature_name'].astype(str).values
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    base=pd.Series([_re.sub(r'-\d+$','',str(v)) for v in a.var_names]); k=~base.duplicated().values
    a=a[:,k].copy(); a.var_names=base[k].values; return a

cov_rows=[]
def score_chunked(a, tag):
    parts=[]
    for i in range(0,a.n_obs,CH):
        s=pct.score_modules(a[i:i+CH].copy(),method='mean_raw')
        if isinstance(s,dict):
            if i==0:
                cv=s['coverage'].copy().reset_index(); cv.insert(0,'dataset',tag); cov_rows.append(cv)
            s=s['scores']
        parts.append(s[mods].values.astype(np.float32)); gc.collect()
    return np.vstack(parts)

rows=[]; BADCOV=[]
for name,p,pf in DATA:
    try:
        a=sc.read_h5ad(p)
        if pf is not None: a=pf(a)
        a=nv(a); ng=np.asarray((a.X>0).sum(1)).flatten(); a=a[ng>=1000].copy()
        if a.n_obs>CAP:
            a=a[rng.choice(a.n_obs,CAP,replace=False)].copy()
        ng=np.asarray((a.X>0).sum(1)).flatten()
        if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
        S=score_chunked(a,name); z=(S-mu)/sd; pc=z@load
        pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        top=idx[np.argmax(pcn@refn.T,axis=1)]           # over all 154, exactly as Fig. 3D
        for g,t in zip(ng,top): rows.append((name,int(g),t))
        if cov_rows and cov_rows[-1]['n_genes_found'].median()<2: BADCOV.append(name)
        print(f'   {name}: {a.n_obs} cells',flush=True)
        del a,S,z,pc,pcn; gc.collect()
    except Exception as ex:
        print(f'   !! {name}: {ex}',flush=True)
if BADCOV: raise SystemExit('catalog did not match the gene names in: '+', '.join(BADCOV))
PC=pd.DataFrame(rows,columns=['cancer','ng','target'])
PC.to_csv(OUT/'E43_percell_target.csv.gz',index=False,compression='gzip')
print(f'\npooled {len(PC):,} cells, {PC.cancer.nunique()} cohorts',flush=True)

if cov_rows:
    cov=pd.concat(cov_rows,ignore_index=True)
    cov['missing']=cov['n_genes_total']-cov['n_genes_found']
    cov.to_csv(OUT/'E43_gene_coverage.csv',index=False)
    bad=cov[cov.missing>0].groupby('module' if 'module' in cov.columns else cov.columns[1]).agg(datasets_incomplete=('dataset','nunique'),
                                                 worst_found=('n_genes_found','min'),total=('n_genes_total','max'))
    print('\n=== gene recovery: modules missing at least one gene in at least one dataset ===')
    print(bad.sort_values('datasets_incomplete',ascending=False).to_string())

TRACK=['cytotrophoblasts','megakaryocytes','gastric progenitor cells','late spermatids','gastric chief cells',
       'esophageal apical cells','pancreatic acinar cells','hofbauer cells','migrating cytotrophoblasts',
       'differentiating spermatogonia','epididymal clear cells','early primary spermatocytes',
       'undifferentiated spermatogonia','erythrocyte progenitors','erythrocytes','platelets']
PC['dec']=pd.qcut(PC.ng,10,labels=False,duplicates='drop')
tab=PC.groupby('dec').apply(lambda d: pd.Series({**{'ng':d.ng.median(),'n':len(d)},
     **{t:100*np.mean(d.target.values==t) for t in TRACK}}))
tab.to_csv(OUT/'E43_decile_by_target.csv')
ver=[]
for t in TRACK:
    y=tab[t].values; rho,pv=spearmanr(tab.index.values,y)
    ver.append(dict(target=t,first_decile=round(y[0],2),last_decile=round(y[-1],2),
                    rho_with_complexity=round(rho,3),p=float(pv),
                    verdict=('rises (complexity-robust)' if rho>0 else 'FALLS with complexity (erythrocyte-like)')))
V=pd.DataFrame(ver); V.to_csv(OUT/'E43_target_verdicts.csv',index=False)
print('\n=== complexity-decile control, every tracked target ===')
print(V.to_string(index=False))

# boundaries restricted to the complexity-robust members, same ranking as E41
rank=pd.read_csv(OUT/'tableS2_ranking154.csv')
rank=rank[rank.qc_status=='candidate'].sort_values(['recurrence_top5','median_depth_pct'],ascending=False)
falls=set(V.loc[V.rho_with_complexity<0,'target'])
NC=int((pd.read_csv(OUT/'tableS2_ranking154.csv').qc_status=='candidate').sum())
out=[]
for nm,k in [('core (top-2)',2),('neighbourhood (top-4)',4),('broader region (top-14)',14)]:
    mem=rank.cell.head(k).tolist()
    clean=[c for c in mem if c not in falls]
    occ_all=PC.target.isin(mem).mean()*100; occ_cl=PC.target.isin(clean).mean()*100
    out.append(dict(boundary=nm,k=k,members=len(mem),members_clean=len(clean),
                    occupancy_all=round(occ_all,1),chance_all=round(100*k/NC,2),fold_all=round(occ_all/(100*k/NC),1),
                    occupancy_clean=round(occ_cl,1),chance_clean=round(100*len(clean)/NC,2),
                    fold_clean=round(occ_cl/(100*len(clean)/NC),1) if clean else np.nan,
                    dropped='; '.join(c for c in mem if c in falls)))
B=pd.DataFrame(out); B.to_csv(OUT/'E43_boundaries_clean.csv',index=False)
print('\n=== boundary enrichment before and after dropping the erythrocyte-like members ===')
print(B.to_string(index=False))
print('\nDONE')
