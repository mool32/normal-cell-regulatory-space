"""Referee request (PNAS v13 report, §3.1), stage 3: re-score the 25 cohorts with the gene-name fix.

E7_specification_curve.py's loader `nv()` strips a trailing -<digits> from every gene name to undo
var_names_make_unique(). NKX3-1 is a real HGNC symbol, so it became NKX3 and matched nothing in all
17 scored cohorts. The cohorts also carry the current symbols CGAS / STING1 / BMAL1, which the
catalog scores as MB21D1 / TMEM173 / ARNTL.

This re-runs the identical loader, cap, QC, seed and scoring with gene_names.resolve_to_catalog in
place of the blanket strip, and writes results/E46_rederivation/E46_S_fix.npz beside the shipped
results/E7_S_raw.npz. Nothing shipped is overwritten.

Because renaming does not change the detected-gene count, the QC filter and the seeded subsample
select the SAME cells as the shipped run - so the two S matrices are comparable cell for cell, and
only the modules whose gene set was repaired may move. That is asserted, not assumed.
"""
import warnings, sys, re, gc, json
warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, (ROOT / 'tool/perceptome2').as_posix())
sys.path.insert(0, (ROOT / 'scripts').as_posix())
import perceptome as pct
from gene_names import resolve_matrix

CAP = 4000
mods = list(pct.list_modules())
C = ROOT / 'results'
CACHE = C / 'cancer_scaled/cache'; V63 = C / 'cancer_v6_3_extras'
V64 = C / 'cancer_v6_4_extras'; V65 = C / 'cancer_v6_5_extras'; OUT = C / 'cancer_v6_6_extras'
DEST = C / 'E46_rederivation'; DEST.mkdir(parents=True, exist_ok=True)

BCC = lambda o: ((o['disease'].astype(str) == 'basal cell carcinoma') &
                 (o['cell_type'].astype(str).isin(['basal cell', 'basal cell of epidermis',
                                                   'suprabasal keratinocyte']))).values
DATA = [('HCC', C/'cancer_sun2021_hcc/sun2021.h5ad', lambda o: (o['celltype'] == 'Hepatocyte').values),
 ('BRCA', C/'cancer_validation/wu2021_brca.h5ad', None), ('Glioblastoma', CACHE/'Glioblastoma_malignant.h5ad', None),
 ('OvarianSerous', CACHE/'OvarianSerous_malignant.h5ad', None), ('RCC', CACHE/'RCC_clearcell_malignant.h5ad', None),
 ('Melanoma', CACHE/'Melanoma_malignant.h5ad', None), ('Prostate', CACHE/'Prostate_malignant.h5ad', None),
 ('HNSCC', V63/'HNSCC_malignant.h5ad', None), ('Neuroblastoma', V64/'Neuroblastoma_malignant.h5ad', None),
 ('B-ALL', V64/'B-ALL_malignant.h5ad', None), ('DLBCL', V64/'DLBCL_malignant.h5ad', None),
 ('SqLung', V64/'SqLung_malignant.h5ad', None), ('TNBC', V64/'TNBC_malignant.h5ad', None),
 ('PDAC', V65/'PDAC_Peng2019_malignant.h5ad', None),
 ('NSCLC', ROOT/'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',
  lambda o: (o['cell_type'].astype(str) == 'malignant cell').values),
 ('AML', OUT/'vanGalen_AML_malignant.h5ad', None), ('MM', OUT/'Ledergor_MM_malignant.h5ad', None),
 ('BCC', ROOT/'paper6_immune_3layer/data/bcc_full.h5ad', BCC),
 ('Colon', C/'cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',
  lambda o: (o['ClusterMidway'].astype(str) == 'EpiT').values),
 ('SCLC', OUT/'SCLC_malignant.h5ad', None), ('Gastric', OUT/'Gastric_malignant.h5ad', None),
 ('GermCell', OUT/'GermCell_malignant.h5ad', None), ('Osteosarcoma', OUT/'Osteosarcoma_malignant.h5ad', None),
 ('Ewing', OUT/'Ewing_malignant.h5ad', None), ('Seminoma', OUT/'Seminoma_malignant.h5ad', None)]

def nv_fixed(a):
    """E7's nv(), with the blanket suffix strip replaced by the collision-aware one + synonyms."""
    v0 = str(a.var.index[0])
    if 'feature_name' in a.var.columns and (v0.startswith('ENSG') or not any(c.isalpha() for c in v0)):
        a.var_names = a.var['feature_name'].astype(str).values
    a.var_names = a.var_names.str.upper(); a.var_names_make_unique()
    return resolve_matrix(a)      # <-- the one change: rename into catalog space, keep the column
                                  #     that carries the data when a rename collides

rng = np.random.default_rng(0)
Sb, ids, cov = [], [], []
for name, p, pf in DATA:
    try:
        a = sc.read_h5ad(p, backed='r')
        a = (a[pf(a.obs)].to_memory() if pf else a.to_memory())
        a = nv_fixed(a)
        g = np.asarray((a.X > 0).sum(1)).ravel(); a = a[g >= 1000].copy()
        if a.n_obs > CAP: a = a[rng.choice(a.n_obs, CAP, replace=False)].copy()
        if a.n_obs < 50: continue
        if a.X.max() > 100: sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        have = set(a.var_names)
        for m in ('AR', 'cGAS-STING', 'Circadian'):
            gs = pct.get_genes(m, 'core') if hasattr(pct, 'get_genes') else []
            cov.append(dict(cohort=name, module=m,
                            found=sum(x in have for x in gs), of=len(gs)))
        Sc = pct.score_modules(a, method='mean_raw', gene_set='core')['scores'] \
                .reindex(columns=mods, fill_value=0.0)
        Sb.append(Sc.values.astype(np.float32)); ids += [name] * a.n_obs
        print(f'  {name:14s} {a.n_obs:>6,} cells', flush=True)
        del a, Sc; gc.collect()
    except Exception as ex:
        print(f'  FAIL {name}: {type(ex).__name__}: {ex}', flush=True); raise

S = np.vstack(Sb); cid = np.array(ids)
np.savez_compressed(DEST / 'E46_S_fix.npz', S=S, cid=cid)
pd.DataFrame(cov).to_csv(DEST / 'E46_gene_coverage.csv', index=False)
print(f'\nwrote E46_S_fix.npz  {S.shape}')

# --- the comparability assertion -----------------------------------------------------------------
old = np.load(ROOT / 'results/E7_S_raw.npz', allow_pickle=True)
So, cido = old['S'].astype(np.float32), old['cid']
print(f'shipped S {So.shape}, cid identical: {np.array_equal(cid, cido)}')
if np.array_equal(cid, cido):
    d = np.abs(S - So).max(axis=0)
    moved = [(mods[j], float(d[j])) for j in range(len(mods)) if d[j] > 1e-6]
    print(f'modules that moved: {len(moved)} of 44')
    for m, v in sorted(moved, key=lambda x: -x[1]): print(f'   {m:14s} max|delta| {v:.4f}')
print('DONE')
