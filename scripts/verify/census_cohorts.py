#!/usr/bin/env python3
"""Resolve every Census-derived panel cohort and print what Table S1 does not say.

Table S1 gives the twelve Census cohorts no identifier, and `E18_purity_audit.csv` records six of
them as "pre-filtered ... author-annotated (high)". Reading the stored objects directly contradicts
that for four, and shows a second thing nothing in the paper states: `is_primary_data` was never
filtered, so a cell deposited in more than one Census dataset is counted more than once.

Three columns are the point:
  malignant label    the cell_type values the cohort was selected on. "malignant cell" and
                     "neoplastic cell" are author malignant annotations; anything else is a
                     heuristic -- and a normal epithelial label is the same heuristic that produced
                     the retracted lung control (Fig. S11-II).
  unique cells       distinct `observation_joinid`. Rows minus this is double counting.
  primary rows       rows with `is_primary_data == True`, the filter a Census query should carry.

Writes results/census_cohorts.csv. Reads obs only, so it is seconds, not minutes.
"""
import sys, os
from pathlib import Path
import anndata as ad, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
C = ROOT / 'results'
CACHE = C / 'cancer_scaled/cache'
V63, V64, OUT = C / 'cancer_v6_3_extras', C / 'cancer_v6_4_extras', C / 'cancer_v6_6_extras'

# every panel cohort whose cells come from the Census (the other 13 come from GEO)
COHORTS = [
    ('Glioblastoma', CACHE / 'Glioblastoma_malignant.h5ad', None),
    ('Ovarian serous', CACHE / 'OvarianSerous_malignant.h5ad', None),
    ('RCC', CACHE / 'RCC_clearcell_malignant.h5ad', None),
    ('Melanoma', CACHE / 'Melanoma_malignant.h5ad', None),
    ('Prostate', CACHE / 'Prostate_malignant.h5ad', None),
    ('HNSCC', V63 / 'HNSCC_malignant.h5ad', None),
    ('Neuroblastoma', V64 / 'Neuroblastoma_malignant.h5ad', None),
    ('B-ALL', V64 / 'B-ALL_malignant.h5ad', None),
    ('DLBCL', V64 / 'DLBCL_malignant.h5ad', None),
    ('Squamous lung', V64 / 'SqLung_malignant.h5ad', None),
    ('TNBC', V64 / 'TNBC_malignant.h5ad', None),
    # the only Census cohort whose malignant selection happens at load time rather than in the
    # stored object (scripts/E2e_module_robustness_fullfid.py), so the same filter is applied here
    ('NSCLC', ROOT / 'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',
     lambda o: (o['cell_type'].astype(str) == 'malignant cell').values),
    ('SCLC', OUT / 'SCLC_malignant.h5ad', None),
    ('Gastric', OUT / 'Gastric_malignant.h5ad', None),
]
AUTHOR = {'malignant cell', 'neoplastic cell'}


def main():
    rows = []
    for name, p, pf in COHORTS:
        if not p.exists():
            print(f'  MISSING {name}: {p}'); continue
        a = ad.read_h5ad(p, backed='r')
        o = a.obs
        if pf is not None: o = o[pf(o)]
        n_rows = len(o)
        lab = sorted(str(x) for x in o['cell_type'].unique()) if 'cell_type' in o else []
        lab = [x for x in lab if (o['cell_type'].astype(str) == x).any()]
        uniq = int(o['observation_joinid'].nunique()) if 'observation_joinid' in o else None
        prim = int((o['is_primary_data'] == True).sum()) if 'is_primary_data' in o else None
        ds = [d for d in (o['dataset_id'].astype(str).unique() if 'dataset_id' in o else [])
              if (o['dataset_id'].astype(str) == d).any()]
        rows.append(dict(cohort=name, rows=n_rows, unique_cells=uniq,
                         duplicated_rows=(n_rows - uniq) if uniq is not None else None,
                         primary_rows=prim,
                         donors=int(o['donor_id'].nunique()) if 'donor_id' in o else None,
                         datasets=len(ds),
                         assays='; '.join(sorted({str(x) for x in o['assay']})) if 'assay' in o else '',
                         malignant_label='; '.join(lab),
                         selection=('author malignant annotation'
                                    if lab and set(lab) <= AUTHOR else 'heuristic cell-type label'),
                         dataset_ids='; '.join(d[:8] for d in ds)))
        a.file.close()
    D = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(D[['cohort', 'rows', 'unique_cells', 'duplicated_rows', 'primary_rows', 'donors',
             'datasets', 'selection']].to_string(index=False))
    print()
    for r in rows:
        if r['selection'] != 'author malignant annotation':
            print(f"  HEURISTIC  {r['cohort']:16s} selected on: {r['malignant_label']}")
    print()
    for r in rows:
        if r['duplicated_rows']:
            print(f"  DOUBLE-COUNTED {r['cohort']:16s} {r['duplicated_rows']:6,d} of "
                  f"{r['rows']:,d} rows ({100 * r['duplicated_rows'] / r['rows']:.1f}%) are a cell "
                  f"already in the cohort")
    D.to_csv(C / 'census_cohorts.csv', index=False)
    print(f'\nwrote results/census_cohorts.csv  ({len(D)} Census cohorts of the 25-cancer panel)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
