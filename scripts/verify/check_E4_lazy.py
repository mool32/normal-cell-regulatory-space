#!/usr/bin/env python3
"""Assert E4's lazy cohort loader returns exactly what the eager one does.

E4's gene-level control is one of the two bounds on the curation prior, and its stored run covered
6 of the 17 listed cohorts because Ledergor_MM (dense, 7.4 GB) could not be held in memory twice.
The loader was split into an eager and a backed path so the run can complete; the two must agree
cell for cell and value for value, or the control is a different analysis than the one described.

Run on every dense cohort small enough to take both paths, forcing the lazy one on.
"""
import sys, os, re, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[2]
src = (ROOT / 'scripts' / 'E4_genelevel_signature.py').read_text()
ns = {'np': np, 'pd': pd, 'sc': sc, 'ad': ad, 're': re, 'Path': Path}
exec(compile(src[src.index('def nv_plan'):src.index('percell=[]')], 'E4', 'exec'), ns)

C = ROOT / 'results'
HEP = lambda o: (o['celltype'] == 'Hepatocyte').values
COHORTS = [('HCC', C / 'cancer_sun2021_hcc/sun2021.h5ad', HEP),
           ('Glioblastoma', C / 'cancer_scaled/cache/Glioblastoma_malignant.h5ad', None),
           ('OvarianSerous', C / 'cancer_scaled/cache/OvarianSerous_malignant.h5ad', None),
           ('Melanoma', C / 'cancer_scaled/cache/Melanoma_malignant.h5ad', None),
           ('AML', C / 'cancer_v6_6_extras/vanGalen_AML_malignant.h5ad', None),
           ('DLBCL', C / 'cancer_v6_4_extras/DLBCL_malignant.h5ad', None),
           ('TNBC', C / 'cancer_v6_4_extras/TNBC_malignant.h5ad', None),
           # a positional mask on a cohort that fits both ways, so the obs-mask branch itself is
           # tested: HCC and NSCLC are the only cohorts that carry a filter and HCC is too large here
           ('AML+mask', C / 'cancer_v6_6_extras/vanGalen_AML_malignant.h5ad',
            lambda o: (np.arange(len(o)) % 2 == 0))]


def arr(a):
    return np.asarray(a.todense() if hasattr(a, 'todense') else a)


def main():
    bad = 0
    for nm, p, pf in COHORTS:
        if not p.exists():
            print(f'  skip {nm:14s} absent'); continue
        gb = ns['dense_gb'](p)
        if not gb:
            print(f'  skip {nm:14s} sparse on disk; the lazy path is dense-only'); continue
        if gb > 2.0:
            # the eager arm of the comparison would have to hold it, which is the thing the lazy
            # path exists to avoid; the code under test is the same for every dense cohort
            print(f'  skip {nm:14s} {gb:.1f} GB dense, too large to run both arms here'); continue
        ns['BACKED_GB'] = 1e9; e, _ = ns['read_cohort'](p, pf)
        if e.n_obs > 1500:
            e = e[np.random.default_rng(0).choice(e.n_obs, 1500, replace=False)].copy()
        ns['BACKED_GB'] = 0.0; l, _ = ns['read_cohort'](p, pf)
        ok = (e.shape == l.shape and list(e.var_names) == list(l.var_names)
              and list(e.obs_names) == list(l.obs_names) and np.array_equal(arr(e.X), arr(l.X)))
        bad += not ok
        print(f'  {"ok  " if ok else "FAIL"} {nm:14s} {e.shape} eager vs {l.shape} lazy')
        del e, l
    print(f'\n{len(COHORTS)} cohorts, {bad} disagreement(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
