#!/usr/bin/env python3
"""Does the unfiltered `is_primary_data` change the panel result? Recompute the baseline both ways.

`scripts/verify/census_cohorts.py` shows six Census cohorts count some cells more than once, because
a Census query without `is_primary_data == True` returns a cell once per dataset it was deposited in.
Prostate is the extreme: 47% of its 30,000 rows are a cell already in the cohort, the atlas plus its
Epithelial and Cycling subsets, so cycling prostate cells carry triple weight in that centroid.

This runs the canonical pipeline of `scripts/E2e_module_robustness_fullfid.py` twice over the same
25 cohorts -- once as the paper does, once dropping duplicate `observation_joinid` within each cohort
before QC -- and prints both baselines. Everything else is identical, including the seed, so the only
difference is the duplicate rows.
"""
import sys, re, gc, json, warnings
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tool/perceptome2'))
import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity

CORE = ['megakaryocytes', 'cytotrophoblasts', 'migrating cytotrophoblasts']
NPC, CAP = 9, 6000
mods = list(pct.list_modules())
Rm = load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
refcells = list(Rm.index)
core_idx = np.array([refcells.index(c) for c in CORE])
EXCLUDE = np.array([refcells.index(c) for c in ['erythrocytes', 'platelets'] if c in refcells])
RM = Rm.values.astype(float)

C = ROOT / 'results'; CACHE = C / 'cancer_scaled/cache'
V63, V64, V65, OUT = (C / 'cancer_v6_3_extras', C / 'cancer_v6_4_extras',
                      C / 'cancer_v6_5_extras', C / 'cancer_v6_6_extras')
BCC = lambda o: ((o['disease'].astype(str) == 'basal cell carcinoma') & (o['cell_type'].astype(str).isin(
    ['basal cell', 'basal cell of epidermis', 'suprabasal keratinocyte']))).values
DATA = [('HCC', C / 'cancer_sun2021_hcc/sun2021.h5ad', lambda o: (o['celltype'] == 'Hepatocyte').values),
 ('BRCA', C / 'cancer_validation/wu2021_brca.h5ad', None), ('Glioblastoma', CACHE / 'Glioblastoma_malignant.h5ad', None),
 ('OvarianSerous', CACHE / 'OvarianSerous_malignant.h5ad', None), ('RCC', CACHE / 'RCC_clearcell_malignant.h5ad', None),
 ('Melanoma', CACHE / 'Melanoma_malignant.h5ad', None), ('Prostate', CACHE / 'Prostate_malignant.h5ad', None),
 ('HNSCC', V63 / 'HNSCC_malignant.h5ad', None), ('Neuroblastoma', V64 / 'Neuroblastoma_malignant.h5ad', None),
 ('B-ALL', V64 / 'B-ALL_malignant.h5ad', None), ('DLBCL', V64 / 'DLBCL_malignant.h5ad', None),
 ('SqLung', V64 / 'SqLung_malignant.h5ad', None), ('TNBC', V64 / 'TNBC_malignant.h5ad', None),
 ('PDAC', V65 / 'PDAC_Peng2019_malignant.h5ad', None),
 ('NSCLC', ROOT / 'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',
  lambda o: (o['cell_type'].astype(str) == 'malignant cell').values),
 ('AML', OUT / 'vanGalen_AML_malignant.h5ad', None), ('MM', OUT / 'Ledergor_MM_malignant.h5ad', None),
 ('BCC', ROOT / 'paper6_immune_3layer/data/bcc_full.h5ad', BCC),
 ('Colon', C / 'cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',
  lambda o: (o['ClusterMidway'].astype(str) == 'EpiT').values),
 ('SCLC', OUT / 'SCLC_malignant.h5ad', None), ('Gastric', OUT / 'Gastric_malignant.h5ad', None),
 ('GermCell', OUT / 'GermCell_malignant.h5ad', None), ('Osteosarcoma', OUT / 'Osteosarcoma_malignant.h5ad', None),
 ('Ewing', OUT / 'Ewing_malignant.h5ad', None), ('Seminoma', OUT / 'Seminoma_malignant.h5ad', None)]


def nv(a):
    v0 = str(a.var.index[0])
    if 'feature_name' in a.var.columns and (v0.startswith('ENSG') or not any(c.isalpha() for c in v0)):
        a.var_names = a.var['feature_name'].astype(str).values
    a.var_names = a.var_names.str.upper(); a.var_names_make_unique()
    base = pd.Series([re.sub(r'-\d+$', '', str(v)) for v in a.var_names]); k = ~base.duplicated().values
    a = a[:, k].copy(); a.var_names = base[k].values; return a


def eigenspace(cols):
    X = RM[:, cols]; mean = X.mean(0); std = X.std(0); std[std < 1e-10] = 1.0
    Z = (X - mean) / std; ev, vec = np.linalg.eigh(np.corrcoef(Z.T))
    idx = np.argsort(ev)[::-1]; L = vec[:, idx[:NPC]]
    refP = Z @ L
    return mean, std, L, refP / (np.linalg.norm(refP, axis=1, keepdims=True) + 1e-12), int((ev > 1).sum())


def dense_gb(path):
    """In-memory GB if X is a dense HDF5 dataset, else 0 (same test as E4's loader)."""
    import h5py
    with h5py.File(path, 'r') as f:
        X = f['X']
        return 0.0 if not isinstance(X, h5py.Dataset) else (
            X.shape[0] * X.shape[1] * X.dtype.itemsize / 2**30)


def score(dedup):
    """Score every cohort on the 44 modules, optionally dropping duplicate Census cells first.

    Backed reading is used only for a cohort whose DENSE matrix will not fit twice (Sun HCC 6.9 GB,
    Ledergor MM 7.4 GB). A backed SPARSE matrix is read into memory first and subset there: indexing
    a CSRDataset by an index array raises inside scipy in this anndata, which silently cost five
    cohorts the first time this ran -- and a comparison that loses cohorts on one side is no
    comparison at all, so the two passes are asserted to cover the same cohorts below."""
    rng = np.random.default_rng(0); blocks, ids, dropped = [], [], {}
    for name, p, pf in DATA:
        try:
            lazy = dense_gb(p) >= 1.0
            a = sc.read_h5ad(p, backed='r') if lazy else sc.read_h5ad(p)
            keep = np.where(np.asarray(pf(a.obs)))[0] if pf else np.arange(a.n_obs)
            if dedup and 'observation_joinid' in a.obs.columns:
                j = a.obs['observation_joinid'].values[keep]
                first = ~pd.Series(j).duplicated().values
                if (~first).any(): dropped[name] = int((~first).sum())
                keep = keep[first]
            if len(keep) == a.n_obs:
                a = a.to_memory() if lazy else a
            else:
                a = a[keep].to_memory() if lazy else a[keep].copy()
            a = nv(a)
            ng = np.asarray((a.X > 0).sum(1)).ravel(); a = a[ng >= 1000].copy()
            if a.n_obs > CAP: a = a[rng.choice(a.n_obs, CAP, replace=False)].copy()
            if a.n_obs < 50: print(f'  !!{name}: only {a.n_obs} cells', flush=True); continue
            if a.X.max() > 100: sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
            S = pct.score_modules(a, method='mean_raw', gene_set='core')['scores'].reindex(
                columns=mods, fill_value=0.0)
            blocks.append(S.values.astype(np.float32)); ids += [name] * a.n_obs
            del a, S; gc.collect()
        except Exception as ex:
            print(f'  !!{name}: {ex}', flush=True)
    return np.vstack(blocks), np.array(ids), dropped


def metrics(Sall, cancer_ids):
    cols = list(range(len(mods)))
    mean, std, L, refPn, nk = eigenspace(cols)
    z = (Sall[:, cols] - mean) / std; pc = z @ L
    pcn = pc / (np.linalg.norm(pc, axis=1, keepdims=True) + 1e-12)
    sim = pcn @ refPn.T; sim[:, EXCLUDE] = -1e9; top1 = np.argmax(sim, axis=1)
    uc = sorted(set(cancer_ids))
    recur = np.zeros(154); cents = []
    for cid in uc:
        m = cancer_ids == cid
        recur[np.argsort(-np.bincount(top1[m], minlength=154))[:5]] += 1
        cents.append(pcn[m].mean(0))
    recur[EXCLUDE] = -1; ro = np.argsort(-recur)
    rrank = {refcells[ro[i]]: i + 1 for i in range(154)}
    CC = np.array(cents); CC = CC / (np.linalg.norm(CC, axis=1, keepdims=True) + 1e-12)
    attr = CC.mean(0); attr /= np.linalg.norm(attr) + 1e-12
    cosF = refPn @ attr; cosF[EXCLUDE] = -1e9; ao = np.argsort(-cosF)
    arank = {refcells[ao[i]]: i + 1 for i in range(154)}
    return dict(cells=len(top1), cohorts=len(uc), nk=nk,
                corefrac=round(100 * float(np.isin(top1, core_idx).mean()), 1),
                cyto_rec=rrank['cytotrophoblasts'], cyto_attr=arank['cytotrophoblasts'],
                mega_rec=rrank['megakaryocytes'], attr_top3=[refcells[ao[i]] for i in range(3)],
                rec_top2=[refcells[ro[0]], refcells[ro[1]]])


def main():
    out, cohorts = {}, {}
    for tag, dedup in [('as published', False), ('deduplicated', True)]:
        S, ids, dropped = score(dedup)
        cohorts[tag] = sorted(set(ids))
        out[tag] = metrics(S, ids)
        if dropped:
            print(f'  dedup dropped: ' + ', '.join(f'{k} {v:,}' for k, v in dropped.items()))
        del S; gc.collect()
    print()
    keys = ['cells', 'cohorts', 'nk', 'corefrac', 'cyto_attr', 'cyto_rec', 'mega_rec']
    print(f'{"":14s} ' + ' '.join(f'{k:>10s}' for k in keys))
    for tag in out:
        print(f'{tag:14s} ' + ' '.join(f'{str(out[tag][k]):>10s}' for k in keys))
    for tag in out:
        print(f'  {tag:14s} attractor top-3: {out[tag]["attr_top3"]}')
    assert cohorts['as published'] == cohorts['deduplicated'], (
        'the two passes cover different cohorts, so nothing can be concluded: '
        + repr(set(cohorts['as published']) ^ set(cohorts['deduplicated'])))
    print(f"both passes cover the same {len(cohorts['as published'])} cohorts")
    same = all(out['as published'][k] == out['deduplicated'][k]
               for k in ['cyto_attr', 'cyto_rec', 'mega_rec'])
    d = abs(out['as published']['corefrac'] - out['deduplicated']['corefrac'])
    print(f'\nranks unchanged: {same} | core fraction moves {d:.1f} points')
    out['harness'] = dict(cap_per_cohort=CAP, npc=NPC, note=(
        'E2e\'s harness: the cap is 6,000 cells per cohort, not the 4,000 of the paper\'s own '
        'figures, because this script reproduces E2e\'s baseline so the two passes differ only '
        'in the duplicate rows.'))
    (ROOT / 'results/census_dedup_effect.json').write_text(json.dumps(out, indent=1))
    print('wrote results/census_dedup_effect.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
