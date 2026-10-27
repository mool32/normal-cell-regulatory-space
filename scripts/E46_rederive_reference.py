"""Referee request (PNAS v13 report, §3.1), stage 1-2: re-derive the HPA reference and the eigenspace
with the gene-name fix of scripts/gene_names.py active.

Note S-Q(ii) found three genes systematically lost and says the corrected catalog "has not been
re-derived end to end". The HPA source carries the CURRENT symbols (CGAS, STING1, BMAL1: 154 rows
each) while the catalog scores the DEPRECATED ones (MB21D1, TMEM173, ARNTL: 0 rows), so the
reference has been built without both cGAS-STING sensors and without one of five Circadian genes.

Nothing shipped is overwritten: the corrected artifacts go to results/E46_rederivation/.

The method here is a re-implementation of tool/perceptome2/scripts/02 and 12, so it is only
trustworthy if it reproduces them. Stage 0 therefore rebuilds BOTH references with the fix OFF and
asserts agreement with the shipped npz before the fixed build is allowed to run.
"""
import csv, json, sys, math
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, (ROOT / 'tool/perceptome2').as_posix())
sys.path.insert(0, (ROOT / 'scripts').as_posix())
from perceptome.catalog import load_catalog, get_genes
from gene_names import SYNONYMS

HPA_TSV = ROOT / 'paper2/data/rna_single_cell_type.tsv'
HPA_CTCLASS = ROOT / 'paper2/data/rna_single_cell_type_cell_types.tsv'
SHIPPED = ROOT / 'tool/perceptome2/perceptome/perceptivity/data'
OUT = ROOT / 'results/E46_rederivation'; OUT.mkdir(parents=True, exist_ok=True)

cat = load_catalog(); assert cat['n_modules'] == 44
MODULES = sorted(cat['modules'].keys())
CORE = {m: set(get_genes(m, 'core', cat)) for m in MODULES}
ACT = {m: set(get_genes(m, 'activity', cat)) for m in MODULES}

def read_hpa(resolve):
    """(cell_type, gene) -> nCPM. `resolve` maps the data's symbol into the catalog's symbol space."""
    rows = []
    with open(HPA_TSV) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            g = r['Gene name'].upper()
            if resolve: g = SYNONYMS.get(g, g)
            rows.append((r['Cell type'], g, float(r['nCPM'])))
    return rows

def build(resolve):
    """Both calibrations, exactly as the shipped builders compute them.

    The two differ in how a (gene, cell type) pair absent from the HPA file is treated, and the
    difference is load-bearing: scripts/02 (pseudobulk) SKIPS the pair, so the module mean is over
    the genes actually measured in that cell type, while scripts/12 (sc-scaled) builds a dense
    matrix and lets the absent pair enter as a zero. Reproducing the shipped artifacts means
    reproducing both conventions, not reconciling them - this re-implementation got it wrong in the
    first pass and stage 0 caught it.
    """
    import scanpy as sc, anndata as ad
    rows = read_hpa(resolve)
    # --- pseudobulk, scripts/02 semantics: dict of dicts, missing pairs skipped -------------------
    expr = {}
    for c, g, v in rows: expr.setdefault(g, {})[c] = v
    cts = sorted({c for d in expr.values() for c in d})
    def mat_pb(sets):
        M = np.full((len(cts), len(MODULES)), np.nan, dtype=np.float32)
        for j, m in enumerate(MODULES):
            gs = [g for g in sets[m] if g in expr]
            for i, ct in enumerate(cts):
                v = [math.log1p(expr[g][ct]) for g in gs if ct in expr[g]]
                if v: M[i, j] = float(np.mean(v))
        return M
    Rp, Ap = mat_pb(CORE), mat_pb(ACT)
    # --- sc-scaled, scripts/12 semantics: dense nCPM matrix, normalize_total + log1p --------------
    genes = sorted({g for _, g, _ in rows})
    ci = {c: i for i, c in enumerate(cts)}; gi = {g: i for i, g in enumerate(genes)}
    X = np.zeros((len(cts), len(genes)), dtype=np.float32)
    for c, g, v in rows: X[ci[c], gi[g]] = v          # a collision keeps the last row, as in 12
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=cts), var=pd.DataFrame(index=genes))
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    Xs = a.X if isinstance(a.X, np.ndarray) else a.X.toarray()
    def mat_sc(sets):
        M = np.zeros((len(cts), len(MODULES)), dtype=np.float32)
        for j, m in enumerate(MODULES):
            idx = [gi[g] for g in sets[m] if g in gi]
            if idx: M[:, j] = Xs[:, idx].mean(axis=1)
        return M
    Rs, As = mat_sc(CORE), mat_sc(ACT)
    return dict(cell_types=cts, modules=MODULES, R=Rp, A=Ap, Rs=Rs, As=As,
                n_genes_core={m: sum(g in expr for g in CORE[m]) for m in MODULES})

def coverage(b, tag):
    print(f'  [{tag}] core-gene coverage of the three repaired modules:')
    for m in ('AR', 'cGAS-STING', 'Circadian'):
        print(f'      {m:12s} {b["n_genes_core"][m]} of {len(CORE[m])} core genes found in HPA')

# ---------------- stage 0: equivalence against the shipped artifacts -----------------------------
print('stage 0  rebuilding with the fix OFF, to test this re-implementation')
off = build(resolve=False)
coverage(off, 'fix off')
ok = True
for tag, key, path in [('pseudobulk', 'R', SHIPPED / 'hpa_perceptivity_v03.npz'),
                       ('sc-scaled', 'Rs', SHIPPED / 'hpa_perceptivity_v03_scscaled.npz')]:
    ship = np.load(path)['R']
    mine = off[key]
    if ship.shape != mine.shape:
        print(f'  FAIL {tag}: shape {mine.shape} vs shipped {ship.shape}'); ok = False; continue
    both = ~(np.isnan(ship) | np.isnan(mine))
    d = np.abs(ship[both] - mine[both])
    nan_same = (np.isnan(ship) == np.isnan(mine)).all()
    print(f'  {tag:11s} max|diff| = {d.max():.3e}   median {np.median(d):.3e}   '
          f'NaN pattern identical: {nan_same}')
    if not (d.max() < 1e-4 and nan_same): ok = False
assert ok, 'the re-implementation does not reproduce the shipped reference; the fixed build is not trustworthy'
print('  EQUIVALENT - the fixed build below is a one-change delta on a reproduced baseline\n')

# ---------------- stage 1: the fixed reference ---------------------------------------------------
print('stage 1  rebuilding with the fix ON')
on = build(resolve=True)
coverage(on, 'fix on')
np.savez_compressed(OUT / 'hpa_perceptivity_FIX.npz', R=on['R'], A=on['A'],
                    C=on['R'] - on['A'], headroom=np.nanmax(on['A'], 0)[None, :] - on['A'])
np.savez_compressed(OUT / 'hpa_perceptivity_FIX_scscaled.npz', R=on['Rs'], A=on['As'],
                    C=on['Rs'] - on['As'], headroom=on['As'].max(0)[None, :] - on['As'])
ctc = {}
with open(HPA_CTCLASS) as f:
    for r in csv.DictReader(f, delimiter='\t'): ctc[r['Cell type']] = r.get('Cell type class', '')
meta = dict(cell_types=on['cell_types'], modules=MODULES,
            cell_type_class=[ctc.get(c, '') for c in on['cell_types']],
            A_max=np.nanmax(on['A'], 0).tolist(), version='0.3.genefix')
json.dump(meta, open(OUT / 'hpa_perceptivity_FIX.json', 'w'), indent=1)
meta_s = dict(meta, A_max=on['As'].max(0).tolist(), mode='single_cell_scaled')
json.dump(meta_s, open(OUT / 'hpa_perceptivity_FIX_scscaled.json', 'w'), indent=1)

# how much did the three modules move?
print('\n  movement in the three repaired modules (sc-scaled R, across 154 cell types):')
for m in ('AR', 'cGAS-STING', 'Circadian'):
    j = MODULES.index(m)
    a_, b_ = off['Rs'][:, j], on['Rs'][:, j]
    print(f'      {m:12s} mean {a_.mean():.4f} -> {b_.mean():.4f}   '
          f'sd {a_.std():.4f} -> {b_.std():.4f}   max|delta| {np.abs(b_-a_).max():.4f}   '
          f'rank-corr {pd.Series(a_).corr(pd.Series(b_), method="spearman"):.3f}')
unchanged = [m for m in MODULES if np.abs(on['Rs'][:, MODULES.index(m)] - off['Rs'][:, MODULES.index(m)]).max() < 1e-6]
print(f'  modules unchanged by the fix: {len(unchanged)} of 44')

# ---------------- stage 2: the eigenspace from the fixed reference -------------------------------
print('\nstage 2  rebuilding the eigenspace (same rebuild() as scripts/03: Kaiser, 100 bootstraps, seed 42)')
from perceptome.eigenspace.rebuild import rebuild
def eig(Rmat, tag):
    Rdf = pd.DataFrame(Rmat, index=on['cell_types'], columns=MODULES).fillna(0.0)
    out = rebuild(Rdf, kaiser=True, bootstrap_n=100, seed=42)
    print(f'  [{tag}] n_kaiser={out["n_kaiser"]}  var PC1-3 = '
          f'{[f"{v:.3f}" for v in out["var_explained"][:3]]}  cum9={sum(out["var_explained"][:9]):.3f}')
    return out
e_off, e_on = eig(off['R'], 'fix off'), eig(on['R'], 'fix on')
ship = json.load(open(ROOT / 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json'))
print(f'  [shipped ] n_kaiser={ship.get("n_kaiser")}  var PC1-3 = '
      f'{[f"{v:.3f}" for v in ship["var_explained"][:3]]}')
# sign-insensitive agreement of the loading matrices
def L(o): return np.array([[o['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in MODULES])
for tag, o in (('fix off vs shipped', e_off), ('fix on vs shipped', e_on)):
    a_, b_ = L(o), L(ship)
    cos = np.abs((a_ * b_).sum(0) / (np.linalg.norm(a_, axis=0) * np.linalg.norm(b_, axis=0)))
    print(f'  {tag:20s} |cos| per PC: {" ".join(f"{c:.3f}" for c in cos)}')
e_on['module_order'] = MODULES
json.dump(e_on, open(OUT / 'reference_FIX.json', 'w'), indent=1, default=float)
print(f'\nwrote {OUT.relative_to(ROOT)}/  (reference_FIX.json + two FIX npz)')
print('DONE')
