"""Referee request (PNAS v13 report, §3.1), item (iii): the colorectal within-patient increase,
raw and shape-only, recomputed on the re-derived catalog.

Same cells, same arms, same three-cell core, same quarantine-before-assignment rule, same per-arm
floor of 20 and the same full-precision Wilcoxon as scripts/E36_colorectal_comparator.py and
scripts/verify/regen_pelka_quarantined.py. The only changes are the gene names (resolve_matrix),
the reference and the eigenspace.

Gate: the arm rebuilt WITHOUT the fix must reproduce E36 row (i) exactly, or the comparison is
measuring this script rather than the fix.
"""
import warnings, sys, gc, json
warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import wilcoxon, binomtest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, (ROOT / 'tool/perceptome2').as_posix())
sys.path.insert(0, (ROOT / 'scripts').as_posix())
import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from gene_names import resolve_matrix
FIX = ROOT / 'results/E46_rederivation'
CH = 12000
demag = lambda X: X - X.mean(1, keepdims=True)
ident = lambda X: X

def frame(fixed):
    """(module order, loadings, REF dict, core indices, quarantine indices) for one catalog."""
    if fixed:
        e = json.load(open(FIX / 'reference_FIX.json'))
        meta = json.load(open(FIX / 'hpa_perceptivity_FIX_scscaled.json'))
        Rm = pd.DataFrame(np.load(FIX / 'hpa_perceptivity_FIX_scscaled.npz')['R'],
                          index=meta['cell_types'], columns=meta['modules'])
    else:
        e = json.load(open(ROOT / 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json'))
        Rm = load_hpa_perceptivity(mode='single_cell_scaled')['R']
    mods = e['module_order']; Rm = Rm[mods]; cells = list(Rm.index)
    load = np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
    def refset(tf):
        R = tf(Rm.values.astype(float)); mu, sd = R.mean(0), R.std(0) + 1e-12
        rp = ((R - mu) / sd) @ load
        return mu, sd, rp / (np.linalg.norm(rp, axis=1, keepdims=True) + 1e-12)
    REF = {'RAW': (refset(ident), ident), 'SHAPE': (refset(demag), demag)}
    CORE3 = [cells.index(c) for c in ('megakaryocytes', 'cytotrophoblasts', 'migrating cytotrophoblasts')]
    EX = [cells.index(c) for c in ('erythrocytes', 'platelets')]
    return mods, load, REF, CORE3, EX

# ---- the cells, and the arm definition E36 used --------------------------------------------------
orig = sc.read_h5ad(ROOT / 'results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',
                    backed='r')
full = orig.obs['ClusterFull'].astype(str).values; del orig
C = ad.read_h5ad(ROOT / 'perceptome_validation/cache/pelka_epi.h5ad', backed='r')
obs = C.obs.copy()
obs['arm2'] = np.where(pd.Series(full[obs['src_row'].astype(int).values]).str.startswith('Tumor').values,
                       'malig', 'normal')
print(f'Pelka epithelium {C.shape}, {obs.patient.nunique()} patients', flush=True)

# ---- pass 1: per-column nonzero counts, so a rename collision keeps the column with the data ----
nz = np.zeros(C.shape[1], dtype=np.int64)
for st in range(0, C.shape[0], CH):
    m = np.zeros(C.shape[0], bool); m[st:st + CH] = True
    a = C[m].to_memory()
    nz += np.asarray((a.X != 0).sum(axis=0)).ravel()
    del a; gc.collect()
print('pass 1 done (nonzero counts for the collision rule)', flush=True)

def keep_mask(names, nz):
    from gene_names import resolve_to_catalog
    s = pd.Series(resolve_to_catalog(list(names)))
    dup = set(s[s.duplicated(keep=False)].tolist())
    keep = np.zeros(len(s), bool); best = {}
    for i, nm in enumerate(s):
        if nm not in dup: keep[i] = True
        elif nm not in best or nz[i] > nz[best[nm]]: best[nm] = i
    for i in best.values(): keep[i] = True
    return keep, s

def run(fixed, tag):
    mods, load, REF, CORE3, EX = frame(fixed)
    if fixed:
        keep, names = keep_mask(C.var_names, nz)
    assign = {k: np.zeros(C.shape[0], int) for k in REF}
    for st in range(0, C.shape[0], CH):
        m = np.zeros(C.shape[0], bool); m[st:st + CH] = True
        a = C[m].to_memory()
        if fixed:
            a = a[:, keep].copy(); a.var_names = names[keep].values
        if a.X.max() > 100: sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        S = pct.score_modules(a, method='mean_raw', gene_set='core')['scores'] \
               .reindex(columns=mods, fill_value=0.0).values.astype(float)
        for k, ((mu, sd, rpn), tf) in REF.items():
            pc = ((tf(S) - mu) / sd) @ load
            pcn = pc / (np.linalg.norm(pc, axis=1, keepdims=True) + 1e-12)
            sim = pcn @ rpn.T; sim[:, EX] = -1e9
            assign[k][st:st + CH] = np.argmax(sim, axis=1)
        del a, S; gc.collect()
    o = obs.copy()
    o['core_RAW'] = np.isin(assign['RAW'], CORE3)
    o['core_SHAPE'] = np.isin(assign['SHAPE'], CORE3)
    rows = []
    for p, g in o.groupby('patient'):
        mi, ni = g[g.arm2 == 'malig'], g[g.arm2 == 'normal']
        if len(mi) >= 20 and len(ni) >= 20:
            rows.append(dict(patient=p, m_raw=100 * mi.core_RAW.mean(), n_raw=100 * ni.core_RAW.mean(),
                             m_shp=100 * mi.core_SHAPE.mean(), n_shp=100 * ni.core_SHAPE.mean()))
    R = pd.DataFrame(rows)
    out = {}
    for lab, mc, nc in (('RAW', 'm_raw', 'n_raw'), ('SHAPE', 'm_shp', 'n_shp')):
        mm, nn = R[mc], R[nc]; d = mm - nn
        out[lab] = dict(n=len(R), malig=mm.median(), normal=nn.median(), delta=d.mean(),
                        higher=int((d > 0).sum()), p=wilcoxon(mm, nn).pvalue,
                        sign_p=binomtest(int((d > 0).sum()), len(R), 0.5).pvalue)
        print(f'  [{tag}] {lab:6s} n={len(R)}  malig {mm.median():.3f}%  normal {nn.median():.3f}%  '
              f'delta {d.mean():+.3f}pp  higher {int((d>0).sum())}/{len(R)}  '
              f'W P={out[lab]["p"]:.3e}', flush=True)
    R.to_csv(FIX / f'E46_colorectal_{"fix" if fixed else "shipped"}.csv', index=False)
    return out

print('\nrebuilding the shipped arm (the gate)')
A = run(False, 'shipped')
T = pd.read_csv(ROOT / 'results/referee2_newdata/E36_tests.csv')
g = T[T.test.str.startswith('(i) all malignant vs all normal')].iloc[0]
for k, v in [('n', A['RAW']['n']), ('malig_median', A['RAW']['malig']),
             ('normal_median', A['RAW']['normal']), ('higher', A['RAW']['higher'])]:
    ok = abs(float(g[k]) - v) < 1e-6
    print(f'  GATE {k:14s} E36={g[k]:<12} rebuilt={v:<12} {"ok" if ok else "MISMATCH"}')
    assert ok, f'{k}: E36 {g[k]} vs rebuilt {v}'
print('  gate passed - the comparison below measures the fix\n')
print('rebuilding with the gene-name fix, the fixed reference and the fixed eigenspace')
B = run(True, 're-derived')
print('\n(iii) colorectal within-patient increase')
for lab in ('RAW', 'SHAPE'):
    a, b = A[lab], B[lab]
    print(f'  {lab:6s} {a["malig"]:.1f}% vs {a["normal"]:.1f}% ({a["higher"]}/{a["n"]}, P={a["p"]:.1e})'
          f'   ->   {b["malig"]:.1f}% vs {b["normal"]:.1f}% ({b["higher"]}/{b["n"]}, P={b["p"]:.1e})')
json.dump({'shipped': A, 're-derived': B}, open(FIX / 'E46_colorectal_compare.json', 'w'),
          indent=1, default=float)
print('DONE')
