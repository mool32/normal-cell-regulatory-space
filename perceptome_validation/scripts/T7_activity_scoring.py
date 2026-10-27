"""T7 — Activity-level scoring (readiness -> activity).
Every module is re-scored by regulatory OUTPUT instead of component abundance:
  TF-output modules with a DoRothEA A/B regulon of >= 5 targets (union over the module's primary TFs; signed weights)
      -> regulon activity (decoupler VIPER; ULM reported as the fast linear alternative)
  pathway-footprint modules without a TF output that PROGENy covers (ERK/MAPK -> MAPK, PI3K/PTEN -> PI3K)
      -> PROGENy footprint (top-100 weights, same estimator)
  everything else -> mean expression of the catalog's documented activity (target) genes
Reference (154 HPA cell types, both calibrations) and the 83,819 cells are scored identically, per cohort on that cohort's present
genes; the eigenspace is rebuilt from the activity-scored 154 x 44 matrix (pseudobulk loadings primary; single_cell_scaled secondary)
and the full pipeline re-run: E1 / E2 / E3 (Pelka) raw and shape-only.
Criteria: PASS if E1_3 >= 11.7 %, cyto rank <= 3, E3 shape-only same sign and >= 5 pp; FAIL if any breaks.
Directional prediction (one-sided): shape-only E1 / raw E1 rises under activity scoring relative to readiness scoring.
"""
import os, sys, json, time, warnings, gc
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prereg_lib import *  # noqa
import decoupler as dc
ESTIMATOR = sys.argv[1] if len(sys.argv) > 1 else 'viper'      # viper | ulm

hpa = hpa_gene_level(); cts = hpa['cts']; genes = hpa['genes']; pan = perceptome_panel(); mods = list(pan)
cache = load_panel_cache(); pelka = load_pelka_cache()
dor = pd.read_csv(os.path.join(PAN, 'dorothea_AB.csv')); dor['target'] = dor.target.str.upper()
prog = pd.read_csv(os.path.join(PAN, 'progeny_top100.csv')); prog['target'] = prog.target.str.upper()
PROGENY_MAP = {'ERK/MAPK': 'MAPK', 'PI3K/PTEN': 'PI3K'}
register_inputs('T7', [os.path.join(PAN, 'dorothea_AB.csv'), os.path.join(PAN, 'progeny_top100.csv'), os.path.join(CACHE, 'panel25.h5ad'), os.path.join(CACHE, 'pelka_epi.h5ad'), os.path.abspath(__file__)])

# ---- build the activity definition of every module ----
net_rows = []; mean_sets = {}; method = {}
for m in mods:
    info = pct.get_module_info(m); tfs = info.get('primary_tf') or []
    reg = dor[dor.source.isin(tfs)] if tfs else dor.iloc[0:0]
    reg = reg.drop_duplicates('target')
    if tfs and len(reg) >= 5:
        for _, r in reg.iterrows(): net_rows.append(dict(source=m, target=r.target, weight=float(r.weight)))
        method[m] = f'dorothea:{"+".join(tfs)} ({len(reg)} targets)'
    elif (not tfs) and m in PROGENY_MAP:
        pr = prog[prog.source == PROGENY_MAP[m]]
        for _, r in pr.iterrows(): net_rows.append(dict(source=m, target=r.target, weight=float(r.weight)))
        method[m] = f'progeny:{PROGENY_MAP[m]} ({len(pr)} genes)'
    else:
        mean_sets[m] = [str(g).upper() for g in (info.get('activity_genes') or [])]; method[m] = f'mean activity_genes ({len(mean_sets[m])})'
net = pd.DataFrame(net_rows)
print("activity definitions:"); [print(f"  {m:14s} {method[m]}") for m in mods]
json.dump(method, open(os.path.join(REG, 'T7_activity_definitions.json'), 'w'), indent=1, ensure_ascii=False)

def estimate(X, obs_names, var_names, present=None):
    """activity matrix (n x 44) in module order: decoupler estimator for regulon/footprint modules, mean for the rest."""
    X = sp.csr_matrix(X, dtype=np.float32) if sp.issparse(X) else np.asarray(X, np.float32)
    vn = np.array(var_names); keep = np.ones(len(vn), bool) if present is None else np.asarray(present, bool)
    Xk = X[:, keep]; vk = vn[keep]
    a = ad.AnnData(X=Xk if sp.issparse(Xk) else sp.csr_matrix(Xk)); a.obs_names = [str(o) for o in obs_names]; a.var_names = list(vk)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if ESTIMATOR == 'viper': dc.mt.viper(a, net, tmin=3, verbose=False)
        else: dc.mt.ulm(a, net, tmin=3, verbose=False)
    key = 'score_viper' if ESTIMATOR == 'viper' else 'score_ulm'
    sc_ = a.obsm[key]; out = np.zeros((a.n_obs, len(mods)), np.float32)
    gi = {g: i for i, g in enumerate(vk)}
    for j, m in enumerate(mods):
        if m in sc_.columns: out[:, j] = sc_[m].values
        elif m in mean_sets:
            ii = [gi[g] for g in mean_sets[m] if g in gi]
            if len(ii) >= 2: out[:, j] = np.asarray(Xk[:, ii].mean(1)).ravel()
    return out

t0 = time.time()
A_pb = estimate(hpa['H_pb'], cts, genes); A_sc = estimate(hpa['H_sc'], cts, genes); print(f"reference scored ({time.time()-t0:.0f}s)", flush=True)
BLOCK = 20000          # VIPER is memory-heavy: score large cohorts in blocks (a 127k-cell cohort OOMs a 16 GB machine)
def score_cache(c, cohort_key='cancer', tag=''):
    ck = os.path.join(CACHE, f'T7_{ESTIMATOR}_{tag}.npy')
    if os.path.exists(ck):
        A = np.load(ck); print(f"   [{tag}] loaded cached activity scores {A.shape}", flush=True); return A
    obs = c['obs']; cid = obs[cohort_key].astype(str).values; out = np.zeros((c['X'].shape[0], len(mods)), np.float32)
    for k in pd.unique(cid):
        rows = np.where(cid == k)[0]; t1 = time.time()
        for b in range(0, len(rows), BLOCK):
            r = rows[b:b + BLOCK]
            out[r] = estimate(c['X'][r], obs.index[r], c['genes'], present=c['present'][k]); gc.collect()
        print(f"   {k}: {len(rows)} cells ({time.time()-t1:.0f}s)", flush=True)
    np.save(ck, out); return out
A_cells = score_cache(cache, tag='panel25'); A_pelka = score_cache(pelka, 'cohort', tag='pelka')
np.savez_compressed(os.path.join(RES, f'T7_activity_scores_{ESTIMATOR}.npz'), A_pb=A_pb, A_sc=A_sc, A_cells=A_cells, A_pelka=A_pelka, mods=np.array(mods))

# ---- pipeline on activity scores ----
cid = cache['obs']['cancer'].astype(str).values; excl = np.array([cts.index(c) for c in EXCL])
def pipeline(A_ref_load, A_ref_sc, A_c, A_p, tag, tf=lambda X: X):
    es = eigenspace(tf(A_ref_load)); mu, sd, refPn = ref_frame(tf(A_ref_sc), es['L'])
    Pn = project(tf(A_c), mu, sd, es['L']); top1, _ = assign_top1(Pn, refPn, excl); ep = endpoints(top1, cid, cts, excl)
    # E3
    core_idx = np.array([cts.index(c) for c in CORE3]); arm = pelka['obs']['arm'].astype(str).values; pat = pelka['obs']['patient'].astype(str).values
    Pp = project(tf(A_p), mu, sd, es['L']); tp, _ = assign_top1(Pp, refPn, excl); isc = np.isin(tp, core_idx); rows = []
    for p in pd.unique(pat):
        mi = (pat == p) & (arm == 'malig'); ni = (pat == p) & (arm == 'normal')
        if mi.sum() < 20 or ni.sum() < 20: continue
        rows.append(dict(patient=p, malig=100 * isc[mi].mean(), normal=100 * isc[ni].mean()))
    d = pd.DataFrame(rows); dif = d.malig - d.normal
    from scipy.stats import wilcoxon
    e3 = dict(delta=float(dif.mean()), deeper=int((dif > 0).sum()), n=len(d), p=float(wilcoxon(d.malig, d.normal).pvalue) if len(d) > 5 else np.nan)
    r = dict(tag=tag, n_kaiser=es['n_kaiser'], E1_two=100 * ep['E1_2'], E1_three=100 * ep['E1_3'], cyto_rank=ep['cyto_rank_dense'], cyto_rec=ep['cyto_rec'], mega_rec=ep['mega_rec'], top5=ep['top5'], E3_delta=e3['delta'], E3_deeper=f"{e3['deeper']}/{e3['n']}", E3_p=e3['p'])
    print(f"[{tag:44s}] Kaiser {es['n_kaiser']:2d} | E1 two {r['E1_two']:5.1f}% three {r['E1_three']:5.1f}% | cyto {ep['cyto_rec']:2d}/25 rank {ep['cyto_rank_dense']:3d} | mega {ep['mega_rec']:2d}/25 | E3 {e3['delta']:+5.1f} ({e3['deeper']}/{e3['n']}) | top5 {ep['top5']}", flush=True)
    return r
demag = lambda X: X - X.mean(1, keepdims=True)
res = []
for calib, Aload in [('pseudobulk', A_pb), ('single_cell_scaled', A_sc)]:
    res.append(pipeline(Aload, A_sc, A_cells, A_pelka, f'ACTIVITY-{ESTIMATOR} | {calib} loadings | raw'))
    res.append(pipeline(Aload, A_sc, A_cells, A_pelka, f'ACTIVITY-{ESTIMATOR} | {calib} loadings | shape-only', tf=demag))
# readiness baseline for the directional prediction (same harness, both transforms)
S_cells = score_cells(cache, pan); S_pelka = score_cells(pelka, pan, cohort_key='cohort'); R_pb = score_ref(hpa['H_pb'], genes, pan); R_sc = score_ref(hpa['H_sc'], genes, pan)
res.append(pipeline(R_pb, R_sc, S_cells, S_pelka, 'READINESS | pseudobulk loadings | raw'))
res.append(pipeline(R_pb, R_sc, S_cells, S_pelka, 'READINESS | pseudobulk loadings | shape-only', tf=demag))
T = pd.DataFrame(res); T.to_csv(os.path.join(RES, f'T7_activity_table_{ESTIMATOR}.csv'), index=False)
# ---- verdict ----
raw = T[T.tag == f'ACTIVITY-{ESTIMATOR} | pseudobulk loadings | raw'].iloc[0]; shp = T[T.tag == f'ACTIVITY-{ESTIMATOR} | pseudobulk loadings | shape-only'].iloc[0]
rr = T[T.tag == 'READINESS | pseudobulk loadings | raw'].iloc[0]; rs = T[T.tag == 'READINESS | pseudobulk loadings | shape-only'].iloc[0]
ok = (raw.E1_three >= 11.7) and (raw.cyto_rank <= 3) and (shp.E3_delta >= 5.0) and (np.sign(shp.E3_delta) == np.sign(rs.E3_delta))
verdict = 'PASS' if ok else 'FAIL'
ratio_act = shp.E1_three / raw.E1_three if raw.E1_three > 0 else np.nan; ratio_read = rs.E1_three / rr.E1_three
pred = 'CONFIRMED' if ratio_act > ratio_read else 'REFUTED'
print("\n" + "=" * 100)
print(f"T7 VERDICT ({ESTIMATOR}; primary: pseudobulk loadings): {verdict}")
print(f"  activity raw: E1_3 {raw.E1_three:.1f}% (floor 11.7) | cyto rank {raw.cyto_rank} ({raw.cyto_rec}/25) | E3 shape-only {shp.E3_delta:+.1f} pp ({shp.E3_deeper}) [readiness shape-only {rs.E3_delta:+.1f}]")
print(f"  directional prediction (shape/raw E1 ratio rises under activity): readiness {ratio_read:.2f} -> activity {ratio_act:.2f}: {pred}")
write_verdict('T7', verdict, f"activity({ESTIMATOR}) E1_3={raw.E1_three:.1f}%, cyto rank {raw.cyto_rank}, E3 shape={shp.E3_delta:+.1f}pp",
              'PASS: E1_3>=11.7 & cyto rank<=3 & E3 shape >=5pp same sign', f"directional prediction {pred} (shape/raw ratio readiness {ratio_read:.2f} -> activity {ratio_act:.2f}); scscaled-loadings rows in table")
json.dump(dict(verdict=verdict, prediction=pred, ratio_readiness=ratio_read, ratio_activity=ratio_act, estimator=ESTIMATOR), open(os.path.join(RES, f'T7_verdict_{ESTIMATOR}.json'), 'w'), indent=1, default=float)
print('DONE')
