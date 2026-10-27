"""Shared harness for the perceptome validation pre-registration (T1–T7).

Everything here reproduces the paper-3 v10 pipeline byte-for-byte for the incumbent 44-module
catalog and generalises it to ANY feature panel {module: [genes]}:

  HPA gene level  ->  panel scores (154 x k) under BOTH calibrations
                      pseudobulk      : mean log1p(nCPM)                      (eigenspace loadings)
                      single_cell_scaled : mean log1p(normalize_total 1e4)    (mu/sd + reference centroids)
  eigenspace      ->  z-score per module -> correlation PCA -> Kaiser (recorded) -> 9 PCs forced
  cancers         ->  per-cell mean_raw over the panel genes present in that cohort (min 2 genes, else 0)
                      z vs single_cell_scaled mu/sd -> project -> cosine top-1 among 152 (erythrocyte/platelet quarantined)
  endpoints       ->  E1 core fraction (two-cell and three-cell), E2 recurrence ranking, E3 Pelka within-patient deepening

The 25-cancer panel is cached ONCE at gene level (cache/panel25.h5ad) with the exact cells of the paper
(seed 0, CAP 4000, cohort order of scripts/E3_referee_controls.py) so that every test scores the same
83,819 cells.  Gate script T0_gate.py verifies the reproduction against results/E7_S_raw.npz and the
bundled eigenspace before any endpoint is computed.
"""
import os, re, json, gc, hashlib, warnings, sys
import numpy as np, pandas as pd, scipy.sparse as sp
warnings.filterwarnings('ignore')

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
VAL = os.path.join(ROOT, 'perceptome_validation')
CACHE = os.path.join(VAL, 'cache'); RES = os.path.join(VAL, 'results'); REG = os.path.join(VAL, 'registry'); PAN = os.path.join(VAL, 'panels')
for d in (CACHE, RES, REG, PAN): os.makedirs(d, exist_ok=True)
sys.path.insert(0, os.path.join(ROOT, 'tool/perceptome2'))
import perceptome as pct  # noqa: E402
from perceptome.perceptivity.reference import load_hpa_perceptivity  # noqa: E402

CORE2 = ['cytotrophoblasts', 'megakaryocytes']
CORE3 = ['cytotrophoblasts', 'megakaryocytes', 'migrating cytotrophoblasts']
EXCL = ['erythrocytes', 'platelets']
NPC = 9; CAP = 4000; MIN_GENES_CELL = 1000; SEED = 0

def stage(m): print(f"[stage] {m}", flush=True)
def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''): h.update(chunk)
    return h.hexdigest()

# ----------------------------------------------------------------------------------------------
# Panels
# ----------------------------------------------------------------------------------------------
def load_panel(name):
    """panels/<name>.json -> dict {module: [GENES]} (uppercased, order preserved, de-duplicated)."""
    p = json.load(open(os.path.join(PAN, f'{name}.json')))
    return {m: list(dict.fromkeys(str(g).upper() for g in gl)) for m, gl in p.items()}

def perceptome_panel():
    return {m: [str(g).upper() for g in pct.get_genes(m)] for m in pct.list_modules()}

# ----------------------------------------------------------------------------------------------
# HPA gene level (both calibrations)
# ----------------------------------------------------------------------------------------------
_HPA = None
def hpa_gene_level():
    """dict(genes, cts, H_pb, H_sc).  cts = the 154 reference cell types in bundled order.
    H_pb = log1p(nCPM); H_sc = log1p(normalize_total(nCPM, 1e4)).  Duplicated gene names averaged."""
    global _HPA
    if _HPA is not None: return _HPA
    path = os.path.join(CACHE, 'hpa_gene_level.npz')
    meta = json.load(open(os.path.join(ROOT, 'tool/perceptome2/perceptome/perceptivity/data/hpa_perceptivity_v03_scscaled.json')))
    cts = list(meta['cell_types'])
    if os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        _HPA = dict(genes=list(d['genes']), cts=list(d['cts']), H_pb=d['H_pb'], H_sc=d['H_sc'])
        assert _HPA['cts'] == cts
        return _HPA
    stage('building HPA gene-level matrices from rna_single_cell_type.tsv')
    df = pd.read_csv(os.path.join(ROOT, 'results/hpa_download/rna_single_cell_type.tsv'), sep='\t', usecols=['Gene name', 'Cell type', 'nCPM'])
    df['Gene name'] = df['Gene name'].astype(str).str.upper()
    H = df.pivot_table(index='Gene name', columns='Cell type', values='nCPM', aggfunc='mean').reindex(columns=cts).fillna(0.0)
    N = H.values.astype(np.float64).T                     # 154 x genes, nCPM
    H_pb = np.log1p(N)
    H_sc = np.log1p(N / (N.sum(1, keepdims=True) + 1e-12) * 1e4)
    genes = list(H.index)
    np.savez_compressed(path, genes=np.array(genes), cts=np.array(cts), H_pb=H_pb.astype(np.float32), H_sc=H_sc.astype(np.float32))
    _HPA = dict(genes=genes, cts=cts, H_pb=H_pb.astype(np.float32), H_sc=H_sc.astype(np.float32))
    return _HPA

def score_ref(H, genes, panel):
    """HPA reference scoring, exactly as scripts/12_build_hpa_scscaled.py: mean over the panel genes present
    in HPA (>=1 gene; 0 if none).  Returns (154 x k) float64 in panel order."""
    gi = {g: i for i, g in enumerate(genes)}
    out = np.zeros((H.shape[0], len(panel)))
    for j, (m, gl) in enumerate(panel.items()):
        ii = [gi[g] for g in gl if g in gi]
        if ii: out[:, j] = H[:, ii].astype(np.float64).mean(1)
    return out

# ----------------------------------------------------------------------------------------------
# Gene universe for the cancer cache
# ----------------------------------------------------------------------------------------------
def gene_universe():
    """HPA genes  U  every gene of every frozen panel/resource (catalog core+activity, hallmark, progeny, reactome,
    kegg, DoRothEA A/B, PROGENy-500, trophoblast markers).  Genes outside the universe are dropped from the
    cache; they can never enter any module score in the tests below."""
    U = set(hpa_gene_level()['genes'])
    for m in pct.list_modules():
        U |= {str(g).upper() for g in pct.get_genes(m)}
        info = pct.get_module_info(m)
        for key in ('activity_genes', 'sensor_genes', 'cascade_genes', 'tf_genes', 'feedback_genes'):
            U |= {str(g).upper() for g in (info.get(key) or [])}
    for f in os.listdir(PAN):
        if f.endswith('.json'):
            p = json.load(open(os.path.join(PAN, f)))
            for gl in p.values(): U |= {str(g).upper() for g in gl}
        elif f.endswith('.csv'):
            d = pd.read_csv(os.path.join(PAN, f))
            for c in ('source', 'target'):
                if c in d.columns: U |= set(d[c].astype(str).str.upper())
    U |= {'MMP2', 'MMP9', 'ITGA5', 'HLA-G', 'CGA', 'CGB3', 'KRT7', 'TFAP2C', 'GATA3', 'PAGE4', 'XAGE2', 'XAGE3', 'MB21D1', 'TMEM173', 'ARNTL'}
    return sorted(U)

# ----------------------------------------------------------------------------------------------
# Cohort definitions (identical to scripts/E3_referee_controls.py, including order -> same rng draws)
# ----------------------------------------------------------------------------------------------
def cohort_table():
    C = os.path.join(ROOT, 'results'); CACHE_ = f'{C}/cancer_scaled/cache'; V63 = f'{C}/cancer_v6_3_extras'; V64 = f'{C}/cancer_v6_4_extras'
    V65 = f'{C}/cancer_v6_5_extras'; OUT = f'{C}/cancer_v6_6_extras'
    BCC = lambda o: ((o['disease'].astype(str) == 'basal cell carcinoma') & (o['cell_type'].astype(str).isin(['basal cell', 'basal cell of epidermis', 'suprabasal keratinocyte']))).values
    return [('HCC', f'{C}/cancer_sun2021_hcc/sun2021.h5ad', lambda o: (o['celltype'] == 'Hepatocyte').values),
            ('BRCA', f'{C}/cancer_validation/wu2021_brca.h5ad', None), ('Glioblastoma', f'{CACHE_}/Glioblastoma_malignant.h5ad', None),
            ('OvarianSerous', f'{CACHE_}/OvarianSerous_malignant.h5ad', None), ('RCC', f'{CACHE_}/RCC_clearcell_malignant.h5ad', None),
            ('Melanoma', f'{CACHE_}/Melanoma_malignant.h5ad', None), ('Prostate', f'{CACHE_}/Prostate_malignant.h5ad', None),
            ('HNSCC', f'{V63}/HNSCC_malignant.h5ad', None), ('Neuroblastoma', f'{V64}/Neuroblastoma_malignant.h5ad', None),
            ('B-ALL', f'{V64}/B-ALL_malignant.h5ad', None), ('DLBCL', f'{V64}/DLBCL_malignant.h5ad', None),
            ('SqLung', f'{V64}/SqLung_malignant.h5ad', None), ('TNBC', f'{V64}/TNBC_malignant.h5ad', None),
            ('PDAC', f'{V65}/PDAC_Peng2019_malignant.h5ad', None),
            ('NSCLC', f'{ROOT}/paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad', lambda o: (o['cell_type'].astype(str) == 'malignant cell').values),
            ('AML', f'{OUT}/vanGalen_AML_malignant.h5ad', None), ('MM', f'{OUT}/Ledergor_MM_malignant.h5ad', None),
            ('BCC', f'{ROOT}/paper6_immune_3layer/data/bcc_full.h5ad', BCC),
            ('Colon', f'{C}/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad', lambda o: (o['ClusterMidway'].astype(str) == 'EpiT').values),
            ('SCLC', f'{OUT}/SCLC_malignant.h5ad', None), ('Gastric', f'{OUT}/Gastric_malignant.h5ad', None),
            ('GermCell', f'{OUT}/GermCell_malignant.h5ad', None), ('Osteosarcoma', f'{OUT}/Osteosarcoma_malignant.h5ad', None),
            ('Ewing', f'{OUT}/Ewing_malignant.h5ad', None), ('Seminoma', f'{OUT}/Seminoma_malignant.h5ad', None)]

def var_keep_mask(a):
    """The paper's nv(): symbol names (feature_name if ENSG index), uppercase, make_unique, drop '-N' suffix duplicates.
    Returns (keep_mask over a.var, base symbol names of the kept columns)."""
    v0 = str(a.var.index[0])
    if 'feature_name' in a.var.columns and (v0.startswith('ENSG') or not any(c.isalpha() for c in v0)):
        names = pd.Index(a.var['feature_name'].astype(str).values)
    else:
        names = pd.Index(a.var_names.astype(str))
    names = names.str.upper()
    names = pd.Index(pd.Series(names).pipe(lambda s: _make_unique(s)))
    base = pd.Series([re.sub(r'-\d+$', '', str(v)) for v in names])
    k = ~base.duplicated().values
    return k, base[k].values

def _make_unique(s):
    """replicate anndata var_names_make_unique (append -1, -2 ... to later duplicates)."""
    seen = {}; out = []
    counts = s.value_counts()
    for v in s:
        if counts[v] > 1:
            n = seen.get(v, 0)
            out.append(v if n == 0 else f'{v}-{n}'); seen[v] = n + 1
        else: out.append(v)
    return pd.Series(out)

# ----------------------------------------------------------------------------------------------
# Cancer cache
# ----------------------------------------------------------------------------------------------
def load_panel_cache():
    """The 25-cancer panel at gene level.  Returns dict(X csr float32 [cells x universe], genes, obs DataFrame
    [cancer, n_genes, src_row], present {cancer: bool array over genes})."""
    import anndata as ad
    a = ad.read_h5ad(os.path.join(CACHE, 'panel25.h5ad'))
    present = {k: np.asarray(v, bool) for k, v in a.uns['present'].items()}
    return dict(X=sp.csr_matrix(a.X, dtype=np.float32), genes=list(a.var_names), obs=a.obs.copy(), present=present)

def score_cells(cache, panel, cohort_key='cancer', min_genes=2):
    """Per-cell mean_raw scoring of a panel, byte-identical to pct.score_modules(method='mean_raw'):
    mean over the panel genes PRESENT in the cohort's var_names (present after nv()); <min_genes present -> 0."""
    X, genes, obs, present = cache['X'], cache['genes'], cache['obs'], cache['present']
    gi = {g: i for i, g in enumerate(genes)}
    k = len(panel); S = np.zeros((X.shape[0], k), np.float32)
    cid = obs[cohort_key].astype(str).values
    for c in pd.unique(cid):
        rows = np.where(cid == c)[0]; pres = present[c]
        r_, c_, v_ = [], [], []
        for j, (m, gl) in enumerate(panel.items()):
            ii = [gi[g] for g in gl if g in gi and pres[gi[g]]]
            if len(ii) < min_genes: continue
            r_ += ii; c_ += [j] * len(ii); v_ += [1.0 / len(ii)] * len(ii)
        G = sp.csr_matrix((v_, (r_, c_)), shape=(len(genes), k), dtype=np.float32)
        S[rows] = np.asarray((X[rows] @ G).todense(), np.float32)
    return S

# ----------------------------------------------------------------------------------------------
# Eigenspace + projection + endpoints
# ----------------------------------------------------------------------------------------------
def eigenspace(R, n_pc=NPC):
    """Canonical rebuild(): z-score (ddof=0; sd<1e-10 -> 1) -> corrcoef -> eigh -> descending; returns
    dict(mean, std, L [k x n_pc], eigvals, n_kaiser, eigvecs_full)."""
    X = np.asarray(R, float); mean = X.mean(0); std = X.std(0); std[std < 1e-10] = 1.0
    Z = (X - mean) / std
    corr = np.corrcoef(Z.T); corr = np.nan_to_num(corr)
    ev, vec = np.linalg.eigh(corr); idx = np.argsort(ev)[::-1]; ev = ev[idx]; vec = vec[:, idx]
    return dict(mean=mean, std=std, L=vec[:, :n_pc], eigvals=ev, n_kaiser=int((ev > 1).sum()), eigvecs=vec, var_explained=ev / ev.sum())

def bundled_eigenspace():
    e = json.load(open(os.path.join(ROOT, 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')))
    mods = e['module_order']; L = np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(NPC)] for m in mods])
    return mods, L, e

def ref_frame(Rsc, L):
    """single_cell_scaled reference frame: mu, sd (ddof=0, +1e-12) and unit centroids in the eigenspace."""
    Rsc = np.asarray(Rsc, float); mu = Rsc.mean(0); sd = Rsc.std(0) + 1e-12
    P = ((Rsc - mu) / sd) @ L; Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    return mu, sd, Pn

def project(S, mu, sd, L):
    P = ((np.asarray(S, float) - mu) / sd) @ L
    return P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)

def assign_top1(Pn, refPn, excl_idx, chunk=20000):
    n = Pn.shape[0]; top1 = np.empty(n, int); dmin = np.empty(n, float)
    for i in range(0, n, chunk):
        sim = Pn[i:i + chunk] @ refPn.T; sim[:, excl_idx] = -1e9
        t = np.argmax(sim, 1); top1[i:i + chunk] = t; dmin[i:i + chunk] = 1 - sim[np.arange(len(t)), t]
    return top1, dmin

def recurrence(top1, cid, n_ref, excl_idx, topn=5):
    """# cancers in which each cell type is among the topn most frequent top-1 targets (excluded cells -> -1)."""
    rec = np.zeros(n_ref)
    for c in pd.unique(cid):
        cnt = np.bincount(top1[cid == c], minlength=n_ref); rec[np.argsort(-cnt, kind='stable')[:topn]] += 1
    rec[excl_idx] = -1
    return rec

def rank_of(rec, refcells, cell):
    """1-based rank in the recurrence ranking; ties broken by stable argsort (same convention as E2e)."""
    order = np.argsort(-rec, kind='stable'); return int(np.where(order == refcells.index(cell))[0][0]) + 1

def rank_dense(rec, refcells, cell):
    """rank counting only cells with STRICTLY higher recurrence + 1 (tie-generous; reported alongside)."""
    v = rec[refcells.index(cell)]; return int((rec > v).sum()) + 1

def endpoints(top1, cid, refcells, excl_idx):
    core2 = [refcells.index(c) for c in CORE2]; core3 = [refcells.index(c) for c in CORE3]
    rec = recurrence(top1, cid, len(refcells), excl_idx)
    order = np.argsort(-rec, kind='stable')
    depth_med = {}
    for i in range(len(refcells)):
        depth_med[i] = float(np.median([np.mean(top1[cid == c] == i) for c in pd.unique(cid)]))
    return dict(E1_2=float(np.isin(top1, core2).mean()), E1_3=float(np.isin(top1, core3).mean()),
                cyto_rank=rank_of(rec, refcells, 'cytotrophoblasts'), cyto_rank_dense=rank_dense(rec, refcells, 'cytotrophoblasts'),
                cyto_rec=int(rec[refcells.index('cytotrophoblasts')]), mega_rec=int(rec[refcells.index('megakaryocytes')]),
                mega_rank=rank_of(rec, refcells, 'megakaryocytes'), migr_rec=int(rec[refcells.index('migrating cytotrophoblasts')]),
                top5=[(refcells[i], int(rec[i]), round(100 * depth_med[i], 2)) for i in order[:5]],
                cyto_depth_median=100 * depth_med[refcells.index('cytotrophoblasts')], rec=rec)

def run_panel(panel, cache, hpa=None, loadings_calibration='pseudobulk', n_pc=NPC, transform=None, return_all=False):
    """Full pipeline for one panel.  loadings_calibration: 'pseudobulk' (paper: eigenspace built on pseudobulk R,
    projection under single_cell_scaled) or 'single_cell_scaled' (E2e-style: everything on single_cell_scaled).
    transform: optional f(matrix)->matrix applied to BOTH reference and cancer module scores (magnitude removal)."""
    hpa = hpa or hpa_gene_level(); cts = hpa['cts']
    R_pb = score_ref(hpa['H_pb'], hpa['genes'], panel); R_sc = score_ref(hpa['H_sc'], hpa['genes'], panel)
    S = score_cells(cache, panel)
    if transform is not None: R_pb, R_sc, S = transform(R_pb), transform(R_sc), transform(S)
    es = eigenspace(R_pb if loadings_calibration == 'pseudobulk' else R_sc, n_pc=n_pc)
    mu, sd, refPn = ref_frame(R_sc, es['L'])
    excl = np.array([cts.index(c) for c in EXCL])
    Pn = project(S, mu, sd, es['L']); top1, dmin = assign_top1(Pn, refPn, excl)
    cid = cache['obs']['cancer'].astype(str).values
    ep = endpoints(top1, cid, cts, excl); ep['n_kaiser'] = es['n_kaiser']; ep['n_cells'] = int(len(top1))
    ep['eigvals_top'] = [round(float(v), 3) for v in es['eigvals'][:12]]
    if return_all: ep.update(dict(top1=top1, dmin=dmin, Pn=Pn, refPn=refPn, mu=mu, sd=sd, L=es['L'], R_pb=R_pb, R_sc=R_sc, S=S, es=es))
    return ep

# ----------------------------------------------------------------------------------------------
# E3 — Pelka within-patient deepening (raw + shape-only), template scripts/E17_deepening_shapeonly.py
# ----------------------------------------------------------------------------------------------
def load_pelka_cache():
    import anndata as ad
    a = ad.read_h5ad(os.path.join(CACHE, 'pelka_epi.h5ad'))
    present = {k: np.asarray(v, bool) for k, v in a.uns['present'].items()}
    return dict(X=sp.csr_matrix(a.X, dtype=np.float32), genes=list(a.var_names), obs=a.obs.copy(), present=present)

def e3_deepening(panel, pelka, R_sc, L, core=CORE3, min_cells=20):
    """Per-patient core fraction, malignant (EpiT) vs matched normal (Epi), RAW and SHAPE-ONLY (within-cell centering of
    both reference and cells), as in E17.  Returns DataFrame [metric, patient, malig, normal] + summary dict."""
    from scipy.stats import wilcoxon, binomtest
    cts = hpa_gene_level()['cts']; excl = np.array([cts.index(c) for c in EXCL]); core_idx = np.array([cts.index(c) for c in core])
    S = score_cells(pelka, panel, cohort_key='cohort')
    arm = pelka['obs']['arm'].astype(str).values; pat = pelka['obs']['patient'].astype(str).values
    rows = []
    for tag, tf in [('RAW', lambda X: X), ('SHAPE-ONLY', lambda X: X - X.mean(1, keepdims=True))]:
        mu, sd, refPn = ref_frame(tf(np.asarray(R_sc, float)), L)
        Pn = project(tf(S.astype(float)), mu, sd, L); top1, _ = assign_top1(Pn, refPn, excl); isc = np.isin(top1, core_idx)
        for p in pd.unique(pat):
            mi = (pat == p) & (arm == 'malig'); ni = (pat == p) & (arm == 'normal')
            if mi.sum() < min_cells or ni.sum() < min_cells: continue
            rows.append(dict(metric=tag, patient=p, malig=100 * isc[mi].mean(), normal=100 * isc[ni].mean()))
    df = pd.DataFrame(rows); summ = {}
    for tag in ['RAW', 'SHAPE-ONLY']:
        d = df[df.metric == tag]; dif = d.malig.values - d.normal.values
        try: pw = float(wilcoxon(d.malig, d.normal).pvalue)
        except Exception: pw = float('nan')
        summ[tag] = dict(delta_pp=float(dif.mean()), delta_median_pp=float(np.median(dif)), deeper=int((dif > 0).sum()), n=int(len(d)),
                         malig_median=float(d.malig.median()), normal_median=float(d.normal.median()), wilcoxon_p=pw,
                         sign_p=float(binomtest(int((dif > 0).sum()), int((dif != 0).sum()), 0.5, alternative='greater').pvalue) if (dif != 0).sum() else float('nan'))
    return df, summ

# ----------------------------------------------------------------------------------------------
# Registry helpers
# ----------------------------------------------------------------------------------------------
def write_verdict(test_id, verdict, primary_statistic, threshold, notes):
    import datetime
    p = os.path.join(REG, 'verdicts.tsv')
    if not os.path.exists(p): open(p, 'w').write('test_id\tdate\tverdict\tprimary_statistic\tthreshold\tnotes\n')
    open(p, 'a').write(f"{test_id}\t{datetime.date.today().isoformat()}\t{verdict}\t{primary_statistic}\t{threshold}\t{notes}\n")

def register_inputs(test_id, files, extra=None):
    d = {os.path.relpath(f, ROOT): sha256(f) for f in files if os.path.exists(f)}
    if extra: d.update(extra)
    json.dump(d, open(os.path.join(REG, f'inputs_{test_id}.json'), 'w'), indent=1)
    return d
