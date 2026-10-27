#!/usr/bin/env python3
"""T9 — human cells transformed in vitro from a KNOWN precursor, projected into the paper's frozen frame.

Pre-registration: perceptome_validation/PREREG_T9_invitro.md (frozen 2026-10-19, before any value was read).
Attractor + region definitions: perceptome_validation/results/T9_attractor.json (frozen first, separately).

Every statistic is a WITHIN-SERIES paired delta (gate G4).  Nothing is refit: eigenspace, loadings, mu, sd
and the 152 reference centroids are the shipped ones, recomputed here and asserted against the freeze.

Usage:  python3 -I scripts/T9_invitro_transformation.py [GSE...]
"""
import os, sys, json, gzip, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'perceptome_validation', 'scripts'))
import prereg_lib as L  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'invitro')
PLAT = os.path.join(DATA, 'platforms')
RES = os.path.join(ROOT, 'perceptome_validation', 'results')
FREEZE = os.path.join(RES, 'T9_attractor.json')
RNG = np.random.default_rng(0)

PROLIF = ['MKI67', 'CCNB1', 'CCNA2', 'TOP2A', 'PLK1', 'BUB1']
SENESC = ['CDKN2A', 'CDKN1A', 'SERPINE1', 'IL6', 'CXCL8', 'GLB1']


# ----------------------------------------------------------------------------- gene-length map
_GC = None
def gencode():
    """ENSG (unversioned) and SYMBOL -> exonic bp length, from the recount3 Gencode v26 gene annotation."""
    global _GC
    if _GC is not None:
        return _GC
    by_ens, by_sym = {}, {}
    with gzip.open(os.path.join(PLAT, 'gencode_G026_genes.gtf.gz'), 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'gene':
                continue
            try:
                ln = float(f[5])
            except ValueError:
                continue
            m = re.search(r'gene_id "([^".]+)', f[8]); s = re.search(r'gene_name "([^"]+)"', f[8])
            if m:
                by_ens[m.group(1)] = ln
            if s:
                by_sym[s.group(1).upper()] = max(by_sym.get(s.group(1).upper(), 0), ln)
    _GC = (by_ens, by_sym)
    return _GC


# ----------------------------------------------------------------------------- loaders
def collapse(df, how='max'):
    """symbol-indexed frame, duplicate symbols collapsed (max; 'mean' is the pre-specified sensitivity)."""
    df = df[~df.index.isna()]
    df.index = pd.Index([str(i).upper() for i in df.index])
    return (df.groupby(level=0).max() if how == 'max' else df.groupby(level=0).mean())


def load_array(gse, gpl, how='max'):
    """GEO series matrix -> (linear symbol x sample frame, sample titles).  Log-scaled matrices are
    exponentiated (rule: max < 30 => log2)."""
    path = os.path.join(DATA, gse, f'{gse}_series_matrix.txt.gz')
    titles, rows, hdr = None, [], None
    with gzip.open(path, 'rt', errors='replace') as fh:
        intab = False
        for line in fh:
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                intab = True
            elif line.startswith('!series_matrix_table_end'):
                break
            elif intab:
                f = line.rstrip('\n').split('\t')
                if hdr is None:
                    hdr = [x.strip('"') for x in f[1:]]
                else:
                    rows.append((f[0].strip('"'), f[1:]))
    X = pd.DataFrame({k: v for k, v in [(p, np.array([np.nan if x in ('', 'null', 'NA') else float(x)
                                                      for x in vals], float)) for p, vals in rows]}).T
    X.columns = hdr
    if np.nanmax(X.values) < 30:
        X = np.power(2.0, X)
    X = X.clip(lower=0).fillna(0.0)
    # probe -> symbol
    if gpl == 'GPL9188':
        ez = pd.read_csv(os.path.join(PLAT, 'GPL9188_table.tsv'), sep='\t', dtype=str,
                         low_memory=False).set_index('ID')['Gene_ID']
        sym = _entrez_to_symbol()
        ann = ez.dropna().map(lambda e: sym.get(str(e).split('.')[0]))
    elif gpl == 'GPL15207':
        ann = pd.read_csv(os.path.join(PLAT, 'GPL15207_table.tsv'), sep='\t', usecols=['ID', 'Gene Symbol'],
                          dtype=str, low_memory=False).set_index('ID')['Gene Symbol']
    else:
        p = os.path.join(PLAT, f'{gpl}.annot.gz')
        with gzip.open(p, 'rt', errors='replace') as fh:
            skip = 0
            for line in fh:
                skip += 1
                if line.startswith('!platform_table_begin'):
                    break
        ann = pd.read_csv(p, sep='\t', skiprows=skip, dtype=str, low_memory=False,
                          usecols=['ID', 'Gene symbol']).set_index('ID')['Gene symbol']
    ann = ann.dropna()
    ann = ann[~ann.astype(str).str.contains('///')]          # ambiguous probes dropped
    X = X.loc[X.index.intersection(ann.index)]
    X.index = ann.reindex(X.index).values
    X = collapse(X, how)
    X.columns = titles if titles and len(titles) == X.shape[1] else hdr
    return X


_E2S = None
def _entrez_to_symbol():
    """Entrez gene id -> HGNC symbol, read off the GPL570 annotation (no extra download)."""
    global _E2S
    if _E2S is not None:
        return _E2S
    p = os.path.join(PLAT, 'GPL570.annot.gz')
    with gzip.open(p, 'rt', errors='replace') as fh:
        skip = 0
        for line in fh:
            skip += 1
            if line.startswith('!platform_table_begin'):
                break
    a = pd.read_csv(p, sep='\t', skiprows=skip, dtype=str, low_memory=False,
                    usecols=['Gene ID', 'Gene symbol']).dropna()
    a = a[~a['Gene ID'].str.contains('///')]
    _E2S = dict(zip(a['Gene ID'].astype(str), a['Gene symbol'].str.upper()))
    return _E2S


def load_exp_dir(gse, pattern='_rna', how='max'):
    """GSE65885-style: one file per sample, columns symbol/refseq/chr/strand/value, no header."""
    d = os.path.join(DATA, gse, 'ext')
    cols = {}
    for f in sorted(os.listdir(d)):
        if pattern not in f or not f.endswith('.gz'):
            continue
        t = pd.read_csv(os.path.join(d, f), sep='\t', header=None,
                        names=['sym', 'refseq', 'chr', 'strand', 'val'], low_memory=False)
        nm = re.sub(r'^GSM\d+_|\.exp\.txt\.gz$', '', f)
        cols[nm] = t.groupby(t['sym'].astype(str).str.upper())['val'].max()
    X = pd.DataFrame(cols).fillna(0.0)
    _, by_sym = gencode()
    ln = np.array([by_sym.get(s, np.nan) for s in X.index], float)
    keep = ~np.isnan(ln)
    X = X.loc[keep].div(ln[keep], axis=0)
    return collapse(X, how)


def load_recount3(gse, how='max'):
    """recount3 gene_sums (coverage counts) -> length-normalised symbol x sample frame, columns = GSM titles."""
    d = os.path.join(DATA, gse)
    X = pd.read_csv(os.path.join(d, 'gene_sums.G026.gz'), sep='\t', comment='#', index_col=0)
    by_ens, _ = gencode()
    base = pd.Index([g.split('.')[0] for g in X.index])
    ln = np.array([by_ens.get(g, np.nan) for g in base], float)
    keep = ~np.isnan(ln)
    X = X.loc[keep]; ln = ln[keep]
    X = X.div(ln, axis=0)                                     # coverage / bp  ~ molar abundance
    meta = pd.read_csv(os.path.join(d, 'meta.MD.gz'), sep='\t', low_memory=False)
    name = {}
    for _, r in meta.iterrows():
        t = str(r.get('experiment_title', ''))
        m = re.match(r'GSM\d+:\s*(.*?)\s+RNA-seq', t)
        name[str(r['external_id'])] = (m.group(1) if m else t).strip()
    X.columns = [name.get(c, c) for c in X.columns]
    sym = pd.Index([_ens_symbol(g) for g in base[keep]])
    X.index = sym
    return collapse(X, how)


def _ens_symbol(ens):
    _, by_sym = gencode()
    if not hasattr(_ens_symbol, 'inv'):
        inv = {}
        with gzip.open(os.path.join(PLAT, 'gencode_G026_genes.gtf.gz'), 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 9 or f[2] != 'gene':
                    continue
                m = re.search(r'gene_id "([^".]+)', f[8]); s = re.search(r'gene_name "([^"]+)"', f[8])
                if m and s:
                    inv[m.group(1)] = s.group(1).upper()
        _ens_symbol.inv = inv
    return _ens_symbol.inv.get(ens, ens)


def load_counts(gse, fname, sep=',', gene_col=None, ens_col=None, how='max'):
    """plain count table -> length-normalised symbol x sample frame."""
    p = os.path.join(DATA, gse, fname)
    X = pd.read_csv(p, sep=sep, low_memory=False)
    by_ens, by_sym = gencode()
    if ens_col is None and X[X.columns[0]].astype(str).str.startswith('ENSG').mean() > 0.5:
        ens_col = X.columns[0]                               # first column is a versioned Ensembl id
    if ens_col and ens_col in X.columns:
        base = pd.Index([str(g).split('.')[0] for g in X[ens_col]])
        ln = np.array([by_ens.get(g, np.nan) for g in base], float)
        sym = pd.Index([str(s).upper() for s in X[gene_col]]) if gene_col else pd.Index(
            [_ens_symbol(g) for g in base])
        drop = [c for c in (ens_col, gene_col) if c and c in X.columns]
        V = X.drop(columns=drop)
    else:
        sym = pd.Index([str(s).upper() for s in X[X.columns[0]]])
        ln = np.array([by_sym.get(s, np.nan) for s in sym], float)
        V = X.drop(columns=[X.columns[0]])
    V = V.apply(pd.to_numeric, errors='coerce')
    keep = ~np.isnan(ln) & V.notna().all(1).values
    V = V.loc[keep]; V.index = sym[keep]
    V = V.div(ln[keep], axis=0)
    V.columns = [re.sub(r'\.\d+$', '', str(c)) for c in V.columns]
    return collapse(V, how)


# ----------------------------------------------------------------------------- series registry (arms)
def arms_118206(cols):
    a = {}
    for c in cols:
        s = c.lower()
        if re.search(r'basal[123]$', s):            a[c] = 'prec_prostate'
        elif re.search(r'nhbe[145]$', s):           a[c] = 'prec_lung'
        elif 'nhbe' in s and 'parcb' in s:          a[c] = 'tr_lung'
        elif re.search(r'parcb\d', s):              a[c] = 'tr_prostate'
        elif re.match(r'.*\b(acb|pacb|arcb)\d', s): a[c] = 'ctrl_dropout'
        else:                                        a[c] = 'ref_line'
    return a


def arms_6885(cols):
    f = lambda c: ('prec_HME' if c.startswith('HME ') else 'prec_BPE' if c.startswith('BPE ')
                   else 'tr_HME' if c.startswith('HMLER') else 'tr_BPE' if c.startswith('BPLER') else 'other')
    return {c: f(c) for c in cols}


def arms_22573(cols):
    f = lambda c: ('ctrl_U0126' if 'U0126' in c else 'tr_cdk4RAS' if 'cdk4 Ras' in c
                   else 'prec_LacZ' if 'LacZ' in c else 'other')
    return {c: f(c) for c in cols}


def arms_112369(cols):
    def f(c):
        s = c.lower()
        if 'engineered control' in s:                       return 'ctrl_eng_colon' if 'colon' in s else 'ctrl_eng_stomach'
        if 'colon' in s and 'knockout' in s:                return 'tr_colon'
        if 'normal stomach' in s and 'knockout' in s:       return 'tr_stomach'
        if s.startswith('normal stomach'):                  return 'prec_stomach'
        if 'normal counterpart' in s:                       return 'prec_counterpart'
        if 'normal like' in s:                              return 'other_normallike'
        if 'intestinal metaplasia' in s:                    return 'other_IM'
        if 'gastric cancer organoid' in s:                  return 'tumour_patient'
        return 'other'
    return {c: f(c) for c in cols}


def arms_57965(cols):
    def f(c):
        s = c.lower()
        if s.startswith('adenoma'):              return 'adenoma'
        if 'modifed adenoma' in s or 'modified adenoma' in s: return 'tr_mod_adenoma'
        if 'genome-engineered' in s:             return 'tr_engineered'
        if s.startswith('crc'):                  return 'tumour_CRC'
        return 'other'
    return {c: f(c) for c in cols}


def arms_3151(cols):
    def f(c):
        s = c.upper()
        if s.startswith('GFP'):  return 'prec_GFP'
        for k in ('RAS', 'MYC', 'SRC', 'E2F3', 'BCAT'):
            if s.startswith(k):  return f'tr_{k}'
        return 'other'
    return {c: f(c) for c in cols}


def arms_172029(cols):
    def f(c):
        s = c.upper()
        if not s.startswith(('WT_', 'KO_')):  return 'other'
        geno, stage, cond = s.split('_')[:3]
        return f"{'prec' if cond == 'CONTROL' else 'tr'}_{geno}{stage}"
    return {c: f(c) for c in cols}


def arms_231776(cols):
    def f(c):
        s = c.upper()
        if s == 'EP':             return 'prec_BJ'
        if s.startswith('HTERT'): return 'tr_hTERT'
        if s.startswith('SV40'):  return 'tr_SV40'
        if s.startswith('HRAS'):  return 'tr_HRAS'
        return 'other'
    return {c: f(c) for c in cols}


def arms_17941(cols):
    def f(c):
        s = c.replace('-', '').upper()
        if 'MCF10A' in s:
            late = ('24 HOUR' in s) or ('36 HOUR' in s)
            if 'ERSRC' in s and 'TAM' in s:    return 'tr_MCF10A' if late else 'other_early'
            if 'PBABE' in s and 'TAM' in s:    return 'prec_MCF10A' if late else 'other_early'
            if 'ERSRC' in s and 'ETOH' in s:   return 'ctrl_MCF10A_noTAM'
            return 'other'
        if 'BJ FIBROBLAST' in s:
            if ', ELR,' in s or s.endswith('ELR'):  return 'tr_BJ'
            if ', EL,' in s or s.endswith('EL'):    return 'prec_BJ'
            if ', EH,' in s or s.endswith('EH'):    return 'ctrl_BJ_EH'
        return 'other'
    return {c: f(c) for c in cols}


def arms_65885(cols):
    def f(c):
        s = c.upper()
        if s.startswith('ELR'):                      return 'tr_BJ'
        if s.startswith('EL'):                       return 'prec_BJ'
        if s.startswith('EH'):                       return 'ctrl_BJ_EH'
        if 'ERSRCTAM' in s:
            return 'tr_MCF10A' if '24HR' in s else 'other_early'
        if 'ERSRCETOH' in s:
            return 'prec_MCF10A' if '24HR' in s else 'ctrl_MCF10A_noTAM'
        return 'other'
    return {c: f(c) for c in cols}


SERIES = {
    'GSE118206': dict(load=lambda: load_recount3('GSE118206'), arms=arms_118206, tech='rnaseq',
                      pairs=[('prec_prostate', 'tr_prostate'), ('prec_lung', 'tr_lung'),
                             ('prec_prostate', 'ctrl_dropout'), ('prec_prostate', 'ref_line')],
                      cross=('prec_prostate', 'prec_lung', 'tr_prostate', 'tr_lung')),
    'GSE6885':   dict(load=lambda: load_array('GSE6885', 'GPL570'), arms=arms_6885, tech='array',
                      pairs=[('prec_HME', 'tr_HME'), ('prec_BPE', 'tr_BPE')],
                      cross=('prec_HME', 'prec_BPE', 'tr_HME', 'tr_BPE')),
    'GSE22573':  dict(load=lambda: load_array('GSE22573', 'GPL571'), arms=arms_22573, tech='array',
                      pairs=[('prec_LacZ', 'tr_cdk4RAS'), ('prec_LacZ', 'ctrl_U0126')], cross=None),
    'GSE112369': dict(load=lambda: load_array('GSE112369', 'GPL15207'), arms=arms_112369, tech='array',
                      pairs=[('prec_stomach', 'tr_stomach'), ('prec_stomach', 'ctrl_eng_stomach'),
                             ('prec_stomach', 'tumour_patient'), ('ctrl_eng_colon', 'tr_colon')], cross=None),
    'GSE57965':  dict(load=lambda: load_array('GSE57965', 'GPL15207'), arms=arms_57965, tech='array',
                      pairs=[('adenoma', 'tr_engineered'), ('adenoma', 'tr_mod_adenoma'),
                             ('adenoma', 'tumour_CRC')], cross=None),
    'GSE3151':   dict(load=lambda: load_array('GSE3151', 'GPL570'), arms=arms_3151, tech='array',
                      pairs=[('prec_GFP', f'tr_{k}') for k in ('RAS', 'MYC', 'SRC', 'E2F3', 'BCAT')], cross=None),
    'GSE172029': dict(load=lambda: load_counts('GSE172029', 'counts.csv.gz', ',', 'symbol', 'Ensembl'),
                      arms=arms_172029, tech='rnaseq',
                      pairs=[('prec_WTNC', 'tr_WTNC'), ('prec_WTMB', 'tr_WTMB'), ('prec_WTMEL', 'tr_WTMEL'),
                             ('prec_KONC', 'tr_KONC'), ('prec_KOMB', 'tr_KOMB'), ('prec_KOMEL', 'tr_KOMEL')],
                      cross=None),
    'GSE17941':  dict(load=lambda: load_array('GSE17941', 'GPL9188'), arms=arms_17941, tech='array',
                      pairs=[('prec_MCF10A', 'tr_MCF10A'), ('prec_BJ', 'tr_BJ'),
                             ('prec_MCF10A', 'ctrl_MCF10A_noTAM'), ('prec_BJ', 'ctrl_BJ_EH')],
                      cross=('prec_MCF10A', 'prec_BJ', 'tr_MCF10A', 'tr_BJ'),
                      cross_ctrl=('prec_MCF10A', 'prec_BJ', 'ctrl_MCF10A_noTAM', 'ctrl_BJ_EH')),
    'GSE65885':  dict(load=lambda: load_exp_dir('GSE65885'), arms=arms_65885, tech='rnaseq',
                      pairs=[('prec_MCF10A', 'tr_MCF10A'), ('prec_BJ', 'tr_BJ'),
                             ('prec_BJ', 'ctrl_BJ_EH')], cross=None),
    'GSE231776': dict(load=lambda: load_counts('GSE231776', 'cnt.txt.gz', '\t'), arms=arms_231776, tech='rnaseq',
                      pairs=[('prec_BJ', 'tr_hTERT'), ('prec_BJ', 'tr_SV40'), ('prec_BJ', 'tr_HRAS')], cross=None),
}


# ----------------------------------------------------------------------------- frame
def build_frame():
    panel = L.perceptome_panel()
    hpa = L.hpa_gene_level()
    R_pb = L.score_ref(hpa['H_pb'], hpa['genes'], panel)
    R_sc = L.score_ref(hpa['H_sc'], hpa['genes'], panel)
    es = L.eigenspace(R_pb, n_pc=L.NPC)
    mu, sd, refPn = L.ref_frame(R_sc, es['L'])
    tf = lambda X: X - X.mean(1, keepdims=True)
    mu_s, sd_s, refPn_s = L.ref_frame(tf(np.asarray(R_sc, float)), es['L'])
    fz = json.load(open(FREEZE))
    A = np.array(fz['attractor_raw']); A_s = np.array(fz['attractor_shape'])
    cts = hpa['cts']
    assert cts == fz['ref_cell_types']
    assert abs(float(refPn @ A @ np.ones(1)) if False else 1.0) == 1.0
    chk = float((refPn[cts.index('cytotrophoblasts')] @ A))
    assert abs(chk - fz['ref_cos_to_attractor'][cts.index('cytotrophoblasts')]) < 1e-9, chk
    Z = (np.asarray(R_sc, float) - mu) / sd
    core_idx = [cts.index(c) for c in L.CORE3]
    rest = [i for i in range(len(cts)) if i not in core_idx]
    core_sig = Z[core_idx].mean(0) - Z[rest].mean(0)      # + = module high in the core, - = brake held off
    return dict(panel=panel, cts=cts, R_sc=R_sc, L=es['L'], mu=mu, sd=sd, refPn=refPn, core_sig=core_sig,
                mu_s=mu_s, sd_s=sd_s, refPn_s=refPn_s, A=A, A_s=A_s, regions=fz['regions'])


def score_samples(X, panel, universe=None):
    """X: linear symbol x sample.  -> (n_samples x 44) module scores in the H_sc calibration, plus coverage.

    `universe` (sensitivity calibration C2, declared post hoc after the pre-registered 'all-genes' run gave
    platform-dependent absolute addresses): restrict the matrix to the genes HPA itself measures BEFORE
    scaling to a total of 1e4, so the denominator is comparable across platforms.  The pre-registered
    primary is universe=None."""
    if universe is not None:
        X = X.loc[X.index.intersection(pd.Index(universe))]
    V = X.values.astype(float)
    V = V / (V.sum(0, keepdims=True) + 1e-12) * 1e4
    H = np.log1p(V)                                        # genes x samples
    gi = {g: i for i, g in enumerate(X.index)}
    S = np.zeros((X.shape[1], len(panel)))
    cov = {}
    for j, (m, gl) in enumerate(panel.items()):
        ii = [gi[g] for g in gl if g in gi]
        cov[m] = (len(ii), len(gl))
        if len(ii) >= 2:
            S[:, j] = H[ii].mean(0)
    return S, cov, H, gi


def gene_z(H, gi, genes, cols):
    """mean within-series z of a gene set, per sample."""
    ii = [gi[g] for g in genes if g in gi]
    if not ii:
        return np.full(len(cols), np.nan)
    M = H[ii]
    Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-12)
    return Z.mean(0)


def run_series(gse, F, how='max', universe=None):
    spec = SERIES[gse]
    X = spec['load']()
    arm = pd.Series(spec['arms'](list(X.columns)), index=X.columns)
    S, cov, H, gi = score_samples(X, F['panel'], universe)

    # projections, raw and shape-only
    Pn = L.project(S, F['mu'], F['sd'], F['L'])
    Ss = S - S.mean(1, keepdims=True)
    Pns = L.project(Ss, F['mu_s'], F['sd_s'], F['L'])
    cts = F['cts']; excl = np.array([cts.index(c) for c in L.EXCL])
    top1, _ = L.assign_top1(Pn, F['refPn'], excl)
    top1s, _ = L.assign_top1(Pns, F['refPn_s'], excl)

    per = pd.DataFrame(dict(
        series=gse, sample=list(X.columns), arm=arm.values,
        cos_A=Pn @ F['A'], cos_A_shape=Pns @ F['A_s'],
        top1=[cts[i] for i in top1], top1_shape=[cts[i] for i in top1s],
        magnitude=S.mean(1),
        prolif_z=gene_z(H, gi, PROLIF, X.columns), senesc_z=gene_z(H, gi, SENESC, X.columns),
        cdkn2a_z=gene_z(H, gi, ['CDKN2A'], X.columns),
    ))
    ngenes = sum(1 for m, (a, b) in cov.items() for _ in range(0))
    present = sum(a for a, _ in cov.values()); total = sum(b for _, b in cov.values())
    per.attrs['coverage'] = (present, total)
    return per, S, Pn, Pns, cov, X


def region_hit(labels, members):
    return float(np.mean([lb in members for lb in labels]))


def main(which=None, calib='all'):
    F = build_frame()
    universe = L.hpa_gene_level()['genes'] if calib == 'hpa' else None
    tag = '' if calib == 'all' else f'_{calib}'
    suf = '' if which is None else '_' + '-'.join(sorted(which))   # never overwrite a full run
    reg = F['regions']
    out_rows, pair_rows, cross_rows = [], [], []
    for gse in (which or list(SERIES)):
        print(f'\n=== {gse}')
        per, S, Pn, Pns, cov, X = run_series(gse, F, universe=universe)
        pres, tot = per.attrs['coverage']
        g1 = pres / tot
        print(f'  G1 module-gene coverage {pres}/{tot} = {g1:.1%}  ({X.shape[0]} symbols, {X.shape[1]} samples)')
        thin = [m for m, (a, b) in cov.items() if a < 2]
        if thin:
            print(f'  modules with <2 genes (scored 0): {thin}')
        per['G1_coverage'] = g1
        out_rows.append(per)
        for a, b in [(k, v) for k, v in per.groupby('arm').size().items()]:
            pass
        print('  arms:', dict(per.groupby('arm').size()))

        spec = SERIES[gse]
        for p, t in spec['pairs']:
            dp = per[per.arm == p]; dt = per[per.arm == t]
            if len(dp) == 0 or len(dt) == 0:
                continue
            d_raw = dt.cos_A.mean() - dp.cos_A.mean()
            d_shp = dt.cos_A_shape.mean() - dp.cos_A_shape.mean()
            try:
                u = float(mannwhitneyu(dt.cos_A, dp.cos_A, alternative='greater').pvalue) if min(len(dp), len(dt)) >= 2 else np.nan
            except ValueError:
                u = np.nan
            row = dict(series=gse, precursor=p, transformed=t, n_prec=len(dp), n_tr=len(dt),
                       cosA_prec=dp.cos_A.mean(), cosA_tr=dt.cos_A.mean(), delta_raw=d_raw,
                       cosAs_prec=dp.cos_A_shape.mean(), cosAs_tr=dt.cos_A_shape.mean(), delta_shape=d_shp,
                       mwu_p=u,
                       prolif_delta=dt.prolif_z.mean() - dp.prolif_z.mean(),
                       cdkn2a_delta=dt.cdkn2a_z.mean() - dp.cdkn2a_z.mean(),
                       top1_prec=dp.top1.mode().iat[0] if len(dp) else '', top1_tr=dt.top1.mode().iat[0] if len(dt) else '')
            for rk in ('R1', 'R2', 'R2prime', 'R3'):
                row[f'{rk}_prec'] = region_hit(dp.top1, reg[rk])
                row[f'{rk}_tr'] = region_hit(dt.top1, reg[rk])
                row[f'{rk}_prec_shape'] = region_hit(dp.top1_shape, reg[rk])
                row[f'{rk}_tr_shape'] = region_hit(dt.top1_shape, reg[rk])
            # permutation P on the attractor delta (exchange arm labels within the pair)
            v = np.concatenate([dp.cos_A.values, dt.cos_A.values]); lab = np.r_[np.zeros(len(dp)), np.ones(len(dt))]
            if min(len(dp), len(dt)) >= 1 and len(v) >= 3:
                obs = v[lab == 1].mean() - v[lab == 0].mean()
                tol = 1e-9 * max(1.0, abs(obs))          # ties are lost to summation order without this
                cnt = sum(1 for _ in range(10000)
                          if (lambda q: q[lab == 1].mean() - q[lab == 0].mean() >= obs - tol)(RNG.permutation(v)))
                row['perm_p'] = (cnt + 1) / 10001
            else:
                row['perm_p'] = np.nan
            row['magnitude_delta'] = dt.magnitude.mean() - dp.magnitude.mean()
            # P4: does the module-level move follow the core's own profile?
            mp = S[(per.arm == p).values].mean(0); mt = S[(per.arm == t).values].mean(0)
            row['P4_rho_core_profile'] = float(spearmanr(mt - mp, F['core_sig']).statistic)
            row['G2_OIS'] = bool(row['prolif_delta'] < 0 and row['cdkn2a_delta'] > 0)
            row['G3_prec_inside_R2'] = bool(row['R2_prec'] > 0.5)
            pair_rows.append(row)
            print(f"  {p:18s} -> {t:18s} dcos {d_raw:+.4f} (shape {d_shp:+.4f})  "
                  f"R2 {row['R2_prec']:.2f}->{row['R2_tr']:.2f}  top1 {row['top1_prec']} -> {row['top1_tr']}"
                  f"{'  [OIS]' if row['G2_OIS'] else ''}")

        for ckey in ('cross', 'cross_ctrl'):
          if spec.get(ckey):
            pa, pb, ta, tb = spec[ckey]
            cen = lambda nm: (lambda v: v / (np.linalg.norm(v) + 1e-12))(Pn[(per.arm == nm).values].mean(0))
            cs = lambda u, v: 1 - float(u @ v)
            d0 = cs(cen(pa), cen(pb)); d1 = cs(cen(ta), cen(tb))
            idx_a = np.where(per.arm.isin([pa, ta]).values)[0]
            idx_b = np.where(per.arm.isin([pb, tb]).values)[0]
            lab_a = (per.arm.values[idx_a] == ta); lab_b = (per.arm.values[idx_b] == tb)
            obs = d0 - d1
            nperm, hit = 10000, 0
            for _ in range(nperm):
                la = RNG.permutation(lab_a); lb = RNG.permutation(lab_b)
                ca0 = Pn[idx_a[~la]].mean(0); ca1 = Pn[idx_a[la]].mean(0)
                cb0 = Pn[idx_b[~lb]].mean(0); cb1 = Pn[idx_b[lb]].mean(0)
                n = lambda v: v / (np.linalg.norm(v) + 1e-12)
                if (cs(n(ca0), n(cb0)) - cs(n(ca1), n(cb1))) >= obs - 1e-9:
                    hit += 1
            p = (hit + 1) / (nperm + 1)
            cen_s = lambda nm: (lambda v: v / (np.linalg.norm(v) + 1e-12))(Pns[(per.arm == nm).values].mean(0))
            d0s = cs(cen_s(pa), cen_s(pb)); d1s = cs(cen_s(ta), cen_s(tb)); obs_s = d0s - d1s
            hit_s = 0
            for _ in range(nperm):
                la = RNG.permutation(lab_a); lb = RNG.permutation(lab_b)
                n = lambda v: v / (np.linalg.norm(v) + 1e-12)
                if (cs(n(Pns[idx_a[~la]].mean(0)), n(Pns[idx_b[~lb]].mean(0)))
                        - cs(n(Pns[idx_a[la]].mean(0)), n(Pns[idx_b[lb]].mean(0)))) >= obs_s - 1e-9:
                    hit_s += 1
            ps = (hit_s + 1) / (nperm + 1)
            cross_rows.append(dict(series=gse, kind=ckey, precursors=f'{pa}|{pb}', transformed=f'{ta}|{tb}',
                                   d_prec=d0, d_tr=d1, delta=obs, perm_p=p,
                                   G5_uninformative=bool(d0 < 0.05),
                                   d_prec_shape=d0s, d_tr_shape=d1s, delta_shape=obs_s, perm_p_shape=ps))
            print(f'  P1 shape-only: {d0s:.4f} -> {d1s:.4f}  Delta {obs_s:+.4f}  perm P = {ps:.4f}')
            print(f'  P1 cross-origin: cosdist precursors {d0:.4f} -> transformed {d1:.4f}  '
                  f'Delta {obs:+.4f}  perm P = {p:.4f}')

    per_all = pd.concat(out_rows, ignore_index=True)
    pairs = pd.DataFrame(pair_rows)
    cross = pd.DataFrame(cross_rows)
    per_all.to_csv(os.path.join(RES, f'T9_per_sample{tag}{suf}.csv'), index=False)
    pairs.to_csv(os.path.join(RES, f'T9_pairs{tag}{suf}.csv'), index=False)
    if len(cross):
        cross.to_csv(os.path.join(RES, f'T9_cross_origin{tag}{suf}.csv'), index=False)
    print(f'\nwrote T9_per_sample.csv ({len(per_all)}), T9_pairs.csv ({len(pairs)}), T9_cross_origin.csv ({len(cross)})')


if __name__ == '__main__':
    av = [a for a in sys.argv[1:] if not a.startswith('--')]
    cal = 'hpa' if '--hpa' in sys.argv else 'all'
    main(av or None, cal)
