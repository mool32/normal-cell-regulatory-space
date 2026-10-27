#!/usr/bin/env python3
"""E26b — is the D4 persister trough in Fig. S26 a transcriptome-complexity artifact of our own QC floor?

E26 keeps cells with >= 500 detected genes.  That retains 42% of D0 but only 9.5% of D4 -- and D4 is the
trough that makes the trajectory non-monotone.  Drug-tolerant persisters are transcriptionally quiescent,
so a complexity floor removes preferentially the cells the timepoint exists to show.

Three checks, none of which needs new data:
  (1) sweep the floor (0, 200, 300, 500, 1000) and ask whether the trough survives;
  (2) within each timepoint, core fraction by detected-gene decile -- the Fig. 3D control applied here;
  (3) complexity-matched resampling: resample every timepoint to D4's detected-gene distribution.

Writes results/referee2_newdata/E26b_qc_sweep.csv, E26b_decile.csv, E26b_matched.csv
"""
import os, re, sys, glob, json, gc, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
warnings.filterwarnings('ignore')
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'tool/perceptome2'))
import perceptome as pct  # noqa: E402
from perceptome.perceptivity.reference import load_hpa_perceptivity  # noqa: E402

OUT = os.path.join(ROOT, 'results', 'referee2_newdata')
CORE = ['megakaryocytes', 'cytotrophoblasts', 'migrating cytotrophoblasts']
FLOORS = [0, 200, 300, 500, 1000]
RNG = np.random.default_rng(0)

e = json.load(open(os.path.join(ROOT, 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')))
mods = e['module_order']
load = np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm = load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
refcells = list(Rm.index)
core_idx = np.array([refcells.index(c) for c in CORE])
EXCL = np.array([refcells.index(c) for c in ['erythrocytes', 'platelets'] if c in refcells])
demag = lambda X: X - X.mean(1, keepdims=True)
TF = dict(RAW=lambda X: X, SHAPE=demag)


def refset(tf):
    R = tf(Rm.values.astype(float)); mu = R.mean(0); sd = R.std(0) + 1e-12
    rp = ((R - mu) / sd) @ load
    return mu, sd, rp / (np.linalg.norm(rp, axis=1, keepdims=True) + 1e-12)


REF = {m: refset(TF[m]) for m in TF}


def core_flags(S):
    """-> dict metric -> boolean array, is this cell's top-1 in the three-cell core."""
    out = {}
    for m, (mu, sd, rpn) in REF.items():
        pc = ((TF[m](S) - mu) / sd) @ load
        pcn = pc / (np.linalg.norm(pc, axis=1, keepdims=True) + 1e-12)
        sim = pcn @ rpn.T; sim[:, EXCL] = -1e9
        out[m] = np.isin(np.argmax(sim, 1), core_idx)
    return out


def load_day(f):
    df = pd.read_csv(f, sep='\t', index_col=0)
    df.index = [str(g).upper() for g in df.index]
    if pd.Index(df.index).has_duplicates:
        df = df.groupby(level=0).sum()
    a = ad.AnnData(df.T.values.astype(np.float32))
    a.var_names = list(df.index); a.var_names_make_unique()
    del df
    ng = np.asarray((a.X > 0).sum(1)).ravel()
    if a.X.max() > 100:
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    S = pct.score_modules(a, method='mean_raw', gene_set='core')['scores'] \
        .reindex(columns=mods, fill_value=0.0).values.astype(float)
    del a; gc.collect()
    return S, ng


def main():
    days, Ss, NGs = [], {}, {}
    for f in sorted(glob.glob(os.path.join(ROOT, 'paper6/data/aissa2021/extracted/*.dge.txt.gz'))):
        d = int(re.search(r'_D(\d+)\.', os.path.basename(f)).group(1))
        S, ng = load_day(f)
        days.append(d); Ss[d] = S; NGs[d] = ng
        print(f'  D{d}: {len(ng)} cells, detected genes median {np.median(ng):.0f} '
              f'[{np.percentile(ng,10):.0f}-{np.percentile(ng,90):.0f}]', flush=True)
    days.sort()

    flags = {d: core_flags(Ss[d]) for d in days}

    # ---- (1) floor sweep -------------------------------------------------------------------
    rows = []
    for fl in FLOORS:
        for d in days:
            k = NGs[d] >= fl
            if k.sum() < 10:
                continue
            rows.append(dict(floor=fl, day=d, n=int(k.sum()), kept=float(k.mean()),
                             RAW=100 * flags[d]['RAW'][k].mean(),
                             SHAPE=100 * flags[d]['SHAPE'][k].mean()))
    sw = pd.DataFrame(rows)
    sw.to_csv(os.path.join(OUT, 'E26b_qc_sweep.csv'), index=False)
    print('\n=== (1) core fraction (SHAPE) by QC floor')
    piv = sw.pivot(index='floor', columns='day', values='SHAPE').round(1)
    npiv = sw.pivot(index='floor', columns='day', values='n')
    print(piv.to_string())
    print('n kept:'); print(npiv.to_string())
    for fl in FLOORS:
        r = sw[sw.floor == fl].sort_values('day')
        if len(r) < 4:
            continue
        v = r.SHAPE.values; d4 = r[r.day == 4].SHAPE.values
        trough = bool(len(d4) and d4[0] == v.min())
        print(f'  floor {fl:4d}: D4 = {d4[0] if len(d4) else float("nan"):.1f}%  '
              f'min of the course = {v.min():.1f}%  -> D4 is the trough: {trough}')

    # ---- (2) core fraction by detected-gene decile, within each timepoint ------------------
    rows = []
    for d in days:
        ng = NGs[d]
        q = np.quantile(ng, np.linspace(0, 1, 11))
        for i in range(10):
            k = (ng >= q[i]) & (ng <= q[i + 1])
            if k.sum() < 5:
                continue
            rows.append(dict(day=d, decile=i + 1, n=int(k.sum()), median_genes=float(np.median(ng[k])),
                             RAW=100 * flags[d]['RAW'][k].mean(), SHAPE=100 * flags[d]['SHAPE'][k].mean()))
    de = pd.DataFrame(rows)
    de.to_csv(os.path.join(OUT, 'E26b_decile.csv'), index=False)
    print('\n=== (2) core fraction (SHAPE) by detected-gene decile, within timepoint')
    print(de.pivot(index='decile', columns='day', values='SHAPE').round(1).to_string())
    from scipy.stats import spearmanr
    for d in days:
        s = de[de.day == d]
        print(f'  D{d}: Spearman(decile, core%) = {spearmanr(s.decile, s.SHAPE).statistic:+.2f}')

    # ---- (2b) common ABSOLUTE complexity window: the fair like-with-like comparison --------
    rows = []
    for lo, hi in [(500, 1500), (1000, 3000), (1500, 6000), (2000, 8000)]:
        for d in days:
            k = (NGs[d] >= lo) & (NGs[d] < hi)
            rows.append(dict(lo=lo, hi=hi, day=d, n=int(k.sum()),
                             RAW=100 * flags[d]['RAW'][k].mean() if k.sum() >= 15 else np.nan,
                             SHAPE=100 * flags[d]['SHAPE'][k].mean() if k.sum() >= 15 else np.nan))
    win = pd.DataFrame(rows)
    win.to_csv(os.path.join(OUT, 'E26b_window.csv'), index=False)
    print('\n=== (2b) core fraction (SHAPE) inside a common ABSOLUTE detected-gene window')
    print(win.pivot(index=['lo', 'hi'], columns='day', values='SHAPE').round(1).to_string())
    print('n per cell:')
    print(win.pivot(index=['lo', 'hi'], columns='day', values='n').to_string())

    # ---- (3) complexity-matched resampling to D4's distribution ----------------------------
    ref_ng = NGs[4]
    q = np.quantile(ref_ng, np.linspace(0, 1, 11))
    want = np.histogram(ref_ng, bins=q)[0]
    rows = []
    for d in days:
        ng = NGs[d]
        reps = []
        for _ in range(200):
            pick = []
            for i in range(10):
                pool = np.where((ng >= q[i]) & (ng <= q[i + 1]))[0]
                if len(pool) == 0 or want[i] == 0:
                    continue
                pick += list(RNG.choice(pool, size=int(want[i]), replace=True))
            if pick:
                reps.append(100 * flags[d]['SHAPE'][pick].mean())
        rows.append(dict(day=d, matched_mean=float(np.mean(reps)), matched_sd=float(np.std(reps)),
                         n_draws=len(reps),
                         unmatched=100 * flags[d]['SHAPE'][NGs[d] >= 500].mean() if (NGs[d] >= 500).sum() else np.nan))
    mt = pd.DataFrame(rows)
    mt.to_csv(os.path.join(OUT, 'E26b_matched.csv'), index=False)
    print('\n=== (3) every timepoint resampled to D4\'s detected-gene distribution (SHAPE core %)')
    print(mt.round(2).to_string(index=False))
    print('\nwrote E26b_qc_sweep.csv, E26b_decile.csv, E26b_matched.csv')


if __name__ == '__main__':
    main()
