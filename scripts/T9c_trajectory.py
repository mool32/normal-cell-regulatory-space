#!/usr/bin/env python3
"""T9c — Addendum A2: does the coordinate read TRANSFORMATION or the CONSTRUCT?

3a GSE240057 — PARCB time course, 10 donors, paired within donor:
      basal (0wks) -> PARCB-transduced organoid (0wks, all five factors, no tumour yet)
      -> early -> transition -> endpoint tumour.
   f = delta_cocktail / delta_endpoint, thresholded in advance (>=0.7 cocktail, <=0.3 progression).
3b GSE148714 — same MLL fusion into cord-blood vs bone-marrow CD34+ HSPC; only cord blood transforms.

Usage: python3 -I scripts/T9c_trajectory.py
"""
import os, sys, json, gzip, re
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, spearmanr

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'perceptome_validation', 'scripts'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import prereg_lib as L  # noqa: E402
import T9_invitro_transformation as T9  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'invitro')
RES = os.path.join(ROOT, 'perceptome_validation', 'results')
STAGE_IDX = dict(prec_basal=0, ctrl_transduced=1, tr_early=2, tr_transition=3, tr_endpoint=4)


def load_240057(universe):
    """symbol x sample counts; R make.unique '.N' suffixes stripped; length-normalised."""
    X = pd.read_csv(os.path.join(DATA, 'GSE240057', 'counts.txt.gz'), sep='\t', index_col=0)
    # R make.unique appends '.N' to LATER duplicates, so the bare base is always present too.
    # Clone-based GENCODE symbols (AC000003.1, RP11-432M8.1) also end in '.N' but have no bare base:
    # stripping those destroys ~23,000 real symbols and perturbs the per-sample denominator.  Strip only
    # when the bare base exists in this same file.
    raw = pd.Index([str(i) for i in X.index])
    present = set(raw)
    X.index = pd.Index([re.sub(r'\.\d+$', '', n) if (re.search(r'\.\d+$', n)
                                                      and re.sub(r'\.\d+$', '', n) in present) else n
                        for n in raw]).str.upper()
    _, by_sym = T9.gencode()
    ln = np.array([by_sym.get(s, np.nan) for s in X.index], float)
    X = X.loc[~np.isnan(ln)].div(ln[~np.isnan(ln)], axis=0)
    X = T9.collapse(X, 'max')
    # arms: per donor, TP1 = basal, TP2 = transduced 0wks, last = endpoint, TP3 = early, rest = transition
    don = pd.Series([c.split('-')[0] for c in X.columns], index=X.columns)
    tp = pd.Series([int(c.split('-TP')[1]) for c in X.columns], index=X.columns)
    last = tp.groupby(don).max()
    arm = {}
    for c in X.columns:
        d, t = don[c], tp[c]
        arm[c] = ('prec_basal' if t == 1 else 'ctrl_transduced' if t == 2
                  else 'tr_endpoint' if t == last[d] else 'tr_early' if t == 3 else 'tr_transition')
    return X, pd.Series(arm), don


def load_148714(universe):
    X = pd.read_csv(os.path.join(DATA, 'GSE148714', 'rawCounts.csv.gz'), index_col=0)
    by_ens, _ = T9.gencode()
    base = pd.Index([str(g).split('.')[0] for g in X.index])
    ln = np.array([by_ens.get(g, np.nan) for g in base], float)
    keep = ~np.isnan(ln)
    X = X.loc[keep].div(ln[keep], axis=0)
    X.index = pd.Index([T9._ens_symbol(g) for g in base[keep]])
    X = T9.collapse(X, 'max')
    title = {'I19R031a01': 'CB_MLL-AF4', 'I19R031a02': 'CB_MLL-AF4', 'I19R031a03': 'CB_MLL-AF9',
             'I19R031a04': 'CB_MLL-AF9', 'I19R031a05': 'BM_MLL-AF4', 'I19R031a13': 'Ctrl',
             'I19R031a14': 'Ctrl', 'I19R031a17': 'BM_MLL-AF9', 'I19R031a23': 'BM_MLL-AF9',
             'I19R031a36': 'BM_MLL-AF4'}
    arm = pd.Series({c: title.get(c, 'other') for c in X.columns})
    assert 'other' not in set(arm), set(arm)
    return X, arm


def project(X, F, A, A_s, universe):
    S, cov, _, _ = T9.score_samples(X, F['panel'], universe)
    Pn = L.project(S, F['mu'], F['sd'], F['L'])
    Pns = L.project(S - S.mean(1, keepdims=True), F['mu_s'], F['sd_s'], F['L'])
    cts = F['cts']; excl = np.array([cts.index(c) for c in L.EXCL])
    t1, _ = L.assign_top1(Pn, F['refPn'], excl)
    pres = sum(a for a, _ in cov.values()); tot = sum(b for _, b in cov.values())
    return pd.DataFrame(dict(sample=list(X.columns), cos_A=Pn @ A, cos_A_shape=Pns @ A_s,
                             top1=[cts[i] for i in t1])), pres / tot


def main():
    F = T9.build_frame()
    fz = json.load(open(os.path.join(RES, 'T9_attractor.json')))
    A = np.array(fz['attractor_raw']); A_s = np.array(fz['attractor_shape'])
    R2 = fz['regions']['R2']
    out = {}

    for calib in ('all', 'hpa'):
        universe = L.hpa_gene_level()['genes'] if calib == 'hpa' else None
        tag = '' if calib == 'all' else '_hpa'
        print(f'\n################ calibration: {calib}')

        # ---------------- 3a GSE240057 ----------------
        X, arm, don = load_240057(universe)
        per, cov = project(X, F, A, A_s, universe)
        per['arm'] = arm.values; per['donor'] = don.values
        per['stage'] = per.arm.map(STAGE_IDX)
        print(f'\n=== GSE240057  coverage {cov:.1%}  {X.shape[1]} samples, {per.donor.nunique()} donors')
        print(per.groupby('arm')[['cos_A', 'cos_A_shape']].agg(['mean', 'size']).round(3).to_string())
        print('  top-1 by stage:', {a: per[per.arm == a].top1.mode().iat[0] for a in STAGE_IDX})
        print('  in region R2 by stage:',
              {a: round(float(per[per.arm == a].top1.isin(R2).mean()), 2) for a in STAGE_IDX})

        w = per.pivot_table(index='donor', columns='arm', values='cos_A')
        ws = per.pivot_table(index='donor', columns='arm', values='cos_A_shape')
        d_end = (w.tr_endpoint - w.prec_basal).dropna()
        d_ckt = (w.ctrl_transduced - w.prec_basal).dropna()
        d_end_s = (ws.tr_endpoint - ws.prec_basal).dropna()
        d_ckt_s = (ws.ctrl_transduced - ws.prec_basal).dropna()
        f_raw = float(d_ckt.mean() / d_end.mean()); f_shape = float(d_ckt_s.mean() / d_end_s.mean())
        verdict = ('COCKTAIL' if f_raw >= 0.7 else 'PROGRESSION' if f_raw <= 0.3 else 'AMBIGUOUS')
        print(f'  delta endpoint   raw {d_end.mean():+.4f} (n={len(d_end)}, {int((d_end>0).sum())} up, '
              f'Wilcoxon P={wilcoxon(d_end).pvalue:.4g})   shape {d_end_s.mean():+.4f} '
              f'({int((d_end_s>0).sum())} up, P={wilcoxon(d_end_s).pvalue:.4g})')
        print(f'  delta cocktail   raw {d_ckt.mean():+.4f} ({int((d_ckt>0).sum())} up)   '
              f'shape {d_ckt_s.mean():+.4f} ({int((d_ckt_s>0).sum())} up)')
        print(f'  >>> f = {f_raw:.3f} raw / {f_shape:.3f} shape   -> {verdict}')
        rhos = [spearmanr(g.stage, g.cos_A).statistic for _, g in per.groupby('donor') if g.stage.nunique() > 2]
        print(f'  stage ordering: mean rho {np.mean(rhos):+.3f}, positive in {sum(r>0 for r in rhos)}/{len(rhos)} donors')
        per.to_csv(os.path.join(RES, f'T9c_GSE240057_per_sample{tag}.csv'), index=False)
        out[f'GSE240057{tag}'] = dict(coverage=cov, f_raw=f_raw, f_shape=f_shape, verdict=verdict,
                                      d_end=float(d_end.mean()), d_ckt=float(d_ckt.mean()),
                                      d_end_shape=float(d_end_s.mean()), d_ckt_shape=float(d_ckt_s.mean()),
                                      n_donor=int(len(d_end)), up=int((d_end > 0).sum()),
                                      wilcoxon_p=float(wilcoxon(d_end).pvalue),
                                      wilcoxon_p_shape=float(wilcoxon(d_end_s).pvalue),
                                      mean_rho=float(np.mean(rhos)),
                                      rho_pos=int(sum(r > 0 for r in rhos)), n_rho=len(rhos),
                                      region_by_stage={a: float(per[per.arm == a].top1.isin(R2).mean())
                                                       for a in STAGE_IDX})

        # ---------------- 3b GSE148714 ----------------
        Y, yarm = load_148714(universe)
        yper, ycov = project(Y, F, A, A_s, universe)
        yper['arm'] = yarm.reindex(yper['sample']).values
        base = yper[yper.arm == 'Ctrl'].cos_A.mean(); base_s = yper[yper.arm == 'Ctrl'].cos_A_shape.mean()
        print(f'\n=== GSE148714  coverage {ycov:.1%}')
        tb = yper.groupby('arm')[['cos_A', 'cos_A_shape']].mean()
        tb['delta_raw'] = tb.cos_A - base; tb['delta_shape'] = tb.cos_A_shape - base_s
        tb['top1'] = yper.groupby('arm').top1.agg(lambda s: s.mode().iat[0])
        print(tb.round(4).to_string())
        d_cb = float(tb.loc['CB_MLL-AF4', 'delta_raw']); d_bm = float(tb.loc['BM_MLL-AF4', 'delta_raw'])
        print(f'  >>> MLL-AF4: cord blood (transforms) {d_cb:+.4f}  vs  bone marrow (does not) {d_bm:+.4f}'
              f'  -> {"SPECIFIC" if d_cb > d_bm else "NOT SPECIFIC"}')
        yper.to_csv(os.path.join(RES, f'T9c_GSE148714_per_sample{tag}.csv'), index=False)
        out[f'GSE148714{tag}'] = dict(coverage=ycov, cb_af4=d_cb, bm_af4=d_bm,
                                      cb_af9=float(tb.loc['CB_MLL-AF9', 'delta_raw']),
                                      bm_af9=float(tb.loc['BM_MLL-AF9', 'delta_raw']),
                                      specific=bool(d_cb > d_bm),
                                      top1=dict(tb['top1']))

    json.dump(out, open(os.path.join(RES, 'T9c_summary.json'), 'w'), indent=1)
    print('\nwrote T9c_summary.json')


if __name__ == '__main__':
    main()
