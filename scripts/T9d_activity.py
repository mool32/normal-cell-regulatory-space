#!/usr/bin/env python3
"""T9d — Addendum A3: activity scoring on the GSE240057 PARCB time course.

Does the week-0 jump (f = 0.699 under readiness scoring) survive when every module is scored by regulatory
OUTPUT instead of component abundance?  Plus control C3: drop the three delivered genes that sit in modules.

Reuses T7's activity definitions verbatim and T7's cached reference/cancer activity matrices.

Usage: python3 -I scripts/T9d_activity.py [viper|ulm]
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy.stats import wilcoxon, spearmanr

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'perceptome_validation', 'scripts'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import prereg_lib as L  # noqa: E402
import T9_invitro_transformation as T9  # noqa: E402
import T9c_trajectory as T9c  # noqa: E402
import decoupler as dc  # noqa: E402

EST = (sys.argv[1] if len(sys.argv) > 1 else 'viper').lower()
RES = os.path.join(ROOT, 'perceptome_validation', 'results')
PAN = os.path.join(ROOT, 'perceptome_validation', 'panels')
STAGE = T9c.STAGE_IDX
DELIVERED = ['TP53', 'RB1', 'AKT1']          # PARCB genes that appear in the catalog (MYC, BCL2 do not)


# ------------------------------------------------------------------ T7's activity definitions, verbatim
def activity_net(mods):
    import perceptome as pct
    dor = pd.read_csv(os.path.join(PAN, 'dorothea_AB.csv')); dor['target'] = dor.target.str.upper()
    prog = pd.read_csv(os.path.join(PAN, 'progeny_top100.csv')); prog['target'] = prog.target.str.upper()
    PROGENY_MAP = {'ERK/MAPK': 'MAPK', 'PI3K/PTEN': 'PI3K'}
    rows, mean_sets, how = [], {}, {}
    for m in mods:
        info = pct.get_module_info(m); tfs = info.get('primary_tf') or []
        reg = dor[dor.source.isin(tfs)].drop_duplicates('target') if tfs else dor.iloc[0:0]
        if tfs and len(reg) >= 5:
            rows += [dict(source=m, target=r.target, weight=float(r.weight)) for _, r in reg.iterrows()]
            how[m] = f'dorothea({len(reg)})'
        elif (not tfs) and m in PROGENY_MAP:
            pr = prog[prog.source == PROGENY_MAP[m]]
            rows += [dict(source=m, target=r.target, weight=float(r.weight)) for _, r in pr.iterrows()]
            how[m] = f'progeny:{PROGENY_MAP[m]}'
        else:
            mean_sets[m] = [str(g).upper() for g in (info.get('activity_genes') or [])]
            how[m] = f'mean_activity({len(mean_sets[m])})'
    return pd.DataFrame(rows), mean_sets, how


def estimate(X, var_names, net, mean_sets, mods):
    X = np.asarray(X, np.float32)
    a = ad.AnnData(X=sp.csr_matrix(X))
    a.obs_names = [str(i) for i in range(X.shape[0])]
    a.var_names = list(var_names)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        (dc.mt.viper if EST == 'viper' else dc.mt.ulm)(a, net, tmin=3, verbose=False)
    sc_ = a.obsm['score_viper' if EST == 'viper' else 'score_ulm']
    out = np.zeros((X.shape[0], len(mods)), np.float32)
    gi = {g: i for i, g in enumerate(var_names)}
    for j, m in enumerate(mods):
        if m in sc_.columns:
            out[:, j] = sc_[m].values
        elif m in mean_sets:
            ii = [gi[g] for g in mean_sets[m] if g in gi]
            if len(ii) >= 2:
                out[:, j] = X[:, ii].mean(1)
    return out


# ------------------------------------------------------------------ paired statistics
def paired(per, col):
    w = per.pivot_table(index='donor', columns='arm', values=col)
    return (w.tr_endpoint - w.prec_basal).dropna(), (w.ctrl_transduced - w.prec_basal).dropna()


def report(per, R2, label):
    d_end, d_ckt = paired(per, 'cos_A')
    d_end_s, d_ckt_s = paired(per, 'cos_A_shape')
    f = float(d_ckt.mean() / d_end.mean())
    rhos = [spearmanr(g.stage, g.cos_A).statistic for _, g in per.groupby('donor') if g.stage.nunique() > 2]
    lad = {a: f"{int(per[(per.arm == a) & per.top1.isin(R2)].shape[0])}/{int((per.arm == a).sum())}"
           for a in STAGE}
    print(f'\n--- {label}')
    print(per.groupby('arm').cos_A.mean().reindex(list(STAGE)).round(3).to_string())
    print(f'  delta endpoint raw {d_end.mean():+.4f} ({int((d_end>0).sum())}/{len(d_end)} up, '
          f'P={wilcoxon(d_end).pvalue:.4g})   shape {d_end_s.mean():+.4f} ({int((d_end_s>0).sum())}/{len(d_end_s)} up, '
          f'P={wilcoxon(d_end_s).pvalue:.4g})')
    print(f'  delta cocktail raw {d_ckt.mean():+.4f} ({int((d_ckt>0).sum())}/{len(d_ckt)} up)   '
          f'shape {d_ckt_s.mean():+.4f}')
    verdict = 'ABUNDANCE ARTEFACT' if f <= 0.3 else 'REAL AT ACTIVITY LEVEL' if f >= 0.7 else 'PARTIAL'
    print(f'  >>> f = {f:.3f}   (readiness was 0.699)  -> {verdict}')
    print(f'  stage ordering mean rho {np.mean(rhos):+.3f}, positive {sum(r>0 for r in rhos)}/{len(rhos)}')
    print(f'  in region R2 by stage: {lad}')
    print(f'  top-1 by stage: { {a: per[per.arm==a].top1.mode().iat[0] for a in STAGE} }')
    return dict(f=f, verdict=verdict, d_end=float(d_end.mean()), d_ckt=float(d_ckt.mean()),
                d_end_shape=float(d_end_s.mean()), d_ckt_shape=float(d_ckt_s.mean()),
                up=int((d_end > 0).sum()), n=int(len(d_end)),
                p=float(wilcoxon(d_end).pvalue), p_shape=float(wilcoxon(d_end_s).pvalue),
                shape_up=int((d_end_s > 0).sum()),
                mean_rho=float(np.mean(rhos)), rho_pos=int(sum(r > 0 for r in rhos)), ladder=lad,
                top1={a: per[per.arm == a].top1.mode().iat[0] for a in STAGE})


def main():
    out = {}
    F = T9.build_frame()
    mods = list(F['panel']); cts = F['cts']
    excl = np.array([cts.index(c) for c in L.EXCL])
    fz = json.load(open(os.path.join(RES, 'T9_attractor.json')))
    R2 = fz['regions']['R2']
    X, arm, don = T9c.load_240057(None)

    # ================= control C3: readiness scoring, delivered genes removed ==================
    print('================ C3: readiness scoring, delivered genes removed', DELIVERED)
    panel_c3 = {m: [g for g in gl if g not in DELIVERED] for m, gl in F['panel'].items()}
    hpa = L.hpa_gene_level()
    for nm, pan in (('readiness (all genes, reproduction)', F['panel']), ('readiness minus delivered', panel_c3)):
        R_pb = L.score_ref(hpa['H_pb'], hpa['genes'], pan)
        R_sc = L.score_ref(hpa['H_sc'], hpa['genes'], pan)
        es = L.eigenspace(R_pb, n_pc=L.NPC)
        mu, sd, refPn = L.ref_frame(R_sc, es['L'])
        tf = lambda M: M - M.mean(1, keepdims=True)
        mu_s, sd_s, refPn_s = L.ref_frame(tf(np.asarray(R_sc, float)), es['L'])
        cache = L.load_panel_cache()
        Sc = L.score_cells(cache, pan)
        Pc = L.project(Sc, mu, sd, es['L'])
        cid = cache['obs']['cancer'].astype(str).values
        A = np.mean([(lambda v: v / np.linalg.norm(v))(Pc[cid == c].mean(0)) for c in sorted(set(cid))], 0)
        A = A / np.linalg.norm(A)
        Pcs = L.project(tf(Sc.astype(float)), mu_s, sd_s, es['L'])
        A_s = np.mean([(lambda v: v / np.linalg.norm(v))(Pcs[cid == c].mean(0)) for c in sorted(set(cid))], 0)
        A_s = A_s / np.linalg.norm(A_s)
        del cache, Sc, Pc, Pcs
        S, _, _, _ = T9.score_samples(X, pan, None)
        Pn = L.project(S, mu, sd, es['L']); Pns = L.project(S - S.mean(1, keepdims=True), mu_s, sd_s, es['L'])
        t1, _ = L.assign_top1(Pn, refPn, excl)
        per = pd.DataFrame(dict(cos_A=Pn @ A, cos_A_shape=Pns @ A_s, top1=[cts[i] for i in t1],
                                arm=arm.values, donor=don.values))
        per['stage'] = per.arm.map(STAGE)
        out[nm] = report(per, R2, nm)

    # ================= activity scoring ==================
    print(f'\n================ activity scoring ({EST})')
    net, mean_sets, how = activity_net(mods)
    npz = np.load(os.path.join(RES, f'T7_activity_scores_{EST}.npz'), allow_pickle=True)
    assert list(npz['mods']) == mods
    A_pb, A_sc, A_cells = npz['A_pb'], npz['A_sc'], npz['A_cells']
    es = L.eigenspace(A_pb, n_pc=L.NPC)
    mu, sd, refPn = L.ref_frame(A_sc, es['L'])
    tf = lambda M: M - M.mean(1, keepdims=True)
    mu_s, sd_s, refPn_s = L.ref_frame(tf(np.asarray(A_sc, float)), es['L'])

    cache_obs = L.load_panel_cache()['obs']
    cid = cache_obs['cancer'].astype(str).values
    Pc = L.project(A_cells, mu, sd, es['L'])
    Aat = np.mean([(lambda v: v / np.linalg.norm(v))(Pc[cid == c].mean(0)) for c in sorted(set(cid))], 0)
    Aat = Aat / np.linalg.norm(Aat)
    Pcs = L.project(tf(A_cells.astype(float)), mu_s, sd_s, es['L'])
    Aat_s = np.mean([(lambda v: v / np.linalg.norm(v))(Pcs[cid == c].mean(0)) for c in sorted(set(cid))], 0)
    Aat_s = Aat_s / np.linalg.norm(Aat_s)
    json.dump(dict(estimator=EST, attractor=Aat.tolist(), attractor_shape=Aat_s.tolist(),
                   ref_cos=(refPn @ Aat).tolist(), cell_types=cts, how=how),
              open(os.path.join(RES, f'T9d_attractor_{EST}.json'), 'w'))
    order = np.argsort(-(refPn @ Aat))
    print('  activity-space top 8 reference cells by cosine to the attractor:',
          [(cts[i], round(float((refPn @ Aat)[i]), 3)) for i in order[:8]])
    del Pc, Pcs, A_cells, cache_obs

    # GSE240057 at activity level, same H_sc calibration
    V = X.values.astype(float)
    H = np.log1p(V / (V.sum(0, keepdims=True) + 1e-12) * 1e4).T          # samples x genes
    Aex = estimate(H, list(X.index), net, mean_sets, mods)
    Pn = L.project(Aex, mu, sd, es['L'])
    Pns = L.project(Aex - Aex.mean(1, keepdims=True), mu_s, sd_s, es['L'])
    t1, _ = L.assign_top1(Pn, refPn, excl)
    per = pd.DataFrame(dict(cos_A=Pn @ Aat, cos_A_shape=Pns @ Aat_s, top1=[cts[i] for i in t1],
                            arm=arm.values, donor=don.values))
    per['stage'] = per.arm.map(STAGE)
    per.to_csv(os.path.join(RES, f'T9d_GSE240057_activity_{EST}.csv'), index=False)
    out[f'activity_{EST}'] = report(per, R2, f'activity ({EST})')

    json.dump(out, open(os.path.join(RES, f'T9d_summary_{EST}.json'), 'w'), indent=1)
    print(f'\nwrote T9d_summary_{EST}.json')


if __name__ == '__main__':
    main()
