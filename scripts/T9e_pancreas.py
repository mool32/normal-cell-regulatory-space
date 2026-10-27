#!/usr/bin/env python3
"""T9e — Addendum A4: GSE222990, two primary human pancreatic origins, same four drivers,
with a non-transforming driver arm (KRAS alone) and a culture-only arm.

Readouts (prereg A4): P1 cross-origin; the acinar ladder with culture held constant (k', g);
all under readiness scoring and, as declared in advance, under activity scoring too.

Usage: python3 -I scripts/T9e_pancreas.py
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'perceptome_validation', 'scripts'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import prereg_lib as L  # noqa: E402
import T9_invitro_transformation as T9  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'invitro', 'GSE222990')
RES = os.path.join(ROOT, 'perceptome_validation', 'results')
RNG = np.random.default_rng(11)
ARMS = ['fresh_acinar', 'fresh_ductal', 'acinar_organoid', 'acinar_K_organoid',
        'acinar_KPTS_organoid', 'ductal_KPTS_cultured']


def load():
    X = pd.read_csv(os.path.join(DATA, 'TPM.txt.gz'), sep='\t', index_col=0)
    X.index = pd.Index([str(i).upper() for i in X.index])
    X = T9.collapse(X, 'max')                      # TPM already length-normalised
    arm = {}
    for c in X.columns:
        hit = [a for a in ARMS if c.startswith(a)]
        assert len(hit) == 1, (c, hit)
        arm[c] = hit[0]
    return X, pd.Series(arm)


def cross_origin(P, arm, pa, pb, ta, tb, nperm=10000):
    cen = lambda nm: (lambda v: v / (np.linalg.norm(v) + 1e-12))(P[(arm == nm).values].mean(0))
    cs = lambda u, v: 1 - float(u @ v)
    d0, d1 = cs(cen(pa), cen(pb)), cs(cen(ta), cen(tb))
    obs = d0 - d1
    ia = np.where(arm.isin([pa, ta]).values)[0]; ib = np.where(arm.isin([pb, tb]).values)[0]
    la = (arm.values[ia] == ta); lb = (arm.values[ib] == tb)
    n = lambda v: v / (np.linalg.norm(v) + 1e-12)
    hit = 0
    for _ in range(nperm):
        a, b = RNG.permutation(la), RNG.permutation(lb)
        if (cs(n(P[ia[~a]].mean(0)), n(P[ib[~b]].mean(0)))
                - cs(n(P[ia[a]].mean(0)), n(P[ib[b]].mean(0)))) >= obs:
            hit += 1
    return dict(d_prec=d0, d_tr=d1, delta=obs, perm_p=(hit + 1) / (nperm + 1),
                G5_uninformative=bool(d0 < 0.05))


def analyse(per, R2, label, out):
    m = per.groupby('arm').cos_A.mean(); ms = per.groupby('arm').cos_A_shape.mean()
    print(f'\n--- {label}')
    print(m.reindex(ARMS).round(4).to_string())
    print('  shape-only:', ms.reindex(ARMS).round(4).to_dict())
    print('  top-1:', {a: per[per.arm == a].top1.mode().iat[0] for a in ARMS})
    print('  in region R2:', {a: f"{int(per[(per.arm == a) & per.top1.isin(R2)].shape[0])}/"
                                 f"{int((per.arm == a).sum())}" for a in ARMS})
    # acinar ladder, culture held constant
    dK = float(m.acinar_KPTS_organoid - m.acinar_organoid)
    dk = float(m.acinar_K_organoid - m.acinar_organoid)
    dKs = float(ms.acinar_KPTS_organoid - ms.acinar_organoid)
    dks = float(ms.acinar_K_organoid - ms.acinar_organoid)
    g = float((m.acinar_organoid - m.fresh_acinar) / (m.acinar_KPTS_organoid - m.fresh_acinar))
    kp = dk / dK if dK != 0 else float('nan')
    kps = dks / dKs if dKs != 0 else float('nan')
    u = lambda a, b: float(mannwhitneyu(per[per.arm == a].cos_A, per[per.arm == b].cos_A,
                                        alternative='greater').pvalue)
    print(f"  delta'(KPTS vs organoid)  raw {dK:+.4f} (P={u('acinar_KPTS_organoid','acinar_organoid'):.4f})"
          f"   shape {dKs:+.4f}")
    print(f"  delta'(K    vs organoid)  raw {dk:+.4f} (P={u('acinar_K_organoid','acinar_organoid'):.4f})"
          f"   shape {dks:+.4f}")
    print(f"  >>> k' = {kp:.3f} raw / {kps:.3f} shape   (>=0.7 KRAS-alone explains it, <=0.3 it does not)")
    print(f"  >>> g  = {g:.3f}  (share of the fresh->KPTS move delivered by culture alone)")
    print(f"  ductal (PLATFORM-CONFOUNDED, G6): raw {float(m.ductal_KPTS_cultured - m.fresh_ductal):+.4f}"
          f"   shape {float(ms.ductal_KPTS_cultured - ms.fresh_ductal):+.4f}")
    out[label] = dict(means=m.reindex(ARMS).to_dict(), means_shape=ms.reindex(ARMS).to_dict(),
                      dK=dK, dk=dk, dK_shape=dKs, dk_shape=dks, k_prime=kp, k_prime_shape=kps, g=g,
                      p_KPTS=u('acinar_KPTS_organoid', 'acinar_organoid'),
                      p_K=u('acinar_K_organoid', 'acinar_organoid'),
                      ductal_raw=float(m.ductal_KPTS_cultured - m.fresh_ductal),
                      ductal_shape=float(ms.ductal_KPTS_cultured - ms.fresh_ductal),
                      top1={a: per[per.arm == a].top1.mode().iat[0] for a in ARMS},
                      R2={a: f"{int(per[(per.arm == a) & per.top1.isin(R2)].shape[0])}/{int((per.arm == a).sum())}"
                          for a in ARMS})
    return out


def arm_geometry(Pn, arm, RES, calib):
    """The arm-centroid distance matrix, within-arm spread, and the twin substitutions that replace the
    transformed arm by each non-transforming control arm.  These are the load-bearing statistics of A4 and
    must be reproducible from the repository, not from an ad-hoc shell heredoc."""
    nrm = lambda v: v / (np.linalg.norm(v) + 1e-12)
    cen = {a: nrm(Pn[(arm == a).values].mean(0)) for a in ARMS}
    cs = lambda u, v: 1 - float(u @ v)
    M = pd.DataFrame({a: {b: cs(cen[a], cen[b]) for b in ARMS} for a in ARMS})
    M.to_csv(os.path.join(RES, f'T9e_arm_matrix_{calib}.csv'))
    spread = {a: float(np.mean(1 - Pn[(arm == a).values] @ cen[a])) for a in ARMS}
    d0 = cs(cen['fresh_acinar'], cen['fresh_ductal'])
    subs, dT = {}, 'ductal_KPTS_cultured'
    for a in ('acinar_KPTS_organoid', 'acinar_K_organoid', 'acinar_organoid'):
        subs[a] = dict(d=cs(cen[a], cen[dT]), delta=d0 - cs(cen[a], cen[dT]))
    # leave-one-out on the 3-sample KPTS arm: does it stay closer to ductal_KPTS than the no-driver arm?
    ii = np.where((arm == 'acinar_KPTS_organoid').values)[0]
    base = subs['acinar_organoid']['d']
    loo = [cs(nrm(Pn[np.delete(ii, i)].mean(0)), cen[dT]) for i in range(len(ii))]
    subs['LOO_KPTS_distance'] = [float(x) for x in loo]
    subs['LOO_all_closer_than_no_driver'] = bool(all(x < base for x in loo))
    subs['d_prec'] = d0
    print(f'  arm matrix + substitutions written (T9e_arm_matrix_{calib}.csv); '
          f"culture step {M.loc['fresh_acinar','acinar_organoid']:.3f}, "
          f"driver step {M.loc['acinar_organoid','acinar_KPTS_organoid']:.3f}; "
          f"LOO all closer than no-driver arm: {subs['LOO_all_closer_than_no_driver']}")
    return dict(matrix=M.round(4).to_dict(), spread=spread), subs


def main():
    F = T9.build_frame()
    fz = json.load(open(os.path.join(RES, 'T9_attractor.json')))
    A, A_s, R2 = np.array(fz['attractor_raw']), np.array(fz['attractor_shape']), fz['regions']['R2']
    cts = F['cts']; excl = np.array([cts.index(c) for c in L.EXCL])
    X, arm = load()
    print(f'{X.shape[0]} symbols, {X.shape[1]} samples; arms {dict(arm.value_counts())}')
    out = {}

    # ---------------- readiness scoring, both calibrations ----------------
    for calib in ('all', 'hpa'):
        universe = L.hpa_gene_level()['genes'] if calib == 'hpa' else None
        S, cov, _, _ = T9.score_samples(X, F['panel'], universe)
        pres = sum(a for a, _ in cov.values()); tot = sum(b for _, b in cov.values())
        Pn = L.project(S, F['mu'], F['sd'], F['L'])
        Pns = L.project(S - S.mean(1, keepdims=True), F['mu_s'], F['sd_s'], F['L'])
        t1, _ = L.assign_top1(Pn, F['refPn'], excl)
        per = pd.DataFrame(dict(sample=list(X.columns), arm=arm.values, cos_A=Pn @ A,
                                cos_A_shape=Pns @ A_s, top1=[cts[i] for i in t1]))
        print(f'\n======== readiness, calibration {calib}  (G1 coverage {pres}/{tot} = {pres/tot:.1%})')
        out = analyse(per, R2, f'readiness_{calib}', out)
        for nm, P, Ab in (('raw', Pn, A), ('shape', Pns, A_s)):
            c = cross_origin(P, arm, 'fresh_acinar', 'fresh_ductal',
                             'acinar_KPTS_organoid', 'ductal_KPTS_cultured')
            print(f"  P1 {nm}: precursors {c['d_prec']:.4f} -> transformed {c['d_tr']:.4f}  "
                  f"Delta {c['delta']:+.4f}  perm P = {c['perm_p']:.4f}"
                  f"{'  [G5 uninformative]' if c['G5_uninformative'] else ''}")
            out[f'P1_readiness_{calib}_{nm}'] = c
        per.to_csv(os.path.join(RES, f'T9e_GSE222990_per_sample_{calib}.csv'), index=False)
        out[f'arm_matrix_{calib}'], out[f'substitution_{calib}'] = arm_geometry(Pn, arm, RES, calib)

    # ---------------- activity scoring, declared in advance in A4 note 2 ----------------
    import T9d_activity as T9d
    for est in ('viper', 'ulm'):
        T9d.EST = est
        net, mean_sets, how = T9d.activity_net(list(F['panel']))
        npz = np.load(os.path.join(RES, f'T7_activity_scores_{est}.npz'), allow_pickle=True)
        Apb, Asc = npz['A_pb'], npz['A_sc']
        es = L.eigenspace(Apb, n_pc=L.NPC)
        mu, sd, refPn = L.ref_frame(Asc, es['L'])
        tf = lambda M: M - M.mean(1, keepdims=True)
        mu_s, sd_s, refPn_s = L.ref_frame(tf(np.asarray(Asc, float)), es['L'])
        fza = json.load(open(os.path.join(RES, f'T9d_attractor_{est}.json')))
        Aa, Aas = np.array(fza['attractor']), np.array(fza['attractor_shape'])
        V = X.values.astype(float)
        H = np.log1p(V / (V.sum(0, keepdims=True) + 1e-12) * 1e4).T
        Aex = T9d.estimate(H, list(X.index), net, mean_sets, list(F['panel']))
        Pn = L.project(Aex, mu, sd, es['L'])
        Pns = L.project(Aex - Aex.mean(1, keepdims=True), mu_s, sd_s, es['L'])
        t1, _ = L.assign_top1(Pn, refPn, excl)
        per = pd.DataFrame(dict(sample=list(X.columns), arm=arm.values, cos_A=Pn @ Aa,
                                cos_A_shape=Pns @ Aas, top1=[cts[i] for i in t1]))
        print(f'\n======== activity scoring ({est})')
        out = analyse(per, R2, f'activity_{est}', out)
        for nm, P in (('raw', Pn), ('shape', Pns)):
            c = cross_origin(P, arm, 'fresh_acinar', 'fresh_ductal',
                             'acinar_KPTS_organoid', 'ductal_KPTS_cultured')
            print(f"  P1 {nm}: precursors {c['d_prec']:.4f} -> transformed {c['d_tr']:.4f}  "
                  f"Delta {c['delta']:+.4f}  perm P = {c['perm_p']:.4f}")
            out[f'P1_activity_{est}_{nm}'] = c
        per.to_csv(os.path.join(RES, f'T9e_GSE222990_activity_{est}.csv'), index=False)

    json.dump(out, open(os.path.join(RES, 'T9e_summary.json'), 'w'), indent=1)
    print('\nwrote T9e_summary.json')


if __name__ == '__main__':
    main()
