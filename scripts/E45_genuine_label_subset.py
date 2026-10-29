"""Referee request (PNAS v13 report, §3.3): the target is named by cross-cancer recurrence over the
25 cohorts, and Note S-U reports that in 8 of the 14 Census cohorts "malignant" is a cell-type label
rather than a malignancy annotation. If a third of the panel's malignant compartments are partly
normal cells of the tissue, the recurrence ranking is partly a ranking over normal-cell projections,
and the direction-of-bias defence is not available (Note S-I tested it and reports it not borne out).

So: recompute the recurrence ranking and the Table 1 boundaries on the subset whose malignant
compartment carries a genuine author/atlas annotation, and report whether the top-2 and the region
membership survive. This is NOT the leave-one-cohort-out analysis (Note S-L(iv)): the concern is a
correlated class of eight, not one cohort at a time.

Identical geometry, standardization, quarantine and top-5 rule as E40_boundary_enrichment.py — the
only change is which cohorts enter the recurrence count. Writes results/cancer_v6_6_extras/E45_*.csv.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, (ROOT / 'tool/perceptome2').as_posix())
from perceptome.perceptivity.reference import load_hpa_perceptivity

# --- the eight heuristically-labelled Census cohorts, from results/census_cohorts.csv -------------
CEN = pd.read_csv(ROOT / 'results/census_cohorts.csv')
# cohort names in census_cohorts.csv vs the ids in E7_S_raw.npz
ALIAS = {'Ovarian serous': 'OvarianSerous', 'Squamous lung': 'SqLung'}
HEUR = {ALIAS.get(c, c) for c in CEN.loc[CEN.selection == 'heuristic cell-type label', 'cohort']}

e = json.load(open(ROOT / 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json'))
mods = e['module_order']
load = np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm = load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; cells = list(Rm.index)
mu, sd = Rm.mean(0).values, Rm.std(0).values + 1e-12
refP = ((Rm.values - mu) / sd) @ load
refn = refP / (np.linalg.norm(refP, axis=1, keepdims=True) + 1e-12)
EX = [cells.index(c) for c in ('erythrocytes', 'platelets')]
NCAND = len(cells) - len(EX)

d = np.load(ROOT / 'results/E7_S_raw.npz', allow_pickle=True)
S, cid = d['S'].astype(float), d['cid']
pc = ((S - mu) / sd) @ load
pcn = pc / (np.linalg.norm(pc, axis=1, keepdims=True) + 1e-12)
sim = pcn @ refn.T; sim[:, EX] = -1e9
top = np.argmax(sim, axis=1)

ALL = sorted(set(cid.tolist()))
KEEP = [c for c in ALL if c not in HEUR]
assert len(HEUR & set(ALL)) == 8, f'expected 8 heuristic cohorts in the panel, got {sorted(HEUR & set(ALL))}'
print(f'panel {len(S):,} cells / {len(ALL)} cohorts   ->   genuine-label subset {len(KEEP)} cohorts')
print('dropped (heuristic cell-type label):', ', '.join(sorted(HEUR & set(ALL))), '\n')

def recurrence(cohorts):
    rec = np.zeros(len(cells))
    for c in cohorts:
        t = top[cid == c]
        b = np.bincount(t, minlength=len(cells)); b[EX] = -1
        rec[np.argsort(-b)[:5]] += 1
    rec[EX] = -1
    return rec

rows = []
rec_all, rec_keep = recurrence(ALL), recurrence(KEEP)
o_all, o_keep = np.argsort(-rec_all), np.argsort(-rec_keep)
rank_keep = {i: r for r, i in enumerate(o_keep, 1)}
print(f'{"cell type":34s} {"n/25":>6} {"rank":>5}   {"n/17":>6} {"rank":>5}   fraction')
for i in o_all[:16]:
    print(f'   {cells[i]:34s} {int(rec_all[i]):2d}/25 {list(o_all).index(i)+1:4d}   '
          f'{int(rec_keep[i]):2d}/{len(KEEP)} {rank_keep[i]:4d}   '
          f'{rec_all[i]/len(ALL):.2f} -> {rec_keep[i]/len(KEEP):.2f}')
    rows.append(dict(cell_type=cells[i], rec_25=int(rec_all[i]), rank_25=list(o_all).index(i) + 1,
                     rec_17=int(rec_keep[i]), rank_17=rank_keep[i],
                     frac_25=round(rec_all[i] / len(ALL), 3), frac_17=round(rec_keep[i] / len(KEEP), 3)))
pd.DataFrame(rows).to_csv(ROOT / 'results/cancer_v6_6_extras/E45_recurrence_genuine.csv', index=False)

# --- Table 1 boundaries, both on the same k (comparability) and on the subset's own majority rule --
def boundaries(cohorts, ks, tag):
    sel = np.isin(cid, list(cohorts)); t = top[sel]
    rec = recurrence(cohorts); order = np.argsort(-rec)
    out = []
    for nm, k in ks:
        idx = list(order[:k])
        occ = 100 * np.isin(t, idx).mean(); ch = 100 * k / NCAND
        out.append(dict(panel=tag, boundary=nm, k=k, occupancy_pct=round(occ, 1),
                        chance_pct=round(ch, 2), fold=round(occ / ch, 1),
                        members='; '.join(cells[i] for i in idx)))
        print(f'   {tag:12s} {nm:22s} k={k:2d}  occupancy {occ:5.1f}%  chance {ch:4.2f}%  '
              f'enrichment {occ/ch:5.1f}x')
    return out

print('\nTable 1 boundaries:')
KS = [('core', 2), ('neighbourhood', 4), ('broader region', 14)]
B = boundaries(ALL, KS, '25 cohorts') + boundaries(KEEP, KS, '17 genuine')
# the subset's own majority rule: recurrent in >= ceil(17/2) = 9 of 17, and >= 5 of 17 for the broad band
rk = recurrence(KEEP)
k_major = int((rk >= np.ceil(len(KEEP) / 2)).sum()); k_broad = int((rk >= 5).sum())
print(f'\n   the subset\'s own recurrence rule: majority (>={int(np.ceil(len(KEEP)/2))}/{len(KEEP)}) '
      f'-> {k_major} cells;  >=5/{len(KEEP)} -> {k_broad} cells')
B += boundaries(KEEP, [('neighbourhood (own majority rule)', k_major),
                       ('broader region (own >=5 rule)', k_broad)], '17 own-rule')
pd.DataFrame(B).to_csv(ROOT / 'results/cancer_v6_6_extras/E45_boundaries_genuine.csv', index=False)

# --- does the region membership survive? ----------------------------------------------------------
m25, m17 = set(o_all[:14]), set(o_keep[:14])
print(f'\nregion membership (top-14): {len(m25 & m17)}/14 shared')
print('   lost  :', ', '.join(cells[i] for i in o_all[:14] if i not in m17) or 'none')
print('   gained:', ', '.join(cells[i] for i in o_keep[:14] if i not in m25) or 'none')
print(f'core (top-2) on 25: {cells[o_all[0]]}, {cells[o_all[1]]}')
print(f'core (top-2) on 17: {cells[o_keep[0]]}, {cells[o_keep[1]]}')

# --- the circularity the referee names specifically -----------------------------------------------
# Ovarian "malignant" cells are selected as fallopian tube secretory epithelial cells, and the region
# the paper names is "the same neighbourhood of high-output secretory and progenitor states". If that
# selection were feeding its own destination, ovarian would sit on the core far harder than the panel.
print('\nper-cohort top-1 target and core occupancy, the eight dropped cohorts vs the panel:')
core25 = list(o_all[:2])
det = []
for c in ALL:
    t = top[cid == c]
    b = np.bincount(t, minlength=len(cells)); b[EX] = -1
    det.append(dict(cohort=c, heuristic=c in HEUR, n=int(len(t)),
                    top1=cells[int(np.argmax(b))],
                    top1_pct=round(100 * b.max() / len(t), 1),
                    core_pct=round(100 * np.isin(t, core25).mean(), 1),
                    label=CEN.loc[CEN.cohort.map(lambda x: ALIAS.get(x, x)) == c, 'malignant_label'].iloc[0]
                          if c in set(CEN.cohort.map(lambda x: ALIAS.get(x, x))) else 'author (GEO/atlas)'))
D = pd.DataFrame(det).sort_values(['heuristic', 'core_pct'], ascending=[False, False])
for r in D.itertuples():
    flag = 'drop' if r.heuristic else 'keep'
    print(f'   {flag}  {r.cohort:14s} n={r.n:>6,}  core {r.core_pct:5.1f}%  '
          f'top1 {r.top1:26s} {r.top1_pct:5.1f}%   [{str(r.label)[:46]}]')
D.to_csv(ROOT / 'results/cancer_v6_6_extras/E45_per_cohort_labels.csv', index=False)
h, k = D[D.heuristic].core_pct, D[~D.heuristic].core_pct
print(f'\n   core occupancy: heuristic-label cohorts median {h.median():.1f}% (n={len(h)}), '
      f'genuine-label median {k.median():.1f}% (n={len(k)})')
ov = D[D.cohort == 'OvarianSerous'].iloc[0]
print(f'   ovarian specifically: core {ov.core_pct}%, top-1 {ov.top1} ({ov.top1_pct}%) '
      f'- panel median {D.core_pct.median():.1f}%')

# --- how the three boundaries were chosen, and the tie at the widest one -------------------------
# Found 2026-10-29: the rule Note S-K states does not produce the k the paper reports, and the
# fourteen-cell membership passes through a five-way tie, which two scripts broke differently
# (E41 -> 76.9%, E40 -> 77.7% for the same quantity). Both facts are now registered claims.
import itertools
rows2 = []
for lab, k in [('core: the two leading both rankings', 2),
               ('neighbourhood: recurrent in a majority (>=13 of 25)', int((rec_all >= 13).sum())),
               ('broader: recurrent in >=5 of 25', int((rec_all >= 5).sum()))]:
    idx = list(o_all[:k]); occ = 100 * np.isin(top, idx).mean(); ch = 100 * k / NCAND
    rows2.append(dict(definition='rule', boundary=lab, k=k, occupancy_pct=round(occ, 1),
                      chance_pct=round(ch, 2), fold=round(occ / ch, 1)))
for lab, k in [('core (top-2)', 2), ('neighbourhood (top-4)', 4), ('broader region (top-14)', 14)]:
    idx = list(o_all[:k]); occ = 100 * np.isin(top, idx).mean(); ch = 100 * k / NCAND
    rows2.append(dict(definition='reported k', boundary=lab, k=k, occupancy_pct=round(occ, 1),
                      chance_pct=round(ch, 2), fold=round(occ / ch, 1)))
hard = [i for i in o_all if rec_all[i] > 2][:12]
tied = [i for i in range(len(cells)) if int(rec_all[i]) == 2]
vals = []
for pair in itertools.combinations(tied, 2):
    idx = hard + list(pair); o = 100 * np.isin(top, idx).mean()
    vals.append((o, o / (100 * 14 / NCAND)))
rows2.append(dict(definition='top-14 tie', boundary=f'{len(tied)}-way tie at rec=2, choose 2',
                  k=len(vals), occupancy_pct=round(min(v[0] for v in vals), 1),
                  chance_pct=round(max(v[0] for v in vals), 1),
                  fold=round(min(v[1] for v in vals), 1)))
rows2.append(dict(definition='top-14 tie (max)', boundary='highest of the ten orderings', k=len(vals),
                  occupancy_pct=round(max(v[0] for v in vals), 1), chance_pct=float('nan'),
                  fold=round(max(v[1] for v in vals), 1)))
pd.DataFrame(rows2).to_csv(ROOT / 'results/cancer_v6_6_extras/E45_boundary_choice.csv', index=False)
print('\nboundary definitions: rule vs reported k, and the tie range')
for r in rows2:
    print(f'   {r["definition"]:18s} {r["boundary"][:46]:48s} k={r["k"]:3d}  '
          f'occ {r["occupancy_pct"]:5.1f}  fold {r["fold"]}')

print('DONE')
