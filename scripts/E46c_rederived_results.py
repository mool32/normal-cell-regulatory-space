"""Referee request (PNAS v13 report, §3.1), stage 4: the three numbers asked for, recomputed on the
re-derived catalog -- (i) core fraction and enrichment, (ii) cross-cancer recurrence top-5.
(iii) the colorectal within-patient increase is in E46d.

Everything is computed twice from the same code: once on the shipped reference, eigenspace and score
matrix, and once on the fixed ones (results/E46_rederivation/). The shipped arm must reproduce the
paper's published values, which is the test that this comparison is measuring the fix and not the
re-implementation.

The paper's frame is nine-dimensional. The fixed reference's Kaiser cutoff is TEN, so the 9-PC
projection is kept for comparability and the tenth component is reported separately rather than
silently changing the frame.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, (ROOT / 'tool/perceptome2').as_posix())
from perceptome.perceptivity.reference import load_hpa_perceptivity
FIX = ROOT / 'results/E46_rederivation'

def arm(tag, eig_path, ref_npz, S_npz, npc=9):
    e = json.load(open(eig_path)); mods = e['module_order']
    load = np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(npc)] for m in mods])
    if ref_npz is None:
        Rm = load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
        cells = list(Rm.index); Rv = Rm.values
    else:
        meta = json.load(open(str(ref_npz).replace('.npz', '.json')))
        cells = meta['cell_types']
        Rv = pd.DataFrame(np.load(ref_npz)['R'], index=cells, columns=meta['modules'])[mods].values
    mu, sd = Rv.mean(0), Rv.std(0) + 1e-12
    refP = ((Rv - mu) / sd) @ load
    refn = refP / (np.linalg.norm(refP, axis=1, keepdims=True) + 1e-12)
    EX = [cells.index(c) for c in ('erythrocytes', 'platelets')]
    d = np.load(S_npz, allow_pickle=True); S, cid = d['S'].astype(float), d['cid']
    pc = ((S - mu) / sd) @ load
    pcn = pc / (np.linalg.norm(pc, axis=1, keepdims=True) + 1e-12)
    sim = pcn @ refn.T; sim[:, EX] = -1e9
    top = np.argmax(sim, axis=1)
    rec = np.zeros(len(cells))
    for c in sorted(set(cid.tolist())):
        t = top[cid == c]; b = np.bincount(t, minlength=len(cells)); b[EX] = -1
        rec[np.argsort(-b)[:5]] += 1
    rec[EX] = -1
    return dict(tag=tag, cells=cells, top=top, cid=cid, rec=rec, EX=EX,
                ncand=len(cells) - len(EX), order=np.argsort(-rec))

A = arm('shipped', ROOT / 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json', None,
        ROOT / 'results/E7_S_raw.npz')
B = arm('re-derived', FIX / 'reference_FIX.json', FIX / 'hpa_perceptivity_FIX_scscaled.npz',
        FIX / 'E46_S_fix.npz')
print(f'cells: shipped {len(A["top"]):,}  re-derived {len(B["top"]):,}  '
      f'(the fix recovers 573 real symbols per cohort, so the detected-gene QC admits one more cell)')

# ---------------- (ii) cross-cancer recurrence ----------------------------------------------------
print('\n(ii) cross-cancer recurrence, top 8')
print(f'    {"cell type":32s} {"shipped":>14s}   {"re-derived":>14s}')
rows = []
for i in A['order'][:8]:
    j = B['cells'].index(A['cells'][i])
    ra = list(A['order']).index(i) + 1; rb = list(B['order']).index(j) + 1
    print(f'    {A["cells"][i]:32s} {int(A["rec"][i]):2d}/25 rank {ra:<3d}   '
          f'{int(B["rec"][j]):2d}/25 rank {rb:<3d}')
    rows.append(dict(cell_type=A['cells'][i], rec_shipped=int(A['rec'][i]), rank_shipped=ra,
                     rec_rederived=int(B['rec'][j]), rank_rederived=rb))
pd.DataFrame(rows).to_csv(FIX / 'E46_recurrence_compare.csv', index=False)

# ---------------- (i) core fraction and enrichment ------------------------------------------------
print('\n(i) occupancy and enrichment at the three boundaries')
out = []
for X in (A, B):
    anchor = [X['order'][0]] if X['cells'][X['order'][0]] == 'cytotrophoblasts' else \
             [int(np.where(np.array(X['cells']) == 'cytotrophoblasts')[0][0])]
    cyto = int(np.where(np.array(X['cells']) == 'cytotrophoblasts')[0][0])
    for nm, idx in [('cytotrophoblast alone', [cyto]), ('core (top-2)', list(X['order'][:2])),
                    ('neighbourhood (top-4)', list(X['order'][:4])),
                    ('broader region (top-14)', list(X['order'][:14]))]:
        occ = 100 * np.isin(X['top'], idx).mean(); ch = 100 * len(idx) / X['ncand']
        out.append(dict(panel=X['tag'], boundary=nm, k=len(idx), occupancy_pct=round(occ, 1),
                        chance_pct=round(ch, 2), fold=round(occ / ch, 1),
                        members='; '.join(X['cells'][i] for i in idx)))
D = pd.DataFrame(out)
for nm in D.boundary.unique():
    s = D[(D.panel == 'shipped') & (D.boundary == nm)].iloc[0]
    r = D[(D.panel == 're-derived') & (D.boundary == nm)].iloc[0]
    print(f'    {nm:26s} {s.occupancy_pct:5.1f}% ({s.fold:4.1f}x)   ->   '
          f'{r.occupancy_pct:5.1f}% ({r.fold:4.1f}x)')
D.to_csv(FIX / 'E46_boundaries_compare.csv', index=False)
same = [A['cells'][i] for i in A['order'][:14]]
newm = [B['cells'][i] for i in B['order'][:14]]
print(f'\n    region membership: {len(set(same) & set(newm))}/14 shared')
print(f'    lost  : {", ".join(c for c in same if c not in newm) or "none"}')
print(f'    gained: {", ".join(c for c in newm if c not in same) or "none"}')
print(f'    core  : {same[:2]}  ->  {newm[:2]}')

# ---------------- the tenth component ------------------------------------------------------------
e = json.load(open(FIX / 'reference_FIX.json'))
print(f'\n    Kaiser cutoff on the re-derived reference: {e["n_kaiser"]} components '
      f'(the paper uses nine; PC10 eigenvalue {e["eigenvalues"][9]:.3f})')
B10 = arm('re-derived, 10-PC', FIX / 'reference_FIX.json',
          FIX / 'hpa_perceptivity_FIX_scscaled.npz', FIX / 'E46_S_fix.npz', npc=10)
c10 = int(np.where(np.array(B10['cells']) == 'cytotrophoblasts')[0][0])
print(f'    on a ten-dimensional frame: cytotrophoblast alone '
      f'{100*np.mean(B10["top"]==c10):.1f}%, core (top-2) '
      f'{100*np.isin(B10["top"], list(B10["order"][:2])).mean():.1f}%, '
      f'top-2 = {[B10["cells"][i] for i in B10["order"][:2]]}')
print('DONE')
