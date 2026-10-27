"""E39 — "why 44 and not 30 or 100?", answered on normal cells only.
Take random subsets of k of the 44 modules, rebuild the eigenspace from the 154 normal cell types alone, and
measure how much of the normal-cell biology the reduced catalog still recovers: (a) nearest-neighbour
biological-family concordance, (b) Mantel correlation of the reduced cell-cell geometry with the full
44-module geometry. If both saturate well before 44, the catalog size is adequate rather than arbitrary;
if they are still climbing at 44, the paper should say the catalog is not yet saturated."""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import prereg_lib as L

hpa = L.hpa_gene_level()
panel = L.perceptome_panel(); mods = list(panel)
R_full = L.score_ref(hpa['H_sc'], hpa['genes'], panel)
CTS = list(hpa['cts'])
fam = json.load(open(os.path.join(L.ROOT, 'results/fam.json')))
keep = [i for i, c in enumerate(CTS) if c in fam and c not in L.EXCL]
famv = np.array([fam[CTS[i]] for i in keep])
Rk = np.asarray(R_full, float)[keep]
print(f"normal reference: {Rk.shape[0]} cell types x {Rk.shape[1]} modules")

def geometry(Rm, npc=9):
    npc = min(npc, Rm.shape[1])
    mu, sd = Rm.mean(0), Rm.std(0) + 1e-12
    Z = (Rm - mu) / sd
    C = np.cov(Z.T); ev, vec = np.linalg.eigh(C)
    o = np.argsort(-ev); Lm = vec[:, o[:npc]]
    P = Z @ Lm
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    return Pn

def nn_family(Pn):
    S = Pn @ Pn.T; np.fill_diagonal(S, -9)
    nn = np.argmax(S, axis=1)
    return float(np.mean(famv[nn] == famv))

def mantel(Pa, Pb):
    da = 1 - Pa @ Pa.T; db = 1 - Pb @ Pb.T
    iu = np.triu_indices(da.shape[0], 1)
    return float(np.corrcoef(da[iu], db[iu])[0, 1])

full = geometry(Rk); base_nn = nn_family(full)
print(f"full 44-module catalog: NN same-family {100*base_nn:.1f}%\n")
rng = np.random.default_rng(0); rows = []
for k in [5, 10, 15, 20, 25, 30, 35, 40, 44]:
    nns, mts = [], []
    draws = 1 if k == 44 else 25
    for _ in range(draws):
        idx = np.sort(rng.choice(len(mods), k, replace=False)) if k < 44 else np.arange(44)
        P = geometry(Rk[:, idx])
        nns.append(nn_family(P)); mts.append(mantel(P, full))
    rows.append(dict(k=k, nn_mean=np.mean(nns), nn_sd=np.std(nns), mantel_mean=np.mean(mts), mantel_sd=np.std(mts)))
    print(f"  k={k:3d}  NN same-family {100*np.mean(nns):5.1f}% (sd {100*np.std(nns):4.1f})"
          f"  | geometry agreement with full catalog r = {np.mean(mts):.3f} (sd {np.std(mts):.3f})")
D = pd.DataFrame(rows); D.to_csv(os.path.join(L.ROOT, 'results/referee2_newdata/E39_saturation.csv'), index=False)
for col, lab in [('nn_mean', 'family recovery'), ('mantel_mean', 'geometry agreement')]:
    v = D.set_index('k')[col]
    reach = next((k for k in v.index if v[k] >= 0.95 * v[44]), None)
    print(f"\n  {lab}: reaches 95% of the full-catalog value at k = {reach}")
print("DONE")
