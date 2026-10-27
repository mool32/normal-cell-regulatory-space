"""T8 (exploratory, not pre-registered) — Gudkov/Teodor question: when the three independent curations (T5)
name GASTRIC PROGENITOR instead of the cytotrophoblast, is the cytotrophoblast FAR, or is it the same
neighbourhood with a different winner? I.e. does the catalog determine the REGION or only the LABEL?
Reports, for the incumbent and for C2/C3/C4:
  (a) the cosine gap structure of the ranking — is rank 27 'far' or inside a flat band?
  (b) where the cytotrophoblast sits by COSINE to the malignant attractor (not just by recurrence rank)
  (c) cosine between that catalog's winner and the cytotrophoblast, in that catalog's own space
  (d) overlap of the top-10 target sets across catalogs (is the region reproducible?)
  (e) the trophoblast FAMILY as a whole (cyto + migrating cyto + syncytio + extravillous) vs the winner."""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import prereg_lib as L

TROPH_KEYS = ('cytotrophoblast', 'trophoblast', 'syncytio', 'extravillous')
def troph_members(cts): return [c for c in cts if any(k in c.lower() for k in TROPH_KEYS)]

def analyse(tag, panel, cache, hpa):
    ep = L.run_panel(panel, cache, hpa=hpa, return_all=True)
    cts = hpa['cts']; refPn = ep['refPn']; Pn = ep['Pn']; top1 = ep['top1']
    excl = [cts.index(c) for c in L.EXCL]
    cid = cache['obs']['cancer'].astype(str).values; uc = sorted(set(cid))
    # malignant attractor = mean unit direction of the 25 cohort centroids
    M = np.vstack([Pn[cid == c].mean(0) for c in uc]); M /= np.linalg.norm(M, axis=1, keepdims=True)
    att = M.mean(0); att /= np.linalg.norm(att)
    cos = refPn @ att
    for i in excl: cos[i] = -9
    order = np.argsort(-cos)
    rank_cos = {cts[i]: 1 + n for n, i in enumerate(order)}
    # recurrence
    rec = np.zeros(len(cts))
    for c in uc:
        t = top1[cid == c]; b = np.bincount(t, minlength=len(cts)); b[excl] = -1
        rec[np.argsort(-b)[:5]] += 1
    occ = 100 * np.bincount(top1, minlength=len(cts)) / len(top1)
    rec_order = np.argsort(-rec)
    rank_rec = {cts[i]: 1 + n for n, i in enumerate(rec_order)}
    winner = cts[rec_order[0]]
    cyto = 'cytotrophoblasts'
    troph = troph_members(cts)
    out = dict(tag=tag, winner=winner, winner_rec=int(rec[cts.index(winner)]),
               cyto_rank_rec=rank_rec[cyto], cyto_rec=int(rec[cts.index(cyto)]),
               cyto_rank_cos=rank_cos[cyto], cyto_cos=float(cos[cts.index(cyto)]),
               winner_cos=float(cos[cts.index(winner)]),
               cos_winner_cyto=float(refPn[cts.index(winner)] @ refPn[cts.index(cyto)]),
               occ_cyto=float(occ[cts.index(cyto)]), occ_winner=float(occ[cts.index(winner)]),
               troph_family_occ=float(sum(occ[cts.index(t)] for t in troph)),
               troph_family_best_rank_cos=min(rank_cos[t] for t in troph),
               top10_cos=[cts[i] for i in order[:10]], top10_rec=[cts[i] for i in rec_order[:10]],
               cos_rank1=float(cos[order[0]]), cos_rank10=float(cos[order[9]]),
               cos_rank30=float(cos[order[29]]), cos_rank76=float(cos[order[75]]), cos_rank152=float(cos[order[-1]]))
    return out, cos, cts, refPn, order

if __name__ == '__main__':
    hpa = L.hpa_gene_level(); cache = L.load_panel_cache()
    panels = {'incumbent': L.perceptome_panel()}
    for c in ['C2', 'C3', 'C4']:
        d = json.load(open(os.path.join(L.REG, f'T5_catalog_{c}.json')))
        panels[c] = {m: [str(g).upper() for g in (gl['genes'] if isinstance(gl, dict) else gl)]
                     for m, gl in d['modules'].items()}
    rows = []; keep = {}
    for tag, p in panels.items():
        L.stage(f'{tag}: {len(p)} modules / {len(set(sum(p.values(),[])))} genes')
        r, cos, cts, refPn, order = analyse(tag, p, cache, hpa); rows.append(r); keep[tag] = (cos, cts, refPn, order)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(L.RES, 'T8_region_vs_label.csv'), index=False)

    print('\n' + '=' * 100)
    print('IS THE CYTOTROPHOBLAST *FAR* UNDER THE INDEPENDENT CURATIONS, OR THE SAME REGION WITH ANOTHER WINNER?')
    print('=' * 100)
    for r in rows:
        print(f"\n[{r['tag']}]  winner = {r['winner']} (recurrence {r['winner_rec']}/25, occupancy {r['occ_winner']:.1f}%)")
        print(f"   cytotrophoblast: recurrence rank {r['cyto_rank_rec']:3d}/152 ({r['cyto_rec']}/25)   "
              f"COSINE rank {r['cyto_rank_cos']:3d}/152 (cos {r['cyto_cos']:+.3f})   occupancy {r['occ_cyto']:.1f}%")
        print(f"   winner cos to malignant attractor {r['winner_cos']:+.3f}   cos(winner, cytotrophoblast) = {r['cos_winner_cyto']:+.3f}")
        print(f"   trophoblast FAMILY: total occupancy {r['troph_family_occ']:.1f}%, best cosine rank {r['troph_family_best_rank_cos']}/152")
        print(f"   cosine band: rank1 {r['cos_rank1']:.3f} | rank10 {r['cos_rank10']:.3f} | rank30 {r['cos_rank30']:.3f} "
              f"| median {r['cos_rank76']:.3f} | last {r['cos_rank152']:.3f}")
        print(f"   top-10 by cosine: {', '.join(r['top10_cos'][:6])} …")

    print('\n' + '=' * 100); print('REGION REPRODUCIBILITY — overlap of the top-10 target sets')
    tags = [r['tag'] for r in rows]
    for a in tags:
        A = set(dict(zip([x['tag'] for x in rows], [x['top10_cos'] for x in rows]))[a])
        line = [f"{len(A & set(dict(zip([x['tag'] for x in rows],[x['top10_cos'] for x in rows]))[b]))}/10" for b in tags]
        print(f"   {a:10s} vs {tags}: {line}")
    inc = set(dict(zip([x['tag'] for x in rows], [x['top10_cos'] for x in rows]))['incumbent'])
    shared = set.intersection(*[set(x['top10_cos']) for x in rows])
    print(f"\n   shared by ALL FOUR catalogs (top-10 by cosine): {sorted(shared) if shared else 'none'}")
    print('DONE')
