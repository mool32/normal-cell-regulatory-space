"""E38 — CANCER-FREE justification of the module catalog, computed on the 154 normal cell types only.
Gudkov's flow forbids cancer from appearing while the modules are introduced, so the robustness controls
(module-drop, leave-one-cancer-out, gene shuffle) cannot be used there: all are measured in cancer units.
This script supplies a justification that uses no cancer data at all.

A. WHY NEGATIVE FEEDBACK (criterion iii), restated mechanistically: a module is admissible only if its
   activity is READABLE FROM RNA. The feedback element is induced by the module's own output, so its
   abundance tracks pathway activity, whereas receptors and kinases are largely structural and act
   post-translationally. Prediction, testable on normal cells: feedback/target genes should vary across
   normal cell types substantially more than sensor/cascade genes of the same modules.
B. WHY THESE 44 (the selection), without cancer: the retained modules should be more coherent internally
   and more discriminative of normal cell identity than off-the-shelf pathway collections put through the
   identical scoring, and than size-matched random gene sets."""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import prereg_lib as L
from scipy.stats import mannwhitneyu

hpa = L.hpa_gene_level()
H, GEN, CTS = hpa['H_sc'], list(hpa['genes']), list(hpa['cts'])
gi = {g: i for i, g in enumerate(GEN)}
X = np.asarray(H, float)
if X.shape[0] == len(GEN): X = X.T          # cells x genes
print(f"HPA normal reference: {X.shape[0]} cell types x {X.shape[1]} genes")
sd_all = X.std(0)

D = json.load(open(os.path.join(L.ROOT, 'tool/perceptome2/perceptome/catalog/data/modules_v03.json')))['modules']
up = lambda l: [str(g).upper() for g in (l or [])]
rows = []
for m, v in D.items():
    for role in ['sensor_genes', 'cascade_genes', 'tf_genes', 'feedback_genes', 'activity_genes', 'core_genes']:
        for g in up(v.get(role)):
            if g in gi: rows.append(dict(module=m, role=role, gene=g, sd=float(sd_all[gi[g]])))
G = pd.DataFrame(rows).drop_duplicates(['module', 'role', 'gene'])
print("\n" + "=" * 92)
print("A. IS THE FEEDBACK/TARGET ARM THE TRANSCRIPTIONALLY READABLE PART? (normal cells only)")
print("=" * 92)
summ = G.groupby('role')['sd'].agg(['count', 'median', 'mean']).round(3).sort_values('median', ascending=False)
print(summ.to_string())
induced = G[G.role.isin(['feedback_genes', 'activity_genes'])]['sd']
structural = G[G.role.isin(['sensor_genes', 'cascade_genes'])]['sd']
u, p = mannwhitneyu(induced, structural, alternative='greater')
print(f"\n  induced arm (feedback + TF targets) n={len(induced)}  median SD {induced.median():.3f}")
print(f"  structural arm (sensors + cascade)   n={len(structural)}  median SD {structural.median():.3f}")
print(f"  induced > structural across normal cell types: Mann-Whitney P = {p:.2e}")
fb = G[G.role == 'feedback_genes']['sd']; cas = G[G.role == 'cascade_genes']['sd']
print(f"  feedback alone vs cascade alone: {fb.median():.3f} vs {cas.median():.3f}, P = {mannwhitneyu(fb, cas, alternative='greater')[1]:.2e}")

print("\n" + "=" * 92)
print("B. ARE THE RETAINED MODULES MORE COHERENT AND MORE DISCRIMINATIVE THAN ALTERNATIVES? (no cancer)")
print("=" * 92)
fam = json.load(open(os.path.join(L.ROOT, 'results/fam.json')))
keep = [i for i, c in enumerate(CTS) if c in fam and c not in L.EXCL]
famv = np.array([fam[CTS[i]] for i in keep])

def coherence(genes):
    idx = [gi[g] for g in genes if g in gi]
    if len(idx) < 3: return np.nan
    A = X[:, idx]
    A = A[:, A.std(0) > 1e-9]
    if A.shape[1] < 3: return np.nan
    C = np.corrcoef(A.T); iu = np.triu_indices(C.shape[0], 1)
    return float(np.nanmean(C[iu]))

def eta2(score):
    s = score[keep]
    gm = s.mean(); sst = ((s - gm) ** 2).sum()
    if sst <= 0: return np.nan
    ssb = sum(((s[famv == f].mean() - gm) ** 2) * (famv == f).sum() for f in set(famv))
    return float(ssb / sst)

def panel_stats(panel, tag):
    coh, et, n = [], [], []
    for mname, genes in panel.items():
        gs = [g for g in up(genes) if g in gi]
        if len(gs) < 3: continue
        c = coherence(gs)
        sc = X[:, [gi[g] for g in gs]].mean(1)
        coh.append(c); et.append(eta2(sc)); n.append(len(gs))
    coh = np.array(coh, float); et = np.array(et, float)
    print(f"  {tag:26s} modules {len(coh):3d}  median genes {int(np.median(n)):4d}  "
          f"median within-module coherence {np.nanmedian(coh):+.3f}  median eta^2 vs family {np.nanmedian(et):.3f}")
    return coh, et

panels = {'perceptome-44 (ours)': L.perceptome_panel()}
for nm in ['kegg-signalling-46', 'reactome-top44', 'hallmark-50', 'progeny-14']:
    try: panels[nm] = L.load_panel(nm)
    except Exception as ex: print(f"  (panel {nm} unavailable: {str(ex)[:50]})")
res = {}
for nm, p in panels.items(): res[nm] = panel_stats(p, nm)

rng = np.random.default_rng(0)
ours = L.perceptome_panel(); sizes = [len([g for g in up(v) if g in gi]) for v in ours.values()]
expressed = [g for g in GEN if sd_all[gi[g]] > np.percentile(sd_all, 25)]
rc, re_ = [], []
for _ in range(200):
    gs = list(rng.choice(expressed, int(np.median(sizes)), replace=False))
    rc.append(coherence(gs)); re_.append(eta2(X[:, [gi[g] for g in gs]].mean(1)))
print(f"  {'random size-matched x200':26s} modules 200  median genes {int(np.median(sizes)):4d}  "
      f"median within-module coherence {np.nanmedian(rc):+.3f}  median eta^2 vs family {np.nanmedian(re_):.3f}")

oc, oe = res['perceptome-44 (ours)']
print(f"\n  ours vs random — coherence P = {mannwhitneyu(oc[~np.isnan(oc)], np.array(rc)[~np.isnan(rc)], alternative='greater')[1]:.1e}"
      f" | eta^2 P = {mannwhitneyu(oe[~np.isnan(oe)], np.array(re_)[~np.isnan(re_)], alternative='greater')[1]:.1e}")
for nm in panels:
    if nm == 'perceptome-44 (ours)': continue
    c, e = res[nm]
    pc = mannwhitneyu(oc[~np.isnan(oc)], c[~np.isnan(c)], alternative='greater')[1]
    pe = mannwhitneyu(oe[~np.isnan(oe)], e[~np.isnan(e)], alternative='greater')[1]
    print(f"  ours vs {nm:22s} coherence P = {pc:.1e} | eta^2 P = {pe:.1e}")
out = os.path.join(L.ROOT, 'results/referee2_newdata')
G.to_csv(os.path.join(out, 'E38_gene_roles_sd.csv'), index=False)
print("\nDONE")
