"""M3 — self-contained justification of criterion (iii) inside paper 3.
Question a referee will ask: why make NEGATIVE FEEDBACK the primary criterion for admitting a module?
Answer, computed here from two public resources and this paper's own catalog, with no dependence on
unpublished companion work: genes that occupy the negative-feedback position of a module are enriched for
cancer-driver status RELATIVE TO THE OTHER GENES OF THE SAME MODULES — i.e. independently of pathway
membership, which is the comparison that matters. Resources: COSMIC Cancer Gene Census v103 (740 genes) and
the module catalog (modules_v03.json). Reports the within-catalog test (primary), a genome background for
scale, a TF-output comparison (both are well-studied named classes, controlling for study bias), and a
per-module breakdown."""
import json, sys
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from pathlib import Path
ROOT = Path('.').resolve()

cgc = set()
for line in open(ROOT / 'paper2/data/cosmic_cgc_full.txt', encoding='utf8'):
    s = line.strip().split('\t')[0].strip().strip('"')
    if s and not s.lower().startswith('gene'):
        cgc.add(s.upper())
print(f"COSMIC CGC v103: {len(cgc)} genes")

D = json.load(open(ROOT / 'tool/perceptome2/perceptome/catalog/data/modules_v03.json'))['modules']
ROLE = ['sensor_genes', 'cascade_genes', 'tf_genes', 'feedback_genes', 'core_genes', 'activity_genes']
up = lambda L: {str(g).upper() for g in (L or [])}

ann = [m for m, v in D.items() if v.get('feedback_genes')]
FB, OTHER, TF = set(), set(), set()
per = []
for m in ann:
    v = D[m]
    fb = up(v.get('feedback_genes'))
    rest = set().union(*[up(v.get(r)) for r in ROLE]) - fb
    FB |= fb; OTHER |= rest; TF |= up(v.get('tf_genes'))
    per.append(dict(module=m, n_fb=len(fb), fb_cgc=len(fb & cgc), n_rest=len(rest), rest_cgc=len(rest & cgc)))
OTHER -= FB; TF -= FB
print(f"modules with an annotated feedback element: {len(ann)}/44")

def test(A, B, la, lb):
    a, b = len(A & cgc), len(A - cgc)
    c, d = len(B & cgc), len(B - cgc)
    orr, p = fisher_exact([[a, b], [c, d]])
    pa, pb = 100 * a / max(len(A), 1), 100 * c / max(len(B), 1)
    print(f"  {la:34s} {a:4d}/{len(A):4d} = {pa:5.1f}%   vs  {lb:30s} {c:4d}/{len(B):5d} = {pb:4.1f}%"
          f"   OR = {orr:6.1f}   P = {p:.2e}")
    return dict(test=f"{la} vs {lb}", a=a, nA=len(A), c=c, nB=len(B), pct_A=pa, pct_B=pb, OR=orr, P=p)

print("\n" + "=" * 118)
print("DOES OCCUPYING THE NEGATIVE-FEEDBACK POSITION PREDICT CANCER-DRIVER STATUS?")
print("=" * 118)
rows = []
rows.append(test(FB, OTHER, "feedback genes of the modules", "all other genes, SAME modules"))
print("     ^ the primary, within-pathway control: both sets come from the same pathways, so pathway membership cannot explain it")
GENOME = 20000
g_cgc = len(cgc)
a, b = len(FB & cgc), len(FB - cgc)
orr, p = fisher_exact([[a, b], [g_cgc - a, GENOME - g_cgc - b]])
print(f"  {'feedback genes':34s} {a:4d}/{len(FB):4d} = {100*a/len(FB):5.1f}%   vs  {'genome background':30s} {g_cgc:4d}/{GENOME:5d} = {100*g_cgc/GENOME:4.1f}%   OR = {orr:6.1f}   P = {p:.2e}")
rows.append(dict(test="feedback vs genome", a=a, nA=len(FB), c=g_cgc, nB=GENOME, pct_A=100*a/len(FB), pct_B=100*g_cgc/GENOME, OR=orr, P=p))
rows.append(test(FB, TF, "feedback genes", "TF-output genes, same modules"))
print("     ^ study-bias control: TF-output genes are an equally named, equally well-studied class")

print("\nper-module breakdown (modules with an annotated feedback element):")
P = pd.DataFrame(per).sort_values('fb_cgc', ascending=False)
P['fb_pct'] = (100 * P.fb_cgc / P.n_fb).round(0)
P['rest_pct'] = (100 * P.rest_cgc / P.n_rest.clip(lower=1)).round(0)
print(P.to_string(index=False))
higher = int((P.fb_pct > P.rest_pct).sum())
print(f"\n  feedback arm richer in CGC than the rest of its own module in {higher}/{len(P)} modules")
out = ROOT / 'results/referee2_newdata'
pd.DataFrame(rows).to_csv(out / 'E37_feedback_cgc.csv', index=False)
P.to_csv(out / 'E37_feedback_cgc_per_module.csv', index=False)
print(f"\nfeedback genes in CGC: {sorted(FB & cgc)}")
print("DONE")
