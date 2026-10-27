# Perceptome validation suite — results (pre-registration v1, frozen 2026-09-19)

Run 2026-09-19 on the paper-3 v10 pipeline. Every number below comes from `results/`; verdicts are in `registry/verdicts.tsv`,
amendments in `registry/amendments.md`, provenance quirks in `registry/provenance_notes.md`. Scripts: `scripts/T*.py`; shared harness
`scripts/prereg_lib.py`. Inputs are SHA-256-registered in `registry/inputs_T*.json`.

## 0. Reproduction gate (T0) — passed

| check | result |
|---|---|
| gene-level HPA → bundled single_cell_scaled reference | exact (max Δ 1e-4) |
| gene-level HPA → bundled pseudobulk reference | 42/44 modules exact; 18 entries (kidney × cGAS-STING, brain × PXR/CAR) differ → older HPA build behind the bundled file (provenance note 1) |
| cancer cache = the paper's 83,819 cells | identical cohort counts and cell order; perceptome-44 scores cell-by-cell max Δ 5e-7 vs `results/E7_S_raw.npz` |
| bundled eigenspace → endpoints | E1 two-cell 23.9 %, three-cell **27.1 %**; cytotrophoblast **18/25**, rank 1; megakaryocyte 18/25 |
| rebuilt eigenspace (the pipeline every panel goes through) | Kaiser 9; E1 21.9 % / 24.7 %; cytotrophoblast 18/25 rank 2; loadings \|cos\| 0.98–1.00 vs bundled |
| E3 Pelka (bundled loadings) | +15.0 pp raw (43/46), +9.9 pp shape-only (43/46) — paper +15.5 / +10.5 (NKX3-1 handling; provenance note 2) |

Specification curve reported everywhere: loadings built on pseudobulk (paper construction, PRIMARY) vs on single_cell_scaled
(scripts/E2e construction); forced 9 PCs vs native Kaiser k; two-cell vs three-cell E1; dense vs stable-argsort rank.

## T1 — Established-pathway baselines → **PASS (modules load-bearing)**

Primary cells (pseudobulk loadings, k = 9, two-cell E1, dense rank). Incumbent: E1 21.9 % (three-cell 24.7 %), cytotrophoblast
18/25 rank 2, E3 shape-only +9.0 pp (44/46).

| panel | axes / genes | Kaiser k | E1 two-cell | E1 ratio | cyto rank (rec) | E3 raw | E3 shape-only |
|---|---|---|---|---|---|---|---|
| perceptome-44 | 44 / 206 | 9 | 21.9 % | 1.00 | 2 (18/25) | +13.9 (43/46) | +9.0 (44/46) |
| hallmark-50 | 50 / 4,384 | 6 | 2.3 % | 0.10 | 17 (1/25) | +19.2 (46/46) | +4.0 (36/46) |
| hallmark-50-trim5 | 50 / 184 | 6 | 1.8 % | 0.08 | 29 (0/25) | +9.2 (45/46) | +6.1 (42/46) |
| progeny-14 | 14 / 1,295 | 3 | 6.4 % | 0.29 | 27 (0/25) | +10.2 (41/46) | +12.9 (40/46) |
| reactome-top44 | 44 / 3,774 | 4 | 1.2 % | 0.05 | 15 (2/25) | +6.5 (38/46) | +2.5 (33/46) |
| reactome-top44-trim5 | 44 / 69 | 3 | 1.1 % | 0.05 | 34 (0/25) | +5.9 (33/46) | +2.3 (28/46) |
| kegg-signalling-46 | 46 / 2,383 | 4 | 1.4 % | 0.06 | 12 (2/25) | +14.9 (46/46) | +15.0 (44/46) |
| kegg-signalling-44 | 44 / 2,352 | 4 | 1.5 % | 0.07 | 12 (2/25) | +17.5 (46/46) | +6.1 (39/46) |

Both PASS conditions hold at once: no competitor puts the cytotrophoblast in the top 3 (best rank 9–12 across all specifications),
and every competitor has E1 < 0.5 × incumbent (max 0.29 ×). Same verdict under single_cell_scaled loadings and at native k.
Competitors converge on pancreatic acinar / gastric chief / oocyte / spermatid / lactating-breast cells — high-output secretory cells.

**Not in the verdict but decisive for the framing:** the within-patient deepening (E3) is *not* catalog-specific. KEGG-46 gives
Pelka malignant 15.5 % vs matched-normal 0.9 % (46/46 raw; +15.0 pp shape-only, 44/46), PROGENy +12.9 pp shape-only (40/46),
Hallmark +19.2 pp raw (46/46) — as large as or larger than the incumbent (+13.9 / +9.0). Cross-lineage convergence onto the
cytotrophoblast is a property of this catalog; the "deeper in the state within patients" result is a property of curated regulatory
biology generally. Paper 3 currently presents E3 as its most defensible claim; it should be presented as catalog-independent.


**Exploratory (not pre-registered), T1 addendum.** Ranking all 152 centroids by the Pelka malignant-minus-normal assignment
shift under each panel: the largest shift is toward a *trophoblast* state under every curated panel — migrating cytotrophoblast #1
under KEGG-46 (15.9 % vs 1.1 %) and Hallmark-50 (23.2 % vs 0.9 %), syncytiotrophoblast #2 under KEGG-46, extravillous / migrating
trophoblast #2–3 under PROGENy (10.4 % / 10.3 % vs 0.8 % / 2.8 %); under the incumbent: gastric progenitor #1, migrating
cytotrophoblast #2, cytotrophoblast #3. So the within-patient trophoblast-ward shift of colorectal malignant epithelium is
catalog-independent; what the incumbent adds is the *cross-lineage generality* of the trophoblast target across the 25 cancers
(competitors: 0–2 of 25 cancers put the cytotrophoblast in their top 5). Script: inline in the session log; to be added as
`T1_addendum_pelka_shift.py` if the manuscript uses it.

## T2 — The correct null for map validity → **REMOVE THE CLAIM**

Nearest-neighbour same-family agreement (152 cell types, cosine, family labels `results/fam.json`, exactly E34-M7):

| representation | NN same-family |
|---|---|
| perceptome-44 (bundled loadings, paper) | 62.5 % |
| perceptome-44 (rebuilt) | 63.8 % (scscaled loadings 65.1 %) |
| 1,000 random 206-gene panels, identical module structure | median 60.5 %, 95th pct 66.4 % (perceptome at 78th pct) |
| 1,000 random panels matched on expression decile × detection quartile | **median 61.2 %**, 95th pct 67.1 % (perceptome at 71st pct) |
| hvg-2000 → PCA-9 (ceiling) | 68.4 % (seurat HVG 71.1 %) |
| hallmark-50 → PCA-9 | 62.5 % |
| label permutation (incumbent null) | mean 14.4 %, P = 1e-4 |

Criterion "median of (3) ≥ 55 %" is met (61.2 %) → the 62 % statistic is unremarkable: any comparable gene panel through the
same recipe gives ~60 %, and HVG-PCA gives 68–71 %. Delete the nearest-neighbour statistic from the validity argument; keep only the
axis-interpretability evidence. Figure: `results/T2_nulls.png`.

## T3 — Off-manifold control → **CURRENT CLAIM STANDS**

d_min = cosine distance to the nearest of 152 normal centroids (bundled pipeline):

| population | n | median d_min | IQR | p95 |
|---|---|---|---|---|
| A malignant (25 cohorts) | 83,819 | 0.113 | 0.070–0.179 | 0.296 |
| B held-out normal single cells (skin/kidney/breast) | 16,394 | 0.076 | 0.045–0.163 | 0.282 |
| B, n_genes ≥ 1000 | 11,324 | 0.083 | 0.047–0.179 | 0.292 |
| C Pelka matched-normal epithelium | 42,013 | 0.114 | 0.078–0.165 | — |

Median(A) lies inside IQR(B) → STANDS (also vs B≥1000, with rebuilt loadings, and shape-only). Off-manifold fraction
(d_min > p99 of B): A 1.8 %, C 1.1 %. Exploratory: core-assigned malignant cells are *less* off-manifold than non-core (1.0 % vs
2.1 %; median d_min 0.118 vs 0.111) — evidence against the "joint exile" reading. Per-cohort table `results/T3_by_cohort_primary.csv`.
Not done: gastric (GSE150290) and lung (GSE131907) matched-normal arms for population C (only Pelka included).

## T4 — Module depth → **FAIL — fragile (identity, not magnitude)**

| depth | rule | E1 three-cell | E1 two-cell | cyto rank (rec) | top-5 |
|---|---|---|---|---|---|
| 3 | drop two lowest-expressed | 9.2 % | 5.1 % | 14 (1/25) | acinar, lactating breast, spermatids… |
| **5** | incumbent | **24.7 %** | 21.9 % | **2 (18/25)** | mega, cyto, gastric progenitor… |
| 10 | DoRothEA A→B, then KEGG/Reactome | 16.1 % | 14.4 % | 8 (5/25) | mega, spermatids, oocytes… |
| 15 | same | 15.1 % | 11.3 % | 9 (5/25) | chief, spermatids, mega… |
| 20 | same | 19.5 % | 15.5 % | 5 (10/25) | acinar, mega, gastric progenitor… |
| 10/15/20 | sensitivity A5: pathway members only (no TF targets) | 9.8 / 18.9 / 22.2 % | | 9 / 17 / 28 | mega + chief take over |

E1 stays ≥ 11.7 % at every depth ≥ 10, but the cytotrophoblast leaves the top 3 at every depth (both expansion rules, both loading
calibrations). Gene jackknife (206 rebuilds): 202/206 keep the cytotrophoblast at rank 1–2, E1 range 15.4–30.7 %; **four genes each
dethrone it alone** — EIF4EBP1 (mTOR; → rank 12, E1 15.4 %), FKBP5 (AR/GR/MR; rank 6), SLC2A1 (HIF; rank 8), CD36 (PPARγ; rank 7);
> 3 flags → FAIL by the pre-registered rule. The convergence *magnitude* is robust to depth; the *cytotrophoblast identity* is a
property of the 4–5-gene sampling and of a handful of genes. The manuscript must name these four genes as fragilities and hold the
identity claim more provisionally than it now does. Files: `results/T4_summary.json`, `results/T4b_gene_jackknife.csv`,
`registry/T4_expansion_catalog_{10,15,20}.json` (with the source of every added gene).

## T6 — Interpretability/fidelity frontier → **RETIRE THE PCA BASIS (adopt varimax)**

| representation | recon err (44-module) | CKA vs HVG-PCA-9 | NN-family | Interp-S | Interp-N (3 blinded raters) | E1 two / three | cyto rank | E3 shape |
|---|---|---|---|---|---|---|---|---|
| PCA (incumbent) | 0.338 | 0.661 | 62.5 % | 0.31 | 2.52 | 23.9 / 27.1 | 1 | +9.9 |
| varimax | 0.338 | 0.661 | 62.5 % | 0.52 | **3.78** | 23.9 / 27.1 | 1 | +9.9 |
| quartimax | 0.338 | 0.661 | 62.5 % | 0.53 | 3.59 | 23.9 / 27.1 | 1 | +9.9 |
| promax (oblique) | 0.338 | 0.625 | 62.5 % | 0.54 | 3.78 | 27.3 / 30.5 | 2 | +12.6 |
| ICA | 0.338 | 0.488 | 57.9 % | 0.43 | 2.85 | 19.5 / 20.8 | 2 | +5.7 |
| sparse PCA α 0.01–0.5 | 0.318–0.321 | 0.62–0.65 | 62–66 % | 0.50–0.60 | 3.67–3.89 | 11–13 / 18–19 | 1 | +11 to +15 |
| NMF-9 | 0.468 | 0.416 | 40.1 % | 0.67 | 3.04 | 5.7 / 8.5 | 34 | +6.4 |

Correctness check passed (all rotations reconstruct identically). Varimax and quartimax satisfy every RETIRE condition
(Interp-N Δ +1.26 / +1.07 ≥ 0.5, Interp-S up, E1 and cytotrophoblast rank unchanged by construction). Raters: three independent
LLM instances, blinded sheets (`results/T6_rating_sheets_blind.json`, key in registry), inter-rater Spearman ρ 0.90, Kendall's W 0.81.
The incumbent PCA basis is Pareto-dominated by varimax on (Interp-N, Fidelity-K) — state this in Limitations if the basis is kept.
Note: the bundled subspace was built on pseudobulk R, so it reconstructs the single_cell_scaled matrix worse (0.338) than a
PCA-9 refit on that matrix would (sparse PCA reaches 0.318) — a small, real fidelity cost of mixing calibrations.
T6c (guided expansion, HVG-2000): frozen 9 seeds + k_free residual PCs reach R² 0.685 at k_free = 15 vs PCA ceiling 0.713 (gap closes
from 0.12 at k_free = 0 to 0.03) — "approach the ceiling by ≈15" confirmed; a seeded non-linear AE drifts to seed CKA 0.62–0.68,
below the predicted ≥ 0.70 — that secondary prediction is refuted. Figure `results/T6_pareto.png`.

## T5 — Adversarial reconstruction of the catalog → **KILLED — curator-determined**

Three independent curations (C2, C3, C4), each given only the four operational criteria verbatim and the size target, with the
incumbent catalog, Table S9, paper 3 and the words "cytotrophoblast"/"trophoblast"/"megakaryocyte" withheld (frozen prompt
`registry/T5_prompts/curation_prompt_frozen.md`, SHA 594e1ce9…; the three catalogs are archived in `registry/T5_catalog_C*.json`).

| catalog | modules | genes | screened / dropped | Jaccard vs incumbent | E1 three-cell | cyto rank (rec) | mega rec | E3 raw | E3 shape-only |
|---|---|---|---|---|---|---|---|---|---|
| incumbent | 44 | 206 | ~75 / — | 1.000 | 24.7 % | **2 (18/25)** | 19/25 | +13.9 | +9.0 (44/46) |
| C2 | 49 | 261 | 151 / 102 | 0.283 | 8.4 % | **29 (0/25)** | 8/25 | +1.1 | +0.7 (29/46) |
| C3 | 48 | 256 | 146 / 98 | 0.305 | 16.7 % | **28 (0/25)** | 18/25 | +10.1 | +4.7 (40/46) |
| C4 | 47 | 239 | 132 / 85 | 0.309 | 10.7 % | **27 (0/25)** | 10/25 | +5.3 | +2.2 (40/46) |

0 of 3 place the cytotrophoblast in the top 5 (ranks 27–29 of 152 under both loading calibrations), mean gene Jaccard 0.299 → KILLED.
The overlap threshold is borderline and read as an aggregate (amendment A6); under a strict per-catalog reading only C2 is below 0.30,
but the pre-reg's SUGGESTIVE branch still does not apply (it needs 1 of 3 in the top 5, or overlap > 0.60), so the test is "not PASS"
on every reading.

**Between-catalog variance — the number the framework did not have:** E1 three-cell 12.0 ± 4.3 % across curations vs 24.7 % for the
incumbent; E3 shape-only +2.5 ± 2.0 pp vs +9.0; E3 raw +5.5 ± 4.5 pp vs +13.9. Representational choice moves the headline endpoint by
about a factor of two and the deepening by a factor of three.

**What the independent curations DO agree on.** Their mutual Jaccard is 0.62–0.69 — far above their 0.28–0.31 with the incumbent —
and all three put **gastric progenitor cells** first by a wide margin (21, 21, 22 of 25 cancers; the incumbent ranks them third,
13/25). C3, whose catalog is closest to the incumbent in spirit, also recovers the megakaryocyte (18/25, rank 2) and the largest E1
(16.7 %). So the four criteria do reproduce *a* cross-lineage convergence onto a high-output progenitor-like normal cell — but not
onto the trophoblast. Caveat stated plainly: the three curators are instances of one model family and share a prior, which explains
their agreement with each other; it does not explain their disagreement with the incumbent.

Each curation contains ~50 % of the incumbent's 206 genes (103, 108, 105) — the divergence is in which regulatory axes are
represented and how densely, not in an exotic gene universe. Of the four genes T4b showed to be load-bearing for the cytotrophoblast
anchor, C2/C4 contain two (FKBP5, CD36) and C3 one — another indication that the anchor depends on specific gene choices rather than
on the criteria.

## T7 — Activity-level scoring → **FAIL** (both estimators; directional prediction CONFIRMED)

Modules re-scored by regulatory output: 34 of 44 by their primary TF's DoRothEA A/B regulon (signed weights; 6–755 targets),
ERK/MAPK and PI3K/PTEN by the PROGENy footprint, and 8 without a usable regulon (AMPK, Autophagy, MR, mTOR, NPAS4, PXR/CAR, TR,
UPR-IRE1) by the mean of the catalog's documented target genes. (Corrected 2026-10-25: this paragraph read "38 ... and 6" while
naming seven cells and omitting mTOR; `registry/T7_activity_definitions.json` is the artifact and gives 34 / 2 / 8.) Definitions frozen in `registry/T7_activity_definitions.json`.
VIPER is the pre-registered estimator (primary); ULM is reported as the linear alternative.

| scoring | Kaiser | E1 two / three | cyto rank (rec) | mega rec | E3 (Pelka) | top-5 targets |
|---|---|---|---|---|---|---|
| readiness, raw (incumbent) | 9 | 21.9 / 24.7 % | **2 (18/25)** | 19/25 | +13.9 (43/46) | mega, **cyto**, gastric progenitor… |
| readiness, shape-only | 12 | 11.0 / 12.7 % | 4 (8/25) | 12/25 | +4.9 (40/46) | gastric progenitor, mega, hofbauer, cyto… |
| **activity VIPER, raw** | 8 | 33.3 / **33.3 %** | **17 (0/25)** | **25/25** | +1.5 (26/46) | **mega (25/25)**, oocytes, tuft, MEP… |
| **activity VIPER, shape-only** | 10 | 32.8 / 32.8 % | 13 (0/25) | 25/25 | **+19.5 (41/46)** | mega, oocytes, spermatids, tuft… |
| activity VIPER, scscaled loadings, raw | 9 | 43.4 / 43.5 % | 19 (0/25) | 25/25 | +11.3 (34/46) | mega (43.9 % depth), MEP, tuft… |
| activity ULM, raw | 8 | 27.6 / 32.8 % | 4 (13/25) | 20/25 | +14.5 (39/46) | oocytes (21/25), mega, spermatids, cyto… |
| activity ULM, shape-only | 9 | 25.2 / 32.5 % | 7 (5/25) | 18/25 | +8.3 (35/46) | oocytes, extravillous trophoblast, mega… |

Two of the three PASS conditions hold under both estimators — E1 three-cell 33.3 % (VIPER) / 32.8 % (ULM), floor 11.7 %; and E3
shape-only +19.5 pp (VIPER, 41/46) / +8.3 pp (ULM), threshold 5 pp, same sign. The test fails on the third, decisively under the
pre-registered estimator: the **cytotrophoblast drops to rank 17 with 0 of 25 cancers** (ULM: rank 4, 13/25). The activity-level
target is the **megakaryocyte in 25 of 25 cancers** (VIPER) or the **oocyte in 21 of 25** (ULM).

Read straight: measured by *regulatory output* rather than by *parts list*, cross-lineage convergence is stronger, within-patient
deepening is stronger in shape-only terms — and the cell it converges on is not the cytotrophoblast. Per the pre-registration, this
FAIL "is informative and must be reported: it would mean the convergence is about *equipment* for a state rather than occupancy of
it, which is a materially different claim and changes the Abelev reading."

Caveat on the VIPER megakaryocyte: the paper's own complexity control (Fig. S7) shows megakaryocyte assignment rising as detected-gene
count falls, and VIPER's rank-based enrichment is sensitive to the number of detected targets, so part of the 25/25 may be the same
artifact at a new measurement level. That does not rescue the cytotrophoblast, which is absent from the top five in every activity
specification tested (6 of 6).

**Directional prediction CONFIRMED** (one-sided, pre-registered): the shape-only fraction of E1 rises relative to raw — readiness
0.51 → activity 0.98 (VIPER) / 0.99 (ULM). Magnitude-carried dilution therefore does explain the gap between raw and shape-only
convergence: at the activity level the two coincide. The manuscript's "the convergence is largely a shared high-output signal"
caveat is a property of readiness-level measurement, not of the phenomenon.

## Kill-list / what changes in the manuscript

1. **KILLED (T5):** "the perceptome is a method" — the four operational criteria do not determine the catalog. Three independent
   curations following the same criteria produce catalogs sharing half the incumbent's genes but ranking the cytotrophoblast
   27th–29th of 152, and all three name a different target (gastric progenitor cells, 21–22/25). The cancer result must be restated
   as conditional on *this* catalog, and the framework is demoted from method to instrument-with-a-provenance. The between-catalog
   variance now exists and must be reported: E1 12.0 ± 4.3 % (incumbent 24.7 %), E3 shape-only +2.5 ± 2.0 pp (incumbent +9.0).
2. **Killed (T2):** "the map is biologically meaningful before any cancer", as evidenced by 62 % vs 14 %. Matched random panels give
   61 %; HVG-PCA gives 68–71 %. Replace with the axis-interpretability evidence, or state the corrected null.
3. **Weakened (T4):** the cytotrophoblast *identity* survives every established-pathway baseline (T1) but not module depth (rank 5–9
   at 10/15/20 genes) and not the gene jackknife (four single genes dethrone it). Name them in the manuscript as fragilities.
4. **Reframed (T1 secondary):** within-patient deepening is reproduced — and exceeded — by KEGG/PROGENy/Hallmark panels, and the
   trophoblast-ward shift in Pelka is the top shift under every curated panel. E3 is a property of curated regulatory panels
   projected onto a normal-cell gallery, not evidence for the 44-module catalog. What the catalog adds is cross-lineage generality.
5. **Retained (T3, T1 primary):** proximity is genuine — malignant cells sit inside the normal cloud, and core-assigned cells are
   *less* off-manifold than the rest. And no off-the-shelf pathway collection reproduces the cross-cancer cytotrophoblast target.
6. **Presentation (T6):** adopt the varimax rotation — strictly better interpretability at identical fidelity and identical results.
7. **Scope-restricted (T7):** the cytotrophoblast anchor is specific to readiness-level (component-abundance) measurement. Scored by
   regulatory output, convergence is *stronger* (E1 33 %, shape-only deepening +19.5 pp) but lands on the megakaryocyte (25/25, VIPER)
   or the oocyte (21/25, ULM); the cytotrophoblast is absent from the top five in all six activity specifications. The paper must say
   that its claim is about *equipment for* the state, not occupancy of it — which weakens the direct Abelev reading. The pre-registered
   directional prediction was confirmed: readiness/activity dilution is what made the raw convergence look magnitude-carried.

## Suite-level reading

The pre-registration says "a suite in which nothing is killed should be treated with suspicion, not satisfaction." Two claims were
killed (T2, T5), two were found fragile or scope-restricted (T4, T7), one was reframed (T1 secondary), two survived (T1 primary, T3),
and one presentational choice was retired (T6).

**What survives, stated exactly.** Malignant cells of 25 lineages, scored on this 44-module catalog at the readiness level and
projected into a cancer-blind normal-cell frame, (i) sit genuinely *inside* the normal manifold rather than in exile from it (T3),
(ii) converge far more strongly than any off-the-shelf pathway collection reproduces (T1, E1 ratio ≤ 0.29 ×), and (iii) sit deeper in
that region than the same patient's matched normal — though (iii) is a general property of curated regulatory panels (T1 secondary).

**What does not survive.** The cytotrophoblast as *the* named target. It is lost under module-depth expansion (T4a), under removal of
any of four single genes (T4b), under three independent curations following the same criteria (T5), and under activity-level scoring
(T7) — four independent perturbations, each moving the address to a different cell (gastric chief/acinar, gastric progenitors,
megakaryocyte, oocyte). The convergence is robust; its address is not. That is the honest headline, and it is a different paper from
the one now drafted — narrower, but one no reviewer can take apart the way this suite just did.

Amendments used: A1–A6 (`registry/amendments.md`); one in the favourable direction (A5), so the HALT rule (two) was not triggered.
