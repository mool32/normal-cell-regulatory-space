#!/usr/bin/env python3
"""The ordered SI figure list, after the v12 consolidation: 30 numbered items -> 16.

`scripts/build_v10_pdf.py` imports SI_FIGURES and injects each entry over the matching bullet in
`manuscript/paper3_v10_SI.md`, so the numbers live in exactly one place and a renumbering is a change
to this file alone. The grouping is decided and composed by `scripts/build_SI_consolidate.py`, which
also prints the `width` values below; keep the two files in step.

Citation style inside the manuscript. A composite figure's blocks carry ROMAN part labels (I, II, III)
drawn by the compositor, because several of the source panels already carry their own (A), (B), (C)
letters baked into the raster and an outer A/B/C would put a literal "A" inside the block labelled
"B". A part is cited as "Fig. S5-I", a sub-panel of that part as "Fig. S5-IB".

Each entry: (bullet label, figure file, caption, pandoc width).
Captions are single-line by construction (adjacent string literals joined by the parser); a literal
newline inside one would break the pandoc image syntax that carries it.
"""

SI_FIGURES = [

    ('Fig. S1', 'v12_figS1_axis_stability.png',
     "**Figure S1. How well each principal axis is determined, and how much biology each one carries.** "
     "**(I)** The 154 cell types are resampled with replacement 500 times; the eigenspace is rebuilt by "
     "the canonical recipe each time and each full-data axis matched to its best bootstrap partner. Bars "
     "are the mean |cos| to the full-data axis, whiskers the 5th–95th percentile: PC1 0.99, PC2 0.97, "
     "PC3 0.93, PC4 0.86, PC5 0.81, PC6 0.75, PC7 0.71, PC8 0.63, PC9 0.53. Stability falls "
     "**monotonically** with the axis, and only the first five sit above 0.80; the ninth axis is the least "
     "determined, with a 5th percentile of 0.11, which is the honest limit on reading PC9 on its own. "
     "*This part is a recomputation.* The figure shipped in earlier versions had neither a builder nor a "
     "stored result anywhere in the repository, and its caption reported a different profile (PC7–PC9 at "
     "0.40–0.54, with PC8 dipping below PC9); that non-monotonicity does not reproduce under the method "
     "stated here, and the original method cannot be checked because no code for it survives. "
     "**(II)** For each PC1–PC6, the fraction of cell position explained by the independent biological "
     "family (η², one-way ANOVA) against a permutation null (families shuffled 2,000×, 95% band). Every "
     "axis explains 31–70% — far above chance (≈ 0.09; all P < 0.001). The number inside each bar is that "
     "axis's share of total variance. `scripts/build_E0_bootstrap.py` (data "
     "`results/E0_bootstrap_stability.csv`), `scripts/E7_eta2_figure.py`.",
     '100%'),

    ('Fig. S2', 'v7_fig_E2_redundancy.png',
     "**Figure S2. Module redundancy.** Pairwise correlation of the 44 modules across the 154 normal cell "
     "types; only 6 of 946 pairs exceed |r| > 0.7 (boxed) — the three UPR arms, JAK-STAT/Type I IFN, "
     "ERK/NF-κB and BMP/Hippo. Four of the six are biological and two are partly definitional, because "
     "the paired modules share a gene (UPR-ATF6 and UPR-IRE1 share XBP1; JAK-STAT and Type I IFN share "
     "STAT1; Note S-L). Mean |r| = 0.21 — the modules are near-independent.",
     '100%'),

    ('Fig. S3', 'v7_fig_axes.png',
     "**Figure S3. Each principal axis is biologically structured (neither panel uses the modules).** For "
     "PC1–PC6, cell position on the axis (x) vs an external quality signature scored from gene-level HPA "
     "(y); the axis is monotonically ordered by that quality (immune signalling PC1 ρ = −0.89, "
     "translation PC2 +0.75, cholesterol PC5 +0.54). Biological family explains 31–70% of position on "
     "every axis (η²; all P < 10⁻⁶; Fig. S1-II). Axis conventions, previously printed inside the figure: "
     "x is the position along each principal axis (standardized, the same scale on every PC); y is an "
     "external quality signature scored from gene-level HPA without using the modules; the trend line "
     "shows the axis is monotonically ordered by that signature; colour is biological family, and "
     "families cluster along x, which is what η² measures. The leading Reactome terms at each end of "
     "every axis are listed in Note S1.",
     '100%'),

    ('Fig. S4', 'v12_figS4_geometry.png',
     "**Figure S4. Where the cancers land, shown three ways.** **(I)** A converging against a "
     "non-converging cancer, seen directly: malignant cells of B-ALL and neuroblastoma projected onto "
     "interpretable 9-D axes — cosine to the core (convergence) against cosine to the neural-program "
     "region (neurons + adrenal medulla); 20/50/80% KDE density contours over the 154 normal cells. B-ALL "
     "mass sits on the core (red stars); neuroblastoma occupies the neural-program region (toward green "
     "triangles). B-ALL is shown because it is the clearest geometric contrast, not as a specimen of "
     "cross-lineage convergence: of its 60.8% two-cell-core occupancy, 59.7 points are the megakaryocyte, "
     "i.e. largely a return to its own bone-marrow lineage, and only 1.1 the cytotrophoblast. Static "
     "positions, not trajectories; the raw PC1–PC2 projection does not resolve the convergence "
     "(Methods). **(II)** UMAP of the 154 normal cell types (position = state similarity, not distance to "
     "the core; colour = mean top-1 convergence votes from the 25 cancers). The 25 malignant states "
     "(◆) pile into one region (dashed hull); large markers = the two-cell core; the manifold's "
     "geometric centre (alveolar type-1 pneumocyte, ×) receives zero votes. Embryonic/germline cells "
     "co-localize (Fisher P = 5.6 × 10⁻⁶). **(III)** Per-cell-type convergence, column-standardized: the "
     "cytotrophoblast row is uniformly warm across solid and liquid cancers (the complexity-robust "
     "anchor); the megakaryocyte row runs hot specifically in the hematopoietic cancers (B-ALL, AML, "
     "DLBCL) that share its bone-marrow lineage.",
     '100%'),

    ('Fig. S5', 'v12_figS5_complexity.png',
     "**Figure S5. The transcriptome-complexity trap, and the control for it.** **(I)** (A) The 22-candidate "
     "screen; five datasets below the 1,000-gene floor are excluded (AML and MM replaced by GEO cohorts). "
     "(B) Per-target complexity: the core sits far above the erythrocyte artifact. **(II)** The "
     "megakaryocyte complexity control (Gudkov), across 312,220 malignant cells — the 17 of the 25 "
     "panel cohorts whose objects retain per-cell gene counts. (A) Assignment to the megakaryocyte rises as a "
     "cell's detected-gene count falls — like the erythrocyte artifact — while the cytotrophoblast rises "
     "*with* complexity. (B) Among high-complexity cells (≥ 4,000 genes) the cytotrophoblast dominates "
     "(20%) and the megakaryocyte decays to the erythrocyte-residual level (5%). (C) Causal: thinning the "
     "*same* cells collapses cytotrophoblast assignment, surges erythrocyte (positive control), and "
     "raises megakaryocyte only modestly. The cytotrophoblast convergence is complexity-robust; the "
     "megakaryocyte's is genuine in blood lineages (full complexity there) and partly a low-complexity "
     "effect in solid tumors.",
     '100%'),

    ('Fig. S6', 'v12_figS6_ranking.png',
     "**Figure S6. Ranking the candidate targets two ways.** **(I)** Target selection on breadth × depth, "
     "every candidate labelled: the SI companion to Fig. 3C, the same two axes with all fourteen "
     "recurrent candidates named rather than the seven the main-figure panel can carry. Only the core "
     "leads on both recurrence breadth and median depth; the candidate third (gastric progenitor) leads "
     "on breadth alone. **(II)** The cross-cancer recurrence ranking that gives the data-driven core: the "
     "cytotrophoblast and megakaryocyte lead with a clear gap (18→13); enriched neighbourhood in orange; "
     "post-cellular bodies (erythrocyte, platelet) excluded. Before that exclusion the erythrocyte would "
     "rank third (16/25, between the core and the gastric progenitor), but it is the documented "
     "low-complexity artifact: Table S2 carries each cell's pre-quarantine occupancy and Note S-K states "
     "the quarantine rule.",
     '100%'),

    ('Fig. S7', 'v12_figS7_catalog.png',
     "**Figure S7. How much of the result depends on this particular module catalog.** **(I)** Two nulls of "
     "matched size, 150 draws each, scored on 17 cohorts against a five-cell target set (the trophoblast "
     "states plus the megakaryocyte) — so these percentages are not in the same currency as the headline "
     "occupancy. (A) Fraction of malignant cells reaching that target under the real 44 modules (red) "
     "against random genome-wide panels (grey) and against the **same curated genes reshuffled into "
     "random modules** (blue). The real catalog clears the random-genome floor (21% against 6%, z = +3.1, "
     "P = 0.013) but is **not** separated from its own reshuffled grouping (21% against 14%, z = +1.3, "
     "P = 0.11). Above the random-genome floor that places roughly half the effect on the gene set, with "
     "the gene-set-versus-grouping contrast unresolved rather than settled. (B) The geometric "
     "concentration ratio: the real catalog concentrates **less** than either null (5.1 against 6.0 and "
     "6.8; P = 0.93 and 0.92), so concentration as such is not a property of this catalog. **(II)** "
     "Robustness to the exact 44-module choice, rebuilding the whole 25-cancer pipeline per draw "
     "(canonical per-cell scoring; eigenspace rebuilt by the correlation-PCA recipe; erythrocyte/platelet "
     "excluded as in the main analysis). (A) Dropping five of the 44 modules at random over 200 draws: "
     "the fraction of malignant cells assigned to the cytotrophoblast/megakaryocyte core is "
     "17.5 ± 2.5% (unperturbed 17.4%), never below 11.7% (≥ 6× over chance); leave-one-out 14.4–21.2%. "
     "(B) The cytotrophoblast anchor stays at the top under the same drops — rank 1 unperturbed, in the "
     "top three (by cosine to the malignant attractor) in 84% of five-module drops and 100% of "
     "single-module drops. The megakaryocyte, the secondary and partly complexity-dependent target, ranks "
     "lower in this rebuilt-eigenspace reconstruction and is not separately assessed here. **(III)** The "
     "signature is gene-level and so independent of the grouping: per-gene cytotrophoblast suppression "
     "against cancer suppression across the 206 catalog genes (202 detectable in both reference and "
     "cancer panel, over the 17 of the 25 panel cohorts whose objects carry gene-level counts; "
     "Spearman ρ = +0.56, with the real ρ above all 2,000 gene-label permutations — "
     "P < 5 × 10⁻⁴, the floor 2,000 draws can reach; the agreement is rank-wise and not linear, "
     "Pearson r = +0.10), holding for genes from internally co-expressed modules (ρ = +0.60) and from "
     "the rest (ρ = +0.54). Read together: the "
     "convergence needs curated biological genes, does not need this particular partition of them into "
     "modules, and does not need this particular set of 44. `scripts/build_E2c_target.py` (data "
     "`results/E2c_target_null.csv`), `scripts/build_E2e_robustness.py`, "
     "`scripts/build_figS10_genelevel.py`.",
     '91%'),

    ('Fig. S8', 'v7_fig_duetprofile.png',
     "**Figure S8. The core against the normal-cell distribution, module by module.** For each of the "
     "twelve most target-defining modules, the distribution of all 154 normal cells (grey) with the core "
     "marked; the core sits in the bottom ~1–16th percentile of every one — an outlier outside the normal "
     "cloud, not merely at its low edge.",
     '100%'),

    ('Fig. S9', 'v7_figS28_reprogramming.png',
     "**Figure S9. The named target is detected in an independent system with no cancer in it.** Mouse "
     "fibroblast-to-iPSC reprogramming (Schiebinger et al. 2019; 236,285 cells, 39 timepoints), scored on "
     "the same 44 modules through orthologs (43 of 44 recovered) and projected through the same fixed "
     "coordinate. The fraction of cells assigned to the cytotrophoblast collapses to 1.8% by day 4, then "
     "rises about ninefold, to between 11% and 20% from day 14 onward (median 16.7% over those nine timepoints; Spearman ρ = +0.76 with day, "
     "P = 3 × 10⁻⁸) — the window in which reprogramming is known to generate extraembryonic and "
     "trophoblast-like off-target fates. The coordinate therefore reports trophoblast-ness where it is "
     "independently expected, in another species and another laboratory, with no cancer in the analysis. "
     "We do not use this figure to establish that claim. The starting fibroblast timepoints (shaded) "
     "score anomalously high and we cannot explain it; if unreprogrammed fibroblasts already score high "
     "on cytotrophoblast assignment, the coordinate may be reading proliferative high-output character "
     "rather than trophoblast identity, and the collapse-then-rise could be a composition effect. The raw "
     "metric is also flatter than the shape-only one, and the object carries no cell-type labels, so we "
     "cannot confirm these are the cells the original study annotated as trophoblast-like. Suggestive "
     "only. `scripts/E28_reprogramming_validation.py`.",
     '88%'),

    ('Fig. S10', 'v7_figS19_speccurve.png',
     "**Figure S10. Specification curve of the core convergence (multiverse).** Core convergence recomputed "
     "under all 36 crossings of four analytic choices — reference calibration (single-cell against "
     "pseudobulk), magnitude handling (raw / within-cell-centered / per-cell rank), core definition "
     "(three-cell / two-cell / cytotrophoblast-only), and erythroid handling (excluded / included). "
     "(A) The specifications sorted by outcome, each against its own chance baseline (ticks). Under the "
     "single-cell calibration the single-cell tumor data require, the convergence exceeds chance in every "
     "one of the 18 specifications (median 7.2%, worst case 2.9× chance); the only specifications that "
     "collapse to the chance baseline all use the pseudobulk calibration, which independently degenerates "
     "(Note S-B). (B) Median convergence by analytic choice — the estimate is stable to core definition "
     "and erythroid handling and is carried by module magnitude (the rank transform retains a smaller "
     "residual), while the single-cell calibration is the load-bearing preprocessing choice. "
     "`scripts/E7_specification_curve.py`.",
     '100%'),

    ('Fig. S11', 'v12_figS11_withinpat.png',
     "**Figure S11. The within-patient increase replicates across organs and survives magnitude removal.** "
     "**(I)** Paired gastric cohort (GSE150290): twenty-three patients (Kim et al., 2022), a "
     "low-baseline organ analysed only after the account was fixed; its tumour samples are also "
     "the panel\u2019s gastric cohort (Note S-U). Per patient, matched adjacent-normal (Pat##-A) and cancer (Pat##-B) "
     "epithelium (EPCAM+/CD45−) were projected and scored; lines join the two arms of each of the 22 "
     "patients meeting the ≥ 20-epithelial-cell floor (red = deeper in cancer, grey = not; black = mean). "
     "Malignant epithelium is deeper than matched normal in 18 of 22 patients (mean 8.8% → 14.2%; "
     "cytotrophoblast-specific 7.0% → 11.9%; paired Wilcoxon *P* = 0.001). The one grey line falling to "
     "zero is Pat01, whose cancer arm yielded no core-assigned epithelium. **(II)** High-baseline lung "
     "deepens too, and the earlier \"no deepening\" there was a dilution artifact. "
     "(A) Lung cohort (GSE131907; Kim et al. 2020), five of whose ten paired patients are also among "
     "the panel\u2019s smallest cohort (Note S-U); within-patient normal (nLung) against tumor "
     "(tLung) epithelium using the authors' own tumor-subtype annotations: malignant epithelium is deeper "
     "than matched normal in 10 of 10 patients (mean 23% → 39%; cytotrophoblast-specific 13% → 29%; "
     "paired Wilcoxon *P* = 0.002). (B) Why the pan-tissue pipeline had scored lung as non-deepening: on "
     "the CELLxGENE Census lung data it used, the broad epithelial filter admitted a \"tumor arm\" that is "
     "only ~9% annotated-malignant cells and ~91% normal-adjacent alveolar epithelium, so its core "
     "fraction (~13%) collapses to the level of a pool of healthy Census donors (~12%); the clean "
     "malignant cells sit at ~29%. Three normal-lung values appear here and they are not interchangeable: "
     "23.4% is the matched normal arm of the GSE131907 cohort, 26.7% is the patient-matched "
     "adjacent normal of the Census cohort, and 11.7% is the healthy-donor pool. **(III)** The deepening "
     "is a genuine module-shape shift rather than transcriptional brightness, and is about constant "
     "across baseline. (A) Within-patient deepening (malignant minus "
     "matched-normal core fraction) recomputed with the overall-magnitude component of each cell removed "
     "(per-cell centering across the 44 modules), in three cohorts: colorectal +10.5 points (44/46; "
     "paired Wilcoxon *P* = 10⁻¹³), corrected lung +12.1 (9/10; *P* = 0.004, sign-test *P* = 0.011), and "
     "gastric +3.4 (15/22; *P* = 0.04 by Wilcoxon but borderline by sign test, *P* = 0.07) — "
     "63–77% of the raw deepening retained. Patient-matched malignant epithelium therefore sits deeper in "
     "the *shape-only* core, not merely brighter; the colorectal and lung results are significant under "
     "either test, the smaller gastric cohort suggestive. (B) Additive deepening (points) against "
     "normal-tissue baseline across five cohorts: no dependence (*r* = +0.4, *P* = 0.5; mean about +10 "
     "points) — the apparent \"more room where the baseline is lower\" holds only in fold-change (a "
     "small-denominator artifact) and is retired. `scripts/E10_heldout_gastric.py`, "
     "`scripts/E12_highbaseline_lung.py`, `scripts/E13_lung_dilution_proof.py`, "
     "`scripts/E17_deepening_shapeonly.py`, `scripts/E17b_shapeonly_gastric_lung.py`.",
     '93%'),

    ('Fig. S12', 'v7_figS24_fibrosis.png',
     "**Figure S12. The deepening is specific to malignancy, not to injury or invasive remodelling "
     "(replacement negative control).** Idiopathic pulmonary fibrosis (Habermann et al. 2020; 33,321 "
     "epithelial cells, 22 donors) is the strictest non-malignant comparator available in the same organ "
     "as the corrected lung result: chronically injured, ECM-remodelling epithelium containing the "
     "aberrant-basaloid and transitional-AT2 populations often described as cancer-like. (A) Deepening "
     "relative to each condition's own control epithelium: fibrosis moves +3.8 points on the raw core and "
     "+1.3 on the shape-only core (neither significant, *P* = 0.12 and 0.21, 22 donors), while malignancy "
     "in the same organ moves +15.7 on the raw core (10 of 10 patients, *P* = 0.002) and +12.1 on the "
     "shape-only core (9 of 10, *P* = 0.004; Fig. S11-IIA and S11-IIIA). (B) Within the fibrotic lung the "
     "ordering is biologically sensible — transitional-AT2 (18.4%) and proliferating epithelium (16.2%) "
     "score highest, ciliated cells lowest (5.2%) — so the state tracks invasive, high-output phenotypes "
     "in general, but only malignancy carries a tissue into it. "
     "`scripts/E24_regeneration_specificity.py`.",
     '100%'),

    ('Fig. S13', 'v12_figS13_progression.png',
     "**Figure S13. Before invasion and after dissemination.** **(I)** Pre-invasive epithelium already "
     "occupies the convergent region (out-of-sample test of the prediction). Within-patient colorectal "
     "progression (GSE161277; Zheng et al. 2022): epithelial cells (EPCAM+/CD45−) from matched normal, "
     "pre-invasive adenoma and carcinoma of the same patient, projected and scored for core convergence. "
     "(A) Full core (cytotrophoblast + megakaryocyte) and (B) the cytotrophoblast-specific fraction, per "
     "patient (thin) and mean (red), for the three patients with the full series. Pre-invasive adenoma "
     "epithelium is already elevated above matched normal in all three patients — mean core 4.1% "
     "(normal) → 11.2% (adenoma) → 16.8% (carcinoma), adenoma reaching two-thirds of the carcinoma level "
     "— consistent with the prediction that the convergent state is entered before invasion. One patient "
     "(Patient 3) is non-monotone at carcinoma (215 epithelial cells, low power). A simple EPCAM+/CD45− "
     "marker gate stands in for author malignancy calls; a larger pre-invasive cohort is the natural "
     "confirmation. **(II)** The same prediction across two further organs: directionally supported in "
     "breast, refuted in esophagus. (A) Breast (GSE195861): epithelial core fraction by stage — normal "
     "breast 7.7% (a single sample), ductal carcinoma *in situ* 11.3% (7 samples), invasive ductal "
     "carcinoma 21.8% (6), lymph-node metastasis 13.4% (6). The pre-invasive lesion sits between normal "
     "and invasive, at about half the invasive level (DCIS against IDC *P* = 0.07); donors are unpaired "
     "and the normal arm is a single sample, so this is directional support, not a test. (B) Esophagus "
     "(GSE291080), within patient: Barrett metaplasia against the same donor's own normal squamous "
     "esophagus — the metaplasia sits *below* normal in four of four patients (2.1% against 7.7%), and "
     "below gastric epithelium in three of three. The prediction fails here. Two readings remain open: "
     "that it is confined to neoplastic precursors (adenoma and in-situ carcinoma are clonal neoplasms; "
     "Barrett metaplasia is a lineage switch), or that it is simply wrong in esophagus. The available "
     "samples are largely non-dysplastic, so esophageal dysplasia — the properly pre-invasive lesion — "
     "is untested. **(III)** The state accompanies invasion and is relinquished on colonization "
     "(suggestive). Within-patient primary tumour against matched lymph-node metastasis (breast, "
     "GSE195861; six pairs): the metastasis sits *lower* in the convergent state than the same patient's "
     "primary — means of 21.8% → 13.4% on the raw core (lower in six of six; paired *P* = 0.031, the "
     "floor a two-sided Wilcoxon can reach at n = 6) and 14.3% → 9.6% shape-only (five of six; *P* = 0.063). "
     "The single evaluable liver pair runs the same way. Of the two portal-vein tumour thrombi available — "
     "tumour actively occupying a vessel — one sits well above its primary (49.5% → 69.9% raw, "
     "31.4% → 46.3% shape-only) while the other does not move on the raw core (54.1% → 54.0%) and falls "
     "with magnitude removed (49.4% → 42.0%); one case is not a direction. Nodal (not distant) metastases and small numbers: reported as suggestive. "
     "`scripts/E9_preinvasive_colon.py`, `scripts/E22_preinvasive_barrett_dcis.py`, "
     "`scripts/E25_primary_vs_metastasis.py`.",
     '91%'),

    ('Fig. S14', 'v7_figS26_dynamics.png',
     "**Figure S14. The coordinate moves during targeted therapy, and the movement is not a "
     "transcriptome-complexity artifact.** Human PC9 lung-adenocarcinoma cells through eleven days of "
     "erlotinib (Aissa et al. 2021); Drop-seq, median 161–556 detected genes per cell. **(A)** Shape-only "
     "core fraction per timepoint. The grey dashed curve is the analysis as first run, at a floor of "
     "≥ 500 detected genes; the three blue curves instead restrict *every* timepoint to one common "
     "absolute detected-gene window, which is the like-with-like comparison. The level is largely set by "
     "that threshold — day 2 reads 18%, 43%, 48%, 53% or 62% at floors of 0, 200, 300, 500 or 1,000 "
     "genes — so no absolute core fraction from this dataset is interpretable without its floor, and we "
     "report none. The shape is not: the rise to day 2 and the day-4 drug-tolerant persister trough "
     "appear in every window above 1,000 genes (n = 67–368 cells per point in the 1,500–6,000 window). "
     "**(B)** Why that control is necessary. Within every timepoint, core assignment tracks detected-gene "
     "count almost monotonically (Spearman +0.77 to +1.00; at day 2 the lowest decile is 2% and the "
     "highest 68%) — the same complexity dependence for which the erythrocyte is quarantined in quality "
     "control (Fig. 3D, Note S-Q), here operating inside the time course itself. Limits, stated plainly: "
     "survivorship thins the later timepoints, different cells are sampled at each one, and no lineage "
     "tracing is available, so movement of cells cannot be separated from selection among them; and the "
     "endpoint is **not** stably above baseline — days 9–11 sit 16–21 points above it in the 1,500–6,000 "
     "window, return to it in the 1,000–3,000 window, and fall below it when every timepoint is resampled "
     "to the untreated cells' own complexity distribution. Read as evidence that the coordinate moves "
     "under a real perturbation, not that the state is entered and held. `scripts/E26_dynamics_persisters.py`, "
     "`scripts/E26b_qc_floor_sweep.py`.",
     '100%'),

    ('Fig. S15', 'v7_fig_detuning.png',
     "**Figure S15. What malignancy adds within the shared state (within-patient, colorectal).** Malignant "
     "against matched-normal core cells compared per module (Pelka; detection-depth controlled). The "
     "malignant occupant runs retained stress/integrity sensors (HIF, NRF2, HSF1, UPR-PERK, autophagy) "
     "and anabolic programs hotter; the p53 module rises through downstream stress/damage targets "
     "(p21/BAX/PUMA/NOXA), not *TP53* itself (transcript unchanged) — a replication-stress footprint, not "
     "retained p53 function. The separation is mostly off the target-defining axis (ρ = +0.36): a "
     "qualitative detuning toward stress-and-growth. Open circles in (A) are the raw difference "
     "before the detection-depth control (balanced 2,640 against 2,496 detected genes, from 4,351 "
     "against 2,395 raw).",
     '100%'),
]
