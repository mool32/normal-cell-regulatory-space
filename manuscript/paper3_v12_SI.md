# Supporting Information — Human malignancies show increased occupancy of a shared region of normal-cell regulatory space

**Theodor Spiro** · Vaika Inc. · Archived at Zenodo: https://doi.org/10.5281/zenodo.20542130

This SI contains extended methods (Notes S1–S5), the control notes S-A to S-M, two notes that are primary evidence about the method — **Note S-N** (the panel evaluated without cancer data) and **Note S-O** (independent reconstructions of the catalog, which is where the catalog-dependence of the anchor is established), plus **Note S-Q** (the complexity control run across the whole convergence region and the gene-symbol recovery audit), **Note S-R** (the pre-registration and per-contrast detail of the five *in vitro* transformation series of Fig. 6), **Note S-T** (the pre-registered activity-scoring test, which is where the measurement-level dependence of the anchor is established), and **Note S-U** (what the twelve CELLxGENE Census cohorts of the panel actually resolve to, with Table S-Q). Notes S-E and S-S were withdrawn and their letters are left unused, and Tables S-P and S-Q carry letters for the same reason. There are fifteen supplementary figures (S1–S15), each composed of one to three parts labelled I, II and III, and twelve supplementary tables (S1–S9, S9b, S-P and S-Q), plus Table 1 in the main text. Tables S1, S2, S3, S4, S6 and S7 are provided as attached CSV files (`tableS1_datasets.csv`, `tableS2_ranking154.csv`, `tableS3_invariance.csv`, `tableS4_modules.csv`, `tableS6_complexity.csv`, `tableS7_core_signature.csv`); all others are inline.

---

## Note S1 — Full methods

**Modules.** A regulatory module is defined by four operational criteria: a sensor input; a transcription-factor output with documented targets; a negative-feedback target enabling cyclical operation; and receptor-field specificity. Screening ~75 candidate signaling and metabolic axes (KEGG, Reactome, primary literature) against these criteria yielded 44 modules. Module activity is scored as the mean log-normalized expression of component genes.

**Eigenspace.** The 154 × 44 Human Protein Atlas single-cell matrix (Uhlén et al., 2015; scRNA-seq data, not protein measurements) was z-scored per module and decomposed by PCA. Nine components were retained by the Kaiser criterion (73.8% of variance). Bootstrap resampling with replacement (n = 154, 500 iterations) gives eigenvector stability that falls monotonically with the axis: 0.75–0.99 for PC1–PC6 and 0.53–0.71 for PC7–PC9 (Fig. S1-I). The 100-iteration field stored with the reference object (0.63–0.99 and 0.40–0.54, with PC8 below PC9) does not reproduce under this recipe and is withdrawn; no code for it survives.

**Two HPA calibrations.** The eigenspace is built on the *pseudobulk* HPA reference. For per-cell tumor projection, the reference is re-derived under the *single_cell_scaled* calibration — the same `normalize_total(1e4)+log1p` preprocessing applied to the tumor data — so per-cell z-scores are commensurate; projecting under the pseudobulk calibration collapses tumor cells to low-complexity targets (scale mismatch). All cancer projections use single_cell_scaled.

**Cancer projection.** Per cancer: load malignant cells; remap ENSEMBL→symbol where needed; uppercase and de-duplicate gene names; quality-filter to n_genes ≥ 1000; `normalize_total(1e4)+log1p` if max > 100; score 44 modules (`mean_raw`); z-score against the single_cell_scaled reference; project through the 44×9 loading matrix; record the top-1 nearest of 154 cell-type centroids (cosine).

**Convergence projection (Fig. S4-I).** Malignant cells of two contrasting cancers (B-ALL, neuroblastoma; ≤2,000 cells each, n_genes ≥ 1000) are projected onto two interpretable 9-D directions — cosine to the core centroid (convergence axis) and cosine to a neural-program centroid (brain excitatory/inhibitory neurons and adrenal medulla cells — the neuronal/sympathoadrenal program neuroblastoma derives from; neural-program axis) — with 20/50/80% KDE density contours over the 154 normal cell types. B-ALL median core-cosine 0.86 (neural-program 0.07); neuroblastoma 0.59 (neural-program 0.47). Positions are static snapshots, not trajectories. These axes are used because the raw PC1–PC2 projection does not resolve the convergence (below).

**Convergence map (Fig. S4-II).** The 154 normal cell types are colored by the mean across 25 cancers of the fraction of malignant cells assigned to each (top-1), painted on the PC1–PC2 layout. Color carries the convergence statistic; position is eigenspace location, shown for context only. A raw projection of per-cell tumor coordinates into PC1–PC2 and re-derivation of nearest-neighbor assignment recovers only 8.2% core-assignment for B-ALL against 60.8% in the full 9-D space, and 6.8% against 23.9% pooled over all 25 cancers; any two cancers' raw clouds overlap across the whole manifold — convergence is a high-dimensional nearest-neighbor property, not a 2-D location. (Both arms of each pair are recomputed on the current 83,819-cell panel under the single quarantine convention; `scripts/E42_2d_vs_9d.py`.)

**Target identification (non-circular).** Cell types are ranked by cross-cancer *recurrence* (number of the 25 cancers in which the cell is among the five most frequent top-1 assignments), not by proximity to a pre-chosen centroid. The two-cell core (megakaryocyte, cytotrophoblast) leads this ranking with a clear margin (cytotrophoblast 18/25, megakaryocyte 18/25) and is independently the only pair leading the *median per-cancer fraction* (depth) ranking. No third cell leads on both metrics (Table S2; Fig. S6-I). The erythrocyte and platelet — post-cellular bodies — are excluded (Note S3); the erythrocyte is also the documented low-complexity artifact (Fig. S5-I/S4).

**Enrichment.** Reported as fold-enrichment over the k/152 chance baseline (152 candidates after the erythrocyte/platelet quarantine; Note S-K) across three principled boundaries (top-2 core, 18.2×; top-4 neighbourhood, 15.9×; top-14 broader region, 8.4×) to show boundary-independence: boundaries are defined by rank in the cross-cancer recurrence ranking, and all three are computed on the same 83,819-cell panel: the top-2 core captures 23.9% of malignant cells (chance 2/152 = 1.32%, 18.2×), the top-4 neighbourhood 41.9% (chance 2.63%, 15.9×) and the top-14 broader region 76.9% (chance 9.21%, 8.4×). Region coherence is the mean pairwise eigenspace cosine among the convergent cells (0.77 vs −0.01 for random pairs).

**Progression.** Within-cancer, patient-level aggregation with n_genes control; HCC clinical stage (Lu 2022) and colorectal grade/T-stage/nodal status (Pelka 2021); Table S5.

**Convergence contraction (Fig. 2).** Each cancer's malignant 9-D centroid (mean over QC-passing malignant cells) is paired with its normal cell of origin in the 154 reference (e.g. HCC→hepatocyte, B-ALL→B-cell, osteosarcoma→fibroblast, seminoma→spermatogonia; full map in Table S8). The contraction is the ratio of the mean pairwise cosine distance among the 25 origins to that among the 25 malignant centroids. Significance is **not** taken against 2,000 random draws of 25 normal cell types (a Stayton-type convergence test): that null is **superseded**, because it ignores the averaging itself, and the fold-contraction it produced is withdrawn. The test that stands is the matched-averaging null of coherent normal micro-populations, which is what the main text and Fig. 2 report — see **Note S-G**, which is the Methods statement for this figure. A replicator should read Note S-G before reproducing this panel. The MDS basin (Fig. S4-II) is a metric-MDS embedding of the origins + malignant centroids on their full 9-D distance matrix (all axes, not leading PCs).

**Panel expansion to 25.** Six datasets were added to span origin classes absent from the original nineteen: small-cell lung (neuroendocrine) and gastric (Census, malignant = epithelial compartment); osteosarcoma (GSE162454) and Ewing sarcoma (GSE130024, PDX) for mesenchymal origin; embryonal germ-cell tumor (GSE228501) and seminoma (GSE197778) for germ-cell origin. For GEO datasets without author malignant annotations, malignant cells were identified by Leiden clustering and lineage markers (osteoblastic SATB2/RUNX2/COL1A1; Ewing CD99/NKX2-2; germ-cell POU5F1/SALL4/KIT) versus immune (PTPRC). Per-dataset detail in Table S1.

**Geometry / coordinate invariance.** Convergent-region compactness was tested in three coordinate systems — PCA (9-D), varimax-rotated (6-D), and diffusion-map (8-D) — against a null of random equal-size cell subsets (Note S4, Table S3).

**Disease negative control and within-patient matched test (Fig. 4; Figs. S11-II, S12).** The identical pipeline (QC n_genes ≥ 1000 → normalize_total(1e4)+log1p → 44-module mean_raw scoring → z vs single_cell_scaled HPA → 9-PC projection → top-1 nearest of 154) was first validated to reproduce the cancer core fractions (Gastric 16.5% vs panel 16.1%; AML 20.4% vs 21.3%), then applied to four chronic-disease atlases with matched-normal arms: rheumatoid arthritis (excluded — median 286 genes/cell, below the QC floor), idiopathic pulmonary fibrosis (Habermann 2020), diabetic kidney disease, and Alzheimer's (Grubman 2019). Core fraction was computed per arm with the erythrocyte low-complexity default reported separately; tissue baselines: kidney 0.5–1.0%, neural 7.6%, lung epithelium 16.1% (Control) / 19.9% (IPF) on the raw core and 7.6% / 8.9% shape-only (22 donors; Fig. S12), breast epithelium 25% (Wu 2021). For the within-patient test, malignant vs patient-matched normal epithelium of the *same* tissue was compared per patient (≥ 20 cells/arm) by paired Wilcoxon plus a sign test, across every tissue with paired tumour-and-adjacent-normal single-cell data. Low-baseline tissues: Pelka colorectal (46 patients, malignant median 18.9% vs normal 4.7%, 43/46, *P* = 9 × 10⁻¹³; separate normal samples, clean); an independent colorectal cohort (Human Tumor Atlas via Census, 23/25, one-sided *P* = 2 × 10⁻⁶); kidney (clear-cell RCC, Census, 6/7, median 11.5% vs 5.3%, *P* = 0.055); and liver (Lu HCC — 8 of its 10 patients have ≥ 20 hepatocytes in *both* the tumour and adjacent-normal arms, distinct from the 10-patient stage-trend analysis in the main text which only requires tumour-side stage metadata — 5/8, 19.4% vs 10.9%, *P* = 0.074). Pooled across these 86 patients, 77 are deeper (sign-test *P* = 7 × 10⁻¹⁵); a mixed-effects model with a dataset random intercept, which credits the between-dataset clustering, gives +13.0 points (SE 1.87, *z* = 6.95, one-sided *P* = 2 × 10⁻¹²).

**Retracted negative control.** An earlier version of this work reported lung adenocarcinoma (Census, 2/8, median 14% for the “tumour” arm against 27% for the patient-matched adjacent-normal arm of the same cohort) as a high-baseline tissue showing no deepening. Three different normal lung values appear in this paper and they are three different things: 26.7% is that cohort's *patient-matched adjacent normal*, 11.7% is a *pool of healthy Census donors* (Fig. S11-IIB), and 23.4% is the *matched normal arm of the GSE131907 cohort* (Fig. S11-IIA). That reading is withdrawn: the Census "tumour" arm assembled by the broad epithelial filter is only ~9% annotated-malignant cells and ~91% normal-adjacent alveolar epithelium, so it was diluted to the normal baseline. With cleanly identified malignant cells (GSE131907) lung deepens like every other tissue tested — 23.4% → 39.0% raw in 10 of 10 patients (*P* = 0.002) and 18.2% → 30.3% shape-only in 9 of 10 (*P* = 0.004) (Fig. S11-II, Note S-J). No baseline-gating claim is made anywhere in this paper. Wu BRCA is uninformative (its "normal epithelial" cells are intratumoral, not separate adjacent-normal), and pooled cross-dataset ratios are vulnerable to Simpson's paradox (Lu pooled 4.6× collapses per patient); only the within-patient paired test is reported as causal.

**Prostate — an unresolved tissue, reported as such.** In a 15-patient paired tumour/adjacent-benign atlas (CELLxGENE) prostate shows no within-patient deepening (7 of 15 patients). We previously set it aside as a saturated high-baseline tissue (~15%); that justification rested on the baseline-gating hypothesis we now retract, so we do not use it. Re-analysing an independent paired prostate cohort (Song 2022, GSE176031, four patients with paired benign and tumour epithelium) does not resolve it either: there the benign prostate baseline is *low* (2.1%, near the chance level) and the tumour arm reaches only 3.8% (paired *P* = 0.38, n.s.). The two prostate datasets therefore disagree about the baseline itself, and neither shows deepening. Prostate is counted as a tested tissue with a negative result, and the claim in the main text is stated accordingly as deepening in every tissue tested *except prostate*, whose cohorts are mutually inconsistent.

---

## Note S2 — NPAS4 leave-one-out

The single neuron-restricted module (NPAS4) does not shape the cross-tissue eigenspace. Rebuilding on 43 modules (pseudobulk basis, matching the canonical reference) reproduces all nine components with per-PC eigenvector cosine ≥ 0.95 (mean 0.99 over PC1–PC6), so all cross-tissue statistics transfer unchanged to the 44-module catalog. NPAS4 is retained for compatibility with neural single-cell datasets; its module score is near zero in non-neural cell types and does not affect cross-tissue analyses.

---

## Note S3 — Low-complexity artifact and quality control

Per-cell module scores depend on transcriptome complexity: cells with few detected genes carry degenerate profiles that map preferentially to low-complexity normal targets (erythrocytes median 890 genes, late spermatids 763 — both *below* the QC floor; Table S6). Screening 22 candidate cancer datasets by median genes/cell, five fall below 1,000 (AML 137, OralSCC 639, MM 665, Follicular 831, Gastric 912; Fig. S5-IA). These are excluded; AML and MM are replaced by GEO cohorts that pass the floor (van Galen 2019, Seq-Well, median ~1,050; Ledergor 2018, InDrops), the other three dropped. The final panel is 25 — nineteen passing the initial screen plus six added (two mesenchymal sarcomas, two germ-cell tumors, one neuroendocrine, one new epithelial tissue) to span all major origin classes (Note S1). The gastric dataset excluded here (median 912) is the original screen cohort; it is **not** the gastric in the final panel — that is a different CELLxGENE Census stomach dataset (the "one new epithelial tissue" added), whose twenty-three tumour donors are however the same twenty-three patients as the paired gastric test of §5, which is why that test is no longer described as held out (Note S-U), with median 1,169 genes/cell, above the floor (54% of its cells pass the per-cell n_genes ≥ 1000 filter). After QC, convergence-fraction shows no association with dataset cell count (Spearman ρ = +0.14, P = 0.60).

The erythrocyte family is excluded from the recurrence ranking and the main-panel figures specifically: before exclusion it ranks third (16 of 25), filling the gap between the core and the rest (Table S2 carries each cell's pre-quarantine occupancy; Note S-K). The core is not a residual artifact — megakaryocyte- and cytotrophoblast-assigned cells carry median 1,396 and 2,430 detected genes, well above the floor (Fig. S5-IB).

**Platform / label-quality.** Convergence-fraction ordering is sensitive to malignant-label quality (three cancers carry heuristic or Census labels; a colon dataset's author-validated relabeling moved it from 13% to 36%) and, despite uniform QC, to residual platform effects (one InDrops dataset; the rest 10x / Census / Seq-Well). Per-dataset provenance in Table S1.

**Marker-annotated additions do not drive the result.** Four of the six added cancers carry marker-based rather than author malignant labels (both sarcomas, both germ-cell tumors). Dropping all four and recomputing on the remaining 21 leaves the central numbers essentially unchanged: the malignant compartment's mean pairwise cosine distance is 0.166 across the full 25 and 0.174 across the 21 (the drop-4 value; not to be confused with the 0.174 obtained when only 25 cells per cancer are averaged, Fig. 2B) (both far below the 0.99 separating the tissues of origin), and the core holds rank 2 and 3 of 154 in the full-9-D cosine to the malignant attractor in both panels. The marker-based cancers neither create the convergence nor inflate the core's standing.

---

## Note S4 — Coordinate invariance attaches to the region, not the bare core

The convergent **region** (the four-cell majority-enriched neighborhood) is an unusually tight subset of the 154 cell types, and that tightness is coordinate-invariant: it is significant in both the linear PCA basis (percentile 0.005) and the varimax rotation (0.001), and directional in the diffusion-map embedding (0.072), confirming the convergence is geometry, not an artifact of the linear PCA basis (Table S3). The two-cell **core** alone is too small for the compactness statistic (a single pairwise distance; percentile 0.07–0.24, not significant) — identity attaches to the core (which two cells), invariance to the region (whether the neighborhood is tight). The two carriers are kept distinct.

**Table S3. Coordinate-invariance compactness, region versus bare core, in three bases** (`tableS3_invariance.csv`).

| Set | Coordinates | Percentile | z |
|---|---|---:|---:|
| neighborhood (4) | PCA (9-D) | **0.005** | −2.09 |
| neighborhood (4) | varimax (6-D) | **0.001** | −2.38 |
| neighborhood (4) | diffusion (8-D) | 0.072 | −0.83 |
| core (2) | PCA (9-D) | 0.113 (ns) | −1.13 |
| core (2) | varimax (6-D) | 0.073 (ns) | −1.32 |
| core (2) | diffusion (8-D) | 0.243 (ns) | −0.53 |

The diffusion embedding is no tighter than PCA, consistent with the manifold being only weakly curved (PCA adequate). Compactness in the pseudobulk calibration is diffuse (percentile ~0.46); the tightness is specific to the single_cell_scaled calibration used for tumor projection, as expected from the calibration's construction. This is not calibration-shopping: Note S-H shows the pseudobulk alternative degenerates independently, on normal single cells with no cancer involved, collapsing 99.6% of them onto erythroid targets — so the diffuseness here is a property of the wrong calibration for single-cell data, not evidence that the tightness was selected for.

---

## Note S5 — Metabolic-enzyme stratification (secondary)

On the normal side, convergent-region cells are glycolytic-oxidative-active and distinct from NAMPT-pathway inflammatory cells; on the cancer side the corresponding stratification is null. Reported as a secondary observation only; it does not bear on the convergence result.

---

## Supplementary figures

Sixteen numbered figures. Several are composites: their blocks carry Roman part labels (I, II, III),
and a block's own sub-panels keep their letters, so a part is cited as Fig. S5-I and a sub-panel of
it as Fig. S5-IB. Grouping and composition: `scripts/build_SI_consolidate.py`; the ordered list with
captions: `scripts/si_figures_v12.py`.



![**Figure S1. How well each principal axis is determined, and how much biology each one carries.** **(I)** The 154 cell types are resampled with replacement 500 times; the eigenspace is rebuilt by the canonical recipe each time and each full-data axis matched to its best bootstrap partner. Bars are the mean |cos| to the full-data axis, whiskers the 5th–95th percentile: PC1 0.99, PC2 0.97, PC3 0.93, PC4 0.86, PC5 0.81, PC6 0.75, PC7 0.71, PC8 0.63, PC9 0.53. Stability falls **monotonically** with the axis, and only the first five sit above 0.80; the ninth axis is the least determined, with a 5th percentile of 0.11, which is the honest limit on reading PC9 on its own. *This part is a recomputation.* The figure shipped in earlier versions had neither a builder nor a stored result anywhere in the repository, and its caption reported a different profile (PC7–PC9 at 0.40–0.54, with PC8 dipping below PC9); that non-monotonicity does not reproduce under the method stated here, and the original method cannot be checked because no code for it survives. **(II)** For each PC1–PC6, the fraction of cell position explained by the independent biological family (η², one-way ANOVA) against a permutation null (families shuffled 2,000×, 95% band). Every axis explains 31–70% — far above chance (≈ 0.09; all P < 0.001). The number inside each bar is that axis's share of total variance. `scripts/build_E0_bootstrap.py` (data `results/E0_bootstrap_stability.csv`), `scripts/E7_eta2_figure.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS1_axis_stability.png){width=100%}




![**Figure S2. Module redundancy.** Pairwise correlation of the 44 modules across the 154 normal cell types; only 6 of 946 pairs exceed |r| > 0.7 (boxed) — the three UPR arms, JAK-STAT/Type I IFN, ERK/NF-κB and BMP/Hippo. Four of the six are biological and two are partly definitional, because the paired modules share a gene (UPR-ATF6 and UPR-IRE1 share XBP1; JAK-STAT and Type I IFN share STAT1; Note S-L). Mean |r| = 0.21 — the modules are near-independent.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_fig_E2_redundancy.png){width=100%}




![**Figure S3. Each principal axis is biologically structured (neither panel uses the modules).** For PC1–PC6, cell position on the axis (x) vs an external quality signature scored from gene-level HPA (y); the axis is monotonically ordered by that quality (immune signalling PC1 ρ = −0.89, translation PC2 +0.75, cholesterol PC5 +0.54). Biological family explains 31–70% of position on every axis (η²; all P < 10⁻⁶; Fig. S1-II). Axis conventions, previously printed inside the figure: x is the position along each principal axis (standardized, the same scale on every PC); y is an external quality signature scored from gene-level HPA without using the modules; the trend line shows the axis is monotonically ordered by that signature; colour is biological family, and families cluster along x, which is what η² measures. The leading Reactome terms at each end of every axis are listed in Note S1.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_fig_axes.png){width=100%}




![**Figure S4. Where the cancers land, shown three ways.** **(I)** A converging against a non-converging cancer, seen directly: malignant cells of B-ALL and neuroblastoma projected onto interpretable 9-D axes — cosine to the core (convergence) against cosine to the neural-program region (neurons + adrenal medulla); 20/50/80% KDE density contours over the 154 normal cells. B-ALL mass sits on the core (red stars); neuroblastoma occupies the neural-program region (toward green triangles). B-ALL is shown because it is the clearest geometric contrast, not as a specimen of cross-lineage convergence: of its 60.8% two-cell-core occupancy, 59.7 points are the megakaryocyte, i.e. largely a return to its own bone-marrow lineage, and only 1.1 the cytotrophoblast. Static positions, not trajectories; the raw PC1–PC2 projection does not resolve the convergence (Methods). **(II)** UMAP of the 154 normal cell types (position = state similarity, not distance to the core; colour = mean top-1 convergence votes from the 25 cancers). The 25 malignant states (◆) pile into one region (dashed hull); large markers = the two-cell core; the manifold's geometric centre (alveolar type-1 pneumocyte, ×) receives zero votes. Embryonic/germline cells co-localize (Fisher P = 5.6 × 10⁻⁶). **(III)** Per-cell-type convergence, column-standardized: the cytotrophoblast row is uniformly warm across solid and liquid cancers (the complexity-robust anchor); the megakaryocyte row runs hot specifically in the hematopoietic cancers (B-ALL, AML, DLBCL) that share its bone-marrow lineage.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS4_geometry.png){width=100%}




![**Figure S5. The transcriptome-complexity trap, and the control for it.** **(I)** (A) The 22-candidate screen; five datasets below the 1,000-gene floor are excluded (AML and MM replaced by GEO cohorts). (B) Per-target complexity: the core sits far above the erythrocyte artifact. **(II)** The megakaryocyte complexity control (Gudkov), across 312,220 malignant cells — the 17 of the 25 panel cohorts whose objects retain per-cell gene counts. (A) Assignment to the megakaryocyte rises as a cell's detected-gene count falls — like the erythrocyte artifact — while the cytotrophoblast rises *with* complexity. (B) Among high-complexity cells (≥ 4,000 genes) the cytotrophoblast dominates (20%) and the megakaryocyte decays to the erythrocyte-residual level (5%). (C) Causal: thinning the *same* cells collapses cytotrophoblast assignment, surges erythrocyte (positive control), and raises megakaryocyte only modestly. The cytotrophoblast convergence is complexity-robust; the megakaryocyte's is genuine in blood lineages (full complexity there) and partly a low-complexity effect in solid tumors.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS5_complexity.png){width=100%}




![**Figure S6. Ranking the candidate targets two ways.** **(I)** Target selection on breadth × depth, every candidate labelled: the SI companion to Fig. 3C, the same two axes with all fourteen recurrent candidates named rather than the seven the main-figure panel can carry. Only the core leads on both recurrence breadth and median depth; the candidate third (gastric progenitor) leads on breadth alone. **(II)** The cross-cancer recurrence ranking that gives the data-driven core: the cytotrophoblast and megakaryocyte lead with a clear gap (18→13); enriched neighbourhood in orange; post-cellular bodies (erythrocyte, platelet) excluded. Before that exclusion the erythrocyte would rank third (16/25, between the core and the gastric progenitor), but it is the documented low-complexity artifact: Table S2 carries each cell's pre-quarantine occupancy and Note S-K states the quarantine rule.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS6_ranking.png){width=100%}




![**Figure S7. How much of the result depends on this particular module catalog.** **(I)** Two nulls of matched size, 150 draws each, scored on 17 cohorts against a five-cell target set (the trophoblast states plus the megakaryocyte) — so these percentages are not in the same currency as the headline occupancy. (A) Fraction of malignant cells reaching that target under the real 44 modules (red) against random genome-wide panels (grey) and against the **same curated genes reshuffled into random modules** (blue). The real catalog clears the random-genome floor (21% against 6%, z = +3.1, P = 0.013) but is **not** separated from its own reshuffled grouping (21% against 14%, z = +1.3, P = 0.11). Above the random-genome floor that places roughly half the effect on the gene set, with the gene-set-versus-grouping contrast unresolved rather than settled. (B) The geometric concentration ratio: the real catalog concentrates **less** than either null (5.1 against 6.0 and 6.8; P = 0.93 and 0.92), so concentration as such is not a property of this catalog. **(II)** Robustness to the exact 44-module choice, rebuilding the whole 25-cancer pipeline per draw (canonical per-cell scoring; eigenspace rebuilt by the correlation-PCA recipe; erythrocyte/platelet excluded as in the main analysis). (A) Dropping five of the 44 modules at random over 200 draws: the fraction of malignant cells assigned to the cytotrophoblast/megakaryocyte core is 17.5 ± 2.5% (unperturbed 17.4%), never below 11.7% (≥ 6× over chance); leave-one-out 14.4–21.2%. (B) The cytotrophoblast anchor stays at the top under the same drops — rank 1 unperturbed, in the top three (by cosine to the malignant attractor) in 84% of five-module drops and 100% of single-module drops. The megakaryocyte, the secondary and partly complexity-dependent target, ranks lower in this rebuilt-eigenspace reconstruction and is not separately assessed here. **(III)** The signature is gene-level and so independent of the grouping: per-gene cytotrophoblast suppression against cancer suppression across the 206 catalog genes (202 detectable in both reference and cancer panel, over the 17 of the 25 panel cohorts whose objects carry gene-level counts; Spearman ρ = +0.56, with the real ρ above all 2,000 gene-label permutations — P < 5 × 10⁻⁴, the floor 2,000 draws can reach; the agreement is rank-wise and not linear, Pearson r = +0.10), holding for genes from internally co-expressed modules (ρ = +0.60) and from the rest (ρ = +0.54). Read together: the convergence needs curated biological genes, does not need this particular partition of them into modules, and does not need this particular set of 44. `scripts/build_E2c_target.py` (data `results/E2c_target_null.csv`), `scripts/build_E2e_robustness.py`, `scripts/build_figS10_genelevel.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS7_catalog.png){width=91%}




![**Figure S8. The core against the normal-cell distribution, module by module.** For each of the twelve most target-defining modules, the distribution of all 154 normal cells (grey) with the core marked; the core sits in the bottom ~1–16th percentile of every one — an outlier outside the normal cloud, not merely at its low edge.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_fig_duetprofile.png){width=100%}




![**Figure S9. The named target is detected in an independent system with no cancer in it.** Mouse fibroblast-to-iPSC reprogramming (Schiebinger et al. 2019; 236,285 cells, 39 timepoints), scored on the same 44 modules through orthologs (43 of 44 recovered) and projected through the same fixed coordinate. The fraction of cells assigned to the cytotrophoblast collapses to 1.8% by day 4, then rises about ninefold, to between 11% and 20% from day 14 onward (median 16.7% over those nine timepoints; Spearman ρ = +0.76 with day, P = 3 × 10⁻⁸) — the window in which reprogramming is known to generate extraembryonic and trophoblast-like off-target fates. The coordinate therefore reports trophoblast-ness where it is independently expected, in another species and another laboratory, with no cancer in the analysis. We do not use this figure to establish that claim. The starting fibroblast timepoints (shaded) score anomalously high and we cannot explain it; if unreprogrammed fibroblasts already score high on cytotrophoblast assignment, the coordinate may be reading proliferative high-output character rather than trophoblast identity, and the collapse-then-rise could be a composition effect. The raw metric is also flatter than the shape-only one, and the object carries no cell-type labels, so we cannot confirm these are the cells the original study annotated as trophoblast-like. Suggestive only. `scripts/E28_reprogramming_validation.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_figS28_reprogramming.png){width=88%}




![**Figure S10. Specification curve of the core convergence (multiverse).** Core convergence recomputed under all 36 crossings of four analytic choices — reference calibration (single-cell against pseudobulk), magnitude handling (raw / within-cell-centered / per-cell rank), core definition (three-cell / two-cell / cytotrophoblast-only), and erythroid handling (excluded / included). (A) The specifications sorted by outcome, each against its own chance baseline (ticks). Under the single-cell calibration the single-cell tumor data require, the convergence exceeds chance in every one of the 18 specifications (median 7.2%, worst case 2.9× chance); the only specifications that collapse to the chance baseline all use the pseudobulk calibration, which independently degenerates (Note S-B). (B) Median convergence by analytic choice — the estimate is stable to core definition and erythroid handling and is carried by module magnitude (the rank transform retains a smaller residual), while the single-cell calibration is the load-bearing preprocessing choice. `scripts/E7_specification_curve.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_figS19_speccurve.png){width=100%}




![**Figure S11. The within-patient increase replicates across organs and survives magnitude removal.** **(I)** Paired gastric cohort (GSE150290): twenty-three patients (Kim et al., 2022), a low-baseline organ analysed only after the account was fixed; its tumour samples are also the panel’s gastric cohort (Note S-U). Per patient, matched adjacent-normal (Pat##-A) and cancer (Pat##-B) epithelium (EPCAM+/CD45−) were projected and scored; lines join the two arms of each of the 22 patients meeting the ≥ 20-epithelial-cell floor (red = deeper in cancer, grey = not; black = mean). Malignant epithelium is deeper than matched normal in 18 of 22 patients (mean 8.8% → 14.2%; cytotrophoblast-specific 7.0% → 11.9%; paired Wilcoxon *P* = 0.001). The one grey line falling to zero is Pat01, whose cancer arm yielded no core-assigned epithelium. **(II)** High-baseline lung deepens too, and the earlier "no deepening" there was a dilution artifact. (A) Lung cohort (GSE131907; Kim et al. 2020), five of whose ten paired patients are also among the panel’s smallest cohort (Note S-U); within-patient normal (nLung) against tumor (tLung) epithelium using the authors' own tumor-subtype annotations: malignant epithelium is deeper than matched normal in 10 of 10 patients (mean 23% → 39%; cytotrophoblast-specific 13% → 29%; paired Wilcoxon *P* = 0.002). (B) Why the pan-tissue pipeline had scored lung as non-deepening: on the CELLxGENE Census lung data it used, the broad epithelial filter admitted a "tumor arm" that is only ~9% annotated-malignant cells and ~91% normal-adjacent alveolar epithelium, so its core fraction (~13%) collapses to the level of a pool of healthy Census donors (~12%); the clean malignant cells sit at ~29%. Three normal-lung values appear here and they are not interchangeable: 23.4% is the matched normal arm of the GSE131907 cohort, 26.7% is the patient-matched adjacent normal of the Census cohort, and 11.7% is the healthy-donor pool. **(III)** The deepening is a genuine module-shape shift rather than transcriptional brightness, and is about constant across baseline. (A) Within-patient deepening (malignant minus matched-normal core fraction) recomputed with the overall-magnitude component of each cell removed (per-cell centering across the 44 modules), in three cohorts: colorectal +10.5 points (44/46; paired Wilcoxon *P* = 10⁻¹³), corrected lung +12.1 (9/10; *P* = 0.004, sign-test *P* = 0.011), and gastric +3.4 (15/22; *P* = 0.04 by Wilcoxon but borderline by sign test, *P* = 0.07) — 63–77% of the raw deepening retained. Patient-matched malignant epithelium therefore sits deeper in the *shape-only* core, not merely brighter; the colorectal and lung results are significant under either test, the smaller gastric cohort suggestive. (B) Additive deepening (points) against normal-tissue baseline across five cohorts: no dependence (*r* = +0.4, *P* = 0.5; mean about +10 points) — the apparent "more room where the baseline is lower" holds only in fold-change (a small-denominator artifact) and is retired. `scripts/E10_heldout_gastric.py`, `scripts/E12_highbaseline_lung.py`, `scripts/E13_lung_dilution_proof.py`, `scripts/E17_deepening_shapeonly.py`, `scripts/E17b_shapeonly_gastric_lung.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS11_withinpat.png){width=93%}




![**Figure S12. The deepening is specific to malignancy, not to injury or invasive remodelling (replacement negative control).** Idiopathic pulmonary fibrosis (Habermann et al. 2020; 33,321 epithelial cells, 22 donors) is the strictest non-malignant comparator available in the same organ as the corrected lung result: chronically injured, ECM-remodelling epithelium containing the aberrant-basaloid and transitional-AT2 populations often described as cancer-like. (A) Deepening relative to each condition's own control epithelium: fibrosis moves +3.8 points on the raw core and +1.3 on the shape-only core (neither significant, *P* = 0.12 and 0.21, 22 donors), while malignancy in the same organ moves +15.7 on the raw core (10 of 10 patients, *P* = 0.002) and +12.1 on the shape-only core (9 of 10, *P* = 0.004; Fig. S11-IIA and S11-IIIA). (B) Within the fibrotic lung the ordering is biologically sensible — transitional-AT2 (18.4%) and proliferating epithelium (16.2%) score highest, ciliated cells lowest (5.2%) — so the state tracks invasive, high-output phenotypes in general, but only malignancy carries a tissue into it. `scripts/E24_regeneration_specificity.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_figS24_fibrosis.png){width=100%}




![**Figure S13. Before invasion and after dissemination.** **(I)** Pre-invasive epithelium already occupies the convergent region (out-of-sample test of the prediction). Within-patient colorectal progression (GSE161277; Zheng et al. 2022): epithelial cells (EPCAM+/CD45−) from matched normal, pre-invasive adenoma and carcinoma of the same patient, projected and scored for core convergence. (A) Full core (cytotrophoblast + megakaryocyte) and (B) the cytotrophoblast-specific fraction, per patient (thin) and mean (red), for the three patients with the full series. Pre-invasive adenoma epithelium is already elevated above matched normal in all three patients — mean core 4.1% (normal) → 11.2% (adenoma) → 16.8% (carcinoma), adenoma reaching two-thirds of the carcinoma level — consistent with the prediction that the convergent state is entered before invasion. One patient (Patient 3) is non-monotone at carcinoma (215 epithelial cells, low power). A simple EPCAM+/CD45− marker gate stands in for author malignancy calls; a larger pre-invasive cohort is the natural confirmation. **(II)** The same prediction across two further organs: directionally supported in breast, refuted in esophagus. (A) Breast (GSE195861): epithelial core fraction by stage — normal breast 7.7% (a single sample), ductal carcinoma *in situ* 11.3% (7 samples), invasive ductal carcinoma 21.8% (6), lymph-node metastasis 13.4% (6). The pre-invasive lesion sits between normal and invasive, at about half the invasive level (DCIS against IDC *P* = 0.07); donors are unpaired and the normal arm is a single sample, so this is directional support, not a test. (B) Esophagus (GSE291080), within patient: Barrett metaplasia against the same donor's own normal squamous esophagus — the metaplasia sits *below* normal in four of four patients (2.1% against 7.7%), and below gastric epithelium in three of three. The prediction fails here. Two readings remain open: that it is confined to neoplastic precursors (adenoma and in-situ carcinoma are clonal neoplasms; Barrett metaplasia is a lineage switch), or that it is simply wrong in esophagus. The available samples are largely non-dysplastic, so esophageal dysplasia — the properly pre-invasive lesion — is untested. **(III)** The state accompanies invasion and is relinquished on colonization (suggestive). Within-patient primary tumour against matched lymph-node metastasis (breast, GSE195861; six pairs): the metastasis sits *lower* in the convergent state than the same patient's primary — means of 21.8% → 13.4% on the raw core (lower in six of six; paired *P* = 0.031, the floor a two-sided Wilcoxon can reach at n = 6) and 14.3% → 9.6% shape-only (five of six; *P* = 0.063). The single evaluable liver pair runs the same way. Of the two portal-vein tumour thrombi available — tumour actively occupying a vessel — one sits well above its primary (49.5% → 69.9% raw, 31.4% → 46.3% shape-only) while the other does not move on the raw core (54.1% → 54.0%) and falls with magnitude removed (49.4% → 42.0%); one case is not a direction. Nodal (not distant) metastases and small numbers: reported as suggestive. `scripts/E9_preinvasive_colon.py`, `scripts/E22_preinvasive_barrett_dcis.py`, `scripts/E25_primary_vs_metastasis.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v12_figS13_progression.png){width=91%}




![**Figure S14. The coordinate moves during targeted therapy, and the movement is not a transcriptome-complexity artifact.** Human PC9 lung-adenocarcinoma cells through eleven days of erlotinib (Aissa et al. 2021); Drop-seq, median 161–556 detected genes per cell. **(A)** Shape-only core fraction per timepoint. The grey dashed curve is the analysis as first run, at a floor of ≥ 500 detected genes; the three blue curves instead restrict *every* timepoint to one common absolute detected-gene window, which is the like-with-like comparison. The level is largely set by that threshold — day 2 reads 18%, 43%, 48%, 53% or 62% at floors of 0, 200, 300, 500 or 1,000 genes — so no absolute core fraction from this dataset is interpretable without its floor, and we report none. The shape is not: the rise to day 2 and the day-4 drug-tolerant persister trough appear in every window above 1,000 genes (n = 67–368 cells per point in the 1,500–6,000 window). **(B)** Why that control is necessary. Within every timepoint, core assignment tracks detected-gene count almost monotonically (Spearman +0.77 to +1.00; at day 2 the lowest decile is 2% and the highest 68%) — the same complexity dependence for which the erythrocyte is quarantined in quality control (Fig. 3D, Note S-Q), here operating inside the time course itself. Limits, stated plainly: survivorship thins the later timepoints, different cells are sampled at each one, and no lineage tracing is available, so movement of cells cannot be separated from selection among them; and the endpoint is **not** stably above baseline — days 9–11 sit 16–21 points above it in the 1,500–6,000 window, return to it in the 1,000–3,000 window, and fall below it when every timepoint is resampled to the untreated cells' own complexity distribution. Read as evidence that the coordinate moves under a real perturbation, not that the state is entered and held. `scripts/E26_dynamics_persisters.py`, `scripts/E26b_qc_floor_sweep.py`.](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_figS26_dynamics.png){width=100%}




![**Figure S15. What malignancy adds within the shared state (within-patient, colorectal).** Malignant against matched-normal core cells compared per module (Pelka; detection-depth controlled). The malignant occupant runs retained stress/integrity sensors (HIF, NRF2, HSF1, UPR-PERK, autophagy) and anabolic programs hotter; the p53 module rises through downstream stress/damage targets (p21/BAX/PUMA/NOXA), not *TP53* itself (transcript unchanged) — a replication-stress footprint, not retained p53 function. The separation is mostly off the target-defining axis (ρ = +0.36): a qualitative detuning toward stress-and-growth. Open circles in (A) are the raw difference before the detection-depth control (balanced 2,640 against 2,496 detected genes, from 4,351 against 2,395 raw).](/Users/teo/Desktop/Desktop - teo’s MacBook Air/research/perceptual_modules/release/manuscript/figures/v7_fig_detuning.png){width=100%}



## Supplementary tables

**Table S1. Datasets (25).**

| cancer                | source         | tissue      | platform     |   n_cells_QC | origin_class   | malignant_label   |
|:----------------------|:---------------|:------------|:-------------|-------------:|:---------------|:------------------|
| HCC                   | Lu 2022        | liver       | 10x          |        54978 | epithelial     | author            |
| BRCA                  | Wu 2021        | breast      | 10x          |        65104 | epithelial     | author            |
| Glioblastoma          | Census         | brain       | 10x          |        25661 | glial          | Census            |
| Ovarian serous        | Census         | ovary       | 10x          |        28950 | epithelial     | Census            |
| RCC                   | Census         | kidney      | 10x          |        26012 | epithelial     | Census            |
| Melanoma              | Census         | skin        | 10x          |        12030 | neural-crest   | Census            |
| Prostate              | Census         | prostate    | 10x          |        26109 | epithelial     | Census            |
| HNSCC                 | Census         | head/neck   | 10x          |        14808 | epithelial     | Census            |
| Neuroblastoma         | Census         | adrenal     | 10x          |        30000 | neural-crest   | Census            |
| B-ALL                 | Census         | bone marrow | 10x          |        29132 | hematopoietic  | Census            |
| DLBCL                 | Census         | lymph node  | 10x          |        12786 | hematopoietic  | Census            |
| Squamous lung         | Census         | lung        | 10x          |        19707 | epithelial     | Census            |
| TNBC                  | Census         | breast      | 10x          |        15419 | epithelial     | Census            |
| PDAC                  | Moncada 2020   | pancreas    | inDrop       |          505 | epithelial     | author            |
| NSCLC                 | Census         | lung        | 10x          |          204 | epithelial     | Census            |
| AML                   | van Galen 2019 | bone marrow | Seq-Well     |         6274 | hematopoietic  | author            |
| MM                    | Ledergor 2018  | bone marrow | InDrops      |         5138 | hematopoietic  | author            |
| BCC                   | Yost 2019      | skin        | 10x          |         5526 | epithelial     | heuristic         |
| Colorectal            | Pelka 2021     | colon       | 10x/Seq-Well |        23712 | epithelial     | author            |
| Small-cell lung       | Census         | lung        | 10x          |        22388 | neuroendocrine | Census            |
| Gastric               | Census         | stomach     | 10x          |         4632 | epithelial     | Census*           |
| Germ-cell (embryonal) | GSE228501      | testis      | 10x          |          443 | germ cell      | marker            |
| Osteosarcoma          | GSE162454      | bone        | 10x          |        16214 | mesenchymal    | marker            |
| Ewing sarcoma         | GSE130024      | bone (PDX)  | 10x          |         2007 | mesenchymal    | marker            |
| Seminoma              | GSE197778      | testis      | 10x          |          660 | germ cell      | marker            |

*malignant_label*: `author` = author/atlas malignant annotation; `Census` = malignant (tumour/epithelial) compartment from CELLxGENE Census; `marker` = marker-based Leiden clustering (Note S1); `heuristic` = heuristic label. `*` (Gastric) marks the Census gastric cohort in the panel — median 1,169 genes/cell, just above the QC floor (a separate, earlier gastric dataset that scored 912 in the 22-candidate screen was excluded; Note S3).

**Table S2. Full 154-cell ranking** (cross-cancer recurrence, median and mean depth, occupancy and QC status; the two quarantined post-cellular bodies carry their pre-quarantine occupancy for reference).



```{=latex}
\begingroup\footnotesize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.06}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{\dimexpr 0.31\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.08\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.11\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.1\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.11\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.18\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.11\linewidth-2\tabcolsep\relax}@{}}
\toprule
\textbf{cell type} & \textbf{recur. (of 25)} & \textbf{median depth \%} & \textbf{mean depth \%} & \textbf{occupancy \%} & \textbf{QC status} & \textbf{pre-quar. occ. \%} \\
\midrule
\endfirsthead
\toprule
\textbf{cell type} & \textbf{recur. (of 25)} & \textbf{median depth \%} & \textbf{mean depth \%} & \textbf{occupancy \%} & \textbf{QC status} & \textbf{pre-quar. occ. \%} \\
\midrule
\endhead
adipocytes & 0 & 0 & 0 & 0 & candidate & 0 \\
adrenal cortex cells & 0 & 0 & 0 & 0 & candidate & 0 \\
adrenal medulla cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
alveolar cells type 1 & 0 & 0 & 0.04 & 0.04 & candidate & 0.04 \\
alveolar cells type 2 & 0 & 0 & 0.04 & 0.01 & candidate & 0.01 \\
astrocytes & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
b-cells & 0 & 0 & 0.04 & 0.04 & candidate & 0.03 \\
basal keratinocytes & 0 & 0.27 & 0.59 & 0.63 & candidate & 0.63 \\
basal prostatic cells & 0 & 0 & 0 & 0 & candidate & 0 \\
bergmann glia & 0 & 0 & 0.03 & 0.04 & candidate & 0.04 \\
brain excitatory neurons & 0 & 0 & 0.04 & 0.05 & candidate & 0.05 \\
brain inhibitory neurons & 0 & 0 & 0.04 & 0.04 & candidate & 0.04 \\
breast hormone-responsive cells & 0 & 0 & 0 & 0 & candidate & 0 \\
breast lactating cells & 0 & 0.07 & 0.12 & 0.15 & candidate & 0.14 \\
breast myoepithelial cells & 1 & 0.5 & 0.94 & 0.98 & candidate & 0.92 \\
breast secretory cells & 0 & 0.38 & 0.92 & 0.91 & candidate & 0.9 \\
cardiomyocytes & 0 & 0 & 0.01 & 0.02 & candidate & 0.02 \\
cdc & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
cholangiocytes & 0 & 0 & 0 & 0 & candidate & 0 \\
choroid plexus epithelial cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
colonocytes & 0 & 0 & 0.24 & 0.28 & candidate & 0.28 \\
cone photoreceptor cells & 0 & 0 & 0.04 & 0.05 & candidate & 0.05 \\
conjunctival goblet cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
corticotrophs & 0 & 0 & 0.07 & 0.08 & candidate & 0.08 \\
cytotrophoblasts & 18 & 10.75 & 11.03 & 11.1 & candidate & 10.89 \\
decidual stromal cells & 0 & 0.03 & 0.26 & 0.16 & candidate & 0.15 \\
differentiating spermatogonia & 5 & 1.12 & 2.53 & 2.42 & candidate & 2.42 \\
distal convoluted tubule cells & 0 & 0 & 0 & 0 & candidate & 0 \\
early primary spermatocytes & 4 & 2.15 & 2.87 & 2.21 & candidate & 2.15 \\
early spermatids & 0 & 0.15 & 0.55 & 0.57 & candidate & 0.54 \\
endometrial ciliated cells & 0 & 0 & 0 & 0 & candidate & 0 \\
endometrial glandular cells & 0 & 0 & 0 & 0 & candidate & 0 \\
endometrial luminal cells & 0 & 0 & 0 & 0 & candidate & 0 \\
endometrial secretory cells & 0 & 0 & 0.21 & 0.08 & candidate & 0.08 \\
endometrial stromal cells & 0 & 0 & 0 & 0 & candidate & 0 \\
enteric stem cells & 2 & 0.53 & 1.43 & 1.57 & candidate & 1.56 \\
enteric transient amplifying cells & 0 & 0 & 0.18 & 0.21 & candidate & 0.21 \\
enterocytes & 0 & 0 & 0.03 & 0.03 & candidate & 0.03 \\
ependymal cells & 0 & 0 & 0.02 & 0.03 & candidate & 0.03 \\
epicardial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
epididymal basal cells & 0 & 0.05 & 0.12 & 0.14 & candidate & 0.14 \\
epididymal clear cells & 4 & 1.12 & 3 & 2.62 & candidate & 2.61 \\
epididymal efferent duct absorptive cells & 0 & 0 & 0.02 & 0.01 & candidate & 0.01 \\
epididymal efferent duct ciliated cells & 0 & 0.22 & 0.32 & 0.33 & candidate & 0.33 \\
epididymal principal cells & 0 & 0.15 & 0.31 & 0.35 & candidate & 0.35 \\
erythrocyte progenitors & 2 & 0.62 & 1.32 & 1.47 & candidate & 1.45 \\
erythrocytes & — & — & — & — & quarantined (post-cellular, low-complexity) & 9.45 \\
esophageal apical cells & 9 & 2.38 & 3.47 & 3.58 & candidate & 2.99 \\
esophageal basal cells & 1 & 1.18 & 1.53 & 1.5 & candidate & 1.49 \\
esophageal suprabasal cells & 1 & 0.4 & 1.36 & 1.2 & candidate & 1.17 \\
extravillous trophoblasts & 0 & 0.03 & 0.57 & 0.59 & candidate & 0.59 \\
fallopian secretory cells & 0 & 0 & 0.15 & 0.04 & candidate & 0.04 \\
fallopian tube ciliated cells & 0 & 0 & 0 & 0 & candidate & 0 \\
fibro-adipogenic progenitors & 0 & 0 & 0 & 0 & candidate & 0 \\
fibroblasts & 0 & 0 & 0 & 0 & candidate & 0 \\
foveolar cells & 0 & 0 & 0.01 & 0 & candidate & 0 \\
gastric chief cells & 10 & 2.73 & 6.52 & 7.33 & candidate & 5.52 \\
gastric progenitor cells & 13 & 4.55 & 10.36 & 9.37 & candidate & 9.15 \\
goblet cells & 0 & 0 & 0.07 & 0.08 & candidate & 0.08 \\
gonadotrophs & 0 & 0 & 0 & 0 & candidate & 0 \\
granulosa cells & 1 & 0.3 & 0.78 & 0.72 & candidate & 0.72 \\
hematopoietic stem cells & 0 & 0 & 0 & 0 & candidate & 0 \\
hepatic stellate cells & 0 & 0 & 0.07 & 0.05 & candidate & 0.05 \\
hepatocytes & 0 & 0 & 0.09 & 0.11 & candidate & 0.11 \\
hofbauer cells & 6 & 1.47 & 2.95 & 3.26 & candidate & 3.09 \\
innate lymphoid cells & 0 & 0 & 0 & 0 & candidate & 0 \\
kupffer cells & 0 & 0.03 & 0.34 & 0.39 & candidate & 0.35 \\
lacrimal acinar cells & 0 & 0 & 0.17 & 0.2 & candidate & 0.2 \\
lactotrophs & 0 & 0 & 0.04 & 0.05 & candidate & 0.05 \\
late primary spermatocytes & 1 & 0.78 & 1.12 & 1 & candidate & 0.89 \\
late spermatids & 11 & 3.82 & 8.34 & 8.62 & candidate & 5.98 \\
leydig cells & 0 & 0 & 0.1 & 0.12 & candidate & 0.12 \\
loop of henle epithelial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
lymphatic endothelial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
macrophages & 0 & 0 & 0.03 & 0.03 & candidate & 0.03 \\
mast cells & 0 & 0 & 0.11 & 0.13 & candidate & 0.12 \\
medullary thymic epithelial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
megakaryocyte progenitors & 0 & 0.03 & 0.06 & 0.06 & candidate & 0.06 \\
megakaryocyte-erythroid progenitors & 0 & 0 & 0.1 & 0.12 & candidate & 0.08 \\
megakaryocytes & 18 & 8.22 & 13.51 & 12.81 & candidate & 8.27 \\
melanocytes & 0 & 0 & 0.17 & 0.2 & candidate & 0.2 \\
mesothelial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
microglia & 0 & 0 & 0.05 & 0.06 & candidate & 0.06 \\
migrating cytotrophoblasts & 5 & 2.45 & 2.99 & 3.18 & candidate & 3.18 \\
monocyte progenitors & 0 & 0 & 0.02 & 0.02 & candidate & 0.01 \\
monocytes & 0 & 0 & 0.03 & 0.04 & candidate & 0.04 \\
mucous neck cells & 1 & 0.12 & 0.43 & 0.42 & candidate & 0.4 \\
myonuclei & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
myosatellite cells & 0 & 0 & 0 & 0 & candidate & 0 \\
müller glia & 0 & 0 & 0.12 & 0.09 & candidate & 0.09 \\
neuroendocrine cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
neutrophil progenitors & 0 & 0 & 0.11 & 0.13 & candidate & 0.1 \\
neutrophils & 0 & 0 & 0.02 & 0.03 & candidate & 0.02 \\
nk-cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
ocular epithelial cells & 0 & 0.03 & 0.31 & 0.27 & candidate & 0.27 \\
oligodendrocyte progenitor cells & 0 & 0 & 0.27 & 0.33 & candidate & 0.33 \\
oligodendrocytes & 0 & 0.03 & 0.33 & 0.37 & candidate & 0.37 \\
oocytes & 0 & 0.23 & 0.49 & 0.45 & candidate & 0.45 \\
other brain neurons & 0 & 0 & 0.08 & 0.09 & candidate & 0.09 \\
ovarian stromal cells & 0 & 0 & 0.14 & 0.14 & candidate & 0.13 \\
pancreatic acinar cells & 6 & 3.33 & 6.97 & 7.45 & candidate & 6.45 \\
pancreatic duct cells & 0 & 0 & 0.11 & 0.09 & candidate & 0.09 \\
pancreatic islet cells & 0 & 0 & 0.07 & 0.08 & candidate & 0.07 \\
paneth cells & 2 & 0.05 & 1.52 & 1.75 & candidate & 1.66 \\
papillary tip epithelial cells & 0 & 0 & 0.02 & 0.01 & candidate & 0.01 \\
parietal cells & 0 & 0.22 & 0.61 & 0.69 & candidate & 0.59 \\
pdcs & 0 & 0 & 0.15 & 0.18 & candidate & 0.15 \\
pericytes & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
peritubular myoid cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
pituicytes/fscs & 0 & 0 & 0 & 0 & candidate & 0 \\
pituitary stem cells & 0 & 0 & 0 & 0 & candidate & 0 \\
plasma cells & 0 & 0 & 0.13 & 0.15 & candidate & 0.12 \\
platelets & — & — & — & — & quarantined (post-cellular, low-complexity) & 2.51 \\
podocytes & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
prostatic club cells & 0 & 0 & 0 & 0 & candidate & 0 \\
prostatic glandular cells & 0 & 0 & 0.02 & 0.03 & candidate & 0.03 \\
prostatic hillock cells & 0 & 0 & 0.1 & 0.04 & candidate & 0.04 \\
proximal tubule cells & 0 & 0 & 0 & 0 & candidate & 0 \\
renal collecting duct intercalated cells & 0 & 0 & 0 & 0 & candidate & 0 \\
renal collecting duct principal cells & 0 & 0 & 0.02 & 0 & candidate & 0 \\
renal connecting tubule cells & 0 & 0 & 0 & 0 & candidate & 0 \\
respiratory basal cells & 0 & 0 & 0.01 & 0 & candidate & 0 \\
respiratory ciliated cells & 0 & 0.03 & 0.26 & 0.22 & candidate & 0.22 \\
respiratory deuterosomal cells & 0 & 0 & 0.06 & 0.02 & candidate & 0.02 \\
respiratory ionocytes & 0 & 0 & 0 & 0 & candidate & 0 \\
respiratory secretory cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
retinal amacrine cells & 0 & 0 & 0.15 & 0.18 & candidate & 0.18 \\
retinal bipolar cells & 0 & 0 & 0.09 & 0.1 & candidate & 0.1 \\
retinal ganglion cells & 0 & 0 & 0.02 & 0.03 & candidate & 0.03 \\
retinal horizontal cells & 0 & 0 & 0.1 & 0.1 & candidate & 0.1 \\
retinal pigment epithelial cells & 0 & 0 & 0.08 & 0.05 & candidate & 0.05 \\
rod photoreceptor cells & 0 & 0 & 0.02 & 0.01 & candidate & 0.01 \\
salivary acinar cells & 0 & 0 & 0.01 & 0.01 & candidate & 0.01 \\
salivary basal cells & 0 & 0 & 0 & 0 & candidate & 0 \\
salivary duct cells & 0 & 0 & 0.06 & 0.04 & candidate & 0.04 \\
salivary ionocytes & 0 & 0 & 0.01 & 0 & candidate & 0 \\
salivary myoepithelial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
schwann cells & 0 & 0 & 0 & 0 & candidate & 0 \\
sertoli cells & 0 & 0 & 0.16 & 0.19 & candidate & 0.19 \\
smooth muscle cells & 0 & 0 & 0.02 & 0 & candidate & 0 \\
somatotrophs & 0 & 0 & 0 & 0 & candidate & 0 \\
submucosal glandular cells & 0 & 0 & 0.03 & 0.04 & candidate & 0.04 \\
suprabasal keratinocytes & 0 & 0.15 & 0.48 & 0.52 & candidate & 0.5 \\
syncytiotrophoblasts & 2 & 1.18 & 2.11 & 2.2 & candidate & 2.19 \\
t-cells & 0 & 0 & 0 & 0 & candidate & 0 \\
thymic myoid cells & 0 & 0 & 0 & 0 & candidate & 0 \\
thymocytes & 0 & 0 & 0 & 0 & candidate & 0 \\
thyrotrophs & 0 & 0 & 0.16 & 0.19 & candidate & 0.19 \\
transitional alveolar cells & 0 & 0 & 0 & 0 & candidate & 0 \\
tuft cells & 0 & 0.15 & 0.41 & 0.44 & candidate & 0.43 \\
undifferentiated spermatogonia & 2 & 0.05 & 1.1 & 1.23 & candidate & 1.23 \\
urothelial cells & 0 & 0 & 0.08 & 0.05 & candidate & 0.05 \\
vascular endothelial cells & 0 & 0 & 0 & 0 & candidate & 0 \\
vascular smooth muscle cells & 0 & 0 & 0 & 0 & candidate & 0 \\
\bottomrule
\end{longtable}
\endgroup
```


- **Table S3** (inline, above; `tableS3_invariance.csv`). Coordinate-invariance compactness, region vs core, three bases.
**Table S4. Per-module cancer→core contributions** (44 modules, all 25 cohorts). Each cancer centroid is standardized against the reference (mean and SD over the 154 normal cell types) and unit-normalized; the two-component signature is the unit-normalized mean of the standardized cytotrophoblast and megakaryocyte profiles; a module contribution is the elementwise term of the cosine between them, so the 44 contributions sum to that cosine — here 0.711, which is the mean over the 25 cancers of the cosine between a cancer centroid and the two-component signature in standardized module space. It is not the quantity plotted in Fig. 3A, which is the cosine of each normal cell type to the mean malignant direction in the projected nine-dimensional space (0.89–0.92 for the core). mean_contribution and sd are over the 25 cancers; pct_positive is the fraction of cancers with a positive term. (An earlier version was computed on a 16-cohort panel, which is why its values fell on a k/16 grid; recomputed by scripts/E33_fix_tables.py.)

| module       |   mean_contribution |     sd |   pct_positive |   n_cohorts |
|:-------------|--------------------:|-------:|---------------:|------------:|
| Circadian    |              0.0831 | 0.008  |           1    |          25 |
| MR           |              0.0576 | 0.0113 |           1    |          25 |
| Insulin/FOXO |              0.0535 | 0.014  |           1    |          25 |
| Wnt          |              0.0416 | 0.0125 |           1    |          25 |
| PI3K/PTEN    |              0.0412 | 0.0204 |           0.96 |          25 |
| JAK-STAT     |              0.0382 | 0.0175 |           0.96 |          25 |
| NFAT         |              0.0342 | 0.0107 |           1    |          25 |
| PR           |              0.0273 | 0.0101 |           0.96 |          25 |
| BMP          |              0.027  | 0.0159 |           0.96 |          25 |
| TR           |              0.0263 | 0.0043 |           1    |          25 |
| GR           |              0.0246 | 0.0169 |           0.92 |          25 |
| Hippo        |              0.0209 | 0.003  |           1    |          25 |
| Notch        |              0.0207 | 0.0184 |           0.88 |          25 |
| Type I IFN   |              0.0184 | 0.0236 |           0.76 |          25 |
| SREBP        |              0.0162 | 0.0099 |           0.92 |          25 |
| NRF2         |              0.0152 | 0.0116 |           0.96 |          25 |
| AMPK         |              0.0149 | 0.0061 |           0.96 |          25 |
| ERK/MAPK     |              0.0145 | 0.0106 |           0.88 |          25 |
| TGF-β        |              0.0136 | 0.0029 |           1    |          25 |
| NF-κB        |              0.0126 | 0.0069 |           0.96 |          25 |
| Hedgehog     |              0.0117 | 0.005  |           0.96 |          25 |
| AR           |              0.0114 | 0.0102 |           0.96 |          25 |
| Autophagy    |              0.0097 | 0.0035 |           0.96 |          25 |
| AhR          |              0.0088 | 0.0082 |           0.8  |          25 |
| PPARα        |              0.0082 | 0.0054 |           0.92 |          25 |
| Calcium      |              0.0077 | 0.0042 |           0.92 |          25 |
| UPR-IRE1     |              0.0074 | 0.0086 |           0.84 |          25 |
| cAMP/CREB    |              0.0071 | 0.0026 |           0.96 |          25 |
| UPR-PERK     |              0.0069 | 0.0087 |           0.76 |          25 |
| mTOR         |              0.0064 | 0.0055 |           0.84 |          25 |
| FXR          |              0.0063 | 0.0062 |           0.92 |          25 |
| LXR          |              0.0062 | 0.0015 |           1    |          25 |
| ER           |              0.0056 | 0.0109 |           0.8  |          25 |
| cGAS-STING   |              0.0041 | 0.022  |           0.48 |          25 |
| VDR          |              0.0031 | 0.0035 |           0.88 |          25 |
| PXR/CAR      |              0.0028 | 0.002  |           0.88 |          25 |
| p53          |              0.0004 | 0.0006 |           0.76 |          25 |
| UPR-ATF6     |              0.0002 | 0.0018 |           0.64 |          25 |
| RAR          |             -0      | 0      |           0.16 |          25 |
| NPAS4        |             -0.0002 | 0.0003 |           0.16 |          25 |
| HIF          |             -0.0007 | 0.0008 |           0.2  |          25 |
| PPARγ        |             -0.0008 | 0.0002 |           0    |          25 |
| Cell Cycle   |             -0.0013 | 0.0058 |           0.48 |          25 |
| HSF1         |             -0.0015 | 0.0023 |           0.24 |          25 |
- **Table S5** (inline, below; progression). HCC stage and colorectal grade/T/N.
**Table S6. Transcriptome-complexity diagnostic** (median detected genes per assigned cell type). The twenty cell types receiving the most assignments, then the five lowest-complexity types below them; all 144 types that receive at least one assignment, of the 154 in the reference, are in the attached `tableS6_complexity.csv`. The QC floor is 1,000 detected genes per cell, and the two quarantined post-cellular bodies sit below it.

| cell type                             |   cells assigned |   median genes |   mean genes |
|:--------------------------------------|-----------------:|---------------:|-------------:|
| erythrocytes                          |            82227 |          890   |         1257 |
| megakaryocytes                        |            68023 |         1396   |         1737 |
| pancreatic acinar cells               |            57283 |          992   |         1202 |
| cytotrophoblasts                      |            45556 |         2430   |         2796 |
| late spermatids                       |            37601 |          763   |         1024 |
| gastric chief cells                   |            35489 |          904   |         1309 |
| gastric progenitor cells              |            27280 |         2506.5 |         2774 |
| platelets                             |            25709 |         1541   |         1732 |
| esophageal apical cells               |            23060 |         1521   |         1876 |
| hofbauer cells                        |            15285 |         2527   |         2860 |
| epididymal clear cells                |            12546 |         3150   |         3493 |
| breast myoepithelial cells            |             9953 |         1630   |         1891 |
| differentiating spermatogonia         |             9423 |         2483   |         2771 |
| early primary spermatocytes           |             8545 |         2384   |         2787 |
| syncytiotrophoblasts                  |             8048 |         2990   |         3367 |
| migrating cytotrophoblasts            |             7555 |         3802   |         3864 |
| esophageal suprabasal cells           |             7275 |         2751   |         2979 |
| paneth cells                          |             7150 |         2456   |         2755 |
| undifferentiated spermatogonia        |             7113 |         2576   |         2762 |
| late primary spermatocytes            |             7003 |         1202   |         1439 |
| hepatic stellate cells                |              215 |         1320   |         1445 |
| early spermatids                      |             3071 |         1449   |         1607 |
| ovarian stromal cells                 |              660 |         1480.5 |         1684 |
| renal collecting duct principal cells |                2 |         1491   |         1491 |
| fallopian tube ciliated cells         |                1 |         1496   |         1496 |
**Table S7. The core's target-defining module signature** (44 modules; target-distinctiveness vs cancer contribution).



```{=latex}
\begingroup\footnotesize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.06}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{\dimexpr 0.23\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.16\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.15\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.15\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.16\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.15\linewidth-2\tabcolsep\relax}@{}}
\toprule
\textbf{module} & \textbf{target distinct.} & \textbf{core suppression} & \textbf{megakaryocyte z} & \textbf{cytotrophoblast z} & \textbf{cancer contrib.} \\
\midrule
\endfirsthead
\toprule
\textbf{module} & \textbf{target distinct.} & \textbf{core suppression} & \textbf{megakaryocyte z} & \textbf{cytotrophoblast z} & \textbf{cancer contrib.} \\
\midrule
\endhead
Circadian & -1.828 & 1.828 & -2.33 & -1.21 & 0.083 \\
Type I IFN & -1.612 & 1.612 & -1.65 & -1.5 & 0.018 \\
MR & -1.581 & 1.581 & -1.25 & -1.84 & 0.058 \\
JAK-STAT & -1.578 & 1.578 & -1.2 & -1.89 & 0.038 \\
PI3K/PTEN & -1.497 & 1.497 & -1.47 & -1.46 & 0.041 \\
Insulin/FOXO & -1.489 & 1.489 & -1.49 & -1.43 & 0.053 \\
Notch & -1.339 & 1.339 & -1.87 & -0.76 & 0.021 \\
PR & -1.226 & 1.226 & -1.22 & -1.19 & 0.027 \\
BMP & -1.141 & 1.141 & -1.58 & -0.66 & 0.027 \\
NFAT & -1.108 & 1.108 & -0.95 & -1.24 & 0.034 \\
Wnt & -1.093 & 1.093 & -1.93 & -0.21 & 0.042 \\
SREBP & -0.981 & 0.981 & -1.28 & -0.66 & 0.016 \\
GR & -0.973 & 0.973 & -0.74 & -1.19 & 0.025 \\
TR & -0.936 & 0.936 & -1.25 & -0.6 & 0.026 \\
UPR-PERK & -0.83 & 0.83 & -1.14 & -0.5 & 0.007 \\
ER & -0.808 & 0.808 & -0.7 & -0.89 & 0.006 \\
ERK/MAPK & -0.725 & 0.725 & -1.44 & 0.01 & 0.015 \\
NF-κB & -0.719 & 0.719 & -0.62 & -0.8 & 0.013 \\
Hippo & -0.691 & 0.691 & -1.46 & 0.1 & 0.021 \\
AR & -0.684 & 0.684 & -0.22 & -1.13 & 0.011 \\
UPR-IRE1 & -0.674 & 0.674 & -1.3 & -0.04 & 0.007 \\
Hedgehog & -0.666 & 0.666 & -0.78 & -0.54 & 0.012 \\
AhR & -0.664 & 0.664 & -0.67 & -0.65 & 0.009 \\
NRF2 & -0.656 & 0.656 & -0.91 & -0.39 & 0.015 \\
AMPK & -0.642 & 0.642 & -0.1 & -1.17 & 0.015 \\
FXR & -0.604 & 0.604 & -0.62 & -0.58 & 0.006 \\
PPARα & -0.478 & 0.478 & -0.5 & -0.45 & 0.008 \\
VDR & -0.398 & 0.398 & -0.24 & -0.55 & 0.003 \\
TGF-β & -0.396 & 0.396 & 0.45 & -1.23 & 0.014 \\
PXR/CAR & -0.369 & 0.369 & -0.45 & -0.28 & 0.003 \\
Calcium & -0.346 & 0.346 & 0.36 & -1.04 & 0.008 \\
Autophagy & -0.31 & 0.31 & -0.48 & -0.13 & 0.01 \\
cAMP/CREB & -0.295 & 0.295 & -0.5 & -0.09 & 0.007 \\
Cell Cycle & -0.291 & 0.291 & 0.23 & -0.81 & -0.001 \\
mTOR & -0.269 & 0.269 & 0.3 & -0.83 & 0.006 \\
LXR & -0.21 & 0.21 & -1.04 & 0.63 & 0.006 \\
HSF1 & -0.071 & 0.071 & -0.24 & 0.1 & -0.001 \\
p53 & -0.046 & 0.046 & -0.44 & 0.35 & 0 \\
HIF & -0.029 & 0.029 & -0.7 & 0.64 & -0.001 \\
RAR & 0.002 & -0.002 & -1.35 & 1.35 & -0 \\
NPAS4 & 0.022 & -0.022 & 0.62 & -0.58 & -0 \\
PPARγ & 0.066 & -0.066 & -0.04 & 0.17 & -0.001 \\
UPR-ATF6 & 0.097 & -0.097 & -0.72 & 0.92 & 0 \\
cGAS-STING & 0.977 & -0.977 & 0.77 & 1.16 & 0.004 \\
\bottomrule
\end{longtable}
\endgroup
```


**Table S8. Cancer → cell-of-origin mapping** (25; nominal annotation — the Fig. 3 contraction is tested against random 25-cell draws and does not depend on the pairings).

| cancer          | cell_of_origin (HPA)            | origin_class   |
|:----------------|:--------------------------------|:---------------|
| HCC             | hepatocytes                     | epithelial     |
| BRCA            | breast hormone-responsive cells | epithelial     |
| Glioblastoma    | astrocytes                      | glial          |
| Ovarian serous  | fallopian secretory cells       | epithelial     |
| RCC             | proximal tubule cells           | epithelial     |
| Melanoma        | melanocytes                     | neural-crest   |
| Prostate        | prostatic glandular cells       | epithelial     |
| HNSCC           | suprabasal keratinocytes        | epithelial     |
| Neuroblastoma   | adrenal medulla cells           | neural-crest   |
| B-ALL           | b-cells                         | hematopoietic  |
| DLBCL           | b-cells                         | hematopoietic  |
| Squamous lung   | respiratory ciliated cells      | epithelial     |
| TNBC            | breast myoepithelial cells      | epithelial     |
| PDAC            | pancreatic duct cells           | epithelial     |
| NSCLC           | alveolar cells type 2           | epithelial     |
| AML             | hematopoietic stem cells        | hematopoietic  |
| MM              | plasma cells                    | hematopoietic  |
| BCC             | basal keratinocytes             | epithelial     |
| Colorectal      | colonocytes                     | epithelial     |
| Small-cell lung | neuroendocrine cells            | neuroendocrine |
| Gastric         | mucous neck cells               | epithelial     |
| Germ-cell       | undifferentiated spermatogonia  | germ cell      |
| Osteosarcoma    | fibroblasts                     | mesenchymal    |
| Ewing sarcoma   | fibroblasts                     | mesenchymal    |
| Seminoma        | undifferentiated spermatogonia  | germ cell      |
**Table S9. The 44 regulatory modules and their component genes** (206 unique core genes — the explicit catalog underlying the eigenspace).

| module       |   n_genes | core_genes                            |
|:-------------|----------:|:--------------------------------------|
| AMPK         |         5 | PRKAA1, PRKAA2, PRKAB1, STK11, ADIPOQ |
| AR           |         5 | AR, KLK3, TMPRSS2, FKBP5, NKX3-1      |
| AhR          |         5 | AHR, ARNT, CYP1A1, CYP1B1, AHRR       |
| Autophagy    |         5 | ULK1, BECN1, MAP1LC3B, TFEB, ATG5     |
| BMP          |         5 | BMPR1A, SMAD1, SMAD5, SMAD6, ID1      |
| Calcium      |         4 | CALM1, CAMK2A, ITPR1, PLCG1           |
| Cell Cycle   |         5 | CCND1, CDK4, RB1, E2F1, MKI67         |
| Circadian    |         5 | CLOCK, ARNTL, PER1, NR1D1, CRY1       |
| ER           |         5 | ESR1, TFF1, GREB1, PGR, CCND1         |
| ERK/MAPK     |         5 | MAPK1, MAPK3, MAP2K1, DUSP1, FOS      |
| FXR          |         5 | NR1H4, NR0B2, FGF19, ABCB11, SLC51A   |
| GR           |         5 | NR3C1, FKBP5, TSC22D3, SGK1, DUSP1    |
| HIF          |         5 | HIF1A, VEGFA, EGLN1, SLC2A1, LDHA     |
| HSF1         |         5 | HSF1, HSPA1A, HSP90AA1, HSPA4, DNAJB1 |
| Hedgehog     |         5 | SMO, PTCH1, GLI1, GLI2, SUFU          |
| Hippo        |         5 | YAP1, WWTR1, LATS1, LATS2, TEAD1      |
| Insulin/FOXO |         5 | INSR, IRS1, AKT2, FOXO1, GSK3B        |
| JAK-STAT     |         5 | JAK1, STAT1, STAT3, SOCS1, SOCS3      |
| LXR          |         5 | NR1H3, ABCA1, ABCG1, MYLIP, NR1H2     |
| MR           |         5 | NR3C2, SGK1, SCNN1A, NDRG2, FKBP5     |
| NF-κB        |         5 | NFKB1, RELA, IKBKB, NFKBIA, TNFAIP3   |
| NFAT         |         5 | NFATC1, NFATC2, PPP3CA, RCAN1, IL2    |
| NPAS4        |         5 | NPAS4, ARNT2, CAMK2A, CAMK2B, CALM1   |
| NRF2         |         5 | NFE2L2, KEAP1, NQO1, HMOX1, TXNRD1    |
| Notch        |         5 | NOTCH1, HES1, HEY1, DLL1, RBPJ        |
| PI3K/PTEN    |         5 | PIK3CA, AKT1, PTEN, FOXO3, PDK1       |
| PPARα        |         5 | PPARA, CYP4A11, ACOX1, CPT1A, FABP1   |
| PPARγ        |         5 | PPARG, FABP4, LPL, ADIPOQ, CD36       |
| PR           |         5 | PGR, TNFSF11, WNT4, HAND2, FOXO1      |
| PXR/CAR      |         5 | NR1I2, CYP3A4, CYP2B6, ABCB1, UGT1A1  |
| RAR          |         5 | RARA, RARB, CYP26A1, CRABP2, TGM2     |
| SREBP        |         5 | SREBF1, SREBF2, HMGCR, FASN, SCD      |
| TGF-β        |         5 | TGFB1, TGFBR1, SMAD2, SMAD3, SMAD7    |
| TR           |         5 | THRA, THRB, DIO1, DIO2, THRSP         |
| Type I IFN   |         5 | IFNAR1, STAT1, STAT2, IRF9, USP18     |
| UPR-ATF6     |         5 | ATF6, HSPA5, CALR, PDIA4, XBP1        |
| UPR-IRE1     |         5 | ERN1, XBP1, DNAJB9, SEC61A1, EDEM1    |
| UPR-PERK     |         5 | EIF2AK3, ATF4, DDIT3, PPP1R15A, ASNS  |
| VDR          |         5 | VDR, CYP24A1, CYP27B1, CAMP, S100A8   |
| Wnt          |         5 | CTNNB1, APC, GSK3B, AXIN2, TCF7       |
| cAMP/CREB    |         5 | PRKACA, CREB1, PDE4B, ADCY6, CREM     |
| cGAS-STING   |         5 | MB21D1, TMEM173, IRF3, IFNB1, CXCL10  |
| mTOR         |         5 | MTOR, RPTOR, RPS6KB1, EIF4EBP1, RHEB  |
| p53          |         5 | TP53, MDM2, CDKN1A, BAX, BBC3         |

**Table S5. Within-cancer progression (no strong gradient).**

| Cohort | Test | Statistic | n | Verdict |
|---|---|---|---:|---|
| HCC (Lu 2022) | convergence vs clinical stage I–IV | Spearman ρ = +0.40, P = 0.26 | 10 patients | weak positive, not significant; component not excluded |
| Colorectal (Pelka 2021) | convergence vs grade | Mann–Whitney, P > 0.05 | 62 patients | null |
| Colorectal (Pelka 2021) | convergence vs T-stage | Mann–Whitney, P > 0.05 | 62 patients | null |
| Colorectal (Pelka 2021) | convergence vs nodal status | Mann–Whitney, P > 0.05 | 62 patients | null |

n_genes-controlled; convergence is substantial already at the earliest stages sampled.


**Terminology and boundaries (glossary).** *Regulatory module* and *perceptual module* denote the same object; the main text uses "regulatory module" throughout. Region boundaries, used consistently: the **two-cell core** = cytotrophoblast + megakaryocyte; the **three-cell core** = those two plus the migrating cytotrophoblast (used in the control notes; see Note S-K); the **four-cell majority-enriched neighbourhood** and the **14-cell broader region** are successively wider boundaries used only in the boundary-independence check (Table S2), which shows the enrichment survives the choice. *The vocabulary, fixed.* The paper uses several near-synonyms and they are not interchangeable, so each is pinned here. **The region** (equivalently *convergent region*, *convergent state*) is whichever of the four boundaries above is under discussion, and is always qualified in context; it is the object of the claim this paper defends. **The core** without a number means the two-cell core. **The target** is the cell a ranking puts nearest the region — a property of the ranking, not an identification of the state (§3). **The anchor** is the target under *this* module catalog and at the readiness level of measurement, i.e. the cytotrophoblast; independently reconstructed catalogs reproduce the region and anchor it elsewhere (Note S-O), and output-level scoring does the same (Note S-T), so "anchor" never carries a claim beyond this catalog. **Occupancy** is the share of cells whose top-1 nearest normal neighbour is a member of the boundary named; **depth** is the per-cancer occupancy; **recurrence** is the number of cancers in which a cell is among the five most frequent nearest neighbours. Where a sentence needs the distinction between region and cell and does not make it, the region is meant.


```{=latex}
\clearpage
```

**Note S-A. Robustness of the convergence to the assignment rule and to transcriptional magnitude (`scripts/E3_referee_controls.py`).** Re-scoring the 25 cancers (83,819 malignant cells) and re-projecting: the winner-take-all top-1 core fraction (27%) is matched by a soft cosine-softmax assignment (21% core mass) and by top-k membership (core in the top 3 of 152 for 57% of cells, top 5 for 75%; median nearest-core rank 3), so the convergence is not an artifact of hard assignment. It is not the low-complexity/magnitude artifact that inflated the erythrocyte: core-assigned cells carry higher transcriptional complexity than non-core (median 2,629 vs 2,360 detected genes), and cells assigned to a centroid do not track that centroid's projection norm (r = +0.14). Applying the same magnitude-removing transform to BOTH reference and cancer (within-cell centering, or per-cell rank transform) reduces the core convergence from 27% to 7–11% (still 4–5× over the ~2% chance baseline): the effect is largely a signal of overall regulatory output (the terminal-high-output biology) with a smaller, real module-shape residual.



**Note S-B. Calibration dependence (`scripts/E4_calibration.py`).** Scoring the 25 cancers once and projecting under each reference calibration: under the single-cell-scaled calibration the core fraction is 17% and the 25 malignant centroids are tighter than every one of 2,000 random 25-cell draws; under a pseudobulk (population-averaged) calibration the core fraction falls to 7% and the projection degenerates, collapsing tumor cells toward low-complexity non-core targets (mean pairwise centroid distance 0.001). Convergence onto the core is present under both calibrations but is properly resolved only under the single-cell calibration — the appropriate one for single-cell tumor data, and the choice we adopt on unit-matching grounds (not because it maximizes the effect).



**Note S-C. The address is specific, not merely a high-output signal (`scripts/E5_specificity.py`, 83,819 malignant cells, 25 cancers).** Four independent tests separate *which cell is named* from *overall transcriptional level*. (i) *Shape-only ranking.* Applying the magnitude-removing rank transform to both reference and cancer and re-ranking the cross-cancer recurrence, the cytotrophoblast remains the single most recurrent normal target (14/25 cancers) and the syncytiotrophoblast — a second trophoblast state — rises to second (9/25); the megakaryocyte falls from tied-first to sixth. Magnitude removal therefore *sharpens* the address onto the trophoblast lineage (the named anchor) and reveals that the residual magnitude dependence lives in the megakaryocyte arm, not in the cytotrophoblast. (ii) *High-output specificity.* The professional high-output secretory normals that a pure-magnitude account predicts should also attract convergence do not: of all cross-cancer top-1 votes, hepatocytes receive 0.11%, plasma cells 0.15%, parietal cells 0.69%, prostatic glandular 0.03%, proximal-tubule 0.00% — against 11.1% (cytotrophoblast) and 12.8% (megakaryocyte) for the core. The only high-baseline normals that draw appreciable convergence, pancreatic acinar (7.5%) and gastric chief (7.3%), are the resident secretory cells of two cancers in the panel (PDAC, gastric) — tissue-of-origin return rather than cross-lineage convergence. (iii) *Library-size control.* Core assignment is essentially uncorrelated with sequencing depth: r(core, log10 total counts) = +0.03, r(core, detected genes) = +0.04. (iv) *Orthogonal markers.* Scored on trophoblast/EVT genes lying entirely outside the 206 module genes, core-assigned malignant cells express the extravillous-trophoblast invasion program more highly than non-core (MMP2 0.140 vs 0.082; MMP9 0.098 vs 0.042; ITGA5 0.090 vs 0.064), though the canonical identity markers (HLA-G, CGA/CGB, KRT7, TFAP2C) are not specifically enriched — a modest, invasion-flavored independent corroboration of trophoblast-ness, not a full identity match.



**Note S-D. Within-patient deepening under a mixed-effects model.** The pooled sign test (77/86 patients deeper) treats patients as exchangeable across the four low-baseline datasets, which risks a Simpson's-paradox confound. A linear mixed-effects model of the per-patient malignant-minus-normal core difference with a dataset random intercept estimates a fixed deepening of +13.0 percentage points (SE 1.87; z = 6.95; one-sided *P* = 1.9 × 10⁻¹²), with all four datasets deepening in the mean (colorectal-Pelka 43/46, colorectal-HTA 23/25, kidney 6/7, liver 5/8). The effect thus survives crediting the between-dataset clustering; we note the coarseness of a four-group random-effect estimate (two of the four cohorts colorectal), so the low-baseline conclusion is carried by colorectal with same-signed support from kidney and liver.



**Note S-F. The deepening survives within-pair complexity matching (`scripts/E15_deepening_complexity_matched.py`).** Malignant epithelium carries more detected genes than matched-normal (colorectal median 3,852 vs 2,272), and core assignment is weakly complexity-dependent, so a within-pair complexity offset could in principle manufacture deepening. Resampling the malignant cells of each patient to the *matched-normal* detected-gene distribution (nearest-gene-count matching, 20 draws) and recomputing the core fraction: the deepening falls only from +15.5 to +12.9 percentage points (43/46 patients, paired *P* = 10⁻⁸), retaining ~80% of the raw effect. The paired complexity offset therefore accounts for at most a fifth of the deepening.



**Note S-G. The contraction is mostly magnitude and averaging, with a real shape residue (`scripts/E14_contraction_matched_null.py`).** The raw 25-cancer clustering (mean pairwise cosine distance ~0.166, vs ~0.99 among the near-orthogonal origins) is largely a matched-averaging and shared-magnitude effect: 25 equal-sized coherent normal micro-populations, formed by the same averaging over a seven-tissue normal pool, cluster nearly as tightly (null mean 0.129, 95% interval [0.074, 0.213]), and the malignant clustering is flat across subsample size (0.174 at 25 cells/cancer, 0.166 at 4,000), so it is a stable biological quantity rather than the shrinkage of large means. Removing the overall-magnitude component isolates a residue that we report as a specification curve rather than a single number, because the comparison is sensitive to how the magnitude-removed scores are rescaled before projection. (i) In the **raw module space**, with no standardization or projection at all, the malignant centroids cluster at 0.326 against a null mean of 0.493 — below every one of 200 draws (and in the same space the raw metric is 0.168 vs 0.287, likewise below every draw). (ii) Under the **reference standardization used throughout the paper**, the shape-only clustering is 0.026 against a null of 0.034 (the observed value falls below 292 of 300 matched-averaging draws; permutation *P* = (8+1)/(300+1) = 0.030). (iii) Under a third standardization — dividing by the standard deviations of the *magnitude-removed* reference, which up-weights modules whose pattern barely varies across the 154 normal types and so amplifies low-variance axes — the advantage vanishes (0.400 vs 0.382, 66th percentile). Two of three specifications, including the simplest and assumption-free one, support a real shape convergence; one does not. We report all three and no longer headline the magnitude-inflated 6.3-fold figure. (The coherent-micro-population null is if anything biased tight, making the surviving results conservative; a cell-type-level null is the natural refinement.)



**Note S-H. Cancer-blind justification of the single-cell calibration (`scripts/E16_calibration_cancerblind.py`).** To show the single-cell-scaled calibration is chosen on principled, not result-maximizing, grounds, we test it on 16,312 independent *normal* single cells (skin, kidney, breast; author cell types; no cancer). Projected and assigned to their nearest normal reference type: the **pseudobulk** calibration degenerates — it collapses these normal single cells onto only 2 of 152 targets (effective assignment diversity 1.0), 99.6% of them forced onto low-complexity erythroid cells — the same degeneration it produces on tumor cells (Note S-B). The **single-cell-scaled** calibration keeps the assignments diverse and complexity-appropriate (85 targets, diversity 14.6). Pseudobulk is thus wrong for single-cell data of any kind, normal or malignant; the single-cell calibration is required by the data type and selected with no cancer input. (A fine-cell-type silhouette metric did not distinguish the calibrations; we report the degeneration test, not the metric that would have favored our choice.)



**Note S-I. Malignant-cell purity and the direction of its bias.** The lung dilution artifact (Fig. S11-II) shows automated malignant selection can be wrong in a specific cohort, so we state the purity basis of each within-patient cohort and the direction of any residual bias. The cohorts carrying the deepening use author-annotated malignant compartments — Pelka `ClusterMidway = EpiT` (tumor epithelial clusters) and lung Kim-2020 tumor subtypes tS1/tS2/tS3 — which are high-purity; the gastric arm uses an EPCAM+/CD45− gate on the cancer sample with no malignant annotation, so it may retain some normal-adjacent epithelium. Crucially the bias is signed: normal-adjacent contamination pulls a tumor arm toward its normal baseline and therefore *understates* deepening (the lung case: a 9%-malignant arm read ~13% ≈ normal, versus ~29% for the clean malignant cells; `E13`). Marker-based selection (germ-cell/sarcoma) biases convergence the other way and is bounded by the drop-4 robustness. We do not, however, claim that contamination makes the reported numbers a lower bound: we tested that argument directly and it was not borne out (below). We tested the contamination argument directly on the two cohorts whose labels allow it — restricting HCC to hepatocytes from tumour sites and BRCA to annotated cancer-epithelial cells — and it is *not* borne out: cleaning lowered rather than raised those cohorts (HCC 24.3% → 18.6%, BRCA 31.0% → 26.7%) and left the panel essentially unchanged (core 27.1% → 26.6% raw, 17.3% → 17.1% shape-only; clustering 0.166 → 0.166). Contamination therefore does not appear to be masking convergence, and we report the direction-of-bias argument as untested-and-unsupported rather than as a defence. A CNV-inference malignant floor across all 25 cohorts is deferred to production, as the convergence is now presented as the modest, magnitude-carried companion (Note S-G). Separately, on the deepening-vs-baseline relationship: across five cohorts the additive deepening does not depend on baseline (Pearson *r* = +0.4, *P* = 0.5; mean ~+10 points), the fold-change inverse being a small-denominator artifact — so we make no baseline-gating claim (Fig. S11-IIIB).



**Note S-J. What moves the state, and what does not (perturbation and immune controls).** *Single genes do not.* Across the genome-wide Perturb-seq resource (Replogle et al. 2022; 9,675 CRISPRi knockdowns in K562 with 75,328 non-targeting controls), knocking down a target-defining brake does not move cells toward the core: the correlation between a gene's brake score and its knockdown's core shift is +0.07 (*P* = 0.4), +0.13 after controlling for the generic brightening that essential-gene knockdown causes (*P* = 0.1). A three-level decomposition shows where the signal stops (`scripts/E21_why_null_diagnostic.py`): CRISPRi works (target mRNA down, median −0.43 SD), the containing module responds strongly (LDHA→HIF −1.84 SD, HSP90AA1→HSF1 −1.65; and, as a positive control that the modules read pathway biology rather than their own transcripts, KEAP1 knockdown *raises* the NRF2 module by +0.96 SD), but among the 23 strong module hits the nine-dimensional position moves only a fifth as much and in no consistent direction (13 of 23 toward the core; *r* = −0.10 between module movement and coordinate movement). We read this as a many-module consensus buffered against single-gene perturbation — consistent with a state defined by a dozen simultaneously suppressed brakes, and with the large movement seen under whole-network pressure (Fig. S14). We label that reading for what it is: a post-hoc account of a failed causal prediction, not a tested one. Two limits bound it. The screen is K562, a single erythroleukemia line and an unfavourable substrate for a state anchored on trophoblast and epithelial biology; we do not report its unperturbed position in the coordinate, so we cannot exclude that the cells had little room to move toward the core in the first place — that baseline should accompany any future use of this screen. And no combinatorial knockdown was run, so the “coordinated pressure” that the account requires is untested rather than supported. *Immune privilege is not inherited.* Of the three properties that motivate the cytotrophoblast — invasion, proliferation, immune evasion — the first is corroborated by orthogonal markers (Note S-C) but the third is not: across 62 colorectal tumors the malignant core fraction is unrelated to immune infiltration (Spearman ρ = +0.13 raw and +0.04 shape-only for T/NK cells, *P* = 0.3 and 0.8; myeloid +0.16 and +0.03), and stratifying by mismatch-repair status — the dominant confounder — leaves it unchanged (MMRd ρ = +0.06, n = 34; MMRp ρ = +0.08, n = 28). Sixty-two tumors exclude a moderate association. This is what the coordinate predicts: Type I interferon is the one target-defining brake cancer does not reproduce, so cancer reaches the address without its immune privilege. Cell-fraction infiltration from dissociated tissue is noisy and infiltration is not exclusion; spatial and functional tests remain open. `scripts/E20_gwps_causal.py`, `scripts/E27_immune_exclusion.py`.



**Note S-K. Definitions, denominators and how the reported percentages relate (reader key).** Three choices govern every percentage in this paper and we fix them here. (i) *Gallery denominator.* The reference holds 154 normal cell types; the two post-cellular bodies (erythrocyte, platelet) are quarantined as a documented low-complexity artifact, leaving **152 candidates** for every assignment. The quarantine is applied *before* the nearest-neighbour step, not after it: the two cells are removed from the candidate pool and the cells that would have gone to them are reassigned to their best surviving neighbour, rather than being dropped from the denominator. Two implementations of this existed in earlier drafts and gave different top-4 boundaries; `scripts/E41_regen_ranking.py` now regenerates the ranking, the boundary table and the per-cancer matrix from one assignment, and Table S2, Table 1, Fig. 3B, 3C, 3E and 3F are all derived from that single artifact. Ranks are quoted out of 154 only where the quarantined cells are deliberately shown. (ii) *Which cells define the measured state, and where each definition is used.* Two definitions appear and we separate them explicitly, because they belong to different sections. The **ranking statistics of §3** — cross-cancer occupancy and its enrichment — are quoted for the single most recurrent cell (cytotrophoblast alone, chance 1/152 = 0.66%) and for the **two-cell core** (cytotrophoblast + megakaryocyte, chance 2/152 = 1.3%). The **within-patient occupancy of §5, which is the paper’s primary assay, uses a three-cell definition** that additionally includes the migrating cytotrophoblast, the trophoblast state adjacent to the first (chance 3/152 = 2.0%); every number in §5 and every control quoting a “~2% chance baseline” is on that definition (`scripts/pelka_perpatient.py`, `E10`, `E12`, `E17`, `E36`). An earlier version of this note stated that the main text used the two-cell definition throughout, which was wrong for §5. Splitting the colorectal effect by member shows why the distinction matters: the trophoblast states carry +16.6 points of the increase and the megakaryocyte −0.8. (iii) *Metric.* The primary statistic is the hard top-1 assignment under the single-cell-scaled calibration. The values a reader will meet are therefore: **23.9%** (top-2 core, main text), **27%** (three-cell core, top-1, the 83,819-cell control panel) and **17%** (three-cell core with the eigenspace rebuilt per calibration, Note S-B; that run has no shape-only variant, so 17% and 17.3% below are *not* the same number). Soft assignment gives 21% core mass for the same cells, and Note S-C reports the two core cells separately (11.1% + 12.8% = 23.9% of cross-cancer votes) rather than as one fraction. These are the same phenomenon under different definitions, not conflicting measurements; the numbers the abstract leads with are the raw ones, and each shape-only value is labelled as such where it appears. (iv) *Shape-only, defined.* “Shape-only” means: subtract from every cell (reference and malignant alike) its own mean across the 44 module scores, `X - X.mean(axis=1, keepdims=True)`, **before** the reference standardization, and recompute the reference mean and SD on the centred reference; projection, normalization and top-1 assignment are otherwise unchanged. This is the convention used in the main text and throughout the figures, and under it the three-cell core falls from 27.1% to **17.3%** (`scripts/E23_purity_fixed_convergence.py`). Note S-A quotes **7–11%** for what reads as the same operation: it is not. There the centring is applied *after* z-standardization, against the unchanged reference statistics (`scripts/E3_referee_controls.py`), which removes more than the magnitude component; on the same 83,819 cells the two orders give 7.4% and 10.7% against 17.3%. Both are honest shape-only variants and both leave the core four- to fivefold over chance, but they are different transforms and we do not average or interchange them. *(v) The two enrichment panels.* The three boundary enrichments are quoted twice in this paper from two different panels, and the difference is the panel, not the result. **Table 1 and Fig. 3F give 18.2-fold, 15.9-fold and 8.4-fold** on the 83,819-cell, 25-cohort panel that every main-text statistic uses. **Note S-Q(i) gives 16.0-fold, 11.5-fold and 6.6-fold**, and 15.8, 12.5 and 6.0 after dropping the complexity-failing members, on the 312,220-cell, 17-cohort pool — that control re-scores every tracked target across the full stored cohorts rather than the capped panel, which is why its baseline differs. Neither supersedes the other: the first is the panel, the second is the control run on the pool where the control can be run. A third 17-cohort figure appears in Note S-V (17.0, 16.3 and 8.7) and is a **different** seventeen — the cohorts with a genuine malignant annotation, on the capped panel — so the two "17 cohorts" in this paper are not the same set and are never compared with each other. *(vi) How the three boundaries were chosen, and what that costs.* The rule stated above — core = the two cells leading both rankings, neighbourhood = recurrent in a majority of the 25 cancers, broader region = recurrent in at least five — does **not** produce the *k* = 2, 4 and 14 that Table 1 and Fig. 3F report. On this panel the rule gives *k* = 2, **3** and **10** (occupancy 23.9%, 33.3% and 69.1%; enrichment 18.2-fold, **16.9-fold** and **10.5-fold**). The reported *k* are round numbers chosen to widen the region in roughly even steps, and we now say so rather than implying they were thresholded; the rule-derived boundaries are *more* enriched, not less, so nothing in the boundary-independence argument depends on the choice. One further caveat belongs with the widest boundary: twelve cell types are recurrent in more than two of the 25 cancers and five are tied at exactly two, so the fourteen-cell membership takes two of those five and **which two is decided by index order, not by the data**. Across all ten ways of breaking the tie the occupancy runs 76.6% to 77.9% and the enrichment 8.3-fold to 8.5-fold; two scripts in this repository broke it differently and reported 76.9% and 77.7% for the same quantity (`E41_regen_ranking.py`, which Table 1 is built from, and `E40_boundary_enrichment.py`). The enrichment is unaffected to one decimal place; the named membership of the broader region is not unique and should not be read as though it were.



**Note S-L. Five further methodological controls.** *(i) Map validity, with the right null.* A cell's nearest neighbour in the nine-dimensional space shares its biological family in **62% of the 152 candidate cell types (95% CI 55–70%)**. The relevant null is not 1/K over families (which would be 7% and is too permissive of the claim) but a family-size-preserving one: permuting the family labels 10,000 times over the same neighbour graph gives a null mean of **14%** (maximum 28%) and a permutation **P = 1 × 10⁻⁴**, which is the floor 10,000 draws can reach rather than a measured value. The 15% previously quoted as "chance" corresponds to this size-preserving rate, not to 1/K. *(ii) Modules share a few genes.* Twelve genes appear in more than one module (ADIPOQ, CALM1, CAMK2A, CCND1, DUSP1, FKBP5, FOXO1, GSK3B, PGR, SGK1, STAT1, XBP1), touching 17 of the 44 modules; 219 gene slots reduce to 206 unique genes. Of the six module pairs with |r| > 0.7, **two are partly definitional** — UPR-ATF6 with UPR-IRE1 (share XBP1) and JAK-STAT with Type I IFN (share STAT1) — while four (ERK/MAPK–NF-κB, UPR-PERK–p53, NF-κB–UPR-PERK, BMP–Hippo) share no genes and are biological. We therefore no longer describe all six as "biologically expected". Recomputing the two definitional correlations on disjoint gene sets requires gene-level reference scores and is listed as a production item. *(iii) Dimensionality.* Three of the nine retained components are the least determined under bootstrap (PC7–PC9, 0.53–0.71 over 500 draws), so we report the result at six as well as nine dimensions: the core fraction is 18.9% (6-D) versus 27.1% (9-D) and the contraction 0.131 versus 0.166 — and in six dimensions the **cytotrophoblast becomes the single top recurrent target**, ahead of the megakaryocyte. The magnitude depends on dimensionality; the anchor does not. *(iv) Leave-one-cancer-out.* Dropping each of the 25 cohorts in turn — not only the four marker-labelled ones — the core fraction stays within **25.4–28.4%** (full panel 27.1%) and the cytotrophoblast's recurrence rank is **2 in every one of the 25 panels** (top-3 in 25/25); the most influential single cohort is B-ALL. No cohort drives the result. *(v) Specification multiplicity.* The shape-only convergence is supported by two of three scalings — in raw module space it falls below every one of 200 null draws, under the paper's reference scaling the permutation P is 0.030, and under the magnitude-removed-reference scaling there is no effect. A Bonferroni correction across the three specifications leaves the raw-module-space result significant and renders the reference-scaled one marginal (P = 0.09); we report the curve rather than a single corrected number. `scripts/E34_methodology_controls.py`.



**Note S-M. Supporting detail for six claims made in the Results.** (i) *Lung correction, full account.* On the reference lung data themselves the cleanly identified malignant cells sit about 18 points above the normal epithelium with which they had been averaged; the diluted arm was ~9% annotated-malignant (Fig. S11-II; scripts/E12_highbaseline_lung.py, scripts/E13_lung_dilution_proof.py). (ii) *Prostate.* A fifteen-patient paired atlas (CELLxGENE) shows no within-patient increase (7/15) at a normal baseline of ~15%; an independent paired cohort (Song et al. 2022, GSE176031, four patients) puts the benign baseline at 2.1% and the tumour arm at 3.8% (paired P = 0.38). The two cohorts disagree about the baseline and neither shows an increase, so prostate is reported as an open negative. (iii) *The fold-change impression.* Across five tissues the added occupancy is flat against baseline on the additive scale (Pearson r = +0.4, P = 0.5) and not significant on a bounded logit scale; the appearance of more room where the baseline is lower comes from dividing by colorectal's 4.7% denominator. (iv) *Fibrotic lung, per cell type.* Within IPF epithelium the shape-only core fraction is highest in transitional AT2 (18.4%) and proliferating epithelium (16.2%) and lowest in ciliated cells (5.2%), all below the malignant lung level of 30.3% (Fig. S12). (v) *Stage independence.* The occupancy increase is already saturated at the earliest stages sampled and is independent of tumour grade and stage in the best-powered cohort. (vi) *Germline alternative.* The germline matches the core on PC1 (cosine 1.00) but diverges across the remaining eight components; by the full recurrence-and-depth ranking the nearest germline cells fall well below the core, stable across 2,000 bootstrap resamples of the 25 cancers.



**Table S9b. The four operational criteria, per module (all 44).** For each module: the signal-sensing component, the transcription-factor output, the negative-feedback element and the receptor field. Modules annotated at construction are marked `recorded`; the remainder were curated from the catalog's own gene lists together with canonical pathway biology by two independent curations that were then reconciled, and carry a confidence grade. A feedback element is named for 38 of 44 modules; six are marked `unclear` rather than filled in — the five for which no induced feedback gene could be named with confidence (AR, PR, PPARgamma, NPAS4, PI3K/PTEN), and Hippo, whose recorded entry (NKD1, NKD2, BIRC2, BIRC3, BIRC5) we withdraw: NKD1 and NKD2 are Wnt antagonists and appear correctly in the Wnt row, so this was a transcription error, and the BIRCs are YAP targets but not feedback inhibitors. **A larger gap follows from this table documenting design intent rather than composition.** Cross-checking against Table S9 — the genes actually scored — the named feedback element is **absent from the scored gene list in 11 of the 38 modules** that name one: AMPK (ULK1), Autophagy (RRAGD/FLCN), Calcium (RCAN1), Cell Cycle (E2F7/E2F8), ER (NRIP1), Hippo, Notch (NUMB/NUMBL/NRARP), SREBP (INSIG1/INSIG2/PCSK9), TR (DIO3), cGAS-STING (TREX1/DNASE2/SAMHD1) and mTOR (DEPTOR/DDIT4/GRB10). For those modules the score cannot be reading the feedback arm at all, which matters because the readability argument in Note S-N(ii) rests on the induced arm. Taken with the two failed validations in Note S-N(iv), the fair description is that the four criteria record how the panel was curated, not a property of the 206 genes as scored — which also helps explain why three curators given the criteria alone built materially different catalogs (Note S-O). The receptor-field column is blank for sixteen modules (AhR, Circadian, ERK/MAPK, HIF, HSF1, Hedgehog, Hippo, NF-kB, Notch, SREBP, UPR-ATF6, VDR, Wnt, cGAS-STING, mTOR, p53), and so is their feedback-mechanism sentence. These are exactly the sixteen whose annotation was `recorded` at catalog construction, where the sensor was entered as a gene list without a prose ligand-class description; all twenty-eight re-curated here carry one. We therefore report receptor-field specificity as documented for 28 of 44 and as asserted at construction but never written down for the other 16 — a gap in the record, which we state rather than fill in retrospectively. Printed in abbreviated form; the complete table, with the mechanism of each feedback loop and the provenance of each annotation, is the attached `tableS9b_module_criteria.csv`.



```{=latex}
\begingroup\scriptsize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.06}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{\dimexpr 0.13\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.19\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.22\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.18\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.21\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.07\linewidth-2\tabcolsep\relax}@{}}
\toprule
\textbf{module} & \textbf{sensor} & \textbf{TF output} & \textbf{feedback element} & \textbf{receptor field} & \textbf{conf.} \\
\midrule
\endfirsthead
\toprule
\textbf{module} & \textbf{sensor} & \textbf{TF output} & \textbf{feedback element} & \textbf{receptor field} & \textbf{conf.} \\
\midrule
\endhead
AMPK & PRKAG1, PRKAG2, PRKAG3 & No direct AMPK transcription factor; output via T… & ULK1 & No ligand class; intracellular AMP/ADP:ATP e… & medium \\
AR & AR & AR & unclear & Androgens: testosterone, 5a-dihydrotestoster… & low \\
AhR & AHR & AHR, ARNT & AHRR, CYP1A1, TIPARP & — & rec. \\
Autophagy & MTOR, RPTOR, RRAGA & TFEB, TFE3, FOXO3 & RRAGD, FLCN & No ligand/receptor; nutrient and energy depr… & medium \\
BMP & BMPR1A, BMPR1B, ACVR1 & SMAD1, SMAD5, SMAD9 with SMAD4 & SMAD6, BAMBI & BMP and GDF ligands of the TGF-beta superfam… & high \\
Calcium & CALM1, CALM2, CALM3 & NFATC1, NFATC2, NFATC3 & RCAN1 & No single ligand class; any cytosolic Ca2+ r… & medium \\
Cell Cycle & CCND1, CDK4, CDK6 & E2F1, E2F2, E2F3 & E2F7, E2F8 & No dedicated ligand class; mitogenic growth … & high \\
Circadian & CSNK1D, CSNK1E, PRKAA1 & CLOCK, BMAL1, NPAS2 & PER1, PER2, PER3 & — & rec. \\
ER & ESR1, ESR2 & ESR1, ESR2 & NRIP1 & Estrogens: 17b-estradiol, estrone, xenoestro… & medium \\
ERK/MAPK & EGFR, FGFR1, FGFR2 & FOS, JUN, JUND & DUSP1, DUSP2, DUSP3 & — & rec. \\
FXR & NR1H4 & NR1H4, RXRA & NR0B2, FGF19 & Bile acids: chenodeoxycholic, lithocholic, d… & high \\
GR & NR3C1 & NR3C1 & FKBP5 & Glucocorticoids: cortisol, corticosterone, d… & high \\
HIF & EGLN1, EGLN2, EGLN3 & HIF1A, EPAS1, ARNT & EGLN1, EGLN2, EGLN3 & — & rec. \\
HSF1 & HSF1 & HSF1, HSF2 & HSP90AA1, HSP90AB1, HSPA1A & — & rec. \\
Hedgehog & PTCH1, PTCH2, SMO & GLI1, GLI2, GLI3 & PTCH1, PTCH2, HHIP & — & rec. \\
Hippo & CDH1, NF2, FRMD6 & YAP1, WWTR1, TEAD1 & unclear & — & low \\
Insulin/FOXO & INSR, IGF1R, IRS1 & FOXO1, FOXO3, FOXO4 & INSR, IRS2 & Insulin and IGF-1/IGF-2 peptide hormones (sy… & medium \\
JAK-STAT & IL6ST, IL6R, OSMR & STAT1, STAT2, STAT3 & SOCS1, SOCS3, CISH & Type I/II cytokines: interleukins, interfero… & high \\
LXR & NR1H3, NR1H2 & NR1H3, NR1H2, RXRA & MYLIP, ABCA1, ABCG1 & Oxysterols and desmosterol (cholesterol exce… & high \\
MR & NR3C2 & NR3C2 & FKBP5 & Mineralocorticoids: aldosterone (cortisol ex… & medium \\
NF-κB & TLR4, IL1R1, TNFRSF1A & NFKB1, NFKB2, RELA & NFKBIA, TNFAIP3, CYLD & — & rec. \\
NFAT & PPP3CA, PPP3CB, PPP3R1 & NFATC1, NFATC2, NFATC3 & RCAN1 & Immunoreceptor (TCR/BCR/Fc) and Gq-GPCR stim… & high \\
NPAS4 & CAMK2A, CAMK2B, CALM1 & NPAS4, ARNT2 & unclear & No ligand class; membrane depolarization and… & low \\
NRF2 & KEAP1 & NFE2L2 heterodimerized with MAFF, MAFG, MAFK & KEAP1 & Electrophiles, oxidants/ROS, heavy metals, t… & high \\
Notch & NOTCH1, NOTCH2, NOTCH3 & RBPJ, RBPJL & NUMB, NUMBL, NRARP & — & rec. \\
PI3K/PTEN & INSR, IGF1R, EGFR & FOXO1, FOXO3, FOXO4 & unclear & Insulin, IGF-1/IGF-2 and growth-factor RTK l… & low \\
PPARα & PPARA & PPARA, RXRA & ACOX1, CPT1A, CYP4A11 & Free fatty acids, eicosanoids, fibrates & medium \\
PPARγ & PPARG & PPARG, RXRA & unclear & Fatty acids, prostanoids (15d-PGJ2), thiazol… & low \\
PR & PGR & PGR & unclear & Progestogens: progesterone, 17-OH-progestero… & low \\
PXR/CAR & NR1I2, NR1I3 & NR1I2, NR1I3, RXRA & CYP3A4, CYP2B6, UGT1A1 & Xenobiotics, drugs, secondary bile acids (li… & high \\
RAR & RARA, RARB, RARG & RARA, RARB, RARG & CYP26A1 & Retinoids: all-trans and 9-cis retinoic acid & high \\
SREBP & SCAP, INSIG1, INSIG2 & SREBF1, SREBF2 & INSIG1, INSIG2, PCSK9 & — & rec. \\
TGF-β & TGFBR1, TGFBR2, TGFBR3 & SMAD2, SMAD3 with SMAD4 & SMAD7, SKIL, PMEPA1 & TGF-beta1/2/3, activins, nodal, myostatin & high \\
TR & THRA, THRB & THRA, THRB, RXRA & DIO3 & Thyroid hormone: T3, generated locally from … & high \\
Type I IFN & IFNAR1, IFNAR2 & STAT1, STAT2, IRF9 & USP18 & Type I interferons: IFN-alpha subtypes, IFN-… & high \\
UPR-ATF6 & ATF6, ATF6B, HSPA5 & ATF6, ATF6B & HSPA5, CALR, PDIA3 & — & rec. \\
UPR-IRE1 & ERN1, ERN2, HSPA5 & XBP1 & DNAJB9, HSPA5 & Unfolded/misfolded proteins in the ER lumen;… & high \\
UPR-PERK & EIF2AK3, HSPA5 & ATF4, DDIT3, ATF3 & PPP1R15A, TRIB3 & Unfolded/misfolded ER-lumenal proteins (ER p… & high \\
VDR & VDR & VDR, RXRA & CYP24A1 & — & rec. \\
Wnt & FZD1, FZD2, FZD3 & TCF7, TCF7L1, TCF7L2 & AXIN2, NKD1, NKD2 & — & rec. \\
cAMP/CREB & PRKAR1A, PRKAR2A, PRKAR2B & CREB1, ATF1, CREM & CREM, PDE4B & Gs-coupled GPCR ligands: catecholamines, pep… & high \\
cGAS-STING & CGAS, IFI16, AIM2 & IRF3, IRF7, NFKB1 & TREX1, DNASE2, SAMHD1 & — & rec. \\
mTOR & INSR, IGF1R, TNFRSF1A & SREBF1, HIF1A & DEPTOR, DDIT4, GRB10 & — & rec. \\
p53 & ATM, ATR, CHEK1 & TP53, TP73 & MDM2, MDM4, COP1 & — & rec. \\
\bottomrule
\end{longtable}
\endgroup
```



**Note S-N. The module panel, evaluated without cancer data.** Everything in this note is computed on the 154 normal cell types; no cancer data enter. *(i) Coherence — the modules are units.* For each module we take the mean pairwise Pearson correlation among its genes across the 154 normal cell types. Median within-module coherence is +0.154 for the 44 modules (median 5 genes each), against +0.014 for 200 random gene sets matched on size and drawn from expressed genes (Mann-Whitney P = 3.2e-15), and against +0.062 for KEGG signalling maps (median 121 genes), +0.058 for the largest Reactome pathways (153), +0.070 for Hallmark (179) and +0.051 for PROGENy (99); the modules beat every collection at P = 1e-4 to 6e-9 while carrying twenty- to thirty-fold fewer genes. *(ii) Readability — why a transcriptional output is required.* Grouping every annotated gene of every module by the role it plays and taking its standard deviation across the 154 normal cell types, the induced arm (feedback elements and TF targets, n = 387) has median SD 0.287 against 0.222 for the structural arm (sensors and cascade components, n = 704), P = 2.8e-05. A pathway is readable from RNA through what it induces, not through the abundance of its receptors and kinases. *(iii) Catalog size.* Random subsets of k of the 44 modules, with the eigenspace rebuilt from the normal cells alone, recover nearest-neighbour biological-family concordance of 38.7% (k=5), 54.5% (15), 61.2% (30), 63.0% (40) and 65.8% (44), and agree with the sampled-catalog cell-to-cell geometry at r = 0.53, 0.78, 0.92, 0.976 and 1.0 — the value at k = 44 being a comparison with itself, not an independent agreement. Ninety-five per cent of the full value is reached at k = 35 on both measures and both are still rising at 44. One caveat on the endpoint: this curve rebuilds the eigenspace from each sampled subset, so its k = 44 value (65.8%, 100 of 152) is not the 62% (95 of 152) quoted in the main text and in Fig. 1C, which is measured in the shipped, fixed nine-dimensional reference (`reference_v03`). The two spaces are near-identical but not identical (cell-to-cell geometry r = 0.964), and the 3.3-point gap is the size of that difference, not a disagreement about the statistic. The catalog therefore sits at the shoulder of saturation. This bounds whether 44 is sufficient; it cannot say whether a larger catalog built on different criteria would do better. *(iv) What we could NOT establish — the negative-feedback criterion.* We tried twice to support criterion (iii) empirically and report both failures. Against the Cancer Gene Census (v103, 740 genes), the genes occupying the feedback position of a module are NOT enriched for cancer-driver status relative to the other genes of the same modules: 17/114 (14.9%) versus 128/638 (20.1%), OR = 0.70, P = 0.25, and significantly in the opposite direction against tighter backgrounds (core genes only, OR = 0.33, P = 0.004). And in (ii) above the feedback genes alone are only marginally more variable than cascade genes (0.242 versus 0.222, P = 0.16) — the induced-arm effect is carried by TF targets. We therefore present negative feedback as the panel's organizing design principle and not as an empirically validated selection rule. `scripts/E37_feedback_cgc_enrichment.py`, `perceptome_validation/scripts/E38_cancerfree_module_justification.py`, `perceptome_validation/scripts/E39_saturation.py`.



**Note S-O. How much of the result depends on this particular catalog, and on which normal cells are compared.** *(i) Independent reconstruction of the catalog.* Three curations were produced by independent large-language-model agents from the four operational criteria and the size target alone, with the incumbent catalog, this paper, and the words trophoblast and megakaryocyte withheld; each shares about half of the incumbent's 206 genes (Jaccard 0.28-0.31) while agreeing with one another far more closely (0.62-0.69). What reproduces is the REGION: the top-ten nearest targets overlap the incumbent's by 5-7 of 10, and five cell types - gastric progenitor, gastric chief, pancreatic acinar, oesophageal apical and megakaryocyte - are in the top ten of all four catalogs. What does not reproduce is the ANCHOR: all three curations place the cytotrophoblast 26th-35th of 152 by cosine to the malignant attractor (0.507-0.639 against 0.955-0.963 for their own winner, the gastric progenitor), and 0 of 3 put it in the top five. The cosine band is flat only at the very top (rank 1 to rank 10 is 0.96 to 0.84) and has fallen to 0.51-0.62 by rank 30, so this is a real separation and not an artifact of ordinal ranking. Between-catalog variance, which the framework previously did not have: core occupancy 12.0 ± 4.3% against 24.7% for the incumbent — the three-cell core on the rebuilt eigenspace this comparison uses, which is the right comparator here and is not the 23.9% two-cell main-text figure (Note S-K). Caveat stated plainly: the three curators are language-model agents of one model family and share its priors (Use of generative AI), which explains their agreement with each other but not their disagreement with the incumbent. *(ii) Established pathway collections.* KEGG (46 signalling maps), Reactome (44 largest pathways), Hallmark (50) and PROGENy (14), each scored through the identical pipeline, reproduce the within-patient occupancy increase at a comparable size (KEGG +15.0 points shape-only, 44/46 patients; PROGENy +12.9; Hallmark +19.2 raw, 46/46; Reactome, the weakest, +6.5 raw and +2.5 shape-only, 33/46, *P* = 0.014; incumbent +9.0 shape-only and +13.9 raw), and under every one of them the largest within-patient shift of colorectal epithelium is toward a trophoblast state. All six numbers in that comparison, the incumbent's included, come from one harness run under a single common setting — pseudobulk loadings with the eigenspace rebuilt per panel — chosen so the collections are scored identically. That is why the incumbent reads +9.0/+13.9 here and +10.5/+15.5 in the main text and Note S-F, where it is scored under the single-cell-matched calibration the tumour data require: same 46 patients, same cells, different calibration (under which raw and shape-only nearly exchange places, 8.81/10.98). The comparison between collections is therefore valid; only the incumbent's absolute value differs from the body of the paper (`perceptome_validation/registry/provenance_notes.md`). None of them reproduces the cross-lineage cytotrophoblast target: their core fractions are 0.05-0.29 times the incumbent's and they place the cytotrophoblast 12th-34th, recurrent in 0-2 of the 25 cancers. The occupancy increase is therefore catalog-independent; the cross-lineage generality of the target is not. *(iii) Which normal cells the comparison uses.* Colorectal cancer arises from crypt-base stem and transit-amplifying cells, not from the tissue at large, so the matched-normal arm must be the right subtype. Scoring whole transcriptomes on 2,000 highly variable genes rather than on the modules, the malignant epithelium of 57 of 62 patients is nearest the proliferating stem/TA population (cE03), and each tumour subtype is nearest its own normal counterpart (Tumor cE03 to cE03, r = 0.74). Restricting the normal arm accordingly: malignant versus normal stem/TA, higher in 35 of 39 patients (+14.8 points); subtype against the same subtype, +14.6 (cE01) and +15.4 (cE03). Normal crypt-base cells occupy the state at 7-9% against 0.1-1% for differentiated enterocytes, so the unrestricted comparison overstates the effect (+19.1 points) and the crypt-base comparison is reported as primary. *(iv) Decoy modules, and what the two robustness percentages are really comparing.* Adding five modules of random genes (44 to 49) requires scoring new modules gene by gene, which the reduced 23-cancer harness does and the full pipeline of Fig. S7-II does not, so no full-fidelity decoy arm exists. Within the reduced harness, where both arms are measured the same way, the decoy arm is not the harsher perturbation: the cytotrophoblast holds the top three by cross-cancer recurrence in 42 of 100 decoy draws against 96 of 200 draws that instead drop five real modules (48%; Fisher *P* = 0.39, rank distributions Mann-Whitney *P* = 0.49), and on the quantitative endpoints the decoy arm is marginally the better of the two (median core fraction 19.0% against 17.7%, enrichment 9.7x against 9.1x). The 42% should therefore not be read against the 84% of Fig. S7-II: those are two rankings of the same perturbation class, not two perturbations. The cytotrophoblast is rank 1 by cosine to the malignant attractor and rank 3 by cross-cancer recurrence when nothing is perturbed, so “stays in the top three” is nearly free in the first ranking and close to a coin-flip in the second: the same 200 five-module drops give 84% by attractor rank and 44% by recurrence rank. `scripts/E2d_module_robustness.py`, `scripts/E2e_module_robustness_fullfid.py`. *(v) Single genes that move the anchor.* Rebuilding the catalog 206 times, once per gene with that gene removed, leaves the cytotrophoblast at rank 1–2 in 202 of the 206 rebuilds (core occupancy 15.4–30.7%). Four genes dethrone it on their own: **EIF4EBP1** (mTOR; cytotrophoblast to rank 12, occupancy 15.4%), **FKBP5** (AR/GR/MR; rank 6), **SLC2A1** (HIF; rank 8) and **CD36** (PPARγ; rank 7). The pre-registered rule treats more than three such genes as a failure of this test, and we report it as one. Deepening every module to 10, 15 or 20 genes tells the same story from the other side: core occupancy stays at or above 11.7% at every depth, while the cytotrophoblast leaves the top three at every depth and under both expansion rules. The magnitude of the convergence is robust to how the modules are populated; the identity of the nearest cell rests on a four-to-five-gene sampling of each module and, within that, on a handful of individual genes. `scripts/E36_colorectal_comparator.py`, `perceptome_validation/scripts/T8_region_vs_label.py`, `perceptome_validation/results/T4b_gene_jackknife.csv`.



**Note S-Q. Two audits: the complexity control applied to the whole convergence region, and gene-symbol recovery per module per cohort.** *(i) The complexity control, every region member.* Fig. 3D ran the transcriptome-complexity decile control for the two core cells and the erythrocyte only. Table S6 shows five further members of the top-14 region at or below the quarantined platelet (1,541 median detected genes) and three at or below the quarantined erythrocyte (890): late spermatids 763, gastric chief 904, pancreatic acinar 992. We therefore ran the same control for every tracked target, on 312,220 malignant cells from 17 cohorts, with the identical QC, cap, seed and assignment rule (`scripts/E43_region_complexity_control.py`). **Six of the fourteen region members fall with complexity exactly as the erythrocyte does** — megakaryocyte (13.2% in the lowest decile to 3.5% in the highest), late spermatids (5.9 to 0.2), gastric chief (7.8 to 1.7), oesophageal apical (5.0 to 2.0), pancreatic acinar (13.4 to 0.5) and erythrocyte progenitors — while eight rise, the cytotrophoblast most steeply (6.4 to 15.7). *What this costs the boundary-independence argument: almost nothing.* Dropping the six and recomputing both the occupancy and the chance baseline, enrichment goes from 16.0× to 15.8× at the core, from 11.5× to **12.5×** at the neighbourhood (it rises, because the two survivors are the two strongest targets) and from 6.6× to 6.0× at the broader region. The region is therefore not an artifact class: widening it does add complexity-suspect cells, but the enrichment does not depend on them. A sharper form of the worry — that four of the five cell types common to all four independent catalogs (Note S-O) are in the quarantine complexity range — stands as stated, and we now say so: the cross-catalog agreement on gastric chief, pancreatic acinar and oesophageal apical is agreement on targets this control rejects. What survives the control **and** reproduces across catalogs is the gastric progenitor (2,506 detected genes, rises with complexity) alongside the cytotrophoblast. *(ii) Gene-symbol recovery.* The pipeline reports per-module coverage but never printed it. Auditing all 44 modules against all 17 cohorts found three genes systematically lost, for two different reasons. **`NKX3-1` was lost in all 17 cohorts** — not a synonym problem but our own de-duplication: the loader stripped any trailing `-<digits>` to undo `var_names_make_unique()`, and `NKX3-1` is a real symbol, so it became `NKX3` and matched nothing. The HPA reference, read from a different file, kept it, so the AR module was scored on five genes in the reference and four in every cancer cohort. **`MB21D1` and `TMEM173` (cGAS-STING) and `ARNTL` (Circadian) are deprecated symbols** (now CGAS, STING1 and BMAL1). The HPA source file uses the current forms, so these genes matched *neither* the reference nor most cohorts: cGAS-STING has in fact been scored throughout on IRF3, IFNB1 and CXCL10 — without either of its sensors — and Circadian on four genes of five. `scripts/gene_names.py` now resolves both classes of name and recovers AR to 84 of 85 module×cohort gene slots (from 67), Circadian and cGAS-STING to 85 of 85 (from 73 and 61); `results/cancer_v6_6_extras/E44_gene_recovery_fix.csv`. *What we withdraw because of this.* cGAS-STING is the one module with inverted sign in Table S7 and the largest cross-cohort variance in Table S4, and the Results sentence that malignant occupants run it hotter is therefore **withdrawn**: a three-gene fragment missing both sensors cannot carry that reading. The rest of the catalog is bounded by the existing module-drop robustness (Fig. S7-II: dropping five modules at random leaves the core fraction at 17.5 ± 2.5% and no single module is load-bearing), so the headline results do not rest on these three modules. **The end-to-end re-derivation has since been run, and it is reported in Note S-W** — including a third name defect that the fix as first written would not have repaired.



**Note S-R. The five isogenic transformation series: pre-registration, per-contrast detail and the gates (Fig. 6).** The design, the two-stage displacement and the cross-origin control are in the main text. Across the nine primary within-series contrasts five reach the pre-stated significance level and four do not, and two further candidate series were set aside by the senescence gate — oncogene-induced senescence in the transduced arm (GSE22573, GSE172029), which the pre-registration defines as uninterpretable rather than negative. The clean contraction of GSE118206 holds under the alternative reference calibration (+0.302, *P* = 0.013 raw; +0.161, *P* = 0.009 magnitude-removed) and is visible in a plain highly-variable-gene view of the same arms (precursors 0.172 apart, transformed 0.023), so the module panel sharpens an effect that is present without it; against 200 expression-matched random gene panels the displacement is specific in raw scores (*P* = 0.005) but not after magnitude removal (*P* = 0.06). Pre-registration `perceptome_validation/PREREG_T9_invitro.md`, frozen before each run; results `T9_RESULTS.md`; scripts `T9_invitro_transformation.py`, `T9c_trajectory.py`, `T9d_activity.py`, `T9e_pancreas.py`.

```{=latex}
\clearpage
```

**Note S-T. The same modules scored by regulatory *output* instead of component *abundance*: the convergence grows, the address does not survive (pre-registered test T7; `perceptome_validation/scripts/T7_activity_scoring.py`).** A module score in this paper is the mean abundance of the module's component transcripts. That measures how *equipped* a cell is for a pathway, not how active the pathway is: activation is largely post-translational and reaches RNA only through induced targets. The distinction is structural rather than a technicality, and a framework named for perception that scores the parts list owes the reader a test at the other level. We pre-registered one, with its three endpoints, its thresholds and a one-sided directional prediction fixed before the data were read.

*Method.* All 44 modules were re-scored by regulatory output. Thirty-four modules whose output is a transcription factor were scored as regulon activity over that factor's DoRothEA A/B regulon (signed weights, 6–755 targets; the union where a module names more than one factor — for example HIF through HIF1A + EPAS1, JAK-STAT through STAT3 + STAT1). ERK/MAPK and PI3K/PTEN, which name no transcription-factor output, were scored by their PROGENy footprints (99 and 100 genes). The remaining eight — AMPK, Autophagy, MR, mTOR, NPAS4, PXR/CAR, TR and UPR-IRE1 — have no usable regulon and were scored by the mean expression of the catalog's own documented target genes, which is the weakest arm of the test and is stated as such. VIPER is the pre-registered primary estimator; ULM is reported as the linear alternative. The 154-type reference and all 83,819 malignant cells were re-scored identically, cohort by cohort on that cohort's present genes; the eigenspace was then rebuilt from the activity-scored 154 × 44 matrix (Kaiser rule as everywhere else) and the whole pipeline re-run under both calibrations, with and without magnitude removal. The erythrocyte/platelet quarantine of Note S-K applies throughout, so the candidate pool is 152 and the three-cell chance baseline is 2.0%. The per-module definitions are frozen in `perceptome_validation/registry/T7_activity_definitions.json`.

| scoring | K | E1 2 / 3-cell | cyto rank | rank-1 cell | E3 Δ (deeper) |
|------------------|---|--------------|-----------|------------------|--------------|
| readiness·pb·raw | 9 | 21.9 / 24.7% | **2** (18/25) | mega 19/25 | +13.9 (43/46) |
| readiness·pb·shape | 12 | 11.0 / 12.7% | 4 (8/25) | gastric prog. 21/25 | +4.9 (40/46) |
| **VIPER·pb·raw** | 8 | 33.3 / **33.3%** | **17** (0/25) | **mega 25/25** | +1.5 (26/46) |
| **VIPER·pb·shape** | 10 | 32.8 / 32.8% | 13 (0/25) | mega 25/25 | **+19.5** (41/46) |
| VIPER·sc·raw | 9 | 43.4 / 43.5% | 19 (0/25) | mega 25/25 | +11.3 (34/46) |
| VIPER·sc·shape | 10 | 21.7 / 21.8% | 17 (0/25) | oocyte 25/25 | +15.7 (41/46) |
| ULM·pb·raw | 8 | 27.6 / 32.8% | 4 (13/25) | oocyte 21/25 | +14.5 (39/46) |
| ULM·pb·shape | 9 | 25.2 / 32.5% | 7 (5/25) | oocyte 22/25 | +8.3 (35/46) |
| ULM·sc·raw | 9 | 14.5 / 23.6% | 5 (6/25) | oocyte 23/25 | +15.3 (43/46) |
| ULM·sc·shape | 8 | 22.2 / 29.7% | 5 (11/25) | oocyte 22/25 | +7.5 (35/46) |

Every endpoint of the test, in the order the pre-registration fixed them. `pb` = pseudobulk loadings, the pre-registered primary; `sc` = the single-cell-matched loadings the tumour data require (Note S-H); `shape` = the per-cell magnitude component removed from reference and cells alike; `mega` = megakaryocyte, `gastric prog.` = gastric progenitor cell. *K* is the number of components the Kaiser rule retains in the rebuilt eigenspace; *E1* is the fraction of the 83,819 malignant cells assigned to the two-cell and three-cell core (chance 1.3% and 2.0%); *cyto rank* is the cytotrophoblast's dense rank among the 152 candidates by cross-cancer recurrence, with that recurrence in parentheses; *rank-1 cell* is the cell that recurrence ranks first, with its own recurrence; *E3* is the within-patient increase in three-cell occupancy across the 46 paired Pelka colorectal patients. Paired Wilcoxon *P* for the ten rows in order: 7 × 10⁻¹³, 1 × 10⁻⁷, **0.20**, 3 × 10⁻¹⁰, 8 × 10⁻⁶, 1 × 10⁻⁸, 3 × 10⁻⁹, 1 × 10⁻⁵, 3 × 10⁻¹² and 4 × 10⁻⁶. Rows as committed: pseudobulk loadings primary, single-cell loadings secondary. The table is `perceptome_validation/results/T7_activity_table_viper.csv` and `..._ulm.csv`.

*Verdict: FAIL, under both estimators.* Two of the three pre-registered conditions hold, with room to spare. Core occupancy at the output level is 33.3% (VIPER) and 32.8% (ULM) against a floor of 11.7% and against 24.7% for readiness scoring in the same harness — the cross-lineage convergence is *stronger*, not weaker, when modules are read through what they induce. The magnitude-removed within-patient deepening is +19.5 points (VIPER, 41 of 46) and +8.3 (ULM, 35 of 46) against a threshold of 5 points and the same sign, and against +4.9 for readiness scoring in this harness — three to four times larger. The test fails on the third condition, and decisively under the primary estimator: the **cytotrophoblast falls to rank 17 of 152 and is a top target in 0 of 25 cancers** (ULM: rank 4, 13 of 25), and the output-level target is the **megakaryocyte in 25 of 25 cancers** (ULM: the oocyte in 21 of 25). The cytotrophoblast is outside the top five in all eight activity specifications. Read straight: measured by regulatory output rather than by parts list, cross-lineage convergence onto one region is stronger and the within-patient deepening is stronger in shape-only terms — and the cell converged upon is not the cytotrophoblast. The pre-registration states what this means and we adopt its words: the convergence is about *equipment* for a state rather than occupancy of it, which is a materially different claim. The identity reported in §3 is therefore scoped twice over — to this catalog (Note S-O) and to readiness-level measurement.

*One raw endpoint collapses, and we do not fold it into the others.* Under VIPER, the within-patient deepening without magnitude removal is only +1.5 points, higher in 26 of 46 patients, *P* = 0.20 — a null, where the same arm reads +13.9 under readiness scoring and +14.5 under ULM. At the output level the Pelka deepening is carried entirely by module shape under one estimator and by both components under the other. The deepening claim of §5 is therefore robust to magnitude removal at the readiness level (where it is primary) but not estimator-independent at the output level.

*The directional prediction was confirmed.* Pre-registered one-sided: the shape-only fraction of E1 should rise relative to raw under activity scoring. It rises from 0.51 under readiness to 0.98 (VIPER) and 0.99 (ULM) — at the output level the raw and shape-only convergences coincide. Magnitude-carried dilution therefore does explain the gap between the two at the readiness level, which means the paper's "largely a shared high-output signal" caveat (Note S-A, Discussion) is a property of readiness-level measurement and not of the phenomenon.

*A caveat that cuts both ways, stated rather than used.* The paper's own complexity control (Fig. S5-II, Fig. S7) shows megakaryocyte assignment rising as detected-gene count falls, and VIPER's rank-based enrichment is sensitive to how many of a regulon's targets are detected; part of the 25-of-25 megakaryocyte may be that same low-complexity artifact appearing at a new measurement level. That is a reason to hold the output-level *target* loosely — it is not a reason to retain the cytotrophoblast, which is absent from the top five in every activity specification tested, under both estimators and both calibrations.

*An independent reproduction.* The *in vitro* series of Note S-R were analysed under the same two estimators as a pre-registered amendment. In the ten-donor paired series (GSE240057) the activity-space attractor's nearest normal cells are the megakaryocyte, the oocyte and spermatids, with the cytotrophoblast absent from the VIPER top eight — the T7 result reproduced in a dataset that contributed nothing to it.

*Which scoring the framework adopts, as the pre-registration requires us to state.* Readiness-level scoring remains primary and is what every number in this paper is computed under: it is the level the catalog was built for, the level at which the four operational criteria were applied, and the level at which the module-level and gene-level orthogonal signatures of §4 were established. Activity-level scoring is adopted as a required robustness test and reported at every endpoint rather than as an alternative primary, for two stated reasons: eight of the 44 modules have no regulon and fall back to mean target expression, so the two levels are not measured with equal quality across the catalog; and the complexity sensitivity above is unresolved at the output level. What changes is the scope of the claims, not the numbers: the increase in occupancy and the cross-lineage convergence are reported as measurement-level-robust, and the identity of the convergent cell is reported as specific to readiness-level measurement, in §3, §4 and the Discussion.

*Provenance of the readiness rows.* They are this test's own internal baseline and are not the main text's figures. T7 rebuilds the eigenspace from the matrix it scores, so under magnitude removal the loadings themselves are rebuilt from the magnitude-removed reference (which is why *K* rises from 9 to 12 in the second row). The main text's shape-only control (+10.5 points; `scripts/E17_deepening_shapeonly.py`, Note S-A) and the panel comparison of Note S-O(ii) (+9.0) instead hold the loadings fixed and remove magnitude at the projection step. The raw rows agree to the decimal (+13.9 here and in Note S-O(ii)); the shape-only rows differ because of that rebuilt-versus-fixed eigenspace, not because of a disagreement about the effect. Every comparison drawn inside this note is between rows of this table, scored in one harness under one setting, which is the only purpose the readiness rows serve. The two-cell readiness figure here, 21.9%, likewise is not the main text's 23.9% for the same reason (Note S-K, `perceptome_validation/registry/provenance_notes.md`).



**Table S-P. Per-cohort malignant-cell purity audit (`scripts/E18_purity_audit.py`).** Selection method and, where the source carries annotations, the annotated-malignant fraction of the arm actually used. Three cohorts are flagged at dilution risk — HCC (68% of the hepatocyte arm from tumor samples), BRCA (24% annotated cancer-epithelial), BCC (12%) — the same failure mode as the lung artifact; we tested directly whether contamination masks convergence, on the two of these three whose labels allow it, and it does not (Note S-I: cleaning HCC and BRCA lowered rather than raised them), so no conservativeness is claimed for this group. NSCLC filters to the census malignant label; the colorectal/gastric/lung *deepening* cohorts use author-annotated tumor compartments (Table above / Fig. S11-II). Seventeen solid/liquid cohorts use source-pipeline malignant files; three germ-cell/sarcoma cohorts are marker-labeled (overstates, bounded by the drop-4 robustness, 0.166→0.174). A CNV-inference floor across all 25 is the one remaining production step.



```{=latex}
\begingroup\scriptsize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.06}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{\dimexpr 0.17\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.43\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.2\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.2\linewidth-2\tabcolsep\relax}@{}}
\toprule
\textbf{cohort} & \textbf{malignant-cell selection} & \textbf{purity class} & \textbf{annotated-malignant \% of arm} \\
\midrule
\endfirsthead
\toprule
\textbf{cohort} & \textbf{malignant-cell selection} & \textbf{purity class} & \textbf{annotated-malignant \% of arm} \\
\midrule
\endhead
HCC & celltype=="Hepatocyte" (all hepatocytes, tumor+normal) & PROXY — AT RISK & 68\% \\
BRCA & pf=None (ALL cells of the breast atlas) & AT RISK (no malignant filter) & 24\% \\
NSCLC & cell\_type=="malignant cell" (census) & author census-malignant & 5\% \\
BCC & disease==BCC \& basal/keratinocyte cell types & PROXY — AT RISK & 12\% \\
Colorectal & ClusterMidway=="EpiT" (author tumor-epithelial) & author-annotated (high) & 64\% \\
Glioblastoma & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Ovarian serous & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
RCC & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Melanoma & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Prostate & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
HNSCC & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Neuroblastoma & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
B-ALL & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
DLBCL & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Squamous lung & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
TNBC & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
PDAC & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
AML & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
MM & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Small-cell lung & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Gastric & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Osteosarcoma & pre-filtered *\_malignant.h5ad (author/inferCNV in source pipeline) & author-annotated (high) & n/a — arm is pre-filtered to malignant in the source \\
Germ-cell (embryonal) & marker-based malignant label & MARKER — overstates (bounded by drop-4) & not available in source \\
Ewing sarcoma & marker-based malignant label & MARKER — overstates (bounded by drop-4) & not available in source \\
Seminoma & marker-based malignant label & MARKER — overstates (bounded by drop-4) & not available in source \\
\bottomrule
\end{longtable}
\endgroup
```



**Note S-U. What the Census cohorts of the panel resolve to (`scripts/verify/census_cohorts.py`, Table S-Q).** Table S1 names fourteen of the twenty-five panel cohorts only as “Census”, which is neither an identifier nor a selection rule. Reading each cohort back out of its stored object gives both, and four things follow that Table S1 does not say. *(i) For eight of the fourteen, “malignant” is a cell-type label and not a malignant annotation.* Six cohorts are selected on the author or atlas annotation `malignant cell`/`neoplastic cell` (glioblastoma, melanoma, dlbcl, squamous lung, tnbc, nsclc). The rest are selected on a label that is not about malignancy: a *normal* lineage label in ovarian serous (`fallopian tube secretory epithelial cell`), prostate (basal, luminal and unspecified prostate epithelium), neuroblastoma (`neuroblast (sensu Vertebrata)`) and B-ALL (B-lineage progenitors); a bare `epithelial cell` in a tumour sample in HNSCC, small-cell lung and gastric; and `abnormal cell` in RCC. Epithelial-in-tumour is the same rule that produced the lung dilution artifact retracted in Fig. S11-II, where a 9%-malignant arm read as normal. `E18_purity_audit.csv` records several of these as “pre-filtered … author-annotated (high)”, which holds for the first group only, and Table S1’s `malignant_label` for HNSCC said `author` and now says `Census`. The direction of the resulting bias is the one argued in Note S-I — normal-adjacent epithelium pulls a tumour arm toward its own normal baseline and so understates deepening — but the purity of these cohorts is a heuristic’s purity, and we say that here instead of claiming an annotation we do not have. *(ii) `is_primary_data` was not filtered, so some cells are counted more than once.* A Census query without that predicate returns a cell once per dataset it was deposited in. Six cohorts carry such rows — prostate 14,095 of 30,000 (47%); squamous lung 6,255 of 30,000 (21%); tnbc 2,980 of 30,000 (10%); neuroblastoma 959 of 30,000 (3%); ovarian serous 898 of 30,000 (3%); glioblastoma 2 cells — and the prostate case is the sharpest: its three dataset identifiers are one atlas plus that atlas’s own Epithelial and Cycling subsets, so cycling prostate cells enter the centroid three times. Recomputing the panel baseline with duplicate rows dropped within each cohort before quality control — everything else identical, including the seed and the same 25 cohorts — moves the core fraction from 17.4% to 17.2% and the anchor from rank 3 to rank 1 by cross-cancer recurrence, while staying rank 1 by cosine to the malignant attractor — that is, towards the paper’s claim and not away from it; the attractor top three is the same list either way. The cell count is unchanged at 121,115 because the 6,000-cell-per-cohort cap of this harness absorbs the removal in every affected cohort, so what changes is which cells the cap draws from. (That cap is the one Fig. S7-II’s robustness harness uses, not the 4,000 of the paper’s own figures: the point here is that the two passes differ in nothing but the duplicate rows.) The double counting is therefore a defect of the Census query and not a driver of the result, and we state it here instead of leaving it implicit in a cell count (`scripts/verify/census_dedup_effect.py`). *(iii) No paired cohort is outside the panel.* The panel’s gastric cohort is Census dataset `19053a82` with donors `Pat01-B`–`Pat24-B`: the twenty-three GSE150290 tumour samples, which are the twenty-three patients of the paired gastric test in §5, so that test is no longer described as held out (Note S3). Checking the other two the same way (`scripts/verify/paired_overlap.py`) gives the same answer: the panel’s colorectal cohort *is* Pelka GSE178341, so the flagship paired test runs on patients whose malignant cells are in the panel by construction; and the panel’s NSCLC cohort — its smallest, 241 malignant cells over 55 donors — carries nine `Kim_Lee_2020` donors, five of them (P0018, P0020, P0028, P0031, P0034) among the ten paired lung patients of E12. The paired design was never an out-of-sample test of *which* normal cell the panel names; what it cancels is lineage and batch, because both arms are one patient on one protocol, and §5 now says that instead. *(iv) Platform and tissue.* Table S1 records “10x” for every Census cohort: glioblastoma pools five datasets over nine assays including microwell-seq, Drop-seq, CEL-seq2 and STRT-seq, and melanoma includes 3,058 Smart-seq2 cells. Ovarian serous spans nineteen tissues including omentum, intestine, peritoneum and ascitic fluid — disseminated sites — where Table S1 says “ovary”; RCC pools two independent studies, one an ipilimumab/nivolumab-resistant series, in a single Table S1 row. Table S-Q carries the identifiers so that any of this can be checked against the Census snapshot the paper pins (2025-11-08).


**Table S-Q. The Census cohorts of the 25-cancer panel, resolved** (`scripts/verify/census_cohorts.py`; `results/census_cohorts.csv`). *rows* is the stored cohort before the paper’s quality control and the 4,000-cell cap; *unique cells* counts distinct `observation_joinid`, so rows minus unique cells is double counting; *primary rows* is how many rows carry `is_primary_data == True`. Blank cells are cohorts whose stored object does not carry the column. Three cohorts (NSCLC, small-cell lung, gastric) were assembled by an earlier pull that kept fewer Census columns.


```{=latex}
\begingroup\scriptsize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.06}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{\dimexpr 0.13\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.06\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.07\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.07\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.05\linewidth-2\tabcolsep\relax}>{\raggedleft\arraybackslash}p{\dimexpr 0.05\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.3\linewidth-2\tabcolsep\relax}>{\raggedright\arraybackslash}p{\dimexpr 0.27\linewidth-2\tabcolsep\relax}@{}}
\toprule
\textbf{cohort} & \textbf{rows} & \textbf{unique cells} & \textbf{primary rows} & \textbf{donors} & \textbf{datasets} & \textbf{Census cell\_type selected on} & \textbf{dataset ids (truncated)} \\
\midrule
\endfirsthead
\toprule
\textbf{cohort} & \textbf{rows} & \textbf{unique cells} & \textbf{primary rows} & \textbf{donors} & \textbf{datasets} & \textbf{Census cell\_type selected on} & \textbf{dataset ids (truncated)} \\
\midrule
\endhead
Glioblastoma & 30000 & 29998 & 22802 & 202 & 5 & malignant cell; neoplastic cell & 774c18c5; 145dcf6a; 21722308; c888b684; 56c4912d \\
Ovarian serous & 30000 & 29102 & 16249 & 41 & 2 & fallopian tube secretory epithelial cell & 44c93f2b; b252b015 \\
RCC & 30000 & 30000 & 30000 & 13 & 2 & abnormal cell & be39785b; bd65a70f \\
Melanoma & 12259 & 12259 & 10992 & 64 & 2 & malignant cell & b617ee1b; 7b20c613 \\
Prostate & 30000 & 15905 & 15040 & 24 & 3 & basal cell of prostate epithelium; epithelial cell of prostate; luminal cell of prostate epithelium & 04557650; b347446f; 881e0e6b \\
HNSCC & 18062 & 18062 & 18062 & 10 & 1 & epithelial cell & b6b5ea88 \\
Neuroblastoma & 30000 & 29041 & 16103 & 22 & 2 & neuroblast (sensu Vertebrata) & 232e432a; afae4675 \\
B-ALL & 30000 & 30000 & 30000 & 18 & 1 & common lymphoid progenitor; early pro-B cell; precursor B cell; pro-B cell & cf83c98a \\
DLBCL & 19136 & 19136 & 19136 & 5 & 1 & malignant cell & a77176ec \\
Squamous lung & 30000 & 23745 & 15458 & 44 & 3 & malignant cell & 1e4214ce; 232f6a5a; 1e6a6ef9 \\
TNBC & 30000 & 27020 & 17314 & 20 & 5 & malignant cell & e500acbf; 5d3fc988; ed880090; 7b20c613; de5416ef \\
NSCLC & 241 &  &  & 55 & 1 & malignant cell & 1e6a6ef9 \\
SCLC & 25000 &  & 25000 &  & 0 & epithelial cell &  \\
Gastric & 8659 &  & 8659 &  & 0 & epithelial cell &  \\
\bottomrule
\end{longtable}
\endgroup
```



**Note S-V. The naming, recomputed without the cohorts whose malignant compartment is a cell-type label (`scripts/E45_genuine_label_subset.py`).** Note S-U reports that in **eight of the fourteen** Census cohorts “malignant” is not a malignancy annotation but a normal lineage label (ovarian serous, prostate, neuroblastoma, B-ALL), a bare *epithelial cell* in a tumour sample (HNSCC, small-cell lung, gastric) or *abnormal cell* (RCC) — the same rule that produced the lung dilution artifact. Because the target is named by cross-cancer recurrence over the panel, a third of the panel's malignant compartments being partly normal cells of the tissue would make the ranking partly a ranking over normal-cell projections, and the direction-of-bias defence is not available: Note S-I tested it and reports it not borne out. One case is sharper still — ovarian “malignant” cells are selected as *fallopian tube secretory epithelial cells* while the region named is a neighbourhood of high-output secretory and progenitor states, so a secretory-epithelial selection could be feeding its own destination. We therefore recomputed the recurrence ranking and the Table 1 boundaries on the **17 cohorts whose malignant compartment carries a genuine author or atlas annotation** (the six Census cohorts labelled *malignant cell* or *neoplastic cell*, plus the eleven author-annotated GEO and atlas cohorts), with identical geometry, standardization, quarantine and top-five rule. This is not the leave-one-cohort-out analysis of Note S-L(iv): the concern is a correlated class of eight, not one cohort at a time. **The naming does not depend on them.** The two-cell core is the same two cells in the same order — megakaryocyte and cytotrophoblast, each a top target in 12 of the 17 (0.71 of cohorts, against 0.72 on the full panel) — and the gastric progenitor stays third, its share rising from 0.52 to 0.65. Thirteen of the fourteen region members are shared (syncytiotrophoblast leaves, Paneth cells enter). On the same boundaries the enrichment is **17.0-fold, 16.3-fold and 8.7-fold** against 18.2, 15.9 and 8.4 on the full panel: the core falls by a little and the two wider boundaries rise. Under the subset's own recurrence rule the majority band is three cells (16.2-fold) and the ≥ 5-of-17 band six cells (13.5-fold). Nor are the heuristic cohorts the ones carrying the convergence: their median core occupancy is 24.0% against 22.5% for the genuine-label cohorts. Ovarian specifically sits above the panel median (core 32.5%, top target the cytotrophoblast at 19.7%), so the circularity is a real contribution rather than an imagined one — but it is removed in the subset above, and the core is unchanged without it. What this does **not** establish is that the heuristic labels are harmless for anything else: they remain a reason to read per-cohort occupancy, and the ordering of cohorts in Table 1, as label-quality-dependent (Note S-I).

```{=latex}
\clearpage
```

**Note S-W. The corrected catalog, re-derived end to end (`scripts/E46_rederive_reference.py`, `E46b_rescore_cohorts.py`, `E46c_rederived_results.py`, `E46d_rederive_colorectal.py`).** Note S-Q(ii) found three genes systematically lost and left the re-derivation undone. It is done here. Every artifact is built beside the shipped one, never over it, so the published numbers stay reproducible; and each stage is gated on reproducing the shipped result with the fix switched **off**, because a re-implementation that cannot reproduce the baseline is measuring itself rather than the fix. *What the fix is.* The HPA source carries the current symbols CGAS, STING1 and BMAL1 (154 rows each) while the catalog scores the deprecated MB21D1, TMEM173 and ARNTL (zero rows), so the **reference** was built without both cGAS-STING sensors and without one of five Circadian genes; after the fix all three modules carry 5 of 5. AR was never affected in the reference — NKX3-1 is present in HPA — which confirms the note's reading that the AR loss was cohort-side only. In the **cohorts** the loader's blanket strip of a trailing `-<digits>` renamed NKX3-1 to NKX3; after the fix all three modules carry 5 of 5 in 24 of the 25 cohorts (HCC has 4 of 5 for AR, a genuinely absent gene). *Two things the audit had not seen.* First, the blanket strip was dropping real symbols far beyond the catalog: in the DLBCL cohort alone it collapsed **575** genuine HGNC symbols (histones H1-0 to H1-8, immunoglobulin V and D segments, snoRNAs), keeping 573 fewer columns. Only one of them, NKX3-1, is in the catalog — we checked all 333 referenced symbols and it is the only one matching the pattern — so module scores are untouched by the other 574; but the detected-gene count that the quality control thresholds on is not, which is why the re-derived panel holds 83,821 cells against 83,819. Second, and more consequential: several cohorts carry **both** symbol forms as separate features with the deprecated one identically zero, and in the Pelka colon object the empty `ARNTL` sits at a lower index than the populated `BMAL1`. Renaming and then dropping duplicates by first occurrence therefore keeps the empty column, so the fix as first written would have scored Circadian on a zero and reported it as repaired. `gene_names.resolve_matrix` now keeps whichever column carries the data. *What moves, and what does not.* Forty-two of the 44 modules are numerically unchanged in the reference. The eigenspace rebuilt with the fix off reproduces the shipped loadings exactly (|cos| = 1.000 on all nine components); with the fix on, PC1 and PC2 are unchanged (1.000, 0.999) and the later components rotate modestly (PC7–PC9: 0.954, 0.957, 0.965). **The Kaiser criterion now selects ten components, not nine** (PC10 eigenvalue 1.011), and we return to that below. On the paper's nine-dimensional frame: occupancy of the cytotrophoblast alone falls from 11.1% to **8.5%** (16.9-fold to **13.0-fold**), the two-cell core from 23.9% to **18.6%** (18.2-fold to **14.1-fold**), the four-cell neighbourhood from 41.9% to 38.4% (15.9-fold to 14.6-fold) and the fourteen-cell region is unchanged at **8.4-fold** (its occupancy is not quoted here because the fourteenth member is tie-dependent, Note S-K(vi)). Thirteen of the fourteen region members are shared (syncytiotrophoblast leaves, erythrocyte progenitors enter). The cross-cancer ranking keeps the **same two leading cells and swaps their order**: the cytotrophoblast is now first at 18 of 25 and the megakaryocyte second at 16 of 25, and the gastric progenitor holds 13 of 25. The one substantial reshuffle is below them — Hofbauer cells rise from 6 of 25 to 14 of 25, from rank 7 to rank 3, which is what moves the gastric progenitor from third to fourth. The colorectal within-patient increase, the paper's primary assay, holds with the same effect: medians **14.1% against 3.5%**, higher in **43 of 46** patients, Wilcoxon *P* = 3 × 10⁻¹² (shipped: 18.9% against 4.7%, 43 of 46, 9 × 10⁻¹³), and the ratio is the same fourfold. Magnitude-removed: 8.8% against 1.6%, 42 of 46, *P* = 6 × 10⁻¹¹. *The one place the anchor does move.* On the frame the corrected catalog's **own** Kaiser criterion selects — ten components rather than nine — the cytotrophoblast falls to 4.3% and the two leading cells become the megakaryocyte and the gastric progenitor: the cytotrophoblast leaves the top two. We keep nine dimensions in this paper for comparability with everything published here, and we report the tenth rather than adopting it silently. This is the **fourth** independent demonstration of the scope already stated in §3: the region and the occupancy are robust, and the name is a property of a catalog, a measurement level and now also a dimension count. *What this costs the paper.* The magnitudes quoted throughout are 3 to 5 points high at the core and the enrichments about a quarter high, while the widest boundary, the within-patient effect and its significance are unchanged; the direction of every reported result is preserved. We have not re-run the paper on the corrected catalog, because doing so would change every number in it for a defect that moves none of its conclusions; what we have done is measure the cost and publish it here.