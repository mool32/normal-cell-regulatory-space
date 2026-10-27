# Perceptome framework — validation pre-registration v1

**Status:** PRE-REGISTERED — no analysis has been run against any endpoint below.
**Author:** T. Spiro
**Date frozen:** 2026-09-19
**Scope:** seven tests (T1–T7) that bound the perceptome's residual priors, establish
competitive baselines, and map the interpretability/fidelity frontier.
**Supersedes:** nothing. Complements `perceptome_generalization_task.md`.

---

## 0. Global rules

These bind every test below. A test that violates any of them is void and must be re-run
from a fresh pre-registration.

### 0.1 Frozen inputs

| Input | Pinned version |
|---|---|
| Normal reference | HPA scRNA-seq, 154 cell types, as used in paper 3 v10 |
| Trophoblast states | Vento-Tormo 2018 |
| CELLxGENE Census | snapshot `2025-11-08` |
| Cancer panel | the 25 cohorts of Table S1 |
| Held-out normal single cells | the 16,312 cells (skin/kidney/breast) of Note S-H |
| Perturb-seq | Replogle 2022 via scPerturb |
| Eigenspace recipe | z-score per module → correlation PCA → Kaiser → 9 PCs |
| Calibration | `single_cell_scaled` for all per-cell projection |

Record SHA-256 of every input file and of the analysis script in
`registry/inputs_<test_id>.json` before the first run.

### 0.2 Verdict vocabulary

- **PASS** — the pre-registered criterion for the hypothesis was met.
- **FAIL** — the criterion was not met; the hypothesis is not supported.
- **SUGGESTIVE** — directionally consistent, does not meet the criterion, explicitly
  under-powered. May not be reported as support.
- **KILLED** — a result that removes a claim from the program. Killed claims are written
  to the kill-list and never silently revived.

Every test appends one line to `registry/verdicts.tsv`:
`test_id | date | verdict | primary_statistic | threshold | notes`

### 0.3 Amendment discipline (Lesson L7)

Thresholds are fixed by this document. An amendment after seeing any endpoint must be
logged in `registry/amendments.md` with timestamp, reason, and direction of effect on the
program's favoured conclusion. **Two amendments in the favourable direction within one
test suite → HALT the suite and escalate for external review.**

### 0.4 Reporting

- Effect sizes lead; P-values follow. Never lead with a P-value whose magnitude is driven
  by feature count (the paper-3 policy).
- Where more than one defensible analytic choice exists, report the full specification
  curve, not the best cell.
- One-sided tests only where the direction was pre-specified in this document.
- Multiplicity: within each test, correct across the specifications that test declares.
  Across tests, no correction (they test different hypotheses) — but this is stated, not
  assumed.

### 0.5 Shared endpoints

Three quantities recur. Defined once here, computed identically everywhere.

- **E1 — core fraction.** Fraction of QC-passing malignant cells whose top-1 nearest of
  152 candidate normal centroids (erythrocyte and platelet quarantined) is cytotrophoblast
  or megakaryocyte. Paper-3 reference values: 21% (two-cell, main text), 27.1%
  (three-cell, referee panel), floor under drop-5 perturbation 11.7%.
- **E2 — target identity.** Rank of the cytotrophoblast in the cross-cancer recurrence
  ranking (number of the 25 cancers in which the cell is among the five most frequent
  top-1 assignments). Paper-3 reference: rank 1–2, recurrence 18/25.
- **E3 — within-patient deepening.** Malignant-minus-matched-normal core fraction,
  per patient, Pelka colorectal (46 patients), paired Wilcoxon + sign test, reported raw
  and shape-only. Paper-3 reference: +15.5 raw, +10.5 shape-only, 43–44/46.

---

## T1 — Established-pathway baselines (decisive)

**Priority: 1. Gates everything else. Run first.**

### Question
Does the curated 44-module catalog do work that a standard, off-the-shelf pathway
collection does not?

### Why it matters
Paper 3 compares against random gene panels (6%) and against reshuffled-module versions of
its own genes (14%). Both are floors, not competitors. No established pathway set has been
run through the pipeline. A reviewer will ask, and the answer decides whether the
framework's methodological novelty survives.

### Inputs
Competing feature sets, each substituted for the 44 modules with the rest of the pipeline
byte-identical:

| Panel | Axes | Genes | Role |
|---|---|---|---|
| `perceptome-44` | 44 | 206 | incumbent |
| `hallmark-50` | 50 | ~4,000 | standard practice baseline |
| `hallmark-50-trim5` | 50 | 250 | size-matched (5 highest-mean-expression genes per set in the HPA reference, selected with no cancer data) |
| `progeny-14` | 14 | ~1,400 | pathway-footprint baseline (activity-flavoured) |
| `reactome-top44` | 44 | variable | 44 largest Reactome signalling pathways |
| `reactome-top44-trim5` | 44 | 220 | size-matched |
| `kegg-signalling-44` | 44 | variable | 44 KEGG signalling pathways |

Trimming rule is fixed here and uses no cancer data: rank a set's genes by mean expression
across the 154 HPA cell types, keep the top 5.

### Procedure
For each panel: score → z-score → correlation PCA → Kaiser (record k; if k ≠ 9, report at
both native k and forced k = 9) → project the 25 cancers → compute E1, E2, E3.

### Pre-registered criteria

Let `X` be any of the six non-incumbent panels.

- **KILLED (modules are not load-bearing)** — some `X` satisfies all three:
  (i) E1(X) ≥ 0.80 × E1(perceptome-44);
  (ii) cytotrophoblast at E2 rank ≤ 2 under `X`;
  (iii) E3(X) within 2.0 percentage points of E3(perceptome-44), shape-only, same sign.
  → The curated catalog is not the source of the result. Paper 3's methodological framing
  must be rewritten around the coordinate-system discipline alone, and the module-selection
  rule is demoted from contribution to implementation detail.

- **PASS (modules load-bearing)** — no `X` reaches cytotrophoblast E2 rank ≤ 3, **or**
  every `X` has E1(X) < 0.50 × E1(perceptome-44).
  → The four-criteria selection rule survives as a contribution and is reported as such,
  with this test as its evidence.

- **REFRAME** — some `X` reproduces convergence (i) but not identity (ii).
  → The convergence is a property of curated regulatory biology generally; the specific
  target is a property of this catalog. Both are reported, separately, and the title claim
  narrows to the target.

- **SUGGESTIVE** — anything else. No claim either way; the test is reported as
  inconclusive with the numbers in full.

### Deliverables
`T1_baseline_panels.py`, a 7 × 3 table (panel × E1/E2/E3), the seven recurrence rankings,
and a verdict line.

### Estimated cost
2–3 agent-days. No new data acquisition.

---

## T2 — The correct null for map validity

**Priority: 2. Cheap. Needed for the paper-3 revision regardless of T1.**

### Question
Is 62% nearest-neighbour family agreement a property of the perceptome, or of any
comparable gene panel?

### Why it matters
The manuscript's first claim ("the map is biologically meaningful before any cancer") is
currently defended against a label-permutation null (14%). That null tests whether the map
has any structure at all, which is guaranteed. It does not test whether *this* map is
special. If a random 206-gene panel reaches ~55%, the claim is uninformative and must be
removed or restated.

### Procedure
Compute family-level nearest-neighbour agreement in the 9-D space for:

1. `perceptome-44` (observed: 62%, 95% CI 55–70%)
2. 1,000 random panels of 206 genes drawn from genes detected in the HPA reference,
   grouped into 44 pseudo-modules of matched size, each rebuilt through the full recipe
3. 1,000 random panels matched additionally on mean expression and detection rate
   (stratified draw — controls for the fact that curated genes are well-expressed)
4. `hvg-2000` → PCA-9 (fidelity ceiling; standard practice)
5. `hallmark-50` → PCA-9
6. Label-permutation null (the incumbent, for continuity): expected ~14%

### Pre-registered criteria

- **REMOVE THE CLAIM** — median of (3) ≥ 55%.
  → 62% is unremarkable. Delete the nearest-neighbour statistic from the manuscript's
  validity argument; retain only the axis-interpretability evidence (external quality
  signature ρ = −0.89 / +0.75, family η² = 0.31–0.70), which no random panel can produce
  a priori.
- **RETAIN WITH CORRECTED NULL** — perceptome-44 exceeds the 95th percentile of (3).
  → Replace the permutation null with the matched-random-panel null in the manuscript.
- **SUGGESTIVE** — perceptome-44 falls between the 50th and 95th percentile of (3).
  → Report as a weak sanity check, explicitly not as evidence.

Report (4) unconditionally as the ceiling, including if it exceeds the perceptome. A map
that loses to HVG-PCA on cell-type recovery is expected and is not a defect — but it must
be stated rather than omitted.

### Deliverables
`T2_map_validity_nulls.py`, a density plot of the three null distributions with the
observed value marked, and the corrected manuscript paragraph.

### Estimated cost
1 agent-day.

---

## T3 — Off-manifold control ("joint exile" vs true proximity)

**Priority: 2. Cheap. Closes the most serious unaddressed criticism of paper 3.**

### Question
Do malignant cells land *near* the cytotrophoblast, or merely *nearest* to it while being
far from every normal cell?

### Why it matters
Assignment is hard top-1 among 152 candidates with no "none of the above" option. Fig. S17
establishes that the core sits in the bottom 1–16th percentile of the normal distribution
on a dozen modules — i.e. outside the normal cloud, in the "all brakes released" corner.
Dysregulation moves cells toward exactly that corner. If malignant cells are systematically
off-manifold, "convergence onto the cytotrophoblast" is better described as "a common exit
direction from the normal manifold whose nearest landmark happens to be the
cytotrophoblast". That is a weaker identity claim but a *more defensible* one, and it does
not depend on trophoblast-ness being real.

### Procedure
For every cell, compute `d_min` = cosine distance to the nearest of the 152 normal
centroids, under the identical pipeline.

Populations:
- **A** — 83,819 QC-passing malignant cells (25 cohorts)
- **B** — 16,312 held-out normal single cells (Note S-H: skin, kidney, breast)
- **C** — matched-normal epithelium from the paired cohorts (Pelka, gastric, lung),
  as the within-patient control

Report: median and IQR of `d_min` per population; the **off-manifold fraction** =
proportion of cells with `d_min` above the 99th percentile of B; and the same split by
cohort and by assigned target (core vs non-core).

### Pre-registered criteria

- **REFRAME REQUIRED** — median `d_min`(A) > 95th percentile of `d_min`(B).
  → Malignant cells are systematically off-manifold. The identity claim is restated as a
  direction-of-exit claim throughout the manuscript. E2 and the recurrence ranking remain
  valid as descriptions of *which landmark is nearest*, not of what the cells resemble.
- **CURRENT CLAIM STANDS** — median `d_min`(A) falls within the IQR of `d_min`(B).
  → Malignant cells are as close to their assigned normal centroid as normal cells are to
  theirs. Proximity is genuine; report this as a new control and cite it where the
  identity claim is made.
- **PARTIAL** — anything between.
  → Report the off-manifold fraction as a quantitative bound on how much of the
  convergence is exile, and qualify the identity claim by that fraction.

Secondary, no threshold pre-registered (exploratory, reported as such): is the off-manifold
fraction higher for core-assigned than non-core-assigned malignant cells? A *lower* value
for core-assigned cells would be positive evidence against the exile reading.

### Deliverables
`T3_offmanifold_control.py`, overlaid `d_min` distributions, per-cohort table, verdict.

### Estimated cost
1 agent-day.

---

## T4 — Module depth: genes per module

**Priority: 3.**

### Question
Does the result depend on the fact that each module is summarised by only 4–5 genes?

### Why it matters
Paper 3 tests robustness to *which modules* (drop-5, leave-one-out) but never to *how
deeply each module is sampled*. 5 genes is a very thin estimator. Two failure modes are
open in opposite directions: too few genes means the module score is noise-dominated; too
many means the regulatory core is diluted by downstream bulk. Neither has been excluded.

### Procedure

**T4a — expansion.** Rebuild each of the 44 modules at 10, 15 and 20 genes, adding
documented pathway members by a fixed rule stated here (TF output modules: highest-confidence
DoRothEA regulon members, confidence A then B; other modules: KEGG/Reactome pathway members
ranked by database annotation confidence, ties broken by mean HPA expression). No cancer
data enters the expansion. Rerun the full pipeline at each depth; compute E1, E2, E3.

**T4b — gene-level jackknife.** Drop each of the 206 genes in turn (206 full pipeline
rebuilds); record E1 and the cytotrophoblast's E2 rank.

**T4c — module-level depth floor.** Rebuild at 3 genes per module (dropping the two
lowest-expressed) as the thin-end control.

### Pre-registered criteria

- **PASS** — at depths 10, 15 and 20: E1 ≥ 11.7% (the established drop-5 floor) **and**
  cytotrophoblast E2 rank ≤ 3 at every depth.
  → The result is not an artifact of thin sampling. Adopt whichever depth maximises E1 as
  the canonical catalog only if pre-registered here — it is not; the catalog stays at 5
  genes and the expansion is reported as robustness.
- **FAIL — fragile** — E1 falls below 11.7% at any depth ≥ 10, or cytotrophoblast leaves
  the top 3.
  → The effect lives in a small number of genes. This must be reported prominently, and
  T4b identifies which genes. A result carried by <10 genes is a different paper.
- **JACKKNIFE FLAG** — any single gene whose removal moves cytotrophoblast out of the
  top 5.
  → Named explicitly in the manuscript as a fragility. More than three such genes →
  treat as FAIL above.

### Deliverables
`T4_module_depth.py`, depth × endpoint table, jackknife influence ranking of all 206 genes.

### Estimated cost
3–4 agent-days (T4b is 206 pipeline rebuilds; parallelise).

---

## T5 — Adversarial reconstruction of the catalog

**Priority: 4. The most expensive and the most decisive about the framework's status.**

### Question
Do the four operational criteria determine the catalog, or does the curator?

### Why it matters
This is the only test that answers "did we take the wrong modules" at the level of the
question rather than by perturbation around the existing answer. Drop-5 explores a tiny
neighbourhood of one point in the space of possible catalogs; it says nothing about the
landscape. If two independent curations following the same criteria converge on the same
target, the criteria carry the result and the perceptome is a method. If they diverge, the
curation carries it and the perceptome is a personal artifact — which is publishable, but
as something else.

### Procedure
Run **three independent curations** (`C2`, `C3`, `C4`). Each is performed by an agent
instance given:

- the four operational criteria verbatim (sensor input; TF output with documented targets;
  negative-feedback element; receptor-field specificity)
- access to KEGG, Reactome and primary literature
- the instruction to screen candidate signalling and metabolic axes and retain those
  passing all four criteria, dropping redundant ones
- the target catalog size range (40–50 modules, 4–6 genes each)

Each is **withheld**: the existing 44-module catalog, Table S9, paper 3 in any version, the
words "cytotrophoblast", "trophoblast", "megakaryocyte", and any statement of the result.
Log the exact prompt in `registry/T5_prompts/`.

Then: for each of C2–C4, compute Jaccard overlap of gene sets and of module names with the
incumbent; rerun the full pipeline; compute E1, E2, E3.

### Pre-registered criteria

- **PASS — criteria-determined** — ≥2 of 3 independent catalogs give cytotrophoblast
  E2 rank ≤ 3 **and** E1 ≥ 11.7%.
  → The selection rule reproduces the finding. This becomes the framework's strongest
  single piece of evidence and should lead the methods section.
- **KILLED — curator-determined** — 0 of 3 place cytotrophoblast in the top 5, while gene
  overlap with the incumbent is < 0.30 (i.e. the criteria genuinely admit very different
  catalogs).
  → The perceptome is one arbitrary point in a large space. The cancer result must be
  restated as conditional on this catalog, and the framework is demoted from method to
  instrument-with-a-provenance.
- **SUGGESTIVE** — 1 of 3, or high overlap (> 0.60) with the incumbent.
  → High overlap means the criteria are tight and the test is weak, not that the framework
  passed; report the overlap and state that the test was uninformative. Consider widening
  the candidate pool and re-running.

Record the **between-catalog variance** of E1 and E3 regardless of verdict. This is the
number that currently does not exist and that the framework most needs: an estimate of how
much the result depends on representational choice.

### Deliverables
Three catalogs as CSV, overlap matrix, three full pipeline runs, a variance estimate for
E1/E2/E3 across catalogs, verdict.

### Estimated cost
5–7 agent-days.

---

## T6 — The interpretability/fidelity frontier

**Priority: 3. Answers a standing conceptual question and may retire the PCA basis.**

### Question
How much fidelity does interpretability actually cost, and is the perceptome currently on
the frontier or under it?

### Why it matters
Two structurally different choices are currently entangled: which 9-dimensional *subspace*
to keep, and which *basis* to use within it. Reconstruction fidelity depends only on the
first. Rotations within a fixed subspace are fidelity-free. If varimax or ICA buys
interpretability at exactly zero fidelity cost — the prior generalization work suggests it
does — then the fixed PCA basis is a presentational choice with no defence, and the
perceptome sits strictly inside the achievable frontier.

### Procedure

**T6a — free rotations (subspace fixed, fidelity invariant by construction).**
Apply to the same 9-D subspace: identity (raw PCA), varimax, quartimax, promax, ICA
(FastICA on the whitened subspace). Verify reconstruction error is numerically identical
across all of them — this is a correctness check on the implementation, not a result.

**T6b — subspace-changing alternatives (fidelity is spent).**
Sparse PCA at λ ∈ {0.01, 0.05, 0.1, 0.5}, NMF-9. Record reconstruction error explicitly;
these do move the subspace.

**T6c — guided expansion.**
Guided autoencoder with the 9 perceptome axes frozen as seed and `k_free` ∈
{0, 5, 10, 15, 20, 25}. Record reconstruction error, seed stability (CKA against the frozen
seeds), and CKA against an HVG-PCA reference at matched total dimensionality.

**Metrics, for every representation above:**

| Metric | Definition |
|---|---|
| Fidelity-R | reconstruction error of the 44-module matrix |
| Fidelity-K | CKA against HVG-2000 PCA at matched k |
| Fidelity-T | family-level nearest-neighbour agreement (as T2) |
| Interp-S | mean across axes of the top-3 modules' share of squared loadings |
| Interp-N | blinded nameability: three independent raters receive only each axis's top-5 signed module loadings, no other context, and return a one-phrase name plus confidence 1–5. Score = mean confidence; report inter-rater agreement. Rubric frozen in `registry/T6_nameability_rubric.md`. |
| Result-survival | E1, E2, E3 under that representation |

### Pre-registered criteria

- **RETIRE THE PCA BASIS** — some rotation in T6a achieves Interp-S and Interp-N both
  above the raw PCA basis (Interp-N by ≥ 0.5 on the 1–5 scale) **and** preserves E1 within
  2 percentage points and cytotrophoblast E2 rank ≤ 3.
  → Adopt that rotation as the canonical presentation of the perceptome axes. Fidelity is
  unchanged by construction, so this is a free gain and there is no argument for keeping
  the unrotated basis.
- **PCA BASIS RETAINED** — no rotation improves both interpretability metrics, or every
  improving rotation breaks the result.
  → Report the frontier as measured and keep the current basis, now with a defence.
- **FRONTIER MAP** (reported regardless of verdict): a Pareto plot of Interp-N against
  Fidelity-K over all representations, with the incumbent marked. If the incumbent is
  strictly dominated, that is stated in the manuscript's limitations, not omitted.

Secondary, pre-registered direction: in T6c, seed CKA is expected to stay ≥ 0.70 and
reconstruction to approach the HVG-PCA ceiling by `k_free` ≈ 15. Deviation in either
direction is reported.

### Deliverables
`T6_frontier.py`, the Pareto plot, the nameability rating sheets, a recommendation on the
canonical basis.

### Estimated cost
3–4 agent-days plus rater time.

---

## T7 — Activity-level scoring (readiness → activity)

**Priority: 5. Exploratory but pre-registered.**

### Question
Does scoring modules by regulatory *output* rather than component *abundance* preserve the
result — and does it strengthen the shape-only residue?

### Why it matters
The readiness/activity dissociation is already documented as structural: module expression
measures how equipped a cell is, while actual activation is post-translational and visible
only through target-gene expression. A framework named for perception currently measures
the parts list, not the percept. This test asks whether the finding survives measurement at
the right level — and predicts, from the dissociation itself, that the magnitude-carried
component should shrink and the shape component should grow.

### Procedure
Rescore all 44 modules by output:
- modules whose output is a TF → regulon activity via VIPER on DoRothEA regulons
  (confidence A/B), or mean expression of documented targets where no regulon exists
- pathway-footprint modules → PROGENy-style footprint scores where available

Rebuild the eigenspace from the activity-scored 154 × 44 matrix. Rerun the full pipeline.

### Pre-registered criteria

- **PASS** — E1 ≥ 11.7%, cytotrophoblast E2 rank ≤ 3, and E3 same sign and ≥ 5 percentage
  points (shape-only).
  → The result is not an artifact of readiness-level measurement. Activity scoring becomes
  a candidate canonical scoring method pending T6-style fidelity comparison.
- **FAIL** — any of the three breaks.
  → The finding is specific to readiness-level measurement. This is informative and must be
  reported: it would mean the convergence is about *equipment* for a state rather than
  occupancy of it, which is a materially different claim and changes the Abelev reading.
- **DIRECTIONAL PREDICTION** (pre-registered, one-sided): the shape-only fraction of E1
  rises relative to the raw fraction, compared with readiness scoring. Confirmed →
  supporting evidence that readiness/activity dilution was suppressing the shape residue.
  Refuted → the dissociation does not explain the magnitude dominance, and that explanation
  is withdrawn from the program.

### Deliverables
`T7_activity_scoring.py`, side-by-side endpoint table (readiness vs activity), verdict,
and an explicit statement of which scoring the framework adopts going forward.

### Estimated cost
4–5 agent-days.

---

## Execution order and gating

```
Week 1     T1 ──────────── decisive; gates the framing of everything downstream
           T2 ─┐  run in parallel, cheap, needed for the paper-3 revision
           T3 ─┘  regardless of T1's outcome

Week 2     T6 ──────────── independent of T1; may change canonical presentation
           T4 ────────────

Week 3-4   T7 ────────────
           T5 ──────────── most expensive; run last, informed by T1/T4 results
```

**Gate after T1.** If T1 returns KILLED, stop and re-plan: T4, T5 and T7 all assume the
catalog is the object of interest, and that assumption would be void. T2, T3 and T6 remain
valid and should still be completed for the manuscript revision.

**Gate after T3.** If T3 returns REFRAME REQUIRED, the manuscript's identity claim changes
before T5 is designed, because T5's E2 criterion would need restating in exit-direction
terms.

---

## What each test can kill

| Test | Claim at risk |
|---|---|
| T1 | that the four-criteria module catalog is a methodological contribution |
| T2 | that "the map is biologically meaningful before any cancer" is evidence rather than a sanity check |
| T3 | that malignant cells *resemble* the cytotrophoblast, as opposed to exiting the normal manifold in its direction |
| T4 | that the result is a property of the modules rather than of a handful of genes |
| T5 | that the perceptome is a method rather than one curator's artifact |
| T6 | that the fixed PCA basis is the right presentation of the axes |
| T7 | that the convergence is about occupancy of a state rather than equipment for one |

A suite in which nothing is killed should be treated with suspicion, not satisfaction.

---

## Registry files created by this suite

```
registry/
  inputs_<test_id>.json        SHA-256 of every input and script
  verdicts.tsv                 one line per test
  amendments.md                timestamped, with direction of effect
  T5_prompts/                  exact curation prompts, for auditability
  T6_nameability_rubric.md     frozen before rating
```
