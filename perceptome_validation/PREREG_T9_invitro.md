# PREREG T9 — transformation *in vitro* with a known precursor

**Frozen 2026-10-19, before any expression value was read.** Written at A. V. Gudkov's request: *«можно найти данные по транскриптомам клеток которые инвитро превратились в опухолевые… и там известны предшественники… это проверка что мы видим те же переходы»*.

Data were downloaded before this file was written (`data/invitro/`), but only **sample titles and platform ids** were inspected, which define the arms. No expression matrix was opened, no score computed, no projection run. Anything below that reads as a prediction is a prediction.

---

## 0. Why this test exists

Every within-patient result in paper 3 compares a malignant cell with a normal cell we *believe* to be its counterpart. The belief is the weak link: the cell of origin is never known with certainty, and a mismatched comparator can manufacture the effect. *In vitro* transformation removes the belief — the precursor is in the dish, and the same dish, with defined genetic hits, produces the transformed state.

It also tests something the patient data cannot: **whether two different normal epithelia given the same oncogenic hits end up in the same place.** That is the paper's headline claim, run as an experiment rather than inferred from a cross-sectional snapshot.

---

## 1. Series, arms, and what each one is for

| id | precursor arm (n) | transformed arm (n) | control / specificity arm | platform | role |
|---|---|---|---|---|---|
| **GSE118206** | prostate basal (3), NHBE bronchial (3) | prostate-PARCB (6), NHBE-PARCB (5) | factor-dropout ACB/PACB/ARCB (13); cancer lines (6) | RNA-seq (recount3, SRP156528) | **primary**: two epithelia, same five hits, both precursors in series |
| **GSE6885** | HME (4), BPE (5) | HMLER (6), BPLER (6) | — | GPL570 array | second two-origin series, mammary |
| **GSE172029** | NC, MB, MEL ±DOX, WT and KO | BRAF-V600E induced (DOX) at each stage | CONTROL (no DOX) at each stage | RNA-seq | one oncogene, three precursor stages, only one competent |
| **GSE231776** | BJ fibroblast EP (**1**) | +hTERT (3), +SV40 (3), +HRAS (3) | — | RNA-seq | ordered ladder; non-epithelial ⇒ is convergence epithelium-specific |
| **GSE22573** | organotypic epidermis LacZ (8) | CDK4+HRAS (8) | CDK4+HRAS **+U0126** (2) — pharmacological reversal | GPL571 array | best-powered isogenic pair; reversal arm |
| **GSE112369** | normal stomach 1a–2c (6), normal-counterpart (2) | stomach CDH1/TP53 DKO (4), CDH1 KO (1); colon CDH1/TP53 DKO (1) | engineered control stomach (1), colon (1); R-spondin ± | GPL15207 array | normal → engineered → **real patient tumour organoid** (36) on one platform |
| **GSE57965** | **absent** | genome-engineered organoid (6), modified adenoma (3) | adenoma (5), CRC organoid (9) | GPL15207 array | **adversarial**: authors report five mutations do *not* reach a CRC profile |
| **GSE3151** | HMEC + GFP (10) | RAS/MYC/E2F3/SRC/β-cat (7–10 each) | GFP is the control | GPL570 array | 18 h activation = first step only, direction not endpoint |

GSE57965 has no normal-organoid arm; it is scored on the adenoma→engineered→CRC axis only and is never used for a precursor-paired statistic.

---

## 2. Pipeline (fixed; no free parameters)

Every series goes through the paper's own frame, unchanged:

1. **Gene level.** Array: probe→symbol via the GEO platform annotation; collapse duplicate symbols by **max** (mean reported as sensitivity). If the matrix is log-scaled (max value < 30) it is exponentiated first.
2. **Calibration.** Each sample is treated exactly as an HPA cell type is treated: linear expression → scale to a total of 1e4 → `log1p`. This is `prereg_lib.hpa_gene_level()['H_sc']` applied to a new row.
3. **Module scores.** Mean over the module's genes present in that series, minimum 2 genes, else 0 — `score_ref` semantics, identical to the reference.
4. **Projection.** `mu, sd, refPn = ref_frame(R_sc, L)` from the 154 HPA types with the bundled loadings `L`; `project(S, mu, sd, L)` → unit 9-D vector. Erythrocyte and platelet stay quarantined (152 targets).
5. **No refitting of anything.** The eigenspace, the loadings, `mu`, `sd` and the reference centroids are the shipped ones.

**Shape-only** variants of every statistic use the same within-sample centering as the paper (`X − X.mean(axis=1)`) applied to reference and samples alike.

---

## 3. Endpoints, in order of weight

### P1 — cross-origin convergence (label-free, no normal cell named)

For a series with two distinct precursors A and B given the same hits:

```
Δ = cosdist(centroid(prec_A), centroid(prec_B)) − cosdist(centroid(tr_A), centroid(tr_B))
```

Positive Δ = the two origins ended closer than they started. Significance by permuting the precursor/transformed label within origin, 10 000 draws, one-sided.

Applies to **GSE118206** (prostate basal vs NHBE) and **GSE6885** (HME vs BPE). This is the paper's Fig. 2 claim done experimentally, and it names no normal cell at all.

### P2 — movement toward the malignant attractor (label-free)

Attractor **A** = mean of the 25 per-cancer malignant centroids (unit 9-D), computed once from `perceptome_validation/cache/panel25.h5ad` under the shipped frame and frozen to `results/T9_attractor.json` before any *in vitro* sample is projected. It is built from cancer cells only and carries no normal-cell label.

```
δ = cos(transformed, A) − cos(precursor, A)        per series, per origin
```

Reported raw and shape-only, with a paired test across samples where n allows and as a centroid difference where it does not.

### P3 — region occupancy (labelled; three frozen definitions, primary declared now)

Top-1 nearest normal among the 152, precursor vs transformed. Region membership, all three reported, **primary = R2**:

- **R2 (primary)** — the 14 most recurrent cell types of Table 1, i.e. the paper's own "region". Chosen because it pre-exists this test.
- **R2′** — R2 minus the six members that fall with transcriptome complexity (Note S-Q). Sensitivity; depth is not limiting in bulk, so R2 is primary.
- **R1** — the five cell types in the top 10 of all four catalogs (gastric progenitor, gastric chief, pancreatic acinar, oesophageal apical, megakaryocyte; T8). The catalog-robust region.
- **R3** — the three-cell core (cytotrophoblast, migrating cytotrophoblast, megakaryocyte), the unit §5 actually uses. Reported so that the *in vitro* test speaks to §5 in §5's own currency, **and** so that the trophoblast-specific share of any movement is visible rather than hidden.

### P4 — direction along the released-brake signature

Module-level Spearman ρ between the transformed−precursor module-z delta and the target-defining suppression signature of Fig. 4. Secondary; descriptive.

---

## 4. Gates, with numbers

| gate | rule | consequence |
|---|---|---|
| **G1 coverage** | ≥ 80% of the 203 module genes present; every module ≥ 2 genes | series excluded and reported; per-module coverage table shipped |
| **G2 OIS** | transformed arm fails if proliferation delta < 0 **and** CDKN2A delta > 0. Proliferation = mean z of MKI67, CCNB1, CCNA2, TOP2A, PLK1, BUB1; senescence = CDKN2A, CDKN1A, SERPINE1, IL6, CXCL8, GLB1 | arm reported as senescent, excluded from the headline, kept as a specificity control |
| **G3 baseline ceiling** | if the precursor's top-1 is already a region member, the series **cannot** show entry | reported as "precursor already inside" — neither support nor refutation. This is the most likely uninteresting failure: culture is proliferative and stressed |
| **G4 platform** | every statistic is a within-series paired delta. Arrays and RNA-seq are never pooled; no absolute cosine is compared across series | cross-series statements are forbidden by construction |

---

## 5. Decision rule, fixed now

- **SUPPORTS the convergence claim** if P1 > 0 with permutation *P* < 0.05 in **GSE118206**, and P2 > 0 in at least **5 of the 7** precursor-paired series that pass G1–G3, while the non-oncogenic arms (U0126 reversal, no-DOX controls, engineered control organoids) do not move.
- **REFUTES it** if P2 ≤ 0 in the majority of gate-passing series, or if the control arms move as much as the transformed arms.
- **BOUNDARY** otherwise: report which lineages move and which do not, and say so in those words.

Pre-committed expectations, so that the result cannot be re-read afterwards:

1. **GSE57965 is expected to be negative or weak.** Matano's own conclusion is that five engineered mutations do not reach a CRC expression profile. If our coordinates nevertheless place the engineered organoids with the CRC organoids, that is a strong result for the module space; if they sit with the adenomas, that is the honest boundary and it will be reported as such.
2. **GSE3151 tests direction only.** Eighteen hours of oncogene activation is not transformation. A positive δ there says the first step points the right way; a negative one does not refute the endpoint claim.
3. **GSE118206's endpoint is small-cell neuroendocrine**, i.e. transdifferentiation driven by a fixed five-factor cocktail. Convergence of the two origins there is partly guaranteed by the shared cocktail. This is stated in advance and will be stated in any write-up: the series proves that two epithelia can be driven to one state, not that spontaneous cancers do so.
4. **GSE231776's precursor arm is n = 1.** No confidence interval is available for its δ; it is a direction, not an estimate.

## 6. What the result licenses

A positive result licenses one sentence in the paper: transformation of human cells with defined hits, from known precursors, moves them toward the convergence region — measured against the precursor itself, so no cell-of-origin assumption is involved. It does **not** license any statement about which normal cell the region "is"; P3/R3 is reported for transparency about the trophoblast share, not as evidence for an anchor.

A negative result is reported in the paper regardless of what the rest of the analysis shows, with the series named.

---

*Scripts: `scripts/T9_invitro_transformation.py` (one file, writes `results/T9_*.csv` and `results/T9_verdict.json`). Inputs hashed into `registry/inputs_T9.json`. Downloaded data live under `data/invitro/<GSE>/` and are read by path, never executed.*

---

## Addendum A1 — the second two-origin series (frozen 2026-10-19, after the first run, before any value of these series was read)

The first run left P1 resting on a single series, GSE118206, whose stated weakness is that both origins
receive **one fixed five-factor cocktail**, so part of the convergence is built into the design. This
addendum adds the series that removes exactly that objection.

**GSE17941** (Hirsch/Iliopoulos/Struhl 2010, GPL9188, 49 samples, one lab, one platform) carries two
different human cell types transformed by **different** hits:

| origin | precursor arm (frozen) | transformed arm (frozen) | non-transforming control |
|---|---|---|---|
| MCF-10A, mammary epithelium | `pBabe, TAM, 24 hour` + `36 hour` (vector + the same drug) | `Er-src, TAM, 24 hour` + `36 hour` | `Er-src, EtOH` (construct, no drug) |
| BJ foreskin fibroblast | `EL` (hTERT + LT, immortalised, not transformed), n = 3 | `ELR` (hTERT + LT + HRAS), n = 3 | `EH` (hTERT + HRAS, does not transform), n = 3 |

If the two arms converge here, the shared-cocktail objection does not apply: the mammary cells are
transformed by Src activation and the fibroblasts by RAS on an LT background.

**Pre-specified statistics**, identical in form to the first run:

- **P1** Δ = cosdist(MCF-10A precursor, BJ EL) − cosdist(MCF-10A Er-src/TAM late, BJ ELR), permutation *P*, raw and shape-only.
- **P1c, the specificity twin**: the same quantity computed on the two *non-transforming* arms
  (`Er-src, EtOH` against `EH`). A convergence that also appears here is a culture/handling effect, not transformation.
- **P2** per origin, and **P3/P4** as before.

**Gate G5, new and numeric.** GSE6885 failed to be informative because its two precursors started 0.024
apart. MCF-10A and BJ are also both immortalised lines. **If d_prec < 0.05 the series is declared
uninformative for P1 before its Δ is interpreted**, exactly as GSE6885 was.

**GSE65885** (the RNA-seq version of the same two systems, GPL11154; `EL_rna`, `EH_rna`, `ELR_rna`, and the
Er-Src TAM/EtOH series) is run with the same arms as a cross-platform replication. It has **n = 1 per arm**,
so it yields a direction, never an estimate, and no permutation test.

**Pre-committed expectation.** MCF-10A at 24–36 h of Src activation is early transformation, and BJ-ELR is a
fully transformed fibroblast; the two systems differ in how far along they are. A null here would therefore
not refute GSE118206, and will be reported as what it is: two immortalised lines, different hits, no
convergence detected.

---

## Addendum A2 — the third series (frozen 2026-10-19, after A1, before any value of these series was read)

A1 removed the replication and left the first series' stated weakness untested: **does the coordinate read
transformation, or does it read the five-factor cocktail?** A2 tests exactly that, and adds a lineage where
the biology supplies its own negative control.

### 3a. GSE240057 — the PARCB time course, 10 donors, paired within donor (PRIMARY)

Each of ten donors contributes, in one series and on one platform:

| arm | sample | what it is |
|---|---|---|
| `prec_basal` | TP1, "human prostate basal cells, 0wks" | the primary precursor, n = 10 |
| `ctrl_transduced` | TP2, "PARCB transformed prostate organoid, 0wks" | **all five factors already delivered, not yet a tumour**, n = 10 |
| `tr_early` | 4–5 wks | n = 8 |
| `tr_transition` | 4.5–9 wks | n = 14 |
| `tr_endpoint` | last timepoint per donor, 8–18 wks | n = 10 |

**The decisive statistic, with its reading fixed now.** Let
δ_endpoint = cos(endpoint, A) − cos(basal, A) and δ_cocktail = cos(transduced 0wks, A) − cos(basal, A),
both paired within donor, and f = δ_cocktail / δ_endpoint.

- **f ≥ 0.7** → the displacement is delivered by the construct, before a tumour exists. The coordinate is
  reading the cocktail, GSE118206's convergence is largely an artefact of it, and the paper must say so.
- **f ≤ 0.3 together with a positive stage ordering** → the coordinate is reading progression, and
  GSE118206 stands as a statement about transformation.
- **0.3 < f < 0.7** → ambiguous; reported as such, with f.

**Supporting statistics:** paired Wilcoxon on δ_endpoint across the ten donors, raw and shape-only; per-donor
Spearman between stage index (basal 0, transduced 1, early 2, transition 3, endpoint 4) and cos to A, reporting
how many of the ten are positive; and the stage at which each donor's top-1 first enters R2.

This series cannot test P1 — one origin — and is not asked to.

### 3b. GSE148714 — the same fusion, two sources, only one transforms (SECONDARY)

CD34+ HSPC transduced with MLL-AF4 or MLL-AF9; the paper's own result is that **only cord-blood HSPC are
transformed by MLL-AF4**. Arms: `Ctrl` (2, CD34+ HSPC, normal karyotype — the two control samples are not
labelled by source, so they serve as one pooled precursor), `CB_MLL-AF4` (2), `BM_MLL-AF4` (2),
`CB_MLL-AF9` (2), `BM_MLL-AF9` (2).

**Pre-specified:** δ(CB-AF4) > δ(BM-AF4). BM-AF4 is the biological non-transforming twin that GSE118206
lacked and GSE17941 failed. If BM-AF4 moves as far as CB-AF4, the coordinate does not read transformation
in this lineage. n = 2 per arm: a direction, never an estimate, no permutation test, and a null is not a
refutation at this n.

### Pre-committed expectations

1. Some displacement at `ctrl_transduced` is expected — MYC and myrAKT1 are delivered at week 0 and act
   immediately. The question is the **fraction**, which is why f is thresholded in advance.
2. GSE240057's endpoint tumours are the same small-cell neuroendocrine endpoint as GSE118206, so this series
   cannot remove the transdifferentiation caveat; it can only separate cocktail from progression.
3. Haematopoietic precursors are far from the convergence region in the published frame, so 3b has room to
   move in either direction.

---

## Addendum A3 — activity scoring on GSE240057 (frozen 2026-10-19, after A2, before any activity score of this series was computed)

A2 found that 63–70% of the continuous displacement in GSE240057 is present at week 0, in transduced cells
that are not yet tumours. Two readings are open: the coordinate is reading the **abundance of module
component transcripts** that oncogene delivery perturbs directly, or it is reading a real change of state
that precedes the tumour. T7 already showed that re-scoring every module by regulatory **output** instead of
component abundance strengthens convergence and raises the shape-only share from 0.51 to 0.99. If the week-0
jump is an abundance artefact, activity scoring should shrink it.

This also answers Gudkov's "вторая производная" objection on the same data, which is why it is run now.

### Scoring

Exactly T7's definitions, unchanged (`T7_activity_definitions.json`): a DoRothEA A/B regulon for every
module with a primary TF and ≥ 5 targets, PROGENy footprints for ERK/MAPK and PI3K/PTEN, mean of the
catalog's documented `activity_genes` for the rest. Estimator **VIPER primary** (T7's pre-registered primary),
ULM reported alongside. The reference and the 25-cancer attractor are rebuilt in activity space from the
**cached** T7 matrices and frozen to `T9d_attractor_<estimator>.json` before GSE240057 is scored.

### Primary statistic and its reading, fixed now

f_act = δ_cocktail / δ_endpoint in activity space, paired within donor, as in A2.

- **f_act ≤ 0.3** → the week-0 jump was component-transcript abundance. The activity coordinate reads
  transformation, and the paper should report the activity version.
- **f_act ≥ 0.7** → the jump is real at the activity level too: oncogene delivery genuinely moves the
  regulatory state before a tumour exists. That is a finding, not an artefact, and it must be stated as the
  limit of what this test can show.
- **0.3 < f_act < 0.7** → partial; report the number against the readiness value of 0.699.

### Secondary, all pre-specified

- Region-membership ladder by stage in activity space (readiness gave 0/10 → 1/10 → 6/10 → 12/15 → 10/10).
- Per-donor stage ordering ρ (readiness gave +0.77, 10/10 positive).
- δ_endpoint shape-only (readiness gave −0.047, 3/10). T7's directional result predicts this should improve.

### Control C3 — drop the delivered genes

PARCB delivers shRNAs against **TP53** and **RB1** and overexpresses **myrAKT1**, MYC and BCL2. Three of
those appear in the catalog's scored genes: RB1 in Cell Cycle, AKT1 in PI3K/PTEN, TP53 in p53 (MYC and BCL2
are in no module). Re-score with those three genes removed, under ordinary readiness scoring, and recompute f.

- **f falls below 0.5** → the week-0 jump was carried by the delivered transcripts themselves.
- **f essentially unchanged** → it is not; three of 219 scored slots cannot explain it, and the jump is a
  broader change of the expression profile.

### Pre-committed expectation

Only 3 of 219 scored gene slots are delivered genes, so I expect C3 to change f very little. The activity
result I genuinely cannot predict: T7 showed activity scoring raises the shape share of the signal, which
argues f_act should fall, but oncogene delivery also acts through real signalling at week 0 (myrAKT1 is a
constitutively active kinase), which argues it should not. Both outcomes are informative and both are
reportable.

---

## Addendum A4 — GSE222990, the fourth series (frozen 2026-10-19, after A3, before any value of this series was read)

Xu et al. 2024, *Nat Commun* 15:818 (PMID 38280869). This is the design T9 has been missing: **two different
PRIMARY, non-immortalised human cell types** given the **same four drivers**, with a **stated
non-transforming driver arm** and a **culture-only arm**. All 28 sample titles were read from the GEO record
before this was written; no expression value was opened.

| arm | n | what it is | platform |
|---|---|---|---|
| `fresh_acinar` | 5 | freshly sorted primary human acinar cells | GPL24676 |
| `fresh_ductal` | 5 | freshly sorted primary human ductal cells | GPL24676 |
| `acinar_organoid` | 6 | the same acinar cells in organoid culture, **no drivers** | GPL24676 |
| `acinar_K_organoid` | 5 | acinar organoid + **KRAS alone**; the authors state KRAS alone does not transform them | GPL24676 |
| `acinar_KPTS_organoid` | 3 | acinar + four drivers; forms PDAC in vivo | GPL24676 |
| `ductal_KPTS_cultured` | 4 | ductal + four drivers; forms PDAC in vivo | **GPL21290** |

### Gate G6 — the platform split, declared before analysis

The ductal transformed arm is on a **different platform and a later submission** (GSM7947043-46) from
everything else, `fresh_ductal` included. Consequences, fixed now:

- **For P1 this is conservative.** Any platform effect pushes the two transformed arms apart, so a positive
  cross-origin Δ is a lower bound, not an inflation.
- **The ductal per-origin delta is confounded** in an unknown direction and is reported as such, never as
  evidence on its own.
- **The acinar ladder is entirely within GPL24676** and is therefore the clean part of this series.

### Primary statistics

**P1 (cross-origin, as in GSE118206).**
Δ = cosdist(fresh_acinar, fresh_ductal) − cosdist(acinar_KPTS, ductal_KPTS), raw and shape-only,
10 000 label permutations within origin. Gate G5 applies: if the two precursors start < 0.05 apart the
series is declared uninformative for P1, as GSE6885 was.

**The acinar ladder, culture held constant — the part GSE118206 and GSE240057 could not supply.**
With `acinar_organoid` as the baseline (everything cultured, only the drivers differ):

- δ′(KPTS) = cos(acinar_KPTS, A) − cos(acinar_organoid, A)
- δ′(K) = cos(acinar_K, A) − cos(acinar_organoid, A)
- **k′ = δ′(K) / δ′(KPTS)** — the share delivered by a driver that does *not* transform.

And with `fresh_acinar` as the baseline, to size the culture step itself:

- **g = δ(acinar_organoid) / δ(acinar_KPTS)** — the share delivered by culture alone.

Thresholds, fixed now and in the same convention as A2 and A3: **≥ 0.7 = that arm explains it;
≤ 0.3 = it does not; between = partial, report the number.**

### Decision rule

- **SUPPORTS**: P1 Δ > 0 at permutation *P* < 0.05, **and** δ′(KPTS) > 0, **and** k′ ≤ 0.7.
- **REFUTES**: δ′(KPTS) ≤ 0, **or** k′ ≥ 0.7, **or** g ≥ 0.7 with δ′(KPTS) null — i.e. the move is culture,
  or a non-transforming driver does it too.
- Everything else is reported as partial, with the numbers.

### Pre-committed expectations

1. Fresh cells versus organoids is a large state change in its own right; I expect *g* to be substantial.
   That is exactly why δ′, measured from the organoid baseline, is the primary and not δ.
2. KRAS alone is a real oncogene with real signalling output, so I do **not** expect k′ ≈ 0. The question is
   whether it reaches the threshold, and A3 is the reason to expect readiness scoring to be generous to it:
   under readiness, oncogene delivery alone bought 70% of the GSE240057 displacement. **If k′ is high under
   readiness, the pre-registered follow-up is to repeat the ladder under activity scoring (A3's estimators),
   where that artefact vanished — and that follow-up is declared here, before the first number.**
3. Pancreatic acinar cells are a region member in the paper's own catalog (R1 and R2 both contain
   "pancreatic acinar cells"). **Gate G3 is therefore at serious risk on the acinar side**: the precursor may
   already be inside the region, which would make region-membership readouts uninformative for that origin —
   though the continuous statistics would still be valid, and the ductal side is unaffected.
