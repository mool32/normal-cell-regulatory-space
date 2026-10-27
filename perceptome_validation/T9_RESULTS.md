# T9 — results: human cells transformed *in vitro* from a known precursor

Run 2026-10-19 against `PREREG_T9_invitro.md` (frozen the same day, before any expression value was read).
Scripts: `scripts/T9_freeze_attractor.py`, `scripts/T9_invitro_transformation.py`, `scripts/T9b_specificity.py`.
Outputs: `results/T9_attractor.json`, `T9_per_sample[_hpa].csv`, `T9_pairs[_hpa].csv`, `T9_cross_origin[_hpa].csv`, `T9b_specificity.json`.

**Verdict: SUPPORTS, with the effect carried by one series and by magnitude.** The pre-registered decision rule is met — P1 positive at *P* = 0.012 in GSE118206 and P2 positive in every gate-passing series — but four of the nine primary contrasts do not reach significance, two whole series fail the senescence gate, and the magnitude-removed convergence does not beat a random-panel null.

---

## 1. The attractor freeze reproduces the paper

`T9_freeze_attractor.py` rebuilt the frame from scratch and the reference cosines to the malignant attractor come out as the published geometry: megakaryocyte +0.899, erythrocyte +0.899, cytotrophoblast **+0.898** (T8 reports +0.8984 for the incumbent catalog), gastric progenitor +0.887. Gate passed; nothing downstream was refit.

## 2. P1 — two different epithelia, the same five hits (GSE118206)

Prostate basal cells and bronchial epithelium, transformed by the same PARCB cocktail, both precursors in the same series:

| | cos-distance between the two origins | |
|---|---|---|
| | **precursors** | **transformed** |
| raw | 0.285 | **0.004** |
| shape-only (magnitude removed) | 0.190 | **0.010** |

Δ raw **+0.281**, permutation *P* = **0.012**; Δ shape-only **+0.179**, *P* = **0.004** (10 000 label permutations within origin). Both hold under the alternative calibration (§6): +0.302 / *P* = 0.013 and +0.161 / *P* = 0.009.

Two human epithelia that start a third of a radian apart end up essentially on top of each other. This is the paper's Fig. 2 claim run as an experiment rather than inferred from a cross-sectional snapshot, and it names no normal cell.

**GSE6885 cannot test it.** HME and BPE, the two mammary precursors, start 0.024 apart — already indistinguishable in this frame — so there is no room to converge (Δ −0.0003, *P* = 0.61). Reported as uninformative, not as a negative.

## 3. P2 — movement toward the malignant attractor

Per-pair Δcos to the frozen attractor, pre-registered calibration (permutation *P*):

| series | contrast | Δ raw | *P* | Δ shape | gate |
|---|---|---|---|---|---|
| GSE118206 | prostate basal → PARCB | **+0.449** | 0.011 | −0.073 | |
| GSE118206 | NHBE → NHBE-PARCB | **+0.169** | 0.010 | +0.061 | |
| GSE6885 | HME → HMLER | **+0.207** | 0.009 | −0.020 | |
| GSE6885 | BPE → BPLER | +0.115 | 0.082 | +0.079 | |
| GSE3151 | HMEC → RAS (18 h) | **+0.190** | 0.0001 | −0.114 | |
| GSE3151 | HMEC → SRC / MYC / β-cat | +0.078 / +0.061 / +0.026 | 0.0001 / 0.0005 / 0.08 | | **OIS** |
| GSE3151 | HMEC → E2F3 | −0.044 | 1.00 | −0.053 | |
| GSE112369 | normal stomach → CDH1/TP53 KO organoid | +0.063 | 0.16 | +0.039 | |
| GSE57965 | adenoma → genome-engineered organoid | +0.055 | 0.17 | +0.019 | |
| GSE231776 | BJ fibroblast → +HRAS | +0.012 | 0.24 | +0.011 | n(prec)=1 |
| GSE22573 | keratinocyte LacZ → CDK4+HRAS | +0.029 | 0.11 | −0.099 | **OIS** |
| GSE172029 | melanocytic lineage → BRAF-V600E | −0.11 … +0.09 | ns | | **OIS 5/6** |

**Specificity controls behave.**

- **U0126** (MEK inhibition of the transformed keratinocytes) moves them **away**: −0.114.
- **Engineered control organoid** (CRISPR machinery, no driver): −0.003.
- **E2F3**, the one HMEC oncogene that does not transform: −0.044.
- **Real patient tumour organoids** in the same series move **+0.142** (*P* = 0.005) — more than twice the engineered ones.
- **Established prostate cancer lines** (LNCaP, VCaP, DU145, 22RV1, H660) move **+0.324** — less than the PARCB cells.

That last one is a caveat, not a confirmation: engineered cells move *further* than real cancer lines, so the displacement is not a pure malignancy signal.

**No dose-response.** The factor-dropout arms (4 of the 5 PARCB hits) move +0.449, i.e. as far as the full cocktail. The series cannot grade the effect.

## 4. The destination: region reproduced, cell inside it not

All eleven transformed PARCB samples leave their precursors' addresses (oesophageal basal, endometrial secretory) and land inside the paper's region, R2 occupancy 0.00 → **1.00**:

- under the **pre-registered** calibration, 10 of 11 land on **gastric progenitor cells**;
- under the **HPA-universe** calibration (§6), 11 of 11 land on **migrating cytotrophoblasts**.

The region is invariant; the cell named inside it flips with a normalisation choice. This is now the **fourth independent demonstration** of the same thing (T4 gene depth, T5 independent curation, T7 activity scoring, T9 here) and the first on an independent data type. It is direct support for demoting the anchor out of Results.

## 5. Specificity: is the module panel doing the work? (T9b)

200 expression-matched random gene panels of the same sizes, through the identical pipeline, on the same samples:

| | observed (44 modules) | random matched panels, median [95%] | *P* |
|---|---|---|---|
| Δ raw | **+0.281** | +0.023 [−0.001, +0.118] | **0.005** |
| Δ shape-only | +0.179 | +0.059 [+0.008, +0.217] | 0.060 |

And the generic view: 2 000 HVG of the series, correlation distance — precursors 0.172 → transformed 0.023, **Δ +0.148**.

Three honest conclusions:

1. The raw convergence is **specific to the module panel** — it beats 199 of 200 matched random panels.
2. The magnitude-removed convergence is **not** better than random panels (*P* = 0.06). Same verdict the paper already reaches for its own shape-only claim.
3. The convergence is **visible without the modules at all**, in a plain highly-variable-gene view. The module space sharpens it; it does not create it. Anyone who runs this series will see the two arms meet.

## 6. One calibration decision taken after seeing data, declared

The pre-registered rule scaled each sample to a total of 1e4 **over whatever genes the platform measures**. That denominator is 55 k genes for recount3 and 12–21 k for the arrays, which made absolute addresses platform-dependent: BJ fibroblasts scored cos 0.94 to the attractor and were assigned "gastric progenitor", HMEC and keratinocytes both landed on "ocular epithelial cells". Gate G4 already forbade using absolute addresses across series, so nothing in the decision rule depended on it — but the labels were plainly wrong.

A second calibration (`--hpa`) restricts every matrix to the genes HPA itself measures before scaling, making the denominator comparable to the reference. It fixes the artefact (BJ fibroblast → 0.47, outside the region) and changes no sign of any reported delta. **Both are reported throughout; the pre-registered one is primary.** The one substantive thing it changes is the *name* of the cell the PARCB samples land on — which is itself §4's result.

## 7. Gates

| gate | outcome |
|---|---|
| G1 coverage | 201–219 of 219 module-gene slots in every series (91.3–100%). None excluded. |
| G2 OIS | fires on GSE22573 (CDK4+HRAS), on 5 of 6 arms of GSE172029, and on 3 of 5 GSE3151 arms. GSE172029 delivers no non-senescent transformed arm and is excluded from the headline, as pre-specified. |
| G3 baseline ceiling | no gate-passing series had a precursor already inside R2 under the HPA calibration. Under the pre-registered one, GSE231776's precursor did — a normalisation artefact (§6), noted. |
| G4 platform | respected: every statistic is a within-series paired delta. |

## 8. What this licenses in the paper

One paragraph, and no more:

> In human cells transformed *in vitro* from a known precursor, the transformed state moves toward the convergence region measured against the precursor itself, so no cell-of-origin assumption is involved. Two different normal epithelia — prostate basal cells and bronchial epithelium — given the same five oncogenic hits end at a cosine distance of 0.004 from each other, against 0.285 for the two precursors (permutation *P* = 0.012; magnitude-removed 0.010 against 0.190, *P* = 0.004), and both land inside the region. The displacement is specific to the module panel against 200 expression-matched random panels in raw scores (*P* = 0.005) but not after magnitude removal (*P* = 0.06), and the convergence is also visible in a plain highly-variable-gene view, so the panel sharpens a real effect rather than creating it.

What it does **not** license:

- Any claim about which normal cell the region is. The identity flipped between two defensible normalisations of the same data.
- A general statement about spontaneous human cancer. The PARCB endpoint is small-cell neuroendocrine transdifferentiation driven by one fixed cocktail, so part of the convergence is built into the design. This was stated in the prereg before the run and must be stated in the paper.
- A dose-response. The 4-of-5-factor dropout arms move as far as the full cocktail.

## 9. Next, in order of value

1. **A second two-origin series.** GSE118206 carries the whole P1 result. Candidates with both precursors profiled and different origins: GSE148714 (cord-blood vs bone-marrow HSPC, same MLL fusion), GSE139822 + GSE118206 (bladder as a third PARCB epithelium, precursor borrowed), GSE240057 (the PARCB time course — does the path pass through the region or jump).
2. **An organoid series with a true wild-type arm.** GSE112369's normal stomach organoids are the only clean normal→engineered pair here and they give +0.063. Chen 2021 (conventional vs serrated colon) would add a second origin inside one organ.
3. **A non-oncogenic stress arm** (hypoxia, GSE109367) as a harder specificity control than U0126.

---

# Addendum A1 — the second two-origin series: **does not replicate**

Run 2026-10-19 against Addendum A1 of the prereg, which was written before any value of these two series
was read. Both calibrations give the same answer; the numbers below are the pre-registered one.

## A1.1 GSE17941 — two cell types, *different* hits, one platform, one lab

| contrast | cos-distance, start → end | Δ | perm *P* |
|---|---|---|---|
| **P1** MCF-10A precursor ↔ BJ-EL **→** Er-Src/TAM ↔ BJ-ELR | 0.285 → 0.159 | +0.126 | 0.0014 |
| **P1c, the specificity twin** — the two *non-transforming* arms (Er-Src+EtOH ↔ BJ-EH) | 0.285 → 0.169 | **+0.115** | 0.0037 |
| P1 shape-only | 0.259 → 0.168 | +0.091 | 0.0095 |
| P1c shape-only | 0.259 → 0.108 | **+0.151** | 0.025 |

**The control converges as much as the transformed arms, and more after magnitude removal.** The prereg
fixed the reading of this in advance: *"A convergence that also appears here is a culture/handling effect,
not transformation."* GSE17941's P1 is therefore void, not positive.

The arm-centroid distance matrix shows why, and it is unambiguous:

```
                   tr_MCF10A  ctrl_MCF10A  prec_MCF10A | ctrl_BJ_EH  prec_BJ  tr_BJ
tr_MCF10A              0.000        0.009        0.025 |     0.137    0.228  0.159
ctrl_MCF10A            0.009        0.000        0.010 |     0.169    0.254  0.189
prec_MCF10A            0.025        0.010        0.000 |     0.184    0.285  0.212
ctrl_BJ_EH             0.137        0.169        0.184 |     0.000    0.049  0.028
prec_BJ                0.228        0.254        0.285 |     0.049    0.000  0.018
tr_BJ                  0.159        0.189        0.212 |     0.028    0.018  0.000
```

Every MCF-10A arm stays within 0.025 of every other MCF-10A arm; every BJ arm within 0.049 of every other
BJ arm; the gap between the two blocks never falls below 0.137. **The two lineages never leave their own
territory.** Neither enters the region (R2 occupancy 0.00 throughout; cos to attractor 0.02–0.38). The MCF-10A
transformed arm in fact moves *away* from the attractor (−0.158), and both transformation arms fire the
senescence gate.

**GSE65885**, the RNA-seq version of the same two systems (n = 1 per arm, direction only), agrees: the BJ
arm moves +0.330 toward the attractor but its non-transforming EH control moves +0.191 of that, and MCF-10A
moves −0.017. Same non-specificity, different platform.

## A1.2 What this does to the first result

It does not void GSE118206 — P1c cannot be computed there, because that series contains no non-transforming
manipulation — but it removes the replication, and it sharpens two qualifications that must now be stated
whenever the first result is quoted.

**(i) In GSE118206 the convergence is mostly one origin travelling.** Arm-centroid distances:

```
               prec_prostate  prec_lung  tr_prostate  tr_lung  ref_line  ctrl_dropout
prec_prostate          0.000      0.285        0.537    0.479     0.452         0.463
prec_lung              0.285      0.000        0.127    0.116     0.108         0.101
tr_prostate            0.537      0.127        0.000    0.004     0.042         0.005
tr_lung                0.479      0.116        0.004    0.000     0.043         0.007
```

Prostate basal cells move 0.537; NHBE moves 0.116, because cultured bronchial epithelium already sits
0.116 from the destination (cos to attractor 0.735 before transformation, against 0.469 for prostate basal).
"Two epithelia converge" is accurate; "both travel" is not.

**(ii) The destination is shared with real cancer lines, and with incomplete cocktails alike.** The
transformed arms sit 0.042 from six established prostate-cancer and small-cell lines (LNCaP, VCaP, DU145,
22RV1, NCI-H660, MSKCC-EF1) that have nothing to do with PARCB — which argues the destination is not an
artefact of the cocktail. But the 4-of-5-factor dropout arms sit 0.005 away, so the coordinate cannot
distinguish a complete transformation from an incomplete one.

## A1.3 Verdict, revised

**Pre-registered P1 holds in one series and fails its specificity twin in the second.** The honest statement
is a boundary, not a confirmation:

> Two primary human epithelia given the same five hits end on top of each other and inside the region, together
> with independently derived cancer lines (GSE118206). A second pair of human cell types, transformed by
> different hits, does not reproduce this: they shrink toward each other, but their non-transformed controls
> shrink just as much, and neither lineage leaves its own territory or enters the region (GSE17941, GSE65885).

What separates the two is not settled by these data. Three candidate explanations, in the order I would test them:

1. **Cell type.** GSE118206 starts from primary tissue; GSE17941 starts from two immortalised lines already
   adapted to culture, which may have spent their room to move (the GSE6885 problem in another form).
2. **Completeness.** MCF-10A at 24–36 h of Src induction is early transformation; a later or fully
   tumorigenic arm might behave differently.
3. **The cocktail.** The pre-stated worry about GSE118206 — that one fixed five-factor NE programme
   guarantees part of the convergence — survives untested, since the series that was supposed to test it
   failed its own control.

Until a third series settles it, the paper should quote GSE118206 with the shared-cocktail caveat **and** the
non-replication, in the same breath. A positive first result whose only replication attempt failed is one
experiment, not a validated claim.

## A1.4 Next

The third series has to supply what both of these lack: **two different primary (not immortalised) human
cell types, a genuinely transformed arm, and a non-transforming manipulation of both**, so P1c exists.
Candidates, in order: a PARCB bladder arm (GSE139822) against GSE118206's prostate and lung — same cocktail,
third epithelium, but no within-series precursor; Chen 2021 conventional vs serrated colon organoids — two
origins inside one organ, both primary; GSE240057, the PARCB time course, which can at least say whether the
path passes through the region or jumps to it.

---

# Addendum A2 — the third series: **the construct delivers most of it, the address does not**

Run 2026-10-19 against Addendum A2, frozen before any value of these series was read. Thresholds on *f* were
fixed in advance. Both calibrations are reported; the pre-registered one is primary.

## A2.1 GSE240057 — PARCB time course, 10 donors, paired within donor

Mean cosine to the frozen attractor by stage (pre-registered calibration):

| stage | cos to A, raw | cos to A, shape-only | donors whose top-1 is in R2 |
|---|---|---|---|
| primary prostate basal, 0 wks | 0.208 | 0.720 | **0 / 10** |
| PARCB-transduced organoid, 0 wks (all five factors, no tumour) | 0.576 | 0.726 | 1 / 10 |
| early tumour, 4–5 wks | 0.512 | 0.724 | 6 / 10 |
| transition tumour | 0.687 | 0.660 | 12 / 15 |
| endpoint tumour, 8–18 wks | **0.734** | 0.673 | **10 / 10** |

Paired within donor:

- δ_endpoint raw **+0.526**, up in **10 of 10** donors, Wilcoxon *P* = 0.002.
- δ_cocktail raw **+0.368**, up in 10 of 10.
- **f = 0.699** (HPA calibration 0.625). The pre-registered thresholds were *f* ≥ 0.7 = cocktail, ≤ 0.3 = progression. The result lands one thousandth below the cocktail line: **AMBIGUOUS, and only by a hair.**
- Stage ordering: mean Spearman **ρ = +0.77**, positive in **10 of 10** donors.
- δ_endpoint **shape-only −0.047**, up in only 3 of 10, *P* = 0.13. δ_cocktail shape-only +0.006.

Two things follow, and they point in different directions depending on which metric is used.

**On the continuous cosine, the oncogenes do most of the work before any tumour exists.** Sixty-three to
seventy per cent of the total displacement is already present at week 0, in cells that carry all five factors
and are not yet tumours. The stated worry about GSE118206 — that the coordinate reads the construct rather
than the transformed state — is therefore largely **confirmed** for the continuous metric.

**On the paper's own metric it is the opposite.** Region membership — the nearest-normal assignment that
Fig. 3 and Fig. 5 actually use — moves late and cleanly: 0 of 10 donors at baseline, 1 of 10 after
transduction, 6 of 10 early, 12 of 15 at transition, **10 of 10 at endpoint**. The construct buys one donor;
the tumour buys the other nine. Under the HPA calibration the ladder is noisier (0 → 3 → 2 → 8 → 8 of 10)
but ends in the same place.

So the two currencies dissociate, and the paper happens to use the one that behaves: *how far along the
attractor direction* a cell sits is mostly oncogene expression; *which normal cell it is nearest* changes
with transformation.

**And the shape claim fails here.** Endpoint tumours are not closer to the attractor than their own donor's
basal cells once magnitude is removed (−0.047, 3 of 10). Prostate basal cells already sit at shape-cosine
0.720, and the tumours drift marginally away.

## A2.2 GSE148714 — the same fusion, two sources, only one transforms

δ to the pooled CD34+ HSPC control:

| arm | δ raw | δ shape-only | top-1 |
|---|---|---|---|
| **CB + MLL-AF4** (transforms) | **+0.079** | **+0.137** | megakaryocytes |
| **BM + MLL-AF4** (does *not* transform) | +0.049 | **+0.017** | late spermatids |
| CB + MLL-AF9 | +0.068 | +0.117 | gastric progenitor cells |
| BM + MLL-AF9 | +0.061 | +0.086 | megakaryocytes |

The sign is right: the source that transforms moves further, by 1.6× raw and **8× after magnitude removal**
(+0.137 against +0.017). That is the sharpest specificity signal anywhere in T9, and it is the one thing
GSE118206 could not supply and GSE17941 failed to supply.

**Two caveats that cost it most of its weight.**

1. **My pre-committed expectation was wrong.** I wrote that "haematopoietic precursors are far from the
   convergence region, so 3b has room to move in either direction". They are not: the CD34+ HSPC control
   already sits at cosine **0.900**, and its top-1 is already a region member (oesophageal apical). Gate G3
   fires — this series measures depth inside the region, not entry into it.
2. **The destination is the cell the paper's own control rejects.** CB-MLL-AF4's top-1 is the
   **megakaryocyte**, which Note S-Q shows falls with transcriptome complexity exactly as the quarantined
   erythrocyte does. A haematopoietic progenitor moving toward the megakaryocyte is close to uninformative,
   because the attractor direction itself sits at cosine 0.899 from the megakaryocyte.

## A2.3 Where T9 now stands, across all three series

| series | design | result |
|---|---|---|
| GSE118206 | two **primary** epithelia, same five hits, no non-transforming control | **positive**: converge to 0.004 from 0.285, raw *P* = 0.012 and shape-only *P* = 0.004; both enter the region; land 0.042 from independent cancer lines |
| GSE17941 / GSE65885 | two **immortalised** lines, **different** hits, with a non-transforming twin | **void**: the twin converges as much (raw) or more (shape); neither lineage leaves its own territory or enters the region |
| GSE240057 | one origin, **10 donors**, paired trajectory with a transduced-but-not-yet-tumour arm | **split**: raw displacement +0.526, 10/10, monotone 10/10 — but 63–70% of it is delivered by the construct; region membership moves late and completely (0/10 → 10/10); shape-only movement absent |
| GSE148714 | same fusion, **two sources, only one transforms** | **right sign, little weight**: 1.6× raw and 8× shape-only in favour of the transforming source, but the precursor already sits at 0.90 inside the region and the destination is the complexity-rejected megakaryocyte |

**Revised verdict: BOUNDARY — one claim survives, two do not.**

Survives, and is new:

> Across ten donors, primary human prostate basal cells progress from no donor inside the convergence region
> to every donor inside it by the endpoint tumour, monotonically in all ten (mean Spearman ρ = +0.77), with a
> paired increase in attractor proximity in 10 of 10 (*P* = 0.002). The comparison is to each donor's own
> precursor, so no cell-of-origin assumption enters.

Does not survive:

- **That the coordinate reads transformation rather than oncogene expression.** On the continuous metric
  it mostly reads the construct: 63–70% of the displacement is present at week 0 in cells that are not
  tumours. Only the discrete region-membership metric separates the two (1 of 10 donors at week 0 against
  10 of 10 at endpoint), and that dissociation must be reported, not chosen.
- **That transformation moves cells toward the region in module *shape*.** One series says yes
  (GSE118206, though not better than matched random panels at *P* = 0.06), one says no (GSE240057, −0.047,
  3 of 10), one says the apparent effect is non-specific (GSE17941).

## A2.4 What I would do next, and what I would not

**Would not:** add a fourth series of the same kind. Three have now been run and the limiting factor is not
sample size, it is that every available human transformation system either starts from an immortalised line,
delivers its oncogenes as one fixed cocktail, or starts already inside the region.

**Would**, in order:

1. **Report the dissociation as a result.** The construct-versus-tumour split in GSE240057 is the most
   informative single number T9 produced and it belongs in the SI whatever else is decided.
2. **Re-run GSE240057 under activity scoring (T7's estimators).** If the week-0 jump is driven by the
   transcript abundance of MYC and AKT1 themselves, footprint-based scoring should shrink it. That is a
   direct test, on data already on disk, and it also answers Gudkov's protein objection in the same move.
3. **Find a system with staged, inducible delivery** so that factor expression and transformed state are
   separable in time by design, rather than separated post hoc.

---

# Addendum A3 — activity scoring: **the week-0 jump is an abundance artefact, and the shape claim flips positive**

Run 2026-10-19 against Addendum A3, frozen before any activity score of this series was computed.
Thresholds on *f* were fixed in advance. Scripts: `scripts/T9d_activity.py`.
Module definitions are T7's, unchanged; the reference and the 25-cancer attractor were rebuilt in activity
space from T7's cached matrices and frozen to `T9d_attractor_{viper,ulm}.json` before GSE240057 was scored.

## A3.1 The headline numbers, GSE240057, paired within 10 donors

| scoring | δ endpoint, raw | δ endpoint, **shape-only** | δ cocktail (week 0) | **f** | stage ρ |
|---|---|---|---|---|---|
| **readiness** (the paper's: mean component abundance) | +0.526, 10/10, *P* = 0.002 | **−0.047, 3/10, *P* = 0.13** | +0.368, 10/10 | **0.699** | +0.77, 10/10 |
| readiness minus the delivered genes (control C3) | +0.503, 10/10, *P* = 0.002 | −0.019, 4/10 | +0.354, 10/10 | 0.705 | +0.77, 10/10 |
| **activity, VIPER** (pre-registered primary) | +0.135, 10/10, *P* = 0.002 | **+0.148, 10/10, *P* = 0.002** | **−0.038, 0/10** | **−0.279** | +0.72, 10/10 |
| **activity, ULM** | +0.495, 10/10, *P* = 0.002 | **+0.324, 10/10, *P* = 0.002** | +0.186, 10/10 | **0.375** | +0.78, 10/10 |

Two things changed, and both were predicted in advance by someone other than me — by T7's own directional
result, on data T7 never touched.

**1. The construct's contribution collapses.** Under readiness scoring 70% of the displacement was present
at week 0, in transduced cells that were not yet tumours. Under VIPER the transduced arm moves *backwards*
(−0.038, 0 of 10 donors), giving *f* = −0.279; under ULM it accounts for 37%. The pre-registered threshold
for "the week-0 jump was component-transcript abundance" was *f* ≤ 0.3. **VIPER meets it; ULM is partial and
in the same direction.** The objection that the coordinate reads the cocktail rather than the transformed
state is therefore an objection to *readiness* scoring specifically, not to the coordinate.

**2. The shape claim flips from failure to success.** Readiness scoring put the endpoint tumours *further*
from the attractor than their own donor's basal cells once magnitude was removed (−0.047, 3 of 10,
*P* = 0.13). Under activity scoring the same comparison is **+0.148 (VIPER) and +0.324 (ULM), 10 of 10
donors, *P* = 0.002 under both**. This is T7's pre-registered directional prediction — that the shape-only
share of the signal rises when modules are read through regulatory output — reproduced on an independent
data type, in a paired within-donor design, with no cell-of-origin assumption anywhere.

**Control C3 behaves as pre-committed.** PARCB delivers shRNAs against TP53 and RB1 and overexpresses
myrAKT1, MYC and BCL2; three of those sit in scored modules (RB1 → Cell Cycle, AKT1 → PI3K/PTEN,
TP53 → p53). Removing all three leaves *f* at 0.705 against 0.699. The week-0 jump is not those three
transcripts; it is the broader abundance profile of transduced, rapidly growing cells — which is exactly
what activity scoring is designed to look past.

## A3.2 What the two estimators do not agree on

| | VIPER | ULM |
|---|---|---|
| nearest reference cells to the activity-space attractor | megakaryocytes 0.923, oocytes 0.907, late spermatids 0.802 | oocytes 0.955, megakaryocytes 0.941, late spermatids 0.904, **cytotrophoblasts 0.900** |
| donors whose endpoint top-1 is in R2 | **10 / 10** | **0 / 10** |
| top-1 at endpoint | pancreatic acinar cells | tuft cells |

So the *continuous* statistic agrees across estimators in sign, size-order and donor count, while the
*discrete* nearest-normal assignment does not agree at all. And the metrics trade places between scorings:

| | continuous cosine | discrete region membership |
|---|---|---|
| **readiness** | contaminated (*f* = 0.699) | clean (1/10 at week 0 → 10/10 at endpoint) |
| **activity** | clean (*f* = −0.28 / 0.38) | contaminated or empty (VIPER 8/10 already at week 0; ULM 0/10 ever) |

In each scoring exactly one of the two metrics separates construct from tumour. That symmetry is
uncomfortable and I am not going to pick the convenient half of it. The defensible reading is narrower than
either: **the continuous, magnitude-removed, activity-level distance from a donor's own basal cells to that
donor's endpoint tumour increases in 10 of 10 donors under both estimators, and the transduced-not-yet-tumour
arm contributes little or none of it.** Everything about *which* normal cell the destination is remains
estimator-dependent, as T7 already established and as A1 and A2 each found again.

## A3.3 What this adds to the paper, and to the two open questions

**To the cocktail objection (Gudkov, point 2):** answered on its own terms. The week-0 displacement that
made A2 ambiguous is a readiness artefact; at the activity level the five factors buy nothing and the tumour
buys everything.

**To the protein objection (Gudkov, point 1, "вторая производная"):** this is now a positive argument rather
than a concession. Scoring modules by regulatory output instead of component abundance does not weaken the
result — it is the only scoring under which the magnitude-removed claim holds at all, in 10 of 10 donors,
on data assembled by someone else for another purpose. T7 found the same on the 25 cancers and on Pelka.
The paper should report the activity-level version alongside the readiness one, which also discharges the
pre-registration's obligation to report T7.

**The cost, stated plainly:** the activity-space attractor's nearest neighbours are the oocyte, the
megakaryocyte and the spermatids — the complexity-suspect cells of Note S-Q, and under VIPER the
cytotrophoblast is not in the top eight at all. Every strengthening of the *state* claim at the activity
level is a further weakening of the *address* claim. That is the same trade A1 and A2 found, now for the
third time, and it is the strongest internal argument for moving the anchor out of Results.

---

# Addendum A4 — GSE222990: a decomposition, not a refutation

**This section replaces an earlier version of A4 that reported a REFUTES verdict. A five-lens adversarial
critique (2026-10-19) found that verdict rested on a statistic one sample decides, that several quoted
figures came from the secondary calibration while the prose claimed the primary, and that my reading of the
substitution test was backwards. Everything below is the corrected account; §A4.4 lists every correction,
including four that reach back into A1–A3.** All figures are the pre-registered primary calibration unless
marked. Script `scripts/T9e_pancreas.py`; artefacts `T9e_summary.json`, `T9e_arm_matrix_{all,hpa}.csv`,
`T9e_GSE222990_per_sample_{all,hpa}.csv`.

## A4.1 The acinar ladder is statistically empty — withdrawn as evidence

With culture held constant (baseline = the same acinar cells in organoid culture, no drivers, same platform):

| | δ′ raw | one-sided *P* | δ′ shape-only |
|---|---|---|---|
| acinar + four drivers (KPTS) | −0.0592 | 0.81 | −0.027 |
| acinar + KRAS alone | −0.0593 | 0.79 | −0.019 |

The pre-registered REFUTES trigger was "δ′(KPTS) ≤ 0". It is met. **It should not be acted on.**

`acinar_KPTS_organoid` has **n = 3**, with cos_A = 0.816, 0.212, 0.285 — SD 0.330, SEM **0.190**. Against the
six-sample organoid arm (SEM 0.041) the difference carries SE ≈ 0.195, *t* = −0.30. **Leave-one-out flips the
sign in 2 of the 3 drops: −0.248, +0.054, +0.017.** Whether this series "refutes" depends on whether one
sample is in the arm. The quoted *P* = 0.81 is also a one-sided test pointing away from the observed effect;
the informative figure is *P*(less) = 0.27.

*k*′ and *g* inherit the same problem and worse: with δ′ negative they are ratios of two negative
displacements, so the pre-registered reading ("≥ 0.7 means that arm explains it") does not apply at all —
the pathology I correctly flagged for ULM's *k*′ = −1.21 and failed to flag for *k*′ = 1.002. Across scorings
*k*′_shape spans **0.087 (VIPER) to 0.798 (ULM)**, i.e. every one of the three pre-registered bins, and
*g* = **0.340 under VIPER**, below the cutoff — a number the earlier write-up omitted while quoting 0.891
and 1.403.

On the ceiling: under the primary calibration the fresh-acinar arm mean is cos **0.884** with top-1
`oesophageal basal cells` (3 of 5 inside R2). The earlier "cos 0.905, above the cytotrophoblast's 0.898, and
top-1 is literally cytotrophoblasts" mixed a centroid cosine with an arm mean and the secondary calibration
with the primary. Gate G3 is still a real concern for this origin, but not in those terms.

## A4.2 What the substitution actually shows — my earlier reading was backwards

P1 is positive: Δ = **+0.477**, permutation *P* = **0.0025**; shape-only **+0.093**, *P* = **0.0186**.
Substituting each non-transforming acinar arm for the transformed one, against the same ductal arm:

| acinar arm paired against ductal_KPTS | distance | Δ from the 0.561 baseline |
|---|---|---|
| + four drivers (KPTS) | **0.084** | +0.477 |
| + KRAS alone | 0.131 | +0.430 |
| culture only, no drivers | 0.137 | +0.424 |

I previously wrote that this means the no-driver arm "reproduces 89% of the convergence", and read that as
non-specificity. **That reading is wrong.** The 89% is a ratio against a very large baseline (0.561); what it
hides is that the drivers close a further **0.053** of distance, and that this residue is *consistent*:

- it survives all three leave-one-out drops of the 3-sample KPTS arm (remaining gaps **+0.029, +0.037,
  +0.043** in favour of KPTS);
- **0 of the 20** three-sample subsets of the six no-driver organoid samples comes as close to the ductal arm
  as the full KPTS arm does (subset range 0.115–0.175 against 0.084).

And the scale argument in the earlier version was wrong too. I compared the 0.063 driver step to a
*within-arm pairwise* spread of up to 0.069. The relevant null is the centroid-level separation,
≈ spread·½(1/n_a + 1/n_b) ≈ **0.017**, so 0.063 is three to four times it. **The drivers do move the cells.**
The culture step is simply 4.5× larger: 0.285 fresh→organoid against 0.063 organoid→KPTS.

## A4.3 Verdict

**Not a refutation. A decomposition, from an underpowered series.**

- The measured cross-origin convergence in GSE222990 is **dominated by culture**: the acinar arm travels 0.285
  on entering organoid culture and 0.063 more on receiving four drivers, and the ductal arm it meets is alone
  on a different platform (gate G6, declared before the run).
- The drivers nonetheless contribute a **small, consistent, leave-one-out-stable** residue.
- The attractor-directed statistic on the acinar side (δ′) is **uninterpretable at n = 3** and is withdrawn.

The pre-registered decision rule gives neither SUPPORTS nor a usable REFUTES here, because both of its
conjuncts run through δ′. The registry row is amended accordingly.

## A4.4 Corrections the critique forced, including into A1–A3

**1. A code bug that changes a pre-registered verdict.** `T9c_trajectory.py` stripped `.N` suffixes from the
GSE240057 gene index as if they were R `make.unique` duplicates. Of 24,878 dotted names, **only 1,869 are**;
the other ~23,000 are genuine clone-based GENCODE symbols (`AC000003.1`, `RP11-432M8.1`) with no bare base,
and stripping them dropped rows and perturbed each sample's normalisation denominator by 1–2%,
non-uniformly. Fixed (strip only when the bare base exists in the same file) and everything re-run.
**Consequence: A2's *f* = 0.699 becomes *f* = 0.708 — the pre-registered COCKTAIL verdict, not "AMBIGUOUS,
and only by a hair".**

**2. The surviving positive is weaker than I stated.** Re-run after the fix, GSE240057 under activity scoring:

| | δ_endpoint raw | δ_endpoint shape-only | *f* |
|---|---|---|---|
| VIPER | +0.052, **9/10**, *P* = 0.0039 | **+0.028, 8/10, *P* = 0.0195** | −0.603 |
| ULM | +0.401, 10/10, *P* = 0.0020 | +0.214, 10/10, *P* = 0.0020 | 0.324 |

**"10 of 10 donors, *P* = 0.002, raw and magnitude-removed, both estimators" is withdrawn.** It is 10/10
under ULM and 8–9/10 under VIPER. The *direction* is positive under both; the unanimity is not. ULM's
*f*_shape is **0.644**, near the 0.7 line, and was not reported before.

**3. Ratio pathology, applied symmetrically.** A3's *f*_act is now −0.603, a ratio whose numerator changed
sign — exactly the defect I called "uninterpretable" for ULM's *k*′ = −1.21 while treating A3's −0.279 as a
clean threshold pass. Both are now reported as sign-flipped and uninterpretable as ratios; what they do say
is that the transduced arm moves *away* from the attractor under VIPER (1 of 10 donors positive).

**4. Every "*P* = 0.002" is one bit.** `wilcoxon()` was called without `alternative`, so scipy's two-sided
exact test, whose minimum at *n* = 10 is 2/2¹⁰ = **0.001953125**. Every reported value equals that floor
exactly, so the statistic carries only "10 of 10 same sign" — raw, shape-only, VIPER and ULM return the
identical number on the same ten donors. Restated honestly as a sign test: 10 of 10, one-sided binomial
*P* = 0.00098. And a T9-sized family (57 stored *P* values) cannot be survived by Bonferroni at that floor;
the surviving claim is FDR-controlled, not family-wise.

**5. A permutation test that halved its own *P* values.** The per-pair test compared `>= obs` after a float
permutation, losing the tied draws to summation order. Exact enumeration for GSE118206's lung pair gives
**1/56 = 0.0179**; the run reported **0.0100**, below the combinatorial floor. Fixed with a tolerance.
GSE231776's 1-vs-3 pair reported 0.24 against a floor of 0.25 and should never have had a *P* at all.

**6. A pre-registered trigger I left unfired.** A3's control C3 returned *f* = 0.705, across the
pre-registered cocktail line. I reported it only as "behaves as pre-committed".

**7. The main verdict restated its own rule.** Prereg §5 requires "P2 > 0 in at least 5 of the 7
precursor-paired series that pass G1–G3"; the write-up said "P2 positive in every gate-passing series". Once
G2 removes GSE22573 and GSE172029 and G3 removes GSE231776, four series remain and the conjunct is
**unsatisfiable**. The T9 row is amended from SUPPORTS to NOT MET.

**8. An inverted pre-commitment.** Prereg §6: *"A negative result is reported in the paper regardless of what
the rest of the analysis shows, with the series named."* My recommendation deleted the cross-origin paragraph
and kept the positive. Corrected in §A4.5.

**9. Overwritten artefacts.** `T9_invitro_transformation.py` wrote fixed filenames and accepts a subset
argument, so the A1 run overwrote the numbers behind §2–§4. Output names are now per-run and the full set has
been re-run. The arm-centroid matrices, within-arm spreads, substitutions and leave-one-out checks — the
load-bearing statistics of A1.2 and A4.2 — were computed in an ad-hoc shell heredoc and are now emitted by
`T9e_pancreas.py` into `T9e_arm_matrix_*.csv` and `T9e_summary.json`.

**10. The surviving claim lives in a failed regime.** T7's own pre-registered verdict is FAIL, and the one
result proposed for the paper exists only under T7's scoring. That has to be said in the same sentence as the
result.

## A4.5 What may be written in the paper

Both halves, because the pre-registration requires the negative to be reported:

> In human cells transformed *in vitro* from a known precursor, comparison is to the precursor itself, so no
> cell-of-origin assumption enters. Across five series the apparent cross-origin convergence is dominated by
> culture rather than by transformation: in primary human pancreatic organoids the acinar arm travels 0.285
> in cosine distance on entering culture and 0.063 more on receiving four PDAC drivers, and two immortalised
> lines converge as much without being transformed. What survives is within-donor and requires scoring
> modules by regulatory output: across ten donors the endpoint tumour sits closer to the malignant attractor
> than that donor's own primary basal cells in 10 of 10 donors under one activity estimator and 8–9 of 10
> under the other (sign-test *P* = 0.001–0.02), with the transduced-but-not-yet-tumour arm contributing none
> of it, whereas under component-abundance scoring 71% of the displacement arrives with the construct before
> any tumour exists.

Numbers and sources: 0.285 and 0.063 → `T9e_arm_matrix_all.csv`; "two immortalised lines" → `T9_cross_origin.csv`
(GSE17941 P1 +0.126 against its non-transforming twin +0.115); 10/10 and 8–9/10 → `T9d_summary_{ulm,viper}.json`;
71% → `T9c_summary.json` (*f* = 0.708).

Not to be written: GSE118206's P1 as evidence that transformation causes convergence (no non-transforming
control exists in that design); any destination, top-1 or region-occupancy figure from GSE222990; the ductal
+0.300, which is alone on GPL21290; and any "*P* = 0.002 in 10 of 10 under both estimators".
