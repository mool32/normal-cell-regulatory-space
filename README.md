# Human malignancies show increased occupancy of a shared region of normal-cell regulatory space

Code and derived data for the paper. 360 files, 93 MB.

## Regenerate everything

```
python3 scripts/make_figures.py
```

About two minutes. It rebuilds every panel, composes the seven main figures and the fifteen
supplementary ones, writes the main text and the Supporting Information as separate documents
(`manuscript/paper3_v12_MAIN.pdf`/`.docx` and `paper3_v12_SI.pdf`/`.docx`), and finishes by running
the verification harnesses. Run it from this directory; the scripts resolve their inputs relative to
the repository root. Two consecutive runs produce byte-identical figures.

## What is here, and what is not

Everything under `results/`, `tool/` and `perceptome_validation/` is a *derived* table that the
figure code reads — the trace that built this bundle recorded each one being opened. The raw
single-cell matrices are not redistributed; they are public, and the accessions are listed in Data
availability. Four stages open one locally and name it before skipping:

  - `build_v7_clouds.py` — Fig. S4-I per-cancer clouds
      - `results/cancer_v6_4_extras/B-ALL_malignant.h5ad` — GSE132509 (B-ALL), Fig. S4-I only
      - `results/cancer_v6_4_extras/Neuroblastoma_malignant.h5ad` — GSE137804 (neuroblastoma), Fig. S4-I only
  - `verify/census_cohorts.py` — the Census cohorts resolved: selection rule and double counting
      - `results/cancer_scaled/cache/Prostate_malignant.h5ad`
  - `verify/paired_overlap.py` — how much of each paired cohort is already in the panel
      - `paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad`
      - `results/cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad`
  - `verify/check_E4_lazy.py` — E4's lazy cohort loader equals its eager one
      - `results/cancer_v6_6_extras/vanGalen_AML_malignant.h5ad`

Everything else builds without them. The three checks among those four establish facts that ship
here as finished tables (`results/census_cohorts.csv`, `results/cohort_provenance.csv`); a check
that cannot run says what it needs rather than passing silently.

## Verification

- `scripts/verify/check_numbers.py` recomputes every load-bearing number in the text from one named
  artifact each, and fails loudly when the artifact and the sentence disagree.
- `scripts/verify/cohort_provenance.py` resolves each of the 25 cohorts to its accession.
- `MANIFEST.txt` gives the sha256 of every file in this bundle.

A non-zero exit from `check_numbers.py` is a report about unresolved claims, not a broken build; the
runner labels those stages GATE.
