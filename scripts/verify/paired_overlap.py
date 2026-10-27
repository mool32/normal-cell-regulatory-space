#!/usr/bin/env python3
"""Are the within-patient cohorts outside the 25-cancer panel? Mostly not, and this says how much.

The paper called the gastric paired test "genuinely held out". Resolving the panel's own cohorts to
their sources shows that none of the three paired cohorts is cleanly outside the panel:

  colorectal  the panel's Colorectal cohort IS Pelka GSE178341, so every paired patient's malignant
              cells are in the panel by construction. This was never claimed otherwise.
  gastric     the panel's Census stomach cohort resolves to dataset 19053a82 with donors Pat01-B..
              Pat24-B -- the 23 GSE150290 tumour samples, i.e. the paired patients themselves. The
              stored cohort object keeps no donor column, so that step needs the Census snapshot and
              is recorded in results/cohort_provenance.csv rather than recomputed here.
  lung        the panel's NSCLC cohort carries nine `Kim_Lee_2020_*` donors. Mapped through the
              GEO sample/patient table, five of them are among the ten paired patients of E12.

What survives is the design, not the independence: both arms come from the same patient, the same
tissue and the same protocol, so lineage and batch cancel. That is what the test is for. It is not
an out-of-sample test of which normal cell the panel names, and the paper no longer says it is.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / 'results'
MAP = R / 'referee2_newdata/GSE131907_sample_patient.csv'      # !Sample_title -> patient id, from
NSCLC = ROOT / 'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad'   # GEO SOFT
E12 = R / 'referee2_newdata/E12_lung_negative.csv'
PELKA = R / 'cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad'
E36 = R / 'referee2_newdata/E36_cells.csv.gz'


def main():
    bad = 0

    # ---- lung: E12's ten paired patients against the panel's NSCLC donors
    m = pd.read_csv(MAP)
    s2p = dict(zip(m['sample'], m.patient))
    e12 = pd.read_csv(E12).patient.astype(str).tolist()
    paired = {s2p.get('LUNG_T' + p[-2:]) for p in e12} - {None}
    a = ad.read_h5ad(NSCLC, backed='r')
    o = a.obs[a.obs['cell_type'].astype(str) == 'malignant cell']
    panel = {d.replace('Kim_Lee_2020_', '') for d in o['donor_id'].astype(str) if 'Kim_Lee' in d}
    a.file.close()
    hit = sorted(paired & panel)
    print(f'lung        E12 paired patients: {len(paired)}  ({", ".join(sorted(paired))})')
    print(f'            panel NSCLC Kim_Lee_2020 donors: {len(panel)}  ({", ".join(sorted(panel))})')
    print(f'            IN BOTH: {len(hit)} of {len(paired)}  ({", ".join(hit) or "none"})')
    print(f'            the NSCLC cohort is the panel\'s smallest: {len(o):,} malignant cells')
    if len(hit) != 5:
        print('            !! expected 5; the registry quotes that number'); bad += 1

    # ---- colorectal: the panel's colorectal cohort is Pelka, so overlap is total by construction
    pe = ad.read_h5ad(PELKA, backed='r')
    col = [c for c in ('PatientTypeID', 'patient', 'PatientID', 'donor_id') if c in pe.obs.columns]
    pat = set(pe.obs[col[0]].astype(str).str.replace('.T$|.N$', '', regex=True)) if col else set()
    pe.file.close()
    e36 = set(pd.read_csv(E36).patient.astype(str))
    print(f'colorectal  paired patients {len(e36)}; the panel cohort is the same object '
          f'({col[0] if col else "no patient column"}, {len(pat)} ids), overlap by construction '
          f'{len(e36 & pat) if pat else "n/a"}')

    # ---- gastric: recorded, not recomputable offline
    cp = pd.read_csv(R / 'cohort_provenance.csv')
    g = cp[cp.cohort == 'Gastric']
    note = str(g['discrepancy against Table S1 / the text'].iloc[0]) if len(g) else ''
    print(f'gastric     {note[:150]}')
    if 'ARE the 23 GSE150290' not in note:
        print('            !! the provenance record no longer states the gastric identity'); bad += 1

    print(f'\nno paired cohort is cleanly outside the panel; {bad} check(s) failed')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
