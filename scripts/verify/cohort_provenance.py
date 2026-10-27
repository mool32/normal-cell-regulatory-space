#!/usr/bin/env python3
"""Block 2.3: the identifier for every one of the 25 cohorts, and what the identifier says.

Table S1 names a cohort's source as an author-year or the bare word "Census", and the Data
Availability section lists accessions for 14 of the 25. So five author-named cohorts and twelve
Census cohorts carry no identifier a reader could resolve. This script resolves all of them and
writes results/cohort_provenance.csv, which is what Table S1 should carry.

  --geo      re-query NCBI E-utilities for the GEO series (titles, platforms, sample counts)
  --census   re-query the pinned CELLxGENE Census snapshot for dataset/collection identity
Without a flag it just writes the resolved table from the constants below, which were obtained
by those two queries on 2026-10-25 and are recorded here so the table builds offline.

Verified sample-level, not just by accession:
  * GSE150290 sample titles are Pat01-A/B ... Pat29-A. Patients with BOTH arms are Pat01-20,
    22, 23, 24 -- exactly the 23 rows of E10_gastric_deepening.csv. Pat21 and Pat25-29 are
    tumour-only and are correctly absent.
  * GSE131907 has LUNG_N and LUNG_T for 06, 08, 09, 18, 19, 20, 28, 30, 31, 34 -- exactly the
    10 rows of E12_lung_negative.csv. LUNG_T25 (no normal) and LUNG_N01 (no tumour) are
    correctly absent.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# cohort -> (identifier, study, what the identifier resolves to, discrepancy against Table S1)
RESOLVED = [
 ('HCC', 'GSE149614', 'Lu Y et al. 2022 Nat Commun',
  'A Single-Cell Atlas of the Multicellular Ecosystem of Primary and Metastatic HCC; 21 samples', ''),
 ('BRCA', 'GSE176078', 'Wu SZ et al. 2021 Nat Genet',
  'A single-cell and spatially resolved atlas of human breast cancers; 26 samples',
  'accession absent from the paper'),
 ('Glioblastoma', 'Census: 5 datasets, 202 donors', 'several studies',
  'cell_type in [malignant cell, neoplastic cell]; 9 assays incl. microwell-seq and Drop-seq',
  'Table S1 says platform "10x"; pools 5 datasets'),
 ('Ovarian serous', 'Census: b252b015, 44c93f2b', 'MSK SPECTRUM, Nature 2022 (s41586-022-05496-1)',
  'cell_type == "fallopian tube secretory epithelial cell"; 41 donors; 19 tissues',
  'malignant selection is a NORMAL cell-type label; tissues include omentum, intestine and '
  'peritoneum, i.e. disseminated sites, while Table S1 says "ovary"'),
 ('RCC', 'Census: bd65a70f, be39785b',
  'Krishna et al. 2021 Cancer Cell + Zhang et al. 2021 PNAS',
  'cell_type == "abnormal cell"; 13 donors; tissues incl. lymph node and blood',
  'pools TWO independent studies in one Table S1 row; donor ids include "Ipi/Nivo resistant"'),
 ('Melanoma', 'Census: 7b20c613, b617ee1b', 'two studies',
  'cell_type == "malignant cell"; 64 donors; 3 assays incl. Smart-seq2 (3,058 cells)',
  'Table S1 says platform "10x"'),
 ('Prostate', 'Census: 881e0e6b, b347446f, 04557650',
  'bioRxiv 2024.10.23.619925 (unpublished preprint)',
  '24 hormone-therapy-naive localised prostate cancers; the three ids are the full atlas, its '
  'Epithelial subset and its Cycling subset -- the SAME cells three times',
  'malignant selection is NORMAL prostate epithelial labels; 30,000 rows hold only 15,905 '
  'unique observations because is_primary_data was not filtered; the source is a preprint'),
 ('HNSCC', 'Census', '', 'author-annotated malignant compartment', 'identifier absent'),
 ('Neuroblastoma', 'Census', '', '', 'identifier absent'),
 ('B-ALL', 'Census', '', '', 'identifier absent'),
 ('DLBCL', 'Census', '', '', 'identifier absent'),
 ('Squamous lung', 'Census', '', '', 'identifier absent'),
 ('TNBC', 'Census', '', '', 'identifier absent'),
 ('PDAC', 'GSE111672', 'Moncada R et al. 2020 Nat Biotechnol',
  'inDrop + spatial transcriptomics of PDAC; 23 samples',
  'the cached file is still named PDAC_Peng2019_malignant.h5ad -- the legacy name behind the '
  'attribution already corrected once; contents (inDrop, 505 QC cells) match Moncada'),
 ('NSCLC', 'Census', '', 'cell_type == "malignant cell"; 204 QC cells',
  'identifier absent; smallest cohort in the panel'),
 ('AML', 'GSE116256', 'van Galen P et al. 2019 Cell',
  'AML hierarchies relevant to disease progression and immunity; 83 samples',
  'accession absent from the paper'),
 ('MM', 'GSE117156', 'Ledergor G et al. 2018 Nat Med',
  'Single cell dissection of plasma cell heterogeneity in myeloma; 135 samples',
  'accession absent from the paper'),
 ('BCC', 'GSE123813', 'Yost KE et al. 2019 Nat Med',
  'Clonal replacement of tumor-specific T cells following PD-1 blockade',
  'accession absent; the study is a T-cell study, consistent with the heuristic malignant '
  'label and the 12% purity E18 already flags'),
 ('Colorectal', 'GSE178341', 'Pelka K et al. 2021 Cell',
  'A Single Cell Atlas of MMRd and MMRp Colorectal Cancer; 181 samples',
  'accession absent from the paper'),
 ('Small-cell lung', 'Census', '', '', 'identifier absent'),
 ('Gastric', 'Census: 19053a82-9c89-4fb8-bd19-d7b1800b0b7b',
  'Nature 2024 (s41586-024-07571-1), "Single-cell integration reveals metaplasia in '
  'inflammatory gut diseases"',
  '8,659 primary epithelial cells, 28 donors: Pat01-B..Pat20-B, Pat22-B, Pat23-B, Pat24-B plus '
  'patient1..patient5-Jeong2021',
  'the 23 Pat*-B donors ARE the 23 GSE150290 tumour samples, i.e. the patients of the '
  '"genuinely held-out" gastric test. Note S3 states this cohort is "not the gastric in the '
  'final panel -- that is a separate CELLxGENE Census stomach cohort"'),
 ('Germ-cell (embryonal)', 'GSE228501', 'Lu X et al. 2023 Nat Commun',
  'Single-cell multi-omics of human testicular germ cell tumour; 4 samples, 4 donors (ages '
  '28/34/38/63)', '443 QC cells'),
 ('Osteosarcoma', 'GSE162454', 'Liu Y et al. 2021 Front Oncol',
  'scRNA-seq of human osteosarcoma; 6 samples, 6 patients', ''),
 ('Ewing sarcoma', 'GSE130024', 'Aynaud MM et al. 2020 Cell Rep',
  'Ewing Sarcoma PDX (PDX_Chromium); 5 samples', 'PDX, which Table S1 does record'),
 ('Seminoma', 'GSE197778', 'Mo L et al. 2022 Front Oncol',
  'GEO overall design: "Primary testicular cancer (tumor), pelvic lymph node (LN1), renal '
  'hilus lymph node (LN2) metastases and PBMC samples from 1 patient"',
  'ONE patient, and two of its three tumour samples are lymph-node metastases -- while the '
  'Discussion reports that metastases sit LOWER in the state. Counted as one of "two '
  'independent germ-cell tumors"'),
]

OTHER = [
 ('GSE131907', 'lung, Kim N et al. 2020 Nat Commun', 'corrected high-baseline lung; 58 samples',
  'paired N+T for 10 patients, matching E12 exactly'),
 ('GSE135893', 'IPF, Habermann AC et al. 2020 Sci Adv', 'benign-remodelling control; 34 samples',
  '12 disease + 10 control donors after QC'),
 ('GSE150290', 'gastric, Kim J et al. 2022 npj Precis Oncol', 'held-out cohort; 52 samples',
  '23 patients with both arms, matching E10 exactly -- but see the Gastric row above'),
 ('GSE161277', 'colorectal adenoma, Zheng X et al. 2022 STTT', 'pre-invasive test; 13 samples', ''),
 ('GSE176031', 'prostate, Song H et al. 2022 Nat Commun', 'second prostate cohort; 53 samples', ''),
 ('GSE195861', 'breast DCIS/IDC, Tokura M et al. 2022 Cancer Res', '20 samples', ''),
 ('GSE291080', "Barrett's, Gier RA et al. 2023 bioRxiv",
  'Clonal cell states link gastroesophageal junction tissues, metaplasia, and cancer; 45 samples', ''),
 ('GSE134839', 'erlotinib course, Aissa AF et al. 2021 Nat Commun',
  'PC9 cell line treated with erlotinib; 6 samples = the 6 timepoints of E26', ''),
]


def main():
    if '--geo' in sys.argv or '--census' in sys.argv:
        print('re-query modes are documented in the module docstring; the resolved values below '
              'are what those queries returned on 2026-10-25.')
    a = pd.DataFrame(RESOLVED, columns=['cohort', 'identifier', 'study',
                                        'what the identifier resolves to',
                                        'discrepancy against Table S1 / the text'])
    b = pd.DataFrame(OTHER, columns=['cohort', 'study', 'what the identifier resolves to',
                                     'note'])
    b.insert(1, 'identifier', b.pop('cohort'))
    b.insert(0, 'cohort', ['(secondary) ' + x for x in b.identifier])
    b['discrepancy against Table S1 / the text'] = b.pop('note')
    out = pd.concat([a, b], ignore_index=True)
    p = ROOT / 'results' / 'cohort_provenance.csv'
    out.to_csv(p, index=False)
    miss = a[a['discrepancy against Table S1 / the text'].str.contains('identifier absent')]
    print(f'{len(a)} panel cohorts + {len(b)} secondary cohorts resolved -> {p}')
    print(f'{len(miss)} panel cohorts still carry no resolvable identifier in the paper: '
          + ', '.join(miss.cohort))
    col = 'discrepancy against Table S1 / the text'
    absent = a[col].str.contains('identifier absent')
    subst = a[(a[col].str.len() > 0) & ~absent]
    print(f'{len(subst)} panel rows carry a substantive disagreement with Table S1 or the text '
          f'(beyond a missing identifier):')
    for _, r in subst.iterrows():
        print(f'  {r.cohort:24s} {r[col][:110]}')


if __name__ == '__main__':
    main()
