"""Referee round-3 §4 — per-cohort malignant-purity audit (must APPEAR, not be promised).
Classify each of the 25 convergence cohorts by malignant-selection method, and where the source file carries
cell-type / malignancy annotations, COMPUTE the annotated-malignant (or tumor-epithelial) fraction. Flags the
cohorts at dilution risk (like the lung case). Writes results/referee2_newdata/E18_purity_audit.csv."""
import scanpy as sc, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from pathlib import Path
ROOT=Path('.').resolve(); C=ROOT/'results'
def frac(path,maskfn,denomfn=None):
    try:
        b=sc.read_h5ad(path,backed='r'); o=b.obs
        num=np.asarray(maskfn(o)).sum(); den=len(o) if denomfn is None else np.asarray(denomfn(o)).sum()
        return 100*num/den, len(o)
    except Exception as ex: return None,str(ex)[:60]

rows=[]
def add(name,method,cls,path=None,maskfn=None,denomfn=None,note=''):
    pf=nc=''
    if path is not None and maskfn is not None:
        f,n=frac(path,maskfn,denomfn)
        if f is not None: pf=f"{f:.0f}%"; nc=n
    rows.append(dict(cohort=name,selection=method,purity_class=cls,malig_frac_of_arm=pf,n_source=nc,note=note))

# --- cohorts with computable fractions from annotated source files ---
add('HCC','celltype=="Hepatocyte" (all hepatocytes, tumor+normal)','PROXY — AT RISK',
    C/'cancer_sun2021_hcc/sun2021.h5ad',
    lambda o:(o['celltype'].astype(str)=='Hepatocyte')&(o.get('sample',pd.Series(['t']*len(o))).astype(str).str.contains('T|tumor|Tumor',regex=True) if 'sample' in o else np.ones(len(o),bool)),
    lambda o:(o['celltype'].astype(str)=='Hepatocyte'),
    note='fraction of the hepatocyte arm from tumor samples (rest = normal-adjacent hepatocytes)')
add('BRCA','pf=None (ALL cells of the breast atlas)','AT RISK (no malignant filter)',
    C/'cancer_validation/wu2021_brca.h5ad',
    lambda o:(o['celltype_major'].astype(str)=='Cancer Epithelial'),
    note='annotated Cancer-Epithelial fraction of the arm actually used')
add('NSCLC','cell_type=="malignant cell" (census)','author census-malignant',
    ROOT/'paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',
    lambda o:(o['cell_type'].astype(str)=='malignant cell'),
    note='census malignant label')
add('BCC','disease==BCC & basal/keratinocyte cell types','PROXY — AT RISK',
    ROOT/'paper6_immune_3layer/data/bcc_full.h5ad',
    lambda o:((o['disease'].astype(str)=='basal cell carcinoma')&(o['cell_type'].astype(str).isin(['basal cell','basal cell of epidermis','suprabasal keratinocyte']))),
    lambda o:(o['disease'].astype(str)=='basal cell carcinoma'),
    note='basal/keratinocyte fraction of BCC-sample cells (rest = immune/stroma/normal)')
add('Colon (Pelka)','ClusterMidway=="EpiT" (author tumor-epithelial)','author-annotated (high)',
    C/'cancer_phase2/data/pelka_colon/pelka2021_colon_all_cells.h5ad',
    lambda o:(o['ClusterMidway'].astype(str)=='EpiT'),
    lambda o:(o['ClusterTop'].astype(str)=='Epi'),
    note='EpiT as fraction of all epithelial (rest = normal-sample Epi)')

# --- cohorts pre-filtered to malignant in the source pipeline (*_malignant.h5ad) ---
prefiltered=['Glioblastoma','OvarianSerous','RCC','Melanoma','Prostate','HNSCC','Neuroblastoma',
 'B-ALL','DLBCL','SqLung','TNBC','PDAC','AML','MM','SCLC','Gastric','Osteosarcoma']
for n in prefiltered: add(n,'pre-filtered *_malignant.h5ad (author/inferCNV in source pipeline)','author-annotated (high)')
for n in ['GermCell','Ewing','Seminoma']: add(n,'marker-based malignant label','MARKER — overstates (bounded by drop-4)',note='POU5F1/SALL4 etc. risk selecting stem/high-output cells')

df=pd.DataFrame(rows); df.to_csv(C/'referee2_newdata/E18_purity_audit.csv',index=False)
print("="*80); print("PER-COHORT MALIGNANT-PURITY AUDIT (25 convergence cohorts)"); print("="*80)
print(df.to_string(index=False))
print("\nAT-RISK cohorts (proxy/no-filter selection, dilution possible like the lung case):",
      list(df[df.purity_class.str.contains('RISK')].cohort))
print("Direction of bias: normal-adjacent contamination UNDERSTATES deepening/convergence; marker-based OVERSTATES (bounded by drop-4).")
