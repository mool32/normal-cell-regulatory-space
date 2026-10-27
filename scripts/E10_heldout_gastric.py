"""Referee round-2 ask #1 (out-of-sample held-out low-baseline deepening), GASTRIC.
GSE150290 (Kim et al. 2022): 24 patients, paired Pat##-A = ADJACENT NORMAL stomach,
Pat##-B = CANCER stomach (verified on GEO GSM4546300/301 — the reverse of a naming guess).
Processed per-sample gene x cell matrices (.txt.gz). Gate epithelial cells (EPCAM+/CD45-),
project into the fixed single-cell eigenspace, and test within-patient deepening
(cancer core fraction > matched-normal core fraction). Floor: >=20 epithelial cells per arm."""
import scanpy as sc, anndata as ad, json, numpy as np, pandas as pd, sys, glob, os, re, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import wilcoxon
ROOT=Path(__file__).resolve().parents[1]
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
refPC=(Rm.values-mu)/sd@load; refPn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)

DDIR=ROOT/'results/referee2_newdata/gastric150290'
def load_txt(path):
    df=pd.read_csv(path,sep='\t',index_col=0)          # genes x cells
    a=ad.AnnData(df.T.values.astype(np.float32))       # cells x genes
    a.var_names=[str(g).upper() for g in df.index]; a.obs_names=list(df.columns); a.var_names_make_unique()
    return a
def core_frac(a):
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    def col(g): return np.asarray(a[:,g].X).ravel() if g in set(a.var_names) else np.zeros(a.n_obs)
    epi=(col('EPCAM')>0)&(col('PTPRC')==0)
    ae=a[epi].copy()
    if ae.n_obs<1: return 0,np.nan,np.nan
    sc.pp.normalize_total(ae,target_sum=1e4); sc.pp.log1p(ae)
    S=pct.score_modules(ae,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values
    Z=(S-mu)/sd; pc=Z@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
    return ae.n_obs, 100*np.isin(t,core_idx).mean(), 100*np.isin(t,cyto_ids).mean()

files={}
for f in glob.glob(str(DDIR/'*.txt.gz')):
    m=re.search(r'Pat(\d+)-([AB])',os.path.basename(f))
    if m: files[(m.group(1),m.group(2))]=f
pats=sorted(set(p for p,_ in files))
rows=[]
for p in pats:
    if (p,'A') not in files or (p,'B') not in files: continue
    try:
        nA,cA,yA=core_frac(load_txt(files[(p,'A')]))   # A = normal
        nB,cB,yB=core_frac(load_txt(files[(p,'B')]))   # B = cancer
        rows.append(dict(patient=p,n_epi_normal=nA,n_epi_cancer=nB,normal_core=cA,cancer_core=cB,normal_cyto=yA,cancer_cyto=yB))
        print(f"  Pat{p}: normal epi {nA} core {cA:.1f}% | cancer epi {nB} core {cB:.1f}%",flush=True)
    except Exception as ex: print(f"  !!Pat{p}: {ex}",flush=True)
df=pd.DataFrame(rows); df.to_csv(ROOT/'results/referee2_newdata/E10_gastric_deepening.csv',index=False)

FLOOR=20; ok=df[(df.n_epi_normal>=FLOOR)&(df.n_epi_cancer>=FLOOR)].copy()
ok['diff']=ok.cancer_core-ok.normal_core
print("\n"+"="*64); print(f"HELD-OUT GASTRIC (GSE150290) — within-patient deepening (A=normal, B=cancer)")
print(f"chance core ~ {100*3/(len(refcells)-len(EXCLUDE)):.1f}%")
print(f"patients with >={FLOOR} epithelial cells in BOTH arms: {len(ok)}/{len(df)}\n")
if len(ok):
    print(ok[['patient','n_epi_normal','n_epi_cancer','normal_core','cancer_core','diff']].round(1).to_string(index=False))
    print(f"\nmean core: normal {ok.normal_core.mean():.1f}%  ->  cancer {ok.cancer_core.mean():.1f}%  (Δ {ok['diff'].mean():+.1f}pp)")
    print(f"cancer deeper than normal in {(ok['diff']>0).sum()}/{len(ok)} patients")
    try:
        w=wilcoxon(ok.cancer_core,ok.normal_core); print(f"paired Wilcoxon p={w.pvalue:.4f}")
    except Exception as ex: print("wilcoxon:",ex)
    print(f"cyto-specific: normal {ok.normal_cyto.mean():.1f}% -> cancer {ok.cancer_cyto.mean():.1f}%")
    print(f"\n[VERDICT] {'out-of-sample deepening REPLICATES (gastric is low-baseline, cancer deepens)' if ok['diff'].mean()>0 and (ok['diff']>0).sum()>len(ok)/2 else 'NUANCED — see per-patient'}")
else:
    print("Too few patients meet the epithelial floor — cohort too sparse for a powered within-patient test.")
print("DONE")
