"""Referee round-2 ask #1 (SYMMETRIC arm) — HIGH-BASELINE held-out NEGATIVE control, LUNG.
GSE131907 (Kim et al. 2020): author cell annotations give Epithelial cells split by Sample_Origin
(nLung = normal lung, tLung = tumor lung); 10 patients have both (N##/T## same number).
Load ONLY those ~11k epithelial columns from the raw UMI matrix (usecols), project, and test
within-patient: does tumor epithelium deepen the core over matched NORMAL lung epithelium?
Lung is the paper's canonical HIGH-baseline tissue -> prediction: normal lung already high, NO deepening."""
import pandas as pd, numpy as np, scanpy as sc, anndata as ad, json, sys, re, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import wilcoxon
ROOT=Path(__file__).resolve().parents[1]; DD=ROOT/'results/referee2_newdata'
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
refPC=(Rm.values-mu)/sd@load; refPn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)

# --- annotation: epithelial cells in nLung/tLung, paired patients only ---
an=pd.read_csv(DD/'lung_annotation.txt.gz',sep='\t')
epi=an[(an.Cell_type=='Epithelial cells')&(an.Sample_Origin.isin(['nLung','tLung']))].copy()
epi['pat']=epi.Sample.str.extract(r'LUNG_[NT](\d+)')[0]
epi['arm']=epi.Sample.str.extract(r'LUNG_([NT])\d+')[0]   # N=normal, T=tumor
paired=[p for p in epi.pat.unique() if {'N','T'}<=set(epi[epi.pat==p].arm)]
epi=epi[epi.pat.isin(paired)].copy()
idcol='Index' if 'Index' in epi.columns else 'Barcode'
targets=epi[idcol].astype(str).tolist()
print(f"paired lung patients: {sorted(paired)} ; epithelial cells to load: {len(targets)}",flush=True)

# --- stream-load only target columns from the big raw UMI matrix ---
mat=DD/'lung_raw_UMI.txt.gz'
hdr=pd.read_csv(mat,sep='\t',nrows=0); allcols=list(hdr.columns); genecol=allcols[0]
present=[c for c in targets if c in set(allcols)]
print(f"matrix cols: {len(allcols)} (gene col '{genecol}'); target present: {len(present)}/{len(targets)}",flush=True)
use=[genecol]+present
df=pd.read_csv(mat,sep='\t',usecols=use,index_col=0)   # genes x cells
df.index=[str(g).upper() for g in df.index]
if pd.Index(df.index).has_duplicates: df=df.groupby(level=0).sum()
X=df.T                                                  # cells x genes
meta=epi.set_index(idcol).loc[X.index]
print(f"loaded matrix {df.shape} (genes x cells)",flush=True)

def core_of(cells):
    sub=X.loc[cells]
    a=ad.AnnData(sub.values.astype(np.float32)); a.var_names=list(df.index); a.obs_names=list(cells); a.var_names_make_unique()
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    if a.n_obs<1: return 0,np.nan,np.nan
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values
    Z=(S-mu)/sd; pc=Z@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
    return a.n_obs,100*np.isin(t,core_idx).mean(),100*np.isin(t,cyto_ids).mean()

rows=[]
for p in sorted(paired):
    nc=meta.index[(meta.pat==p)&(meta.arm=='N')]; tc=meta.index[(meta.pat==p)&(meta.arm=='T')]
    nN,cN,yN=core_of(nc); nT,cT,yT=core_of(tc)
    rows.append(dict(patient='LUNG'+p,n_epi_normal=nN,n_epi_tumor=nT,normal_core=cN,tumor_core=cT,normal_cyto=yN,tumor_cyto=yT))
    print(f"  LUNG{p}: normal(nLung) epi {nN} core {cN:.1f}% | tumor(tLung) epi {nT} core {cT:.1f}%",flush=True)
df2=pd.DataFrame(rows); df2.to_csv(DD/'E12_lung_negative.csv',index=False)
FLOOR=20; ok=df2[(df2.n_epi_normal>=FLOOR)&(df2.n_epi_tumor>=FLOOR)].copy(); ok['diff']=ok.tumor_core-ok.normal_core
print("\n"+"="*64); print("HIGH-BASELINE NEGATIVE — lung (GSE131907), within-patient nLung(normal) vs tLung(tumor)")
print(f"chance core ~ {100*3/(len(refcells)-len(EXCLUDE)):.1f}%")
print(f"patients with >={FLOOR} epithelial cells in BOTH arms: {len(ok)}/{len(df2)}\n")
if len(ok):
    print(ok[['patient','n_epi_normal','n_epi_tumor','normal_core','tumor_core','diff']].round(1).to_string(index=False))
    print(f"\nmean core: NORMAL lung {ok.normal_core.mean():.1f}%  ->  tumor {ok.tumor_core.mean():.1f}%  (Δ {ok['diff'].mean():+.1f}pp)")
    print(f"tumor deeper than normal in {(ok['diff']>0).sum()}/{len(ok)} patients; paired Wilcoxon p={wilcoxon(ok.tumor_core,ok.normal_core).pvalue:.3f}")
    hi=ok.normal_core.mean()>=10
    print(f"\n[VERDICT] normal-lung baseline {'HIGH' if hi else 'LOW'} ({ok.normal_core.mean():.1f}%); "
          f"{'NEGATIVE CONTROL HOLDS — high-baseline lung, no significant deepening (mirror of the low-baseline gastric result)' if hi and (abs(ok['diff'].mean())<3 or (ok['diff']>0).sum()<=len(ok)/2) else 'see per-patient'}")
print("DONE")
