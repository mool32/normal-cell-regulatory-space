"""Referee 3.5, out-of-sample: does the SHAPE-ONLY within-patient deepening also hold in the held-out
gastric cohort and in the corrected high-baseline lung cohort? (magnitude removed = per-cell centering)."""
import scanpy as sc, anndata as ad, json, numpy as np, pandas as pd, sys, glob, os, re, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import wilcoxon
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
core_idx=np.array([refcells.index(c) for c in CORE]); EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
demag=lambda X:X-X.mean(1,keepdims=True); ident=lambda X:X
def refset(tf):
    R=tf(Rm.values.astype(float)); mu=R.mean(0); sd=R.std(0)+1e-12
    rp=((R-mu)/sd)@load; return mu,sd,rp/(np.linalg.norm(rp,axis=1,keepdims=True)+1e-12)
REF={'RAW':refset(ident),'SHAPE-ONLY':refset(demag)}; TF={'RAW':ident,'SHAPE-ONLY':demag}
def corefrac(a,tag):
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    def col(g): return np.asarray(a[:,g].X.todense()).ravel() if (g in set(a.var_names) and hasattr(a.X,'todense')) else (np.asarray(a[:,g].X).ravel() if g in set(a.var_names) else np.zeros(a.n_obs))
    if tag=='epcam': epi=(col('EPCAM')>0)&(col('PTPRC')==0); a=a[epi].copy()
    if a.n_obs<20: return None
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
    out={}
    for m in ['RAW','SHAPE-ONLY']:
        mu,sd,rpn=REF[m]; pc=((TF[m](S)-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        sim=pcn@rpn.T; sim[:,EXCLUDE]=-1e9; out[m]=100*np.isin(np.argmax(sim,axis=1),core_idx).mean()
    out['n']=a.n_obs; return out

def report(name,rows):
    df=pd.DataFrame(rows); print(f"\n=== {name} ({len(df)} patients, both arms >=20 epi) ===")
    for m in ['RAW','SHAPE-ONLY']:
        dif=df[f'malig_{m}'].values-df[f'normal_{m}'].values
        print(f"  [{m:10s}] malig {df[f'malig_{m}'].mean():5.1f}% vs normal {df[f'normal_{m}'].mean():5.1f}% -> Δ {dif.mean():+5.1f}pp, deeper {int((dif>0).sum())}/{len(df)}, p={wilcoxon(df[f'malig_{m}'],df[f'normal_{m}']).pvalue:.2e}")
    rd=(df['malig_RAW']-df['normal_RAW']).mean(); sd_=(df['malig_SHAPE-ONLY']-df['normal_SHAPE-ONLY']).mean()
    print(f"  shape-only retains {100*sd_/rd:.0f}% -> {'SURVIVES' if sd_>0 and ((df['malig_SHAPE-ONLY']>df['normal_SHAPE-ONLY']).mean()>0.6) else 'collapses'}")
    return df

# ---- GASTRIC (A=normal, B=cancer) ----
def load_txt(p):
    d=pd.read_csv(p,sep='\t',index_col=0); x=ad.AnnData(d.T.values.astype(np.float32)); x.var_names=[str(g).upper() for g in d.index]; x.var_names_make_unique(); return x
gf={}
for f in glob.glob(str(DD/'gastric150290/*_Pat*-[AB].txt.gz')):
    mm=re.search(r'Pat(\d+)-([AB])',os.path.basename(f))
    if mm: gf[(mm.group(1),mm.group(2))]=f
grows=[]
for p in sorted(set(k for k,_ in gf)):
    if (p,'A') not in gf or (p,'B') not in gf: continue
    n=corefrac(load_txt(gf[(p,'A')]),'epcam'); t=corefrac(load_txt(gf[(p,'B')]),'epcam')
    if n and t: grows.append({'patient':p,'normal_RAW':n['RAW'],'malig_RAW':t['RAW'],'normal_SHAPE-ONLY':n['SHAPE-ONLY'],'malig_SHAPE-ONLY':t['SHAPE-ONLY']})
report('GASTRIC held-out (GSE150290)',grows)

# ---- LUNG (nLung=normal, tLung=tumor), from cached matrix via annotation ----
an=pd.read_csv(DD/'lung_annotation.txt.gz',sep='\t'); epi=an[(an.Cell_type=='Epithelial cells')&(an.Sample_Origin.isin(['nLung','tLung']))].copy()
epi['pat']=epi.Sample.str.extract(r'LUNG_[NT](\d+)')[0]; epi['ar']=epi.Sample.str.extract(r'LUNG_([NT])\d+')[0]
paired=[p for p in epi.pat.unique() if {'N','T'}<=set(epi[epi.pat==p].ar)]; epi=epi[epi.pat.isin(paired)]
idc='Index' if 'Index' in epi.columns else 'Barcode'; targ=set(epi[idc].astype(str))
hdr=pd.read_csv(DD/'lung_raw_UMI.txt.gz',sep='\t',nrows=0); gc0=list(hdr.columns)[0]
use=[gc0]+[c for c in list(hdr.columns) if c in targ]
M=pd.read_csv(DD/'lung_raw_UMI.txt.gz',sep='\t',usecols=use,index_col=0); M.index=[str(g).upper() for g in M.index]
if pd.Index(M.index).has_duplicates: M=M.groupby(level=0).sum()
X=M.T; meta=epi.set_index(idc).loc[[c for c in X.index]]
lrows=[]
for p in sorted(paired):
    nc=meta.index[(meta.pat==p)&(meta.ar=='N')]; tc=meta.index[(meta.pat==p)&(meta.ar=='T')]
    def mk(cells):
        s=X.loc[cells]; x=ad.AnnData(s.values.astype(np.float32)); x.var_names=list(M.index); x.var_names_make_unique(); return x
    n=corefrac(mk(nc),'all'); t=corefrac(mk(tc),'all')
    if n and t: lrows.append({'patient':p,'normal_RAW':n['RAW'],'malig_RAW':t['RAW'],'normal_SHAPE-ONLY':n['SHAPE-ONLY'],'malig_SHAPE-ONLY':t['SHAPE-ONLY']})
report('LUNG corrected high-baseline (GSE131907)',lrows)
print("\nDONE")
