"""Axis-1 — extend the pre-invasive prediction to two NEW organs (esophagus, breast).
The paper predicts pre-invasive lesions already occupy the convergent region before invasion; tested so far
only in colorectal adenoma (n=3, Fig. S20). Here:
  ESOPHAGUS  GSE291080 (Gier 2025): WITHIN-PATIENT Barrett's metaplasia (be) vs the same patient's normal
             squamous esophagus (es) and stomach (st). Barrett's is columnar metaplasia, so squamous is the
             tissue it replaced and gastric columnar is what it resembles -> we report BOTH comparators.
  BREAST     GSE195861 (Tokura 2022): cross-sectional normal breast -> DCIS -> invasive IDC -> LN metastasis.
Epithelial gate EPCAM+/CD45-, same fixed projection; reports full core and the cytotrophoblast-specific part."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, os, glob, re, tarfile, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
from scipy.stats import wilcoxon, mannwhitneyu
ROOT=Path('.').resolve(); DD=ROOT/'results/referee2_newdata'
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12; refcells=list(Rm.index)
core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
refPC=(Rm.values-mu)/sd@load; refPn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)

def unpack(tar,outdir,pat):
    outdir=Path(outdir); outdir.mkdir(exist_ok=True)
    if any(outdir.iterdir()): return outdir
    with tarfile.open(tar) as t: t.extractall(outdir,members=[m for m in t.getmembers() if pat in m.name or True])
    # organise <sample>/{barcodes,features,matrix}
    for f in glob.glob(str(outdir/'*barcodes.tsv.gz')):
        base=os.path.basename(f); s=re.sub(r'_barcodes.*','',base)
        d=outdir/s; d.mkdir(exist_ok=True)
        for kind in ['barcodes','features','matrix']:
            src=[g for g in glob.glob(str(outdir/f'{s}*{kind}*')) if os.path.isfile(g)]
            if src: os.rename(src[0], d/(f'{kind}.tsv.gz' if kind!='matrix' else 'matrix.mtx.gz'))
    return outdir

def corefrac(d):
    try: a=sc.read_10x_mtx(str(d),var_names='gene_symbols',make_unique=True)
    except Exception as ex: return None
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=500].copy()
    if a.n_obs<20: return None
    def col(g):
        if g not in set(a.var_names): return np.zeros(a.n_obs)
        v=a[:,g].X; return np.asarray(v.todense()).ravel() if hasattr(v,'todense') else np.asarray(v).ravel()
    epi=(col('EPCAM')>0)&(col('PTPRC')==0); ae=a[epi].copy(); del a; gc.collect()
    if ae.n_obs<20: return None
    sc.pp.normalize_total(ae,target_sum=1e4); sc.pp.log1p(ae)
    S=pct.score_modules(ae,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values
    Z=(S-mu)/sd; pc=Z@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
    sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
    n=ae.n_obs; del ae; gc.collect()
    return dict(n_epi=n,core=100*np.isin(t,core_idx).mean(),cyto=100*np.isin(t,cyto_ids).mean())

# ---------------- ESOPHAGUS (within-patient) ----------------
od=unpack(DD/'GSE291080_BARRETT.tar',DD/'barrett291080','')
rows=[]
for d in sorted(glob.glob(str(od/'GSM*'))):
    nm=os.path.basename(d); m=re.search(r'GSM\d+_(a\d+)_([a-z]+)\d*',nm)
    if not m: continue
    r=corefrac(d)
    if r: rows.append(dict(patient=m.group(1),tissue=m.group(2),sample=nm,**r)); print(f"  {nm}: {r['n_epi']} epi, core {r['core']:.1f}%",flush=True)
be=pd.DataFrame(rows); be.to_csv(DD/'E22_barrett.csv',index=False)
print("\n"+"="*66); print("ESOPHAGUS — Barrett's metaplasia (pre-invasive) vs same patient's normal (GSE291080)")
agg=be.groupby(['patient','tissue']).core.mean().unstack()
print(agg.round(1).to_string())
for comp,label in [('es','normal squamous esophagus'),('st','stomach (columnar comparator)')]:
    if comp in agg.columns and 'be' in agg.columns:
        ok=agg[['be',comp]].dropna()
        if len(ok):
            print(f"\n  Barrett vs {label}: {len(ok)} paired patients; mean {ok['be'].mean():.1f}% vs {ok[comp].mean():.1f}% (Δ {ok['be'].mean()-ok[comp].mean():+.1f}pp); higher in {(ok['be']>ok[comp]).sum()}/{len(ok)}")
            if len(ok)>=5:
                try: print(f"     paired Wilcoxon p={wilcoxon(ok['be'],ok[comp]).pvalue:.3f}")
                except Exception as ex: print("     ",ex)

# ---------------- BREAST (cross-sectional ladder) ----------------
od2=unpack(DD/'GSE195861_DCIS.tar',DD/'dcis195861','')
STAGE=lambda s: ('normal' if 'NCCBC13' in s else 'DCIS' if 'NCCBC' in s else 'IDC' if re.search(r'_P\d',s) else 'LN-met' if re.search(r'_M\d',s) else '?')
rows=[]
for d in sorted(glob.glob(str(od2/'GSM*'))):
    nm=os.path.basename(d); r=corefrac(d)
    if r: rows.append(dict(sample=nm,stage=STAGE(nm),**r)); print(f"  {nm} [{STAGE(nm)}]: {r['n_epi']} epi, core {r['core']:.1f}%",flush=True)
br=pd.DataFrame(rows); br.to_csv(DD/'E22_dcis.csv',index=False)
print("\n"+"="*66); print("BREAST — normal -> DCIS -> IDC -> metastasis (GSE195861, cross-sectional)")
for st in ['normal','DCIS','IDC','LN-met']:
    g=br[br.stage==st]
    if len(g): print(f"  {st:8s} n={len(g)} samples: core {g.core.mean():5.1f}% (cyto {g.cyto.mean():4.1f}%)  per-sample {[round(x,1) for x in g.core]}")
if len(br[br.stage=='DCIS']) and len(br[br.stage=='IDC']):
    dq=br[br.stage=='DCIS'].core; iq=br[br.stage=='IDC'].core
    print(f"\n  DCIS reaches {100*dq.mean()/iq.mean():.0f}% of the invasive (IDC) level; Mann-Whitney DCIS vs IDC p={mannwhitneyu(dq,iq).pvalue:.3f}")
    if len(br[br.stage=='normal']): print(f"  normal breast reference: {br[br.stage=='normal'].core.mean():.1f}%")
print("\nDONE")
