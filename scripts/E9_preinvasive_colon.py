"""Referee round-2 ask #4 (pre-invasive prediction), colorectal.
GSE161277 (Zhang et al.): within-patient normal -> adenoma -> carcinoma scRNA (10x).
Gate epithelial cells (EPCAM+ / CD45-), project into the fixed single-cell eigenspace, and ask:
does ADENOMA (pre-invasive) epithelium already occupy the cytotrophoblast/megakaryocyte core
region at a fraction close to carcinoma and above matched normal?  Patients 1,2,3 have the full
normal/adenoma/carcinoma series.  Decomposes core into cytotrophoblast (specific) vs megakaryocyte."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, glob, os, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CORE=['megakaryocytes','cytotrophoblasts','migrating cytotrophoblasts']
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refcells=list(Rm.index); core_idx=np.array([refcells.index(c) for c in CORE])
cyto_ids=[refcells.index('cytotrophoblasts'),refcells.index('migrating cytotrophoblasts')]; mega_id=refcells.index('megakaryocytes')
EXCLUDE=np.array([refcells.index(c) for c in ['erythrocytes','platelets'] if c in refcells])
refPC=(Rm.values-mu)/sd@load; refPn=refPC/(np.linalg.norm(refPC,axis=1,keepdims=True)+1e-12)

DDIR=ROOT/'results/referee2_newdata/colon161277'
rows=[]
for d in sorted(glob.glob(str(DDIR/'GSM*'))):
    name=os.path.basename(d)
    parts=name.split('_'); pat=parts[1]; stage='_'.join(parts[2:])
    if stage=='blood': continue
    try:
        a=sc.read_10x_mtx(d,var_names='gene_symbols',make_unique=True)
        a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
        ng=np.asarray((a.X>0).sum(1)).ravel(); a=a[ng>=1000].copy()
        # epithelial gate: EPCAM detected, CD45/PTPRC not detected
        def col(g):
            return np.asarray(a[:,g].X.todense()).ravel() if g in a.var_names else np.zeros(a.n_obs)
        epcam=col('EPCAM'); ptprc=col('PTPRC')
        epi=(epcam>0)&(ptprc==0)
        ae=a[epi].copy()
        if ae.n_obs<20:
            rows.append(dict(patient=pat,stage=stage,n_epi=int(epi.sum()),core=np.nan,cyto=np.nan,mega=np.nan)); continue
        sc.pp.normalize_total(ae,target_sum=1e4); sc.pp.log1p(ae)
        S=pct.score_modules(ae,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values
        Z=(S-mu)/sd; pc=Z@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
        sim=pcn@refPn.T; sim[:,EXCLUDE]=-1e9; t=np.argmax(sim,axis=1)
        rows.append(dict(patient=pat,stage=stage,n_epi=int(ae.n_obs),
                         core=100*np.isin(t,core_idx).mean(),
                         cyto=100*np.isin(t,cyto_ids).mean(),
                         mega=100*(t==mega_id).mean()))
        print(f"  {name}: {epi.sum()} epithelial / {a.n_obs} QC cells",flush=True)
    except Exception as ex: print(f"  !!{name}: {ex}",flush=True)
df=pd.DataFrame(rows); df.to_csv(ROOT/'results/referee2_newdata/E9_colon_preinvasive.csv',index=False)
print("\n"+"="*64); print("PRE-INVASIVE (GSE161277) — epithelial core fraction by patient x stage")
print(f"chance core ~ {100*3/(len(refcells)-len(EXCLUDE)):.1f}%\n")
piv=df.pivot_table(index='patient',columns='stage',values='core')
order=[c for c in ['normal','para-cancer','adenoma','adenoma_1','adenoma_2','carcinoma'] if c in piv.columns]
print("[CORE %]"); print(piv[order].round(1).to_string())
print("\n[CYTO %]"); print(df.pivot_table(index='patient',columns='stage',values='cyto')[order].round(1).to_string())
# within-patient normal->adenoma->carcinoma for patients with all three
print("\n[within-patient progression, core %]:")
for pat,g in df.groupby('patient'):
    gg={r.stage:r.core for r in g.itertuples()}
    ad=np.nanmean([gg.get('adenoma',np.nan),gg.get('adenoma_1',np.nan),gg.get('adenoma_2',np.nan)])
    if 'normal' in gg and not np.isnan(ad) and 'carcinoma' in gg:
        print(f"   {pat}: normal {gg['normal']:.1f}%  ->  adenoma {ad:.1f}%  ->  carcinoma {gg['carcinoma']:.1f}%")
# aggregate test: adenoma vs normal (paired within patient)
pats=[p for p,g in df.groupby('patient') if 'normal' in set(g.stage) and (set(g.stage)&{'adenoma','adenoma_1','adenoma_2'})]
nb=[]; ab=[]; cb=[]
for pat in pats:
    g=df[df.patient==pat]; d={r.stage:r.core for r in g.itertuples()}
    nb.append(d['normal']); ab.append(np.nanmean([d.get(s,np.nan) for s in ['adenoma','adenoma_1','adenoma_2']]))
    cb.append(d.get('carcinoma',np.nan))
nb,ab,cb=np.array(nb),np.array(ab),np.array(cb)
print(f"\n[SUMMARY over {len(pats)} full-progression patients]  mean core: normal {np.nanmean(nb):.1f}%  adenoma {np.nanmean(ab):.1f}%  carcinoma {np.nanmean(cb):.1f}%")
print(f"   adenoma>normal in {(ab>nb).sum()}/{len(pats)} patients; adenoma reaches {100*np.nanmean(ab)/np.nanmean(cb):.0f}% of the carcinoma level")
print(f"   => prediction (pre-invasive adenoma already occupies the region) {'SUPPORTED' if np.nanmean(ab)>1.5*np.nanmean(nb) and (ab>nb).sum()>=len(pats)-1 else 'NUANCED — see per-patient'}")
print("DONE")
