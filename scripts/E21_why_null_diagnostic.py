"""WHY did the causal Perturb-seq test (E20) return null? Decompose where the signal dies, in 3 levels:
  L1  Did CRISPRi work?          -> is the TARGET GENE's own mRNA reduced in its KD cells vs control?
  L2  Did the MODULE move?       -> is the score of the module CONTAINING that gene reduced?
                                    (the perceptome averages 5-20 genes per module: a 1-gene KD is diluted ~1/n)
  L3  Did the COORDINATE move?   -> shape-only core proximity (E20 said: no)
If L1 yes, L2 no  -> INSTRUMENT-LIMITED: the perceptome is a coarse-grainer, blind to single-gene perturbation
                     by construction (its virtue as a state sensor is its vice as a perturbation readout).
If L1 yes, L2 yes, L3 no -> BIOLOGY: the state is buffered; single brakes do not move the cell.
Uses the genome-wide Replogle K562 gwps file already on disk."""
import scanpy as sc, json, numpy as np, pandas as pd, sys, warnings, gc
warnings.filterwarnings('ignore'); sys.path.insert(0,'tool/perceptome2'); import perceptome as pct
from perceptome.perceptivity.reference import load_hpa_perceptivity
from pathlib import Path
ROOT=Path('.').resolve()
e=json.load(open('tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; refcells=list(Rm.index)
cyto=refcells.index('cytotrophoblasts'); mega=refcells.index('megakaryocytes')
demag=lambda X:X-X.mean(1,keepdims=True)
Rd=demag(Rm.values.astype(float)); mu=Rd.mean(0); sd=Rd.std(0)+1e-12
refP=((Rd-mu)/sd)@load; refPn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
core_dir=(refPn[cyto]+refPn[mega]); core_dir/=np.linalg.norm(core_dir)
# gene -> modules containing it, and module sizes
gene2mod={}; modsize={}
for m in mods:
    gs=[str(g).upper() for g in pct.get_genes(m)]; modsize[m]=len(gs)
    for g in gs: gene2mod.setdefault(g,[]).append(m)

H5=ROOT/'paper6/data/perturbseq/ReplogleWeissman2022_K562_gwps.h5ad'
b=sc.read_h5ad(H5,backed='r')
pert=b.obs['perturbation'].astype(str).values; nper=b.obs['nperts'].astype(int).values
targets=sorted(set(pert[nper==1]) & set(gene2mod.keys()))
print(f"gwps {b.n_obs} cells; KD genes that are perceptome module members: {len(targets)}",flush=True)
# take those KD cells + a control sample
mask=np.isin(pert,targets)&(nper==1)
ci=np.where(nper==0)[0]; rng=np.random.default_rng(0)
cs=rng.choice(ci,min(20000,len(ci)),replace=False); mask[cs]=True
print(f"loading {int(mask.sum())} cells ({len(cs)} controls)",flush=True)
a=b[mask].to_memory(); del b; gc.collect()
a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
if a.X.max()>50: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
S=pct.score_modules(a,method='mean_raw',gene_set='core')['scores'].reindex(columns=mods,fill_value=0.0).values.astype(float)
Z=demag(S); pc=((Z-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12); prox=pcn@core_dir
p2=a.obs['perturbation'].astype(str).values; n2=a.obs['nperts'].astype(int).values
ctrl=(n2==0); X=a.X
def gene_expr(g):
    if g not in set(a.var_names): return None
    v=X[:,list(a.var_names).index(g)]
    return np.asarray(v.todense()).ravel() if hasattr(v,'todense') else np.asarray(v).ravel()

rows=[]
for g in targets:
    km=(p2==g)&(n2==1)
    if km.sum()<25: continue
    ge=gene_expr(g)
    l1=np.nan
    if ge is not None:
        c_mu,c_sd=ge[ctrl].mean(),ge[ctrl].std()+1e-9
        l1=(ge[km].mean()-c_mu)/c_sd          # target-gene knockdown, in control SD (expect strongly negative)
    for m in gene2mod[g]:
        j=mods.index(m); ms=S[:,j]
        l2=(ms[km].mean()-ms[ctrl].mean())/(ms[ctrl].std()+1e-9)   # module score shift, control SD
        l3=(prox[km].mean()-prox[ctrl].mean())/(prox[ctrl].std()+1e-9)
        rows.append(dict(gene=g,module=m,n_cells=int(km.sum()),module_size=modsize[m],
                         L1_target_mrna_z=l1,L2_module_z=l2,L3_coord_z=l3))
df=pd.DataFrame(rows); df.to_csv(ROOT/'results/referee2_newdata/E21_why_null.csv',index=False)
print("\n"+"="*74); print("WHERE DOES THE CAUSAL SIGNAL DIE? (levels in control-SD units)")
d1=df.drop_duplicates('gene')
print(f"L1  CRISPRi knocked the target gene down:   median z = {d1.L1_target_mrna_z.median():+.2f}  "
      f"({100*(d1.L1_target_mrna_z<-0.2).mean():.0f}% of genes below -0.2 SD)   [expect strongly negative]")
print(f"L2  the module containing it moved:         median z = {df.L2_module_z.median():+.3f}  "
      f"(|z|>0.1 in {100*(df.L2_module_z.abs()>0.1).mean():.0f}% of gene-module pairs)")
print(f"L3  the shape-only coordinate moved:        median z = {df.L3_coord_z.median():+.3f}  "
      f"(|z|>0.1 in {100*(df.L3_coord_z.abs()>0.1).mean():.0f}%)")
att=df.L2_module_z.abs().median()/(d1.L1_target_mrna_z.abs().median()+1e-9)
print(f"\nATTENUATION L1->L2 (module averaging): |L2|/|L1| = {att:.3f}  "
      f"(median module size {df.module_size.median():.0f} genes -> naive dilution ~{1/df.module_size.median():.3f})")
print(f"ATTENUATION L2->L3 (coordinate):        |L3|/|L2| = {df.L3_coord_z.abs().median()/(df.L2_module_z.abs().median()+1e-9):.3f}")
big=df.reindex(df.L2_module_z.abs().sort_values(ascending=False).index).head(8)
print("\nlargest module responses (gene KD -> its own module):")
for r in big.itertuples(): print(f"   {r.gene:9s} -> {r.module:16s} L1 {r.L1_target_mrna_z:+.2f} | L2 {r.L2_module_z:+.3f} | L3 {r.L3_coord_z:+.3f} (module {r.module_size} genes, n={r.n_cells})")
print("\nVERDICT:", "INSTRUMENT-LIMITED — CRISPRi works but module averaging absorbs it (perceptome is a coarse-grainer)"
      if (d1.L1_target_mrna_z.median()<-0.2 and df.L2_module_z.abs().median()<0.1)
      else "modules do respond -> the null is biological (state buffered), see numbers")
print("DONE")
