import sys, json, re
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at HALF the printed width: 0.50 x 6.5 in (sits beside its partner in Fig. S5)
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix()
sys.path.insert(0,ROOT+'/tool/perceptome2'); import perceptome as pct
e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json'))
mods=e['module_order']; load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
hpa=pct.load_hpa_perceptivity(mode='single_cell_scaled'); R=hpa['R'][mods]
mu=R.mean(0).values; sd=R.std(0).values+1e-12; ref=(R.values-mu)/sd@load; idx=R.index.tolist()
refn=ref/(np.linalg.norm(ref,axis=1,keepdims=True)+1e-12)
DUET=['megakaryocytes','cytotrophoblasts']
NEURAL=[c for c in ['brain excitatory neurons','brain inhibitory neurons','adrenal medulla cells'] if c in idx]
print('neural ref:',NEURAL)
duet_c=ref[[idx.index(c) for c in DUET]].mean(0); duet_u=duet_c/np.linalg.norm(duet_c)
neur_c=ref[[idx.index(c) for c in NEURAL]].mean(0); neur_u=neur_c/np.linalg.norm(neur_c)
def cos_axes(P):
    Pn=P/(np.linalg.norm(P,axis=1,keepdims=True)+1e-12)
    return Pn@duet_u, Pn@neur_u
def nv(a):
    if 'feature_name' in a.var.columns and not any(c.isalpha() for c in str(a.var.index[0])): a.var_names=a.var['feature_name'].astype(str).values
    a.var_names=a.var_names.str.upper(); a.var_names_make_unique()
    base=pd.Series([re.sub(r'-\d+$','',str(v)) for v in a.var_names]); k=~base.duplicated().values
    a=a[:,k].copy(); a.var_names=base[k].values; return a
def proj(path,n=2000):
    a=sc.read_h5ad(path); a=nv(a); ng=np.asarray((a.X>0).sum(1)).flatten(); a=a[ng>=1000].copy()
    if a.n_obs>n: a=a[np.random.default_rng(0).choice(a.n_obs,n,replace=False)].copy()
    if a.X.max()>100: sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    s=pct.score_modules(a,method='mean_raw'); s=s['scores'] if isinstance(s,dict) else s
    z=(s[mods].values-mu)/sd; return z@load
bx,by=cos_axes(proj(ROOT+'/results/cancer_v6_4_extras/B-ALL_malignant.h5ad'))
nx,ny=cos_axes(proj(ROOT+'/results/cancer_v6_4_extras/Neuroblastoma_malignant.h5ad'))
rx,ry=cos_axes(ref)
print(f'B-ALL  duet-cos median={np.median(bx):.2f}  neural-cos median={np.median(by):.2f}')
print(f'Neurobl duet-cos median={np.median(nx):.2f}  neural-cos median={np.median(ny):.2f}')
# RENDER
fig,ax=plt.subplots(figsize=F.panel_size(0.5, 1.05))
ax.scatter(rx,ry,s=30,c='#d9d9d9',edgecolor='#bbb',lw=.3,zorder=1,label='154 normal cell types')
ei=[idx.index(c) for c in DUET]; ax.scatter(rx[ei],ry[ei],s=180,marker='*',c='#c0392b',edgecolor='k',lw=.5,zorder=6,label='convergence target')
ni=[idx.index(c) for c in NEURAL]; ax.scatter(rx[ni],ry[ni],s=60,marker='^',c='#16a085',edgecolor='k',lw=.4,zorder=6,label='neural-program cells')
ax.scatter(bx,by,s=7,c='#2980b9',alpha=.30,zorder=3,label='B-ALL malignant (converges, 61%)')
ax.scatter(nx,ny,s=7,c='#8e44ad',alpha=.30,zorder=2,label='neuroblastoma malignant (7%)')
ax.set_xlabel('cosine to the convergent state (9-D)',fontsize=BASE*0.95)
ax.set_ylabel('cosine to the neural region (9-D)',fontsize=BASE*0.95)
ax.set_title('',fontsize=BASE*1.121,loc='left',weight='bold')
ax.legend(fontsize=BASE*0.916,loc='lower left',framealpha=.95)
ax.text(0.99,0.02,'Axes are interpretable 9-D directions (not lossy PC1–PC2): B-ALL shifts onto the\ncore target; neuroblastoma retains its neural-crest identity. Shows convergence directly.',
        transform=ax.transAxes,ha='right',va='bottom',fontsize=BASE*0.825,style='italic',color='#555')
plt.tight_layout(); plt.savefig(ROOT+'/manuscript/figures/v7_fig1A_clouds.png',dpi=200); plt.close(); print('saved clouds fig')

# ---- v2: KDE density contours instead of overlapping scatter ----
from scipy.stats import gaussian_kde
fig,ax=plt.subplots(figsize=F.panel_size(0.5, 1.05))
ax.scatter(rx,ry,s=24,c='#dddddd',edgecolor='#bbb',lw=.25,zorder=1,label='154 normal cell types')
def contour(xx,yy,color,label,zo):
    k=gaussian_kde(np.vstack([xx,yy]))
    xg,yg=np.mgrid[-1:1:120j,-1:1:120j]; dd=k(np.vstack([xg.ravel(),yg.ravel()])).reshape(xg.shape)
    levels=np.percentile(k(np.vstack([xx,yy])),[20,50,80])
    ax.contourf(xg,yg,dd,levels=[levels[0],dd.max()],colors=[color],alpha=.18,zorder=zo)
    ax.contour(xg,yg,dd,levels=levels,colors=color,linewidths=[.6,1.0,1.4],zorder=zo+1)
    ax.plot([],[],color=color,lw=2,label=label)
contour(bx,by,'#2980b9','B-ALL (61% on target)',2)
contour(nx,ny,'#8e44ad','neuroblastoma (7%, stays neural)',2)
ei=[idx.index(c) for c in DUET]; ax.scatter(rx[ei],ry[ei],s=200,marker='*',c='#c0392b',edgecolor='k',lw=.6,zorder=8,label='convergence target')
ni=[idx.index(c) for c in NEURAL]; ax.scatter(rx[ni],ry[ni],s=70,marker='^',c='#16a085',edgecolor='k',lw=.5,zorder=8,label='neural-program cells')
ax.set_xlabel('cosine to the convergent state (9-D)',fontsize=BASE*0.95)
ax.set_ylabel('cosine to the neural region (9-D)',fontsize=BASE*0.95)
ax.set_title('',fontsize=BASE*1.121,loc='left',weight='bold')
ax.legend(fontsize=BASE*0.78,loc='lower left',framealpha=.95,handlelength=1.2,borderpad=0.35,labelspacing=0.35)
ax.set_xlim(-1,1.05); ax.set_ylim(-1,1.05)
plt.tight_layout(); plt.savefig(ROOT+'/manuscript/figures/v7_fig1A_clouds.png',dpi=200); plt.close(); print('saved clouds v2 (contours)')
