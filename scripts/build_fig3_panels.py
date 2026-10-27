"""Fig. 3 panels (v11, house style). Per the figure plan:
 B — the cross-cancer recurrence ranking WITHOUT the erythrocyte (a quarantined post-cellular artifact)
     and WITHOUT the megakaryocyte, both of which the complexity control rejects (panel D); the full
     ranking including them stays in the SI.
 D — the transcriptome-complexity control (previously SI Fig. S7A), promoted to the main figure because it
     is what separates the cytotrophoblast from the two low-complexity candidates."""
import sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
import matplotlib.pyplot as plt
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
sys.path.insert(0, ROOT+'/tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
BASE=F.use()

e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]; cells=list(Rm.index)
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
refP=((Rm.values-mu)/sd)@load; refn=refP/(np.linalg.norm(refP,axis=1,keepdims=True)+1e-12)
EX=[cells.index(c) for c in ('erythrocytes','platelets')]
d=np.load(ROOT+'/results/E7_S_raw.npz',allow_pickle=True); S=d['S'].astype(float); cid=d['cid']
uc=sorted(set(cid.tolist()))
pc=((S-mu)/sd)@load; pcn=pc/(np.linalg.norm(pc,axis=1,keepdims=True)+1e-12)
sim=pcn@refn.T; sim[:,EX]=-1e9; top=np.argmax(sim,axis=1)
rec=np.zeros(len(cells))
for c in uc:
    t=top[cid==c]; b=np.bincount(t,minlength=len(cells)); b[EX]=-1
    rec[np.argsort(-b)[:5]]+=1
rec[EX]=-1

# ---------------- B: recurrence ranking, artifacts removed ----------------
DROP={'erythrocytes','platelets','megakaryocytes'}
keep=[i for i in np.argsort(-rec) if cells[i] not in DROP and rec[i]>0][:13]
names=[cells[i] for i in keep]; vals=[rec[i] for i in keep]
SHORT={'cytotrophoblasts':'cytotrophoblast','gastric progenitor cells':'gastric progenitor',
 'late spermatids':'late spermatid','gastric chief cells':'gastric chief',
 'esophageal apical cells':'esophageal apical','hofbauer cells':'hofbauer',
 'pancreatic acinar cells':'pancreatic acinar','migrating cytotrophoblasts':'migrating cytotroph.',
 'differentiating spermatogonia':'diff. spermatogonia','epididymal clear cells':'epididymal clear',
 'early primary spermatocytes':'early spermatocyte','enteric stem cells':'enteric stem',
 'syncytiotrophoblasts':'syncytiotrophoblast'}
FSB=F.panel_size(0.44,1.10)
fig,ax=plt.subplots(figsize=FSB)
col=['#8e1b1b' if n=='cytotrophoblasts' else '#e8a87c' for n in names]
ax.bar(range(len(keep)),vals,color=col,edgecolor='black',lw=1.0,zorder=3)
ax.set_xticks(range(len(keep)))
ax.set_xticklabels([SHORT.get(n,n.replace(' cells','')) for n in names],
                   rotation=42,ha='right',fontsize=BASE*0.80)
# the explanation now lives in the caption; as in-panel text it sat on top of the bars
ax.set_ylabel('cross-cancer recurrence\n(of the 25 cancers)',fontsize=BASE*0.92)
ax.set_ylim(0,max(vals)*1.08)
# recurrence is a count of cancers, so the ticks must be integers (referee item 14)
import matplotlib.ticker as _tk; ax.yaxis.set_major_locator(_tk.MaxNLocator(integer=True))
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(axis='y',alpha=.22,zorder=0)
plt.savefig(FIG+'/v11_fig3B_recurrence.png',dpi=200,bbox_inches='tight'); plt.close()
F.report('Fig 3B recurrence',FSB,0.44)
print('wrote v11_fig3B_recurrence.png |',names[:4])

# ---------------- D: complexity control (REAL data, same source as Fig. S7A) ----------------
PCd=pd.read_csv(ROOT+'/results/E1_percell_sub.csv.gz')
PCd['dec']=pd.qcut(PCd.ng,10,labels=False,duplicates='drop')
g=PCd.groupby('dec').apply(lambda d: pd.Series({
    'ng':d.ng.median(),
    'mega':100*np.mean(d.code==1),
    'cyto':100*np.mean(d.code.isin([2,3])),
    'eryth':100*np.mean(d.code==4)}))
FSD=F.panel_size(0.53,1.42)
fig,ax=plt.subplots(figsize=FSD)
ax.plot(g.ng,g.cyto,'-o',color='#8e1b1b',lw=3.0,ms=10,label='cytotrophoblast',zorder=4)
ax.plot(g.ng,g.mega,'-o',color='#e67e22',lw=2.4,ms=9,label='megakaryocyte',zorder=3)
ax.plot(g.ng,g.eryth,'--s',color='#7f8c8d',lw=2.2,ms=8,label='erythrocyte (quarantined artifact)',zorder=3)
ax.set_xlabel('transcriptome complexity of the malignant cell\n(median detected genes, by decile)',fontsize=BASE*0.92)
ax.set_ylabel('% of malignant cells assigned',fontsize=BASE*0.92)
ymax=float(max(g.cyto.max(),g.mega.max(),g.eryth.max()))
ax.set_ylim(0,ymax*1.30)
xspan=float(g.ng.iloc[-1]-g.ng.iloc[0]); ax.set_xlim(float(g.ng.iloc[0])-0.10*xspan,float(g.ng.iloc[-1])+0.34*xspan)
# a legend in the empty lower-left corner: outside the axes it sat on the x-axis label, and
# labelling the curves directly put text on the neighbouring curve
ax.legend(fontsize=BASE*0.80,frameon=False,loc='lower left',borderaxespad=0.4,handlelength=1.8)
ax.annotate('rises with complexity \u2192\nnot a low-complexity artifact',
            xy=(float(g.ng.iloc[-3]),float(g.cyto.iloc[-3])),
            xytext=(float(g.ng.iloc[2]),ymax*1.27),ha='left',va='top',
            fontsize=BASE*0.82,color='#8e1b1b',arrowprops=dict(arrowstyle='->',color='#8e1b1b',lw=1.8))
# the second annotation ('both fall with complexity') is left to the caption: in-panel it
# crossed both the megakaryocyte and the erythrocyte curve
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(alpha=.22,zorder=0)
plt.savefig(FIG+'/v11_fig3D_complexity.png',dpi=200,bbox_inches='tight'); plt.close()
F.report('Fig 3D complexity',FSD,0.53)
print('wrote v11_fig3D_complexity.png (real deciles)')
print('   cyto %.1f -> %.1f | mega %.1f -> %.1f | eryth %.1f -> %.1f'%(
      g.cyto.iloc[0],g.cyto.iloc[-1],g.mega.iloc[0],g.mega.iloc[-1],g.eryth.iloc[0],g.eryth.iloc[-1]))
