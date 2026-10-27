"""Fig. 1B (v11) — trimmed to the three strongest axes, house style (large fonts, no suptitle).
Original: each principal axis tells a biological story.
One row per stable PC (PC1-PC6). All 154 normal cell types placed on the 1-D number line of that PC
(vertical jitter only), coloured by biological family; the most extreme cells at each end are labelled;
the modules that load each end are annotated. Shows the axes separate biologically coherent groups —
the coordinate system is meaningful, not arbitrary. Run from repo root.
"""
import sys, json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix()); import figstyle as F
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
sys.path.insert(0, __import__('pathlib').Path(__file__).resolve().parents[1].as_posix()+'/tool/perceptome2')
from perceptome.perceptivity.reference import load_hpa_perceptivity
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
e=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')); mods=e['module_order']
load=np.array([[e['loadings'][f'PC{i+1}'][m] for i in range(9)] for m in mods])
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
mu=Rm.mean(0).values; sd=Rm.std(0).values+1e-12
PC=(Rm.values-mu)/sd@load; cells=Rm.index.tolist()
POST={'erythrocytes','platelets','reticulocytes'}   # post-cellular bodies: coordinate is a low-complexity artifact, excluded
_keep=[i for i in range(len(cells)) if cells[i] not in POST]
PC=PC[_keep]; cells=[cells[i] for i in _keep]; N=len(cells)
fam=json.load(open(ROOT+'/results/fam.json'))
FAMS=['germline / embryonic','stem / progenitor','erythroid / megakaryocyte','blood / immune',
 'epithelial — secretory/glandular','epithelial — barrier/squamous','respiratory','hepatic','renal',
 'endocrine','neural','stromal / mesenchymal','endothelial','muscle','melanocyte']
palette=dict(zip(FAMS,[cm.tab20(i/20) for i in range(len(FAMS))]))
col=[palette[fam[c]] for c in cells]
varexp=PC.var(0)/PC.var(0).sum()*100
def axmods(i):
    L={m:e['loadings'][f'PC{i+1}'][m] for m in mods}; s=sorted(L.items(),key=lambda x:x[1])
    return ', '.join(m for m,_ in s[-3:][::-1]), ', '.join(m for m,_ in s[:3])
def short(c):
    return c.replace(' cells','').replace('undifferentiated ','undiff. ').replace('differentiating ','diff. ')[:26]

import re
enr=json.load(open(ROOT+'/results/E5_pc_enrichment.json'))
q6=json.load(open(ROOT+'/results/E6_axis_quality.json'))
# The share of TOTAL variance, from the shipped eigenspace -- NOT E6_axis_quality.json's `var_pct`,
# which is the share of the retained nine-dimensional variance (31%, 21%, 11% ...). Quoting that as
# "% of variance" is the error block 2 withdrew from the main text; it survived here, inside the
# raster, where no text search of the PDF can reach it.
_VAR = [100 * v for v in json.load(open(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json')))['var_explained']]
assert abs(_VAR[0] - 29.17) < 0.05 and abs(_VAR[1] - 13.72) < 0.05, _VAR[:2]

SHORT={  # legible stems for long enrichment terms (no mid-word ellipsis)
 'diseases of signal transduction by growth factor receptors and second messengers':'signal-transduction disease genes',
 'runx1 regulates transcription of genes involved in interleukin signaling':'RUNX1 transcription',
 'pou5f1 (oct4), sox2, nanog activate genes related to proliferation':'OCT4/SOX2/NANOG pluripotency',
 'srp-dependent cotranslational protein targeting to membrane':'SRP cotranslational targeting'}
def clean_term(t):
    t=re.sub(r'\s*R-HSA-\d+$','',t); t=re.sub(r'\s*\(GO:\d+\)$','',t); t=re.sub(r'\s*Homo sapiens.*$','',t)
    if t.lower() in SHORT: return SHORT[t.lower()]
    return (t[:38]+'…') if len(t)>39 else t
GENERIC={'disease','metabolism','immune system','signal transduction','gene expression (transcription)',
 'generic transcription pathway','metabolism of proteins','metabolism of rna','developmental biology',
 'vesicle-mediated transport','hemostasis','cell cycle','metabolism of lipids','post-translational protein modification',
 'innate immune system','adaptive immune system','cytokine signaling in immune system'}
def topterm(pc,pole):
    L=enr.get(f'PC{pc+1}',{}).get(pole,[])
    for it in L:
        t=clean_term(it['term'])
        if t.lower().rstrip('…') not in GENERIC: return t
    return clean_term(L[0]['term']) if L else '?'
qdf=pd.read_csv(ROOT+'/results/E6_quality_scores.csv').set_index('cell').reindex(cells)
GX=3.6; NP=6
NP=2
SHARE,ASPECT=0.63,1.45
B=F.use(); FS=F.panel_size(SHARE,ASPECT)
fig,axes=plt.subplots(NP,1,figsize=FS)
for i in range(NP):
    ax=axes[i]; qd=q6[f'PC{i+1}']
    x=np.clip(qdf[f'PC{i+1}_pos'].values,-GX,GX); y=qdf[f'PC{i+1}_qual'].values
    ax.scatter(x,y,c=col,s=44,edgecolor='white',lw=.3,zorder=3)
    b1,b0=np.polyfit(x,y,1); xs=np.array([x.min(),x.max()])        # trend line = the monotonic "decline/rise"
    ax.plot(xs,b0+b1*xs,color='#111',lw=2.4,zorder=6,solid_capstyle='round')
    ax.axvline(0,color='#ddd',lw=.6,zorder=1)
    order=np.argsort(qdf[f'PC{i+1}_pos'].values)
    # Fixed offsets collided: the three spermatogenic labels of PC1 and the three secretory labels of
    # PC2 sit on nearly the same point, so a constant (-3, 7) stacked them on top of one another. Same
    # treatment as Fig. S3: ax.text artists (declutter cannot move an offset Annotation), a top band
    # of headroom to move them into, and the stat box handed over as static obstacles.
    _ta=[ax.text(x[k],y[k],short(cells[k]),fontsize=B*0.52,color='#555',zorder=7,ha='center',va='bottom')
         for k in list(order[:3])+list(order[-3:])]
    ax.set_xlim(-GX-0.5,GX+0.5)
    _y0,_y1=y.min(),y.max(); _r=max(_y1-_y0,1e-9)
    _up = qd['spearman_r'] <= 0                                  # the box sits low when quality rises
    ax.set_ylim(_y0-(0.10 if _up else 0.42)*_r, _y1+(0.42 if _up else 0.10)*_r)
    _bx=np.linspace(-GX-0.5,GX+0.5,9)[4:]
    _by=(np.linspace(_y1+0.18*_r,_y1+0.42*_r,3) if _up else np.linspace(_y0-0.42*_r,_y0-0.18*_r,3))
    F.declutter(ax,_ta,avoid=np.vstack([np.c_[x,y],np.array([[a,b] for a in _bx for b in _by])]),
                color='#999',lw=0.5,expand=(1.35,1.6),force_text=(0.8,1.0),force_static=(0.4,0.5))
    ax.set_ylabel(f'{qd["quality"]}\n(signature score, z)',fontsize=B*0.62)
    # A5: no titles above the axes — the axis descriptor belongs to the axis itself
    ax.set_xlabel(f'PC{i+1}  ({_VAR[i]:.0f}% of variance):   ◄ {topterm(i,"neg")}'
                  f'      {topterm(i,"pos")} ►', fontsize=B*0.62, weight='bold')
    arrow='falls' if qd['spearman_r']<0 else 'rises'
    _rise = qd['spearman_r'] > 0
    ax.text(0.985, (0.06 if _rise else 0.93),
            f"families η² = {qd['eta2_family']:.2f}  (p = {qd['p_anova']:.0e})\nquality {arrow} along axis:  r = {qd['spearman_r']:+.2f}",
            transform=ax.transAxes,ha='right',va=('bottom' if _rise else 'top'),fontsize=B*0.60,color='#222',
            bbox=dict(boxstyle='round,pad=0.3',fc='#f6f6f4',ec='#bbb',lw=.6,alpha=.93))
    for s in ['top','right']: ax.spines[s].set_visible(False)
# (each axis now carries its own descriptor; the shared label would duplicate it)
handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=palette[f],markersize=8,label=f) for f in FAMS]
fig.legend(handles=handles,fontsize=B*0.56,loc='lower center',bbox_to_anchor=(0.5,-0.115),ncol=5,framealpha=.95,title='biological family',title_fontsize=B*0.60)
plt.tight_layout(rect=(0,0.085,1,1))
plt.savefig(FIG+'/v11_fig1B_axes.png',bbox_inches='tight'); plt.close()
F.report('Fig 1B axes',FS,SHARE)
print('saved v11_fig1B_axes.png  | printed (share of TOTAL variance):',[f'{v:.0f}%' for v in _VAR[:6]],
      '| within the retained 9-D (NOT what the axes print):',[f'{v:.0f}%' for v in varexp[:6]])
for i in range(NP): print(f'PC{i+1}  neg[{topterm(i,"neg")}]  pos[{topterm(i,"pos")}]')
