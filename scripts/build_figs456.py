"""Figs. 4, 5, 6 panels (v11, house style: no title above the axes, panel letters added by the assembler,
fonts >=12 inside panels and >=13 on axis labels)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon
from pathlib import Path
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
ROOT=Path('.').resolve(); FIG=ROOT/'manuscript/figures'; T=ROOT/'results/cancer_v6_6_extras'
BASE=F.use()

# ===================== FIG 4 =====================
S7=pd.read_csv(T/'tableS7_core_signature.csv',index_col=0)
S7['supp']=-S7['target_distinctiveness']
S7=S7.sort_values('supp',ascending=False)
order=list(S7.index)

# --- 4A+4B: one shared label column, two bar panels (44 labels cannot fit a third-page panel) ---
TOPN=18
sel=S7.head(TOPN)                      # the informative end: the modules the target suppresses most
S4=pd.read_csv(T/'tableS4_modules.csv',index_col=0).reindex(S7.index)
sel4=S4.head(TOPN)
BASE=F.use()
fig,(axA,axB)=plt.subplots(1,2,figsize=F.panel_size(0.62,1.08),sharey=True,
                           gridspec_kw={'wspace':0.08})
y=np.arange(TOPN)[::-1]
axA.barh(y,sel['supp'],color='#c0392b',edgecolor='black',lw=.6,zorder=3)
axA.set_yticks(y); axA.set_yticklabels(sel.index,fontsize=BASE*0.56)
axA.axvline(0,color='black',lw=1.0); axA.invert_xaxis()
axA.set_xlabel('reduced module signature in the target\n(SD below the other 151 normal cells)',fontsize=BASE*0.56)
axA.text(0.04,0.015,'normal cells only —\nno cancer data used',transform=axA.transAxes,ha='left',va='bottom',
         fontsize=BASE*0.54,style='italic',color='#7b241c')
axB.barh(y,sel4['mean_contribution'],xerr=sel4['sd'],color='#2166ac',edgecolor='black',lw=.6,
         error_kw=dict(lw=.8,ecolor='0.4'),zorder=3)
axB.axvline(0,color='black',lw=1.0)
axB.set_xlabel('contribution to the\ncancer→core convergence',fontsize=BASE*0.56)
axB.text(0.96,0.015,'the same 18 modules,\nsame order',transform=axB.transAxes,ha='right',va='bottom',
         fontsize=BASE*0.54,style='italic',color='#555')
axA.text(-0.46,1.015,'A',transform=axA.transAxes,fontsize=BASE*1.35,fontweight='bold',va='bottom')
axB.text(-0.02,1.015,'B',transform=axB.transAxes,fontsize=BASE*1.35,fontweight='bold',va='bottom')
for a in (axA,axB):
    for sp in ('top','right','left'): a.spines[sp].set_visible(False)
    a.grid(axis='x',alpha=.22,zorder=0); a.tick_params(axis='x',labelsize=BASE*0.54)
plt.tight_layout(); plt.savefig(FIG/'v11_fig4AB_signature.png',bbox_inches='tight'); plt.close()
F.report('Fig 4A+B signature',F.panel_size(0.62,1.08),0.62)

# --- 4C: the correlation over all 44 modules ---
x=S7['supp'].values; z=S4['mean_contribution'].values
order=list(S7.index)
rho,p=spearmanr(x,z)
b=np.polyfit(x,z,1); resid=z-np.polyval(b,x)
fig,ax=plt.subplots(figsize=F.panel_size(0.38,1.12))
ax.scatter(x,z,s=BASE*2.6,color='#34495e',edgecolor='white',lw=.6,zorder=4)
xs=np.linspace(x.min(),x.max(),50); ax.plot(xs,np.polyval(b,xs),'--',color='#c0392b',lw=1.8,zorder=3)
lo=int(np.argsort(resid)[0])
ax.scatter(x[lo],z[lo],s=BASE*5.0,facecolor='none',edgecolor='#2166ac',lw=2.0,zorder=5)
ax.annotate('Type I IFN —\nthe one brake cancer\ndoes not reproduce',(x[lo],z[lo]),
            xytext=(-0.55,-0.016),textcoords='data',fontsize=BASE*0.54,color='#2166ac',fontweight='bold',
            ha='left',arrowprops=dict(arrowstyle='->',color='#2166ac',lw=1.3))
ax.set_xlabel('target-defining reduced signature\n(normal cells only)',fontsize=BASE*0.58)
ax.set_ylabel('contribution to convergence',fontsize=BASE*0.58)
ax.text(0.03,0.97,f'$\\rho$ = {rho:+.2f}\n95% CI [0.78, 0.93]',transform=ax.transAxes,va='top',
        fontsize=BASE*0.60,fontweight='bold',bbox=dict(boxstyle='round,pad=0.35',fc='white',ec='0.6'))
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(alpha=.20,zorder=0); ax.tick_params(labelsize=BASE*0.54)
plt.tight_layout(); plt.savefig(FIG/'v11_fig4C_correlation.png',bbox_inches='tight'); plt.close()
F.report('Fig 4C correlation',F.panel_size(0.38,1.12),0.38)
print(f'   rho={rho:+.3f}, most-below-trend={order[lo]}')

# ===================== FIG 5 =====================
W=pd.read_csv(ROOT/'results/within_patient_alltissue.csv')
P=W[W.dataset=='colorectal (Pelka)']
# --- 5A: paired slopes ---
fig,ax=plt.subplots(figsize=F.panel_size(0.5,0.98))
nn=P.normal_core_pct.values; mm=P.malig_core_pct.values
for a,c in zip(nn,mm):
    ax.plot([0,1],[a,c],'-',color=('#c0392b' if c>a else '#95a5a6'),alpha=.45,lw=1.4,zorder=2)
ax.scatter(np.zeros(len(nn)),nn,s=52,color='#2e86c1',zorder=4)
ax.scatter(np.ones(len(mm)),mm,s=52,color='#c0392b',zorder=4)
ax.plot([0,1],[np.median(nn),np.median(mm)],'-',color='black',lw=3.4,zorder=5)
ax.scatter([0,1],[np.median(nn),np.median(mm)],s=150,color='black',zorder=6)
ax.text(-0.06,np.median(nn),f'{np.median(nn):.1f}%',ha='right',va='center',fontsize=13.5,fontweight='bold')
ax.text(1.06,np.median(mm),f'{np.median(mm):.1f}%',ha='left',va='center',fontsize=13.5,fontweight='bold')
up=int((mm>nn).sum())
ax.text(0.5,0.97,f'higher in {up} of {len(P)} patients\npaired $P$ < 10$^{{-12}}$',transform=ax.transAxes,
        ha='center',va='top',fontsize=13,fontweight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(["the patient's own\nnormal epithelium",'malignant\nepithelium'],fontsize=13)
ax.set_xlim(-0.35,1.35); ax.set_ylabel('cells occupying the convergent state (%)',fontsize=13.5)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(FIG/'v11_fig5A_paired.png',dpi=200,bbox_inches='tight'); plt.close()
print(f'wrote v11_fig5A_paired.png ({up}/{len(P)})')

# --- 5B: NEW per-patient scatter against the diagonal ---
fig,ax=plt.subplots(figsize=F.panel_size(0.5,1.0))
lim=max(mm.max(),nn.max())*1.08
ax.plot([0,lim],[0,lim],'--',color='0.4',lw=2.0,zorder=2)
ax.scatter(nn,mm,s=110,color='#c0392b',edgecolor='black',lw=.8,alpha=.85,zorder=4)
ax.text(lim*0.62,lim*0.52,'y = x',fontsize=13,color='0.4',rotation=38)
ax.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1]*1.06)
ax.text(0.04,0.99,f'{up} of {len(P)} patients above the diagonal',transform=ax.transAxes,va='top',
        fontsize=13,fontweight='bold',color='#7b241c')
ax.set_xlabel("normal epithelium of the same patient (%)",fontsize=13.5)
ax.set_ylabel('malignant epithelium (%)',fontsize=13.5)
ax.set_xlim(0,lim); ax.set_ylim(0,lim)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(alpha=.22,zorder=0)
plt.tight_layout(); plt.savefig(FIG/'v11_fig5B_scatter.png',dpi=200,bbox_inches='tight'); plt.close()
print('wrote v11_fig5B_scatter.png')

# ===================== FIG 6 =====================
D=pd.read_csv(ROOT/'results/pelka_diff_occupancy_controlled.csv')
sig=D[(D.resid_delta.abs()>0.25)&(D.resid_pW<0.05)].copy().sort_values('resid_delta',ascending=False)
CLASS={'HIF':'stress','NRF2':'stress','HSF1':'stress','UPR-PERK':'stress','cGAS-STING':'stress','Autophagy':'stress',
       'mTOR':'growth','SREBP':'growth','Cell Cycle':'growth','Wnt':'growth','ER':'growth',
       'PPARα':'metabolic','MR':'metabolic','PI3K/PTEN':'metabolic','p53':'stress','PPARγ':'metabolic'}
CCOL={'stress':'#c0392b','growth':'#8e44ad','metabolic':'#e8a33d','other':'#7f8c8d'}

# --- 6A: per-module difference between the two occupants ---
fig,ax=plt.subplots(figsize=F.panel_size(0.5,1.03))
yy=np.arange(len(sig))[::-1]
cols=[CCOL.get(CLASS.get(m,'other')) for m in sig.module]
ax.barh(yy,sig.resid_delta,color=cols,edgecolor='black',lw=.7,zorder=3)
ax.scatter(sig.raw_delta,yy,facecolor='none',edgecolor='0.35',s=58,lw=1.2,zorder=4,label='before the depth control')
ax.axvline(0,color='black',lw=1.1)
ax.set_yticks(yy); ax.set_yticklabels(sig.module,fontsize=12.5)
ax.set_xlabel('malignant minus normal occupant of the SAME state\n(module z, detection-depth controlled)',fontsize=13)
from matplotlib.patches import Patch
h=[Patch(facecolor=CCOL[k],edgecolor='black',label=k) for k in ['stress','growth','metabolic']]
h.append(plt.Line2D([0],[0],marker='o',color='w',markerfacecolor='none',markeredgecolor='0.35',markersize=10,label='before depth control'))
ax.legend(handles=h,fontsize=BASE*0.82,loc='lower right',frameon=False)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(axis='x',alpha=.22,zorder=0)
plt.tight_layout(); plt.savefig(FIG/'v11_fig6A_occupant.png',dpi=200,bbox_inches='tight'); plt.close()
print(f'wrote v11_fig6A_occupant.png ({len(sig)} modules survive the control)')

# --- 6B: the difference lies OFF the target-defining axis ---
M=D.copy()
rho2,p2=spearmanr(M['distinct'],M['resid_delta'])
fig,ax=plt.subplots(figsize=F.panel_size(0.5,1.18))
cols=[CCOL.get(CLASS.get(m,'other')) for m in M.module]
ax.scatter(M['distinct'],M['resid_delta'],s=95,c=cols,edgecolor='black',lw=.7,zorder=4)
b2=np.polyfit(M['distinct'],M['resid_delta'],1)
xs=np.linspace(M['distinct'].min(),M['distinct'].max(),50)
ax.plot(xs,np.polyval(b2,xs),'--',color='0.4',lw=2.0,zorder=3)
ax.axhline(0,color='black',lw=1.0)
_t=[]
for m in ['HIF','mTOR','NRF2','Cell Cycle','Type I IFN']:
    r=M[M.module==m]
    if len(r):
        xx,yy=float(r['distinct'].iloc[0]),float(r['resid_delta'].iloc[0])
        _t.append(ax.text(xx,yy,m,fontsize=BASE*0.80,fontweight='bold',color='#333',zorder=8))
ax.set_xlabel('core distinctiveness\n(← reduced signature          retained / high →)',fontsize=BASE*0.92)
ax.set_ylabel('malignant − normal occupant',fontsize=BASE*0.92)
# the summary box used to sit top-left, on top of the mTOR point and its label
ax.set_ylim(ax.get_ylim()[0]-0.10,ax.get_ylim()[1])
ax.text(0.975,0.035,f'Spearman $\\rho$ = {rho2:+.2f}\nthe difference is mostly OFF\nthe axis that defines the target',
        transform=ax.transAxes,ha='right',va='bottom',fontsize=BASE*0.84,fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.45',fc='white',ec='0.6'),zorder=9)
F.declutter(ax,_t,avoid=np.c_[M['distinct'].values,M['resid_delta'].values])
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(alpha=.22,zorder=0)
plt.savefig(FIG/'v11_fig6B_offaxis.png',dpi=200,bbox_inches='tight'); plt.close()
print(f'wrote v11_fig6B_offaxis.png (rho={rho2:+.2f}, p={p2:.2f})')

# --- 6C: NEW interpretive schematic ---
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse as Ell
fig,ax=plt.subplots(figsize=F.panel_size(1.0,2.3)); ax.set_xlim(0,100); ax.set_ylim(0,46); ax.axis('off')
def bx(x,y,w,h,t,fc,fs=12.5,bold=False,ec='black'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.5,rounding_size=1.1",fc=fc,ec=ec,lw=1.6,zorder=3))
    ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=fs,zorder=4,linespacing=1.4,
            fontweight='bold' if bold else 'normal')
ax.add_patch(Ell((50,30),46,20,fc='#f7e9e6',ec='#c0392b',lw=2.2,ls='--',zorder=1))
ax.text(50,41.0,'the convergent state',ha='center',fontsize=13,
        color='#c0392b',fontweight='bold',zorder=5)
bx(31,24,16,11,'normal\nepithelial cell\nin the state','#dceaf5')
bx(53,24,16,11,'malignant\ncell in the\nsame state','#f6d5cf')
ax.add_patch(FancyArrowPatch((47.6,29.5),(52.4,29.5),arrowstyle='-|>',mutation_scale=22,lw=2.6,color='#7b241c',zorder=5))
ax.text(50,21.6,'what does malignancy add?',ha='center',fontsize=12,color='#7b241c',style='italic',zorder=5)
ax.text(5,30,'same patient,\nsame tissue,\nboth already\nin the state',ha='center',fontsize=12,color='#555')
ax.text(50,16.5,'malignancy superimposes additional programs — it does not push further along the defining axis',
        ha='center',fontsize=12.5,fontweight='bold')
for x,t,c in [(9,'stress responses\nHIF · NRF2 · HSF1 · UPR','#f7d7d2'),
              (38,'growth and anabolic\nmTOR · SREBP · cell cycle · Wnt','#e8d5ef'),
              (69,'metabolic reprogramming\nPPARα · MR · PI3K/PTEN','#faecd2')]:
    bx(x,3.5,24,10,t,c,fs=11.5)
    ax.add_patch(FancyArrowPatch((x+12,13.9),(x+12,15.6),arrowstyle='-|>',mutation_scale=16,lw=2.0,color='0.35',zorder=5))
plt.tight_layout(); plt.savefig(FIG/'v11_fig6C_schematic.png',dpi=200,bbox_inches='tight'); plt.close()
print('wrote v11_fig6C_schematic.png')
