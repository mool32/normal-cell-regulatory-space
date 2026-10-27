"""Fig. 1A — the analytical workflow plus the anatomy of one module.
Sized through scripts/figstyle.py so the printed in-figure text is 8 pt."""
import sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figstyle as F
FIG=Path('.').resolve()/'manuscript/figures'
SHARE, ASPECT = 1.0, 3.25
B=F.use(); FS=F.panel_size(SHARE,ASPECT)
fig,ax=plt.subplots(figsize=FS); ax.set_xlim(0,100); ax.set_ylim(-3,41); ax.axis('off')
s=lambda m: B*m/2.0                       # font helper, relative to the base

def box(x,y,w,h,txt,fc,fs=None,ec='black'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.45,rounding_size=1.0",fc=fc,ec=ec,lw=1.5,zorder=3))
    ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs or s(1.0),zorder=4,linespacing=1.4)
def arrow(x1,y1,x2,y2,color='0.25',lw=2.0,style=None):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=16,lw=lw,color=color,
                                 zorder=5,connectionstyle=style) if style else
                 FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=16,lw=lw,color=color,zorder=5))

Y,H=27.5,10.0
steps=[(1.0,17.0,"154 normal\nhuman cell types\n(single-cell RNA)","#dceaf5"),
       (22.0,17.0,"44 regulatory\nmodules scored\nper cell type","#cfe3d4"),
       (43.0,17.0,"154 × 44\nmodule-activity\nmatrix","#f2e6cc"),
       (64.0,17.0,"principal\ncomponent\ndecomposition","#f2e6cc"),
       (85.0,14.0,"9-D reference\nspace (fixed)","#f6d5cf")]
for x,w,t,c in steps: box(x,Y,w,H,t,c)
for i in range(len(steps)-1):
    arrow(steps[i][0]+steps[i][1]+0.5,Y+H/2,steps[i+1][0]-0.5,Y+H/2)
ax.text(50,23.6,"built from normal cells only — no cancer data enter at any step",
        ha='center',fontsize=s(0.98),style='italic',color='#7b241c')
ax.text(1.0,19.2,"each module is a closed regulatory circuit:",fontsize=s(1.0),fontweight='bold')
Y2,H2=7.0,8.6
mods=[(1.0,18.0,"sensor\n(receptor or\nintracellular)","#dceaf5"),
      (23.5,18.0,"transcription-\nfactor output","#cfe3d4"),
      (46.0,18.0,"target genes\n(the readout)","#f2e6cc"),
      (68.5,18.0,"negative-feedback\nelement","#e8d5e8")]
for x,w,t,c in mods: box(x,Y2,w,H2,t,c,fs=s(0.95))
for i in range(len(mods)-1):
    arrow(mods[i][0]+mods[i][1]+0.5,Y2+H2/2,mods[i+1][0]-0.5,Y2+H2/2)
arrow(77.5,Y2,10.0,Y2,color='#7d3c98',lw=2.2,style="arc3,rad=-0.28")
ax.text(44,-0.6,"the feedback element shuts the circuit off",ha='center',fontsize=s(0.98),
        color='#7d3c98',style='italic')
ax.text(88.0,Y2+H2/2,"e.g. GR:\nNR3C1 → GR →\nFKBP5",ha='left',va='center',fontsize=s(0.88),color='0.3')
plt.tight_layout(); plt.savefig(FIG/'v11_fig1A_workflow.png',bbox_inches='tight'); plt.close()
F.report('Fig 1A workflow',FS,SHARE)
