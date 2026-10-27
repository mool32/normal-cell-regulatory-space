"""v11 figure assembler — compose multi-panel figures from panel PNGs with bold panel letters.
House style (Gudkov): nothing written above the figure; panels labelled A/B/C...; large fonts inside panels.
layout = list of rows; each row = list of (png, letter) or (png, letter, width_fraction).
width_fraction is the share of the page width the panel occupies (they need not sum to 1)."""
import matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path
FIG=Path('.').resolve()/'manuscript/figures'

def compose(out, layout, width=15.0, letter_size=22, gap=0.012):
    rows=[]
    for row in layout:
        items=[]
        for it in row:
            png,letter = it[0],it[1]; rel = it[2] if len(it)>2 else None
            img=mpimg.imread(FIG/png); items.append([img,letter,rel])
        auto=[i for i in items if i[2] is None]
        used=sum(i[2] for i in items if i[2] is not None)
        if auto:
            share=max(0.0,(1.0-used))/len(auto) if used<1 else 1.0/len(items)
            for i in auto: i[2]=share
        rows.append(items)
    row_h=[]
    for items in rows:
        row_h.append(max((im.shape[0]/im.shape[1])*rel*width for im,_,rel in items))
    H=sum(row_h)+gap*width*(len(rows)-1)
    fig=plt.figure(figsize=(width,H))
    y=1.0
    for items,h in zip(rows,row_h):
        y-= h/H
        x=0.0
        for im,letter,rel in items:
            w=rel
            hh=(im.shape[0]/im.shape[1])*rel*width/H
            ax=fig.add_axes([x, y+(h/H-hh)/2, w, hh]); ax.imshow(im); ax.axis('off')
            fig.text(x+0.002, y+h/H-0.002, letter, fontsize=letter_size, fontweight='bold', va='top', ha='left')
            x+=w
        y-= gap*width/H
    fig.savefig(FIG/out, dpi=200, bbox_inches='tight', facecolor='white'); plt.close()
    print('wrote', out)

if __name__=='__main__':
    import sys
    which=sys.argv[1] if len(sys.argv)>1 else 'all'
    if which in ('1','all'):
        compose('v11_fig1.png',
                [[('v11_fig1A_workflow.png','A',1.0)],
                 [('v11_fig1B_axes.png','B',0.63), ('v11_fig1D_nnfamily.png','C',0.37)]], width=13.0)
    if which in ('3','all'):
        compose('v11_fig3.png',
                [[('v11_fig3A_map.png','A',0.56), ('v11_fig3B_recurrence.png','B',0.44)],
                 [('v11_fig3C_breadthdepth.png','C',0.47), ('v11_fig3D_complexity.png','D',0.53)],
                 [('v11_fig3E_heatmap.png','E',0.60), ('v11_fig3F_boundaries.png','F',0.40)]], width=16.0)
    if which in ('4','all'):
        compose('v11_fig4.png',
                [[('v11_fig4AB_signature.png','',0.62), ('v11_fig4C_correlation.png','C',0.38)]], width=13.0)
    if which in ('5','all'):
        compose('v11_fig5.png',
                [[('v11_fig5A_paired.png','A',0.5), ('v11_fig5B_scatter.png','B',0.5)]], width=12.5)
    if which in ('6','all'):
        compose('v11_fig6.png',
                [[('v11_fig6A_occupant.png','A',0.5), ('v11_fig6B_offaxis.png','B',0.5)],
                 [('v11_fig6C_schematic.png','C',1.0)]], width=14.5)
    if which in ('2','all'):
        compose('v11_fig2.png',
                [[('v11_fig2A_normalmap.png','A',0.5), ('v11_fig2B_malignant.png','B',0.5)],
                 [('v11_fig2C_contraction.png','C',0.5), ('v11_fig2D_null.png','D',0.5)]], width=15.0)
