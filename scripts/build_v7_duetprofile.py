"""Gudkov's request: draw the duet itself against the cloud of all 154 normal cells, per module.
For each top target-defining module, the distribution of all 154 cells' z-scored module signature,
with megakaryocyte + cytotrophoblast marked — showing how far the duet falls out of the normal cloud."""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0,__import__('pathlib').Path(__file__).resolve().parent.as_posix())
import figstyle as F
BASE=F.use()   # authored at the printed width: 1.00 x 6.5 in
sys.path.insert(0,'tool/perceptome2'); from perceptome.perceptivity.reference import load_hpa_perceptivity
ROOT=__import__('pathlib').Path(__file__).resolve().parents[1].as_posix(); FIG=ROOT+'/manuscript/figures'
mods=json.load(open(ROOT+'/tool/perceptome2/perceptome/eigenspace/data/reference_v03.json'))['module_order']
Rm=load_hpa_perceptivity(mode='single_cell_scaled')['R'][mods]
z=(Rm-Rm.mean(0))/(Rm.std(0)+1e-9)
d=pd.read_csv('results/cancer_v6_6_extras/tableS7_core_signature.csv',index_col=0)
top=d.sort_values('target_distinctiveness').head(12).index.tolist()  # most core-suppressed (target-defining)
def ordq(n):
    n=int(round(n)); suf='th' if (n%100 in (11,12,13)) else {1:'st',2:'nd',3:'rd'}.get(n%10,'th'); return f'{n}{suf}'
rng=np.random.RandomState(0)
fig,ax=plt.subplots(figsize=F.panel_size(1.0, 1.194))
for i,m in enumerate(top[::-1]):
    y=len(top)-1-i
    vals=z[m].values; jit=rng.uniform(-0.28,0.28,len(vals))
    ax.scatter(vals, y+jit, s=12, color='#c9cdd2', edgecolor='none', zorder=2)
    mk=z.loc['megakaryocytes',m]; ct=z.loc['cytotrophoblasts',m]
    ax.scatter(mk, y+0.0, s=85, color='#c0392b', marker='D', edgecolor='white', lw=.6, zorder=5)
    ax.scatter(ct, y+0.0, s=85, color='#7b241c', marker='o', edgecolor='white', lw=.6, zorder=5)
    # percentile of the duet mean within the 154
    dm=(mk+ct)/2; pct=100*(vals<dm).mean()
    ax.text(z.values.min()-0.3, y, m, ha='right', va='center', fontsize=BASE*0.943)
    ax.text(min(mk,ct)-0.3, y, f'{ordq(pct)} pct', ha='right', va='center', fontsize=BASE*0.825, color='#7b241c', weight='bold')
ax.axvline(0, color='#bbb', lw=.8, zorder=1)
ax.set_yticks([]); ax.set_ylim(-0.7,len(top)-0.3)
ax.set_xlim(z.values.min()-2.6, z.values.max()+1.4)
ax.set_xlabel('module signature (z-score across the 154 normal cell types)', fontsize=BASE*0.976)
# the long two-line ax title duplicated the caption and set the tight bbox (and so the printed type size)
plt.tight_layout(); plt.savefig(FIG+'/v7_fig_duetprofile.png',dpi=190,bbox_inches='tight'); plt.close()
print('saved v7_fig_duetprofile.png  | modules:',top)
