"""v11 figure style — guarantees the printed point size of in-figure text.

The problem this solves: a panel drawn at figsize W_src inches with font F pt, placed at a share S of a
composite that is finally printed P inches wide, appears at  F * S * P / W_src  points. Authoring panels at
13 pt in 7-15 inch figures and printing them at 4.5-6 inches gave 2.4-4.6 pt — below the 6 pt floor.

The rule here: set the panel's width to  S * P * K  and its base font to  FINAL_PT * K.  The K cancels, so
the printed size is exactly FINAL_PT whatever K is; K just keeps the authoring canvas a comfortable size.
"""
import matplotlib.pyplot as plt

PRINT_W  = 6.5    # inches — a full-width figure in the journal text block
FINAL_PT = 8.0    # the point size in-figure body text should have on paper
K        = 2.0    # authoring magnification

def panel_size(share, aspect):
    """figsize for a panel occupying `share` of the printed width, at width:height = `aspect`."""
    w = share * PRINT_W * K
    return (w, w / aspect)

def use(scale=1.0):
    """Apply the house rcParams. `scale` nudges one panel up or down relative to FINAL_PT."""
    b = FINAL_PT * K * scale
    plt.rcParams.update({
        'font.size':        b,
        'axes.labelsize':   b * 1.05,
        'axes.titlesize':   b * 1.05,
        'xtick.labelsize':  b * 0.92,
        'ytick.labelsize':  b * 0.92,
        'legend.fontsize':  b * 0.92,
        'axes.linewidth':   max(0.8, b * 0.07),
        'xtick.major.width': max(0.8, b * 0.07),
        'ytick.major.width': max(0.8, b * 0.07),
        'lines.linewidth':  max(1.2, b * 0.11),
        'lines.markersize': max(4.0, b * 0.42),
        'savefig.dpi':      200,
        'figure.dpi':       200,
    })
    return b

def report(name, figsize, share, printed_w=PRINT_W):
    pt = plt.rcParams['font.size'] * share * printed_w / figsize[0]
    flag = 'OK' if pt >= 6.5 else 'TOO SMALL'
    print(f"   {name:34s} base {plt.rcParams['font.size']:4.1f}pt -> printed {pt:4.1f}pt  [{flag}]")
    return pt


def declutter(ax, texts, avoid=None, color='0.35', lw=0.7, **kw):
    """Push a set of already-placed Text artists apart until none of their rendered boxes overlap,
    drawing a leader line from each back to where it started. `avoid` is an (n, 2) array of data
    points the labels should also keep clear of (e.g. the scatter they annotate).

    Pass plain ax.text artists, not ax.annotate(..., textcoords='offset points') ones: an Annotation's
    position is its offset, not a data coordinate, so they all collapse onto one spot.

    Hand-tuned offsets do not survive a change of figure size or of the underlying embedding; this does.
    """
    import numpy as _np
    from adjustText import adjust_text
    # Determinism. adjustText uses the global RNG to break ties, and -- the part that actually bit --
    # defaults to a one-SECOND wall-clock budget (`time_lim`), so the iteration count, and with it the
    # final label positions, depend on how busy the machine is. Three figures came back with ~1% of
    # their pixels moved on a re-run with no input change. Seeding the RNG and trading the time budget
    # for a fixed `iter_lim` makes the same data give the same bytes. The marks never moved either way.
    _np.random.seed(0)
    opts = dict(ax=ax, expand=(1.12, 1.25), force_text=(0.45, 0.6), force_static=(0.22, 0.3),
                iter_lim=3000,
                arrowprops=dict(arrowstyle='-', color=color, lw=lw, shrinkA=1, shrinkB=2))
    opts.update(kw)
    if avoid is not None:
        a = _np.asarray(avoid, dtype=float)
        opts['x'], opts['y'] = a[:, 0], a[:, 1]
    return adjust_text(texts, **opts)
