#!/usr/bin/env python3
"""Compose the consolidated SI figures: 30 numbered items (S1-S29 plus S9b) -> 16.

Three sources are dropped outright (DROPPED below) and the remaining 27 are grouped into 16 numbered
figures. This script is the single place the grouping lives; `scripts/build_v10_pdf.py` injects the
result and MUST use the same numbers.

Two constraints dictate the method.

1. Printed font. The sources are authored by `scripts/figstyle.py`, where a panel authored at
   `panel_size(share, ...)` prints at 8 pt *when its printed width is `share` x 6.5 in*. So a
   half-width panel (1300 px at 200 dpi) must occupy half of the composite's width, and a full-width
   panel (2600 px) must occupy all of it. Narrow panels are therefore PADDED with white to the
   common width, never rescaled: rescaling is what changes the printed font size.

2. Page height. A composite printed 6.5 in wide is `h * 6.5 / w` inches tall. The usable text block
   is about 8.4 in, and the caption shares the float with the image, so anything over ~7.3 in is
   injected at a reduced `width=` in the builder; `recommend_width()` prints the value to use.

Panel labels. Several sources already carry their own (A), (B), (C) sub-panel letters, baked into the
raster. Adding outer A/B/C letters would put a literal "A" inside the block labelled "B". The blocks
of a composite therefore carry ROMAN part labels (I, II, III) drawn here, and a sub-panel is cited as
"I-A". Sources that are injected on their own keep whatever letters their own builder drew.

Writes manuscript/figures/v12_figS*.png for the merged groups only; a group with a single source is
injected straight from that source (PASSTHROUGH), so no pointless copy is made.
"""
import os
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
FIG = os.path.join(ROOT, 'manuscript', 'figures')
DPI = 200
PRINT_W_IN = 6.5
USABLE_H_IN = 7.3          # image height above which the builder should reduce `width=`
FONT_FLOOR_PT = 6.5        # scripts/figstyle.py: printed font must not go below this
BASE_PT = 8.0              # what the sources are authored to print at
GAP_PX = 18
PAD_PX = 10
ROMAN = ['I', 'II', 'III', 'IV']

# ---------------------------------------------------------------- the plan
# new number -> (output stem, rows, one-line purpose)
# `rows` is a list of rows; a row is a list of source PNGs placed side by side. One source per row
# = a vertical stack. Part labels are assigned in reading order across all rows.
GROUPS = {
    'S1':  ('v12_figS1_axis_stability',
            [['v7_figS1_bootstrap.png'], ['v7_fig_eta2.png']],
            'the nine axes: how well each is determined, and how biological each is'),
    'S2':  ('v12_figS2_redundancy', [['v7_fig_E2_redundancy.png']],
            'the 44 modules are not restating one another'),
    'S3':  ('v12_figS3_axes_biology', [['v7_fig_axes.png']],
            'each axis against an external, module-free signature'),
    'S4':  ('v12_figS4_geometry',
            [['v7_fig1A_clouds.png', 'v7_fig_map.png'], ['v7_fig1C_heatmap_z.png']],
            'where the cancers land, shown three ways'),
    'S5':  ('v12_figS5_complexity',
            [['v7_figS3_qctrap.png'], ['v7_fig_E1_complexity.png']],
            'the transcriptome-complexity trap and the control for it'),
    'S6':  ('v12_figS6_ranking',
            [['v7_figS2_2axis.png'], ['v7_fig1B_recurrence.png']],
            'ranking the candidate targets two ways'),
    'S7':  ('v12_figS7_catalog',
            [['v7_fig_E2c_target.png'], ['v7_figS9b_module_robustness.png'],
             ['v7_fig_E4_genelevel.png']],
            'how much of the result depends on this particular module catalog'),
    'S8':  ('v12_figS8_target', [['v7_fig_duetprofile.png']],
            'what distinguishes the target among the 154 normal cells'),
    'S9':  ('v12_figS9_reprogramming', [['v7_figS28_reprogramming.png']],
            'the named target in an independent system with no cancer in it'),
    'S10': ('v12_figS10_speccurve', [['v7_figS19_speccurve.png']],
            'specification curve over 36 analytic crossings'),
    'S11': ('v12_figS11_withinpat',
            [['v7_figS21_gastric.png'], ['v7_figS22_lung.png'], ['v7_figS23_shapeonly.png']],
            'the within-patient increase replicates and survives magnitude removal'),
    'S12': ('v12_figS12_fibrosis', [['v7_figS24_fibrosis.png']],
            'the benign-remodelling negative control'),
    'S13': ('v12_figS13_progression',
            [['v7_figS20_preinvasive.png'], ['v7_figS27_preinvasive2.png'],
             ['v7_figS25_metastasis.png']],
            'before invasion and after dissemination'),
    'S14': ('v12_figS14_dynamics', [['v7_figS26_dynamics.png']],
            'the coordinate under targeted therapy, with its complexity control'),
    'S15': ('v12_figS15_detuning', [['v7_fig_detuning.png']],
            'what malignancy adds within the shared state'),
}

# Dropped outright, with the reason recorded:
DROPPED = {
    'v7_figS4_preQC_erythrocyte.png': 'one sentence of content (erythrocyte ranks 3rd before '
                                      'quarantine); kept as the pre-quarantine occupancy column of '
                                      'Table S2 and a clause in Note S-K',
    'v7_fig_landscape.png':           'the same geometry as the normal-cell map (now Fig. S4-II), '
                                      'less legibly',
    'v7_fig_commitment.png':          'interpretive rather than a control; its claim is made in the '
                                      'text, and its per-module numbers are Table S7',
}


def load_font(px):
    for p in ('/System/Library/Fonts/Supplemental/Arial Bold.ttf',
              '/System/Library/Fonts/Helvetica.ttc'):
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, px)
            except Exception:
                pass
    return ImageFont.load_default()


def trim_white_rows(im, keep=10):
    """Drop fully white rows from the top and bottom of a panel.

    Vertical only, deliberately. Trimming the SIDES would narrow the panel, and a narrower panel
    printed to the same 6.5 in prints its font LARGER — the one thing the whole layout is built to
    hold fixed. Several sources carry a band of white where a descriptive title used to sit (removed
    under house rule A5), which is pure wasted page height between stacked blocks.
    """
    g = im.convert('L')
    w, h = g.size
    px = g.load()
    def blank(y):
        return all(px[x, y] > 250 for x in range(0, w, 7))   # sampled; exact scan is 10x slower
    top = 0
    while top < h - 1 and blank(top):
        top += 1
    bot = h - 1
    while bot > top and blank(bot):
        bot -= 1
    top, bot = max(0, top - keep), min(h - 1, bot + keep)
    return im if (top == 0 and bot == h - 1) else im.crop((0, top, w, bot + 1))


def edge_ink(path, band=3):
    """Warn when a panel has ink in its outermost pixels: a label that ran off its own canvas.

    Two figures shipped for months with axis labels clipped mid-word this way (the module-robustness
    x-labels, which also overlapped each other, and the fibrosis y-label cut to \"(pt\"), because
    matplotlib draws nothing outside the canvas and nobody looked at the margin.
    """
    g = Image.open(path).convert('L')
    w, h = g.size
    px = g.load()
    sides = []
    if any(px[x, y] < 240 for x in range(band) for y in range(0, h, 3)):          sides.append('left')
    if any(px[w - 1 - x, y] < 240 for x in range(band) for y in range(0, h, 3)):  sides.append('right')
    if any(px[x, y] < 240 for y in range(band) for x in range(0, w, 3)):          sides.append('top')
    if any(px[x, h - 1 - y] < 240 for y in range(band) for x in range(0, w, 3)):  sides.append('bottom')
    return sides


def recommend_width(h_in):
    """The `width=` to inject at, so a tall composite still fits the float with its caption."""
    if h_in <= USABLE_H_IN:
        return 100
    w = int(100 * USABLE_H_IN / h_in)          # scale the image down to the usable height
    floor = int(100 * FONT_FLOOR_PT / BASE_PT)  # ... but never past the printed-font floor
    return max(w, floor)


def compose(stem, rows, label_parts=True):
    """Lay `rows` out top to bottom; within a row, panels sit side by side at native pixel width."""
    loaded = []
    for row in rows:
        ims = []
        for s in row:
            p = os.path.join(FIG, s)
            if not os.path.exists(p):
                print(f'    MISSING {s}')
                return None
            ims.append(trim_white_rows(Image.open(p).convert('RGB')))
        loaded.append(ims)

    # common width = the widest row; narrow rows are padded, never rescaled
    W = max(sum(i.width for i in ims) + GAP_PX * (len(ims) - 1) for ims in loaded)
    letter_px = round(W * 0.028)
    font = load_font(letter_px)
    lead = round(letter_px * 1.55) if label_parts else 0

    H = sum(max(i.height for i in ims) + lead for ims in loaded) \
        + GAP_PX * (len(loaded) - 1) + 2 * PAD_PX
    out = Image.new('RGB', (W + 2 * PAD_PX, H), 'white')
    d = ImageDraw.Draw(out)
    y, n = PAD_PX, 0
    for ims in loaded:
        rw = sum(i.width for i in ims) + GAP_PX * (len(ims) - 1)
        x0 = PAD_PX + (W - rw) // 2                 # centre a narrow row in the white it is padded to
        if lead:
            x = x0
            for i in ims:
                d.text((x, y), ROMAN[n], fill='black', font=font)
                n += 1
                x += i.width + GAP_PX
            y += lead
        x = x0
        for i in ims:
            out.paste(i, (x, y))
            x += i.width + GAP_PX
        y += max(i.height for i in ims) + GAP_PX
    path = os.path.join(FIG, stem + '.png')
    out.save(path, dpi=(DPI, DPI))
    return path, out.size


def main():
    n_src = sum(len(r) for _, rows, _ in GROUPS.values() for r in rows)
    clipped = []
    for num, (_, rows, _) in GROUPS.items():
        for s in (s for r in rows for s in r):
            p = os.path.join(FIG, s)
            if os.path.exists(p):
                sides = edge_ink(p)
                if sides:
                    clipped.append((num, s, sides))
    print(f'CONSOLIDATION: {n_src + len(DROPPED)} sources -> {len(GROUPS)} numbered SI figures '
          f'({len(DROPPED)} dropped)\n')
    print(f'{"new":5s} {"file":34s} {"panels":>6s} {"printed H":>10s} {"width=":>7s}  sources')
    for num, (stem, rows, purpose) in GROUPS.items():
        flat = [s for r in rows for s in r]
        if len(flat) == 1:                          # a lone source needs no composite
            p = os.path.join(FIG, flat[0])
            if not os.path.exists(p):
                print(f'{num:5s} MISSING {flat[0]}')
                continue
            w, h = Image.open(p).size
            h_in = h * PRINT_W_IN / w
            print(f'{num:5s} {flat[0]:34s} {1:6d} {h_in:9.2f}i {recommend_width(h_in):6d}%  '
                  f'PASSTHROUGH')
        else:
            r = compose(stem, rows)
            if r is None:
                continue
            path, (w, h) = r
            h_in = h * PRINT_W_IN / w
            print(f'{num:5s} {os.path.basename(path):34s} {len(flat):6d} {h_in:9.2f}i '
                  f'{recommend_width(h_in):6d}%  ' + ' + '.join(flat))
        print(f'      {purpose}')
    print()
    if clipped:
        print('INK ON THE CANVAS EDGE (a label that ran off its own panel and was clipped):')
        for num, s, sides in clipped:
            print(f'  {num:5s} {s:36s} {", ".join(sides)}')
        print()
    for k, v in DROPPED.items():
        print(f'  DROP {k:34s} {v}')


if __name__ == '__main__':
    main()
