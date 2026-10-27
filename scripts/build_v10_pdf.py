"""Assemble full v7 long-form PDF: main text + 4 figures + Table 1 + full SI (figs S1-S4 + tables).
Run from repo root: python3 scripts/build_v7_pdf.py  (then pandoc call printed at end is executed)
"""
import re, json, subprocess, sys, pandas as pd
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]  # repo root, machine-independent (was a hardcoded stale path)
sys.path.insert(0, str(ROOT/'scripts'))
from si_figures_v12 import SI_FIGURES
MAN=ROOT/'manuscript'; FIG=MAN/'figures'; TAB=ROOT/'results/cancer_v6_6_extras'

def img(fname, caption, width='80%', land=False):
    body=f'\n\n![{caption}]({FIG/fname}){{width={width}}}\n\n'
    if not land: return body
    # Landscape was tried and abandoned: pdflscape rotates the sheet but the figure stays inside a
    # float, which is typeset after the group that widened \\linewidth has closed, so the image kept
    # its portrait size on a half-empty page. These six figures reach the 6 pt floor by font scale
    # alone (scripts/build_*.py, _SIF), which needs no LaTeX plumbing.
    return body

def csv_md(fname, floatfmt='g', idxname='cell type'):
    df=pd.read_csv(TAB/fname)
    df.columns=[idxname if str(c).startswith('Unnamed') else c for c in df.columns]
    return df.to_markdown(index=False, floatfmt=floatfmt)

_TEX_ESC=[('\\','\\textbackslash{}'),('&','\\&'),('%','\\%'),('$','\\$'),('#','\\#'),
          ('_','\\_'),('{','\\{'),('}','\\}'),('~','\\textasciitilde{}'),('^','\\textasciicircum{}')]
def _tex(x):
    if x is None or (isinstance(x,float) and pd.isna(x)): return ''
    s=('%g'%round(x,3)) if isinstance(x,float) else str(x)
    if s=='nan': return ''
    for a,b in _TEX_ESC: s=s.replace(a,b)
    return s

def csv_tex_df(df, widths, size='footnotesize', caption=None):
    """csv_tex for a frame already in hand (Table S-Q is assembled, not read from one CSV)."""
    return _tex_longtable(df, widths, size, caption)


def csv_tex(fname, widths, idxname='cell type', rename=None, drop=None, size='footnotesize',
            caption=None):
    """Wide tables as a raw longtable with p{} columns summing to \\linewidth.

    pandoc\'s pipe tables gave LaTeX `l`/`r` columns here, which cannot wrap: Table S2, S7 and S9b ran
    1-4 cm past the right margin on eleven pages. `widths` are fractions of the text block and must
    sum to <= 0.97 (the rest goes to intercolumn glue)."""
    df=pd.read_csv(TAB/fname)
    df.columns=[idxname if str(c).startswith('Unnamed') else c for c in df.columns]
    if drop: df=df.drop(columns=[c for c in drop if c in df.columns])
    if rename: df=df.rename(columns=rename)
    return _tex_longtable(df, widths, size, caption, fname)


def _tex_longtable(df, widths, size='footnotesize', caption=None, who='table'):
    """The longtable body, shared by csv_tex (one CSV) and csv_tex_df (a frame assembled here)."""
    assert len(widths)==df.shape[1], f'{who}: {len(widths)} widths for {df.shape[1]} columns'
    assert abs(sum(widths)-1.0)<=0.02, f'{who}: widths sum to {sum(widths):.2f}, expected 1.0'
    align=['>{\\raggedright\\arraybackslash}' if df[c].dtype==object else '>{\\raggedleft\\arraybackslash}'
           for c in df.columns]
    # subtract the intercolumn glue inside the column width, so the widths can sum to 1.0 exactly
    spec=''.join(f'{a}p{{\\dimexpr {w}\\linewidth-2\\tabcolsep\\relax}}' for a,w in zip(align,widths))
    head=' & '.join('\\textbf{%s}'%_tex(c) for c in df.columns)+' \\\\'
    body='\n'.join(' & '.join(_tex(v) for v in r)+' \\\\' for r in df.itertuples(index=False))
    cap = (f'\\multicolumn{{{df.shape[1]}}}{{@{{}}p{{\\linewidth}}@{{}}}}{{{caption}}} \\\\[4pt]\n'
           if caption else '')
    return ('\n\n```{=latex}\n'
            f'\\begingroup\\{size}\\setlength{{\\tabcolsep}}{{3pt}}\\renewcommand{{\\arraystretch}}{{1.06}}\n'
            f'\\begin{{longtable}}{{@{{}}{spec}@{{}}}}\n'
            f'{cap}\\toprule\n{head}\n\\midrule\n\\endfirsthead\n'
            f'\\toprule\n{head}\n\\midrule\n\\endhead\n'
            f'{body}\n\\bottomrule\n\\end{{longtable}}\n\\endgroup\n```\n\n')


def sup(s):  # keep Unicode literal; LaTeX maps each glyph via newunicodechar (header below)
    return s

UNI_HEADER = r'''% long `script/paths.py` in typewriter type were the last two overfull lines:
% allow them to hyphenate, and give TeX a little slack on the rest.
\usepackage[htt]{hyphenat}
\emergencystretch=3em
% let a short page end short rather than push a paragraph past the folio (Note S-A did)
\raggedbottom
\usepackage{longtable}
\usepackage{booktabs}
\usepackage{array}
\usepackage{caption}
\captionsetup{labelformat=empty}
\usepackage{newunicodechar}
\newunicodechar{≥}{\ensuremath{\geq}}
\newunicodechar{≤}{\ensuremath{\leq}}
\newunicodechar{×}{\ensuremath{\times}}
\newunicodechar{−}{\textminus}
\newunicodechar{→}{\ensuremath{\rightarrow}}
\newunicodechar{↔}{\ensuremath{\leftrightarrow}}
\newunicodechar{◆}{\ensuremath{\blacklozenge}}
\newunicodechar{●}{\ensuremath{\bullet}}
\newunicodechar{≈}{\ensuremath{\approx}}
\newunicodechar{ρ}{\ensuremath{\rho}}
\newunicodechar{β}{\ensuremath{\beta}}
\newunicodechar{α}{\ensuremath{\alpha}}
\newunicodechar{κ}{\ensuremath{\kappa}}
\newunicodechar{γ}{\ensuremath{\gamma}}
\newunicodechar{η}{\ensuremath{\eta}}
\newunicodechar{⁻}{\textsuperscript{\textminus}}
\newunicodechar{⁰}{\textsuperscript{0}}
\newunicodechar{¹}{\textsuperscript{1}}
\newunicodechar{²}{\textsuperscript{2}}
\newunicodechar{³}{\textsuperscript{3}}
\newunicodechar{⁴}{\textsuperscript{4}}
\newunicodechar{⁵}{\textsuperscript{5}}
\newunicodechar{⁶}{\textsuperscript{6}}
\newunicodechar{⁷}{\textsuperscript{7}}
\newunicodechar{⁸}{\textsuperscript{8}}
\newunicodechar{⁹}{\textsuperscript{9}}
\newunicodechar{₀}{\textsubscript{0}}
\newunicodechar{₁}{\textsubscript{1}}
\newunicodechar{₂}{\textsubscript{2}}
\newunicodechar{₃}{\textsubscript{3}}
\newunicodechar{₄}{\textsubscript{4}}
\newunicodechar{₅}{\textsubscript{5}}
\newunicodechar{₆}{\textsubscript{6}}
\newunicodechar{₇}{\textsubscript{7}}
\newunicodechar{₈}{\textsubscript{8}}
\newunicodechar{₉}{\textsubscript{9}}
'''

main=(MAN/'paper3_v10_longform.md').read_text()
# --- graphical abstract: placed just above the Abstract heading; EMPTY alt so pandoc
#     does NOT count it as a numbered figure (keeps every manual 'Fig. N' correct) ---
# The graphical abstract is retired, not merely unwired (decision of 25.10.2026): it was built for the
# v8 title, four of the values it printed are withdrawn, and a one-panel summary naming the
# cytotrophoblast re-centres exactly the claim §3, §4 and the Discussion now scope to this catalog and
# to readiness-level measurement. `scripts/build_v7_graphical_abstract.py` still builds the old one.

# --- main-figure injections --------------------------------------------------------------
# Every figure is injected by replacing an anchor sentence of the manuscript. An anchor that has
# since been reworded makes the replacement a silent no-op, and a stale COPY of the same caption
# sitting on a dead anchor is worse: an edit lands on the copy and never reaches the document.
# The twelve dead v9 anchors that used to sit here have been removed, and main_sub() below now
# fails the build if any remaining anchor stops matching.
_T1_ANCHOR = ('the convergent cells form a genuine neighbourhood '
              '(mutual cosine 0.77 versus −0.01 for random pairs).')
_MAIN_MISSES=[]
def main_sub(anchor, injection):
    """Append `injection` after `anchor`; record the anchor if it is not in the text."""
    global main
    if anchor not in main: _MAIN_MISSES.append(anchor[:60]); return
    main=main.replace(anchor, anchor+injection, 1)

main_sub('the fact on which everything below rests.',img('v11_fig1.png',"**Fig. 1. A cancer-independent coordinate system of cellular state, and what it recovers before any disease data are introduced.** (A) The analytical workflow: 44 regulatory modules are scored across the 154 normal human cell types of the Human Protein Atlas, and the resulting 154 x 44 activity matrix is decomposed into a nine-dimensional reference space that is then held fixed. Below, the anatomy of one module: a sensor, a transcription-factor output, its target genes, and the negative-feedback element induced by the module's own activity. No cancer data enter any step. (B) The three leading axes are individually interpretable. Each of the 154 normal cell types is placed on one principal axis (x) against an external quality signature scored from gene-level HPA data without using the modules (y); the axis is monotonically ordered by that independent signature (Spearman r = -0.89, +0.75, +0.50 for PC1-PC3), and biological families, an annotation not used to build the space, separate along each axis (eta^2 = 0.45, 0.58, 0.70). (C) Global layout: a cell's nearest neighbour in the full nine-dimensional space shares its biological family 62% of the time (95% CI 55-70%), against 14% when family labels are permuted. The comparison that matters for the panel, however, is the one in grey: size- and expression-matched random gene panels put through the identical pipeline reach 61%, and a generic highly-variable-gene decomposition reaches 68%. The space is organized, but this statistic does not by itself distinguish these modules from other panels of comparable size; the evidence for the panel is its coherence and interpretability (Note S-N).",'100%'))
main_sub('and is present in 24 of 25 cancers.',img('v11_fig2.png','**Fig. 2. Malignant cells of twenty-five lineages converge onto one region of the normal-cell space.** (A) The 154 normal cell types in the fixed nine-dimensional reference space (UMAP), coloured by biological family. (B) The same map with the 25 malignant states projected into it (diamonds, coloured by lineage of origin; normal cell types in grey). Their normal tissues of origin are spread across the whole map, while the malignant states fall into the one outlined region; the median pairwise separation on the map is 4.1 among normal cell types against 1.5 among the malignant states — illustrative only, since UMAP inter-point distances are not metrically meaningful; the quantitative claim is the nine-dimensional comparison in (C). (C) The same contraction measured in the full nine dimensions, against the null that matters: 25 equal-sized aggregates of unrelated normal cells formed by the same averaging. With magnitude included the malignant states are NOT tighter than that null (0.166 against 0.129), so most of the raw contraction is shared transcriptional output and the averaging itself. With the overall-magnitude component removed, the residual convergence in module *pattern* is real (0.026 against 0.034). Their normal tissues of origin sit at 0.99, the value expected for orthogonal states. (D) The shape-only comparison as a distribution: the observed malignant value falls below 292 of 300 matched-averaging draws, a permutation *P* = (8+1)/(300+1) = 0.03. Specification sensitivity of this test is reported in Note S-G and Fig. S10.','100%'))
main_sub('in naming one convergent *normal* target.',img('v11_fig3.png',"**Fig. 3. Where the convergence lands: a reproducible neighbourhood of normal states, and the cell nearest it under this catalog.** (A) The 25 malignant states (diamonds) on the normal-cell map, cell types coloured by the convergence they receive; the manifold's own centre (an alveolar type-1 pneumocyte, cross) receives none. (B) Cross-cancer recurrence — the number of the 25 cancers in which a cell type is a top-5 nearest neighbour; the erythrocyte and the megakaryocyte are omitted because the control in (D) rejects both; the megakaryocyte is still named as the second member of the two-cell core wherever that denominator is used, and carried no further (full ranking including both, Fig. S6-II). (C) Breadth against depth: only the two core cells lead on both axes, and of those only the cytotrophoblast survives (D). (D) The control that separates them: over deciles of transcriptome complexity, cytotrophoblast assignment rises (6.8% to 21.6%) while megakaryocyte assignment falls (13.2% to 3.5%), tracking the quarantined erythrocyte (16.6% to 3.5%); Note S-Q runs the same control for every region member. (E) Convergence per cancer and per target (column-standardized), cancers grouped by lineage of origin; a white cross is a single-cancer resonance (>25% of that cancer's malignant cells). (F) The enrichment survives widening the region (18.2-fold, 15.9-fold, 8.4-fold); occupancy inside each bar.",'90%'))
main_sub('The event also appears early and essentially all-or-none, independent of grade and stage (SI).',img('v11_fig5.png',"**Fig. 5. Malignancy increases occupancy of a pre-existing physiological state.** (A) Within-patient, cell-type-matched paired test in colorectal cancer (Pelka et al. 2021): each line is one patient, joining the fraction of that patient's own normal epithelium occupying the convergent state to the fraction of their malignant epithelium that does. Medians 4.7% and 18.9%; higher in 43 of 46 patients, paired Wilcoxon P = 9 × 10^-13. (B) The same data per patient against the diagonal: 43 of 46 points lie above y = x. Because both arms come from the same patient, the same tissue and the same protocol, lineage and batch cancel. The effect survives matching on transcriptome complexity, survives removal of the magnitude component, replicates in a second, gastric cohort, and holds against the crypt-base progenitors from which these tumours arise (SI). Baseline occupancy of normal epithelium tissue by tissue is in the SI.",'100%'))
main_sub('Details of the module signature and how the malignant occupant differs from the normal occupant of the same state are in SI.)',img('v11_fig4.png','**Fig. 4. The regulatory signature that defines the target is the one cancer reproduces.** (A) What pre-defines the target, using normal cells only: how far the convergence target sits below the other 151 normal cell types on each module, for the eighteen modules it suppresses most (the full 44 are in Table S7). No cancer data enter this panel. (B) What cancer contributes: the mean per-module contribution to the cancer-to-core convergence across the 25 cancers (+/- SD), for the same eighteen modules in the same order, so each pair can be read across. (C) The two against each other across all 44 modules: Spearman rho = +0.88 (bootstrap 95% CI 0.78-0.93, positive in 100% of 2,000 resamples). The programs that make the normal target unusual are the programs cancers reproduce as they converge on it. Type I interferon is the single most below-trend by OLS residual (third by rank residual) module - the one target-defining brake cancer does not reproduce - which is a statement about the residual: it ranks 14th of 44 on raw contribution but 1st of 44 on how far that contribution falls short of what its suppression predicts.','100%'))
# Fig. 6 was Fig. S15 until the co-author review of 2026-10-27: the isogenic series are the only
# experimental validation in the paper and were unfindable in the SI. The PNG keeps its generated
# name (scripts/build_figS29_invitro.py writes it); only its number and its place changed.
main_sub('Both scorings are reported: they date the same displacement differently, and that difference is itself the result.',
         img('v11_figS29_invitro.png', '**Fig. 6. Transformation of human cells from a known precursor: a two-stage displacement, and what the cross-origin convergence does once a non-transforming control exists.** **(A)** Ten donors of primary human prostate basal cells carrying five defined oncogenic factors (PARCB; Chen et al. 2023, GSE240057), paired within donor. Thin lines are donors, bold is the mean; the shaded band marks the arm that carries all five factors but has not yet formed a tumour. Under component-abundance scoring the cells are already 71% of the way from their own untreated basal cells to the endpoint tumour before any tumour exists (+0.360 of +0.509 in cosine to the malignant attractor, 10 of 10 donors in both; sign-test *P* = 0.001), and the stage ordering is monotone in all ten (mean Spearman rho = +0.80). **(B)** The same donors scored by regulatory output instead of component abundance (DoRothEA/PROGENy footprints, VIPER). Here the transduced-not-yet-tumour arm contributes nothing (-0.031; 1 of 10 donors positive) while the endpoint tumour remains closer to the attractor than its own precursor (+0.052, 9 of 10, *P* = 0.011; magnitude-removed +0.028, 8 of 10, *P* = 0.055; under the ULM estimator +0.401 and +0.214, 10 of 10 for both). The two scorings date the move differently: the abundance-level displacement arrives with the genetic lesion, the output-level displacement accrues with the tumour. **(C)** The discrete readout agrees with the latter - nearest-normal assignment places no donor inside the convergence region at baseline, one after transduction, and all ten at the endpoint. **(D)** Cross-origin contraction in the four series that transform two distinct origins, each beside its non-transforming control twin: in GSE222990 (primary acinar and ductal pancreas, same four drivers) substituting an acinar arm that received *no drivers at all* for the transformed one reproduces +0.42 of the +0.48 contraction; in GSE17941 the twin converges as much in raw scores and more after magnitude removal (+0.15 against +0.09); GSE118206 is the only clean positive and the only design in which no such control can be built; GSE6885 is uninformative because its two precursors start 0.02 apart. Pre-registration `perceptome_validation/PREREG_T9_invitro.md` (frozen before each run), results `T9_RESULTS.md`, `scripts/T9_invitro_transformation.py`, `T9c_trajectory.py`, `T9d_activity.py`, `T9e_pancreas.py`, figure `scripts/build_figS29_invitro.py`.', '100%'))
main_sub('awaits replication wherever paired malignant and normal cells both occupy the state in sufficient numbers.',img('v11_fig6.png','**Fig. 7. Having converged, malignant cells modify the state rather than moving further along it.** Malignant and matched-normal epithelial cells that BOTH occupy the convergent state are compared module by module within each of 30 colorectal patients (>=10 core-assigned cells per arm; paired Wilcoxon, detection depth balanced to median 2,640 versus 2,496 detected genes from a raw 4,351 versus 2,395). (A) The sixteen modules that survive the depth control, coloured by functional class; open circles give the uncontrolled value, so the shrinkage imposed by the control is visible. The malignant occupant carries higher transcriptional signatures of the stress responses (HIF, NRF2, HSF1, UPR-PERK, autophagy) and of anabolic programs (mTOR, SREBP, cell cycle, Wnt) and shifts metabolically (PPARalpha, MR, PI3K/PTEN down). Every one of these is a transcriptional signature and not a measured pathway activity; mTOR and PI3K/PTEN in particular are scored almost entirely from component transcripts, so "higher mTOR" means higher mTOR-component message. The p53 module rises through its downstream damage targets rather than TP53 itself, whose transcript is unchanged - a replication-stress footprint, not retained p53 function. (B) The same differences against the axis that defines the target: the correlation is weak (Spearman rho = +0.36 across all 44 modules), so the difference lies largely OFF the defining axis. (C) The reading: at the level of the population malignancy increases occupancy of a pre-existing state (Fig. 5), and among the cells already occupying it malignancy superimposes stress, growth and metabolic programs rather than intensifying the features that make the state special. Colorectal only; this awaits replication wherever paired malignant and normal cells both occupy the state in numbers.','100%'))

# --- Table 1 -----------------------------------------------------------------------------
# Two sentences cite "Table 1" -- the widened-region enrichment in Results 2 and the single-artifact
# statement of Note S-K -- and no such table was in the document: the v7 draft had one, the v8/v9/v10
# builders dropped it, and the pointers stayed. It is injected from the same artifact Fig. 3F is drawn
# from (`table1_boundaries.csv`, written by E41_regen_ranking.py), so the panel and the table cannot
# drift apart, and it is what the main text otherwise only counts: the names of the four and the
# fourteen. `members` carries fourteen cell names, hence the longtable rather than a pipe table.
main_sub(_T1_ANCHOR,
  csv_tex('table1_boundaries.csv', [0.17, 0.05, 0.11, 0.09, 0.07, 0.51],
          rename={'occupancy_pct': 'occupancy (%)', 'chance_pct': 'chance (%)'},
          size='scriptsize', caption=(
            r'\textbf{Table 1. The convergent region at three boundaries.} Occupancy is the share of '
            r'the 83,819 malignant cells whose nearest normal neighbour is a member; chance is the '
            r'membership size over the 152 candidates that remain after the erythrocyte-family '
            r'quarantine (Note S-K). The same three boundaries are plotted in Fig.~3F, and the '
            r'complexity control applied to every member rather than only to the two leading ones '
            r'is Note S-Q.')))


si=(MAN/'paper3_v10_SI.md').read_text()
import re as _re
_SI_MISSES=[]
def si_sub(label, repl):
    """Replace the whole SI bullet line whose label is `label` (e.g. 'Fig. S2', 'Table S7').
    Anchored on the label only, so editing a caption in the SI can never silently disable the swap."""
    global si
    pat=_re.compile(r'^- \*\*'+_re.escape(label)+r'\.?\*\*.*$', _re.M)
    si2,n=pat.subn(lambda m: repl, si, count=1)   # a callable replacement is taken literally:
    # the old repl.replace('\\','\\\\') doubled every backslash and mangled raw-LaTeX blocks
    if n!=1: _SI_MISSES.append(label)
    else: si=si2
   # the live SI (renamed from paper3_v9_SI.md — it was never a v9 file)
# --- embed SI figures replacing the bullet list lines ---
# The sixteen numbered SI figures, in order, live in scripts/si_figures_v12.py: one list, one place
# to renumber. The grouping that produced them (30 sources -> 16 figures, 3 dropped) is composed by
# scripts/build_SI_consolidate.py, which also prints the width each one should be injected at.
for _lbl, _fn, _cap, _w in SI_FIGURES:
    si_sub(_lbl, img(_fn, _cap, _w))
# --- embed CSV-backed tables ---
si_sub('Table S1', '**Table S1. Datasets (25).**\n\n'+csv_md('tableS1_datasets.csv')+
    '\n\n*malignant_label*: `author` = author/atlas malignant annotation; `Census` = malignant (tumour/epithelial) compartment from CELLxGENE Census; `marker` = marker-based Leiden clustering (Note S1); `heuristic` = heuristic label. `*` (Gastric) marks the Census gastric cohort in the panel — median 1,169 genes/cell, just above the QC floor (a separate, earlier gastric dataset that scored 912 in the 22-candidate screen was excluded; Note S3).\n')
si_sub('Table S2', '**Table S2. Full 154-cell ranking** (cross-cancer recurrence, median and mean depth, occupancy and QC status; the two quarantined post-cellular bodies carry their pre-quarantine occupancy for reference).\n\n'
 +csv_tex('tableS2_ranking154.csv',[0.31,0.08,0.11,0.10,0.11,0.18,0.11],
          rename={'cell':'cell type','recurrence_top5':'recur. (of 25)','median_depth_pct':'median depth %',
                  'mean_depth_pct':'mean depth %','occupancy_pct':'occupancy %','qc_status':'QC status',
                  'prequarantine_occupancy_pct':'pre-quar. occ. %'}))
si_sub('Table S4', '**Table S4. Per-module cancer→core contributions** (44 modules, all 25 cohorts). Each cancer centroid is standardized against the reference (mean and SD over the 154 normal cell types) and unit-normalized; the two-component signature is the unit-normalized mean of the standardized cytotrophoblast and megakaryocyte profiles; a module contribution is the elementwise term of the cosine between them, so the 44 contributions sum to that cosine — here 0.711, which is the mean over the 25 cancers of the cosine between a cancer centroid and the two-component signature in standardized module space. It is not the quantity plotted in Fig. 3A, which is the cosine of each normal cell type to the mean malignant direction in the projected nine-dimensional space (0.89–0.92 for the core). mean_contribution and sd are over the 25 cancers; pct_positive is the fraction of cancers with a positive term. (An earlier version was computed on a 16-cohort panel, which is why its values fell on a k/16 grid; recomputed by scripts/E33_fix_tables.py.)\n\n'+csv_md('tableS4_modules.csv', idxname='module'))
si_sub('Table S8', '**Table S8. Cancer → cell-of-origin mapping** (25; nominal annotation — the Fig. 3 contraction is tested against random 25-cell draws and does not depend on the pairings).\n\n'+csv_md('tableS8_origins.csv'))
si_sub('Table S9', '**Table S9. The 44 regulatory modules and their component genes** (206 unique core genes — the explicit catalog underlying the eigenspace).\n\n'+csv_md('tableS9_module_genes.csv', idxname=None))
si_sub('Table S7', '**Table S7. The core\'s target-defining module signature** (44 modules; target-distinctiveness vs cancer contribution).\n\n'+csv_tex('tableS7_core_signature.csv',[0.23,0.16,0.15,0.15,0.16,0.15],idxname='module',
          drop=['cancer_contribution_16cohort_retracted'],
          rename={'target_distinctiveness':'target distinct.','core_suppression':'core suppression',
                  'megakaryocyte_z':'megakaryocyte z','cytotrophoblast_z':'cytotrophoblast z',
                  'cancer_contribution':'cancer contrib.'}))
# Table S6 used to print all 144 rows, three pages of a machine-readable diagnostic that also
# ships as tableS6_complexity.csv. Only its two ends carry an argument -- the low-complexity floor
# the quarantine rests on, and the complexity of the cells the core is made of -- so the document
# keeps the twenty most-assigned targets and the five lowest-complexity ones, and names the file.
_s6 = pd.read_csv(TAB / 'tableS6_complexity.csv')
_s6top = _s6.nlargest(20, 'n_cells_assigned')
_s6low = _s6[~_s6.cell.isin(_s6top.cell)].nsmallest(5, 'median_ngenes')
_s6sel = pd.concat([_s6top, _s6low]).round({'mean_ngenes': 0}).astype({'mean_ngenes': int})
_s6sel.columns = ['cell type', 'cells assigned', 'median genes', 'mean genes']
si_sub('Table S6', '**Table S6. Transcriptome-complexity diagnostic** (median detected genes per '
    'assigned cell type). The twenty cell types receiving the most assignments, then the five '
    'lowest-complexity types below them; all %d types that receive at least one assignment, of the '
    '154 in the reference, are in the attached `tableS6_complexity.csv`. The QC floor is 1,000 '
    'detected genes per cell, and the two quarantined post-cellular bodies sit below it.\n\n'
    % len(_s6) + _s6sel.to_markdown(index=False))
si=si+'\n\n```{=latex}\n\\clearpage\n```\n\n**Note S-A. Robustness of the convergence to the assignment rule and to transcriptional magnitude (`scripts/E3_referee_controls.py`).** Re-scoring the 25 cancers (83,819 malignant cells) and re-projecting: the winner-take-all top-1 core fraction (27%) is matched by a soft cosine-softmax assignment (21% core mass) and by top-k membership (core in the top 3 of 152 for 57% of cells, top 5 for 75%; median nearest-core rank 3), so the convergence is not an artifact of hard assignment. It is not the low-complexity/magnitude artifact that inflated the erythrocyte: core-assigned cells carry higher transcriptional complexity than non-core (median 2,629 vs 2,360 detected genes), and cells assigned to a centroid do not track that centroid\'s projection norm (r = +0.14). Applying the same magnitude-removing transform to BOTH reference and cancer (within-cell centering, or per-cell rank transform) reduces the core convergence from 27% to 7–11% (still 4–5× over the ~2% chance baseline): the effect is largely a signal of overall regulatory output (the terminal-high-output biology) with a smaller, real module-shape residual.\n\n'
si=si+'\n\n**Note S-B. Calibration dependence (`scripts/E4_calibration.py`).** Scoring the 25 cancers once and projecting under each reference calibration: under the single-cell-scaled calibration the core fraction is 17% and the 25 malignant centroids are tighter than every one of 2,000 random 25-cell draws; under a pseudobulk (population-averaged) calibration the core fraction falls to 7% and the projection degenerates, collapsing tumor cells toward low-complexity non-core targets (mean pairwise centroid distance 0.001). Convergence onto the core is present under both calibrations but is properly resolved only under the single-cell calibration — the appropriate one for single-cell tumor data, and the choice we adopt on unit-matching grounds (not because it maximizes the effect).\n\n'
si=si+'\n\n**Note S-C. The address is specific, not merely a high-output signal (`scripts/E5_specificity.py`, 83,819 malignant cells, 25 cancers).** Four independent tests separate *which cell is named* from *overall transcriptional level*. (i) *Shape-only ranking.* Applying the magnitude-removing rank transform to both reference and cancer and re-ranking the cross-cancer recurrence, the cytotrophoblast remains the single most recurrent normal target (14/25 cancers) and the syncytiotrophoblast — a second trophoblast state — rises to second (9/25); the megakaryocyte falls from tied-first to sixth. Magnitude removal therefore *sharpens* the address onto the trophoblast lineage (the named anchor) and reveals that the residual magnitude dependence lives in the megakaryocyte arm, not in the cytotrophoblast. (ii) *High-output specificity.* The professional high-output secretory normals that a pure-magnitude account predicts should also attract convergence do not: of all cross-cancer top-1 votes, hepatocytes receive 0.11%, plasma cells 0.15%, parietal cells 0.69%, prostatic glandular 0.03%, proximal-tubule 0.00% — against 11.1% (cytotrophoblast) and 12.8% (megakaryocyte) for the core. The only high-baseline normals that draw appreciable convergence, pancreatic acinar (7.5%) and gastric chief (7.3%), are the resident secretory cells of two cancers in the panel (PDAC, gastric) — tissue-of-origin return rather than cross-lineage convergence. (iii) *Library-size control.* Core assignment is essentially uncorrelated with sequencing depth: r(core, log10 total counts) = +0.03, r(core, detected genes) = +0.04. (iv) *Orthogonal markers.* Scored on trophoblast/EVT genes lying entirely outside the 206 module genes, core-assigned malignant cells express the extravillous-trophoblast invasion program more highly than non-core (MMP2 0.140 vs 0.082; MMP9 0.098 vs 0.042; ITGA5 0.090 vs 0.064), though the canonical identity markers (HLA-G, CGA/CGB, KRT7, TFAP2C) are not specifically enriched — a modest, invasion-flavored independent corroboration of trophoblast-ness, not a full identity match.\n\n'
si=si+'\n\n**Note S-D. Within-patient deepening under a mixed-effects model.** The pooled sign test (77/86 patients deeper) treats patients as exchangeable across the four low-baseline datasets, which risks a Simpson\'s-paradox confound. A linear mixed-effects model of the per-patient malignant-minus-normal core difference with a dataset random intercept estimates a fixed deepening of +13.0 percentage points (SE 1.87; z = 6.95; one-sided *P* = 1.9 × 10⁻¹²), with all four datasets deepening in the mean (colorectal-Pelka 43/46, colorectal-HTA 23/25, kidney 6/7, liver 5/8). The effect thus survives crediting the between-dataset clustering; we note the coarseness of a four-group random-effect estimate (two of the four cohorts colorectal), so the low-baseline conclusion is carried by colorectal with same-signed support from kidney and liver.\n\n'
si=si+'\n\n**Note S-F. The deepening survives within-pair complexity matching (`scripts/E15_deepening_complexity_matched.py`).** Malignant epithelium carries more detected genes than matched-normal (colorectal median 3,852 vs 2,272), and core assignment is weakly complexity-dependent, so a within-pair complexity offset could in principle manufacture deepening. Resampling the malignant cells of each patient to the *matched-normal* detected-gene distribution (nearest-gene-count matching, 20 draws) and recomputing the core fraction: the deepening falls only from +15.5 to +12.9 percentage points (43/46 patients, paired *P* = 10⁻⁸), retaining ~80% of the raw effect. The paired complexity offset therefore accounts for at most a fifth of the deepening.\n\n'
si=si+'\n\n**Note S-G. The contraction is mostly magnitude and averaging, with a real shape residue (`scripts/E14_contraction_matched_null.py`).** The raw 25-cancer clustering (mean pairwise cosine distance ~0.166, vs ~0.99 among the near-orthogonal origins) is largely a matched-averaging and shared-magnitude effect: 25 equal-sized coherent normal micro-populations, formed by the same averaging over a seven-tissue normal pool, cluster nearly as tightly (null mean 0.129, 95% interval [0.074, 0.213]), and the malignant clustering is flat across subsample size (0.174 at 25 cells/cancer, 0.166 at 4,000), so it is a stable biological quantity rather than the shrinkage of large means. Removing the overall-magnitude component isolates a residue that we report as a specification curve rather than a single number, because the comparison is sensitive to how the magnitude-removed scores are rescaled before projection. (i) In the **raw module space**, with no standardization or projection at all, the malignant centroids cluster at 0.326 against a null mean of 0.493 — below every one of 200 draws (and in the same space the raw metric is 0.168 vs 0.287, likewise below every draw). (ii) Under the **reference standardization used throughout the paper**, the shape-only clustering is 0.026 against a null of 0.034 (the observed value falls below 292 of 300 matched-averaging draws; permutation *P* = (8+1)/(300+1) = 0.030). (iii) Under a third standardization — dividing by the standard deviations of the *magnitude-removed* reference, which up-weights modules whose pattern barely varies across the 154 normal types and so amplifies low-variance axes — the advantage vanishes (0.400 vs 0.382, 66th percentile). Two of three specifications, including the simplest and assumption-free one, support a real shape convergence; one does not. We report all three and no longer headline the magnitude-inflated 6.3-fold figure. (The coherent-micro-population null is if anything biased tight, making the surviving results conservative; a cell-type-level null is the natural refinement.)\n\n'
si=si+'\n\n**Note S-H. Cancer-blind justification of the single-cell calibration (`scripts/E16_calibration_cancerblind.py`).** To show the single-cell-scaled calibration is chosen on principled, not result-maximizing, grounds, we test it on 16,312 independent *normal* single cells (skin, kidney, breast; author cell types; no cancer). Projected and assigned to their nearest normal reference type: the **pseudobulk** calibration degenerates — it collapses these normal single cells onto only 2 of 152 targets (effective assignment diversity 1.0), 99.6% of them forced onto low-complexity erythroid cells — the same degeneration it produces on tumor cells (Note S-B). The **single-cell-scaled** calibration keeps the assignments diverse and complexity-appropriate (85 targets, diversity 14.6). Pseudobulk is thus wrong for single-cell data of any kind, normal or malignant; the single-cell calibration is required by the data type and selected with no cancer input. (A fine-cell-type silhouette metric did not distinguish the calibrations; we report the degeneration test, not the metric that would have favored our choice.)\n\n'
si=si+'\n\n**Note S-I. Malignant-cell purity and the direction of its bias.** The lung dilution artifact (Fig. S11-II) shows automated malignant selection can be wrong in a specific cohort, so we state the purity basis of each within-patient cohort and the direction of any residual bias. The cohorts carrying the deepening use author-annotated malignant compartments — Pelka `ClusterMidway = EpiT` (tumor epithelial clusters) and lung Kim-2020 tumor subtypes tS1/tS2/tS3 — which are high-purity; the gastric arm uses an EPCAM+/CD45− gate on the cancer sample with no malignant annotation, so it may retain some normal-adjacent epithelium. Crucially the bias is signed: normal-adjacent contamination pulls a tumor arm toward its normal baseline and therefore *understates* deepening (the lung case: a 9%-malignant arm read ~13% ≈ normal, versus ~29% for the clean malignant cells; `E13`). Marker-based selection (germ-cell/sarcoma) biases convergence the other way and is bounded by the drop-4 robustness. We do not, however, claim that contamination makes the reported numbers a lower bound: we tested that argument directly and it was not borne out (below). We tested the contamination argument directly on the two cohorts whose labels allow it — restricting HCC to hepatocytes from tumour sites and BRCA to annotated cancer-epithelial cells — and it is *not* borne out: cleaning lowered rather than raised those cohorts (HCC 24.3% → 18.6%, BRCA 31.0% → 26.7%) and left the panel essentially unchanged (core 27.1% → 26.6% raw, 17.3% → 17.1% shape-only; clustering 0.166 → 0.166). Contamination therefore does not appear to be masking convergence, and we report the direction-of-bias argument as untested-and-unsupported rather than as a defence. A CNV-inference malignant floor across all 25 cohorts is deferred to production, as the convergence is now presented as the modest, magnitude-carried companion (Note S-G). Separately, on the deepening-vs-baseline relationship: across five cohorts the additive deepening does not depend on baseline (Pearson *r* = +0.4, *P* = 0.5; mean ~+10 points), the fold-change inverse being a small-denominator artifact — so we make no baseline-gating claim (Fig. S11-IIIB).\n\n'
si=si+'\n\n**Note S-J. What moves the state, and what does not (perturbation and immune controls).** *Single genes do not.* Across the genome-wide Perturb-seq resource (Replogle et al. 2022; 9,675 CRISPRi knockdowns in K562 with 75,328 non-targeting controls), knocking down a target-defining brake does not move cells toward the core: the correlation between a gene\'s brake score and its knockdown\'s core shift is +0.07 (*P* = 0.4), +0.13 after controlling for the generic brightening that essential-gene knockdown causes (*P* = 0.1). A three-level decomposition shows where the signal stops (`scripts/E21_why_null_diagnostic.py`): CRISPRi works (target mRNA down, median −0.43 SD), the containing module responds strongly (LDHA→HIF −1.84 SD, HSP90AA1→HSF1 −1.65; and, as a positive control that the modules read pathway biology rather than their own transcripts, KEAP1 knockdown *raises* the NRF2 module by +0.96 SD), but among the 23 strong module hits the nine-dimensional position moves only a fifth as much and in no consistent direction (13 of 23 toward the core; *r* = −0.10 between module movement and coordinate movement). We read this as a many-module consensus buffered against single-gene perturbation — consistent with a state defined by a dozen simultaneously suppressed brakes, and with the large movement seen under whole-network pressure (Fig. S14). We label that reading for what it is: a post-hoc account of a failed causal prediction, not a tested one. Two limits bound it. The screen is K562, a single erythroleukemia line and an unfavourable substrate for a state anchored on trophoblast and epithelial biology; we do not report its unperturbed position in the coordinate, so we cannot exclude that the cells had little room to move toward the core in the first place — that baseline should accompany any future use of this screen. And no combinatorial knockdown was run, so the “coordinated pressure” that the account requires is untested rather than supported. *Immune privilege is not inherited.* Of the three properties that motivate the cytotrophoblast — invasion, proliferation, immune evasion — the first is corroborated by orthogonal markers (Note S-C) but the third is not: across 62 colorectal tumors the malignant core fraction is unrelated to immune infiltration (Spearman ρ = +0.13 raw and +0.04 shape-only for T/NK cells, *P* = 0.3 and 0.8; myeloid +0.16 and +0.03), and stratifying by mismatch-repair status — the dominant confounder — leaves it unchanged (MMRd ρ = +0.06, n = 34; MMRp ρ = +0.08, n = 28). Sixty-two tumors exclude a moderate association. This is what the coordinate predicts: Type I interferon is the one target-defining brake cancer does not reproduce, so cancer reaches the address without its immune privilege. Cell-fraction infiltration from dissociated tissue is noisy and infiltration is not exclusion; spatial and functional tests remain open. `scripts/E20_gwps_causal.py`, `scripts/E27_immune_exclusion.py`.\n\n'
si=si+'\n\n'+'**Note S-K. Definitions, denominators and how the reported percentages relate (reader key).** Three choices govern every percentage in this paper and we fix them here. (i) *Gallery denominator.* The reference holds 154 normal cell types; the two post-cellular bodies (erythrocyte, platelet) are quarantined as a documented low-complexity artifact, leaving **152 candidates** for every assignment. The quarantine is applied *before* the nearest-neighbour step, not after it: the two cells are removed from the candidate pool and the cells that would have gone to them are reassigned to their best surviving neighbour, rather than being dropped from the denominator. Two implementations of this existed in earlier drafts and gave different top-4 boundaries; `scripts/E41_regen_ranking.py` now regenerates the ranking, the boundary table and the per-cancer matrix from one assignment, and Table S2, Table 1, Fig. 3B, 3C, 3E and 3F are all derived from that single artifact. Ranks are quoted out of 154 only where the quarantined cells are deliberately shown. (ii) *Which cells define the measured state, and where each definition is used.* Two definitions appear and we separate them explicitly, because they belong to different sections. The **ranking statistics of §3** — cross-cancer occupancy and its enrichment — are quoted for the single most recurrent cell (cytotrophoblast alone, chance 1/152 = 0.66%) and for the **two-cell core** (cytotrophoblast + megakaryocyte, chance 2/152 = 1.3%). The **within-patient occupancy of §5, which is the paper’s primary assay, uses a three-cell definition** that additionally includes the migrating cytotrophoblast, the trophoblast state adjacent to the first (chance 3/152 = 2.0%); every number in §5 and every control quoting a “~2% chance baseline” is on that definition (`scripts/pelka_perpatient.py`, `E10`, `E12`, `E17`, `E36`). An earlier version of this note stated that the main text used the two-cell definition throughout, which was wrong for §5. Splitting the colorectal effect by member shows why the distinction matters: the trophoblast states carry +16.6 points of the increase and the megakaryocyte −0.8. (iii) *Metric.* The primary statistic is the hard top-1 assignment under the single-cell-scaled calibration. The values a reader will meet are therefore: **23.9%** (top-2 core, main text), **27%** (three-cell core, top-1, the 83,819-cell control panel) and **17%** (three-cell core with the eigenspace rebuilt per calibration, Note S-B; that run has no shape-only variant, so 17% and 17.3% below are *not* the same number). Soft assignment gives 21% core mass for the same cells, and Note S-C reports the two core cells separately (11.1% + 12.8% = 23.9% of cross-cancer votes) rather than as one fraction. These are the same phenomenon under different definitions, not conflicting measurements; the numbers the abstract leads with are the raw ones, and each shape-only value is labelled as such where it appears. (iv) *Shape-only, defined.* “Shape-only” means: subtract from every cell (reference and malignant alike) its own mean across the 44 module scores, `X - X.mean(axis=1, keepdims=True)`, **before** the reference standardization, and recompute the reference mean and SD on the centred reference; projection, normalization and top-1 assignment are otherwise unchanged. This is the convention used in the main text and throughout the figures, and under it the three-cell core falls from 27.1% to **17.3%** (`scripts/E23_purity_fixed_convergence.py`). Note S-A quotes **7–11%** for what reads as the same operation: it is not. There the centring is applied *after* z-standardization, against the unchanged reference statistics (`scripts/E3_referee_controls.py`), which removes more than the magnitude component; on the same 83,819 cells the two orders give 7.4% and 10.7% against 17.3%. Both are honest shape-only variants and both leave the core four- to fivefold over chance, but they are different transforms and we do not average or interchange them.'+'\n\n'
si=si+'\n\n'+'**Note S-L. Five further methodological controls.** *(i) Map validity, with the right null.* A cell\'s nearest neighbour in the nine-dimensional space shares its biological family in **62% of the 152 candidate cell types (95% CI 55–70%)**. The relevant null is not 1/K over families (which would be 7% and is too permissive of the claim) but a family-size-preserving one: permuting the family labels 10,000 times over the same neighbour graph gives a null mean of **14%** (maximum 28%) and a permutation **P = 1 × 10⁻⁴**, which is the floor 10,000 draws can reach rather than a measured value. The 15% previously quoted as "chance" corresponds to this size-preserving rate, not to 1/K. *(ii) Modules share a few genes.* Twelve genes appear in more than one module (ADIPOQ, CALM1, CAMK2A, CCND1, DUSP1, FKBP5, FOXO1, GSK3B, PGR, SGK1, STAT1, XBP1), touching 17 of the 44 modules; 219 gene slots reduce to 206 unique genes. Of the six module pairs with |r| > 0.7, **two are partly definitional** — UPR-ATF6 with UPR-IRE1 (share XBP1) and JAK-STAT with Type I IFN (share STAT1) — while four (ERK/MAPK–NF-κB, UPR-PERK–p53, NF-κB–UPR-PERK, BMP–Hippo) share no genes and are biological. We therefore no longer describe all six as "biologically expected". Recomputing the two definitional correlations on disjoint gene sets requires gene-level reference scores and is listed as a production item. *(iii) Dimensionality.* Three of the nine retained components are the least determined under bootstrap (PC7–PC9, 0.53–0.71 over 500 draws), so we report the result at six as well as nine dimensions: the core fraction is 18.9% (6-D) versus 27.1% (9-D) and the contraction 0.131 versus 0.166 — and in six dimensions the **cytotrophoblast becomes the single top recurrent target**, ahead of the megakaryocyte. The magnitude depends on dimensionality; the anchor does not. *(iv) Leave-one-cancer-out.* Dropping each of the 25 cohorts in turn — not only the four marker-labelled ones — the core fraction stays within **25.4–28.4%** (full panel 27.1%) and the cytotrophoblast\'s recurrence rank is **2 in every one of the 25 panels** (top-3 in 25/25); the most influential single cohort is B-ALL. No cohort drives the result. *(v) Specification multiplicity.* The shape-only convergence is supported by two of three scalings — in raw module space it falls below every one of 200 null draws, under the paper\'s reference scaling the permutation P is 0.030, and under the magnitude-removed-reference scaling there is no effect. A Bonferroni correction across the three specifications leaves the raw-module-space result significant and renders the reference-scaled one marginal (P = 0.09); we report the curve rather than a single corrected number. `scripts/E34_methodology_controls.py`.'+'\n\n'
si=si+'\n\n'+"**Note S-M. Supporting detail for six claims made in the Results.** (i) *Lung correction, full account.* On the reference lung data themselves the cleanly identified malignant cells sit about 18 points above the normal epithelium with which they had been averaged; the diluted arm was ~9% annotated-malignant (Fig. S11-II; scripts/E12_highbaseline_lung.py, scripts/E13_lung_dilution_proof.py). (ii) *Prostate.* A fifteen-patient paired atlas (CELLxGENE) shows no within-patient increase (7/15) at a normal baseline of ~15%; an independent paired cohort (Song et al. 2022, GSE176031, four patients) puts the benign baseline at 2.1% and the tumour arm at 3.8% (paired P = 0.38). The two cohorts disagree about the baseline and neither shows an increase, so prostate is reported as an open negative. (iii) *The fold-change impression.* Across five tissues the added occupancy is flat against baseline on the additive scale (Pearson r = +0.4, P = 0.5) and not significant on a bounded logit scale; the appearance of more room where the baseline is lower comes from dividing by colorectal's 4.7% denominator. (iv) *Fibrotic lung, per cell type.* Within IPF epithelium the shape-only core fraction is highest in transitional AT2 (18.4%) and proliferating epithelium (16.2%) and lowest in ciliated cells (5.2%), all below the malignant lung level of 30.3% (Fig. S12). (v) *Stage independence.* The occupancy increase is already saturated at the earliest stages sampled and is independent of tumour grade and stage in the best-powered cohort. (vi) *Germline alternative.* The germline matches the core on PC1 (cosine 1.00) but diverges across the remaining eight components; by the full recurrence-and-depth ranking the nearest germline cells fall well below the core, stable across 2,000 bootstrap resamples of the 25 cancers."+'\n\n'
si=si+'\n\n'+"**Table S9b. The four operational criteria, per module (all 44).** For each module: the signal-sensing component, the transcription-factor output, the negative-feedback element and the receptor field. Modules annotated at construction are marked `recorded`; the remainder were curated from the catalog's own gene lists together with canonical pathway biology by two independent curations that were then reconciled, and carry a confidence grade. A feedback element is named for 38 of 44 modules; six are marked `unclear` rather than filled in — the five for which no induced feedback gene could be named with confidence (AR, PR, PPARgamma, NPAS4, PI3K/PTEN), and Hippo, whose recorded entry (NKD1, NKD2, BIRC2, BIRC3, BIRC5) we withdraw: NKD1 and NKD2 are Wnt antagonists and appear correctly in the Wnt row, so this was a transcription error, and the BIRCs are YAP targets but not feedback inhibitors. **A larger gap, which this table documents design intent rather than composition.** Cross-checking against Table S9 — the genes actually scored — the named feedback element is **absent from the scored gene list in 11 of the 38 modules** that name one: AMPK (ULK1), Autophagy (RRAGD/FLCN), Calcium (RCAN1), Cell Cycle (E2F7/E2F8), ER (NRIP1), Hippo, Notch (NUMB/NUMBL/NRARP), SREBP (INSIG1/INSIG2/PCSK9), TR (DIO3), cGAS-STING (TREX1/DNASE2/SAMHD1) and mTOR (DEPTOR/DDIT4/GRB10). For those modules the score cannot be reading the feedback arm at all, which matters because the readability argument in Note S-N(ii) rests on the induced arm. Taken with the two failed validations in Note S-N(iv), the fair description is that the four criteria record how the panel was curated, not a property of the 206 genes as scored — which also helps explain why three curators given the criteria alone built materially different catalogs (Note S-O). The receptor-field column is blank for sixteen modules (AhR, Circadian, ERK/MAPK, HIF, HSF1, Hedgehog, Hippo, NF-kB, Notch, SREBP, UPR-ATF6, VDR, Wnt, cGAS-STING, mTOR, p53), and so is their feedback-mechanism sentence. These are exactly the sixteen whose annotation was `recorded` at catalog construction, where the sensor was entered as a gene list without a prose ligand-class description; all twenty-eight re-curated here carry one. We therefore report receptor-field specificity as documented for 28 of 44 and as asserted at construction but never written down for the other 16 — a gap in the record, which we state rather than fill in retrospectively. Printed in abbreviated form; the complete table, with the mechanism of each feedback loop and the provenance of each annotation, is the attached `tableS9b_module_criteria.csv`."+'\n\n'+csv_tex('tableS9b_print.csv',[0.13,0.19,0.22,0.18,0.21,0.07],size='scriptsize')
si=si+'\n\n'+"**Note S-N. The module panel, evaluated without cancer data.** Everything in this note is computed on the 154 normal cell types; no cancer data enter. *(i) Coherence — the modules are units.* For each module we take the mean pairwise Pearson correlation among its genes across the 154 normal cell types. Median within-module coherence is +0.154 for the 44 modules (median 5 genes each), against +0.014 for 200 random gene sets matched on size and drawn from expressed genes (Mann-Whitney P = 3.2e-15), and against +0.062 for KEGG signalling maps (median 121 genes), +0.058 for the largest Reactome pathways (153), +0.070 for Hallmark (179) and +0.051 for PROGENy (99); the modules beat every collection at P = 1e-4 to 6e-9 while carrying twenty- to thirty-fold fewer genes. *(ii) Readability — why a transcriptional output is required.* Grouping every annotated gene of every module by the role it plays and taking its standard deviation across the 154 normal cell types, the induced arm (feedback elements and TF targets, n = 387) has median SD 0.287 against 0.222 for the structural arm (sensors and cascade components, n = 704), P = 2.8e-05. A pathway is readable from RNA through what it induces, not through the abundance of its receptors and kinases. *(iii) Catalog size.* Random subsets of k of the 44 modules, with the eigenspace rebuilt from the normal cells alone, recover nearest-neighbour biological-family concordance of 38.7% (k=5), 54.5% (15), 61.2% (30), 63.0% (40) and 65.8% (44), and agree with the sampled-catalog cell-to-cell geometry at r = 0.53, 0.78, 0.92, 0.976 and 1.0 — the value at k = 44 being a comparison with itself, not an independent agreement. Ninety-five per cent of the full value is reached at k = 35 on both measures and both are still rising at 44. One caveat on the endpoint: this curve rebuilds the eigenspace from each sampled subset, so its k = 44 value (65.8%, 100 of 152) is not the 62% (95 of 152) quoted in the main text and in Fig. 1C, which is measured in the shipped, fixed nine-dimensional reference (`reference_v03`). The two spaces are near-identical but not identical (cell-to-cell geometry r = 0.964), and the 3.3-point gap is the size of that difference, not a disagreement about the statistic. The catalog therefore sits at the shoulder of saturation. This bounds whether 44 is sufficient; it cannot say whether a larger catalog built on different criteria would do better. *(iv) What we could NOT establish — the negative-feedback criterion.* We tried twice to support criterion (iii) empirically and report both failures. Against the Cancer Gene Census (v103, 740 genes), the genes occupying the feedback position of a module are NOT enriched for cancer-driver status relative to the other genes of the same modules: 17/114 (14.9%) versus 128/638 (20.1%), OR = 0.70, P = 0.25, and significantly in the opposite direction against tighter backgrounds (core genes only, OR = 0.33, P = 0.004). And in (ii) above the feedback genes alone are only marginally more variable than cascade genes (0.242 versus 0.222, P = 0.16) — the induced-arm effect is carried by TF targets. We therefore present negative feedback as the panel's organizing design principle and not as an empirically validated selection rule. `scripts/E37_feedback_cgc_enrichment.py`, `perceptome_validation/scripts/E38_cancerfree_module_justification.py`, `perceptome_validation/scripts/E39_saturation.py`."+'\n\n'
si=si+'\n\n'+"**Note S-Q. Two audits: the complexity control applied to the whole convergence region, and gene-symbol recovery per module per cohort.** *(i) The complexity control, every region member.* Fig. 3D ran the transcriptome-complexity decile control for the two core cells and the erythrocyte only. Table S6 shows five further members of the top-14 region at or below the quarantined platelet (1,541 median detected genes) and three at or below the quarantined erythrocyte (890): late spermatids 763, gastric chief 904, pancreatic acinar 992. We therefore ran the same control for every tracked target, on 312,220 malignant cells from 17 cohorts, with the identical QC, cap, seed and assignment rule (`scripts/E43_region_complexity_control.py`). **Six of the fourteen region members fall with complexity exactly as the erythrocyte does** — megakaryocyte (13.2% in the lowest decile to 3.5% in the highest), late spermatids (5.9 to 0.2), gastric chief (7.8 to 1.7), oesophageal apical (5.0 to 2.0), pancreatic acinar (13.4 to 0.5) and erythrocyte progenitors — while eight rise, the cytotrophoblast most steeply (6.4 to 15.7). *What this costs the boundary-independence argument: almost nothing.* Dropping the six and recomputing both the occupancy and the chance baseline, enrichment goes from 16.0× to 15.8× at the core, from 11.5× to **12.5×** at the neighbourhood (it rises, because the two survivors are the two strongest targets) and from 6.6× to 6.0× at the broader region. The region is therefore not an artifact class: widening it does add complexity-suspect cells, but the enrichment does not depend on them. A sharper form of the worry — that four of the five cell types common to all four independent catalogs (Note S-O) are in the quarantine complexity range — stands as stated, and we now say so: the cross-catalog agreement on gastric chief, pancreatic acinar and oesophageal apical is agreement on targets this control rejects. What survives the control **and** reproduces across catalogs is the gastric progenitor (2,506 detected genes, rises with complexity) alongside the cytotrophoblast. *(ii) Gene-symbol recovery.* The pipeline reports per-module coverage but never printed it. Auditing all 44 modules against all 17 cohorts found three genes systematically lost, for two different reasons. **`NKX3-1` was lost in all 17 cohorts** — not a synonym problem but our own de-duplication: the loader stripped any trailing `-<digits>` to undo `var_names_make_unique()`, and `NKX3-1` is a real symbol, so it became `NKX3` and matched nothing. The HPA reference, read from a different file, kept it, so the AR module was scored on five genes in the reference and four in every cancer cohort. **`MB21D1` and `TMEM173` (cGAS-STING) and `ARNTL` (Circadian) are deprecated symbols** (now CGAS, STING1 and BMAL1). The HPA source file uses the current forms, so these genes matched *neither* the reference nor most cohorts: cGAS-STING has in fact been scored throughout on IRF3, IFNB1 and CXCL10 — without either of its sensors — and Circadian on four genes of five. `scripts/gene_names.py` now resolves both classes of name and recovers AR to 84 of 85 module×cohort gene slots (from 67), Circadian and cGAS-STING to 85 of 85 (from 73 and 61); `results/cancer_v6_6_extras/E44_gene_recovery_fix.csv`. *What we withdraw because of this.* cGAS-STING is the one module with inverted sign in Table S7 and the largest cross-cohort variance in Table S4, and the Results sentence that malignant occupants run it hotter is therefore **withdrawn**: a three-gene fragment missing both sensors cannot carry that reading. The rest of the catalog is bounded by the existing module-drop robustness (Fig. S7-II: dropping five modules at random leaves the core fraction at 17.5 ± 2.5% and no single module is load-bearing), so the headline results do not rest on these three modules — but the corrected catalog has not been re-derived end to end, and that re-derivation is the one substantive step this revision does not complete."+'\n\n'
si=si+'\n\n'+"**Note S-R. The five isogenic transformation series: pre-registration, per-contrast detail and the gates (Fig. 6).** The design, the two-stage displacement and the cross-origin control are in the main text. Across the nine primary within-series contrasts five reach the pre-stated significance level and four do not, and two further candidate series were set aside by the senescence gate — oncogene-induced senescence in the transduced arm (GSE22573, GSE172029), which the pre-registration defines as uninterpretable rather than negative. The clean contraction of GSE118206 holds under the alternative reference calibration (+0.302, *P* = 0.013 raw; +0.161, *P* = 0.009 magnitude-removed) and is visible in a plain highly-variable-gene view of the same arms (precursors 0.172 apart, transformed 0.023), so the module panel sharpens an effect that is present without it; against 200 expression-matched random gene panels the displacement is specific in raw scores (*P* = 0.005) but not after magnitude removal (*P* = 0.06). Pre-registration `perceptome_validation/PREREG_T9_invitro.md`, frozen before each run; results `T9_RESULTS.md`; scripts `T9_invitro_transformation.py`, `T9c_trajectory.py`, `T9d_activity.py`, `T9e_pancreas.py`."
si=si+'\n\n'+"**Note S-O. How much of the result depends on this particular catalog, and on which normal cells are compared.** *(i) Independent reconstruction of the catalog.* Three curations were produced by independent large-language-model agents from the four operational criteria and the size target alone, with the incumbent catalog, this paper, and the words trophoblast and megakaryocyte withheld; each shares about half of the incumbent's 206 genes (Jaccard 0.28-0.31) while agreeing with one another far more closely (0.62-0.69). What reproduces is the REGION: the top-ten nearest targets overlap the incumbent's by 5-7 of 10, and five cell types - gastric progenitor, gastric chief, pancreatic acinar, oesophageal apical and megakaryocyte - are in the top ten of all four catalogs. What does not reproduce is the ANCHOR: all three curations place the cytotrophoblast 26th-35th of 152 by cosine to the malignant attractor (0.507-0.639 against 0.955-0.963 for their own winner, the gastric progenitor), and 0 of 3 put it in the top five. The cosine band is flat only at the very top (rank 1 to rank 10 is 0.96 to 0.84) and has fallen to 0.51-0.62 by rank 30, so this is a real separation and not an artifact of ordinal ranking. Between-catalog variance, which the framework previously did not have: core occupancy 12.0 +/- 4.3% against 24.7% for the incumbent — the three-cell core on the rebuilt eigenspace this comparison uses, which is the right comparator here and is not the 23.9% two-cell main-text figure (Note S-K). Caveat stated plainly: the three curators are language-model agents of one model family and share its priors (Use of generative AI), which explains their agreement with each other but not their disagreement with the incumbent. *(ii) Established pathway collections.* KEGG (46 signalling maps), Reactome (44 largest pathways), Hallmark (50) and PROGENy (14), each scored through the identical pipeline, reproduce the within-patient occupancy increase at a comparable size (KEGG +15.0 points shape-only, 44/46 patients; PROGENy +12.9; Hallmark +19.2 raw, 46/46; Reactome, the weakest, +6.5 raw and +2.5 shape-only, 33/46, *P* = 0.014; incumbent +9.0 shape-only and +13.9 raw), and under every one of them the largest within-patient shift of colorectal epithelium is toward a trophoblast state. All six numbers in that comparison, the incumbent's included, come from one harness run under a single common setting — pseudobulk loadings with the eigenspace rebuilt per panel — chosen so the collections are scored identically. That is why the incumbent reads +9.0/+13.9 here and +10.5/+15.5 in the main text and Note S-F, where it is scored under the single-cell-matched calibration the tumour data require: same 46 patients, same cells, different calibration (under which raw and shape-only nearly exchange places, 8.81/10.98). The comparison between collections is therefore valid; only the incumbent's absolute value differs from the body of the paper (`perceptome_validation/registry/provenance_notes.md`). None of them reproduces the cross-lineage cytotrophoblast target: their core fractions are 0.05-0.29 times the incumbent's and they place the cytotrophoblast 12th-34th, recurrent in 0-2 of the 25 cancers. The occupancy increase is therefore catalog-independent; the cross-lineage generality of the target is not. *(iii) Which normal cells the comparison uses.* Colorectal cancer arises from crypt-base stem and transit-amplifying cells, not from the tissue at large, so the matched-normal arm must be the right subtype. Scoring whole transcriptomes on 2,000 highly variable genes rather than on the modules, the malignant epithelium of 57 of 62 patients is nearest the proliferating stem/TA population (cE03), and each tumour subtype is nearest its own normal counterpart (Tumor cE03 to cE03, r = 0.74). Restricting the normal arm accordingly: malignant versus normal stem/TA, higher in 35 of 39 patients (+14.8 points); subtype against the same subtype, +14.6 (cE01) and +15.4 (cE03). Normal crypt-base cells occupy the state at 7-9% against 0.1-1% for differentiated enterocytes, so the unrestricted comparison overstates the effect (+19.1 points) and the crypt-base comparison is reported as primary. *(iv) Decoy modules, and what the two robustness percentages are really comparing.* Adding five modules of random genes (44 to 49) requires scoring new modules gene by gene, which the reduced 23-cancer harness does and the full pipeline of Fig. S7-II does not, so no full-fidelity decoy arm exists. Within the reduced harness, where both arms are measured the same way, the decoy arm is not the harsher perturbation: the cytotrophoblast holds the top three by cross-cancer recurrence in 42 of 100 decoy draws against 96 of 200 draws that instead drop five real modules (48%; Fisher *P* = 0.39, rank distributions Mann-Whitney *P* = 0.49), and on the quantitative endpoints the decoy arm is marginally the better of the two (median core fraction 19.0% against 17.7%, enrichment 9.7x against 9.1x). The 42% should therefore not be read against the 84% of Fig. S7-II: those are two rankings of the same perturbation class, not two perturbations. The cytotrophoblast is rank 1 by cosine to the malignant attractor and rank 3 by cross-cancer recurrence when nothing is perturbed, so “stays in the top three” is nearly free in the first ranking and close to a coin-flip in the second: the same 200 five-module drops give 84% by attractor rank and 44% by recurrence rank. `scripts/E2d_module_robustness.py`, `scripts/E2e_module_robustness_fullfid.py`. *(v) Single genes that move the anchor.* Rebuilding the catalog 206 times, once per gene with that gene removed, leaves the cytotrophoblast at rank 1–2 in 202 of the 206 rebuilds (core occupancy 15.4–30.7%). Four genes dethrone it on their own: **EIF4EBP1** (mTOR; cytotrophoblast to rank 12, occupancy 15.4%), **FKBP5** (AR/GR/MR; rank 6), **SLC2A1** (HIF; rank 8) and **CD36** (PPARγ; rank 7). The pre-registered rule treats more than three such genes as a failure of this test, and we report it as one. Deepening every module to 10, 15 or 20 genes tells the same story from the other side: core occupancy stays at or above 11.7% at every depth, while the cytotrophoblast leaves the top three at every depth and under both expansion rules. The magnitude of the convergence is robust to how the modules are populated; the identity of the nearest cell rests on a four-to-five-gene sampling of each module and, within that, on a handful of individual genes. `scripts/E36_colorectal_comparator.py`, `perceptome_validation/scripts/T8_region_vs_label.py`, `perceptome_validation/results/T4b_gene_jackknife.csv`."+'\n\n'

# Note S-T is held in its own file: the pre-registered activity-scoring test (T7), whose table is
# a markdown table and whose verdict is a FAIL the pre-registration obliges us to print in full.
si=si+'\n\n```{=latex}\n\\clearpage\n```\n\n'+(MAN/'note_ST_activity.md').read_text().strip()+'\n\n'

# Note S-S is RETIRED (27.10.2026). It was the controls paragraph moved out of the Discussion in
# block 1, and by the time the control notes were written it restated all of them: every sentence
# of it is in the Discussion or in Notes S-A, S-B, S-C and S-H, which was checked phrase by phrase
# before dropping it. `manuscript/note_SS_controls.md` is kept on disk, marked retired, and is no
# longer injected or shipped. Its letter stays unused, like S-E's.

# ===== v10 control notes (E3 assignment/magnitude + E4 calibration) =====
try:
    _pa=pd.read_csv(ROOT/'results/referee2_newdata/E18_purity_audit.csv').fillna('')
    _pa=_pa[['cohort','selection','purity_class','malig_frac_of_arm']]
    _pa.columns=['cohort','malignant-cell selection','purity class','annotated-malignant % of arm']
    si=si+'\n\n**Table S-P. Per-cohort malignant-cell purity audit (`scripts/E18_purity_audit.py`).** Selection method and, where the source carries annotations, the annotated-malignant fraction of the arm actually used. Three cohorts are flagged at dilution risk — HCC (68% of the hepatocyte arm from tumor samples), BRCA (24% annotated cancer-epithelial), BCC (12%) — the same failure mode as the lung artifact; we tested directly whether contamination masks convergence, on the two of these three whose labels allow it, and it does not (Note S-I: cleaning HCC and BRCA lowered rather than raised them), so no conservativeness is claimed for this group. NSCLC filters to the census malignant label; the colorectal/gastric/lung *deepening* cohorts use author-annotated tumor compartments (Table above / Fig. S11-II). Seventeen solid/liquid cohorts use source-pipeline malignant files; three germ-cell/sarcoma cohorts are marker-labeled (overstates, bounded by the drop-4 robustness, 0.166→0.174). A CNV-inference floor across all 25 is the one remaining production step.\n\n'+_pa.to_markdown(index=False)+'\n\n'
except Exception as _ex: si=si+f'\n\n(Table S-P purity audit: {_ex})\n\n'

# ===== Note S-U + Table S-Q: what the Census cohorts actually are =====
# Every number in the prose below is formatted from results/census_cohorts.csv and
# results/census_dedup_effect.json, so the note cannot drift from the audit that produced it.
_cc = pd.read_csv(ROOT / 'results/census_cohorts.csv')
_dd = json.loads((ROOT / 'results/census_dedup_effect.json').read_text())
_W = ('zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen '
      'fifteen sixteen seventeen eighteen nineteen twenty').split()
_w = lambda n: _W[n] if n < len(_W) else str(n)
_AUTH = _cc[_cc.selection == 'author malignant annotation'].cohort.tolist()
_HEUR = _cc[_cc.selection != 'author malignant annotation']
_DUP = _cc[_cc.duplicated_rows.fillna(0) > 0].sort_values('duplicated_rows', ascending=False)
_dup_txt = '; '.join(
    f"{r.cohort.lower()} {int(r.duplicated_rows):,} of {int(r.rows):,} "
    f"({100 * r.duplicated_rows / r.rows:.0f}%)" if r.duplicated_rows >= 100 else
    f"{r.cohort.lower()} {int(r.duplicated_rows):,} cells"
    for r in _DUP.itertuples())
_pub, _ded = _dd['as published'], _dd['deduplicated']
# Whether the ranking moves is a fact about the artifact, not an expectation: the sentence is
# assembled from whatever census_dedup_effect.json says, so it cannot assert a null it did not get.
_rk = [('by cosine to the malignant attractor', 'cyto_attr'), ('by cross-cancer recurrence', 'cyto_rec')]
_same = [n for n, k in _rk if _pub[k] == _ded[k]]
_moved = [(n, k) for n, k in _rk if _pub[k] != _ded[k]]
_dedv = (f"Recomputing the panel baseline with duplicate rows dropped within each cohort before "
         f"quality control \u2014 everything else identical, including the seed and the same "
         f"{_pub['cohorts']} cohorts \u2014 moves the core fraction from {_pub['corefrac']}% to "
         f"{_ded['corefrac']}% and the anchor ")
if not _moved:
    _dedv += ("not at all: it is rank " + ' and rank '.join(
        f"{_pub[k]} {n}" for n, k in _rk) + " either way.")
else:
    _dedv += ('; '.join(
        f"from rank {_pub[k]} to rank {_ded[k]} {n}" for n, k in _moved)
        + (', while staying rank ' + ' and '.join(f"{_pub[k]} {n}" for n, k in _rk
                                                  if n in _same) if _same else '')
        + ' \u2014 that is, towards the paper\u2019s claim and not away from it'
        + ('; the attractor top three is the same list either way.'
           if _pub['attr_top3'] == _ded['attr_top3'] else
           ', and the attractor top three becomes ' + ', '.join(_ded['attr_top3'])
           + ' against ' + ', '.join(_pub['attr_top3']) + ' as published.'))
_dedv += (f" The cell count is unchanged at {_pub['cells']:,} because the "
          f"{_dd['harness']['cap_per_cohort']:,}-cell-per-cohort cap of this harness absorbs the "
          f"removal in every affected cohort, so what changes is which cells the cap draws from. "
          f"(That cap is the one Fig. S7-II\u2019s robustness harness uses, not the 4,000 of the "
          f"paper\u2019s own figures: the point here is that the two passes differ in nothing but "
          f"the duplicate rows.) The double counting is therefore a defect of the Census query and not a driver of "
          f"the result, and we state it here instead of leaving it implicit in a cell count "
          f"(`scripts/verify/census_dedup_effect.py`).")

si = si + '\n\n' + (
 "**Note S-U. What the Census cohorts of the panel resolve to (`scripts/verify/census_cohorts.py`, "
 "Table S-Q).** Table S1 names fourteen of the twenty-five panel cohorts only as “Census”, "
 "which is neither an identifier nor a selection rule. Reading each cohort back out of its stored "
 "object gives both, and four things follow that Table S1 does not say. "
 "*(i) For " + _w(len(_HEUR)) + " of the " + _w(len(_cc)) + ", “malignant” is a "
 "cell-type label and not a malignant annotation.* " + _w(len(_AUTH)).capitalize() + " cohorts are selected on "
 "the author or atlas annotation `malignant cell`/`neoplastic cell` (" + ', '.join(
    c.lower() for c in _AUTH) + "). The rest are selected on a label that is not about malignancy: a "
 "*normal* lineage label in ovarian serous (`fallopian tube secretory epithelial cell`), prostate "
 "(basal, luminal and unspecified prostate epithelium), neuroblastoma (`neuroblast (sensu "
 "Vertebrata)`) and B-ALL (B-lineage progenitors); a bare `epithelial cell` in a tumour sample in "
 "HNSCC, small-cell lung and gastric; and `abnormal cell` in RCC. Epithelial-in-tumour is the same "
 "rule that produced the lung dilution artifact retracted in Fig. S11-II, where a 9%-malignant arm "
 "read as normal. `E18_purity_audit.csv` records several of these as “pre-filtered … "
 "author-annotated (high)”, which holds for the first group only, and Table S1’s "
 "`malignant_label` for HNSCC said `author` and now says `Census`. The direction of the resulting "
 "bias is the one argued in Note S-I — normal-adjacent epithelium pulls a tumour arm toward its "
 "own normal baseline and so understates deepening — but the purity of these cohorts is a "
 "heuristic’s purity, and we say that here instead of claiming an annotation we do not have. "
 "*(ii) `is_primary_data` was not filtered, so some cells are counted more than once.* A Census "
 "query without that predicate returns a cell once per dataset it was deposited in. "
 + _w(len(_DUP)).capitalize() + " cohorts carry such rows — " + _dup_txt + " — and the prostate case "
 "is the sharpest: its three dataset identifiers are one atlas plus that atlas’s own Epithelial "
 "and Cycling subsets, so cycling prostate cells enter the centroid three times. " + _dedv + " "
 "*(iii) No paired cohort is outside the panel.* The panel’s gastric cohort is Census dataset "
 "`19053a82` with donors `Pat01-B`–`Pat24-B`: the twenty-three GSE150290 tumour samples, "
 "which are the twenty-three patients of the paired gastric test in §5, so that test is no "
 "longer described as held out (Note S3). Checking the other two the same way "
 "(`scripts/verify/paired_overlap.py`) gives the same answer: the panel’s colorectal cohort "
 "*is* Pelka GSE178341, so the flagship paired test runs on patients whose malignant cells are "
 "in the panel by construction; and the panel’s NSCLC cohort — its smallest, 241 malignant "
 "cells over 55 donors — carries nine `Kim_Lee_2020` donors, five of them (P0018, P0020, "
 "P0028, P0031, P0034) among the ten paired lung patients of E12. The paired design was never "
 "an out-of-sample test of *which* normal cell the panel names; what it cancels is lineage and "
 "batch, because both arms are one patient on one protocol, and §5 now says that instead. *(iv) Platform and tissue.* Table S1 "
 "records “10x” for every Census cohort: glioblastoma pools five datasets over nine assays "
 "including microwell-seq, Drop-seq, CEL-seq2 and STRT-seq, and melanoma includes 3,058 Smart-seq2 "
 "cells. Ovarian serous spans nineteen tissues including omentum, intestine, peritoneum and ascitic "
 "fluid — disseminated sites — where Table S1 says “ovary”; RCC pools two "
 "independent studies, one an ipilimumab/nivolumab-resistant series, in a single Table S1 row. "
 "Table S-Q carries the identifiers so that any of this can be checked against the Census snapshot "
 "the paper pins (2025-11-08).\n")

_q = _cc[['cohort', 'rows', 'unique_cells', 'primary_rows', 'donors', 'datasets',
          'malignant_label', 'dataset_ids']].copy()
_q.columns = ['cohort', 'rows', 'unique cells', 'primary rows', 'donors', 'datasets',
              'Census cell_type selected on', 'dataset ids (truncated)']
si = si + '\n\n**Table S-Q. The Census cohorts of the 25-cancer panel, resolved** ' \
    '(`scripts/verify/census_cohorts.py`; `results/census_cohorts.csv`). *rows* is the stored ' \
    'cohort before the paper’s quality control and the 4,000-cell cap; *unique cells* counts ' \
    'distinct `observation_joinid`, so rows minus unique cells is double counting; *primary rows* ' \
    'is how many rows carry `is_primary_data == True`. Blank cells are cohorts whose stored object ' \
    'does not carry the column. Three cohorts (NSCLC, small-cell lung, gastric) were assembled by ' \
    'an earlier pull that kept fewer Census columns.\n' \
    + csv_tex_df(_q, [0.13, 0.06, 0.07, 0.07, 0.05, 0.05, 0.30, 0.27], size='scriptsize')


if _MAIN_MISSES: raise SystemExit('main-figure anchors no longer match the manuscript: '+repr(_MAIN_MISSES))
for _f in ['v11_fig1.png','v11_fig2.png','v11_fig3.png','v11_fig4.png','v11_fig5.png',
           'v11_figS29_invitro.png','v11_fig6.png']:
    _c=main.count(_f)
    if _c!=1: raise SystemExit(f'{_f} injected {_c} times, expected exactly 1')
# --- three sources, two documents --------------------------------------------------------
# The co-author review of 2026-10-27 asked for the main text and the SI as separate files, with no
# table of contents. The COMBINED source is still written, because every verifier reads it
# (check_numbers, check_citations, numlib, make_release) and a claim must be checkable across both
# halves at once; it is no longer typeset.
si_body = si.split('# Supporting Information', 1)[1]
combined = main + '\n\n\\newpage\n\n# Supporting Information\n' + si_body
out = MAN / 'paper3_v10_full_src.md'; out.write_text(combined)
main_src = MAN / 'paper3_v12_MAIN.md'; main_src.write_text(main)
# standalone: the body opens with '## Supporting Information — <title>', so the separator
# heading the combined document needs would print the words twice on the SI's own first page
si_src = MAN / 'paper3_v12_SI.md'
si_src.write_text(si_body.lstrip().replace('## Supporting Information — ',
                                           '# Supporting Information — ', 1))
hdr = MAN / '_uni_header.tex'; hdr.write_text(UNI_HEADER)
print('combined source:', out, len(combined), 'chars')
print('main source:', main_src, len(main), 'chars | SI source:', si_src, len(si_body), 'chars')

# --- pandoc -> PDF via tectonic ---
def build_pdf(src, pdf):
    cmd = ['pandoc', str(src), '-s', '-o', str(pdf), '--pdf-engine=tectonic', '-H', str(hdr),
           '-V', 'geometry:margin=2.3cm', '-V', 'fontsize=10pt', '-V', 'linkcolor=blue',
           '-V', 'colorlinks=true', '-V', 'documentclass=article']
    print('RUN:', ' '.join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    print('STDOUT', r.stdout[-2000:]); print('STDERR', r.stderr[-3000:])
    print('PDF exists:', pdf.exists(), pdf.stat().st_size if pdf.exists() else 0)
    return pdf.exists()

ok_main = build_pdf(main_src, MAN / 'paper3_v12_MAIN.pdf')
ok_si = build_pdf(si_src, MAN / 'paper3_v12_SI.pdf')
if _SI_MISSES: print('!! SI ANCHORS NOT FOUND:', _SI_MISSES)
if _MAIN_MISSES: print('!! MAIN ANCHORS NOT FOUND:', _MAIN_MISSES)
if not (ok_main and ok_si): raise SystemExit('a document failed to build')
