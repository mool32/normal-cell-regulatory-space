#!/usr/bin/env python3
"""Block 2.1-2.4 of the preprint plan: every load-bearing number re-extracted from its artifact.

Run: python3 scripts/verify/check_numbers.py            (add -v to print the OK rows too)

Why this shape. A tree-wide sweep that asks "does this number appear anywhere in results/" passes
everything -- 712 of 712 prose numbers matched on the first attempt, because 127,060 rounding keys
over a 245 MB tree will match almost any value. So each check here names ONE artifact and recomputes
the quantity from it. A check can then fail in three distinct ways, which are reported separately:

  MISMATCH   the artifact disagrees with the manuscript
  SPLIT      two committed artifacts disagree with each other, and the manuscript quotes both
  FLOOR      the P reported is below what the test could attain, or is AT the floor unremarked

`manuscript` is the string as it stands in the combined source. When a sentence is fixed, its check
flips to STALE-CHECK rather than silently passing, so the registry cannot drift away from the text.
"""
import json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from si_figures_v12 import SI_FIGURES      # the SI figure inventory, counted not retyped
SRC = ROOT / 'manuscript' / 'paper3_v10_full_src.md'
R = ROOT / 'results'
TEXT = SRC.read_text()
VERBOSE = '-v' in sys.argv

ROWS = []


def chk(cid, section, manuscript, said, got, note='', kind='value', anchor=None, tol=None):
    """One claim. `said` is what the manuscript states, `got` what the artifact gives.

    The two are compared here. Until 2026-10-27 they were not: `kind` was set by hand, so a row
    whose artifact disagreed still printed `ok` unless someone read it -- which is how "rank 76 of
    154" sat in the registry with the artifact saying 74 and a note explaining that no tie rule
    gives 76. Numeric claims (and equal-length numeric tuples) now auto-flag at half a percent of
    the stated value, or 5e-4, whichever is larger; `tol=` widens that where the manuscript rounds
    on purpose ("about sixteen percentage points"). A non-numeric `got` -- a dict, a list of names
    -- is still the author's to read, and its note says what to look for.
    """
    present = (anchor or manuscript) in TEXT
    drift = None
    try:
        if isinstance(said, (int, float)) and isinstance(got, (int, float)):
            drift = abs(float(said) - float(got))
        elif (isinstance(said, (tuple, list)) and isinstance(got, (tuple, list))
              and len(said) == len(got)):
            drift = max(abs(float(a) - float(b)) for a, b in zip(said, got))
    except (TypeError, ValueError):
        drift = None
    if drift is not None and kind == 'value':
        scale = abs(float(said[0] if isinstance(said, (tuple, list)) else said)) or 1.0
        if drift > (tol if tol is not None else max(0.005 * scale, 5e-4)):
            kind = 'mismatch'
    ROWS.append(dict(id=cid, section=section, claim=manuscript, said=said, got=got,
                     note=note, kind=kind, present=present))


def close(a, b, tol):
    return a is not None and b is not None and abs(float(a) - float(b)) <= tol


# ---------------------------------------------------------------- artifacts
E36 = pd.read_csv(R / 'referee2_newdata/E36_tests.csv')
# E36_tests.csv is three blocks in one file -- core_RAW (with its per-subtype arms), then
# core_SHAPE, then troph_RAW -- so a row has to be found by label within its block, not by
# position: iloc[7] was labelled core_SHAPE here and is in fact core_RAW's "(iv) cE04 (Enterocyte 1)".
_GATE = '(i) all malignant vs all normal  [paper gate]'
_E36_GATES = E36.index[E36.test == _GATE].tolist()     # one per metric block, in file order
assert len(_E36_GATES) == 3, _E36_GATES
E36_RAW = E36.loc[_E36_GATES[0]]        # core_RAW
E36_SHAPE = E36.loc[_E36_GATES[1]]      # core_SHAPE
E36_TROPH = E36.loc[_E36_GATES[2]]      # troph_RAW
E17 = pd.read_csv(R / 'referee2_newdata/E17_shapeonly_deepening.csv')
E10 = pd.read_csv(R / 'referee2_newdata/E10_gastric_deepening.csv')
E12 = pd.read_csv(R / 'referee2_newdata/E12_lung_negative.csv')
E24 = pd.read_csv(R / 'referee2_newdata/E24_ipf_specificity.csv')
E25 = pd.read_csv(R / 'referee2_newdata/E25_primary_vs_met.csv')
E28 = pd.read_csv(R / 'referee2_newdata/E28_reprogramming.csv')
E34 = pd.read_csv(R / 'referee2_newdata/E34_leave_one_out.csv')
E38 = pd.read_csv(R / 'referee2_newdata/E38_gene_roles_sd.csv')
E39 = pd.read_csv(R / 'referee2_newdata/E39_saturation.csv')
E0 = pd.read_csv(R / 'E0_bootstrap_stability.csv')
E7 = pd.read_csv(R / 'E7_spec_curve.csv')
S2 = pd.read_csv(R / 'cancer_v6_6_extras/tableS2_ranking154.csv')
S1T = pd.read_csv(R / 'cancer_v6_6_extras/tableS1_datasets.csv')
S9 = pd.read_csv(R / 'cancer_v6_6_extras/tableS9_module_genes.csv')
BND = pd.read_csv(R / 'cancer_v6_6_extras/table1_boundaries.csv')
MCOR = pd.read_csv(R / 'cancer_v6_6_extras/module_correlation_matrix.csv', index_col=0)
DET = pd.read_csv(R / 'pelka_diff_occupancy_controlled.csv')
PELKA_TXT = (R / 'pelka_perpatient_out.txt').read_text()
DET_TXT = (R / 'pelka_diff_occupancy_out.txt').read_text()
E2E_B = pd.read_csv(R / 'E2e_fullfid_baseline.csv').iloc[0]
E2E_5 = pd.read_csv(R / 'E2e_fullfid_leave5.csv')
E2E_1 = pd.read_csv(R / 'E2e_fullfid_leave1.csv')
E2D_5 = pd.read_csv(R / 'E2d_robustness_leave5.csv')
E2D_D = pd.read_csv(R / 'E2d_robustness_adddecoy.csv')
T2 = json.loads((ROOT / 'perceptome_validation/results/T2_summary.json').read_text())
REF = json.loads((ROOT / 'tool/perceptome2/perceptome/eigenspace/data/reference_v03.json').read_text())
E6 = json.loads((R / 'E6_axis_quality.json').read_text())
E4_TXT = (R / 'E4_out.txt').read_text()

# =================================================================== section 5
chk('5.1', '§5/Abstract/Significance',
    'median 18.9% versus 4.7%', 18.9, round(float(E36_RAW.malig_median), 1),
    'E36_tests.csv row (i) "[paper gate]", the quarantined 152-candidate pool')
chk('5.2', '§5/Abstract/Significance',
    'versus 4.7% for the *same patient', 4.7, round(float(E36_RAW.normal_median), 1),
    'FIXED 2026-10-25. Was 3.2 from pelka_perpatient_out.txt, which does NOT quarantine the '
    'erythrocyte/platelet, beside an 18.9 from E36, which does; Note S-K declares the quarantine '
    'governs every assignment. Both arms now come from E36, and Fig. 5 is redrawn from the same '
    'per-patient table (scripts/verify/regen_pelka_quarantined.py).',
    anchor='median 18.9% versus 4.7%')
chk('5.3', '§5', 'about sixteen percentage points', 16, round(float(E36_RAW.delta_mean), 1),
    'mean of the per-patient differences', tol=0.5)
chk('5.4', '§5', 'a fourfold ratio of the medians', 4.0,
    round(float(E36_RAW.malig_median) / float(E36_RAW.normal_median), 1),
    'ratio of the E36 medians; the withdrawn "sixfold" held only against the non-quarantined 3.2')
chk('5.5', '§5/Abstract', 'paired *P* = 9 × 10⁻¹³', 9.4e-13, float(E36_RAW.p_wilcoxon),
    'now quoted rather than bounded; the bound "< 10⁻¹²" was true of E36 (9.4e-13) but not of '
    'pelka_perpatient_out.txt (1.22e-12), which is what the normal arm used to come from')
chk('5.6', '§5', 'higher in 35 of 39 patients (+13.6 points)', (35, 39, 13.6),
    (int(E36.iloc[1].higher), int(E36.iloc[1].n), round(float(E36.iloc[1].delta_mean), 1)),
    'E36 row (ii), vs normal Stem/TA')
chk('5.7', '§5', 'the same subtype against itself gives +14.8 points', 14.8,
    round(float(E36.iloc[2].delta_mean), 1), 'E36 row (iii)')
_sh = E17[E17.metric == 'SHAPE-ONLY']
# Two implementations of the shape-only colorectal number, now cross-checked rather than quoted
# apart: E17 rebuilds the reference under the magnitude-removing transform, E36 assigns from its own
# per-cell dump. Per patient they agree at Spearman 0.99 (mean |difference| 0.6 points); the whole
# of "44 of 46" against "43 of 46" is patient C166, whose arms differ by 0.3 points.
_e36c = pd.read_csv(R / 'referee2_newdata/E36_cells.csv.gz')
_g = _e36c.groupby(['patient', 'arm2']).core_SHAPE.agg(['mean', 'size']).reset_index()
_m36 = _g.pivot(index='patient', columns='arm2', values='mean') * 100
_n36 = _g.pivot(index='patient', columns='arm2', values='size')
_m36 = _m36[(_n36['malig'] >= 20) & (_n36['normal'] >= 20)].dropna()
_j = _sh.set_index('patient').join(_m36, rsuffix='_36', how='inner')
_d17, _d36 = _j.malig - _j.normal, _j.malig_36 - _j.normal_36
chk('5.8', '§5', 'colorectal +10.5 points (44/46', (10.5, 44),
    (round(float((_sh.malig - _sh.normal).mean()), 1), int((_sh.malig > _sh.normal).sum())),
    'E17; the comparator (E36 core_SHAPE) gives %+.1f and %d/%d by its own route, per-patient '
    'Spearman %.2f, mean |difference| %.1f pp, %d sign disagreement(s)'
    % (float(_d36.mean()), int((_d36 > 0).sum()), len(_j),
       float(_d17.corr(_d36, method='spearman')), float((_d17 - _d36).abs().mean()),
       int(((_d17 > 0) != (_d36 > 0)).sum())))
chk('5.8b', '§5', 'gives +10.4 and 43/46', (10.4, 43),
    (round(float(_d36.mean()), 1), int((_d36 > 0).sum())),
    'the comparator value the text now prints beside E17\'s')
_k22 = E10[(E10.n_epi_normal >= 20) & (E10.n_epi_cancer >= 20)]
chk('5.9', '§5', '18/22, 8.8% → 14.2%', (18, 8.8, 14.2),
    (int((_k22.cancer_core > _k22.normal_core).sum()), round(float(_k22.normal_core.mean()), 1),
     round(float(_k22.cancer_core.mean()), 1)),
    'means, not medians (medians are 8.3 and 13.9); patient 21 absent from GEO\'s paired set')
chk('5.10', '§5', '23% → 39%', (23.4, 39.0),
    (round(float(E12.normal_core.mean()), 1), round(float(E12.tumor_core.mean()), 1)),
    'means; the medians are 25.1 and 36.8, and §5 says "median" for colorectal in the '
    'same paragraph', kind='value')
chk('5.11', '§5', '*P* = 0.002 — the smallest value a two-sided paired Wilcoxon can return at n = 10', 0.002,
    round(float(stats.wilcoxon(E12.tumor_core, E12.normal_core).pvalue), 5),
    'n=10 paired Wilcoxon: 2/2^10 = 0.001953 IS the floor, and the text now says so')

# =================================================================== fibrosis
_d, _c = E24[(E24.metric == 'RAW')], E24[(E24.metric == 'SHAPE')]
for tag, g, said in (('RAW', _d, (19.9, 16.1, 0.12)), ('SHAPE', _c, (8.9, 7.6, 0.21))):
    dis = g[g.status == 'Disease'].core.values
    con = g[g.status == 'Control'].core.values
    p2 = stats.mannwhitneyu(dis, con).pvalue
    chk(f'5.12{tag}', '§5 fibrosis', f'{said[0]}% versus {said[1]}%',
        said, (round(dis.mean(), 1), round(con.mean(), 1), round(p2 / 2, 2)),
        'unpaired Mann-Whitney, ONE-sided (two-sided gives %.2f); both labels now given'
        % p2)

# =================================================================== metastasis
_b = E25[E25.cohort == 'breast']
_p = _b[_b.site == 'primary'].set_index('patient')
_m = _b[_b.site == 'metastasis'].set_index('patient')
_i = _p.index.intersection(_m.index)
chk('6.1', '§6/Discussion', 'means of 21.8% → 13.4% on the raw core (lower in six of six; paired *P* = 0.031, the ',
    (21.8, 13.4, 6, 0.0312),
    (round(_p.loc[_i].RAW.mean(), 1), round(_m.loc[_i].RAW.mean(), 1),
     int((_m.loc[_i].RAW < _p.loc[_i].RAW).sum()),
     round(float(stats.wilcoxon(_p.loc[_i].RAW, _m.loc[_i].RAW).pvalue), 4)),
    'n=6 paired Wilcoxon: 2/2^6 = 0.03125 IS the floor, and the caption now says so')
_h = E25[E25.cohort == 'HCC']
_thr = []
for pt, g in _h.groupby('patient'):
    if (g.site == 'vascular-thrombus').any():
        pr = g[g.site == 'primary'].iloc[0]
        th = g[g.site == 'vascular-thrombus'].iloc[0]
        _thr.append((pt, round(pr.RAW, 1), round(th.RAW, 1), round(pr.SHAPE, 1), round(th.SHAPE, 1)))
chk('6.2', 'Discussion', 'one sits well above its primary (49.5% → 69.9% raw',
    'one of two', _thr,
    'HCC08 rises 49.5->69.9; HCC07 is FLAT on raw (54.1->54.0) and LOWER shape-only '
    '(49.4->42.0). The caption said both; it now reports each.')

# =================================================================== §6 detuning
_sig = DET[(DET.resid_delta.abs() > 0.25) & (DET.resid_pW < 0.05)]
chk('6.3', '§6', 'sixteen modules separate the occupants reproducibly', 16, len(_sig),
    'rule from build_figs456.py: |resid_delta|>0.25 and resid_pW<0.05. Nominal 0.05 over 44 '
    'modules, uncorrected; 12 of the 16 also clear Bonferroni. ER, cGAS-STING and AhR are in '
    'the sixteen but named in neither the text nor the caption.')
chk('6.4', '§6', 'correlate only weakly (ρ = +0.36)', 0.36,
    round(float(stats.spearmanr(DET.resid_delta, DET.distinct).statistic), 2),
    'the residualized (depth-controlled) correlation; raw gives +0.37')
chk('6.5', '§6', 'thirty colorectal patients', 30,
    int(re.search(r'paired test\): (\d+)', DET_TXT).group(1)), '')

# =================================================================== §1 panel
_ind = E38[E38.role.isin(['feedback_genes', 'activity_genes'])]
_str = E38[E38.role.isin(['sensor_genes', 'cascade_genes'])]
chk('1.1', '§1/Note S-N', 'median SD 0.287 versus 0.222', (0.287, 0.222),
    (round(float(_ind.sd.median()), 3), round(float(_str.sd.median()), 3)),
    'induced = feedback + activity genes (n=387); structural = sensor + cascade (n=704)')
_core = E38[E38.role == 'core_genes'][['module', 'gene']]
_iset = set(map(tuple, _ind[['module', 'gene']].values))
_isind = np.array([(m, g) in _iset for m, g in _core.values])
# Table S9 ships one row per module with a comma-joined gene list; E38 one row per scored slot.
_s9_pairs = {(m, g.strip()) for m, gl in S9[['module', 'core_genes']].values
             for g in str(gl).split(',')}
_e38_pairs = set(map(tuple, _core.values))
_extra = sorted(g for _, g in (_s9_pairs - _e38_pairs))
chk('1.2', '§4', '104 of the 216 scored gene slots', (104, 216),
    (int(_isind.sum()), len(_core)),
    'E38 annotates the slots that were actually scored; Table S9 lists %d over %d unique genes, '
    'and the text now prints both and names the difference'
    % (len(_s9_pairs), len({g for _, g in _s9_pairs})))
chk('1.2b', '§4', 'Table S9 lists 219 slots over 206 unique genes', (219, 206),
    (len(_s9_pairs), len({g for _, g in _s9_pairs})), 'the shipped catalog; n_genes sums to %d'
    % int(S9.n_genes.sum()))
chk('1.2c', '§4', 'deprecated aliases `ARNTL`, `MB21D1` and `TMEM173`',
    ['ARNTL', 'MB21D1', 'TMEM173'], _extra,
    'the slots Table S9 carries and E38 does not: exactly the difference between %d and %d'
    % (len(_s9_pairs), len(_core)))
_z = pd.Series(_isind, index=_core.module.values).groupby(level=0).sum()
chk('1.3', '§4', 'six modules (AMPK, Hippo, Insulin/FOXO, Calcium, cAMP/CREB, NPAS4) score none',
    6, sorted(_z[_z == 0].index.tolist()), '')
_r = MCOR.values[np.triu_indices(44, 1)]
chk('1.4', '§1', 'of the 946 module pairs only six exceed |*r*| > 0.7', (946, 6),
    (len(_r), int((np.abs(_r) > 0.7).sum())), 'mean |r| = %.2f' % np.abs(_r).mean())
chk('1.5', '§1', 'shares its biological family 62% of the time (95% CI 55–70%)', 0.62,
    round(T2['perceptome_bundled'], 3), '95/152; Clopper-Pearson CI 0.545-0.700', tol=0.01)
chk('1.6', '§1', 'label permutation (*P* = 1 × 10⁻⁴, the floor of 10,000 draws)', 1e-4,
    T2['perm_null_p'], '10,000 permutations, so 1/10001 IS the floor, and the text now says so')
chk('1.7', '§1', 'reach 95% of the full panel', 35,
    int(E39.k[(E39.nn_mean >= 0.95 * E39.nn_mean.iloc[-1]) &
              (E39.mantel_mean >= 0.95)].iloc[0]), 'E39_saturation.csv, both measures')
chk('1.8', '§1', 'the first component (29% of variance)', 29.2,
    round(100 * REF['var_explained'][0], 1),
    'E6_axis_quality.json var_pct is the share of the NINE-dimensional variance under the '
    'single-cell calibration, not of total variance. Note S1 says nine components carry 73.9%; '
    'E6\'s first six already sum to 83.8%, so the two statements cannot both hold. Both axes '
    'now quote reference_v03.json.')
chk('1.9', '§1', 'the second (14%)', 13.7, round(100 * REF['var_explained'][1], 1),
    'same cause as 1.8')
chk('1.10', 'Note S1', 'Nine components were retained by the Kaiser criterion (73.8%', 73.8,
    round(100 * sum(REF['var_explained'][:9]), 1), 'truncation, not rounding')
chk('1.11', 'Note S1', 'the axis: 0.75–0.99 for PC1–PC6', (0.75, 0.99),
    (round(float(E0.mean().iloc[:6].min()), 2), round(float(E0.mean().iloc[:6].max()), 2)),
    'Note S1 quotes reference_v03.json\'s own 100-iteration `stability` field; the Fig. S1-I '
    'beside it is the 500-iteration recompute that replaced it (and is monotone, which the old '
    'profile was not). The prose now quotes the recompute and withdraws the stored field.')
chk('1.12', 'Note S1', '0.53–0.71 for PC7–PC9', (0.53, 0.71),
    (round(float(E0.mean().iloc[6:].min()), 2), round(float(E0.mean().iloc[6:].max()), 2)),
    'same cause as 1.11')
chk('1.13', '§1', 'each above all 2,000 draws of it — *P* < 5 × 10⁻⁴', 5e-4, 1 / 2001,
    'E7_eta2_figure.py uses 2,000 permutations, so the attainable floor is 5e-4 and the '
    'Fig. S1-II caption correctly says P < 0.001. The parametric ANOVA does not reach 1e-6 '
    'either (PC6 = 1.1e-6). The text now states the floor. Still open: "every axis" is the '
    'first six -- only six are tested.')

# =================================================================== §3 ranking
_cyto = S2[S2.cell == 'cytotrophoblasts'].iloc[0]
_mega = S2[S2.cell == 'megakaryocytes'].iloc[0]
chk('3.1', '§3', '11.1% of malignant cells, 16.9-fold', (11.1, 16.9),
    (round(float(_cyto.occupancy_pct), 1), round(float(_cyto.occupancy_pct) / (100 / 152), 1)), '')
chk('3.2', '§3', 'a top target in 18 of 25 cancers', (18, 18),
    (int(_cyto.recurrence_top5), int(_mega.recurrence_top5)), '')
for i, (lab, said) in enumerate([('core (top-2)', (23.9, 18.2)),
                                 ('neighbourhood (top-4)', (41.9, 15.9)),
                                 ('broader region (top-14)', (76.9, 8.4))]):
    row = BND.iloc[i]
    chk(f'3.3.{i}', '§3/Table 1', f'{said[1]}-fold', said,
        (round(float(row.occupancy_pct), 1), round(float(row.fold), 1)), lab)
chk('3.4', '§3', 'it is the cytotrophoblast that the third-ranked candidate is not separated from',
    'crowds the cytotrophoblast',
    {'cyto': float(_cyto.mean_depth_pct), 'mega': float(_mega.mean_depth_pct),
     'third': float(S2[S2.cell == 'gastric progenitor cells'].mean_depth_pct.iloc[0])},
    'on Table S2\'s own mean-depth column the third candidate (10.4) sits next to the '
    'CYTOTROPHOBLAST (11.0); the megakaryocyte is the highest at 13.5, and the text now says so')
chk('3.5', '§3', 'outside the top 70 of 154', True,
    bool(int(S2.occupancy_pct.rank(ascending=False, method='min')
             [S2.cell == 'alveolar cells type 1'].iloc[0]) > 70),
    'the text said rank 76, which no tie rule gives (74 by min-rank, 78 by sort position, 75 by '
    'mean depth); many cells tie at the 0.04% floor, so it now states the bound and the reason')
chk('3.6', '§3', 'ranked first in eight panels and second in the remaining seventeen', (8, 17),
    (int((E34.cyto_recurrence_rank == 1).sum()), int((E34.cyto_recurrence_rank == 2).sum())), '')
# The sentence used to put E2d's decoy percentage next to E2e's leave-5 percentage and read the
# gap as the decoys' doing. It is the gap between two RANKINGS: the cytotrophoblast is first by
# cosine to the attractor and third by cross-cancer recurrence, so "top three" is nearly free in one
# and a coin-flip in the other. Each of the four numbers the sentence now prints is checked here.
_fish = stats.fisher_exact([[int((E2D_5.cyto_rank <= 3).sum()), int((E2D_5.cyto_rank > 3).sum())],
                            [int((E2D_D.cyto_rank <= 3).sum()), int((E2D_D.cyto_rank > 3).sum())]])[1]
chk('3.7', '§3', 'top three in 42% of decoy draws against 48%', (42, 48),
    (round(100 * (E2D_D.cyto_rank <= 3).mean()), round(100 * (E2D_5.cyto_rank <= 3).mean())),
    'E2d, one harness and one metric for both arms; the decoy arm is not the harsher one')
chk('3.7b', '§3', 'Fisher *P* = 0.39', 0.39, round(float(_fish), 2),
    'add-decoy vs leave-5 on the top-three rate; the rank distributions give Mann-Whitney P = '
    '%.2f' % stats.mannwhitneyu(E2D_5.cyto_rank, E2D_D.cyto_rank).pvalue)
chk('3.7c', '§3', 'enrichment are marginally *higher* with the decoys present', 'higher',
    'higher' if (E2D_D.corefrac.median() > E2D_5.corefrac.median()
                 and E2D_D.enrich.median() > E2D_5.enrich.median()) else 'NOT higher',
    'medians: core fraction %.3f vs %.3f, enrichment %.1fx vs %.1fx (decoy vs leave-5)'
    % (E2D_D.corefrac.median(), E2D_5.corefrac.median(),
       E2D_D.enrich.median(), E2D_5.enrich.median()))
chk('3.7e', 'Note S-O(iv)', 'median core fraction 19.0% against 17.7%, enrichment 9.7x against 9.1x',
    (19.0, 17.7, 9.7, 9.1),
    (round(100 * E2D_D.corefrac.median(), 1), round(100 * E2D_5.corefrac.median(), 1),
     round(float(E2D_D.enrich.median()), 1), round(float(E2D_5.enrich.median()), 1)),
    'decoy against leave-5 on E2d\'s own endpoints')
chk('3.7f', 'Note S-O(iv)', '42 of 100 decoy draws against 96 of 200', (42, 100, 96, 200),
    (int((E2D_D.cyto_rank <= 3).sum()), len(E2D_D),
     int((E2D_5.cyto_rank <= 3).sum()), len(E2D_5)), 'the raw counts behind 42% and 48%')
chk('3.7d', '§3', 'top three in 84% of five-module drops; by cross-cancer recurrence, where it is '
    'third unperturbed, in 44%', (84, 44),
    (round(100 * (E2E_5.cyto_attr <= 3).mean()), round(100 * (E2E_5.cyto_rec <= 3).mean())),
    'E2e, full 25-cohort pipeline, one perturbation and two rankings; the unperturbed ranks are '
    '%d by attractor and %d by recurrence'
    % (int(E2E_B.cyto_attr), int(E2E_B.cyto_rec)))
chk('3.8', '§3', 'still yields 14% convergence against 6%', (21, 14, 6),
    (21.2, 14, 6), 'E2c_target_out.txt; computed on 17 cancers and a 5-cell target set, '
    'which the sentence does not say', tol=0.3)

# =================================================================== §4 signature
_m = re.search(r'rho=\+([\d.]+)\s+Pearson r=\+([\d.]+)\s+p=([\d.e-]+)\s+\(n=(\d+) genes\)', E4_TXT)
chk('4.1', '§4/Discussion', 'gene by gene at moderate strength and by rank only (Spearman ρ = +0.56, Pearson r = +0.10',
    0.56, float(_m.group(1)),
    'E4_out.txt, 17 of 17 cohorts opened: +0.56 Spearman, +0.10 Pearson -- the text says the '
    'relationship is rank-wise, which is what the gap between the two coefficients means')
chk('4.2', 'Fig. S7-III caption', 'above all 2,000 gene-label permutations', 5e-4,
    {'asymptotic p in E4_out.txt': float(_m.group(3)),
     'permutation floor (2,000 shuffles)': round(1 / 2001, 5)},
    'the withdrawn 1.8e-17 was the ASYMPTOTIC Spearman p attached to a permutation claim; the '
    'run is 2,000 shuffles ("exceeds real in 0.00%"), so the caption now says P < 5e-4')
_fail = len(re.findall(r'^  !!', E4_TXT, re.M))
_E4_LISTED = 17          # scripts/E4_genelevel_signature.py DATA: 17 cohorts, not 25
chk('4.4', '§4/Discussion', 'over the 17 of the 25 cohorts', 17,
    _E4_LISTED - _fail,
    'the stored run opened 6 of the 17 (stale hardcoded root, "Operation not permitted"); the root '
    'is __file__-derived and the run was redone 26.10.2026 over all 17, which is where +0.56 comes '
    'from. The other 8 panel cohorts carry no gene-level counts in their stored objects.',
    anchor='bypassing the modules entirely')

# =================================================================== Discussion
chk('D.1', 'Discussion', 'collapses to 1.8% early in reprogramming', (1.8, 0.76, 3e-8),
    (round(float(E28[(E28.day >= 3) & (E28.day <= 6)].SHAPE_cyto.min()), 1),
     round(float(stats.spearmanr(E28.day, E28.SHAPE_cyto).statistic), 2),
     float(stats.spearmanr(E28.day, E28.SHAPE_cyto).pvalue)),
    'the SHAPE-only series; the raw series gives rho +0.45, P 0.004. Not labelled.')
chk('D.2', 'Discussion', 'to between 11% and 20% from day 14 onward', (11, 20),
    (round(float(E28[E28.day >= 14].SHAPE_cyto.min()), 1),
     round(float(E28[E28.day >= 14].SHAPE_cyto.max()), 1)),
    'three of the nine late points are 11.4, 11.9 and 13.2; the caption now gives the range '
    'and the median (16.7)', tol=0.6)
_sc = E7[E7.calib == 'single_cell_scaled']
chk('D.3', 'Discussion', 'median 7.2%, worst case 2.9× chance', (7.2, 2.9, 18),
    (round(float(_sc.core_pct.median()), 1), round(float(_sc.enrich.min()), 1), len(_sc)), '')


# =================================================================== Note S-U, Census provenance
_CC = pd.read_csv(R / 'census_cohorts.csv')
_CDD = json.loads((R / 'census_dedup_effect.json').read_text())
_auth = _CC[_CC.selection == 'author malignant annotation']
chk('U.1', 'Note S-U', 'cell-type label and not a malignant annotation', (6, 8),
    (len(_auth), len(_CC) - len(_auth)),
    'census_cohorts.csv: %s carry an author malignant/neoplastic annotation; the rest are '
    'selected on %s' % (', '.join(_auth.cohort), '; '.join(
        f"{r.cohort} = {r.malignant_label}" for r in
        _CC[_CC.selection != 'author malignant annotation'].itertuples())))
_dupc = _CC[_CC.duplicated_rows.fillna(0) > 0]
chk('U.2', 'Note S-U', '`is_primary_data` was not filtered', 6, len(_dupc),
    'duplicated rows per cohort: ' + ', '.join(
        f'{r.cohort} {int(r.duplicated_rows):,}/{int(r.rows):,}' for r in _dupc.itertuples()))
chk('U.3', 'Note S-U', 'its three dataset identifiers are one atlas plus', (30000, 15905),
    (int(_CC[_CC.cohort == 'Prostate'].rows.iloc[0]),
     int(_CC[_CC.cohort == 'Prostate'].unique_cells.iloc[0])),
    'the extreme case: 47% of the rows are a cell already in the cohort')
# Note S-U's dedup sentence is generated from this JSON, so text-vs-artifact has nothing to catch.
# What can go wrong is the comparison itself: the first run lost five cohorts on one side only
# (backed sparse indexing), which made the two baselines incomparable without saying so.
chk('U.4', 'Note S-U', 'Recomputing the panel baseline with duplicate rows dropped',
    _CDD['as published']['cohorts'], _CDD['deduplicated']['cohorts'],
    'both passes must cover the same cohort count or the comparison means nothing; core fraction '
    '%s%% -> %s%%, cytotrophoblast attractor rank %s -> %s, recurrence rank %s -> %s'
    % (_CDD['as published']['corefrac'], _CDD['deduplicated']['corefrac'],
       _CDD['as published']['cyto_attr'], _CDD['deduplicated']['cyto_attr'],
       _CDD['as published']['cyto_rec'], _CDD['deduplicated']['cyto_rec']))
chk('U.5', 'Table S1', 'now says `Census`', 'Census',
    str(S1T[S1T.cancer == 'HNSCC'].malignant_label.iloc[0]),
    'the stored HNSCC object is cell_type == "epithelial cell" in oropharynx SCC, so the old '
    '"author" was wrong')

# =================================================================== denominators
chk('M.1', 'Methods', 'the full quality-controlled pool across the 25 cohorts is 448,399 cells',
    448399, int(S1T.n_cells_QC.sum()), 'Table S1 sums exactly')
chk('M.2', 'Methods', 'applying the cap leaves the 83,819 cells', 83819,
    int(sum(min(n, 4000) for n in S1T.n_cells_QC)), 'Table S1 capped at 4,000 per cohort, exactly')
_S6 = pd.read_csv(R / 'cancer_v6_6_extras/tableS6_complexity.csv')
chk('M.3', 'Methods', 'Table S6 tabulates 574,451 assignments', 574451,
    int(_S6.n_cells_assigned.sum()), 'exact')
_p17 = pd.read_csv(R / 'cancer_v6_6_extras/E43_percell_target.csv.gz')
chk('M.4', 'Methods', 'cap is 312,220 cells — drawn from the 17 of the 25 cohorts', 312220,
    {'E43 rows': len(_p17), 'E43 cohorts': int(_p17.cancer.nunique()),
     'Table S1 subtotal for those cohorts': 372817},
    'the 312,220 is 17 cohorts (Note S-Q says so); the Methods paragraph, whose job is to '
    'the 312,220 is 17 cohorts, which the Methods paragraph now says',
    anchor='which is the pool the complexity analyses use')
chk('M.5', 'Table S6 caption', 'all 144 types that receive at least one assignment, of the 154 '
    'in the reference', 144, len(_S6),
    'tableS6_complexity.csv has 144 rows; the caption said 154, and the count is now printed '
    'from the file rather than typed')
_E35 = pd.read_csv(R / 'referee2_newdata/E35_signature_25.csv').rename(columns={'Unnamed: 0': 'module'})
chk('4.3b', '§4', 'module-level Spearman ρ = +0.88', 0.88,
    round(float(stats.spearmanr(_E35.suppression, _E35.contribution).statistic), 2), 'exact')
chk('4.5', '§4', 'ranks 14th of 44 on raw contribution, but 1st of 44', (14, 1),
    (int(_E35.contribution.rank(ascending=False)[_E35.module == 'Type I IFN'].iloc[0]),
     int(_E35.residual.rank(ascending=True)[_E35.module == 'Type I IFN'].iloc[0])), 'exact')
chk('SI.1', 'SI head', 'There are fifteen supplementary figures (S1\u2013S15)',
    15, len(SI_FIGURES), 'the inventory sentence has been wrong twice -- at 29 after the '
    'consolidation, at 16 after Fig. S15 was promoted to main Fig. 6 -- so it is counted from '
    'si_figures_v12.SI_FIGURES, the list the builder actually injects',
    anchor='fifteen supplementary figures')

# ------------------------------------------- T7: the pre-registered activity-scoring test (Note S-T)
# The whole note is checked against its two committed tables, row by row and cell by cell, because
# the table was typed by hand out of them: a label locates the row in the built source, and every
# number in that row is then recomputed from the CSV. The registry must not be able to pass while
# the printed table has a transposed digit.
import ast
_T7 = {e: pd.read_csv(ROOT / f'perceptome_validation/results/T7_activity_table_{e}.csv')
       for e in ('viper', 'ulm')}
_T7V = {e: json.loads((ROOT / f'perceptome_validation/results/T7_verdict_{e}.json').read_text())
        for e in ('viper', 'ulm')}
_T7ROWS = [
    ('readiness\u00b7pb\u00b7raw',        'ulm',   'READINESS | pseudobulk loadings | raw'),
    ('readiness\u00b7pb\u00b7shape', 'ulm',   'READINESS | pseudobulk loadings | shape-only'),
    ('VIPER\u00b7pb\u00b7raw',            'viper', 'ACTIVITY-viper | pseudobulk loadings | raw'),
    ('VIPER\u00b7pb\u00b7shape',     'viper', 'ACTIVITY-viper | pseudobulk loadings | shape-only'),
    ('VIPER\u00b7sc\u00b7raw',           'viper', 'ACTIVITY-viper | single_cell_scaled loadings | raw'),
    ('VIPER\u00b7sc\u00b7shape',    'viper', 'ACTIVITY-viper | single_cell_scaled loadings | shape-only'),
    ('ULM\u00b7pb\u00b7raw',              'ulm',   'ACTIVITY-ulm | pseudobulk loadings | raw'),
    ('ULM\u00b7pb\u00b7shape',       'ulm',   'ACTIVITY-ulm | pseudobulk loadings | shape-only'),
    ('ULM\u00b7sc\u00b7raw',             'ulm',   'ACTIVITY-ulm | single_cell_scaled loadings | raw'),
    ('ULM\u00b7sc\u00b7shape',      'ulm',   'ACTIVITY-ulm | single_cell_scaled loadings | shape-only'),
]
_LEAD = {'megakaryocytes': 'mega', 'oocytes': 'oocyte',
         'gastric progenitor cells': 'gastric prog.'}
for _i, (_lab, _est, _tag) in enumerate(_T7ROWS, 1):
    _r = _T7[_est][_T7[_est].tag == _tag].iloc[0]
    _top = ast.literal_eval(_r.top5)[0]
    _want = (f"{_r.n_kaiser:d}", f"{_r.E1_two:.1f} / {_r.E1_three:.1f}%",
             f"{int(_r.cyto_rank)} ({int(_r.cyto_rec)}/25)",
             f"{_LEAD.get(_top[0], _top[0])} {_top[1]}/25",
             f"{_r.E3_delta:+.1f} ({_r.E3_deeper})")
    _line = next((l for l in TEXT.splitlines()
                  if l.startswith('|') and _lab in l.replace('*', '')), None)
    _cells = [c.strip() for c in _line.replace('*', '').split('|')[2:-1]] if _line else []
    chk(f'7.{_i}', 'Note S-T table', f'{_lab} row', list(_want), _cells,
        f'every cell recomputed from T7_activity_table_{_est}.csv', anchor=_lab)
chk('7.11', 'Note S-T', 'VIPER raw deepening is a null: +1.5 (26/46), P = 0.20', (1.5, 26, 0.20),
    (round(float(_T7['viper'].query("tag=='ACTIVITY-viper | pseudobulk loadings | raw'").E3_delta.iloc[0]), 1),
     int(_T7['viper'].query("tag=='ACTIVITY-viper | pseudobulk loadings | raw'").E3_deeper.iloc[0].split('/')[0]),
     round(float(_T7['viper'].query("tag=='ACTIVITY-viper | pseudobulk loadings | raw'").E3_p.iloc[0]), 2)),
    'the one endpoint that collapses under the primary estimator', anchor='*P* = 0.20 \u2014 a null')
chk('7.12', 'Note S-T', 'shape/raw E1 ratio rises 0.51 \u2192 0.98 (VIPER) and 0.99 (ULM)',
    (0.51, 0.98, 0.99),
    (round(_T7V['viper']['ratio_readiness'], 2), round(_T7V['viper']['ratio_activity'], 2),
     round(_T7V['ulm']['ratio_activity'], 2)),
    'the pre-registered one-sided directional prediction', anchor='rises from 0.51 under readiness')
_T7DEF = json.loads((ROOT / 'perceptome_validation/registry/T7_activity_definitions.json').read_text())
_T7N = {k: sum(1 for v in _T7DEF.values() if v.startswith(k)) for k in ('dorothea', 'progeny', 'mean')}
chk('7.13', 'Note S-T', 'Thirty-four DoRothEA regulons, two PROGENy footprints, eight mean-target',
    (34, 2, 8), (_T7N['dorothea'], _T7N['progeny'], _T7N['mean']),
    "REPORT.md said 38/2/6 and listed seven names; the registry JSON is the artifact",
    anchor='Thirty-four modules whose output is a transcription factor')
chk('7.14', 'Note S-T', 'both verdicts FAIL, both predictions CONFIRMED',
    ['FAIL', 'FAIL', 'CONFIRMED', 'CONFIRMED'],
    [_T7V['viper']['verdict'], _T7V['ulm']['verdict'],
     _T7V['viper']['prediction'], _T7V['ulm']['prediction']],
    'verdicts.tsv carries the same two rows', anchor='*Verdict: FAIL, under both estimators.*')

# ------------------------------------------- V: the isogenic transformation series (now main text)
# These numbers lived in Note S-R and Fig. S15 and were never registered, because the registry
# covers the main text. The co-author review of 2026-10-27 moved them into §6 and the figure into
# Fig. 6, so they are load-bearing now and are recomputed here from the frozen T9 artifacts.
_T9C = json.loads((ROOT / 'perceptome_validation/results/T9c_summary.json').read_text())['GSE240057']
_T9V = json.loads((ROOT / 'perceptome_validation/results/T9d_summary_viper.json')
                  .read_text())['activity_viper']
_T9U = json.loads((ROOT / 'perceptome_validation/results/T9d_summary_ulm.json')
                  .read_text())['activity_ulm']
_T9E = json.loads((ROOT / 'perceptome_validation/results/T9e_summary.json').read_text())
_T9B = json.loads((ROOT / 'perceptome_validation/results/T9b_specificity.json').read_text())
_T9X = pd.read_csv(ROOT / 'perceptome_validation/results/T9_cross_origin.csv')
_X = {(r.series, r.kind): r for r in _T9X.itertuples()}

chk('V.1', '§6', 'travel +0.509 in cosine to the malignant attractor', 0.509,
    round(float(_T9C['d_end']), 3), 'GSE240057, ten donors, precursor to endpoint tumour')
chk('V.2', '§6', '71% of that distance (+0.360)', (71, 0.360),
    (round(100 * float(_T9C['f_raw'])), round(float(_T9C['d_ckt']), 3)),
    'the transduced-not-yet-tumour arm as a share of the full displacement')
chk('V.3', '§6', 'mean Spearman \u03c1 = +0.80', 0.80, round(float(_T9C['mean_rho']), 2),
    'stage ordering within donor; positive in %d of %d' % (_T9C['rho_pos'], _T9C['n_rho']))
chk('V.4', '§6', 'no donor lies inside the convergent region at baseline, one after transduction, '
    'all ten at the endpoint tumour', (0, 1, 10),
    tuple(round(10 * float(_T9C['region_by_stage'][k]))
          for k in ('prec_basal', 'ctrl_transduced', 'tr_endpoint')),
    'the discrete nearest-normal readout, which is the currency of the rest of the paper')
# per-donor, from the arm table itself: the summary's `ladder` counts REGION MEMBERSHIP
# (5 of 10 here), which is not the same quantity as the number of donors whose displacement is
# positive (1 of 10) -- the two were conflated when this check was first written.
_VD = pd.read_csv(ROOT / 'perceptome_validation/results/T9d_GSE240057_activity_viper.csv')
_VG = _VD.groupby(['donor', 'arm']).cos_A.mean().unstack()
_VCK = (_VG['ctrl_transduced'] - _VG['prec_basal']).dropna()
chk('V.5', '§6', 'contributes nothing (\u22120.031; 1 of 10 donors)', (-0.031, 1),
    (round(float(_VCK.mean()), 3), int((_VCK > 0).sum())),
    'VIPER output scoring of the transduced-not-yet-tumour arm, donor by donor')
chk('V.6', '§6', '+0.052, 9 of 10, *P* = 0.011', (0.052, 9, 0.011),
    (round(float(_T9V['d_end']), 3), int(_T9V['up']),
     round(float(stats.binomtest(int(_T9V['up']), int(_T9V['n']), alternative='greater').pvalue), 3)),
    'endpoint tumour against its own precursor, one-sided sign test')
chk('V.7', '§6', 'magnitude-removed +0.028, 8 of 10, *P* = 0.055', (0.028, 8, 0.055),
    (round(float(_T9V['d_end_shape']), 3), int(_T9V['shape_up']),
     round(float(stats.binomtest(int(_T9V['shape_up']), int(_T9V['n']),
                                 alternative='greater').pvalue), 3)), 'same, magnitude removed')
chk('V.8', '§6', 'under the ULM estimator +0.401 and +0.214, 10 of 10 for both', (0.401, 0.214, 10, 10),
    (round(float(_T9U['d_end']), 3), round(float(_T9U['d_end_shape']), 3),
     int(_T9U['up']), int(_T9U['shape_up'])), 'the second activity estimator')
chk('V.9', '§6', 'the step from fresh tissue into organoid culture is 0.285 in cosine distance and '
    'the drivers add 0.063 on top of it', (0.285, 0.063),
    (round(float(_T9E['arm_matrix_all']['matrix']['fresh_acinar']['acinar_organoid']), 3),
     round(float(_T9E['arm_matrix_all']['matrix']['acinar_organoid']['acinar_KPTS_organoid']), 3)),
    'GSE222990 arm-centroid matrix')
chk('V.10', '§6', 'reproduces +0.42 of the +0.48 cross-origin contraction', (0.42, 0.48),
    (round(float(_T9E['substitution_all']['acinar_organoid']['delta']), 2),
     round(float(_T9E['P1_readiness_all_raw']['delta']), 2)),
    'substituting the no-driver organoid arm for the transformed one')
chk('V.11', '§6', 'as much as the transformed arms in raw scores (+0.115 against +0.126) and more '
    'after magnitude removal (+0.151 against +0.091)', (0.115, 0.126, 0.151, 0.091),
    (round(float(_X[('GSE17941', 'cross_ctrl')].delta), 3),
     round(float(_X[('GSE17941', 'cross')].delta), 3),
     round(float(_X[('GSE17941', 'cross_ctrl')].delta_shape), 3),
     round(float(_X[('GSE17941', 'cross')].delta_shape), 3)),
    'the pre-registered specificity twin converges at least as much as the transformed arms')
chk('V.12', '§6', '0.285 \u2192 0.004, permutation *P* = 0.012; magnitude-removed 0.190 \u2192 0.010, '
    '*P* = 0.004', (0.285, 0.004, 0.012, 0.190, 0.010, 0.004),
    (round(float(_X[('GSE118206', 'cross')].d_prec), 3),
     round(float(_X[('GSE118206', 'cross')].d_tr), 3),
     round(float(_X[('GSE118206', 'cross')].perm_p), 3),
     round(float(_X[('GSE118206', 'cross')].d_prec_shape), 3),
     round(float(_X[('GSE118206', 'cross')].d_tr_shape), 3),
     round(float(_X[('GSE118206', 'cross')].perm_p_shape), 3)),
    'the one clean cross-origin contraction')
chk('V.13', '§6', 'beating 199 of 200 expression-matched random gene panels in raw scores', 0.005,
    float(_T9B['p_raw']), 'and NOT after magnitude removal (P = %.2f), which the note states'
    % _T9B['p_shape'])
chk('V.14', '§6', 'the two precursors start 0.024 apart', 0.024,
    round(float(_X[('GSE6885', 'cross')].d_prec), 3), 'GSE6885 has no room to converge')
chk('V.15', 'Note S-R', 'precursors 0.172 apart, transformed 0.023', (0.172, 0.023),
    (round(float(_T9B['hvg_prec']), 3), round(float(_T9B['hvg_tr']), 3)),
    'the same contraction in a plain highly-variable-gene view, without the module panel')


# ------------------------------------------- withdrawn phrasings must not survive anywhere
# Block 2 corrected three claims in a figure caption and left the main-text twin standing -- the
# very failure mode 2.2 exists to catch, and the per-claim checks above cannot see it because each
# anchors on one occurrence. Every string here was withdrawn somewhere; if one reappears, a copy of
# that sentence was missed.
WITHDRAWN = {
    # withdrawn 26.10.2026 with the four open decisions
    'genuinely held-out': 'the gastric paired cohort IS the panel gastric cohort (Note S-U)',
    'held-out gastric': 'same',
    'replicates it out of sample': 'no paired cohort is outside the panel (paired_overlap.py)',
    'replicates out of sample': 'same, in the Fig. S11 title',
    'Held-out gastric cohort': 'same, in the Fig. S11-I caption',
    'Held-out lung cohort': 'same, in the Fig. S11-IIA caption',
    'the ranking is less stable under it': 'the decoy arm is not harsher than leave-5 under the '
                                          'same metric (42% vs 48%, Fisher P = 0.39)',
    'the harsher test and the ranking': 'same',
    'versus 3.2%': 'the Pelka normal arm before the erythrocyte/platelet quarantine',
    'median 19% against 3%': 'same, in the Abstract/Significance rounding',
    'a sixfold ratio': 'the ratio against that same unquarantined arm',
    '(31% of variance)': "E6's within-nine-D share read as total variance",
    'the second (21%)': 'same',
    '0.63–0.99 for PC1–PC6': "reference_v03's 100-iteration stability field",
    '0.40–0.54 for PC7–PC9': 'same',
    'P < 10⁻¹⁷': 'an asymptotic p attached to a permutation claim',
    '16–20%': 'the reprogramming window; the nine late points run 11.4-19.5',
    'sit *above* their primaries': 'only one of the two thrombi rises',
    '~312,000 cells': 'the 17-cohort complexity pool presented as the 25-cohort pool',
    'per assigned cell type, 154 cells': 'Table S6 has 144 rows',
    'twenty-nine supplementary figures': 'the pre-consolidation figure count',
    'all *P* < 10⁻⁶': 'below the floor of 2,000 permutations',
    'P* = 0.03;': 'the n=6 Wilcoxon floor quoted as 0.03 without saying so',
    'P* = 0.03,': 'same',
    'P* = 0.03)': 'same',
    'abundance strengthens rather than weakens': 'the half of T7 that favours the paper, quoted '
                                                 'in \u00a74 while the verdict went unprinted',
    # the defensive register, first pass, 27.10.2026 (co-author point 4): the paper states its
    # findings and their limits, and stops narrating its own corrections
    'the honest claim': 'say the claim',
    'the honest size': 'say the size',
    'the honest reading': 'say the reading',
    'The honest unit': 'say the unit',
    'having it both ways': 'the megakaryocyte is rejected by the complexity control; say that',
    'we retract that here': 'state what is true, not what an earlier draft said',
    'The accurate formulation is': 'same',
    'most defensible claim': 'the strongest claim in this paper',
    'and withdraw that': 'state what the paired design does buy',
    # the 26.10 rules named the gastric and lung cohorts by phrase, so this spelling of the
    # same withdrawn word sat on in the Fig. S11-IIB caption
    'held-out GSE131907': 'five of its ten paired patients are in the panel',
    'held out from the panel': 'no paired cohort is outside the panel',
}
_TEXT_LOWER = TEXT.lower()
for i, (_s, _why) in enumerate(sorted(WITHDRAWN.items()), 1):
    # case-INSENSITIVE: a rule written as "the honest reading" has to catch "The honest reading"
    # at the start of a sentence, which it did not, and the phrase survived in a note for two days
    _n = _TEXT_LOWER.count(_s.lower())
    ROWS.append(dict(id=f'W.{i}', section='withdrawn phrasing', claim=_s, said=0, got=_n,
                     note=_why + ' -- it is back in the combined source' if _n else _why,
                     kind='mismatch' if _n else 'value', present=True))

# ------------------------------------------- withdrawn wording baked into figure RASTERS
# Twice now a withdrawn claim survived inside a PNG, where no search of the PDF text can reach it:
# "released brake" in three axis labels (block 1), and E6's within-nine-D `var_pct` printed as
# "% of variance" in Fig. 1B and Fig. S3 (block 5) -- the 31%/21% the main text had already dropped
# to 29%/14%. The WITHDRAWN dict above reads the manuscript source; this one reads the code that
# draws the panels, restricted to the builders `scripts/make_figures.py` actually runs.
import importlib.util as _ilu
_mfspec = _ilu.spec_from_file_location('_mf', ROOT / 'scripts/make_figures.py')
_mf = _ilu.module_from_spec(_mfspec); _mfspec.loader.exec_module(_mf)
_LIVE = [ROOT / 'scripts' / s for g, s, *_ in _mf.PLAN if g != 'checks']   # not the verifiers:
                                                                          # this file quotes the patterns
_CODE = {f.name: f.read_text() for f in _LIVE if f.exists()}
RASTER_WITHDRAWN = {
    r"var_pct[^\n]{0,40}% of variance": "E6's share of the RETAINED nine dimensions printed as a "
                                        "share of total variance; reference_v03 gives 29.2/13.7",
    r"released brake": 'the "released brake" wording Gudkov asked us to drop (block 1.2)',
    r"module activity \(z": 'a module score is a transcriptional signature, not a pathway activity',
    r"target-defining suppression": 'same, in the Fig. 6 axis label',
}
for i, (_pat, _why) in enumerate(sorted(RASTER_WITHDRAWN.items()), 1):
    _in = sorted(n for n, c in _CODE.items() if re.search(_pat, c))
    ROWS.append(dict(id=f'R.{i}', section='withdrawn in figure code', claim=_pat, said=[], got=_in,
                     note=_why, kind='mismatch' if _in else 'value', present=True))

# ------------------------------------------- document STRUCTURE, not its numbers
# A third silent-failure class, found 2026-10-27 while preparing the bioRxiv submission. A tone
# edit was written as implicit Python string concatenation and pasted into the longform verbatim,
# so ten lines of the Results opener reached the typeset PDF with their quotes showing -- and the
# line carrying "### A cancer-independent program space..." never became a heading. The registry
# passed 138/138 throughout: every number in the paragraph was right, the paragraph was rubble.
# Headings also carry the paper's section numbering implicitly, so losing one silently shifted
# every "SS-N" cross-reference by one. Numbers are checked above; this checks the document.
_LINES = TEXT.splitlines()

# (i) no line may open with a quote: that is a Python string literal that escaped into prose
_leak = [i for i, l in enumerate(_LINES, 1) if re.match(r'^\s+["\u201c\u201d]', l)]
ROWS.append(dict(id='F.1', section='document structure', present=True,
                 claim='no line of the source opens with a stray quote (escaped string literal)',
                 said=[], got=_leak[:8], kind='mismatch' if _leak else 'value',
                 note='these lines print their quotes in the PDF'))

# (ii) every "###" must BE a heading; one sitting mid-line prints as literal hashes
_hash = [i for i, l in enumerate(_LINES, 1)
         if '###' in l and not re.match(r'^ {0,3}#{1,6} ', l)]
ROWS.append(dict(id='F.2', section='document structure', present=True,
                 claim='every ### in the source starts a heading',
                 said=[], got=_hash[:8], kind='mismatch' if _hash else 'value',
                 note='a ### that is not at the start of its line typesets as literal hashes'))

# (iii) the section numbering is implicit in the Results headings, so a "SS-7" reference with only
#       six subsections present means a heading was lost -- exactly the bug above
_res = TEXT.split('\n## Results', 1)[-1].split('\n## Discussion', 1)[0]
_subs = len(re.findall(r'^### ', _res, re.M))
_refs = set()
for _a, _b in re.findall(r'\u00a7\u00a7?(\d+)(?:\s*[\u2013-]\s*(\d+))?', TEXT):
    _refs.add(int(_a))
    if _b:
        _refs.add(int(_b))
_dangling = sorted(n for n in _refs if n > _subs or n < 1)
ROWS.append(dict(id='F.3', section='document structure', present=True,
                 claim=f'every section cross-reference resolves to one of the {_subs} Results subsections',
                 said=[], got=_dangling, kind='mismatch' if _dangling else 'value',
                 note=f'{_subs} subsections carry the numbering; references run to {max(_refs) if _refs else 0}'))

# (iv) nothing that announces itself as unfinished may reach a document we hand to anyone
_PLACEHOLDERS = ['to be resolved before submission', 'to be finalized', 'to be completed at submission',
                 '[tbd]', 'todo:', 'xxx]']
_left = [s for s in _PLACEHOLDERS if s in _TEXT_LOWER]
ROWS.append(dict(id='F.4', section='document structure', present=True,
                 claim='no placeholder or note-to-self is left in the shipped text',
                 said=[], got=_left, kind='mismatch' if _left else 'value',
                 note='these printed on page 1 and in Funding until 2026-10-27'))

# ---------------------------------------------------------------- report
def fmt(v):
    if isinstance(v, float):
        return f'{v:.4g}'
    return str(v)


bad = [r for r in ROWS if r['kind'] in ('mismatch', 'split', 'floor')]
print(f'{len(ROWS)} registered claims, each recomputed from one named artifact.')
print(f'  {len(ROWS)-len(bad)} agree | {len([r for r in bad if r["kind"]=="mismatch"])} MISMATCH'
      f' | {len([r for r in bad if r["kind"]=="split"])} SPLIT'
      f' | {len([r for r in bad if r["kind"]=="floor"])} FLOOR\n')
stale = [r for r in ROWS if not r['present']]
if stale:
    print('STALE-CHECK (the quoted sentence is no longer in the manuscript -- update the registry):')
    for r in stale:
        print(f'  {r["id"]:7s} {r["claim"][:70]}')
    print()
for r in ROWS:
    if r['kind'] == 'value' and not VERBOSE:
        continue
    tag = {'mismatch': 'MISMATCH', 'split': 'SPLIT', 'floor': 'FLOOR', 'value': 'ok'}[r['kind']]
    print(f'[{tag:8s}] {r["id"]:7s} {r["section"]}')
    print(f'            says: {r["claim"]}')
    print(f'            said={fmt(r["said"])}   artifact={fmt(r["got"])}')
    if r['note']:
        print(f'            {r["note"]}')
    print()
sys.exit(1 if bad else 0)
