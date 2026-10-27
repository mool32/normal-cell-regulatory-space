#!/usr/bin/env python3
"""End-to-end regeneration: every figure the paper prints, then the PDF, the Word file and the checks.

Code availability promises "a single end-to-end script that regenerates the main and supplementary
figures from public inputs". This is it. Everything below runs from the committed derived tables in
`results/`; only ONE stage (v7_clouds) opens raw .h5ad, and it is marked so a reader without the
single-cell downloads can see exactly what they lose.

    python3 scripts/make_figures.py            run everything, stop at the first failure
    python3 scripts/make_figures.py --list     print the plan and the per-stage inputs, run nothing
    python3 scripts/make_figures.py --keep-going
    python3 scripts/make_figures.py --only E1  run only the stages whose script name contains "E1"

Order matters in three places and nowhere else: the panel builders must precede `assemble_figures`,
which composes six of the seven main figures (Fig. 6 is `build_figS29_invitro.py`'s own panel);
`build_SI_consolidate` must follow every panel it groups; and the
document builders must come last, since they read the PNGs. `si_figures_v12.py` is NOT a stage -- it
is the caption/numbering table that `build_v10_pdf.py` imports.
"""
import subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / 'scripts'

# The one stage that needs raw single-cell matrices rather than a committed derived table. They are
# public but large and are not redistributed with the code, so when they are absent the stage is
# SKIPPED and named -- a release that silently dropped a figure would be worse than one that says so.
RAW_INPUTS = {
    # the end-to-end re-derivation (Note S-W): the HPA source table, the 25 cohort matrices and
    # the Pelka epithelium cache. Its conclusions ship as results/E46_rederivation/*.
    'E46_rederive_reference.py': ['paper2/data/rna_single_cell_type.tsv'],
    'E46b_rescore_cohorts.py': ['results/cancer_v6_6_extras/Ledergor_MM_malignant.h5ad'],
    'E46d_rederive_colorectal.py': ['perceptome_validation/cache/pelka_epi.h5ad'],

    'build_v7_clouds.py': ['results/cancer_v6_4_extras/B-ALL_malignant.h5ad',
                           'results/cancer_v6_4_extras/Neuroblastoma_malignant.h5ad'],
    # These two read cohort matrices to establish provenance, which is exactly what cannot be
    # redistributed: in the release bundle they SKIP by name rather than fail, and the facts they
    # establish ship as results/census_cohorts.csv and results/cohort_provenance.csv.
    'verify/census_cohorts.py': ['results/cancer_scaled/cache/Prostate_malignant.h5ad'],
    'verify/paired_overlap.py': ['paper1/oscilatory/data/census_full/census_cancer_lung_tumor.h5ad',
                                 'results/cancer_phase2/data/pelka_colon/'
                                 'pelka2021_colon_all_cells.h5ad'],
    'verify/check_E4_lazy.py': ['results/cancer_v6_6_extras/vanGalen_AML_malignant.h5ad'],
}

# (stage group, script, what it writes / why it is here).  `raw` marks a stage needing .h5ad.
PLAN = [
    ('panels: Fig. 1',  'build_fig1A_workflow.py',   'Fig. 1A, the schematic'),
    ('panels: Fig. 1',  'build_fig1B_axes.py',       'Fig. 1B-C, the nine axes and the normal-cell map'),
    ('panels: Fig. 1',  'build_fig3F_1D.py',         'Fig. 1D and Fig. 3F (one builder, two figures)'),
    ('panels: Fig. 2',  'build_fig2_panels.py',      'Fig. 2C-D, contraction and its null'),
    ('panels: Fig. 2',  'build_fig2_maps.py',        'Fig. 2A-B, normal and malignant maps'),
    ('panels: Fig. 3',  'build_fig3A_map.py',        'Fig. 3A, where the cancers land'),
    ('panels: Fig. 3',  'build_fig3_panels.py',      'Fig. 3B, 3D'),
    ('panels: Fig. 3',  'build_fig3CE.py',           'Fig. 3C, 3E (+ the two SI ranking panels)'),
    ('panels: Fig. 4-6', 'build_figs456.py',         'Fig. 4, 5 and 6 panels'),
    ('assemble',        'assemble_figures.py',       'composes v11_fig1-5 and v11_fig6 (printed as Fig. 7) from the panels above'),

    ('SI panels',       'build_E0_bootstrap.py',     'Fig. S1-I axis stability'),
    ('SI panels',       'E7_eta2_figure.py',         'Fig. S1-II eta-squared'),
    ('SI panels',       'E2ab_modules_redundancy.py','Fig. S2 module redundancy'),
    ('SI panels',       'build_v7_axes.py',          'Fig. S3 axes against external signatures'),
    ('SI panels',       'build_v7_clouds.py',        'Fig. S4-I per-cancer clouds', 'raw'),
    ('SI panels',       'build_v7_map.py',           'Fig. S4-I normal-cell map'),
    ('SI panels',       'build_v7_SI.py',            'Fig. S5-I the complexity trap'),
    ('SI panels',       'build_E1_figure.py',        'Fig. S5-II the complexity control'),
    ('SI panels',       'build_v7_fig1.py',          'Fig. S6 the two ranking views'),
    ('SI panels',       'build_E2c_target.py',       'Fig. S7-I catalog dependence'),
    ('SI panels',       'build_E2d_robustness.py',   'Fig. S7-II decoy modules'),
    ('SI panels',       'build_E2e_robustness.py',   'Fig. S7-II leave-five modules'),
    ('SI panels',       'build_figS10_genelevel.py', 'Fig. S7-III the gene-level control'),
    ('SI panels',       'build_v7_duetprofile.py',   'Fig. S8 the target among 154 normal cells'),
    ('SI panels',       'build_E28_reprog.py',       'Fig. S9 reprogramming, no cancer in the data'),
    ('SI panels',       'build_E7_speccurve.py',     'Fig. S10 the 36-way specification curve'),
    ('SI panels',       'build_E10_gastric.py',      'Fig. S11-I the paired gastric cohort'),
    ('SI panels',       'build_E12_lung.py',         'Fig. S11-II corrected lung'),
    ('SI panels',       'build_E17_shapeonly.py',    'Fig. S11-III magnitude removed'),
    ('SI panels',       'build_E24_E25.py',          'Fig. S12 fibrosis + Fig. S13-III metastasis'),
    ('SI panels',       'build_E9_preinvasive.py',   'Fig. S13-I pre-invasive colon'),
    ('SI panels',       'build_E22_preinvasive2.py', 'Fig. S13-II Barrett / DCIS'),
    ('SI panels',       'build_E26_dynamics.py',     'Fig. S14 the coordinate under therapy'),
    ('SI panels',       'build_figS29_invitro.py',   'Fig. 6 isogenic transformation in vitro'),
    ('SI panels',       'build_v7_detuning.py',      'Fig. S15 what malignancy adds'),
    ('assemble',        'build_SI_consolidate.py',   'groups 26 panels into the 15 numbered SI figures'),

    ('document',        'build_v10_pdf.py',          'main + SI sources -> paper3_v12_MAIN.pdf and paper3_v12_SI.pdf'),
    ('document',        'build_v10_docx.py',         'the same two sources -> paper3_v12_MAIN.docx and paper3_v12_SI.docx'),
    ('re-derivation',   'E46_rederive_reference.py',
     'Note S-W stage 1-2: the corrected HPA reference and eigenspace'),
    ('re-derivation',   'E46b_rescore_cohorts.py',
     'Note S-W stage 3: the 25 cohorts re-scored on corrected gene names'),
    ('re-derivation',   'E46d_rederive_colorectal.py',
     'Note S-W stage 4: the colorectal within-patient arm re-derived'),
    ('checks',          'E45_genuine_label_subset.py',
     'Note S-V: the naming recomputed without the heuristically-labelled cohorts'),
    ('checks',          'E46c_rederived_results.py',
     'Note S-W: the re-derived catalog, occupancy and recurrence'),
    ('checks',          'verify/check_numbers.py',   'every registered claim recomputed from its artifact'),
    ('checks',          'verify/cohort_provenance.py','every cohort resolved to its accession'),
    ('checks',          'verify/check_citations.py', 'every citation, cross-reference and figure part resolved'),
    ('checks',          'verify/census_cohorts.py',  'the Census cohorts resolved: selection rule and double counting'),
    ('checks',          'verify/check_E4_lazy.py',   "E4's lazy cohort loader equals its eager one"),
    ('checks',          'verify/paired_overlap.py',  'how much of each paired cohort is already in the panel'),
]


def main():
    args = sys.argv[1:]
    only = args[args.index('--only') + 1] if '--only' in args else None
    plan = [p for p in PLAN if not only or only in p[1]]
    # The 're-derivation' group re-scores all 25 cohort matrices and the Pelka epithelium from raw
    # .h5ad: about half an hour against three minutes for everything else, and its outputs are
    # committed (results/E46_rederivation/). It is therefore opt-in -- `--rederive` -- while
    # E46c_rederived_results.py re-checks the published Note S-W numbers from those outputs on
    # every run. Without this the one-command build stops being one command anyone will run.
    if '--rederive' not in args:
        heavy = [s for g, s, *_ in plan if g == 're-derivation']
        plan = [p for p in plan if p[0] != 're-derivation']
        if heavy and '--list' not in args:
            print(f"  (skipping {len(heavy)} re-derivation stage(s); pass --rederive to run them: "
                  f"{', '.join(heavy)})", flush=True)
    if '--list' in args:
        for g, s, why, *flag in plan:
            print(f"  {g:17s} {s:30s} {'[needs .h5ad] ' if flag else ''}{why}")
        print(f"\n{len(plan)} stages. Missing scripts: "
              f"{[s for _, s, *_ in plan if not (S / s).exists()] or 'none'}")
        return 0
    rows, failed, gates, skipped = [], [], [], []
    for g, s, why, *flag in plan:
        p = S / s
        if not p.exists():
            rows.append((g, s, 'MISSING', 0.0)); failed.append(s); continue
        _miss = [f for f in RAW_INPUTS.get(s, []) if not (ROOT / f).exists()]
        if _miss:
            rows.append((g, s, 'skipped', 0.0)); skipped.append(s)
            print(f"  SKIP {s:30s}    -    {why}\n       needs {', '.join(_miss)}", flush=True)
            continue
        t0 = time.time()
        r = subprocess.run([sys.executable, str(p)], cwd=ROOT, capture_output=True, text=True)
        dt = time.time() - t0
        ok = r.returncode == 0
        # A verifier exits non-zero when claims are unresolved. That is a finding about the paper,
        # not a broken build, so it is labelled GATE and never stops the run.
        gate = (not ok) and g == 'checks'
        rows.append((g, s, 'ok' if ok else f'rc={r.returncode}', dt))
        print(f"  {'ok  ' if ok else ('GATE' if gate else 'FAIL')} {s:30s} {dt:6.1f}s  {why}", flush=True)
        if not ok:
            (gates if gate else failed).append(s)
            print('\n'.join('      | ' + l for l in (r.stderr or r.stdout).strip().splitlines()[-12:]))
            if not gate and '--keep-going' not in args:
                break
    total = sum(d for *_, d in rows)
    print(f"\n{len(rows)} stages, {len(failed)} failed, {len(skipped)} skipped for missing raw "
          f"inputs, {len(gates)} gate(s) reporting, {total:.0f}s total")
    if failed:
        print('failed: ' + ', '.join(failed))
    if skipped:
        print('skipped (raw inputs absent): ' + ', '.join(skipped))
    if gates:
        print('open findings reported by: ' + ', '.join(gates))
    return 1 if (failed or gates) else 0


if __name__ == '__main__':
    sys.exit(main())
