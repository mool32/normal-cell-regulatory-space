#!/usr/bin/env python3
"""Assemble the release bundle: exactly the files the paper's own pipeline reads, and nothing else.

    python3 scripts/make_release.py [--out DIR] [--no-verify]

How the file list is decided. Not by hand and not by a checked-in list that goes stale: this script
RUNS `scripts/make_figures.py` under a Python audit hook that records every file opened inside the
repository, and copies that set. Whatever the pipeline actually touched is in the bundle; whatever it
did not is out. It then re-runs the whole pipeline INSIDE the copy, so "regenerates from public
inputs" is a statement that was tested rather than hoped for.

Two things the trace cannot see and the bundle therefore names explicitly:
  * the two single-cell matrices Fig. S4-I opens through h5py's C layer (RAW_EXTRAS below). They are
    public, large, and listed by accession in the README rather than copied.
  * the scripts cited in captions that no figure needs (the analysis that produced a committed table).
    Those are collected from the manuscript source so the citation and the bundle cannot disagree.

NOTE on `cancer_paper_release/`: that directory is the PREVIOUS paper's release ("Convergent
perceptual identity in cancers across organ systems", 11 cancers, TME/CAF analysis). None of its
sixteen step*.py scripts is cited by this manuscript. It is not what this script builds and it must
not be the Zenodo deposit for this paper.
"""
import hashlib, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_EXTRAS = {  # opened below the Python layer, so the trace never sees them
    'results/cancer_v6_4_extras/B-ALL_malignant.h5ad': 'GSE132509 (B-ALL), Fig. S4-I only',
    'results/cancer_v6_4_extras/Neuroblastoma_malignant.h5ad': 'GSE137804 (neuroblastoma), Fig. S4-I only',
}
ALWAYS = ['manuscript/paper3_v10_longform.md', 'manuscript/paper3_v10_SI.md',
          'manuscript/note_ST_activity.md',
          'manuscript/paper3_v12_MAIN.pdf', 'manuscript/paper3_v12_SI.pdf',
          'manuscript/paper3_v12_MAIN.docx', 'manuscript/paper3_v12_SI.docx',
          'scripts/make_figures.py', 'scripts/make_release.py', 'scripts/figstyle.py',
          'scripts/si_figures_v12.py', 'LICENSE', 'requirements.txt',
          'CITATION.cff', '.gitignore']

# The audit hook sees data files (io.open). It does NOT see imports -- the import machinery uses
# io.open_code, which raises no `open` event -- so a first version of this bundle shipped every CSV
# and no `prereg_lib.py`, and twelve stages died on ModuleNotFoundError inside the copy. The atexit
# sweep over sys.modules adds the imported sources.
SHIM = """import sys, os, atexit
_log = os.environ.get('AUDIT_LOG')
if _log:
    _f = open(_log, 'a', buffering=1)
    _root = os.environ.get('AUDIT_ROOT', '')
    def _hook(event, args):
        if event == 'open' and args and isinstance(args[0], str):
            try: p = os.path.abspath(args[0])
            except Exception: return
            if p.startswith(_root): _f.write(p + '\\n')
    sys.addaudithook(_hook)
    @atexit.register
    def _mods():
        for m in list(sys.modules.values()):
            f = getattr(m, '__file__', None)
            if not f: continue
            try: f = os.path.abspath(f)
            except Exception: continue
            if f.startswith(_root): _f.write(f + '\\n')
"""


def traced_inputs():
    """Every file inside the repo that a full pipeline run opens."""
    with tempfile.TemporaryDirectory() as td:
        (Path(td) / 'usercustomize.py').write_text(SHIM)
        log = Path(td) / 'opened.txt'; log.touch()
        env = dict(os.environ, AUDIT_LOG=str(log), AUDIT_ROOT=str(ROOT) + os.sep,
                   PYTHONPATH=td + os.pathsep + os.environ.get('PYTHONPATH', ''))
        print('tracing a full pipeline run ...', flush=True)
        subprocess.run([sys.executable, str(ROOT / 'scripts/make_figures.py'), '--keep-going'],
                       cwd=ROOT, env=env, capture_output=True, text=True)
        seen = {l.strip() for l in log.read_text().splitlines() if l.strip()}
    rel = set()
    for p in seen:
        f = Path(p)
        if f.is_file():
            rel.add(str(f.relative_to(ROOT)))
    return rel


def cited_scripts():
    """Script filenames the manuscript cites, resolved to their real location in the tree."""
    src = (ROOT / 'manuscript/paper3_v10_full_src.md').read_text()
    names = set(re.findall(r'[A-Za-z0-9_]+\.py', src))
    out, missing = set(), []
    for n in sorted(names):
        hits = [str(Path(dp).relative_to(ROOT) / n) for dp, _, fs in os.walk(ROOT)
                if n in fs and '/.git' not in dp]
        if hits: out.add(sorted(hits, key=len)[0])
        else: missing.append(n)
    if missing:
        print('  cited but NOT FOUND: ' + ', '.join(missing))
    return out



def _raw_block():
    """The raw inputs each stage needs, read from make_figures.RAW_INPUTS so the two lists cannot
    drift apart; RAW_EXTRAS supplies the accession note where there is one."""
    import importlib.util as ilu
    spec = ilu.spec_from_file_location('_mf', ROOT / 'scripts/make_figures.py')
    mf = ilu.module_from_spec(spec); spec.loader.exec_module(mf)
    why = {s: w for g, s, w, *_ in mf.PLAN}
    out = []
    for stage, needs in mf.RAW_INPUTS.items():
        out.append(f'  - `{stage}` \u2014 {why.get(stage, "")}')
        for f in needs:
            note = RAW_EXTRAS.get(f)
            out.append(f'      - `{f}`' + (f' \u2014 {note}' if note else ''))
    return '\n'.join(out)



def write_manifest(out, files):
    """sha256 of everything that actually shipped, README included, MANIFEST itself excluded."""
    man = ['# MANIFEST  (sha256  size  path)']
    for f in sorted(set(files) | {'README.md'}):
        q = out / f
        if not q.is_file():
            continue
        man.append(f'{hashlib.sha256(q.read_bytes()).hexdigest()}  {q.stat().st_size:>10d}  {f}')
    (out / 'MANIFEST.txt').write_text('\n'.join(man) + '\n')


def main():
    out = Path(sys.argv[sys.argv.index('--out') + 1]) if '--out' in sys.argv else ROOT / 'release'
    files = traced_inputs() | cited_scripts() | {f for f in ALWAYS if (ROOT / f).exists()}
    files |= {str(p.relative_to(ROOT)) for p in (ROOT / 'manuscript/figures').glob('*.png')}
    # The pre-registration record. No script OPENS these, so the trace never sees them -- and the
    # bundle shipped the code that runs T1-T9 without the frozen prompts, the amendments log or the
    # verdicts that make those tests pre-registered. The paper cites them by path.
    files |= {str(q.relative_to(ROOT)) for pat in ('perceptome_validation/registry/**/*',
                                                   'perceptome_validation/*.md')
              for q in ROOT.glob(pat) if q.is_file()}
    # bytecode caches get traced like any other read, but shipping them makes MANIFEST.txt
    # disagree with a fresh clone (.gitignore drops them) and they regenerate anyway
    files = {f for f in sorted(files)
             if (ROOT / f).is_file() and '__pycache__' not in f and not f.endswith('.pyc')}
    if out.exists():
        shutil.rmtree(out)
    nbytes = 0
    for f in sorted(files):
        dst = out / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / f, dst)
        nbytes += dst.stat().st_size
    (out / 'README.md').write_text(README.format(
        n=len(files), mb=nbytes / 1e6,
        raw=_raw_block()))
    print(f'\n{len(files)} files, {nbytes/1e6:.1f} MB -> {out}')

    if '--no-verify' in sys.argv:
        write_manifest(out, files)
        return 0
    print('\nre-running the pipeline inside the copy (this is the actual test) ...', flush=True)
    r = subprocess.run([sys.executable, 'scripts/make_figures.py', '--keep-going'],
                       cwd=out, capture_output=True, text=True)
    tail = (r.stdout or '').strip().splitlines()
    print('\n'.join('  ' + l for l in tail[-6:]))
    bad = [l for l in tail if l.strip().startswith('FAIL')]
    # AFTER the re-run, not before: the re-run regenerates the documents, and a PDF carries a
    # creation timestamp, so a manifest written at copy time disagreed with the seven shipped
    # documents in every bundle this script has ever produced. A checksum list that does not
    # describe what shipped is worse than none.
    write_manifest(out, files)
    print('\nRELEASE VERIFIED: the bundle rebuilds every figure and both documents on its own'
          if not bad else f'\n{len(bad)} stage(s) failed inside the bundle:\n' + '\n'.join(bad))
    return 1 if bad else 0


README = """# Human malignancies show increased occupancy of a shared region of normal-cell regulatory space

Code and derived data for the paper. {n} files, {mb:.0f} MB.

## Regenerate everything

```
python3 scripts/make_figures.py
```

About two minutes. It rebuilds every panel, composes the seven main figures and the fifteen
supplementary ones, writes the main text and the Supporting Information as separate documents
(`manuscript/paper3_v12_MAIN.pdf`/`.docx` and `paper3_v12_SI.pdf`/`.docx`), and finishes by running
the verification harnesses. Run it from this directory; the scripts resolve their inputs relative to
the repository root. Two consecutive runs produce byte-identical figures.

## What is here, and what is not

Everything under `results/`, `tool/` and `perceptome_validation/` is a *derived* table that the
figure code reads — the trace that built this bundle recorded each one being opened. The raw
single-cell matrices are not redistributed; they are public, and the accessions are listed in Data
availability. Four stages open one locally and name it before skipping:

{raw}

Everything else builds without them. The three checks among those four establish facts that ship
here as finished tables (`results/census_cohorts.csv`, `results/cohort_provenance.csv`); a check
that cannot run says what it needs rather than passing silently.

## Verification

- `scripts/verify/check_numbers.py` recomputes every load-bearing number in the text from one named
  artifact each, and fails loudly when the artifact and the sentence disagree.
- `scripts/verify/cohort_provenance.py` resolves each of the 25 cohorts to its accession.
- `MANIFEST.txt` gives the sha256 of every file in this bundle.

A non-zero exit from `check_numbers.py` is a report about unresolved claims, not a broken build; the
runner labels those stages GATE.
"""

if __name__ == '__main__':
    sys.exit(main())
