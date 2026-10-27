"""Build an editable Word version of the manuscript from the same generated source as the PDF.
Run scripts/build_v10_pdf.py first (it writes the two sources), then this. Since the co-author
review of 2026-10-27 the main text and the SI are separate documents with no table of contents, so
this writes paper3_v12_MAIN.docx and paper3_v12_SI.docx.

The source carries three wide tables as raw `{=latex}` longtables (Table S2, S7 and S9b): pandoc's
pipe tables gave uncontrollable column widths in LaTeX, so the PDF builder writes them by hand. A
raw LaTeX block is silently DROPPED when the target is docx, which is how the Word file shipped for
weeks with those three tables simply absent -- including Table S9b, which exists nowhere else. They
are therefore parsed back out of the LaTeX here and re-emitted as pipe tables.
"""
import re, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; MAN=ROOT/'manuscript'

PAGEBREAK='```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```'
# the escapes scripts/build_v10_pdf.py::_tex applies, undone in reverse (backslash last)
_UNESC=[('\\textasciicircum{}','^'),('\\textasciitilde{}','~'),('\\}','}'),('\\{','{'),
        ('\\_','_'),('\\#','#'),('\\$','$'),('\\%','%'),('\\&','&'),('\\textbackslash{}','\\')]
def _untex(s):
    for a,b in _UNESC: s=s.replace(a,b)
    return s.strip()

def _rows(tex):
    """Every `a & b & c \\\\` line of a longtable body, in order, as lists of cells."""
    out=[]
    for ln in tex.splitlines():
        ln=ln.strip()
        if not ln.endswith('\\\\') or ln.startswith('\\'): continue
        out.append([_untex(c) for c in ln[:-2].split(' & ')])
    return out

def latex_longtable_to_md(block):
    head=re.search(r'\\midrule\s*\n(.*?)\\\\', block, re.S)       # the \textbf{..} header line
    hdr=[_untex(c.replace('\\textbf{','').rstrip('}'))
         for c in re.search(r'\\toprule\s*\n(.*?)\s*\\\\', block, re.S).group(1).split(' & ')]
    body=_rows(block.split('\\endhead',1)[-1])
    body=[r for r in body if len(r)==len(hdr)]
    if not body: return ''
    md=['| '+' | '.join(hdr)+' |', '|'+'---|'*len(hdr)]
    md+= ['| '+' | '.join(c.replace('|','\\|') for c in r)+' |' for r in body]
    return '\n'.join(md)

def convert_raw(m):
    block=m.group(1)
    if '\\begin{longtable}' in block: return '\n\n'+latex_longtable_to_md(block)+'\n\n'
    if '\\clearpage' in block or '\\newpage' in block: return '\n\n'+PAGEBREAK+'\n\n'
    return ''                                                     # anything else is PDF-only styling

def to_docx(src_path, out_path):
    src = src_path.read_text()
    n_raw = len(re.findall(r'```\{=latex\}\n(.*?)\n```', src, re.S))
    src = re.sub(r'```\{=latex\}\n(.*?)\n```', convert_raw, src, flags=re.S)
    src = src.replace('\\newpage', PAGEBREAK)
    assert '\\begin{longtable}' not in src, 'a longtable survived the conversion'
    print(f'{src_path.name}: {n_raw} raw LaTeX blocks converted; '
          f'{src.count("|---|")} tables now in pipe form')
    tmp = MAN / ('_docx_src_' + src_path.stem + '.md'); tmp.write_text(src)
    cmd = ['pandoc', str(tmp), '-f', 'markdown+pipe_tables+raw_attribute', '-o', str(out_path),
           '--resource-path', str(MAN)]
    ref = MAN / '_docx_reference.docx'
    if ref.exists(): cmd += ['--reference-doc', str(ref)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout or '', r.stderr or '')
    tmp.unlink(missing_ok=True)
    print(f'wrote {out_path} ({out_path.stat().st_size/1e6:.2f} MB)' if out_path.exists() else 'FAILED')
    return out_path.exists()

ok = [to_docx(MAN / 'paper3_v12_MAIN.md', MAN / 'paper3_v12_MAIN.docx'),
      to_docx(MAN / 'paper3_v12_SI.md', MAN / 'paper3_v12_SI.docx')]
if not all(ok): raise SystemExit('a .docx failed to build')
