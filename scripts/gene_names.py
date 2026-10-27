"""Gene-name resolution for the module catalog, and the two defects it fixes.

Found by the round-4 gene-recovery audit (scripts/E43_region_complexity_control.py):

1. `NKX3-1` was lost in **all 17** scored cohorts. Not a synonym problem: the pipeline's own
   de-duplication stripped any trailing `-<digits>` from every gene name, to undo the suffixes
   `var_names_make_unique()` adds. `NKX3-1` is a real HGNC symbol, so it became `NKX3` and matched
   nothing. `collapse_unique_suffix` below strips the suffix only when the bare name is actually
   present elsewhere in the matrix, which is the only case the suffix can have come from make_unique.

2. Deprecated symbols in the catalog itself: `MB21D1` (now CGAS) and `TMEM173` (now STING1) in
   cGAS-STING, and `ARNTL` (now BMAL1) in Circadian. Datasets annotated against a recent GENCODE
   carry the current symbols, older ones the deprecated ones, and the two colon cohorts carry a mix —
   so the module's effective gene set varied by dataset. `resolve_to_catalog` maps whichever form the
   data uses onto the symbol the catalog scores.
"""
import re

# deprecated/alternate -> the symbol the catalog uses
SYNONYMS = {
    'CGAS': 'MB21D1', 'C6ORF150': 'MB21D1',
    'STING1': 'TMEM173', 'ERIS': 'TMEM173', 'MITA': 'TMEM173',
    'BMAL1': 'ARNTL', 'MOP3': 'ARNTL',
    'NKX3A': 'NKX3-1',
}

def collapse_unique_suffix(names):
    """Undo `var_names_make_unique()` suffixes WITHOUT damaging real symbols such as NKX3-1."""
    present = set(names)
    out = []
    for v in names:
        base = re.sub(r'-\d+$', '', v)
        out.append(base if (base != v and base in present) else v)
    return out

def resolve_to_catalog(names):
    """Rewrite a dataset's gene names into the catalog's symbol space."""
    return [SYNONYMS.get(v, v) for v in collapse_unique_suffix(list(names))]

def resolve_matrix(adata):
    """Rewrite an AnnData's gene names into catalog space, keeping the column that carries the data.

    The third defect, found 2026-10-27 while re-deriving the catalog (PNAS v13 report 3.1): several
    cohorts carry BOTH the current and the deprecated symbol as separate features, and the
    deprecated one is identically zero. In the Pelka colon object CGAS is nonzero in 12% of cells
    while MB21D1 is empty, STING1 in 18% while TMEM173 is empty, BMAL1 in 5% while ARNTL is empty --
    and ARNTL sits at var index 1077 against BMAL1 at 1597. Renaming and then dropping duplicates by
    first occurrence therefore KEEPS THE EMPTY COLUMN for Circadian, so the rename alone would have
    scored the module on a zero and reported it as repaired.

    When a rename collides, keep whichever column has more nonzero entries. These are two
    annotations of one gene, never two genes, so nothing is lost; summing would be defensible too
    but risks double counting if a dataset ever carries real counts in both.
    """
    import numpy as np, pandas as pd
    names = resolve_to_catalog(list(adata.var_names))
    s = pd.Series(names)
    dup = set(s[s.duplicated(keep=False)].tolist())
    if not dup:
        keep = np.ones(len(s), bool)
    else:
        X = adata.X
        nz = np.asarray((X != 0).sum(axis=0)).ravel()
        keep = np.zeros(len(s), bool)
        best = {}
        for i, nm in enumerate(names):
            if nm not in dup:
                keep[i] = True
            elif nm not in best or nz[i] > nz[best[nm]]:
                best[nm] = i
        for i in best.values():
            keep[i] = True
    out = adata[:, keep].copy()
    out.var_names = s[keep].values
    return out
