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
