#!/usr/bin/env python3
"""Check every cross-reference in the built manuscript resolves, and every reference is used.

Four classes, each a different way a reference can be wrong, and none of them visible to a
spell-checker or to `check_numbers.py` (which recomputes values, not pointers):

  CITE   an author-year citation in the text with no entry in the reference list
  ORPHAN an entry in the reference list that nothing cites
  XREF   a "Fig. 3B" / "Fig. S7-II" / "Note S-K" / "Table S9b" pointer with no such caption
  DUP    two captions claiming the same number

Author matching is on surname + year, because the in-text form ("Vento-Tormo et al., 2018") and
the list form ("Vento-Tormo R et al. (2018)") agree on nothing else. Collective in-text citations
("Hanahan & Weinberg, 2000, 2011; Hanahan, 2022") are split on ';' and then on ',' so that a bare
trailing year attaches to the surname it follows.

Run against manuscript/paper3_v10_full_src.md, the GENERATED file the PDF is built from, so that
SI notes injected by the builder are checked too.
"""
import os, re, sys, collections

ROMAN = ['I', 'II', 'III', 'IV']

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
SRC = os.path.join(ROOT, 'manuscript', 'paper3_v10_full_src.md')

# A surname may carry spaces, hyphens or an accent; it never carries a digit.
# A surname may lead with a lowercase particle (van Galen) and may be two words.
NAME = r"(?:(?:van|von|de|del|della|di|da|le|la) )?[A-ZÅÄÖÜÉ][\w’.-]*(?:[ -][A-ZÅÄÖÜÉ][\w’.-]*)?"
YEAR = r"(?:19|20)\d\d[a-z]?"


def _groups():
    sys.path.insert(0, os.path.join(ROOT, 'scripts'))
    import build_SI_consolidate as b
    return b.GROUPS


def load():
    src = open(SRC).read()
    # The reference list itself must not be scanned for citations: every entry looks like one.
    m = re.search(r'\n## References\n(.*?)\n## ', src, re.S)
    if not m:
        sys.exit('no reference list found')
    reflist, body = m.group(1), src[:m.start()] + src[m.end():]
    return src, body, reflist


def parse_reflist(reflist):
    """-> {(surname, year): raw entry}, plus the raw list in order."""
    out, order = {}, []
    for e in (x.strip() for x in reflist.replace('\n', ' ').split('·')):
        if not e:
            continue
        m = re.match(rf'({NAME})[^(]*\(({YEAR})\)', e)
        if not m:
            print(f'  [UNPARSED] {e[:70]}')
            continue
        key = (m.group(1), m.group(2))
        if key in out:
            print(f'  [DUP REF ] {key} twice in the list')
        out[key] = e
        order.append(key)
    return out, order


# One in-text citation: a surname (optionally "et al." or "& Other"), then one or more years.
# The years group is greedy over ", 2011" repeats so that "Hanahan & Weinberg, 2000, 2011" yields
# both. A year must not sit inside a longer token, or GSE131907 reads as "1907".
CITE_RE = re.compile(
    rf'({NAME})'
    rf'(?: et al\.| &amp; {NAME}| & {NAME}|,? and {NAME})?'
    rf'[,.]? \(?((?:19|20)\d\d[a-z]?(?:, ?(?:19|20)\d\d[a-z]?)*)\)?(?![\d\w])')
# "GSE131907" and "bioRxiv 2023.01.26" are not citations; neither is a figure/table/note label.
NOT_A_CITE = re.compile(r'(?:GSE|GSM|SRP|PRJ|CRA|E-MTAB-|bioRxiv |version |doi|zenodo)\s*\d*$', re.I)
# This paper's own Zenodo deposit, cited on the title block; it is not a reference-list entry.
SELF = {('Spiro', '2026'), ('T.', '2026')}   # "Spiro, T. (2026)", this paper's own Zenodo deposit


def intext(body):
    """Every (surname, year) the text cites, with the fragment it sits in."""
    found = []
    for m in CITE_RE.finditer(body):
        before = body[max(0, m.start() - 12):m.start()]
        if NOT_A_CITE.search(before) or re.search(r'[\w.]$', before):
            continue
        for y in re.findall(YEAR, m.group(2)):
            found.append((m.group(1), y, m.group(0)[:70]))
    return found


# A caption carries the number followed by a period or an em-dash title; the SI head lists the
# same notes as bare **Note S-N** cross-references, and the SI table index as **Table S3**.
# `\textbf{...}` because Table 1's caption is a row of its own longtable head (scripts/build_v10_pdf.py),
# which is where a caption has to live to be unable to drift onto another page from its table.
CAPTION = re.compile(r'(?:\*\*|^#+ |\\textbf\{)(?:(Fig\.|Figure) (S?\d+[a-z]?)|(Table) (S?\d+[a-z]?|S-[A-Z])'
                     r'|(Note) (S\d+|S-[A-Z]))(?=\.| —)', re.M)
# Letters deliberately retired: the document says so where the series is listed.
RETIRED = {('Note', 'S-E')}
POINTER = re.compile(r'\b(?:(Figs?\.|Figures?) (S?\d+[a-z]?)(?:-([IVX]+))?(?:([A-Z])\b)?'
                     r'|(Tables?) (S?\d+[a-z]?|S-[A-Z])'
                     r'|(Notes?) (S\d+|S-[A-Z]))')


def main():
    src, body, reflist = load()
    print(f'{os.path.basename(SRC)}  {len(src):,} chars\n')
    refs, order = parse_reflist(reflist)
    cites = intext(body)

    cited = collections.Counter((s, y) for s, y, _ in cites)
    bad = [(s, y, ctx) for s, y, ctx in cites if (s, y) not in refs]
    # Tolerate a surname the list spells with more parts (e.g. "Hanahan & Weinberg" -> "Hanahan D")
    bad = [b for b in bad if not any(k[0].split()[0] == b[0].split()[0] and k[1] == b[1] for k in refs)
           and (b[0], b[1]) not in SELF]
    orphan = [k for k in order if cited[k] == 0
              and not any(c[0].split()[0] == k[0].split()[0] and c[1] == k[1] for c in cited)]

    # ---- captions present, and pointers into them
    have, dup = {}, []
    for m in CAPTION.finditer(src):
        kind = (m.group(1) and 'Fig') or (m.group(3) and 'Table') or 'Note'
        num = m.group(2) or m.group(4) or m.group(6)
        if (kind, num) in have:
            dup.append((kind, num))
        have[(kind, num)] = m.start()
    # A composite SI figure's parts are drawn by scripts/build_SI_consolidate.py, so how many a
    # figure has is knowable; a main figure's panel letters are listed in its own caption.
    parts = {f'S{n[1:]}' if n.startswith('S') else n: sum(len(r) for r in rows)
             for n, (_, rows, _) in _groups().items()}
    panels = {}
    for (kind, num), pos in have.items():
        if kind == 'Fig':
            cap = src[pos:pos + 4000].split('\n')[0]
            panels[num] = set(re.findall(r'\*\*\(([A-Z])\)\*\*|(?<![\w(])\(([A-Z])\)(?= )', cap))
            panels[num] = {a or b for a, b in panels[num]}

    # "Fig. 3B, 3C, 3E and 3F": every letter after the first carries no "Fig.", so POINTER sees
    # only 3B. Expand the list into explicit pointers before scanning.
    body = re.sub(r'(Figs?\. )(S?\d+)([A-Z])((?:,? (?:and )?\d*[A-Z]\b)+)',
                  lambda m: m.group(1) + m.group(2) + m.group(3) + ''.join(
                      f', Fig. {m.group(2)}{x[-1]}' for x in re.findall(r'\d*[A-Z]\b', m.group(4))),
                  body)

    missing, badpart = collections.Counter(), collections.Counter()
    for m in POINTER.finditer(body):
        if m.group(1):
            kind, num, part, letter = 'Fig', m.group(2), m.group(3), m.group(4)
        elif m.group(5):
            kind, num, part, letter = 'Table', m.group(6), None, None
        else:
            kind, num, part, letter = 'Note', m.group(8), None, None
        if (kind, num) not in have and (kind, num) not in RETIRED:
            missing[(kind, num)] += 1
            continue
        if part is not None:                       # "Fig. S13-II"
            n = parts.get(num)
            k = ROMAN.index(part) + 1 if part in ROMAN else 99
            if n is None or k > n:
                badpart[(f'Fig. {num}-{part}', f'{n} part(s)' if n else 'not a composite')] += 1
        elif letter and num.startswith('S') is False and panels.get(num) and letter not in panels[num]:
            badpart[(f'Fig. {num}{letter}', 'caption lists ' + ''.join(sorted(panels[num])))] += 1

    print(f'{len(order)} reference entries | {len(set(cited))} distinct citations | '
          f'{len(have)} captions ({sum(1 for k in have if k[0]=="Fig")} figures, '
          f'{sum(1 for k in have if k[0]=="Table")} tables, {sum(1 for k in have if k[0]=="Note")} notes)')
    print(f'  {len(bad)} CITE | {len(orphan)} ORPHAN | {len(missing)} XREF | '
          f'{len(badpart)} PART | {len(dup)} DUP\n')

    for s, y, ctx in bad:
        print(f'[CITE  ] {s} {y}   cited in: {ctx}')
    for k in orphan:
        kind = 'ORPHAN' if k[0].split()[0] not in body else 'NOYEAR'
        print(f'[{kind}] {k[0]} {k[1]}   {"never cited" if kind == "ORPHAN" else "named without a year"}'
              f': {refs[k][:62]}')
    for (kind, num), n in sorted(missing.items()):
        print(f'[XREF  ] {kind} {num}   pointed at {n}x, no such caption')
    for (ref, why), n in sorted(badpart.items()):
        print(f'[PART  ] {ref:14s} pointed at {n}x, but {why}')
    for kind, num in dup:
        print(f'[DUP   ] {kind} {num}   two captions')
    hard = [k for k in orphan if k[0].split()[0] not in body]
    missing.update(badpart)
    return 1 if (bad or hard or missing or dup) else 0


if __name__ == '__main__':
    sys.exit(main())
