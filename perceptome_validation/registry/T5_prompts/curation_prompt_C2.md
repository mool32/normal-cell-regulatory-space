# T5 curation prompt (frozen 2026-09-19; identical for C2, C3, C4 except the output path and the seed word)

You are an independent curator building a catalog of cellular "regulatory modules" for a pre-registered study. Work ONLY from your
own knowledge of cell signalling and from public databases (KEGG, Reactome, primary literature — you may use web search/fetch).
Do NOT read, list or search any file on this machine (no local files at all), do NOT look for any existing catalog, and do NOT try to
infer what result the study is testing. Your catalog will be evaluated blind.

## Definition (verbatim)
A regulatory module is a self-contained unit that reads one signal and produces a regulated output. It is specified by FOUR
operational criteria, all of which must hold:
  1. a **sensor input** — a receptor or sensor protein that detects the signal (ligand, metabolite, physical stress, nucleic acid…);
  2. a **transcription-factor output with documented targets** — the pathway terminates in a TF (or TF complex) with a documented
     target-gene set;
  3. a **negative-feedback element enabling cyclical operation** — a documented transcriptional or post-translational negative-
     feedback component that returns the module to baseline (an inhibitor induced by the output, a degrader of the sensor, etc.);
  4. **receptor-field specificity** — the module responds to a distinct class of signal (its "receptor field") not already covered by
     another retained module.

## Task
Screen candidate signalling AND metabolic/homeostatic axes of human cells (aim to consider at least ~70 candidates spanning growth
factor/cytokine/hormone/nutrient/stress/nucleic-acid/mechanical/circadian/xenobiotic sensing) against the four criteria. RETAIN every
axis that passes all four; DROP any that fails a criterion or is REDUNDANT with an axis already retained (same receptor field or
same TF output). Target catalog size: **40–50 modules, 4–6 genes each**. For each retained module list its **component genes**
(HGNC symbols, current nomenclature): the sensor, key cascade element(s), the TF, and the negative-feedback element — the genes whose
expression measures how well equipped a cell is to run the module (NOT the TF's target genes). Use only protein-coding genes.

## Output
Write a single JSON file to: /private/tmp/claude-501/-Users-teo-Desktop-Desktop---teo-s-MacBook-Air-research-perceptual-modules/108e5a57-d677-4b97-b4ad-48765c5cebd6/scratchpad/T5/catalog_C2.json
Format:
{
  "curator": "C2",
  "screened_candidates": ["<name>", ...],              // every axis you considered (>= 70)
  "dropped": {"<name>": "<which criterion failed / redundant with which module>", ...},
  "modules": {
     "<Module name>": {"genes": ["SYM1","SYM2","SYM3","SYM4","SYM5"],
                       "sensor": "SYM", "tf": "SYM", "feedback": "SYM", "receptor_field": "<one phrase>"},
     ...
  }
}
Then reply with ONLY the text: DONE <path> <number of modules>.
Do not include any prose in the JSON; do not write any other file.
