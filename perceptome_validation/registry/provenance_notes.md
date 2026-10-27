# Provenance notes — things found while reproducing the paper-3 pipeline (not amendments; recorded for the record)

1. **Bundled pseudobulk reference vs current HPA TSV.** Re-deriving the 154×44 pseudobulk R from
   `results/hpa_download/rna_single_cell_type.tsv` reproduces the bundled `hpa_perceptivity_v03.npz` in 42/44 modules exactly
   and differs in 18 of 6,776 entries: cGAS-STING for 8 kidney cell types (bundled ≈ 1.5× re-derived) and PXR/CAR for 10 brain
   cell types (bundled ≈ 1.25× re-derived). The bundled single_cell_scaled reference reproduces exactly (max |Δ| 1e-4). The bundled
   pseudobulk file therefore predates a small HPA/catalog change. Consequence: the REBUILT incumbent eigenspace (|cos| to the bundled
   loadings 0.98–1.00 per PC) gives E1 three-cell 24.7 % / two-cell 21.9 % and cytotrophoblast rank 2 (18/25), versus 27.1 % / 23.9 %
   and rank 1 (18/25) under the bundled loadings. Every T1 panel — the incumbent included — goes through the rebuilt pipeline, so
   the comparison is fair; the bundled numbers are reported as the paper reference (T0 gate G3).
2. **NKX3-1.** The paper's cancer loader `nv()` strips a trailing "-<digits>" from gene symbols (to drop make_unique suffixes), which
   turns NKX3-1 into NKX3; the AR module is therefore scored on 4 genes in every cancer cohort but on 5 in the HPA reference.
   `scripts/E17` (the E3 endpoint) does not apply this stripping. The harness applies the paper's main pipeline (nv()) uniformly to
   E1/E2/E3, which makes the Pelka baseline +15.0 raw (43/46) / +9.9 shape-only (43/46) instead of the paper's +15.5 / +10.5
   (43–44/46). Internal comparisons (panel vs incumbent) are unaffected.
3. **Gene symbols absent from HPA.** MB21D1, TMEM173 and ARNTL are not in the HPA TSV (HPA uses CGAS, STING1, BMAL1), so in the
   reference cGAS-STING is scored on 3 genes and Circadian on 4; in cancer cohorts annotated with the older symbols they are scored
   on 5 / 5. This is a property of the incumbent pipeline and is reproduced, not corrected.
4. **Held-out normals (T3 population B).** Re-selecting the Note S-H cells with the paper's `nv()` gives 16,394 cells (paper: 16,312;
   scripts/E16 counts n_genes before gene de-duplication).
5. **Two-cell vs three-cell E1 on the 83,819-cell panel, bundled loadings:** 23.9 % / 27.1 %. The main-text "21 %" two-cell figure was
   computed on a different cell pool; it is not reproduced here and not needed (all T1 criteria are ratios to the incumbent).
