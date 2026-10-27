# T6 — blinded nameability rubric (frozen 2026-09-19, before any rating)

**Raters.** Three independent rater instances (LLM agents, fresh context, no access to this repository, the manuscript, or any
other axis). Each rater receives ONE sheet: for every axis of ONE representation, the top-5 modules by |loading| with their
signed loadings, and nothing else (no representation name, no PC index semantics, no cancer context). Sheets are shuffled and
axis labels anonymised (`axis_01 … axis_09`); representation identity is hidden behind a random code.

**Task per axis.** Return (a) a one-phrase biological name for the axis (e.g. "inflammatory ↔ steroid-metabolic"),
(b) a confidence 1–5 that the name is a faithful, specific summary of the listed modules:
  1 = no coherent theme, name is a guess · 2 = weak theme · 3 = theme present but with clear outliers ·
  4 = clear theme, minor outliers · 5 = unambiguous, every listed module fits the name.

**Score.** Interp-N(representation) = mean confidence over its axes and raters. Inter-rater agreement = Kendall's W over axes
(confidence) and mean pairwise Spearman ρ across raters. Raters may not confer.

**Comparison.** RETIRE requires Interp-N(rotation) − Interp-N(raw PCA) ≥ 0.5 (pre-reg T6), with Interp-S also above raw PCA,
and E1 within 2 pp / cytotrophoblast E2 rank ≤ 3 under that rotation.
