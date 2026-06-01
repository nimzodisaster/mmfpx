# mmfpx

**Fractional Polynomial Mixed-Effects Utilities**

[![Release](https://img.shields.io/badge/release-0.6.0-blue.svg)](https://github.com/nimzodisaster/mmfpx/releases/tag/0.6.0)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.0-276DC3.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> Tools for exploring fractional polynomial relationships in neuroimaging and related longitudinal data.

`mmfpx` helps you find and compare flexible non-linear relationships between a predictor (typically age or total brain volume) and one or more outcomes, using the fractional polynomial framework of Royston & Altman. It is built for the realities of neuroimaging data: repeated measures within subjects, many regional outcomes to model at once, and the need to defend a chosen model structure against the charge of overfitting.

---

## What it does

The package provides two complementary workflows:

**Mixed-effects model selection (`mmfp`).** Fits a grid of first- and second-order fractional polynomial models across a set of candidate powers, for one or more outcomes, and returns fit statistics that let you select an appropriate functional form. It supports three fitting backends and reports diagnostics that tell you not just *which* model won, but *how clearly* it won.

---

## Key features

### Multiple fitting engines, one interface

Choose the backend that fits your design:

| Engine | Use case | Notes |
|---|---|---|
| `"nlme"` *(default)* | Longitudinal / repeated-measures data | Supports within-subject correlation structures (`corSymm`) |
| `"lme4"` | Longitudinal data, alternative optimizer | Reports singular-fit and convergence diagnostics |
| `"lm"` | Cross-sectional / non-clustered data | Ordinary least squares; no random effects |

Arguments that don't apply to a chosen engine raise clear errors rather than being silently ignored — so you always know what you're actually fitting.

### Honest model-selection diagnostics

For every outcome, `mmfp` returns more than a single "best" model. Alongside `logLik`, `AIC`, and `BIC`, each candidate carries:

- **`delta_AIC`** — distance from the best model in the set
- **`evidence_ratio`** — how likely each model is relative to the best
- **`weight`** — Akaike weight (relative support within the candidate set)

A compact `selection_summary` reports the best model, its weight, and how many models fall within ΔAIC ≤ 2 of the top — a direct read on whether your selection is *decisive* or whether several forms fit comparably well.

> **Note:** Akaike weights describe relative support *within the fitted candidate set*. They are not probabilities that a model is correct, and they do not measure out-of-sample performance.

### Fractional polynomial flexibility

- First-order (**FP1**) and second-order (**FP2**) models, or restrict to either class
- Configurable candidate powers, including the `log` transform and correct handling of repeated-power FP2 terms
- Interactions between the fractional polynomial terms and grouping variables via a concise `static_formula` shorthand
- Optional age scaling: raw, ratio (divide by a chosen value or the mean), or standard-deviation scaling

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("nimzodisaster/mmfpx")
```

---

## Quick start

```r
library(mmfpx)

# Mixed-effects FP selection across one or more outcomes
res <- mmfp(
  data         = demo_data,
  outcome_vars = c("roi_volume"),
  age_var      = "age",
  id_var       = "subj_id",
  visit_var    = "visit"
)

# Inspect the ranked candidate models for an outcome
res$results$roi_volume$summary

# See how decisive the selection was
res$results$roi_volume$selection_summary
```

### Fit only FP1 or only FP2 models

```r
mmfp(demo_data, "roi_volume", age_var = "age", fp_models = "fp1",
     id_var = "subj_id", visit_var = "visit")
```

### Interact FP terms with a grouping variable

```r
# "*diagnosis" expands to .fp1*diagnosis (FP1) and
# .fp1*diagnosis + .fp2*diagnosis (FP2)
mmfp(demo_data, "roi_volume", age_var = "age",
     static_formula = "*diagnosis",
     id_var = "subj_id", visit_var = "visit")
```

### Cross-sectional data (no random effects)

```r
mmfp(cross_sectional_data, "roi_volume", age_var = "age", engine = "lm")
```


---

## Choosing an engine

`mmfp` defaults to `nlme` because most neuroimaging data are longitudinal, and treating repeated measures as independent invalidates inference. Use `lm` **only** for genuinely cross-sectional or non-clustered data — it treats every row as independent. Use `lme4` if you prefer its optimizer or want its singular-fit reporting; for most selection tasks `nlme` and `lme4` will agree closely on fit statistics.

All models in a selection grid are fit by maximum likelihood, so their AIC/BIC values are directly comparable across the different fractional polynomial structures.

---
## Experimental

The following functions are included but should be considered **experimental** —
their interfaces and outputs may change, and they have seen less testing than `mmfp`.

### `normbytcv` — bagged fractional-polynomial volume normalization

Adjusts regional brain volumes for their non-linear association with total brain
volume by averaging the adjustments from a set of well-supported fractional
polynomial models, rather than committing to a single fit. Models within an
evidence threshold of the best are retained and their adjusted values combined.
This is a functional but rough implementation: the model-averaging is currently
equal-weighted across the retained set rather than weighted by relative support,
and the approach assumes the volume–TBV relationship being removed is nuisance
rather than signal — an assumption worth scrutinizing for your data, since brain
regions and total volume share biologically meaningful variance.

### `brms_mmfp` — multivariate Bayesian FP model selection

Fits joint multivariate Bayesian mixed-effects models over two or more outcomes
with fractional polynomial time terms, comparing candidate bases by PSIS-LOO and
WAIC. Users supply templated fixed and random formulas containing a `{fp_terms}`
token that is expanded per candidate. This is the most experimental function in
the package: it currently selects a single winning model rather than averaging
over candidates (LOO stacking weights are the natural next step), the default
priors assume a correlated random-effects structure, and fitting is expensive —
each candidate is a full MCMC run. Treat it as a starting point for Bayesian FP
workflows rather than a finished tool.

## Roadmap

Planned for a future release: a cluster (subject-level) bootstrap diagnostic that re-runs the full selection across resamples, reporting how stable the chosen powers are and how much the fitted trajectories move — the quantitative complement to the within-sample Akaike weights. Also continued development of normbytcv and brms_mmfp, the most experimental function under development.

---

## Citation

If you use `mmfpx` in published work, please cite the package and the underlying fractional polynomial methodology (Royston & Altman, 1994).

## License

MIT © Josh Lee

---

*Version 0.6.0 — [release notes](https://github.com/nimzodisaster/mmfpx/releases/tag/0.6.0)*