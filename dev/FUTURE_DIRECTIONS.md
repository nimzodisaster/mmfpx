# Future Directions

## Bootstrap diagnostic for FP model selection (`mmfp_bootstrap`)

A separate function (NOT a flag inside `mmfp`) to quantify how stable the
fractional-polynomial model selection is, addressing the overfitting/selection-
instability risk inherent in choosing one model from a large grid.

### Motivation

`mmfp` fits ~54 models and selects by AIC. The Akaike-weight columns already added
to `mmfp` show whether the selection is *identified* within a single sample (do the
data distinguish the candidates). They do NOT show whether selection is *stable*
across samples. The bootstrap answers the second question.

Key prior insight: FP candidate models are highly collinear, so the winning powers
jump around between near-identical curves across resamples. A selection-frequency
table alone is therefore misleading -- it looks alarmingly unstable even when the
fitted trajectory is rock-solid. The fix is to report the curve envelope as the
headline and the frequency table as the explanation.

### Design decisions (settled)

- **Separate function**, expensive and opt-in. Not a default-on `mmfp` feature.
- **Resample unit: subject (`id_var`), not row.** Cluster bootstrap -- resample whole
  subjects with replacement, keep all their visits. Row resampling destroys the
  longitudinal correlation structure and invalidates the results. This is a
  correctness requirement, not a preference.
- **Marginal prediction** for the curve (population trajectory, random effects
  integrated out / zeroed). Conditional-on-subject prediction is meaningless since
  bootstrap subjects differ from the original.
- **Envelope** = pointwise 2.5/97.5 quantiles of fitted curves over an age grid,
  across resamples.

### Design decisions (unresolved -- decide before building)

- **Full re-selection vs fixed-winner.** Full re-selection (re-fit the whole grid and
  re-pick by AIC on each resample) is the statistically correct answer to the
  overfitting question and is the only version that can demonstrate the
  "powers unstable but curves stable" insight. Cost is B x grid-size x per-fit.
  Fixed-winner (bootstrap only the originally selected model's curve) is much cheaper
  but measures only parameter uncertainty, not selection uncertainty. Leaning toward
  full re-selection; the cheap version may be a useful first cut for many-outcome /
  large-data use.
- **Engine scope.** Re-fit goes through existing `fit_single_model` logic, so both
  engines come nearly free. For lme4, singular-fit frequency will be high on
  resamples -- decide whether to keep singular fits in the envelope (leaning yes) and
  report the singular rate as output, consistent with how `mmfp` surfaces it.

### Proposed output (per outcome)

- Power-pair selection frequencies (the explanatory table).
- Curve envelope: tibble of `age, median, lower, upper` (the headline).
- Summary scalars: singular-fit rate, count of resamples where fitting failed entirely.

### Proposed signature (sketch)

`mmfp_bootstrap(data, outcome_vars, ..., n_boot = 200, age_grid = NULL, ...)`
passing through to the same fitting machinery as `mmfp`.

### Cost warning

Full re-selection is minutes-to-tens-of-minutes per outcome at B=1000 and must
parallelize over resamples. Default `n_boot` low (~200) with docs advising higher
for publication. Opt-in only.

## Port weighted model averaging into `normbytcv`

Replace the hard AIC-threshold equal-weighting with Akaike-weighted averaging
(`weighted.mean` over models using `exp(-0.5 * delta_AIC)` normalized weights).
Guard the single-model case (current `simplify2array %>% rowMeans` crashes when
exactly one model passes the filter). Re-document the `proportion` argument: it is
an evidence-ratio threshold, not a probability. Also port `mmfp`'s `tolower`-aware
power handling into `transform_TBV`.

## Shared FP internals

The `engine = "lm"` flag in `mmfp` is a deliberate interface tradeoff, but it adds
real interface debt. The clean long-term resolution is to extract shared
fractional-polynomial setup, transformation, fitting-grid, and summary machinery
across `mmfp`, the ordinary-`lm` path, and `normbytcv` instead of letting parallel
implementations drift.
