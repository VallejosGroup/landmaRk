# Changelog

## landmaRk (development version)

### New Features

- [`predict_longitudinal()`](https://vallejosgroup.github.io/landmaRk/reference/predict_longitudinal.md)
  now supports custom “summary measure” functions for `method` (in the
  same spirit as `"locf"`), computed directly from the raw longitudinal
  data without requiring a prior
  [`fit_longitudinal()`](https://vallejosgroup.github.io/landmaRk/reference/fit_longitudinal.md)
  call. Such functions must have arguments `data`, `id`, `time`, `value`
  and `landmark`.

### Bug Fixes

- [`predict_longitudinal()`](https://vallejosgroup.github.io/landmaRk/reference/predict_longitudinal.md)
  no longer requires
  [`fit_longitudinal()`](https://vallejosgroup.github.io/landmaRk/reference/fit_longitudinal.md)
  to have been called when `method` is a summary measure (e.g. `"locf"`
  or a custom function); the check for an existing model fit is now only
  applied to fit-based prediction methods (e.g. `"lcmm"`, `"lme4"`).
  \`“locf”’ is implemented as a special case of summary measure

### Other Changes

- Removed `timeROC` from `Imports` as it is no longer used.
- Resolved an R CMD check NOTE about a missing global binding for `AUC`
  in
  [`performance_metrics()`](https://vallejosgroup.github.io/landmaRk/reference/performance_metrics.md).

## landmaRk 0.1.1

CRAN release: 2026-05-27

- Package is made publicly available via CRAN.

### New Features

- AUC(t) metric for evaluating landmark models.
- Cross-validation support for model evaluation.
- Subject-specific predictions with LCMM.
- Cluster assignment from LCMM fits can be used in the survival step.
- [`summary()`](https://rdrr.io/r/base/summary.html) method for
  `LandmarkAnalysis` class.
- [`prune()`](https://vallejosgroup.github.io/landmaRk/reference/prune.md)
  method for `LandmarkAnalysis` class.
- Parallelised longitudinal model fitting.
- Survival plots via
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.
- Convergence messages for fitted models.
- Print number of subjects in each risk set for
  [`show()`](https://rdrr.io/r/methods/show.html) method.

### Bug Fixes

- Fixed cross-validation not working correctly with LCMM.
- Fixed LOCF (last observation carried forward) issues.
- Fixed bug where static covariates were not used in Cox PH models.
- Fixed bug selecting IDs for survival analysis.
- Fixed bugs in working out risk sets.
- Fixed Windows parallelism issues.
- Better error handling for
  [`fit_survival()`](https://vallejosgroup.github.io/landmaRk/reference/fit_survival.md)
  when predictions are not available.
- Character covariates are now converted to factors.

### Other Changes

- Renamed package from landmarkR to landmaRk.
- Renamed main class to `LandmarkAnalysis`.
- Standardised horizon parameter handling.
- Now depends on release version of lcmm (\>= 2.2.2).
- Baseline set to 0 for all survival models.
- Internal functions refactored into smaller functions.
- Code formatted using Air.
- Added CODE_OF_CONDUCT.md.

## landmaRk 0.0.0.9000

- Package is made publicly available via r-universe and GitHub.

### Features

- Support for lmer, and hlme is added for longitudinal sub-models.
- Support for last observation carried forward (LOCF).
- Support for cox proportional hazards survival sub-model.
- C-index and Brier score metrics are implemented for evaluating
  landmark models.
- Support for arbitrary longitudinal sub-models.
