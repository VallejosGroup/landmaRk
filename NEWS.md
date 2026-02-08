# landmaRk

# landmaRk 0.1.1

* Package is made publicly available via CRAN.

## New Features

* AUC(t) metric for evaluating landmark models.
* Cross-validation support for model evaluation.
* Subject-specific predictions with LCMM.
* Cluster assignment from LCMM fits can be used in the survival step.
* `summary()` method for `LandmarkAnalysis` class.
* `prune()` method for `LandmarkAnalysis` class.
* Parallelised longitudinal model fitting.
* Survival plots via `plot()` method.
* Convergence messages for fitted models.
* Print number of subjects in each risk set for `show()` method.

## Bug Fixes

* Fixed cross-validation not working correctly with LCMM.
* Fixed LOCF (last observation carried forward) issues.
* Fixed bug where static covariates were not used in Cox PH models.
* Fixed bug selecting IDs for survival analysis.
* Fixed bugs in working out risk sets.
* Fixed Windows parallelism issues.
* Better error handling for `fit_survival()` when predictions are not available.
* Character covariates are now converted to factors.

## Other Changes

* Renamed package from landmarkR to landmaRk.
* Renamed main class to `LandmarkAnalysis`.
* Standardised horizon parameter handling.
* Now depends on release version of lcmm (>= 2.2.2).
* Baseline set to 0 for all survival models.
* Internal functions refactored into smaller functions.
* Code formatted using Air.
* Added CODE_OF_CONDUCT.md.

# landmaRk 0.0.0.9000

* Package is made publicly available via r-universe and GitHub. 

## Features

* Support for lmer, and hlme is added for longitudinal sub-models. 
* Support for last observation carried forward (LOCF). 
* Support for cox proportional hazards survival sub-model. 
* C-index and Brier score metrics are implemented for evaluating landmark models.
* Support for arbitrary longitudinal sub-models.
