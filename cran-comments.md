## R CMD check results

0 errors | 0 warnings | 1 note

* The only NOTE ("checking for future file timestamps ... unable to verify
  current time") is a local network artifact and unrelated to the package.

## This submission

This is an update (version 0.1.2). Changes since the previous CRAN release:

* `predict_longitudinal()` now supports custom "summary measure" functions for
  `method` (in the same spirit as `"locf"`), computed directly from the raw
  longitudinal data without requiring a prior `fit_longitudinal()` call.
* `predict_longitudinal()` no longer requires `fit_longitudinal()` to have been
  called when `method` is a summary measure.
* Fixed the new implementation of LOCF (last observation carried forward).
* Removed `timeROC` from `Imports` as it is no longer used.
* Resolved an R CMD check NOTE about a missing global binding for `AUC` in
  `performance_metrics()`.

## Test environments

- windows-latest (R release)
- macOS-latest (R release)
- macOS-latest (R devel)
- ubuntu-latest (R release)
- ubuntu-latest (R oldrel-1)
- ubuntu-latest (R oldrel-2)

via GitHub actions https://github.com/vallejosgroup/landmaRk/actions
