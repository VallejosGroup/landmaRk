## R CMD check results

0 errors | 0 warnings | 3 notes

* "checking CRAN incoming feasibility ... Days since last update: 5". This
  update is submitted shortly after 0.1.2 because it fixes an important bug:
  landmaRk can now correctly handle more than one dynamic (longitudinal)
  covariate. The previous release did not handle this case correctly.
* "checking for future file timestamps ... unable to verify current time" is a
  local network artifact and unrelated to the package.
* "checking HTML version of manual" reports `<main> is not recognized` errors.
  These come from an outdated local HTML Tidy that does not know the HTML5
  `<main>` element inserted into help pages; they do not appear on CRAN's own
  check machines.

## This submission

This is an update (version 0.1.3). Key change since the previous CRAN release:

* Fixed an important bug so that landmaRk can now handle more than one dynamic
  (longitudinal) covariate.

## Test environments

- windows-latest (R release)
- macOS-latest (R release)
- macOS-latest (R devel)
- ubuntu-latest (R release)
- ubuntu-latest (R oldrel-1)
- ubuntu-latest (R oldrel-2)

via GitHub actions https://github.com/vallejosgroup/landmaRk/actions
