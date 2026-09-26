# CRAN comments

## New submission

This is the first CRAN submission of tulpaObs. It fits occupancy, abundance and
detection models by supplying their observation likelihoods to the 'tulpa'
engine (tulpa 0.6.0, on CRAN since 2026-09-24), which it imports and links
against.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
* INLA (Suggests) is not on CRAN; Additional_repositories names its
  repository, which the check confirms is reachable.
* The incoming check on win-builder r-devel also lists "MacKenzie", "Royle", "et"
  and "al" as possibly misspelled: these are the author surnames and the
  standard abbreviation from the two method references in the Description
  field, not software names, so they are left unquoted.

## Test environments

* local: Windows 11, R 4.6.1, `R CMD check --as-cran` including the PDF manual
* win-builder: R-devel (Status: 1 NOTE, the one above)
* GitHub Actions: ubuntu-latest (R-release), windows-latest and macos-latest
  on dispatch

## Notes

* INLA is in Suggests and is not on CRAN; Additional_repositories names its
  repository. Two vignettes compare against it, and every chunk that calls it
  is evaluated only when it is installed; tests that need it call
  `skip_if_not_installed("INLA")`.

* Model fits in examples are wrapped in \donttest{}; the family constructors,
  simulators and formatting helpers run unwrapped. No example uses \dontrun{}.

* Recovery, coverage and sampler tests are skipped on CRAN with
  `skip_on_cran()`; they run in the package's own CI tiers.

* The testthat run under `--as-cran` takes about two minutes and reports no
  failures; recovery, coverage and sampler tests are skipped on CRAN.

* Thread counts default to at most two under `_R_CHECK_LIMIT_CORES_`.

* Internal helpers that seed the RNG restore the caller's `.Random.seed` on
  exit.

* The package contains a large compiled codebase (C++ likelihood kernels), so
  the installed size may exceed the default threshold on some platforms.

## Downstream dependencies

None; this is a new package.
