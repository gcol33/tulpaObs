# CRAN comments

## Resubmission

This is a resubmission of tulpaObs (first submitted as 0.3.0 on 2026-09-26).
In reply to the review of 2026-10-07:

* "Please write TRUE and FALSE instead of T and F. Please don't use "T" or "F"
  as vector names." -- The argument `T` of `dyn_abun_laplace()` is renamed
  `n_seasons`, and the same argument of the exported `simulate_dyn_abun()` is
  renamed `n.seasons`. No variable, argument or list element in `R/`,
  `tests/`, the examples or the vignettes is named `T` or `F`, and logical
  values are written `TRUE`/`FALSE`.

* "You have examples for unexported functions." -- `simulate_dyn_int_occu()`
  is now exported and documented. A scan of every example in `man/` finds no
  call to an unexported function.

The version is 0.6.2 rather than 0.3.0 because development continued while
the first submission was under review; NEWS.md lists the changes.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
* INLA (Suggests) is not on CRAN; Additional_repositories names its
  repository, which the check confirms is reachable.
* The incoming check also lists "MacKenzie", "Royle", "et" and "al" as
  possibly misspelled: these are the author surnames and the standard
  abbreviation from the two method references in the Description field.

## Test environments

* local: Windows 11, R 4.6.1, `R CMD check --as-cran`
* local: Windows 11, R-devel (2026-10-26 r90766 ucrt), `R CMD check
  --as-cran` on the submitted tarball: 2 NOTEs (new submission / INLA in
  Additional_repositories; HTML manual math rendering skipped because V8 is
  not installed locally), check time 301 s, tests 2248 passed.
  win-builder's R-devel queue was not processing on 2026-10-28: the tarball
  was uploaded at 12:09 CET and the queue held the same packages, none
  processed, for over two hours, so the local R-devel check stands in for it.
* GitHub Actions: ubuntu-latest (R-release) on every push; a CRAN-like job
  without `NOT_CRAN` on x86_64 and aarch64

## Notes

* tulpaObs imports and links against 'tulpa' (on CRAN, 0.7.5).

* INLA is in Suggests; two vignettes compare against it, and every chunk that
  calls it is evaluated only when it is installed; tests that need it call
  `skip_if_not_installed("INLA")`.

* Model fits in examples are wrapped in \donttest{}; the family constructors,
  simulators and formatting helpers run unwrapped. No example uses \dontrun{}.

* Recovery, coverage and sampler tests are skipped on CRAN with
  `skip_on_cran()`; they run in the package's own CI tiers.

* Thread counts default to at most two under `_R_CHECK_LIMIT_CORES_`.

* Internal helpers that seed the RNG restore the caller's `.Random.seed` on
  exit.

* The package contains a large compiled codebase (C++ likelihood kernels), so
  the installed size may exceed the default threshold on some platforms.

## Downstream dependencies

None; this is a new package.
