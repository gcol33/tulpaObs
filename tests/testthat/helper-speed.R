# Speed-tier control for the test suite.
#
# The parameter-recovery and coverage tests fit 15-20 models per seed across
# many seeds (and the NUTS tests run the sampler), which dominates the suite's
# wall-clock. For fast dev iteration, set the TULPAOBS_FAST environment variable
# to "1" -- the heavy statistical loops then skip (reported as skips, never
# silently dropped) and only the structural / dispatch / closed-form unit tests
# run, finishing in well under a second.
#
#   Sys.setenv(TULPAOBS_FAST = "1"); devtools::test()   # fast: structure only
#   Sys.unsetenv("TULPAOBS_FAST");   devtools::test()   # full: + recovery loops
#
# The default (variable unset) runs everything, so CRAN and CI always see the
# full recovery suite. Place skip_if_fast() as the first line of any test_that()
# block whose cost is a multi-seed fitting loop or a NUTS sample.
skip_if_fast <- function() {
  if (identical(Sys.getenv("TULPAOBS_FAST"), "1")) {
    testthat::skip("TULPAOBS_FAST set: skipping slow recovery/coverage loop")
  }
}
