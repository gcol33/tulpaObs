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
#
# RELEASE GATE. A release verifies what it changed; the schedule verifies
# everything. Green before tagging:
#
#   1. R CMD check on the built tarball -- a missing NAMESPACE export surfaces
#      there and nowhere else, since load_all() resolves internals regardless.
#   2. Whole-suite smoke against the INSTALLED package, which is seconds:
#        Sys.setenv(TULPAOBS_FAST = "1", NOT_CRAN = "true")
#        testthat::test_dir("tests/testthat", package = "tulpaObs",
#                           load_package = "installed")
#   3. The recovery files covering the families the diff touches, with
#      TULPAOBS_FAST unset -- a change under R/cover_*.R runs test-cover-*.R.
#      Cost is bounded by the diff, not by the suite.
#
# The full recovery tier -- every multi-seed recovery and coverage loop, all
# seeds -- does NOT gate a release. It carries the calibration evidence
# (estimators recover simulated truth; intervals cover at the nominal rate) and
# runs on the schedule in .github/workflows/full-recovery.yaml. The calibration
# claim rests on that workflow's last green run, recorded with the version it
# ran against; a release states which run it inherits. Untouched families are
# covered by the schedule, touched ones by step 3.
#
# The recovery paths call the tulpa engine, so verify the shared tulpa dep /
# Remotes pin from a clean library before tagging
# (.github/scripts/check-engine-pin.R).
skip_if_fast <- function() {
  if (identical(Sys.getenv("TULPAOBS_FAST"), "1")) {
    testthat::skip("TULPAOBS_FAST set: skipping slow recovery/coverage loop")
  }
}

# Gate the SPDE recovery suite on tulpaMesh (GitHub-only, Additional_repositories).
# tulpaMesh present -> proceed. Absent -> skip, EXCEPT when TULPAOBS_REQUIRE_SPDE
# is "1": then the whole SPDE surface is expected to run and a missing tulpaMesh
# is a loud failure, not an invisible skip. Set the flag on any CI runner that is
# supposed to exercise SPDE so a mesh-install regression cannot green-wash the
# suite by silently skipping every SPDE recovery block.
skip_if_no_tulpamesh <- function() {
  if (requireNamespace("tulpaMesh", quietly = TRUE)) return(invisible())
  if (identical(Sys.getenv("TULPAOBS_REQUIRE_SPDE"), "1")) {
    testthat::fail(paste0(
      "TULPAOBS_REQUIRE_SPDE=1 but tulpaMesh is not installed: the SPDE recovery ",
      "suite would silently skip. Install tulpaMesh (Additional_repositories) or ",
      "unset TULPAOBS_REQUIRE_SPDE."))
  }
  testthat::skip("tulpaMesh not installed: skipping SPDE recovery block")
}
