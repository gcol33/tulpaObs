# Every family in .tobs_family_methods must be reachable from a normal install.
#
# These run against the NAMESPACE rather than the loaded namespace, because that
# is exactly where the gap opens: tests resolve internal functions whether or not
# they are exported, so a family can be missing from NAMESPACE while its whole
# test file passes. abun() shipped that way -- a roxygen block whose `@export`
# bound to the next definition it saw, an internal helper interposed between the
# block and the constructor. `library(tulpaObs); abun()` failed for users while
# `test-abun.R` was green.

test_that("every family constructor is exported", {
  fams <- names(tulpaObs:::.tobs_family_methods)
  expect_gt(length(fams), 0)
  exports <- getNamespaceExports("tulpaObs")
  expect_setequal(intersect(fams, exports), fams)
})

test_that("no internal (dot-prefixed) name is exported", {
  # A dot-prefixed export is the signature of a roxygen block that attached to
  # the wrong definition: the block's @export lands on whatever internal helper
  # follows it. The internals are documented as .tobs_* / .dispatch_* and are
  # deliberately not part of the public surface.
  exports <- getNamespaceExports("tulpaObs")
  expect_equal(grep("^\\.", exports, value = TRUE), character(0))
})

test_that("every exported family constructor returns a tobs_family", {
  # Catches the inverse: a name exported but bound to something that is not the
  # constructor (which is how #147 presented -- .tobs_check_K_max, a validator
  # returning invisible(NULL), sat where the family object was expected).
  fams <- names(tulpaObs:::.tobs_family_methods)
  for (f in fams) {
    fn <- get(f, envir = asNamespace("tulpaObs"))
    expect_true(is.function(fn), info = f)
    obj <- tryCatch(fn(), error = function(e) e)
    if (inherits(obj, "error")) {
      # A family with a required argument (e.g. `species`) legitimately errors
      # on a bare call; it must still be a function, checked above.
      next
    }
    expect_s3_class(obj, "tobs_family")
  }
})
