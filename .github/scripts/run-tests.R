#!/usr/bin/env Rscript

# Runs the test suite against the INSTALLED package and reports the tier it
# actually ran, not the tier it was asked for.
#
# Both CI test workflows call this: the smoke job with TULPAOBS_FAST=1, the
# weekly job with it unset. Which tier ran is the difference between "the
# plumbing holds" and "the estimators recover simulated truth", so the counts
# are logged every run and the tier is stated up front -- a smoke run reporting
# a few thousand assertions and several hundred skips looks identical to a
# broken full run unless the tier is on the record.

library(testthat)
library(tulpaObs)

fast <- identical(Sys.getenv("TULPAOBS_FAST"), "1")
tier <- if (fast) {
  "smoke (TULPAOBS_FAST=1: recovery loops skipped)"
} else {
  "full recovery (all seeds, NUTS, spatial)"
}

cat("tier      :", tier, "\n")
cat("tulpaObs  :", as.character(utils::packageVersion("tulpaObs")), "\n")
cat("tulpa     :", as.character(utils::packageVersion("tulpa")), "\n")
cat("parallel  :", Sys.getenv("TESTTHAT_PARALLEL", "<unset>"), "\n\n")

# load_package = "installed" (#151): under Config/testthat/parallel each
# worker is a separate process that must load the package itself, and without
# this argument testthat has every worker call pkgload::load_all() on the
# source tree instead -- safe for a pure-R package, but this one compiles a
# large C++ backend, so N workers each independently (re)compiling into the
# same src/ race on the shared build artifacts and corrupt the DLL. This
# script already requires the package to be installed (the library() call
# above), so telling every worker to just library(tulpaObs) -- a read of one
# already-built DLL -- costs nothing and is safe for any worker count.
res <- test_dir("tests/testthat", package = "tulpaObs",
                reporter = "summary", stop_on_failure = FALSE,
                load_package = "installed")

df <- as.data.frame(res)
failed <- sum(df$failed)
errors <- sum(df$error)

cat(sprintf(
  "\n%s\nassertions %d | skipped %d | failed %d | errors %d\n",
  tier, sum(df$passed), sum(df$skipped), failed, errors))

if (failed > 0 || errors > 0) {
  # One line per failing test rather than a printed data frame: the frame wraps
  # or drops the file column in a narrow log viewer, which is where this is read.
  bad <- df[df$failed > 0 | df$error, ]
  cat("\nFailing tests:\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  %s :: %s  (failed %d%s)\n",
                bad$file[i], bad$test[i], bad$failed[i],
                if (bad$error[i]) ", errored" else ""))
  }
  quit(status = 1L)
}

cat("Suite green.\n")
