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
#
# Two environment variables shape a run:
#
#   TULPAOBS_SHARD     one shard id from .github/scripts/shard-tests.R. Unset
#                      (the smoke job, and any local run) means the whole
#                      directory, exactly as before.
#   TULPAOBS_TEST_OUT  directory to write the machine-readable results into,
#                      for the aggregate job to collect. Unset writes nothing.

library(testthat)
library(tulpaObs)

fast <- identical(Sys.getenv("TULPAOBS_FAST"), "1")
tier <- if (fast) {
  "smoke (TULPAOBS_FAST=1: recovery loops skipped)"
} else {
  "full recovery (all seeds, NUTS, spatial)"
}

shard <- Sys.getenv("TULPAOBS_SHARD", "")
out_dir <- Sys.getenv("TULPAOBS_TEST_OUT", "")

filter <- NULL
assigned <- NULL
if (nzchar(shard)) {
  source(".github/scripts/shard-tests.R")
  plan <- tier3_plan()
  assigned <- tier3_shard_files(shard, plan)
  filter <- tier3_filter_regex(assigned)
  # Confirm the pattern selects exactly the assigned files BEFORE spending the
  # run on the assumption. A filter that quietly selects the wrong set is the
  # one failure mode of sharding that looks green.
  tier3_check_filter(filter, assigned, plan$file)
}

cat("tier      :", tier, "\n")
cat("tulpaObs  :", as.character(utils::packageVersion("tulpaObs")), "\n")
cat("tulpa     :", as.character(utils::packageVersion("tulpa")), "\n")
cat("parallel  :", Sys.getenv("TESTTHAT_PARALLEL", "<unset>"), "\n")
if (nzchar(shard)) {
  cat(sprintf("shard     : %s (%d of %d test files)\n", shard,
              length(assigned), length(plan$file)))
  cat("files     :", paste(assigned, collapse = ", "), "\n")
} else {
  cat("shard     : <none: whole directory>\n")
}
cat("\n")

started <- Sys.time()

# load_package = "installed" (#151): under Config/testthat/parallel each
# worker is a separate process that must load the package itself, and without
# this argument testthat has every worker call pkgload::load_all() on the
# source tree instead -- safe for a pure-R package, but this one compiles a
# large C++ backend, so N workers each independently (re)compiling into the
# same src/ race on the shared build artifacts and corrupt the DLL. This
# script already requires the package to be installed (the library() call
# above), so telling every worker to just library(tulpaObs) -- a read of one
# already-built DLL -- costs nothing and is safe for any worker count.
res <- test_dir("tests/testthat", package = "tulpaObs", filter = filter,
                reporter = "summary", stop_on_failure = FALSE,
                load_package = "installed")

elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

df <- as.data.frame(res)
failed <- sum(df$failed)
errors <- sum(df$error)

# Per-file wall time, summed over the file's tests. Under parallel workers this
# is each test's own elapsed time rather than the file's share of the job, which
# is what the shard planner wants: a file's serial cost, independent of how many
# workers happened to run alongside it.
by_file <- data.frame(file = character(0), seconds = numeric(0),
                      assertions = integer(0), skipped = integer(0),
                      failed = integer(0), errors = integer(0),
                      stringsAsFactors = FALSE)
ran <- character(0)
if (nrow(df)) {
  df$file <- basename(as.character(df$file))
  ran <- sort(unique(df$file))
  split_by <- factor(df$file, levels = ran)
  num <- function(col) as.numeric(tapply(col, split_by, sum))
  by_file <- data.frame(
    file = ran,
    seconds = round(num(df$real), 3),
    assertions = as.integer(num(df$passed)),
    skipped = as.integer(num(as.integer(df$skipped))),
    failed = as.integer(num(df$failed)),
    errors = as.integer(num(as.integer(df$error))),
    stringsAsFactors = FALSE)
  by_file <- by_file[order(-by_file$seconds), , drop = FALSE]
}

if (nzchar(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(by_file, file.path(out_dir, "timings.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    shard = if (nzchar(shard)) shard else "all",
    tier = tier,
    files_assigned = if (is.null(assigned)) NA_integer_ else length(assigned),
    files_ran = length(ran),
    assertions = sum(df$passed),
    skipped = sum(df$skipped),
    failed = failed,
    errors = errors,
    seconds = round(elapsed, 1),
    stringsAsFactors = FALSE),
    file.path(out_dir, "summary.csv"), row.names = FALSE)
}

cat(sprintf(
  "\n%s\nassertions %d | skipped %d | failed %d | errors %d | %.1f min\n",
  tier, sum(df$passed), sum(df$skipped), failed, errors, elapsed / 60))

if (nrow(by_file)) {
  cat("\nSlowest files:\n")
  for (i in seq_len(min(10L, nrow(by_file)))) {
    cat(sprintf("  %8.1fs  %s\n", by_file$seconds[i], by_file$file[i]))
  }
}

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
}

# A shard must run the files it was given, all of them and nothing else. The
# pre-flight above checks the pattern; this checks the outcome, so a file that
# was selected but produced no results at all still surfaces instead of being
# counted as covered by a green shard.
if (!is.null(assigned)) {
  dropped <- setdiff(assigned, ran)
  extra <- setdiff(ran, assigned)
  if (length(dropped) || length(extra)) {
    cat("\nShard coverage FAILED for '", shard, "':\n", sep = "")
    if (length(dropped)) {
      cat("  assigned but produced no results: ",
          paste(dropped, collapse = ", "), "\n", sep = "")
    }
    if (length(extra)) {
      cat("  ran but not assigned: ", paste(extra, collapse = ", "), "\n", sep = "")
    }
    quit(status = 1L)
  }
}

if (failed > 0 || errors > 0) quit(status = 1L)

cat("Suite green.\n")
