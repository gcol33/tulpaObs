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
# Three environment variables shape a run:
#
#   TULPAOBS_SHARD     one shard id from .github/scripts/shard-tests.R. Unset
#                      (the smoke job, and any local run) means the whole
#                      directory, exactly as before.
#   TULPAOBS_FILES     the subset of the suite to run: file names, testthat
#                      stems or globs, comma or space separated. Read whether
#                      or not a shard is set, so a local run and a dispatched
#                      one select the same way. A pattern matching nothing is
#                      an error, not an empty run.
#   TULPAOBS_TEST_OUT  directory to write the machine-readable results into,
#                      for the aggregate job to collect. Unset writes nothing.
#   TULPAOBS_WORKERS   worker count. Unset takes every core under CI, and a
#                      fraction of them elsewhere so a shared workstation keeps
#                      a share of itself.

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

test_path <- "tests/testthat"
assigned <- NULL
if (nzchar(shard)) {
  source(".github/scripts/shard-tests.R")
  plan <- tier3_plan()
  assigned <- tier3_shard_files(shard, plan)
  # Longest known cost first. A pool's wall time is set by when its longest file
  # starts, so a heavy file picked up last leaves every other worker idle while
  # it finishes.
  files <- assigned[order(-plan$seconds[match(assigned, plan$file)], assigned)]
} else {
  files <- sort(basename(list.files(test_path, pattern = "^test-.*[.][Rr]$")))
  # A shardless run honours the same selection a dispatched one does; the
  # planner owns the matching so the two cannot read a pattern differently.
  if (nzchar(Sys.getenv("TULPAOBS_FILES", ""))) {
    source(".github/scripts/shard-tests.R")
    files <- tier3_selection(files)
  }
}

# One test file per worker, and the worker times its own file. testthat's own
# parallel runner is not used here, and neither is its per-test `real` column:
# that column is populated only for files run in this process, so under workers
# it comes back 0 for every test and a shard that ran for four hours publishes
# per-file costs summing to one second. Those costs feed the weights table the
# shard planner packs against, so a tier that cannot be measured is a tier that
# can only ever be packed against guesses. Owning the dispatch puts a clock
# around each file rather than around each test, which is both measurable under
# parallelism and the number the planner actually wants: a file's cost on its
# own, not its share of a job it happened to share with others.
n_workers <- local({
  env <- suppressWarnings(as.integer(Sys.getenv("TULPAOBS_WORKERS", "")))
  n <- if (!is.na(env) && env >= 1L) {
    env
  } else {
    cores <- parallel::detectCores()
    # A worker is not one thread: the fits inside it run OpenMP, so worker
    # count and load are not proportional. Measured on the 32-processor box,
    # 16 workers took it to 94% busy and 14 left about 30% free, which is the
    # share the machine's owner is meant to keep -- hence a fraction well under
    # one half rather than one core short of all of them. CI has no owner to
    # leave room for.
    if (is.na(cores)) 1L else if (nzchar(Sys.getenv("CI"))) cores else
      max(1L, floor(0.45 * cores))
  }
  max(1L, min(n, length(files)))
})

cat("tier      :", tier, "\n")
cat("tulpaObs  :", as.character(utils::packageVersion("tulpaObs")), "\n")
cat("tulpa     :", as.character(utils::packageVersion("tulpa")), "\n")
cat("workers   :", n_workers, "\n")
if (nzchar(shard)) {
  cat(sprintf("shard     : %s (%d of %d test files)\n", shard,
              length(assigned), length(plan$file)))
  cat("files     :", paste(sort(assigned), collapse = ", "), "\n")
} else {
  cat(sprintf("shard     : <none: whole directory, %d files>\n", length(files)))
}
cat("\n")

# load_package = "installed" (#151): a worker that loads the source tree instead
# runs pkgload::load_all(), and this package compiles a large C++ backend, so N
# workers each independently (re)compiling into the same src/ race on the shared
# build artifacts and corrupt the DLL. This script already requires the package
# to be installed (the library() call above), so telling every worker to just
# library(tulpaObs) -- a read of one already-built DLL -- costs nothing and is
# safe at any worker count.
run_one <- function(file) {
  started <- Sys.time()
  cat(sprintf("[start] %s\n", file))
  utils::flush.console()

  crashed <- NULL
  res <- NULL
  text <- utils::capture.output(
    res <- tryCatch(
      testthat::test_file(file.path("tests/testthat", file),
                          package = "tulpaObs",
                          load_package = "installed",
                          reporter = "summary"),
      error = function(e) {
        crashed <<- conditionMessage(e)
        NULL
      }))

  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  df <- if (is.null(res)) NULL else as.data.frame(res)
  if (!is.null(df) && nrow(df)) df$file <- file

  cat(sprintf("[done ] %8.1fs  %s%s\n", seconds, file,
              if (is.null(crashed)) "" else "  (ERRORED OUT)"))
  utils::flush.console()

  list(file = file, seconds = seconds, df = df, text = text, crashed = crashed)
}

started <- Sys.time()

# outfile = "" forwards each worker's own output to this console, so the
# [start]/[done] lines above arrive live rather than at the end of a job that
# can run for hours. The test output itself is captured inside the worker and
# replayed below, which keeps one file's failures in one block instead of
# interleaved with whatever else was running at the time.
if (n_workers > 1L) {
  cl <- parallel::makeCluster(n_workers, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, "run_one", envir = environment())
  parallel::clusterCall(cl, function(wd) {
    setwd(wd)
    library(testthat)
    library(tulpaObs)
    invisible(NULL)
  }, wd = getwd())
  results <- parallel::clusterApplyLB(cl, files, function(f) run_one(f))
} else {
  results <- lapply(files, run_one)
}

elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

crashed <- vapply(results, function(r) !is.null(r$crashed), logical(1))
df <- do.call(rbind, lapply(results, `[[`, "df"))
if (is.null(df)) {
  df <- data.frame(file = character(0), test = character(0), passed = numeric(0),
                   skipped = logical(0), failed = numeric(0), error = logical(0),
                   stringsAsFactors = FALSE)
}

failed <- sum(df$failed)
errors <- sum(df$error) + sum(crashed)

# Per-file wall time, measured by the worker that ran the file, alone.
by_file <- data.frame(
  file = vapply(results, `[[`, character(1), "file"),
  seconds = round(vapply(results, `[[`, numeric(1), "seconds"), 3),
  assertions = 0L, skipped = 0L, failed = 0L, errors = 0L,
  stringsAsFactors = FALSE)
if (nrow(df)) {
  tally <- function(col) {
    v <- tapply(as.numeric(col), factor(df$file, levels = by_file$file), sum)
    as.integer(ifelse(is.na(v), 0L, v))
  }
  by_file$assertions <- tally(df$passed)
  by_file$skipped <- tally(as.integer(df$skipped))
  by_file$failed <- tally(df$failed)
  by_file$errors <- tally(as.integer(df$error))
}
by_file$errors <- by_file$errors + as.integer(crashed)
ran <- sort(by_file$file[by_file$assertions > 0L | by_file$skipped > 0L |
                           by_file$failed > 0L | by_file$errors > 0L])
by_file <- by_file[order(-by_file$seconds), , drop = FALSE]

if (nzchar(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(by_file, file.path(out_dir, "timings.csv"), row.names = FALSE)
  # Per-test durations as well as per-file. A worker runs its file in this
  # process, so testthat's `real` column is populated here even though it is
  # not under testthat's own parallel runner -- and an expensive file is worth
  # nothing to the planner until you can see which block inside it is the cost.
  if (nrow(df)) {
    utils::write.csv(
      data.frame(file = df$file, test = df$test,
                 seconds = round(as.numeric(df$real), 3),
                 assertions = as.integer(df$passed),
                 skipped = as.integer(df$skipped),
                 failed = as.integer(df$failed),
                 errors = as.integer(df$error),
                 stringsAsFactors = FALSE),
      file.path(out_dir, "tests.csv"), row.names = FALSE)
  }
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

# The captured output, one block per file, for the files that need reading. A
# green file's summary block says nothing its row in the table below does not.
for (r in results) {
  bad <- !is.null(r$crashed) ||
    (!is.null(r$df) && nrow(r$df) &&
       (sum(r$df$failed) > 0 || any(r$df$error)))
  if (!bad) next
  cat("\n", strrep("-", 70), "\n", r$file, "\n", strrep("-", 70), "\n", sep = "")
  cat(r$text, sep = "\n")
  if (!is.null(r$crashed)) {
    cat("\nfile did not return results: ", r$crashed, "\n", sep = "")
  }
}

cat(sprintf(paste0("\n%s\nassertions %d | skipped %d | failed %d | errors %d",
                   " | %.1f min wall, %.1f min summed over files\n"),
            tier, sum(df$passed), sum(df$skipped), failed, errors,
            elapsed / 60, sum(by_file$seconds) / 60))

if (nrow(by_file)) {
  cat("\nSlowest files:\n")
  for (i in seq_len(min(10L, nrow(by_file)))) {
    cat(sprintf("  %8.1fs  %s\n", by_file$seconds[i], by_file$file[i]))
  }
}

if (failed > 0 || errors > 0) {
  # One line per failing test rather than a printed data frame: the frame wraps
  # or drops the file column in a narrow log viewer, which is where this is read.
  bad <- df[df$failed > 0 | df$error, , drop = FALSE]
  cat("\nFailing tests:\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  %s :: %s  (failed %d%s)\n",
                bad$file[i], bad$test[i], bad$failed[i],
                if (bad$error[i]) ", errored" else ""))
  }
  for (r in results[crashed]) {
    cat(sprintf("  %s :: <whole file>  (errored: %s)\n", r$file, r$crashed))
  }
}

# A shard must run the files it was given, all of them and nothing else. A file
# that was dispatched but produced no results at all surfaces here instead of
# being counted as covered by a green shard.
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

# A file that stops RUNNING its tests reports the same clean tail as one that
# ran them: testthat prints a skip only where a skip happened, and a block
# gated by both skip_on_cran() and skip_if_fast() is reached by neither routine
# tier, so an assertion can go stale or a whole block can stop executing with
# nothing to notice (gcol33/tulpaObs#302). The counts a file HAS produced are
# recorded in tests/expected-counts.csv; this compares against them so a run
# has to have run, rather than merely not failed.
#
# A smoke run's skip counts differ from a full run's by construction, so the
# manifest carries one row per (file, tier) and a run reads its own tier's
# rows. Files absent from the manifest are not checked -- it records what has
# actually been measured, and the full tier has never completed once, so it
# starts partial and grows.
drift <- NULL
manifest_path <- "tests/expected-counts.csv"
if (file.exists(manifest_path) && nrow(by_file)) {
  want <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  want <- want[want$tier == if (fast) "smoke" else "full", , drop = FALSE]
  got  <- by_file[match(want$file, by_file$file), , drop = FALSE]
  seen <- !is.na(got$file)
  want <- want[seen, , drop = FALSE]
  got  <- got[seen, , drop = FALSE]
  if (nrow(want)) {
    lost    <- got$assertions < want$assertions
    skipped <- got$skipped    > want$skipped
    if (any(lost | skipped)) {
      drift <- data.frame(
        file = want$file[lost | skipped],
        assertions = sprintf("%d (expected %d)", got$assertions[lost | skipped],
                             want$assertions[lost | skipped]),
        skipped = sprintf("%d (expected %d)", got$skipped[lost | skipped],
                          want$skipped[lost | skipped]),
        stringsAsFactors = FALSE)
      cat("\nCount check FAILED -- a file ran fewer assertions, or skipped more,\n",
          "than tests/expected-counts.csv records for the ",
          if (fast) "smoke" else "full", " tier:\n", sep = "")
      for (i in seq_len(nrow(drift))) {
        cat(sprintf("  %-46s assertions %-18s skipped %s\n",
                    drift$file[i], drift$assertions[i], drift$skipped[i]))
      }
      cat("\nIf the drop is intended (a block deleted or deliberately gated),",
          "update the manifest in the same commit.\n")
    }
    gained <- got$assertions > want$assertions
    if (any(gained)) {
      cat("\nCount check: ", sum(gained), " file(s) ran MORE assertions than ",
          "recorded (manifest is behind, not a failure):\n", sep = "")
      for (i in which(gained)) {
        cat(sprintf("  %-46s %d, recorded %d\n", want$file[i],
                    got$assertions[i], want$assertions[i]))
      }
    }
    cat(sprintf("\ncount check: %d of %d %s-tier manifest files ran here\n",
                nrow(want), sum(utils::read.csv(manifest_path)$tier ==
                                  (if (fast) "smoke" else "full")),
                if (fast) "smoke" else "full"))
  }
}

if (failed > 0 || errors > 0 || !is.null(drift)) quit(status = 1L)

cat("Suite green.\n")
