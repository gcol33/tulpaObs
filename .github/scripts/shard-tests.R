#!/usr/bin/env Rscript

# Shard planner for the full-recovery (tier 3) workflow.
#
# The tier's cost is concentrated in a handful of files -- one measured at 5.12h
# against a 350-minute job cap -- so a single job carrying the whole suite is
# cancelled before it reports anything, and every recovery and coverage number
# in the package stays unverified by CI. Splitting the suite across jobs removes
# that: no job carries the cap, and a job that does overrun costs only the files
# assigned to it.
#
# Two kinds of shard:
#
#   solo-<stem>   one file, chosen because its measured cost is a large fraction
#                 of the cap or because it has never completed a timed run. On
#                 its own it either fits or it does not, and either way the rest
#                 of the tier still reports.
#   pool-<n>      everything else, bin-packed longest-first across N buckets.
#
# Assignment is computed here and nowhere else: the workflow's matrix, each
# shard job's file filter, and the aggregate job's coverage check all call this
# script, so they cannot disagree about which file belongs where.
#
# NOTHING MAY BE DROPPED. A shard scheme that quietly skips files is worse than
# a slow one, so this script fails rather than proceeds when the assignment does
# not account for every file, and run-tests.R separately checks that the files
# that actually ran are the files the shard was assigned.
#
# Usage:
#   Rscript .github/scripts/shard-tests.R ids           JSON array of shard ids
#   Rscript .github/scripts/shard-tests.R summary       human-readable plan
#   Rscript .github/scripts/shard-tests.R plan          shard<TAB>file per line
#   Rscript .github/scripts/shard-tests.R files <id>    the files of one shard
#
# Sourcing it with no arguments defines the functions and runs nothing.

# Number of buckets the non-isolated files are packed into. Higher costs only
# runner concurrency (the repository is public, so job minutes are free) and
# buys resilience: an unmeasured expensive file that overruns takes out 1/N of
# the pooled evidence rather than all of it.
#
# Ten, because the pooled shards overran what their weights said they held. At
# six pools each carried 195.5 min of weight and ran on several workers, so no
# pool should have come near a 350-minute cap; pool-6 reached it and reported
# nothing, and pool-1 came back at 81% of it. The weights cannot explain that,
# which is the point: roughly 200 files carry the placeholder below rather than
# a measurement, and the 21 SBC acceptance files carry an equal share of one
# total rather than a per-file cost. Until a completed run replaces those, the
# only lever that lowers the load a single shard carries is the number of
# shards, and it costs nothing but concurrency.
TIER3_POOLS <- 10L

# Weight for a file with no measurement, used only to spread the unmeasured
# majority evenly. Order of magnitude, from the numbers already on record: the
# tier is documented at roughly 9h serial and the eight measured files account
# for 7.98h of that, leaving about an hour across the remaining ~200 files. It
# is a placeholder for real per-file timings, which the aggregate job publishes.
TIER3_DEFAULT_SECONDS <- 20

# testthat's own pattern for what counts as a test file, so this planner and
# testthat enumerate the same set.
TIER3_TEST_PATTERN <- "^test.*\\.[rR]$"

tier3_test_files <- function(test_path = "tests/testthat") {
  if (!dir.exists(test_path)) {
    stop("test directory not found: ", test_path, call. = FALSE)
  }
  files <- sort(list.files(test_path, pattern = TIER3_TEST_PATTERN))
  if (!length(files)) stop("no test files under ", test_path, call. = FALSE)
  files
}

# testthat matches a filter against the file name with the leading "test-" and
# the extension removed (testthat:::context_name).
tier3_stem <- function(files) {
  sub("[.][Rr]$", "", sub("^test[-_]", "", files))
}

tier3_weights <- function(weights_path = ".github/tier3-weights.csv") {
  if (!file.exists(weights_path)) {
    stop("weights table not found: ", weights_path, call. = FALSE)
  }
  w <- utils::read.csv(weights_path, comment.char = "#", stringsAsFactors = FALSE)
  need <- c("file", "seconds", "isolate")
  miss <- setdiff(need, names(w))
  if (length(miss)) {
    stop("weights table is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  w$file <- trimws(w$file)
  w$seconds <- suppressWarnings(as.numeric(w$seconds))
  iso <- tolower(trimws(as.character(w$isolate)))
  if (!all(iso %in% c("yes", "no"))) {
    stop("weights table column 'isolate' takes yes or no, got: ",
         paste(unique(iso[!iso %in% c("yes", "no")]), collapse = ", "),
         call. = FALSE)
  }
  w$isolate <- iso == "yes"
  if (anyDuplicated(w$file)) {
    stop("weights table lists a file twice: ",
         paste(unique(w$file[duplicated(w$file)]), collapse = ", "), call. = FALSE)
  }
  bad <- w$file[!w$isolate & is.na(w$seconds)]
  if (length(bad)) {
    stop("weights table gives no seconds and does not isolate: ",
         paste(bad, collapse = ", "),
         ". An unmeasured file either takes the default weight (drop the row) ",
         "or gets a shard to itself (isolate = yes).", call. = FALSE)
  }
  w
}

# Longest-processing-time-first packing: the standard makespan heuristic, and
# the only sensible one here because the cost distribution is extremely skewed.
# Ties break on file name so a plan is reproducible across runs.
tier3_plan <- function(test_path = "tests/testthat",
                       weights_path = ".github/tier3-weights.csv",
                       n_pools = TIER3_POOLS) {
  files <- tier3_test_files(test_path)
  w <- tier3_weights(weights_path)

  stale <- setdiff(w$file, files)
  if (length(stale)) {
    stop("weights table names file(s) that do not exist: ",
         paste(stale, collapse = ", "),
         ". Update ", weights_path, " -- a stale weight silently mis-packs the ",
         "shard it was meant to describe.", call. = FALSE)
  }

  seconds <- rep(TIER3_DEFAULT_SECONDS, length(files))
  measured <- rep(FALSE, length(files))
  isolate <- rep(FALSE, length(files))
  names(seconds) <- names(measured) <- names(isolate) <- files

  hit <- match(w$file, files)
  seconds[hit] <- ifelse(is.na(w$seconds), TIER3_DEFAULT_SECONDS, w$seconds)
  measured[hit] <- !is.na(w$seconds)
  isolate[hit] <- w$isolate

  shard <- rep(NA_character_, length(files))
  names(shard) <- files

  solo <- files[isolate]
  shard[solo] <- paste0("solo-", tier3_stem(solo))

  pooled <- files[!isolate]
  pools <- paste0("pool-", seq_len(n_pools))
  load <- stats::setNames(numeric(length(pools)), pools)
  ord <- pooled[order(-seconds[pooled], pooled)]
  for (f in ord) {
    k <- which.min(load)
    shard[f] <- pools[k]
    load[k] <- load[k] + seconds[f]
  }

  out <- data.frame(file = files, shard = unname(shard[files]),
                    seconds = unname(seconds[files]),
                    measured = unname(measured[files]),
                    isolate = unname(isolate[files]),
                    stringsAsFactors = FALSE)

  # Every file assigned exactly once, to a shard that exists.
  if (anyNA(out$shard) || !all(nzchar(out$shard))) {
    stop("planner left ", sum(is.na(out$shard) | !nzchar(out$shard)),
         " file(s) unassigned.", call. = FALSE)
  }
  if (!setequal(out$file, files) || nrow(out) != length(files)) {
    stop("planner assignment does not cover the test directory exactly.",
         call. = FALSE)
  }
  out
}

# Shard ids in a stable order, empty pools dropped: an empty filter makes
# test_dir abort with "No test files found", which would read as a broken shard
# rather than an unused one.
tier3_shard_ids <- function(plan = tier3_plan()) {
  ids <- unique(plan$shard)
  solo <- sort(ids[startsWith(ids, "solo-")])
  pool <- ids[startsWith(ids, "pool-")]
  pool <- pool[order(as.integer(sub("^pool-", "", pool)))]
  c(solo, pool)
}

tier3_shard_files <- function(shard, plan = tier3_plan()) {
  files <- plan$file[plan$shard == shard]
  if (!length(files)) {
    stop("no files assigned to shard '", shard, "'. Known shards: ",
         paste(tier3_shard_ids(plan), collapse = ", "), call. = FALSE)
  }
  files
}


tier3_json_array <- function(x) {
  paste0("[", paste0("\"", x, "\"", collapse = ","), "]")
}

tier3_summary <- function(plan = tier3_plan()) {
  ids <- tier3_shard_ids(plan)
  cat(sprintf("test files : %d\n", nrow(plan)))
  cat(sprintf("measured   : %d\n", sum(plan$measured)))
  cat(sprintf("shards     : %d (%d solo, %d pooled)\n\n", length(ids),
              sum(startsWith(ids, "solo-")), sum(startsWith(ids, "pool-"))))
  cat(sprintf("%-28s %6s %10s  %s\n", "shard", "files", "est (min)", "note"))
  for (id in ids) {
    sub <- plan[plan$shard == id, , drop = FALSE]
    note <- if (any(!sub$measured)) {
      sprintf("%d of %d unmeasured", sum(!sub$measured), nrow(sub))
    } else {
      "all measured"
    }
    # An estimate built entirely out of the placeholder weight is not an
    # estimate; print nothing rather than a number that reads as one.
    est <- if (any(sub$measured)) sprintf("%10.1f", sum(sub$seconds) / 60) else
      sprintf("%10s", "-")
    cat(sprintf("%-28s %6d %s  %s\n", id, nrow(sub), est, note))
  }
  invisible(plan)
}

# Command-line dispatch, only when this file is the script being run. The
# consumers source it, and they take arguments of their own -- keying off
# trailing arguments alone would make aggregate-tier3.R's artifact directory
# look like a mode name.
local({
  invoked <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(invoked) || !identical(basename(invoked[[1]]), "shard-tests.R")) {
    return(invisible(NULL))
  }
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    stop("usage: shard-tests.R ids | summary | plan | files <shard>", call. = FALSE)
  }
  mode <- args[[1]]
  if (identical(mode, "ids")) {
    cat(tier3_json_array(tier3_shard_ids()), "\n", sep = "")
  } else if (identical(mode, "summary")) {
    tier3_summary()
  } else if (identical(mode, "plan")) {
    p <- tier3_plan()
    cat(paste0(p$shard, "\t", p$file), sep = "\n")
  } else if (identical(mode, "files")) {
    if (length(args) < 2L) stop("usage: shard-tests.R files <shard>", call. = FALSE)
    cat(tier3_shard_files(args[[2]]), sep = "\n")
  } else {
    stop("unknown mode '", mode, "'. Use ids, summary, plan or files <shard>.",
         call. = FALSE)
  }
})
