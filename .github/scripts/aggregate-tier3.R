#!/usr/bin/env Rscript

# Collects the shard results of the full-recovery tier into one verdict.
#
# The point of sharding is that a shard which overruns the job cap costs only
# its own files. That is only true if the surviving shards' evidence is put on
# the record, and only useful if a partial tier is never mistaken for a complete
# one -- so this reports what every shard did, names the files no shard verified,
# and exits non-zero unless the whole plan reported green.
#
# A shard cancelled at the cap uploads nothing (its steps do not run), so the
# expected shard list is recomputed from the planner rather than read from the
# artifacts. What is missing is therefore always nameable.
#
# Usage: Rscript .github/scripts/aggregate-tier3.R [artifact-dir]
#   Expects <artifact-dir>/tier3-<shard>/{summary,timings}.csv, the layout
#   actions/download-artifact produces for a name pattern.

source(".github/scripts/shard-tests.R")

args <- commandArgs(trailingOnly = TRUE)
art_dir <- if (length(args)) args[[1]] else "tier3-artifacts"

plan <- tier3_plan()
ids <- tier3_shard_ids(plan)

read_csv_or_null <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) NULL else out
}

summaries <- list()
timings <- list()
for (id in ids) {
  base <- file.path(art_dir, paste0("tier3-", id))
  summaries[[id]] <- read_csv_or_null(file.path(base, "summary.csv"))
  timings[[id]] <- read_csv_or_null(file.path(base, "timings.csv"))
}

reported <- ids[!vapply(summaries[ids], is.null, logical(1))]
silent <- setdiff(ids, reported)

cat("full-recovery tier: shard report\n")
cat(sprintf("plan: %d test files across %d shards\n\n", nrow(plan), length(ids)))
cat(sprintf("%-28s %8s %8s %10s %8s %8s %9s\n",
            "shard", "files", "assert", "skipped", "failed", "errors", "min"))

tot_assert <- 0L; tot_skip <- 0L; tot_fail <- 0L; tot_err <- 0L
for (id in ids) {
  s <- summaries[[id]]
  if (is.null(s)) {
    cat(sprintf("%-28s %8s %8s %10s %8s %8s %9s   NO REPORT\n",
                id, nrow(plan[plan$shard == id, ]), "-", "-", "-", "-", "-"))
    next
  }
  tot_assert <- tot_assert + sum(s$assertions)
  tot_skip <- tot_skip + sum(s$skipped)
  tot_fail <- tot_fail + sum(s$failed)
  tot_err <- tot_err + sum(s$errors)
  cat(sprintf("%-28s %8d %8d %10d %8d %8d %9.1f\n",
              id, s$files_ran[1], s$assertions[1], s$skipped[1],
              s$failed[1], s$errors[1], s$seconds[1] / 60))
}
cat(sprintf("\n%-28s %8s %8d %10d %8d %8d\n", "TOTAL (reporting shards)", "",
            tot_assert, tot_skip, tot_fail, tot_err))

# Coverage is asserted against files that produced results, not against the
# plan: a shard is only credited with the files it demonstrably ran.
verified <- unlist(lapply(timings[reported], function(t) t$file), use.names = FALSE)
verified <- unique(as.character(verified))
unverified <- setdiff(plan$file, verified)

# One timing table for the whole tier: this is what refreshes
# .github/tier3-weights.csv, and the only way the placeholder weights ever
# become measured ones.
timing_parts <- lapply(reported, function(id) {
  t <- timings[[id]]
  if (is.null(t)) return(NULL)
  data.frame(shard = id, t[, c("file", "seconds", "assertions", "skipped",
                               "failed", "errors")],
             stringsAsFactors = FALSE)
})
timing_parts <- timing_parts[!vapply(timing_parts, is.null, logical(1))]
all_timings <- if (length(timing_parts)) {
  do.call(rbind, timing_parts)
} else {
  NULL
}

if (!is.null(all_timings) && nrow(all_timings)) {
  all_timings <- all_timings[order(-all_timings$seconds), , drop = FALSE]
  utils::write.csv(all_timings, "tier3-timings.csv", row.names = FALSE)
  cat("\nSlowest files this run (full table in tier3-timings.csv):\n")
  for (i in seq_len(min(20L, nrow(all_timings)))) {
    cat(sprintf("  %9.1fs  %-46s %s\n", all_timings$seconds[i],
                all_timings$file[i], all_timings$shard[i]))
  }
}

problems <- character(0)
if (length(silent)) {
  problems <- c(problems, sprintf(
    "%d shard(s) did not report: %s", length(silent), paste(silent, collapse = ", ")))
}
if (length(unverified)) {
  cat(sprintf("\n%d test file(s) unverified by this run:\n", length(unverified)))
  cat(paste0("  ", unverified, collapse = "\n"), "\n")
  problems <- c(problems, sprintf("%d test file(s) produced no results",
                                  length(unverified)))
}
if (tot_fail > 0 || tot_err > 0) {
  problems <- c(problems, sprintf("%d failure(s) and %d error(s) across reporting shards",
                                  tot_fail, tot_err))
}

if (length(problems)) {
  cat("\nTier NOT green:\n")
  cat(paste0("  - ", problems, collapse = "\n"), "\n")
  cat("\nThe reporting shards' results above still stand: their recovery and\n")
  cat("coverage blocks ran and are on the record for this commit.\n")
  quit(status = 1L)
}

cat("\nTier green: every shard reported and every test file ran.\n")
