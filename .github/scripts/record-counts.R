#!/usr/bin/env Rscript

# Write measured per-file counts into tests/expected-counts.csv.
#
# The manifest is what lets a run be verified to have RUN rather than read as
# green, but until now nothing could write it: run-tests.R only ever reads it,
# so every row was placed by hand. That is affordable for the smoke tier, which
# a single ten-minute sweep measures end to end, and it is not affordable for
# the full tier, which is roughly a day of compute spread over eleven shards.
# The full side therefore sat at a few dozen rows while the smoke side was
# complete -- not because the full rows were unwanted, but because there was no
# path from a completed run to them.
#
# This is that path. It takes the per-file table a run already produces and
# merges it into the manifest, so a tier that has finished once can record what
# it measured instead of having it transcribed.
#
# Usage:
#   Rscript .github/scripts/record-counts.R <timings.csv> --tier=full [options]
#
#   <timings.csv>   either <TULPAOBS_TEST_OUT>/timings.csv from a single run,
#                   or tier3-timings.csv from aggregate-tier3.R (which is the
#                   same table with a shard column). Several may be given and
#                   are concatenated.
#   --tier=         "smoke" or "full". Required: the manifest carries one row
#                   per (file, tier) and a run only ever measures its own.
#   --manifest=     default tests/expected-counts.csv
#   --dry-run       report what would change and write nothing.
#
# What it refuses to do, and why each one matters:
#
#   * Record from a run carrying any failure or error. The recorded count is a
#     floor the guard holds later runs to, so taking it from a red run pins
#     whatever the breakage happened to produce -- the one way this file can
#     make things worse rather than better. A red run is reported and nothing
#     is written, not even its green files: a failure early in a file truncates
#     the assertions after it, and that file's neighbours in the same shard are
#     not obviously unaffected either.
#   * Record a file that is not in tests/testthat. A manifest row naming a file
#     that does not exist fails the run by design, so writing one would break
#     the next run rather than the next rename.
#   * Record a file that neither asserted nor skipped anything. A 0/0 row
#     cannot fail the comparison -- got is never below want -- so it reads as
#     coverage while checking nothing.
#
# Merging is per (file, tier): a row already present is updated, a row for a
# file this run did not cover is left alone, and every other tier is untouched.
# A shard therefore records its own files without disturbing the rest, which is
# what makes the full tier fillable one completed shard at a time rather than
# only by a whole green tier.

args <- commandArgs(trailingOnly = TRUE)

opt_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[length(hit)]])
}

tier <- opt_value("--tier")
manifest_path <- opt_value("--manifest", "tests/expected-counts.csv")
dry_run <- "--dry-run" %in% args
inputs <- args[!grepl("^--", args)]

die <- function(...) {
  cat("record-counts: ", ..., "\n", sep = "")
  quit(status = 1L)
}

if (!length(inputs)) {
  die("give at least one timings CSV. ",
      "Usage: record-counts.R <timings.csv> --tier=full")
}
if (is.null(tier) || !tier %in% c("smoke", "full")) {
  die("--tier= must be \"smoke\" or \"full\"; a run only measures its own tier.")
}
missing_in <- inputs[!file.exists(inputs)]
if (length(missing_in)) {
  die("no such file: ", paste(missing_in, collapse = ", "))
}
if (!file.exists(manifest_path)) die("no manifest at ", manifest_path)

test_dir <- "tests/testthat"
if (!dir.exists(test_dir)) {
  die("run this from the package root: no ", test_dir)
}

need <- c("file", "assertions", "skipped", "failed", "errors", "seconds")
parts <- lapply(inputs, function(p) {
  d <- utils::read.csv(p, stringsAsFactors = FALSE)
  gap <- setdiff(need, names(d))
  if (length(gap)) {
    die(p, " is missing column(s): ", paste(gap, collapse = ", "),
        ". Expected the per-file table run-tests.R writes.")
  }
  d[, need, drop = FALSE]
})
got <- do.call(rbind, parts)
if (!nrow(got)) die("the input carries no rows.")

# A file split across inputs (one row per shard, or a file rerun) is summed for
# counts and takes the longest observed time: the guard compares a whole file's
# assertions against the row, so a partial row would understate it.
got <- do.call(rbind, lapply(split(got, got$file), function(d) {
  data.frame(file = d$file[1],
             assertions = sum(as.integer(d$assertions)),
             skipped = sum(as.integer(d$skipped)),
             failed = sum(as.integer(d$failed)),
             errors = sum(as.integer(d$errors)),
             seconds = max(as.numeric(d$seconds)),
             stringsAsFactors = FALSE)
}))
rownames(got) <- NULL

red <- got[got$failed > 0L | got$errors > 0L, , drop = FALSE]
if (nrow(red)) {
  cat("record-counts: refusing to record from a run that is not green.\n\n")
  for (i in seq_len(nrow(red))) {
    cat(sprintf("  %-46s failed %d  errors %d\n",
                red$file[i], red$failed[i], red$errors[i]))
  }
  cat("\nThe recorded count becomes the floor later runs are held to, so a red\n",
      "run would pin whatever the breakage produced. Fix the failures, rerun,\n",
      "and record that.\n", sep = "")
  quit(status = 1L)
}

on_disk <- basename(list.files(test_dir, pattern = "^test-.*[.][Rr]$"))
absent <- setdiff(got$file, on_disk)
if (length(absent)) {
  cat("record-counts: skipping ", length(absent),
      " file(s) not in ", test_dir, " (a row naming a missing file fails the run):\n",
      sep = "")
  for (f in absent) cat("  ", f, "\n", sep = "")
  got <- got[!got$file %in% absent, , drop = FALSE]
}

idle <- got[got$assertions == 0L & got$skipped == 0L, , drop = FALSE]
if (nrow(idle)) {
  cat("record-counts: skipping ", nrow(idle),
      " file(s) that neither asserted nor skipped (a 0/0 row cannot fail):\n",
      sep = "")
  for (f in idle$file) cat("  ", f, "\n", sep = "")
  got <- got[!(got$assertions == 0L & got$skipped == 0L), , drop = FALSE]
}

if (!nrow(got)) die("nothing left to record after those exclusions.")

man <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
man_need <- c("file", "tier", "assertions", "skipped", "seconds")
gap <- setdiff(man_need, names(man))
if (length(gap)) {
  die(manifest_path, " is missing column(s): ", paste(gap, collapse = ", "))
}

new_rows <- data.frame(file = got$file, tier = tier,
                       assertions = as.integer(got$assertions),
                       skipped = as.integer(got$skipped),
                       seconds = round(as.numeric(got$seconds), 1),
                       stringsAsFactors = FALSE)

key <- function(d) paste(d$file, d$tier, sep = "\r")
existing <- man[key(man) %in% key(new_rows), , drop = FALSE]
added <- new_rows[!key(new_rows) %in% key(man), , drop = FALSE]

changed <- 0L; same <- 0L
if (nrow(existing)) {
  m <- match(key(existing), key(new_rows))
  for (i in seq_len(nrow(existing))) {
    n <- new_rows[m[i], ]
    if (identical(as.integer(existing$assertions[i]), n$assertions) &&
        identical(as.integer(existing$skipped[i]), n$skipped)) {
      same <- same + 1L
    } else {
      changed <- changed + 1L
      cat(sprintf("  update %-46s assertions %d -> %d   skipped %d -> %d\n",
                  existing$file[i], as.integer(existing$assertions[i]),
                  n$assertions, as.integer(existing$skipped[i]), n$skipped))
    }
  }
}

out <- rbind(man[!key(man) %in% key(new_rows), man_need, drop = FALSE], new_rows)
out <- out[order(out$file, out$tier), , drop = FALSE]
rownames(out) <- NULL

cat(sprintf("\n%s tier: %d row(s) added, %d updated, %d already exact\n",
            tier, nrow(added), changed, same))
cat(sprintf("manifest: %d rows -> %d rows (%d smoke, %d full)\n",
            nrow(man), nrow(out), sum(out$tier == "smoke"),
            sum(out$tier == "full")))

uncovered <- setdiff(on_disk, out$file[out$tier == tier])
if (length(uncovered)) {
  cat(sprintf("%d test file(s) still have no %s row.\n",
              length(uncovered), tier))
}

if (dry_run) {
  cat("\n--dry-run: nothing written.\n")
  quit(status = 0L)
}

utils::write.csv(out, manifest_path, row.names = FALSE)
cat("\nwrote ", manifest_path, "\n", sep = "")
