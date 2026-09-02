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
#   --update        also rewrite rows that already exist (see below). Without
#                   it, existing rows are compared and reported but never
#                   changed.
#   --dry-run       report what would change and write nothing.
#
# ADD-ONLY BY DEFAULT, and that is the whole safety argument.
#
# The two recorded numbers are not estimates of a true value, they are BOUNDS
# that run-tests.R holds a later run to, in opposite directions:
#
#     assertions   a FLOOR   (fails when a run asserts FEWER)
#     skipped      a CEILING (fails when a run skips MORE)
#
# So a row is only sound if it holds in EVERY environment the suite is expected
# to pass in, and one run measures exactly one environment. A richer box -- one
# with more Suggests installed -- asserts more and skips less, so recording
# from it raises the floor and lowers the ceiling together, and every leaner
# box then fails on a suite that is perfectly healthy. Measured instance:
# test-cover-hurdle-beta.R reports 25 assertions where betareg is installed and
# 20 where it is not, and the manifest correctly holds 20. Recording a CI run
# over that row would write 25 and break every box without betareg.
#
# Lowering is not the safe direction either, and this is the sharper of the
# two. A block that stops executing makes a run assert fewer -- which is the
# defect this manifest exists to catch, and the guard does catch it. A recorder
# that then wrote the smaller number would erase the finding and hand back a
# green suite. Automatically taking the minimum would be a tool for hiding
# exactly what the manifest is for.
#
# Since neither direction can be resolved from counts alone -- "fewer
# assertions" is indistinguishable between a leaner environment and a block
# that died -- the default changes no row that already exists. It ADDS rows for
# files that have none, which is the whole of the full tier's backlog and needs
# no judgement, and it REPORTS every difference on the rest for a person to
# read. Pass --update when you know why a count moved and are recording it in
# the same commit as the change that moved it.
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
do_update <- "--update" %in% args
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
added <- new_rows[!key(new_rows) %in% key(man), , drop = FALSE]
existing <- man[key(man) %in% key(new_rows), , drop = FALSE]

# Every difference on a row that already exists is reported whether or not it
# is going to be written, and labelled with which bound it moves and what that
# would do to a run. "assertions below" is the guard's own failing case, so it
# is a finding to read rather than a number to record.
same <- 0L
drift <- list()
if (nrow(existing)) {
  m <- match(key(existing), key(new_rows))
  for (i in seq_len(nrow(existing))) {
    n <- new_rows[m[i], ]
    was_a <- as.integer(existing$assertions[i])
    was_s <- as.integer(existing$skipped[i])
    if (identical(was_a, n$assertions) && identical(was_s, n$skipped)) {
      same <- same + 1L
      next
    }
    note <- if (n$assertions < was_a || n$skipped > was_s) {
      "run is BELOW the recorded bound -- the guard fails on this; read it"
    } else {
      "manifest is behind, or this box has more Suggests than the leanest one"
    }
    drift[[length(drift) + 1L]] <- data.frame(
      file = existing$file[i], was_a = was_a, now_a = n$assertions,
      was_s = was_s, now_s = n$skipped, note = note, stringsAsFactors = FALSE)
  }
}
drift <- if (length(drift)) do.call(rbind, drift) else NULL

if (!is.null(drift)) {
  cat(sprintf("\n%d existing row(s) differ from this run%s:\n", nrow(drift),
              if (do_update) " (--update: rewriting them)" else
                " (not changed; pass --update to rewrite)"))
  for (i in seq_len(nrow(drift))) {
    cat(sprintf("  %-46s assertions %d -> %d   skipped %d -> %d\n    %s\n",
                drift$file[i], drift$was_a[i], drift$now_a[i],
                drift$was_s[i], drift$now_s[i], drift$note[i]))
  }
}

if (do_update) {
  keep <- man[!key(man) %in% key(new_rows), man_need, drop = FALSE]
  out <- rbind(keep, new_rows)
} else {
  # Existing rows survive verbatim; only files with no row for this tier are
  # added.
  out <- rbind(man[, man_need, drop = FALSE], added)
}
out <- out[order(out$file, out$tier), , drop = FALSE]
rownames(out) <- NULL

cat(sprintf("\n%s tier: %d row(s) added, %d rewritten, %d already exact\n",
            tier, nrow(added),
            if (do_update && !is.null(drift)) nrow(drift) else 0L, same))
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
