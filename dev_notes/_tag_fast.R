# One-time migration: pair skip_if_fast() with every skip_on_cran() in the test
# suite, so TULPAOBS_FAST=1 skips the expensive fitting/sampling blocks (the
# skip_on_cran set) while leaving the structural / dispatch / unit tests live.
# Idempotent: re-running does not duplicate the insert. Indentation is copied
# from each skip_on_cran() line (handles 2- and 4-space files alike).
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
files <- list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
total <- 0L
changed <- character(0)
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  hit <- grepl("^\\s*skip_on_cran\\(\\)\\s*$", lines)
  if (!any(hit)) next
  out <- character(0)
  n_ins <- 0L
  for (i in seq_along(lines)) {
    out <- c(out, lines[i])
    if (hit[i]) {
      nxt <- if (i < length(lines)) lines[i + 1L] else ""
      if (!grepl("^\\s*skip_if_fast\\(\\)\\s*$", nxt)) {
        indent <- sub("skip_on_cran.*$", "", lines[i])
        out <- c(out, paste0(indent, "skip_if_fast()"))
        n_ins <- n_ins + 1L
      }
    }
  }
  if (n_ins > 0L) {
    writeLines(out, f)
    total <- total + n_ins
    changed <- c(changed, sprintf("%s (+%d)", basename(f), n_ins))
  }
}
cat("Inserted", total, "skip_if_fast() calls across", length(changed), "files:\n")
cat(paste0("  ", changed, collapse = "\n"), "\n")
