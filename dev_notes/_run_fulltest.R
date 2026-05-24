# Run every test file under load_all in one session, flushing per-file results
# to dev_notes/_fulltest_results.txt BEFORE and AFTER each file so a segfault
# (Windows R, complex C++) pinpoints the offending file. Not inline -e (that
# segfaults on Windows per CLAUDE.md).
#
#   "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/_run_fulltest.R

Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all("."))
library(testthat)

out <- "dev_notes/_fulltest_results.txt"
cat("", file = out)
log_line <- function(s) { cat(s, "\n", file = out, append = TRUE); flush.console() }

files <- list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
tot_fail <- 0L
for (f in files) {
  log_line(sprintf("START %s", basename(f)))
  r <- tryCatch(
    as.data.frame(test_file(f, reporter = "silent")),
    error = function(e) { log_line(sprintf("  ERROR %s", conditionMessage(e))); NULL }
  )
  if (!is.null(r)) {
    fa <- sum(r$failed); wa <- sum(r$warning); sk <- sum(r$skipped)
    tot_fail <- tot_fail + fa
    log_line(sprintf("DONE  %s  fail=%d warn=%d skip=%d", basename(f), fa, wa, sk))
  }
}
log_line(sprintf("=== TOTAL FAILURES: %d ===", tot_fail))
