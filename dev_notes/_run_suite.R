# Scratch: run the tulpaObs suite file-by-file with NOT_CRAN=true, appending a
# one-line result per file to _suite_results.tsv so a SIGKILL on one heavy file
# does not lose results for the rest (see global memory on full-suite SIGKILL).
# Not committed.  Usage: Rscript dev_notes/_run_suite.R [grep-substring-for-files]
Sys.setenv(NOT_CRAN = "true")
pkg  <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
suppressWarnings(suppressMessages(devtools::load_all(pkg, quiet = TRUE)))
cat("LOAD_OK tulpa=", as.character(packageVersion("tulpa")), "\n", sep = "")

tdir  <- file.path(pkg, "tests", "testthat")
files <- list.files(tdir, pattern = "^test-.*\\.R$", full.names = FALSE)
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[[1]])) files <- grep(args[[1]], files, value = TRUE)
# Run known-heavy recovery/coverage files LAST so a possible full-suite SIGKILL
# (global memory) costs the fewest already-computed results.
heavy <- c("test-cover-hurdle-nested-joint-recovery.R", "test-re-laplace-recovery.R",
           "test-re-bar-recovery.R", "test-method-gibbs-recovery.R",
           "test-cover-hurdle-adaptive-grid.R", "test-occu-laplace-se.R",
           "test-spde-occ.R", "test-abun.R", "test-sla-cover-joint.R")
files <- c(setdiff(files, heavy), intersect(heavy, files))

out <- file.path(pkg, "dev_notes", "_suite_results.tsv")
if (length(args) < 2L || args[[2]] != "append") {
  cat("file\tpass\tfail\twarn\tskip\tseconds\n", file = out)
}
library(testthat)
for (f in files) {
  t0 <- Sys.time()
  r  <- tryCatch(
    as.data.frame(test_file(file.path(tdir, f), reporter = SilentReporter$new())),
    error = function(e) {
      cat(sprintf("  [ERROR loading %s]: %s\n", f, conditionMessage(e)))
      NULL
    })
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (is.null(r)) {
    line <- sprintf("%s\tNA\tERROR\tNA\tNA\t%s\n", f, secs)
  } else {
    line <- sprintf("%s\t%d\t%d\t%d\t%d\t%s\n",
                    f, sum(r$passed), sum(r$failed), sum(r$warning),
                    sum(r$skipped), secs)
    if (sum(r$failed) > 0L) {
      bad <- r[r$failed > 0L, c("test", "failed")]
      cat(sprintf("  FAIL in %s:\n", f)); print(bad)
    }
  }
  cat(line)
  cat(line, file = out, append = TRUE)
}
cat("SUITE_DONE\n")
