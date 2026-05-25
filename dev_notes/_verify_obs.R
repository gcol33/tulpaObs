# Scratch verification harness for the tulpaObs Phase-6 migration. Not committed.
# Recompiles tulpaObs C++ against the freshly-installed current tulpa headers,
# then runs a filtered slice of the suite with NOT_CRAN=true so the
# recovery / CI-coverage assertions (skip_on_cran) actually execute.
#
# Usage: Rscript dev_notes/_verify_obs.R [test-file-substring]
Sys.setenv(NOT_CRAN = "true")
pkg <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
suppressWarnings(suppressMessages({
  Rcpp::compileAttributes(pkg)
  devtools::load_all(pkg, quiet = TRUE)
}))
cat("LOAD_OK  tulpa=", as.character(packageVersion("tulpa")),
    " tulpaObs=", as.character(packageVersion("tulpaObs")), "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
filt <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else NULL
res <- testthat::test_local(
  pkg,
  filter   = filt,
  reporter = testthat::SummaryReporter$new(),
  stop_on_failure = FALSE
)
df <- as.data.frame(res)
cat("\n==== SUMMARY ====\n")
cat("files:", length(unique(df$file)),
    " fail:", sum(df$failed),
    " warn:", sum(df$warning),
    " skip:", sum(df$skipped),
    " pass:", sum(df$passed), "\n")
if (sum(df$failed) > 0L) {
  bad <- df[df$failed > 0L, c("file", "context", "test", "failed")]
  cat("---- FAILURES ----\n"); print(bad)
}
cat("TESTS_DONE\n")
