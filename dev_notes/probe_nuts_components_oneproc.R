# Confirm all 5 NUTS-component "errors" were ABI-mismatch artifacts: run the
# whole file in one process at matched ABI 24 (install tulpa, recompile+load
# tulpaObs, test_file -- no between-process window for a concurrent reinstall).
suppressMessages({
  try(remove.packages("tulpa"), silent = TRUE)
  devtools::install("C:/Users/Gilles Colling/Documents/dev/_tulpa_bisect",
                    quick = TRUE, build = FALSE, upgrade = FALSE, quiet = TRUE)
})
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", recompile = TRUE, quiet = TRUE))
cat("in-process tulpa ABI:",
    trimws(grep("ABI_VERSION =", readLines(system.file(
      "include/tulpa/model_data.h", package = "tulpa")), value = TRUE)), "\n")
res <- as.data.frame(testthat::test_file(
  "tests/testthat/test-nuts-components.R",
  reporter = testthat::SilentReporter$new()))
cat(sprintf("NUTS components: failed=%d error=%d skipped=%d passed=%d\n",
            sum(res$failed), sum(res$error), sum(res$skipped),
            sum(res$nb) - sum(res$failed) - sum(res$error)))
