devtools::load_all(".", quiet = TRUE)
Sys.setenv(NOT_CRAN = "true")
for (f in c("tests/testthat/test-formula-terms.R",
            "tests/testthat/test-nuts-components.R")) {
  cat("\n==== ", f, " ====\n")
  res <- testthat::test_file(f, reporter = testthat::SummaryReporter$new())
  df <- as.data.frame(res)
  cat(sprintf("PASS=%d FAIL=%d WARN=%d SKIP=%d\n",
              sum(df$passed), sum(df$failed), sum(df$warning), sum(df$skipped)))
}
