suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})
res <- testthat::test_file("tests/testthat/test-visit-data.R",
                           reporter = testthat::SummaryReporter$new())
invisible(res)
