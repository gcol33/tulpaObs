suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(testthat)
})

# Run only the structural + gate tests (skip_on_cran ones are gated).
Sys.setenv(NOT_CRAN = "false")
testthat::test_file("tests/testthat/test-abun-re.R", reporter = "summary")
