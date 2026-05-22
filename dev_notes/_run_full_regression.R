suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})
tr <- testthat::test_dir("tests/testthat", reporter = "summary")
invisible(tr)
