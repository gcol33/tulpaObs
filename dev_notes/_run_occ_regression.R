suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})
files <- c(
  "tests/testthat/test-occ.R",
  "tests/testthat/test-occu-prior.R",
  "tests/testthat/test-occu-laplace-se.R",
  "tests/testthat/test-spatial-occ.R",
  "tests/testthat/test-spde-occ.R",
  "tests/testthat/test-nested-laplace-occu.R",
  "tests/testthat/test-methods.R",
  "tests/testthat/test-new-features.R",
  "tests/testthat/test-simplified-laplace.R"
)
for (f in files) {
  cat("\n==== ", f, " ====\n", sep = "")
  testthat::test_file(f, reporter = testthat::SummaryReporter$new())
}
