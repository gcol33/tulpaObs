suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})
Sys.setenv(NOT_CRAN = "true")

# Just the S3 test — cheap, fast feedback.
testthat::test_file(
  "C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat/test-ms-abun.R",
  desc = "ms_abun S3 methods work",
  reporter = testthat::ProgressReporter$new()
)
