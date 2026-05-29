suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})
Sys.setenv(NOT_CRAN = "true")
testthat::test_file(
  "C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat/test-ms-abun.R",
  reporter = testthat::ProgressReporter$new()
)
