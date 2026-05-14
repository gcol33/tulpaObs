# Smoke check the tulpaObs joint nested-Laplace dispatch.
suppressMessages({
  setwd("C:/Users/Gilles Colling/Documents/dev/tulpa")
  Rcpp::compileAttributes()
  devtools::document(quiet = TRUE)
  devtools::install(quiet = TRUE)
  setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
  devtools::load_all(quiet = TRUE)
})

library(testthat)
res <- test_file("tests/testthat/test-cover-hurdle-nested-joint.R",
                 reporter = SummaryReporter$new())
print(res)
