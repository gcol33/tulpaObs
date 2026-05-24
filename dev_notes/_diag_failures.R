setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
for (f in c("test-occu-laplace-se.R", "test-occu-prior.R",
            "test-issue8-visit-detection.R")) {
  cat("\n\n##########", f, "##########\n")
  res <- testthat::test_file(file.path("tests/testthat", f),
                             reporter = testthat::SilentReporter$new())
  df <- as.data.frame(res)
  for (i in seq_len(nrow(df))) {
    if (df$failed[i] > 0 || df$error[i] > 0) {
      cat("\n--", df$test[i], "--\n")
      rr <- res[[i]]$results
      for (r in rr) if (inherits(r, c("expectation_failure", "expectation_error")))
        cat(conditionMessage(r), "\n")
    }
  }
}
