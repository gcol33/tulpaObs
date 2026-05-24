setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
for (f in c("test-occu-laplace-se.R", "test-occu-prior.R",
            "test-issue8-visit-detection.R")) {
  res <- as.data.frame(testthat::test_file(
    file.path("tests/testthat", f), reporter = testthat::SilentReporter$new()))
  cat(sprintf("\n%-32s failed=%d error=%d passed=%d\n", f,
              sum(res$failed), sum(res$error),
              sum(res$nb) - sum(res$failed) - sum(res$error)))
  for (i in seq_len(nrow(res))) {
    if (res$failed[i] > 0 || res$error[i] > 0) {
      cat("   FAIL:", res$test[i], "\n")
      rr <- res$result[[i]]
      for (r in rr) {
        if (inherits(r, c("expectation_failure", "expectation_error"))) {
          msg <- conditionMessage(r)
          msg <- gsub("\n", " | ", msg)
          if (nchar(msg) > 300) msg <- paste0(substr(msg, 1, 300), "...")
          cat("        ->", msg, "\n")
        }
      }
    }
  }
}
