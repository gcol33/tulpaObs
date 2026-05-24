setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
cat("tulpa ABI:", tryCatch(as.character(utils::packageVersion("tulpa")),
                           error = function(e) "??"), "\n")
for (f in c("test-nested-laplace-families.R", "test-nested-laplace-occu.R",
            "test-nested-laplace-prediction.R")) {
  res <- as.data.frame(testthat::test_file(
    file.path("tests/testthat", f), reporter = testthat::SilentReporter$new()))
  cat(sprintf("\n%-36s failed=%d error=%d passed=%d\n", f,
              sum(res$failed), sum(res$error),
              sum(res$nb) - sum(res$failed) - sum(res$error)))
  for (i in seq_len(nrow(res))) {
    if (res$failed[i] > 0 || res$error[i] > 0) {
      cat("   FAIL:", res$test[i], "\n")
      for (r in res$result[[i]]) {
        if (inherits(r, c("expectation_failure", "expectation_error"))) {
          msg <- gsub("\n", " | ", conditionMessage(r))
          if (nchar(msg) > 280) msg <- paste0(substr(msg, 1, 280), "...")
          cat("        ->", msg, "\n")
        }
      }
    }
  }
}
