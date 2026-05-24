setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(TULPAOBS_FAST = "1")
suppressMessages(devtools::load_all(".", quiet = TRUE))
t0 <- Sys.time()
for (f in c("test-occu-laplace-se.R", "test-occu-prior.R",
            "test-issue8-visit-detection.R")) {
  res <- as.data.frame(testthat::test_file(
    file.path("tests/testthat", f), reporter = testthat::SilentReporter$new()))
  cat(sprintf("%-32s skipped=%d failed=%d error=%d passed=%d\n", f,
              sum(res$skipped), sum(res$failed), sum(res$error),
              sum(res$nb) - sum(res$failed) - sum(res$error) - sum(res$skipped)))
}
cat(sprintf("FAST-mode wall-clock for 3 files: %.2fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
