Sys.setenv(NOT_CRAN = "true")
pkgbuild::clean_dll("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs",
                                    quiet = TRUE, recompile = TRUE))
library(testthat)
files <- c("test-formula-terms.R",
           "test-issue8-visit-detection.R",
           "test-re-bar-recovery.R")
base <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat"
for (f in files) {
  cat("\n##### ", f, " #####\n", sep = "")
  testthat::test_file(file.path(base, f), reporter = "summary")
}
cat("\n=== all issue tests done ===\n")
