# Scratch: surface the 3 warnings in test-abun.R. Not committed.
Sys.setenv(NOT_CRAN = "true")
pkg <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
suppressWarnings(suppressMessages(devtools::load_all(pkg, quiet = TRUE)))
library(testthat)
res <- as.data.frame(test_file(
  file.path(pkg, "tests/testthat/test-abun.R"),
  reporter = ProgressReporter$new(show_praise = FALSE)))
cat("\n==== rows with warnings ====\n")
w <- res[res$warning > 0L, c("test", "warning")]
print(w)
cat("ABUN_DONE\n")
