# Capture the exact error text of the NUTS-component failures. If it is the
# "tulpa ABI mismatch" guard (a clean stop(), by design "a clear error instead
# of segfault"), the 5 errors are install-race artifacts, not a code crash.
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
res <- testthat::test_file("tests/testthat/test-nuts-components.R",
                           reporter = testthat::SilentReporter$new())
df <- as.data.frame(res)
for (i in seq_len(nrow(df))) if (df$error[i] > 0) {
  for (r in res[[i]]$results) if (inherits(r, "expectation_error")) {
    cat("[", df$test[i], "]\n   ", conditionMessage(r), "\n"); break
  }
}
