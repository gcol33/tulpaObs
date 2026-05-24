setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(TULPAOBS_FAST = "1")
suppressMessages(devtools::load_all(".", quiet = TRUE))
t0 <- Sys.time()
files <- list.files("tests/testthat", pattern = "^test-.*\\.R$")
tot <- c(skipped = 0L, failed = 0L, error = 0L, passed = 0L)
for (f in files) {
  res <- tryCatch(as.data.frame(testthat::test_file(
    file.path("tests/testthat", f), reporter = testthat::SilentReporter$new())),
    error = function(e) { cat("ERROR loading", f, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  sk <- sum(res$skipped); fl <- sum(res$failed); er <- sum(res$error)
  ps <- sum(res$nb) - fl - er - sk
  tot <- tot + c(sk, fl, er, ps)
  if (fl > 0 || er > 0) cat(sprintf("  %-40s failed=%d error=%d\n", f, fl, er))
}
cat(sprintf("\nFULL SUITE in FAST mode: skipped=%d failed=%d error=%d passed=%d  (%.1fs)\n",
            tot["skipped"], tot["failed"], tot["error"], tot["passed"],
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
