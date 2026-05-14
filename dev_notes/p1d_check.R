# Phase 1d cross-check: within_between() helper + beta-positive joint engine.
suppressMessages({
  pkgload::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  pkgload::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

library(testthat)
files <- c(
  "C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat/test-within-between.R",
  "C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat/test-cover-hurdle-nested-joint.R"
)
total_pass <- 0L
total_fail <- 0L
for (f in files) {
  cat("\n## ", basename(f), "\n", sep = "")
  res <- test_file(f, reporter = SilentReporter$new())
  s <- as.data.frame(res)
  total_pass <- total_pass + sum(s$passed)
  total_fail <- total_fail + sum(s$failed)
  cat("  passed:", sum(s$passed), " failed:", sum(s$failed), "\n")
  if (sum(s$failed) > 0) {
    print(s[s$failed > 0, c("test", "failed", "warning", "error")])
  }
}
cat("\n== TOTAL ==  passed:", total_pass, " failed:", total_fail, "\n")
