suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})
Sys.setenv(NOT_CRAN = "true")
res <- testthat::test_file(
  "C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat/test-ms-abun.R",
  reporter = testthat::SilentReporter$new()
)
df <- as.data.frame(res)
cat(sprintf("PASS %d  FAIL %d  WARN %d  SKIP %d\n",
            sum(df$nb), sum(df$failed), sum(df$warning), sum(df$skipped)))
if (sum(df$failed) > 0) {
  cat("\n=== failures ===\n")
  print(subset(df, failed > 0))
}
