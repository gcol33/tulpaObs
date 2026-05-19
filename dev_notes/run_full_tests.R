suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa", quiet = TRUE, export_all = FALSE)
  load_all(".",        quiet = TRUE, export_all = FALSE)
  library(testthat)
})
res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nTotals: pass=%d  fail=%d  warn=%d  skip=%d\n",
            sum(df$passed), sum(df$failed), sum(df$warning), sum(df$skipped)))
if (any(df$failed > 0)) {
  cat("\nFailures:\n")
  print(df[df$failed > 0, c("file", "test", "failed")])
}
