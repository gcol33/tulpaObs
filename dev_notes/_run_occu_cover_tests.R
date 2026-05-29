suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})
files <- c(
  "tests/testthat/test-occu-cover.R",
  "tests/testthat/test-occu-cover-coupling.R",
  "tests/testthat/test-occu-cover-joint-coupled.R",
  "tests/testthat/test-occu-cover-spatial.R"
)
for (f in files) {
  cat("== ", basename(f), " ==\n", sep = "")
  res <- tryCatch(
    testthat::test_file(file.path("C:/Users/Gilles Colling/Documents/dev/tulpaObs", f),
                        reporter = testthat::SilentReporter$new()),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res)) next
  df <- as.data.frame(res)
  cat(sprintf("  PASS %d  FAIL %d  WARN %d  SKIP %d\n",
              sum(df$nb), sum(df$failed), sum(df$warning), sum(df$skipped)))
}
