setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all("."))
res <- testthat::test_dir("tests/testthat",
                          reporter = testthat::SummaryReporter$new(),
                          stop_on_failure = FALSE,
                          stop_on_warning = FALSE)
df <- as.data.frame(res)
cat("\n--- TEST RESULTS ---\n")
cat(sprintf("Tests: %d  Failed: %d  Errors: %d  Warnings: %d  Skipped: %d\n",
            nrow(df), sum(df$failed), sum(df$error), sum(df$warning),
            sum(df$skipped)))
fails <- df[df$failed > 0L | df$error, c("file", "test", "failed", "error")]
if (nrow(fails) > 0L) {
  cat("\n--- FAILURES / ERRORS ---\n")
  print(fails, row.names = FALSE)
}
