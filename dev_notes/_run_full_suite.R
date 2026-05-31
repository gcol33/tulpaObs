Sys.setenv(NOT_CRAN = "true")
Sys.setenv(TESTTHAT_PARALLEL = "false")
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE))
res <- testthat::test_dir("C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat",
                          reporter = testthat::SummaryReporter$new(),
                          stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nTOTAL  PASS=%d FAIL=%d WARN=%d SKIP=%d\n",
            sum(df$passed), sum(df$failed), sum(df$warning), sum(df$skipped)))
fails <- df[df$failed > 0, c("file", "test", "failed")]
if (nrow(fails)) { cat("\nFAILURES:\n"); print(fails, row.names = FALSE) }
