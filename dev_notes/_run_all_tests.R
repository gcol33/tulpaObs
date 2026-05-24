setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
res <- as.data.frame(devtools::test(reporter = testthat::SummaryReporter$new()))
cat(sprintf("\n==> failed=%d  error=%d  skipped=%d  passed=%d\n",
            sum(res$failed), sum(res$error), sum(res$skipped),
            sum(res$nb) - sum(res$failed) - sum(res$error)))
bad <- res[res$failed > 0 | res$error > 0, c("file", "test", "failed", "error")]
if (nrow(bad)) { cat("\nNON-GREEN:\n"); print(bad) } else cat("ALL GREEN\n")
