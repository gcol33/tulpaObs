setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", recompile = TRUE, quiet = TRUE))
res <- testthat::test_file("tests/testthat/test-re-laplace-recovery.R",
                           reporter = testthat::SummaryReporter$new())
df <- as.data.frame(res)
cat(sprintf("\nfailed=%d  error=%d  skipped=%d  passed=%d\n",
            sum(df$failed), sum(df$error), sum(df$skipped),
            sum(df$nb) - sum(df$failed) - sum(df$error)))
