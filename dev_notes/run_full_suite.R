Sys.setenv(NOT_CRAN = "true")
res <- as.data.frame(devtools::test(
    "C:/Users/Gilles Colling/Documents/dev/tulpaObs",
    reporter = "silent"
))
cat("FILES:", length(unique(res$file)), "\n")
cat("TESTS:", nrow(res), "\n")
cat("PASS :", sum(res$nb), "\n")
cat("FAIL :", sum(res$failed), "\n")
cat("WARN :", sum(res$warning), "\n")
cat("SKIP :", sum(res$skipped), "\n")
fails <- subset(res, failed > 0)
if (nrow(fails) > 0L) {
    cat("\nFAILING TESTS:\n")
    for (i in seq_len(nrow(fails))) {
        cat(sprintf("  %s :: %s\n",
                    basename(fails$file[i]), fails$test[i]))
    }
} else {
    cat("\nALL TESTS PASSED.\n")
}
