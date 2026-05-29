devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
fs <- list.files("C:/Users/Gilles Colling/Documents/dev/tulpaObs/tests/testthat",
                  pattern = "^test-.*\\.R$", full.names = TRUE)
totals <- list(pass = 0L, fail = 0L, warn = 0L, skip = 0L)
for (f in fs) {
  cat("== ", basename(f), " ==\n", sep = "")
  res <- tryCatch(testthat::test_file(f, reporter = testthat::SilentReporter$new()),
                  error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  df <- as.data.frame(res)
  totals$pass <- totals$pass + sum(df$nb)   # NB. nb is n passes? actually check
  totals$fail <- totals$fail + sum(df$failed)
  totals$warn <- totals$warn + sum(df$warning)
  totals$skip <- totals$skip + sum(df$skipped)
  cat(sprintf("   PASS %d  FAIL %d  WARN %d  SKIP %d\n",
              sum(df$nb), sum(df$failed), sum(df$warning), sum(df$skipped)))
}
cat("\n==== TOTALS ====\n")
cat(sprintf("PASS %d  FAIL %d  WARN %d  SKIP %d\n",
            totals$pass, totals$fail, totals$warn, totals$skip))
