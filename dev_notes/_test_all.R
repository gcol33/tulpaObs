# Full tulpaObs test suite: confirm the ms_abun wiring (dispatch switch,
# methods.R model_type branches, obs_families status flip) left the existing
# families intact.
res <- devtools::test("C:/Users/Gilles Colling/Documents/dev/tulpaObs",
                      stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n=== FULL SUITE: %d tests, %d failed, %d skipped, %d warnings ===\n",
            nrow(df), sum(df$failed), sum(df$skipped), sum(df$warning)))
if (sum(df$failed) > 0) print(df[df$failed > 0, c("file", "test", "failed")])
