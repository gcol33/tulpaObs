# Standalone smoke for the progress reporter (no R fit).

devtools::load_all(quiet = TRUE)

# throttle = 0 so every call prints.
cat("--- rate only (no max_calls) ---\n")
r1 <- tulpaObs:::.tobs_progress_reporter("smoke", throttle = 0)
for (k in 1:4) { r1(value = -100 + k * 0.5); Sys.sleep(0.05) }

cat("\n--- with max_calls = 100 (ETA shrinks each call) ---\n")
r2 <- tulpaObs:::.tobs_progress_reporter("smoke-eta", throttle = 0,
                                          max_calls = 100L)
for (k in 1:4) { r2(value = -100 + k * 0.5); Sys.sleep(0.1) }

cat("\n--- minutes / hours formatting ---\n")
r3 <- tulpaObs:::.tobs_progress_reporter("smoke-long", throttle = 0,
                                          max_calls = 10000L)
for (k in 1:3) { r3(value = -50 + k); Sys.sleep(0.5) }
cat("DONE\n")
