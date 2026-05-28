# Standalone smoke for the progress reporter (no R fit).

devtools::load_all(quiet = TRUE)

# throttle = 0 so every call prints.
r <- tulpaObs:::.tobs_progress_reporter("smoke", throttle = 0)
for (k in 1:5) {
  r(value = -100 + k * 0.5, extra = sprintf("step %d", k))
  Sys.sleep(0.05)
}
cat("DONE\n")
