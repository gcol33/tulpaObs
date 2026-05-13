source("C:/GillesC/Documents/dev/tulpaObs/R/obs_families.R")
source("C:/GillesC/Documents/dev/tulpaObs/R/tulpa_obs.R")

cat("=== family constructors ===\n")
print(occ())
cat("\n")
print(nmixture(K_max = 50))
cat("\n")
print(cover_hurdle(positive = "beta"))

cat("\n=== validation: missing family ===\n")
tryCatch(
  tulpa_obs(~1, data.frame(x = 1)),
  error = function(e) cat("OK error:", conditionMessage(e), "\n")
)

cat("\n=== validation: planned family ===\n")
tryCatch(
  tulpa_obs(~1, data.frame(x = 1), family = cover_hurdle()),
  error = function(e) cat("OK error:", conditionMessage(e), "\n")
)

cat("\n=== validation: bad family argument ===\n")
tryCatch(
  tulpa_obs(~1, data.frame(x = 1), family = "occ"),
  error = function(e) cat("OK error:", conditionMessage(e), "\n")
)

cat("\nDONE\n")
