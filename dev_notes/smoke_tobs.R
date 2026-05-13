devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs")

cat("=== family constructors ===\n")
print(occu())
cat("\n")
print(abun(K_max = 50))
cat("\n")
print(cover(positive = "beta"))

cat("\n=== validation: missing family ===\n")
tryCatch(
  tobs(~1, data.frame(x = 1)),
  error = function(e) cat("OK error:", conditionMessage(e), "\n")
)

cat("\n=== validation: planned family ===\n")
tryCatch(
  tobs(~1, data.frame(x = 1), family = cover()),
  error = function(e) cat("OK error:", conditionMessage(e), "\n")
)

cat("\n=== validation: bad family argument ===\n")
tryCatch(
  tobs(~1, data.frame(x = 1), family = "occu"),
  error = function(e) cat("OK error:", conditionMessage(e), "\n")
)

cat("\nDONE\n")
