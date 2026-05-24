# Build + targeted test runner for the engine/approx -> method= refactor.
suppressMessages({
  library(devtools)
})

cat("== document() ==\n")
devtools::document(quiet = TRUE)

cat("== load_all() ==\n")
devtools::load_all(quiet = TRUE)

# Smoke: the public method surface resolves correctly.
cat("== method resolution ==\n")
for (m in c("auto", "laplace", "laplace_sla", "laplace_gibbs",
            "laplace_mi", "nested_laplace", "nested_laplace_sla", "nuts")) {
  r <- tulpaObs:::.tobs_resolve_method(m, occu())
  cat(sprintf("  %-20s -> engine=%-15s approx=%-18s correction=%s\n",
              m, r$engine, r$approx, r$correction))
}

cat("\n== targeted tests ==\n")
files <- c(
  "test-method-gibbs-recovery.R",
  "test-occu-prior.R",
  "test-re-laplace-recovery.R",
  "test-cover-hurdle-lognormal.R",
  "test-cover-hurdle-multi-block.R",
  "test-simplified-laplace.R"
)
for (f in files) {
  cat("---", f, "---\n")
  res <- as.data.frame(testthat::test_file(
    file.path("tests/testthat", f), reporter = "summary"))
  fails <- sum(res$failed, na.rm = TRUE)
  warns <- sum(res$warning, na.rm = TRUE)
  cat(sprintf("  -> failed=%d warning=%d\n", fails, warns))
}
cat("\nDONE\n")
