Sys.setenv(NOT_CRAN = "false")  # match CRAN-style skip behavior
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa", quiet = TRUE, export_all = FALSE)
  load_all(".",         quiet = TRUE, export_all = FALSE)
  library(testthat)
})
# Run the full test suite via test() — this exercises the same path
# devtools::check uses but skips the cran-style ERRORs (no INLAocc deps,
# missing test files wrapper, etc) that aren't related to the SLA work.
res <- devtools::test(filter = NULL, reporter = "summary")
cat("\n\n--- summary ---\n")
print(res)
