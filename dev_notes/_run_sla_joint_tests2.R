Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa", quiet = TRUE, export_all = FALSE)
  load_all(".",         quiet = TRUE, export_all = FALSE)
  library(testthat)
})
res <- test_file("tests/testthat/test-sla-cover-joint.R", reporter = "summary")
print(res)
