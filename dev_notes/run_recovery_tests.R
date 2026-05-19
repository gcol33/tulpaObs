suppressPackageStartupMessages({
  library(devtools)
  library(testthat)
})
load_all(".", quiet = TRUE)
test_file("tests/testthat/test-cover-hurdle-nested-joint-recovery.R",
          reporter = "summary")
