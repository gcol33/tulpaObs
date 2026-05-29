suppressMessages({ devtools::load_all(".", quiet = TRUE); library(testthat) })
Sys.setenv(NOT_CRAN = "false")
testthat::test_file("tests/testthat/test-abun.R", reporter = "summary")
