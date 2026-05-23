# Runner: load_all + run the formula-term parser tests in isolation.
devtools::load_all(".", quiet = TRUE)
testthat::test_file("tests/testthat/test-formula-terms.R")
