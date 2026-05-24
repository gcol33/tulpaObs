# Runner: verify the spatial() umbrella parses + dispatches correctly.
devtools::load_all(".", quiet = TRUE)
testthat::test_file("tests/testthat/test-formula-terms.R")
