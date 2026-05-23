# Task 4 verification: builders now parse formulas + attach structured_terms.
# Confirm (a) the parser tests still pass, (b) the plain-formula occupancy
# path still builds identical designs and fits (spatial= arg path untouched).
devtools::load_all(".", quiet = TRUE)

cat("== formula-term parser tests ==\n")
testthat::test_file("tests/testthat/test-formula-terms.R")

cat("\n== single-season occupancy (plain formula) ==\n")
testthat::test_file("tests/testthat/test-occ.R")

cat("\n== spatial occupancy (spatial= arg path, builder unchanged consumer) ==\n")
testthat::test_file("tests/testthat/test-spatial-occ.R")
