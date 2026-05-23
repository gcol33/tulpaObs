# Verify the formula-native migration: load_all + the migrated test files
# (structural / smoke tiers; the multi-seed recovery tiers run separately).
suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))
cat("load_all clean\n\n")

files <- c(
  "test-formula-terms.R",
  "test-spatial-occ.R",
  "test-spde-occ.R",
  "test-nuts-components.R",
  "test-nested-laplace-occu.R",
  "test-cover-hurdle-lognormal.R",
  "test-cover-hurdle-multi-block.R",
  "test-cover-hurdle-nested-joint.R"
)
for (f in files) {
  cat("==== ", f, " ====\n", sep = "")
  testthat::test_file(file.path("tests/testthat", f), reporter = "summary")
  cat("\n")
}
