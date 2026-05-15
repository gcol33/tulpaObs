# Install tulpaObs against the option-B tulpa + run the existing joint-beta
# recovery test. Then a sparse-n_pos recovery probe at n_pos ~ 46 (D7 cell B
# regime) to verify the phi-on-outer-grid fix.
suppressMessages({
  setwd("C:/Users/GillesC/Documents/dev/tulpaObs")
  cat("--- devtools::install() tulpaObs ---\n")
  devtools::install(".", quiet = TRUE, upgrade = FALSE, quick = TRUE,
                    args = "--no-multiarch")
  library(tulpaObs)
  library(tulpa)
})

Sys.setenv(NOT_CRAN = "true")
cat("--- run joint-beta dense recovery (n_pos median ~325, 10 seeds) ---\n")
testthat::test_dir("tests/testthat",
                    filter = "cover-hurdle-nested-joint-recovery",
                    reporter = "summary",
                    stop_on_failure = FALSE)
