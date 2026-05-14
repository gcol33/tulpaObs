# Smoke check the tulpaObs joint nested-Laplace dispatch end-to-end. We
# load tulpa from source (its compiled .dll is already current at the
# install location, so pkgload reuses it) so the new R driver is visible.
suppressMessages({
  pkgload::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  pkgload::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

library(testthat)
res <- test_file("tests/testthat/test-cover-hurdle-nested-joint.R",
                 reporter = SummaryReporter$new())
print(res)
