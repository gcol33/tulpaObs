# Bisect helper: install tulpa from the detached bisect worktree (whatever
# commit it's currently checked out at), then run tulpaObs's no-RE Laplace
# failing files in the SAME process (so a concurrent reinstall can't race in
# between). The Laplace path never calls cpp_occu_fit, so the tulpaObs<->tulpa
# ABI check is not triggered -- we can test any tulpa commit against the
# already-compiled tulpaObs.dll.
WT <- "C:/Users/Gilles Colling/Documents/dev/_tulpa_bisect"
suppressMessages(devtools::install(WT, quick = TRUE, build = FALSE,
                                   upgrade = FALSE, quiet = TRUE))
h <- system.file("include/tulpa/model_data.h", package = "tulpa")
abi <- grep("TULPA_ABI_VERSION =", readLines(h), value = TRUE)
cat("== installed tulpa:", trimws(abi),
    "| desc:", as.character(utils::packageDescription("tulpa")$Version), "==\n")

setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
for (f in c("test-occu-laplace-se.R", "test-occu-prior.R",
            "test-issue8-visit-detection.R")) {
  res <- as.data.frame(testthat::test_file(
    file.path("tests/testthat", f), reporter = testthat::SilentReporter$new()))
  cat(sprintf("  %-32s failed=%d error=%d passed=%d\n", f,
              sum(res$failed), sum(res$error),
              sum(res$nb) - sum(res$failed) - sum(res$error)))
}
