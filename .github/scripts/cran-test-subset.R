# Writes tests/cran-files.txt: the test files CRAN runs (tests/testthat.R reads
# it when NOT_CRAN is unset). Each file is timed under CRAN conditions against
# the INSTALLED package, and files are taken cheapest first until the budget is
# spent, so the check's test step stays inside CRAN's time limit.
#
#   Rscript .github/scripts/cran-test-subset.R [--budget=40] [--timings=file.csv]
#
# --timings reuses a CSV with columns file,elapsed instead of measuring.
# The budget is local seconds; win-builder r-devel ran this suite 1.4x slower.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
budget <- as.numeric(opt("budget", "40"))
timings <- opt("timings", NA_character_)

if (is.na(timings)) {
  Sys.unsetenv(c("NOT_CRAN", "TULPAOBS_FAST", "TULPAOBS_CRAN_ALL"))
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  files <- list.files("tests/testthat", "^test-.*\\.R$", full.names = TRUE)
  elapsed <- vapply(files, function(f) {
    t0 <- proc.time()[["elapsed"]]
    try(testthat::test_file(f, package = "tulpaObs",
                            load_package = "installed", reporter = "silent"),
        silent = TRUE)
    proc.time()[["elapsed"]] - t0
  }, numeric(1))
  tab <- data.frame(file = basename(files), elapsed = unname(elapsed))
} else {
  tab <- utils::read.csv(timings)
}

tab <- tab[order(tab$elapsed, tab$file), ]
keep <- sort(tab$file[cumsum(tab$elapsed) <= budget])
writeLines(c(
  "# Test files run on CRAN; written by .github/scripts/cran-test-subset.R.",
  sprintf("# Budget %g s local, %d of %d files, %.1f s measured.", budget,
          length(keep), nrow(tab), sum(tab$elapsed[tab$file %in% keep])),
  keep
), "tests/cran-files.txt")
cat(sprintf("%d of %d files, %.1f s\n", length(keep), nrow(tab),
            sum(tab$elapsed[tab$file %in% keep])))
