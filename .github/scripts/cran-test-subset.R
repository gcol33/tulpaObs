# Writes tests/cran-files.txt: the test files CRAN runs (tests/testthat.R reads
# it when NOT_CRAN is unset). Each file is run under CRAN conditions against the
# INSTALLED package, recording its wall time and the assertions it makes there.
# Files that assert nothing on CRAN are dropped; the rest are taken in order of
# assertions per second until the budget is spent, so the check's test step
# stays inside CRAN's time limit while exercising as much as it can.
#
#   Rscript .github/scripts/cran-test-subset.R [--budget=15] [--timings=file.csv]
#                                              [--out=file.csv]
#
# --timings reuses a CSV with columns file,elapsed,assertions,failed instead of
# measuring; --out saves the measured table in that layout.
# The budget is local seconds; win-builder r-devel ran this suite 1.4x slower.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
budget <- as.numeric(opt("budget", "15"))
timings <- opt("timings", NA_character_)
out <- opt("out", NA_character_)

if (is.na(timings)) {
  Sys.unsetenv(c("NOT_CRAN", "TULPAOBS_FAST", "TULPAOBS_CRAN_ALL"))
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  files <- list.files("tests/testthat", "^test-.*\\.R$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(
      as.data.frame(testthat::test_file(f, package = "tulpaObs",
                                        load_package = "installed",
                                        reporter = "silent")),
      error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(res)) {
      return(data.frame(file = basename(f), elapsed = elapsed,
                        assertions = 0L, failed = 1L))
    }
    data.frame(file = basename(f), elapsed = elapsed,
               assertions = sum(res$passed),
               failed = sum(res$failed) + sum(res$error))
  })
  tab <- do.call(rbind, rows)
  if (!is.na(out)) utils::write.csv(tab, out, row.names = FALSE)
} else {
  tab <- utils::read.csv(timings)
}

tab <- tab[tab$assertions > 0 & tab$failed == 0, ]
tab <- tab[order(-tab$assertions / tab$elapsed, tab$file), ]
take <- cumsum(tab$elapsed) <= budget
keep <- sort(tab$file[take])
writeLines(c(
  "# Test files run on CRAN; written by .github/scripts/cran-test-subset.R.",
  sprintf("# Budget %g s local, %d of %d asserting files, %.1f s, %d assertions measured.",
          budget, length(keep), nrow(tab), sum(tab$elapsed[take]),
          sum(tab$assertions[take])),
  keep
), "tests/cran-files.txt")
cat(sprintf("%d of %d asserting files, %.1f s, %d of %d assertions\n",
            length(keep), nrow(tab), sum(tab$elapsed[take]),
            sum(tab$assertions[take]), sum(tab$assertions)))
