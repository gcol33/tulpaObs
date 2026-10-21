library(testthat)
library(tulpaObs)

# On CRAN (NOT_CRAN unset) only the files listed in cran-files.txt run, to keep
# the check inside CRAN's time limit. TULPAOBS_CRAN_ALL=1 runs every file under
# the same skip_on_cran() gating (cran-like.yaml).
filter <- NULL
if (!identical(Sys.getenv("NOT_CRAN"), "true") &&
    !identical(Sys.getenv("TULPAOBS_CRAN_ALL"), "1")) {
  files <- readLines("cran-files.txt")
  files <- files[nzchar(files) & !startsWith(files, "#")]
  stems <- gsub(".", "\\.", sub("^test[-_]", "", sub("\\.[rR]$", "", files)),
                fixed = TRUE)
  filter <- paste0("^(", paste(stems, collapse = "|"), ")$")
}

test_check("tulpaObs", filter = filter)
