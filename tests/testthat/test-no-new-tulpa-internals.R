# =============================================================================
# test-no-new-tulpa-internals.R
# -- the unexported-tulpa-internal surface stays empty (#357)
#
# An unexported engine symbol has no stability contract: a rename or a
# signature change in tulpa breaks tulpaObs at run time, and no version
# constraint can catch it (#324 was exactly this -- tulpaObs called the
# post-0.4.0 `.nl_grid_log_quad(refining =)` while DESCRIPTION pinned v0.4.0,
# 15 smoke errors). `.tulpa_iter_progress` then went the other way: tulpa
# renamed it to the exported `tulpa_iter_progress()` and every reach here
# would have errored had it not been swept.
#
# Every computation `R/` needs from the engine now has an exported door
# so the audited ceiling this file used to carry is
# gone and the budget is zero. A reach that reappears is either a door that
# exists and was not used, or a door that has to be filed in gcol33/tulpa
# before the call site can land.
#
# Tests are out of scope: they drive engine internals deliberately (compiled
# fixtures such as `cpp_test_log_prior_icar`, AGHQ objective gradients), which
# is what a fixture is for.
# =============================================================================

# Per-file count of code lines (excludes comment-only lines) reaching a
# tulpa internal via `tulpa:::`, `asNamespace("tulpa")` or
# `getFromNamespace(*, "tulpa")`.
.tobs_scan_tulpa_internals <- function(pkg_root) {
  pat <- "tulpa:::|asNamespace\\(\"tulpa\"\\)|getFromNamespace\\([^)]*\"tulpa\""
  files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                      full.names = TRUE)
  counts <- integer(0)
  for (f in files) {
    ln <- readLines(f, warn = FALSE)
    hit <- grepl(pat, ln) & !grepl("^\\s*#", ln)
    n <- sum(hit)
    if (n > 0) counts[basename(f)] <- n
  }
  counts
}

test_that("R/ reaches no unexported tulpa internal", {
  pkg_root <- system.file(package = "tulpaObs")
  # Under devtools::load_all() the installed tree is not this checkout's R/;
  # fall back to the package source root testthat runs from.
  if (!dir.exists(file.path(pkg_root, "R")) &&
      !file.exists(file.path(pkg_root, "R", "sysdata.rda"))) {
    pkg_root <- NULL
  }
  src_root <- testthat::test_path("..", "..")
  if (dir.exists(file.path(src_root, "R"))) pkg_root <- src_root
  skip_if(is.null(pkg_root), "package R/ source tree not found from this runner")

  observed <- .tobs_scan_tulpa_internals(pkg_root)

  expect_identical(
    sprintf("%s (%d)", names(observed), observed), character(0),
    info = paste0(
      "File(s) in R/ calling a tulpa internal (tulpa:::, ",
      "asNamespace(\"tulpa\"), or getFromNamespace(*, \"tulpa\")): ",
      paste(sprintf("%s (%d)", names(observed), observed), collapse = ", "),
      ". #357 emptied this surface -- every computation the engine supplies ",
      "has an exported door. If the one needed here has none, file it in ",
      "gcol33/tulpa and use the door, rather than reaching past it."))
})
