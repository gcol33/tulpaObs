# =============================================================================
# test-no-new-tulpa-internals.R
# -- the unexported-tulpa-internal surface cannot regrow (#357)
#
# tulpaObs reaches tulpa internals it has no exported door for (`tulpa:::x`,
# `asNamespace("tulpa")`, `getFromNamespace(x, "tulpa")`) at a known, audited
# set of sites -- gcol33/tulpaObs#357's table. Those have no stability
# contract: a signature change in the engine breaks tulpaObs at run time with
# no version constraint able to catch it (#324 was exactly this). The fix is
# tracked upstream in gcol33/tulpa (exported doors / C-callables for each);
# most are not exported yet (`tulpa_hyper_draws()` and, since #810,
# `tulpa_pit(log_lik = )` for the LOO-PIT kernel are), so the remaining
# audited sites stay until each further door lands.
#
# What this guards is the ONE thing available right now: the surface cannot
# grow past what is already audited. A new call site -- a new file reaching
# an internal for the first time, or an existing file reaching MORE of them
# -- fails here instead of shipping unnoticed. It is deliberately keyed by
# FILE, not by line number: reformatting an existing call, or moving it a few
# lines, is not a regrowth and should not fail this test.
#
# When a door in gcol33/tulpa is exported and a call site here switches to
# it (#357's actual fix), LOWER that file's ceiling (or drop the row) in the
# SAME commit -- this is a ceiling, not a floor, and nothing should raise it
# except a genuinely new, justified internal reach.
# =============================================================================

# Per-file ceiling on unexported-tulpa-internal call SITES (audited count as
# of #357, commit fed61c1). Comment-only lines (roxygen docs, prose) are
# excluded, same as the issue's own `grep -v ':[0-9]*: *#'` reproduce recipe.
.TOBS_TULPA_INTERNAL_CEILING <- c(
  "areal_bfgs.R"                = 2L,
  "ccd_outer.R"                 = 1L,
  "community_em.R"              = 1L,
  "cover_hurdle_joint.R"        = 1L,
  "cover_hurdle_joint_decode.R" = 1L,
  "dyn_abun.R"                  = 2L,
  "em_laplace_re.R"             = 1L,
  "fp_occu.R"                   = 1L,
  "joint_postprocess_shared.R"  = 2L,
  "joint_substrate.R"           = 3L,
  "ms_count_pg_gibbs.R"         = 1L,
  "ms_dyn_occu_pg_gibbs.R"      = 1L,
  "ms_int_occu_pg_gibbs.R"      = 1L,
  "ms_occu_pg_gibbs.R"          = 1L,
  "ms_occu_spatial.R"           = 2L,
  "nmix_laplace_re_spatial.R"   = 2L,
  "nmix_laplace_spde.R"         = 2L,
  "occu_cover_batch.R"          = 1L,
  "occu_cover_nuts.R"           = 6L,
  "occu_multiscale_cover.R"     = 1L,
  "occu_pg_gibbs.R"             = 2L,
  "t_occu.R"                    = 1L)

# Per-file count of code lines (excludes comment-only lines) reaching a
# tulpa internal via `tulpa:::`, `asNamespace("tulpa")` or
# `getFromNamespace(*, "tulpa")`.
.tobs_scan_tulpa_internals <- function(pkg_root) {
  pat <- "tulpa:::|asNamespace\\(\"tulpa\"\\)|getFromNamespace"
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

test_that("the tulpa-internal call surface has not grown past #357's audit", {
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
  ceiling  <- .TOBS_TULPA_INTERNAL_CEILING

  # A file reaching a tulpa internal for the FIRST time, not in the audited
  # table at all.
  new_files <- setdiff(names(observed), names(ceiling))
  expect_identical(
    new_files, character(0),
    info = paste0(
      "New file(s) calling a tulpa internal (tulpa:::, asNamespace(\"tulpa\"), ",
      "or getFromNamespace(*, \"tulpa\")) with no exported door and no entry ",
      "in #357's audited ceiling: ", paste(new_files, collapse = ", "),
      ". If tulpa now exports a door for this, use it instead. If not, this ",
      "grows the internal surface #357 tracks -- add the file to ",
      "gcol33/tulpaObs#357's table and to .TOBS_TULPA_INTERNAL_CEILING here, ",
      "with the reason."))

  # An audited file reaching MORE internals than its recorded ceiling.
  over <- observed[names(observed) %in% names(ceiling) &
                     observed > ceiling[names(observed)]]
  expect_identical(
    unname(over), integer(0),
    info = paste0(
      "File(s) with more tulpa-internal call sites than #357's audited ",
      "ceiling: ", paste(sprintf("%s (%d > %d)", names(over), over,
                                 ceiling[names(over)]), collapse = ", "),
      ". Raise the ceiling here only for a genuinely new, justified internal ",
      "reach -- not to silence this test."))

  # A file that no longer reaches ANY tulpa internal (the door it needed got
  # exported, or the code was removed) is exactly #357's fix landing -- not a
  # failure. Lower or drop its row here in the same commit as that change.
})
