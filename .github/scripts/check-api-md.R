#!/usr/bin/env Rscript

# Assert that the machine-derivable parts of API.md still describe the package.
#
# README.md calls API.md "the full reference", and the parts of it that go stale
# mechanically are the ones nothing recomputes: the export list, the simulator
# table, the `tobs()` signature, and the sampler-control defaults. Each of those
# has one authority in the package -- NAMESPACE, the exported `simulate_*`
# names, `formals(tobs)` and `.TOBS_ENGINE_DEFAULTS` -- so the check reads that
# authority rather than a second copy of it.
#
# Prose, family notes and worked examples are deliberately not checked: they
# carry judgement a diff cannot supply.

api <- readLines("API.md", warn = FALSE)
txt <- paste(api, collapse = "\n")
fail <- character()

ns <- readLines("NAMESPACE", warn = FALSE)
exports <- grep("^export[(]", ns, value = TRUE)
exports <- sort(substr(exports, 8L, nchar(exports) - 1L))

# ---- 1. every export is named, and nothing named is not an export ----------
block <- txt
start <- regexpr("[*][*]Fitter[*][*]", block)
if (start < 0) {
  fail <- c(fail, "API.md has no export list (looked for a **Fitter** group)")
} else {
  stop_at <- regexpr("> Structured terms", substring(block, start))
  listing <- substring(block, start,
                       start + (if (stop_at > 0) stop_at - 2L else nchar(block)))
  named <- unique(regmatches(listing,
                             gregexpr("`[^`]+`", listing))[[1L]])
  named <- gsub("`", "", named)
  absent <- setdiff(exports, named)
  extra  <- setdiff(named, exports)
  if (length(absent))
    fail <- c(fail, sprintf("export list is missing %d export(s): %s",
                            length(absent), paste(absent, collapse = ", ")))
  if (length(extra))
    fail <- c(fail, sprintf("export list names %d non-export(s): %s",
                            length(extra), paste(extra, collapse = ", ")))
}

# ---- 2. every simulator has a row in the simulator table ------------------
sims <- grep("^simulate_", exports, value = TRUE)
no_row <- sims[!vapply(sims, function(s)
  grepl(paste0("`", s, "()`"), txt, fixed = TRUE), logical(1))]
if (length(no_row))
  fail <- c(fail, sprintf("%d simulator(s) have no table row: %s",
                          length(no_row), paste(no_row, collapse = ", ")))

# ---- 3. the tobs() signature block lists the real formals -----------------
# Checks 3 and 4 read the package itself, so they run only where it is
# installed: the smoke job, which installs before this step. The tarball job
# runs this beside the other source-only checks and gets 1 and 2.
if (requireNamespace("tulpaObs", quietly = TRUE)) {
  fm <- names(formals(tulpaObs::tobs))
  sig_start <- regexpr("tobs(formula", txt, fixed = TRUE)
  if (sig_start > 0) {
    sig <- substring(txt, sig_start, sig_start + 500L)
    sig <- substring(sig, 1L, regexpr("```", sig, fixed = TRUE) - 1L)
    gone <- fm[!vapply(fm, function(a)
      grepl(a, sig, fixed = TRUE), logical(1))]
    if (length(gone))
      fail <- c(fail, sprintf("tobs() signature block omits %s",
                              paste(gone, collapse = ", ")))
  }

  # ---- 4. the sampler-control defaults match the table --------------------
  d <- getFromNamespace(".tobs_engine_defaults", "tulpaObs")("nuts")
  for (k in names(d)) {
    row <- grep(sprintf("^[|] `%s`", k), api, value = TRUE)
    if (!length(row)) {
      fail <- c(fail, sprintf("control table has no row for `%s`", k))
      next
    }
    v <- format(d[[k]])
    if (!grepl(v, row[1L], fixed = TRUE))
      fail <- c(fail, sprintf("control table row for `%s` does not carry the resolved default %s: %s",
                              k, v, trimws(row[1L])))
  }
}

if (length(fail)) {
  cat("API.md is out of step with the package:\n\n")
  cat(paste0("  - ", fail), sep = "\n")
  cat("\nRegenerate the affected block, or add the row by hand.\n")
  quit(status = 1L)
}
cat(sprintf("API.md matches the package (%d exports, %d simulators).\n",
            length(exports), length(sims)))
