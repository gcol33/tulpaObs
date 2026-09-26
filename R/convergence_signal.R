# One convergence record per fit, and one signal when an optimiser stops short.
#
# Every family reaches the user through `tobs()`, and its fitters recorded
# convergence in two layouts: a `fit$convergence` list on all of them and a flat
# top-level `fit$converged` on six. `convergence()` / `converged()` are the
# public accessors and read the record, so the record is completed here from
# whichever layout the fitter wrote and the flat slot is removed: a fit carries
# one answer.
#
# The signal is written here for the same reason. A Laplace-family fit that
# ended without meeting its criterion used to warn on four families, print a
# line on two and say nothing on the rest; one warning at the fit tail reaches
# every family once, where a per-fitter warning fired once per species on the
# community routes that loop a single-species fitter.
#
# Only the optimiser routes are signalled. A sampler's `converged` is an Rhat
# verdict with its own diagnostics (`check_model()`, `convergence()`), not a
# stopped iteration budget, so the advice here does not apply to it.
.tobs_finalize_convergence <- function(fit, family = NULL) {
  rec <- fit$convergence
  if (!is.list(rec)) rec <- list()
  rec$converged <- rec$converged %||% fit$converged %||% NA
  rec$n_iter    <- rec$n_iter    %||% fit$n_iter    %||% NA_integer_
  fit$convergence <- rec
  fit$converged   <- NULL

  if (identical(rec$converged, FALSE) &&
      isTRUE(fit$method %in% c("laplace", "nested_laplace"))) {
    fam <- if (!is.null(family$name)) paste0(family$name, "()") else "the fit"
    n_it <- if (is.finite(suppressWarnings(as.numeric(rec$n_iter))))
      sprintf(" after %d iterations", as.integer(rec$n_iter)) else ""
    warning(sprintf(paste0(
      "%s: the %s fit stopped%s without meeting its convergence criterion, ",
      "so the estimates are not at the mode. Raise control$max.iter (or ",
      "loosen control$tol); convergence(fit) has the record."),
      fam, fit$method, n_it), call. = FALSE)
  }
  fit
}
