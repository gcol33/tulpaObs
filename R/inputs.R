# =============================================================================
# inputs.R — Single source of truth for tobs() response / site / visit inputs
#
# Every model family hand-derives its own design matrices, but two pieces of
# input handling must not drift between families:
#
#   * the policy that the site dimension of `y` matches the rows of the
#     site-level `data` -- the cross-check that occu / abun / removal / distance
#     / fp_occu / dyn_abun / the community and cover families each used to
#     hand-roll with a near-identical stop(); and
#   * the canonical (sites x visits, plus sources for the integrated families)
#     totals tobs() reports under control$verbose and stores on the fit.
#
# A `tobs_data` frame (R/data.R) is unpacked here into the (data, y, visits)
# triple the fitter consumes, so a pre-validated frame routes through exactly the
# same pipeline as raw arguments rather than a parallel path.
# =============================================================================


# Cross-check the site dimension of `y` against the rows of the site-level
# `data`. The single home for the check each family binder used to hand-roll.
# `y_unit` keeps the message faithful to the response shape: a 2D response counts
# "rows", a 3D array / per-source list counts "sites", and the cover hurdle's
# response vector counts "values".
.tobs_check_site_count <- function(n_y, n_data, y_unit = "rows") {
  if (n_y != n_data) {
    stop(sprintf("y has %d %s but data has %d rows", n_y, y_unit, n_data),
         call. = FALSE)
  }
  invisible(NULL)
}


# Canonical response totals for a fit, read from the response shape alone so one
# reader backs the verbose fit message and the stored `fit$dims`. `n_sites` is
# the site count, `max_visits` the replicate dimension, and `n_sources` the
# number of per-source response matrices for the integrated families (NA when
# the shape does not carry that dimension).
.tobs_input_dims <- function(y) {
  out <- list(n_sites = NA_integer_, max_visits = NA_integer_,
              n_sources = NA_integer_)
  if (is.null(y)) return(out)

  # Integrated families: a list of per-source responses (matrix or 3D array).
  if (is.list(y) && !is.array(y)) {
    out$n_sources <- length(y)
    first <- if (length(y)) y[[1]] else NULL
    d <- dim(first)
    if (!is.null(d)) {
      out$n_sites    <- d[1]
      out$max_visits <- d[2]
    } else if (!is.null(first)) {
      out$n_sites <- length(first)
    }
    return(out)
  }

  d <- dim(y)
  if (is.null(d)) {            # bare response vector (cover hurdle)
    out$n_sites <- length(y)
    return(out)
  }
  out$n_sites    <- d[1]
  out$max_visits <- d[2]
  out
}


# Build the one-line "fitting <family> on N sites x J visits" message emitted
# before dispatch when control$verbose is TRUE. Reports only the dimensions that
# are unambiguous from the response shape (sites, replicate visits, sources) --
# it does not guess whether a third array dimension is species or seasons.
.tobs_input_message <- function(family_name, n_sites, max_visits, n_sources) {
  parts <- sprintf("%d site%s", n_sites, if (identical(n_sites, 1L)) "" else "s")
  if (!is.null(max_visits) && !is.na(max_visits)) {
    parts <- c(parts, sprintf("%d visit%s", max_visits,
                              if (identical(as.integer(max_visits), 1L)) "" else "s"))
  }
  body <- paste(parts, collapse = " x ")
  if (!is.null(n_sources) && !is.na(n_sources)) {
    body <- sprintf("%s, %d data source%s", body, n_sources,
                    if (identical(as.integer(n_sources), 1L)) "" else "s")
  }
  sprintf("tobs(): fitting %s on %s.", family_name, body)
}


# Unpack a `tobs_data` frame into the (data, y, visits) triple tobs() fits from.
# The frame's site-level `occ.covs` becomes `data`, its response `y` becomes `y`,
# and its visit-level `det.covs` (a named list of [n_sites x max_visits]
# matrices) becomes `visits` -- the exact shape .normalize_visits() already
# consumes. Passing `y =` / `visits =` alongside a frame is ambiguous and errors.
.tobs_unpack_frame <- function(frame, y, visits) {
  if (!is.null(y)) {
    stop("Pass the response through the `tobs_data` frame, not `y =`, when ",
         "`data` is a frame.", call. = FALSE)
  }
  if (!is.null(visits)) {
    stop("Pass visit-level covariates through the frame's `det.covs`, not ",
         "`visits =`, when `data` is a frame.", call. = FALSE)
  }

  fy <- frame$y
  if (is.null(fy)) {
    stop("`tobs_data` frame carries no response `y`.", call. = FALSE)
  }

  d <- dim(fy)
  n_sites <- if (is.null(d)) length(fy) else d[1]

  data <- frame$occ.covs
  if (is.null(data)) {
    # No site-level covariates: a placeholder frame with the right row count so
    # an intercept-only design (model.matrix(~ 1, data)) binds correctly.
    data <- data.frame(.tobs_site = seq_len(n_sites))
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  list(data = data, y = fy, visits = frame$det.covs)
}
