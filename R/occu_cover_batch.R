# Batched multi-response occu_cover (gcol33/tulpa#66).
#
# Fitting occu_cover per species at EVA scale re-pays a cost that is identical
# across species: the site x visit structure (cells, sites, visits, detection
# design) is species-invariant; only the response y / y_pos differs. The B
# species' latent blocks are independent, so the joint Hessian is block-diagonal
# and each species' fit is statistically independent of the others.
#
# This file implements the *batched-independent* path: one tobs() call fits B
# species, each with exactly the per-species model it would get from a separate
# fit. It is distinct from ms_occu_cover() (the community model, which borrows
# strength via Gaussian community priors): here there is no cross-species term.
#
# Backend staging (see tulpa/dev_notes/design_batched_multiresponse_joint.md):
#   * Stage 1 (this file): per-species loop reusing the full single-species
#     pipeline. Correct + per-species bit-identical to independent fits by
#     construction; it is the validation oracle and the production driver
#     skeleton (parallelisable across machines, checkpointable). It does NOT
#     amortise the C++ occupancy-mixture scatter -- that is the bandwidth win of
#     the fused backend, which slots behind this same public API.
#   * Later stages swap in the fused block-diagonal C++ solver with no API
#     change; the 2-species equivalence test gates that swap.

# Number of responses B encoded in `y`, or NULL when `y` is an ordinary
# single-species response. Multi-response is either a list of >= 2 matrices or a
# 3D array [n_sites x max_visits x B]. A data.frame (is.list TRUE) or a plain
# matrix is single-response.
.tobs_multiresponse_n <- function(y) {
  if (is.null(y)) return(NULL)
  if (is.array(y) && length(dim(y)) == 3L) {
    return(dim(y)[3L])
  }
  if (is.list(y) && !is.data.frame(y)) {
    if (length(y) < 1L) return(NULL)
    if (!all(vapply(y, function(e) is.matrix(e) || is.data.frame(e),
                    logical(1)))) {
      return(NULL)
    }
    return(length(y))
  }
  NULL
}

# Slice response `s` (1-based) out of a multi-response `y` (list or 3D array) as
# a single n_sites x max_visits matrix.
.tobs_response_slice <- function(y, s) {
  if (is.array(y) && length(dim(y)) == 3L) {
    return(y[, , s, drop = TRUE])
  }
  as.matrix(y[[s]])
}

# Resolve species labels for a batch of B responses: explicit `species`, else
# the names carried on a `y` list, else sp1..spB. Validated to length B + unique.
.tobs_batch_species_labels <- function(species, y, B) {
  if (!is.null(species)) {
    labs <- as.character(species)
  } else if (is.list(y) && !is.null(names(y)) && all(nzchar(names(y)))) {
    labs <- names(y)
  } else {
    labs <- paste0("sp", seq_len(B))
  }
  if (length(labs) != B) {
    stop(sprintf(paste0("Batched occu_cover: `species` has length %d but `y` ",
                        "carries %d responses."), length(labs), B),
         call. = FALSE)
  }
  if (anyDuplicated(labs)) {
    stop("Batched occu_cover: `species` labels must be unique.", call. = FALSE)
  }
  labs
}

# Fit B species, each with the per-species occu_cover model, by replaying the
# full single-species tobs() pipeline once per response. Returns a `tobs_batch`.
#
# `tobs_args` is the captured argument list of the originating tobs() call
# (formula/data/family/detection/visits/method/priors/control plus `...`),
# already shorn of `y`; this driver inserts the per-species `y` / `y_pos`.
.tobs_fit_occu_cover_batch <- function(tobs_args, y, B) {
  dots     <- tobs_args$dots
  y_pos    <- dots$y_pos
  species  <- dots$species

  B_pos <- .tobs_multiresponse_n(y_pos)
  if (is.null(B_pos)) {
    stop(paste0("Batched occu_cover: `y` is multi-response (", B,
                " species) so `y_pos` must be a matching list/3D array of the ",
                "same length."), call. = FALSE)
  }
  if (B_pos != B) {
    stop(sprintf(paste0("Batched occu_cover: `y` carries %d responses but ",
                        "`y_pos` carries %d."), B, B_pos), call. = FALSE)
  }

  labels <- .tobs_batch_species_labels(species, y, B)

  # Per-species `...`: drop the batch-only `species`, override `y_pos` with the
  # species slice. Everything else (positive, etc.) flows through unchanged.
  base_dots <- dots
  base_dots$species <- NULL

  fits <- vector("list", B)
  for (s in seq_len(B)) {
    sp_dots <- base_dots
    sp_dots$y_pos <- .tobs_response_slice(y_pos, s)
    call_args <- c(
      list(
        formula   = tobs_args$formula,
        data      = tobs_args$data,
        family    = tobs_args$family,
        detection = tobs_args$detection,
        y         = .tobs_response_slice(y, s),
        visits    = tobs_args$visits,
        method    = tobs_args$method,
        priors    = tobs_args$priors,
        control   = tobs_args$control
      ),
      sp_dots
    )
    fits[[s]] <- do.call(tobs, call_args)
  }
  names(fits) <- labels

  structure(
    list(
      fits      = fits,
      species   = labels,
      n_species = B,
      family    = tobs_args$family,
      method    = tobs_args$method
    ),
    class = "tobs_batch"
  )
}


# ---------------------------------------------------------------------------
# S3 surface for tobs_batch
# ---------------------------------------------------------------------------

#' Print a batched multi-response occu_cover fit.
#'
#' @param x A `tobs_batch` returned by [tobs()] on multi-response `y`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.tobs_batch <- function(x, ...) {
  cat(sprintf("<tobs_batch> %d species, family = %s\n",
              x$n_species, x$family$name))
  cat(sprintf("  species: %s\n",
              paste(x$species, collapse = ", ")))
  cat("  per-species fits in $fits[[<species>]]\n")
  invisible(x)
}

#' Per-species coefficients from a batched occu_cover fit.
#'
#' @param object A `tobs_batch` returned by [tobs()] on multi-response `y`.
#' @param ... Passed to each per-species `coef()`.
#' @return A named list of per-species coefficient vectors.
#' @export
coef.tobs_batch <- function(object, ...) {
  lapply(object$fits, stats::coef, ...)
}

#' Extract one species' fit from a batched occu_cover fit.
#'
#' @param x A `tobs_batch`.
#' @param species Species label (character) or index (integer).
#' @return The single-species `tobs_fit` for that species.
#' @export
tobs_batch_fit <- function(x, species) {
  if (!inherits(x, "tobs_batch")) {
    stop("tobs_batch_fit(): `x` must be a tobs_batch.", call. = FALSE)
  }
  x$fits[[species]]
}
