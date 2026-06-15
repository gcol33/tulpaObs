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

  # Backend selection. The looped backend (B independent single-species fits) is
  # the DEFAULT: it is correct by construction and as fast as possible. The fused
  # block-diagonal backend (gcol33/tulpa#66) is correct (bit-identical per
  # species, gated by test-occu-cover-batch.R) but delivers no measured speed
  # benefit for occu_cover -- the per-species sparse factorization dominates and
  # is not amortizable, so it is at best parity and slower than looped at large
  # fields (dev_notes/_probe_batch_bsweep.R). It is reachable via
  # control$batch.backend = "fused" for experimentation; the sparse-native
  # variant that would be needed to make it competitive is tracked as an open
  # issue (gcol33/tulpa#69). `.tobs_fit_occu_cover_batch_fused` returns NULL when
  # the configuration is not fused-eligible, falling through to the looped path.
  backend <- tobs_args$control[["batch.backend"]] %||% "looped"
  if (identical(backend, "fused")) {
    fused <- .tobs_fit_occu_cover_batch_fused(tobs_args, y, y_pos, B, labels)
    if (!is.null(fused)) return(fused)
  }

  # Looped backend. `batch.backend` is a batch-orchestration knob, not a
  # single-species control key, so strip it before the per-species fits (it would
  # otherwise be rejected by .tobs_validate_control). Per-species `...`: drop the
  # batch-only `species`, override `y_pos` with the species slice; everything else
  # (positive, etc.) flows through unchanged.
  sp_control <- tobs_args$control
  sp_control[["batch.backend"]] <- NULL
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
        control   = sp_control
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
      method    = tobs_args$method,
      backend   = "looped"
    ),
    class = "tobs_batch"
  )
}


# Fused block-diagonal backend (gcol33/tulpa#66). Runs B species through ONE
# multi-block nested-Laplace solve: the species share the design + sparsity
# pattern, their latent systems are block-diagonal, and the fused cell-coupling
# scatter loads each design row once and loops species inner. Per-species
# trajectory is bit-identical to an independent single-species fit (the fused
# path only reorganises the work), so each species post-processes to the same
# tobs_fit a looped fit produces.
#
# Returns a `tobs_batch` (backend = "fused"), or NULL when the configuration is
# not fused-eligible -- the caller then falls back to the looped path. Eligible:
# spatial nested-Laplace on the default joint_coupled engine with a FIXED pos-arm
# dispersion (no latent cover RE, no phi.grid.pos), i.e. the common occu_cover
# spatial fit. The fused driver integrates a single shared FIXED outer grid
# across species (per-species adaptive refinement is inherently not shareable),
# so a fused fit equals an adaptive single-species fit only with adaptive grid
# off; the equivalence gate fixes both sides' grid.
.tobs_fit_occu_cover_batch_fused <- function(tobs_args, y, y_pos, B, labels) {
  if (!identical(tobs_args$method, "nested_laplace")) return(NULL)
  engine_pick <- tobs_args$control[["engine"]] %||% "joint_coupled"
  if (!identical(engine_pick, "joint_coupled")) return(NULL)

  dots <- tobs_args$dots

  # Collect per-species prep by replaying the dispatch in collect mode. This
  # reuses ALL of .dispatch_occu_cover + the joint_coupled Part-A builder (model
  # construction, field resolution, arm priors, sigma_pos pre-fit, grids); no
  # model-building logic is duplicated here. A species whose dispatch does not
  # return an `occu_cover_jc_prep` (non-spatial, v2/v3, an error) is ineligible.
  preps <- vector("list", B)
  for (s in seq_len(B)) {
    sp_dots          <- dots
    sp_dots$species  <- NULL
    sp_dots$y_pos    <- .tobs_response_slice(y_pos, s)
    ctrl_s           <- tobs_args$control
    ctrl_s$.batch_collect <- TRUE
    prep <- tryCatch(
      do.call(.dispatch_occu_cover, c(list(
        formula   = tobs_args$formula, data = tobs_args$data,
        family    = tobs_args$family,  detection = tobs_args$detection,
        y         = .tobs_response_slice(y, s), visits = tobs_args$visits,
        engine    = "nested_laplace", priors = tobs_args$priors,
        control   = ctrl_s), sp_dots)),
      error = function(e) e)
    if (!inherits(prep, "occu_cover_jc_prep")) return(NULL)
    preps[[s]] <- prep
  }

  # Fused-eligible only with a fixed pos-arm dispersion: a latent cover RE or an
  # explicit phi.grid.pos puts sigma on the outer grid as a per-arm phi axis,
  # which the batched driver does not carry.
  ineligible <- vapply(preps, function(p)
    !is.null(p$fit_call$phi_grid) || isTRUE(p$is_latent), logical(1))
  if (any(ineligible)) return(NULL)

  fc1       <- preps[[1L]]$fit_call
  arms1     <- fc1$responses
  n_arms    <- length(arms1)
  spec_name <- preps[[1L]]$spec_name
  has_trend <- isTRUE(preps[[1L]]$has_trend)

  # Per-data-arm species-column response matrix; per-arm per-species dispersion.
  y_batch <- vector("list", n_arms)
  for (k in seq_len(n_arms)) {
    yk <- arms1[[k]]$y
    if (is.null(yk) || length(yk) == 0L) next
    y_batch[[k]] <- do.call(cbind, lapply(preps, function(p)
      as.numeric(p$fit_call$responses[[k]]$y)))
  }
  phi_batch <- matrix(0, n_arms, B)
  for (k in seq_len(n_arms)) {
    for (s in seq_len(B)) {
      phi_batch[k, s] <- preps[[s]]$fit_call$responses[[k]]$phi %||% 1
    }
  }

  bat <- tulpa:::tulpa_nl_joint_batch(
    responses     = arms1, prior = fc1$prior, copy = fc1$copy,
    n_batch       = B, y_batch = y_batch, phi_batch = phi_batch,
    max_iter      = as.integer(fc1$control$max_iter %||% 200L),
    tol           = as.numeric(fc1$control$tol %||% 1e-6),
    cell_coupling = spec_name, store_Q = TRUE)

  arm_layout <- bat$arm_layout
  theta_grid <- bat$theta_grid
  # The single-field multi-block grid carries b1.-prefixed axis names; Part B's
  # no-trend branch reads bare "sigma"/"alpha" (the single-block convention).
  # Strip the single block's prefix so the hyperparameter summary resolves.
  if (!has_trend && !is.null(colnames(theta_grid))) {
    colnames(theta_grid) <- sub("^b1\\.", "", colnames(theta_grid))
  }

  fits <- lapply(seq_len(B), function(s) {
    ps <- bat$per_species[[s]]
    engine_fit <- list(
      arm_layout       = arm_layout,
      theta_grid       = theta_grid,
      log_marginal     = ps$log_marginal,
      weights          = ps$weights,
      modes            = ps$modes,
      Q_csc_p_per_grid = ps$Q_csc_p_per_grid,
      Q_csc_i_per_grid = ps$Q_csc_i_per_grid,
      Q_csc_x_per_grid = ps$Q_csc_x_per_grid,
      Q_csc_n          = ps$Q_csc_n
    )
    .occu_cover_jc_postprocess(engine_fit, preps[[s]]$ctx)
  })
  names(fits) <- labels

  structure(
    list(
      fits      = fits,
      species   = labels,
      n_species = B,
      family    = tobs_args$family,
      method    = tobs_args$method,
      backend   = "fused"
    ),
    class = "tobs_batch"
  )
}


# ---------------------------------------------------------------------------
# by = "<species_col>": per-species batched fit from a long plot-level frame
# ---------------------------------------------------------------------------

# Build each species' response from a long / plot-level `data` and route the
# per-species responses through the batched-independent driver, returning a
# `tobs_batch`. The species share the design (cell-level covariates, the spatial
# graph, the detection / cover formulas, the visit grid); only the per-species
# response differs -- exactly the batched-independent contract.
#
# REUSE: the long -> response pivot is `tobs_data()`'s, called once per species
# onto ONE shared site x visit grid (sites = first-appearance order over the FULL
# frame, visits = sorted unique). The site-level covariate frame the dispatchers
# need is the first row per site, taken the way tobs_data(occ.covs = ) does. The
# fitting itself is the existing batch driver (occu_cover) or per-species
# single-species tobs() calls (cover) -- no model-building logic is duplicated.
#
# `dots` carries the originating tobs() `...`. The long -> matrix column names
# are read from it: `site`, `visit`, `response` (the detection 0/1 column for
# occu_cover, the cover column for cover), `y_pos` (the cover column, occu_cover
# only), and `det.covs` (visit-level covariate columns carried to `visits`).
.tobs_fit_by_species <- function(formula, data, family, detection, visits,
                                 method, priors, control, by, dots) {
  fam <- family$name
  if (!fam %in% c("occu_cover", "cover")) {
    stop(sprintf(paste0(
      "tobs(by = ): per-species batched fitting is wired for occu_cover() and ",
      "cover() (a per-plot response a species column splits), not %s(). ",
      "Fit %s() one species at a time, or for a community model that pools ",
      "across species use ms_occu_cover()."), fam, fam), call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("tobs(by = ): `data` must be a long / plot-level data frame.",
         call. = FALSE)
  }
  if (length(by) != 1L || !is.character(by) || !by %in% names(data)) {
    stop(sprintf("tobs(by = '%s'): `by` must name a single column of `data`.",
                 by), call. = FALSE)
  }
  # `visits` (the tobs() formal) carries the long-format VISIT column name in
  # by = mode for occu_cover: the user writes `visit = "<col>"`, which partial-
  # matches the `visits` formal, so it arrives here as a length-1 string. (The
  # visit-level covariate matrices are built internally from the long `data` via
  # `det.covs = `, not handed in pre-pivoted.) A non-string `visits` here is a
  # mis-supplied pre-built visit grid.
  visit <- NULL
  if (!is.null(visits)) {
    if (is.character(visits) && length(visits) == 1L) {
      visit <- visits
    } else {
      stop("tobs(by = ): pass the long-format visit column as `visit = ",
           "\"<col>\"` (a column name); visit-level covariates are built from ",
           "`data` via `det.covs = `, not handed in as a pre-pivoted grid.",
           call. = FALSE)
    }
  }

  site     <- dots$site
  response <- dots$response
  det.covs <- dots$det.covs
  if (is.null(site) || is.null(response)) {
    stop("tobs(by = ): supply `site = ` and `response = ` (the long-format ",
         "site identifier and response column names).", call. = FALSE)
  }
  for (col in c(site, response, det.covs)) {
    if (!col %in% names(data))
      stop(sprintf("tobs(by = ): column '%s' not found in `data`.", col),
           call. = FALSE)
  }

  sp_raw <- data[[by]]
  labels <- as.character(unique(sp_raw))
  B <- length(labels)
  if (B < 1L) stop("tobs(by = ): no species found in `data`.", call. = FALSE)

  if (identical(fam, "cover")) {
    return(.tobs_fit_cover_by(
      formula = formula, data = data, family = family, method = method,
      priors = priors, control = control,
      by = by, site = site, response = response, labels = labels,
      sp_raw = sp_raw, dots = dots))
  }

  # occu_cover: a site x visit response. Detection arm = `response` (0/1), cover
  # arm = `y_pos` (continuous), both pivoted onto the shared grid. Visit-level
  # covariates flow to `visits` as the det.covs matrices tobs_data() returns.
  y_pos <- dots$y_pos
  if (is.null(visit) || is.null(y_pos)) {
    stop("tobs(by = ) for occu_cover(): supply `visit = ` (the long-format ",
         "visit/replicate column) and `y_pos = ` (the cover column).",
         call. = FALSE)
  }
  for (col in c(visit, y_pos)) {
    if (!col %in% names(data))
      stop(sprintf("tobs(by = ): column '%s' not found in `data`.", col),
           call. = FALSE)
  }

  # One canonical site x visit grid over the FULL frame, so every species'
  # response rows align with each other AND with the shared cell-level design.
  sites_lvl  <- unique(data[[site]])
  visits_lvl <- sort(unique(data[[visit]]))

  # Shared cell-level covariate frame: first row per site over the full frame.
  # The dispatchers read X_occ / X_det / X_pos from this; it is species-invariant
  # because the design covariates do not vary by species (only the response does).
  cell_rows <- match(sites_lvl, data[[site]])
  cell_data <- data[cell_rows, , drop = FALSE]
  rownames(cell_data) <- NULL

  # Shared visit-level covariate grid (the det.covs matrices). The visit-level
  # covariates are plot attributes (one value per site x visit cell), shared
  # across the species observed at that plot, so the grid is built from the FULL
  # frame -- one record per (site, visit) cell -- NOT from a single species,
  # which would leave a species' unvisited cells NA in the shared design.
  visits_sp <- NULL
  if (!is.null(det.covs)) {
    cell_first <- !duplicated(data[c(site, visit)])
    vdf <- data[cell_first, , drop = FALSE]
    visits_sp <- tobs_data(vdf, y = response, site = site, visit = visit,
                           type = "occurrence", det.covs = det.covs,
                           sites = sites_lvl, visits = visits_lvl)$det.covs
  }

  y_list    <- vector("list", B)
  ypos_list <- vector("list", B)
  for (b in seq_len(B)) {
    sp_df <- data[sp_raw == labels[b], , drop = FALSE]
    od <- tobs_data(sp_df, y = response, site = site, visit = visit,
                    type = "occurrence",
                    sites = sites_lvl, visits = visits_lvl)
    op <- tobs_data(sp_df, y = y_pos, site = site, visit = visit,
                    type = "cover",
                    sites = sites_lvl, visits = visits_lvl)
    yp <- op$y
    yp[is.na(yp)] <- 0                 # cover meaningful only where y == 1
    y_list[[b]]    <- od$y
    ypos_list[[b]] <- yp
  }
  names(y_list) <- labels

  # The detection / cover formulas reference visit-level covariates by name; the
  # batch driver consumes `visits` (the det.covs matrices) the same way the
  # single-species path does. dots$y_pos is replaced by the per-species list;
  # site / visit / response / det.covs are batch-build keys, not fitter args.
  fit_dots <- dots
  fit_dots$site <- NULL; fit_dots$visit <- NULL
  fit_dots$response <- NULL; fit_dots$det.covs <- NULL
  fit_dots$y_pos <- ypos_list

  .tobs_fit_occu_cover_batch(
    tobs_args = list(formula = formula, data = cell_data, family = family,
                     detection = detection, visits = visits_sp,
                     method = method, priors = priors, control = control,
                     dots = fit_dots),
    y = y_list, B = B)
}

# Per-species cover() batch: each species is an independent single-species
# cover() fit on the shared cell-level design, looped through tobs() itself (so
# each fit is byte-identical to a separate cover() call). Returns a `tobs_batch`
# with the same shape as the occu_cover looped backend.
.tobs_fit_cover_by <- function(formula, data, family, method, priors, control,
                               by, site, response, labels, sp_raw, dots) {
  B <- length(labels)
  sites_lvl <- unique(data[[site]])
  cell_rows <- match(sites_lvl, data[[site]])
  cell_data <- data[cell_rows, , drop = FALSE]
  rownames(cell_data) <- NULL

  # Cover is single-replicate: one response value per site. Align each species'
  # cover vector to the shared site set; a site the species never reports is NA
  # (dropped from both arms by the cover encoder, i.e. not observed there).
  fit_dots <- dots
  fit_dots$site <- NULL; fit_dots$response <- NULL
  fit_dots$visit <- NULL; fit_dots$det.covs <- NULL; fit_dots$y_pos <- NULL

  fits <- vector("list", B)
  for (b in seq_len(B)) {
    sp_df <- data[sp_raw == labels[b], , drop = FALSE]
    yv <- sp_df[[response]][match(sites_lvl, sp_df[[site]])]
    call_args <- c(
      list(formula = formula, data = cell_data, family = family,
           y = as.numeric(yv), method = method, priors = priors,
           control = control),
      fit_dots)
    fits[[b]] <- do.call(tobs, call_args)
  }
  names(fits) <- labels

  structure(
    list(fits = fits, species = labels, n_species = B,
         family = family, method = method, backend = "looped"),
    class = "tobs_batch")
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
tobs_get <- function(x, species) {
  if (!inherits(x, "tobs_batch")) {
    stop("tobs_get(): `x` must be a tobs_batch.", call. = FALSE)
  }
  x$fits[[species]]
}
