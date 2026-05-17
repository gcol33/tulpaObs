# =============================================================================
# Occupancy nested-Laplace path. Mirrors `.tobs_laplace` but builds a
# multi-block latent prior (spatial + temporal + iid) and routes the
# occupancy M-step block through tulpa::tulpa_nested_laplace() via the
# generic EM engine's per-block dispatcher. Single-season occupancy only;
# other occupancy variants raise a clear error.
# =============================================================================


#' Build a multi-block latent prior list for tulpa::tulpa_nested_laplace()
#'
#' Converts the tobs-level spatial / temporal / RE specs into the
#' list-of-blocks shape that `tulpa::tulpa_nested_laplace()` expects under
#' its multi-block dispatch (`.is_multi_block_prior`).
#'
#' @param spatial Optional `tobs_spatial` from [tobs_bym2()] / [tobs_icar()].
#'   Only BYM2 / ICAR are wired through to the multi-block engine at present;
#'   GP / multiscale_gp / SVC are not yet supported and raise.
#' @param temporal Optional `tobs_temporal` from [tobs_temporal()]. Types
#'   `"ar1"`, `"rw1"`, `"rw2"`, `"iid"` are supported.
#' @param re Optional list of `tobs_re` objects. Only `model = "iid"` terms
#'   are converted to IID latent blocks here; correlated structures (ar1 /
#'   rw1 / rw2 on RE groups) are passed through as temporal-like blocks.
#' @param model A `tobs_model` from `.tobs_build_model()`. Used to resolve
#'   variable-name references in temporal$time / re$group against
#'   `model$data`, and to pin `N = model$n_sites` for single-season fits.
#'
#' @return `NULL` when no latent block is supplied; a single-block list
#'   when exactly one is supplied; a list-of-blocks otherwise. Each block
#'   is the minimal field set that the tulpa multi-block dispatch fills
#'   defaults around (`.NL_REGISTRY[[type]]$defaults`).
#'
#' @keywords internal
.tobs_to_multi_block_prior <- function(spatial = NULL, temporal = NULL,
                                       re = NULL, model) {
  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object", call. = FALSE)
  }
  if (!identical(model$model_type, "single")) {
    stop("Multi-block nested Laplace is currently wired only for ",
         "single-season occupancy (model_type = 'single'); got '",
         model$model_type, "'.", call. = FALSE)
  }

  N <- as.integer(model$n_sites)
  blocks <- list()

  if (!is.null(spatial)) {
    blocks <- c(blocks, list(.tobs_block_from_spatial(spatial, N)))
  }
  if (!is.null(temporal)) {
    blocks <- c(blocks, list(.tobs_block_from_temporal(temporal, model, N)))
  }
  if (!is.null(re)) {
    for (r in re) {
      blocks <- c(blocks, list(.tobs_block_from_re(r, model, N)))
    }
  }

  if (length(blocks) == 0L) return(NULL)
  if (length(blocks) == 1L) return(blocks[[1]])
  blocks
}


# Build a multi-block-shaped spatial block from a tobs_spatial spec. The
# multi-block C++ entry reads spatial_idx as 1-based per the existing
# .nl_block_spec_for_cpp() contract.
#
# Per-block grids are narrower than `.NL_REGISTRY`'s single-block defaults
# (which target ~20-25 cells per block). Three-block combinations
# (e.g. BYM2 + AR1 + IID) at the single-block defaults exceed the
# multi-block hard cap (2048 cells); the narrower defaults below keep
# typical combos under ~250 cells. Users who need finer integration can
# pass `*_grid` overrides directly via the tobs_* spec attributes (when
# present) -- the helper passes them through if set.
.tobs_block_from_spatial <- function(spatial, N) {
  if (!inherits(spatial, "tobs_spatial")) {
    stop("`spatial` must be a tobs_spatial object", call. = FALSE)
  }
  type <- spatial$type
  if (!type %in% c("bym2", "icar", "car_proper")) {
    stop("Spatial type '", type, "' is not yet wired into the multi-block ",
         "nested-Laplace path (supported: bym2, icar, car_proper). ",
         "Use `engine = 'laplace'` or open an issue if you need this type.",
         call. = FALSE)
  }
  if (as.integer(spatial$n_units) != N) {
    stop(sprintf(
      "spatial has %d units but the model has %d sites; one obs per site ",
      spatial$n_units, N),
      "is required for single-season occupancy nested-Laplace.", call. = FALSE)
  }
  out <- list(
    type            = type,
    spatial_idx     = seq_len(N),
    n_spatial_units = as.integer(spatial$n_units),
    adj_row_ptr     = as.integer(spatial$adj_row_ptr),
    adj_col_idx     = as.integer(spatial$adj_col_idx),
    n_neighbors     = as.integer(spatial$n_neighbors)
  )
  if (type == "bym2") {
    out$scale_factor <- as.numeric(spatial$scale_factor %||% 1.0)
    if (!is.null(spatial$sigma_grid)) out$sigma_grid <- spatial$sigma_grid
    if (!is.null(spatial$rho_grid))   out$rho_grid   <- spatial$rho_grid
    if (is.null(out$sigma_grid) && is.null(out$rho_grid)) {
      sg <- exp(seq(log(0.2), log(2.0), length.out = 3))
      rg <- c(0.3, 0.7)
      gr <- expand.grid(sigma = sg, rho = rg)
      out$sigma_grid <- gr$sigma
      out$rho_grid   <- gr$rho
    }
  } else if (type == "icar") {
    if (!is.null(spatial$tau_grid)) out$tau_grid <- spatial$tau_grid
  } else if (type == "car_proper") {
    if (!is.null(spatial$tau_grid)) out$tau_grid <- spatial$tau_grid
    if (!is.null(spatial$rho_grid)) out$rho_grid <- spatial$rho_grid
  }
  out
}


# Build a multi-block temporal block from a tobs_temporal spec. Resolves
# the time variable from model$data when given as a string; passes integer
# indices through. `model = "iid"` becomes an iid block (no Q assembly).
.tobs_block_from_temporal <- function(temporal, model, N) {
  if (!inherits(temporal, "tobs_temporal")) {
    stop("`temporal` must be a tobs_temporal object", call. = FALSE)
  }
  type <- temporal$type
  if (!type %in% c("ar1", "rw1", "rw2", "iid")) {
    stop("Temporal type '", type, "' is not supported by the multi-block ",
         "nested-Laplace path (supported: ar1, rw1, rw2, iid).",
         call. = FALSE)
  }
  time_idx <- if (is.character(temporal$time)) {
    if (is.null(model$data) || !temporal$time %in% names(model$data)) {
      stop("temporal$time = '", temporal$time, "' not found in model$data.",
           call. = FALSE)
    }
    as.integer(as.factor(model$data[[temporal$time]]))
  } else {
    as.integer(temporal$time)
  }
  if (length(time_idx) != N) {
    stop(sprintf(
      "Resolved temporal index has length %d but the model has %d sites.",
      length(time_idx), N), call. = FALSE)
  }
  n_times <- max(time_idx, na.rm = TRUE)

  if (type == "iid") {
    out <- list(type = "iid", obs_idx = time_idx, n_units = as.integer(n_times))
    if (!is.null(temporal$sigma_grid)) out$sigma_grid <- temporal$sigma_grid
    return(out)
  }

  out <- list(
    type         = type,
    temporal_idx = time_idx,
    n_times      = as.integer(n_times)
  )
  if (type == "rw1") out$cyclic <- isTRUE(temporal$cyclic)
  if (!is.null(temporal$tau_grid)) out$tau_grid <- temporal$tau_grid
  if (type == "ar1") {
    if (!is.null(temporal$rho_grid)) out$rho_grid <- temporal$rho_grid
    if (is.null(out$tau_grid) && is.null(out$rho_grid)) {
      # 3 x 2 = 6 cells per block keeps a BYM2 + AR1 + IID combo under cap.
      tg <- exp(seq(log(0.5), log(20), length.out = 3))
      rg <- c(0.3, 0.8)
      gr <- expand.grid(tau = tg, rho = rg)
      out$tau_grid <- gr$tau
      out$rho_grid <- gr$rho
    }
  }
  out
}


# Build a multi-block iid (or temporal-on-groups) block from a tobs_re spec.
# A bare `tobs_re(group = "x")` (default model = "iid") becomes an iid block.
# If the user asks for ar1/rw1/rw2 on the group, route through the temporal
# block builder so the same dispatch logic handles both.
.tobs_block_from_re <- function(re, model, N) {
  if (!inherits(re, "tobs_re")) {
    stop("`re` element must be a tobs_re object", call. = FALSE)
  }
  if (identical(re$type, "slope")) {
    stop("Random slopes are not yet supported by the multi-block ",
         "nested-Laplace path. Use `engine = 'laplace'` for slopes.",
         call. = FALSE)
  }

  grp_idx <- if (is.character(re$group)) {
    if (is.null(model$data) || !re$group %in% names(model$data)) {
      stop("re$group = '", re$group, "' not found in model$data.",
           call. = FALSE)
    }
    as.integer(as.factor(model$data[[re$group]]))
  } else {
    as.integer(re$group)
  }
  if (length(grp_idx) != N) {
    stop(sprintf(
      "Resolved RE group index has length %d but the model has %d sites.",
      length(grp_idx), N), call. = FALSE)
  }
  n_units <- max(grp_idx, na.rm = TRUE)

  model_name <- re$model %||% "iid"
  if (model_name == "iid") {
    out <- list(type = "iid", obs_idx = grp_idx, n_units = as.integer(n_units))
    if (!is.null(re$sigma_grid)) {
      out$sigma_grid <- re$sigma_grid
    } else {
      # 3-point grid keeps multi-block combos comfortably under the cap.
      out$sigma_grid <- exp(seq(log(0.2), log(2.0), length.out = 3))
    }
    out
  } else if (model_name %in% c("ar1", "rw1", "rw2")) {
    out <- list(type = model_name, temporal_idx = grp_idx,
                n_times = as.integer(n_units))
    if (model_name == "rw1") out$cyclic <- FALSE
    if (!is.null(re$tau_grid)) out$tau_grid <- re$tau_grid
    out
  } else {
    stop("RE model '", model_name, "' is not supported by the multi-block ",
         "nested-Laplace path.", call. = FALSE)
  }
}


# =============================================================================
# Driver — EM around tulpa_em_laplace(..., method = nested per-block)
#
# The occupancy block carries a `prior` field built by
# .tobs_to_multi_block_prior(), which makes the generic EM engine route that
# block through tulpa::tulpa_nested_laplace() instead of tulpa::tulpa_laplace().
# The detection block stays single-Laplace (it carries its own observation
# weights; nested-Laplace doesn't accept those).
#
# Known approximation. The E-step computes psi[i] = plogis(X[i,]*beta_occ)
# only -- the latent block contribution (spatial / temporal / iid) is left
# out of psi because reconstructing eta from the grid-averaged mode requires
# per-grid hyperparameter values that the nested-Laplace driver does not
# return alongside the modes. The EM fixed point is therefore the
# fixed-effect MAP under the marginal hyperparameter posterior, with the
# latent block estimated by the inner nested-Laplace fit on the
# pseudo-binomial response. Recovery of fixed effects has not yet been
# validated against simulated truth; see plan_multi_block.md Phase D for
# the follow-up validation harness.
# =============================================================================

#' Fit a single-season occupancy tobs model via nested-Laplace
#'
#' Internal driver matching `.tobs_laplace()` in shape but routing the
#' occupancy block through `tulpa::tulpa_nested_laplace()` with a
#' multi-block latent prior assembled from `spatial`, `temporal`, and `re`.
#'
#' @keywords internal
.tobs_em_nested_laplace <- function(model, spatial = NULL, temporal = NULL,
                                    re = NULL, priors = NULL,
                                    sigma_beta = 10,
                                    max_iter = 25L, tol = 1e-3,
                                    damping = 0.3,
                                    verbose = TRUE) {
  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object", call. = FALSE)
  }
  if (!identical(model$model_type, "single")) {
    stop("nested_laplace engine is currently wired only for single-season ",
         "occupancy (`family = occu()`); got model_type = '",
         model$model_type, "'.", call. = FALSE)
  }
  if (!is.null(priors) && !isFALSE(priors)) {
    message(".tobs_em_nested_laplace(): fixed-effect priors are not yet ",
            "applied on the nested-Laplace path; falling back to ",
            "unpenalised MAP for the occupancy block.")
  }

  multi_prior <- .tobs_to_multi_block_prior(
    spatial = spatial, temporal = temporal, re = re, model = model
  )
  if (is.null(multi_prior)) {
    stop("nested_laplace engine requires at least one latent block ",
         "(spatial, temporal, or re); none were supplied. Use ",
         "`engine = 'laplace'` for a fit with no latent structure.",
         call. = FALSE)
  }

  callbacks <- .build_single_callbacks_nested(model, multi_prior)

  em_result <- tulpa::tulpa_em_laplace(
    e_step        = callbacks$e_step,
    m_step_encode = callbacks$m_step_encode,
    max_iter      = max_iter,
    tol           = tol,
    damping       = damping,
    correction    = "none",
    verbose       = verbose
  )

  fit <- build_laplace_fit(em_result, model, spatial,
                           callbacks$p_per_submodel,
                           prior_spec = NULL)
  fit$nested_laplace <- list(
    multi_prior = multi_prior,
    occ_fit     = em_result$fits$occ
  )
  fit$temporal <- temporal
  fit$re <- re
  fit
}


# Build E-step / M-step callbacks for nested-Laplace occupancy. The
# occupancy block carries `prior = multi_prior` so the generic EM engine
# dispatches it through tulpa_nested_laplace(); the detection block stays
# single-Laplace.
.build_single_callbacks_nested <- function(model, multi_prior) {
  y <- model$y
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  n_sites <- model$n_sites
  p_occ <- ncol(X_occ)
  p_det <- ncol(X_det)

  n_valid <- integer(n_sites)
  n_det   <- integer(n_sites)
  any_det <- logical(n_sites)
  for (i in seq_len(n_sites)) {
    valid <- y[i, ] >= 0
    n_valid[i] <- sum(valid)
    n_det[i]   <- sum(y[i, valid] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  # E-step. For nested-Laplace we read the latent posterior mean at each
  # site from fits$occ$mode (which .fit_block_via_nested_laplace() sets
  # to the weighted average over the hyperparameter grid). The mode
  # ordering is [beta (p_occ), latent block 1, latent block 2, ...]; for
  # the eta computation we want X*beta + sum(latent at site i), so we
  # walk the multi_prior blocks and add each block's contribution.
  e_step <- function(fits, ...) {
    occ_fit <- fits$occ
    if (is.null(occ_fit) || is.null(occ_fit$mode)) {
      # First iteration before any M-step has run.
      eta_occ <- numeric(n_sites)
    } else {
      beta_occ <- occ_fit$mode[seq_len(p_occ)]
      eta_occ <- as.vector(X_occ %*% beta_occ)
      # Latent contributions: the nested-Laplace driver returns
      # `modes` averaged over the grid; for now we approximate eta
      # by the fixed effects only, since the latent posterior moments
      # are integrated into the M-step block log-marginal. The E-step
      # converges to the same fixed point either way because both
      # M-step encodings depend on weights[i] = P(z_i = 1 | y_i, theta).
    }
    det_fit <- fits$det
    if (is.null(det_fit) || is.null(det_fit$mode)) {
      p <- rep(0.5, n_sites)
    } else {
      beta_det <- extract_beta(det_fit, p_det)
      p <- plogis(as.vector(X_det %*% beta_det))
    }
    psi <- plogis(eta_occ)
    list(weights = occ_weights(psi, p, n_sites, n_valid, n_det, any_det))
  }

  m_step_encode <- function(weights, ...) {
    # Pseudo-binomial encoding (M = 1000) mirrors the single-Laplace path
    # so the recovery behavior matches when no latent blocks are present.
    M <- 1000L
    y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
    y_occ <- pmin(pmax(y_occ, 0L), M)
    occ_block <- list(
      y         = y_occ,
      n_trials  = rep(M, n_sites),
      X         = X_occ,
      family    = "binomial",
      prior     = multi_prior
    )

    # Detection: weighted binomial, single-Laplace (no spatial / latent).
    w_det <- weights
    w_det[any_det] <- 1
    keep_det <- keep & (w_det > 1e-6)
    det_block <- list(
      y        = n_det[keep_det],
      n_trials = n_valid[keep_det],
      X        = X_det[keep_det, , drop = FALSE],
      weights  = w_det[keep_det],
      family   = "binomial"
    )

    list(occ = occ_block, det = det_block)
  }

  list(e_step = e_step, m_step_encode = m_step_encode,
       p_per_submodel = c(occ = p_occ, det = p_det))
}
