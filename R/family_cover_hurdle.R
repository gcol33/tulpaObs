# =============================================================================
# family_cover_hurdle.R — Vegetation cover hurdle on the tulpa Laplace backend
#
# Phase 1a: lognormal-positive variant via two independent tulpa_laplace()
# calls (binomial on occur, gaussian on log(cover[occur==1])).
# Phase 1c: joint shared-field model via tulpa_nested_laplace_joint(),
# lognormal-positive on BYM2/ICAR/CAR_proper.
# Phase 1d: beta-positive on the joint engine. phi is profiled (pre-fit on
# the positive subset alone via tulpa_laplace_beta(), then plugged into the
# joint as fixed dispersion). Mirrors the sigma_pos handling for lognormal.
# Full posterior integration over phi is scheduled for Phase 3.
# =============================================================================


# ---------------------------------------------------------------------------
# Dispatcher (called from tobs())
# ---------------------------------------------------------------------------

.dispatch_cover <- function(formula, data, family, detection, y,
                            visit_data, spatial, temporal, engine,
                            priors, control, ...) {
  positive <- family$params$positive
  if (!positive %in% c("lognormal", "beta")) {
    stop("cover(positive = '", positive, "') is not supported. ",
         "Use 'lognormal' or 'beta'.", call. = FALSE)
  }
  if (!is.null(temporal)) {
    stop("`temporal = ` is scheduled for Phase 1d (Mundlak helper). ",
         "Add a year covariate to `formula` for Phase 1a.", call. = FALSE)
  }
  if (!is.null(detection)) {
    stop("`cover()` does not use a detection formula ",
         "(replicates = 'single'). Drop the `detection` argument.",
         call. = FALSE)
  }
  if (is.null(y)) {
    stop("`cover()` requires `y` (a length-N numeric vector of cover ",
         "in [0, 1]).", call. = FALSE)
  }
  enc  <- encode_cover_hurdle(formula, data, y, spatial, positive = positive)

  if (identical(engine, "nested_laplace")) {
    return(decode_cover_hurdle_joint(
      fit_cover_hurdle_joint_nested(enc, data, positive, control), enc, family
    ))
  }

  fits <- fit_cover_hurdle(enc, positive, engine, priors, control)
  decode_cover_hurdle(fits, enc, family)
}


# ---------------------------------------------------------------------------
# Encode
# ---------------------------------------------------------------------------

#' Encode cover-hurdle data for the two-Laplace fit
#'
#' Splits `y` into a binomial occurrence indicator and a positive-cover
#' subset, builds design matrices for each arm using the same `formula`,
#' and packages the spatial spec for both arms.
#'
#' For `positive = "lognormal"` the positive arm's response is
#' `log(y[occur == 1])`. For `positive = "beta"` it is `y[occur == 1]` on
#' the natural (0, 1) scale; an additional eps-clip is applied so the
#' Laplace engine does not see exact 1's introduced upstream.
#'
#' @param formula State-process formula (no LHS); used for both occurrence
#'   and positive-cover arms.
#' @param data Data frame with `nrow(data) == length(y)`.
#' @param y Length-N numeric vector of cover in `[0, 1]`. NAs are dropped
#'   from both arms (treated as missing, not as zero cover).
#' @param spatial Optional `tulpa_spatial` spec, passed to both arms.
#' @param positive One of `"lognormal"` or `"beta"`.
#' @return A list with: `occ_data`, `pos_data`, `spatial_spec`, `N`,
#'   `idx_pos` (row indices of the positive subset within `data`),
#'   `formula`, `positive`.
#' @keywords internal
encode_cover_hurdle <- function(formula, data, y, spatial = NULL,
                                positive = c("lognormal", "beta")) {
  positive <- match.arg(positive)
  if (!is.numeric(y)) stop("`y` must be numeric.", call. = FALSE)
  if (length(y) != nrow(data)) {
    stop(sprintf("length(y) (%d) must equal nrow(data) (%d).",
                 length(y), nrow(data)), call. = FALSE)
  }
  rng <- range(y, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    stop("`y` must be in [0, 1] (got range [", rng[1], ", ", rng[2], "]).",
         call. = FALSE)
  }

  obs_keep <- !is.na(y)
  y_obs    <- y[obs_keep]
  data_obs <- data[obs_keep, , drop = FALSE]
  occur    <- as.integer(y_obs > 0)

  X_occ <- stats::model.matrix(formula, data_obs)

  is_pos <- occur == 1L
  data_pos <- data_obs[is_pos, , drop = FALSE]
  y_pos    <- y_obs[is_pos]
  if (positive == "lognormal") {
    y_pos_resp <- log(y_pos)
  } else {
    # Beta arm needs y strictly in (0, 1). Cap at 1 - 1e-6; lower bound is
    # already guaranteed by occur == 1 + the range check above.
    y_pos_resp <- pmin(y_pos, 1 - 1e-6)
  }
  X_pos <- stats::model.matrix(formula, data_pos)

  list(
    occ_data = list(y = occur, n_trials = rep(1L, length(occur)), X = X_occ),
    pos_data = list(y = y_pos_resp, X = X_pos),
    spatial_spec = spatial,
    N            = length(occur),
    idx_pos      = which(is_pos),
    formula      = formula,
    positive     = positive,
    obs_keep     = obs_keep
  )
}


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

#' Fit the two arms of a cover hurdle
#'
#' Two independent `tulpa::tulpa_laplace()` calls. The joint shared-field
#' fit is Phase 1c. For `positive = "lognormal"` the positive arm is a
#' Gaussian fit on `log(cover)` with sigma estimated post-hoc as the
#' residual standard error. For `positive = "beta"` the positive arm uses
#' `tulpa::tulpa_laplace_beta()` which estimates the precision `phi` via
#' an outer 1-D optimisation and weights the Hessian accordingly.
#'
#' @param enc Output of [encode_cover_hurdle()].
#' @param positive `"lognormal"` or `"beta"` (taken from `enc$positive`).
#' @param engine `"laplace"` (default) or `"nested_laplace"`. The latter is
#'   routed through [fit_cover_hurdle_joint_nested()].
#' @param priors Currently ignored — passed through for forward compat.
#' @param control List with optional `max_iter`, `tol`, `n_threads`.
#' @return List with `m_occ`, `m_pos`, `positive`, `pos_fit_n`, `pos_fit_p`,
#'   plus one of `sigma_pos` (lognormal) or `phi_pos` (beta).
#' @keywords internal
fit_cover_hurdle <- function(enc, positive = enc$positive,
                             engine = "laplace",
                             priors = NULL, control = list()) {
  if (!engine %in% c("laplace", "auto")) {
    stop("cover() currently supports only engine = 'laplace' or ",
         "'nested_laplace' (got '", engine,
         "'). nuts lands in later phases.", call. = FALSE)
  }
  max_iter  <- control$max_iter  %||% 100L
  tol       <- control$tol       %||% 1e-6
  n_threads <- control$n_threads %||% 1L

  m_occ <- tulpa::tulpa_laplace(
    y        = enc$occ_data$y,
    n_trials = enc$occ_data$n_trials,
    X        = enc$occ_data$X,
    family   = "binomial",
    spatial  = enc$spatial_spec,
    max_iter = max_iter, tol = tol, n_threads = n_threads
  )

  if (length(enc$pos_data$y) < ncol(enc$pos_data$X) + 1L) {
    stop("Too few positive-cover sites (", length(enc$pos_data$y),
         ") for the requested formula (", ncol(enc$pos_data$X),
         " coefficients). Need at least ncol(X) + 1.", call. = FALSE)
  }

  n_pos <- length(enc$pos_data$y)
  p_pos <- ncol(enc$pos_data$X)

  if (positive == "lognormal") {
    m_pos <- tulpa::tulpa_laplace(
      y        = enc$pos_data$y,
      n_trials = rep(1L, n_pos),
      X        = enc$pos_data$X,
      family   = "gaussian",
      spatial  = enc$spatial_spec,
      max_iter = max_iter, tol = tol, n_threads = n_threads
    )
    # Gaussian Laplace runs with phi = 1; estimate residual SD post-hoc.
    beta_pos <- m_pos$mode[seq_len(p_pos)]
    eta_pos  <- as.numeric(enc$pos_data$X %*% beta_pos)
    resid    <- enc$pos_data$y - eta_pos
    sigma_pos <- sqrt(sum(resid^2) / max(n_pos - p_pos, 1L))
    return(list(
      m_occ     = m_occ,
      m_pos     = m_pos,
      positive  = "lognormal",
      sigma_pos = sigma_pos,
      pos_fit_n = n_pos,
      pos_fit_p = p_pos
    ))
  }

  # positive == "beta"
  m_pos <- tulpa::tulpa_laplace_beta(
    y         = enc$pos_data$y,
    X         = enc$pos_data$X,
    spatial   = enc$spatial_spec,
    max_iter  = max_iter, tol = tol, n_threads = n_threads
  )
  list(
    m_occ     = m_occ,
    m_pos     = m_pos,
    positive  = "beta",
    phi_pos   = m_pos$phi,
    pos_fit_n = n_pos,
    pos_fit_p = p_pos
  )
}


# ---------------------------------------------------------------------------
# Decode
# ---------------------------------------------------------------------------

#' Decode the two-arm fit into a cover_fit object
#'
#' Extracts beta vectors and SEs for each arm. SEs are scaled to match
#' each arm's dispersion convention:
#'
#' * lognormal arm: `tulpa_laplace(family = "gaussian")` computes the
#'   Hessian assuming phi = 1, so SEs are rescaled by `sigma_pos^2`.
#' * beta arm: `tulpa_laplace_beta()` already weights the Hessian by phi
#'   (Fisher information), so SEs are returned at scale 1.
#'
#' @keywords internal
decode_cover_hurdle <- function(fits, enc, family) {
  beta_occ <- fits$m_occ$mode[seq_len(ncol(enc$occ_data$X))]
  beta_pos <- fits$m_pos$mode[seq_len(ncol(enc$pos_data$X))]
  names(beta_occ) <- colnames(enc$occ_data$X)
  names(beta_pos) <- colnames(enc$pos_data$X)

  se_occ <- .se_from_hessian(fits$m_occ$H_beta, scale = 1)
  se_pos_scale <- if (fits$positive == "lognormal") fits$sigma_pos^2 else 1
  se_pos <- .se_from_hessian(fits$m_pos$H_beta, scale = se_pos_scale)
  if (length(se_occ)) names(se_occ) <- names(beta_occ)
  if (length(se_pos)) names(se_pos) <- names(beta_pos)

  hyperpar <- list(
    occ = .extract_spatial_hyperpar(fits$m_occ, enc$spatial_spec),
    pos = .extract_spatial_hyperpar(fits$m_pos, enc$spatial_spec)
  )
  if (fits$positive == "lognormal") {
    hyperpar$sigma_pos <- fits$sigma_pos
  } else {
    hyperpar$phi_pos <- fits$phi_pos
  }

  out <- structure(
    list(
      occ          = fits$m_occ,
      pos          = fits$m_pos,
      beta_occ     = beta_occ,
      beta_pos     = beta_pos,
      se_occ       = se_occ,
      se_pos       = se_pos,
      positive     = fits$positive,
      sigma_pos    = if (fits$positive == "lognormal") fits$sigma_pos else NA_real_,
      phi_pos      = if (fits$positive == "beta")      fits$phi_pos    else NA_real_,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = length(enc$idx_pos),
      converged    = isTRUE(fits$m_occ$converged) && isTRUE(fits$m_pos$converged),
      log_marginal = c(occ = fits$m_occ$log_marginal,
                       pos = fits$m_pos$log_marginal)
    ),
    class = c("cover_fit", "tobs_fit", "tulpa_fit")
  )
  out
}


# ---------------------------------------------------------------------------
# Predict
# ---------------------------------------------------------------------------

#' Predict cover from a cover_fit
#'
#' Occurrence probability is always `p = plogis(X * beta_occ)`. The
#' conditional positive cover `mu` depends on the positive-part family:
#'
#' * `positive = "lognormal"`: `mu = exp(eta_pos + sigma_pos^2 / 2)`
#'   (lognormal back-transform on log-cover).
#' * `positive = "beta"`: `mu = plogis(eta_pos)` (mean of the beta on
#'   the natural cover scale with logit link).
#'
#' Expected cover is `E[y] = p * mu` under both positive parts.
#'
#' Spatial random effects are not yet projected at new locations
#' (Phase 1a is fixed-effects-only for prediction; spatial projection
#' lands in 1c).
#'
#' @param object A `cover_fit`.
#' @param newdata A data frame of covariates matching the original formula.
#' @param type One of `"expected"`, `"occupancy"`, `"conditional"`.
#' @param include_RE Currently ignored (no spatial projection in 1a).
#' @param ... Unused.
#' @return Numeric vector of predictions.
#' @export
predict.cover_fit <- function(object, newdata,
                                     type = c("expected", "occupancy",
                                              "conditional"),
                                     include_RE = FALSE, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` is required.", call. = FALSE)
  }
  if (isTRUE(include_RE) && !is.null(object$encoding$spatial_spec)) {
    message("predict.cover_fit(): spatial RE projection at new ",
            "locations is not implemented in Phase 1a; returning ",
            "fixed-effects-only predictions.")
  }

  X <- stats::model.matrix(object$encoding$formula, newdata)
  if (ncol(X) != length(object$beta_occ)) {
    stop("Design-matrix column count (", ncol(X), ") does not match the ",
         "fitted model (", length(object$beta_occ), "). Check `newdata`.",
         call. = FALSE)
  }

  eta_occ <- as.numeric(X %*% object$beta_occ)
  eta_pos <- as.numeric(X %*% object$beta_pos)
  p  <- stats::plogis(eta_occ)
  positive <- object$positive %||% "lognormal"
  mu <- if (positive == "beta") {
    stats::plogis(eta_pos)
  } else {
    exp(eta_pos + object$sigma_pos^2 / 2)
  }

  switch(
    type,
    expected    = p * mu,
    occupancy   = p,
    conditional = mu
  )
}


# ---------------------------------------------------------------------------
# Print / summary
# ---------------------------------------------------------------------------

#' @export
print.cover_fit <- function(x, ...) {
  positive <- x$positive %||% "lognormal"
  cat(sprintf("<cover_fit (%s positive part)>\n", positive))
  cat(sprintf("  N total      : %d\n", x$n_total))
  cat(sprintf("  N positive   : %d (%.1f%%)\n",
              x$n_positive, 100 * x$n_positive / x$n_total))
  if (positive == "lognormal") {
    cat(sprintf("  sigma_pos    : %.4f\n", x$sigma_pos))
  } else {
    cat(sprintf("  phi_pos      : %.4f\n", x$phi_pos))
  }
  cat(sprintf("  converged    : %s\n",
              if (isTRUE(x$converged)) "yes" else "no"))
  cat("\nOccurrence (binomial logit):\n")
  print(.coef_table(x$beta_occ, x$se_occ))
  pos_header <- if (positive == "beta") {
    "Cover (beta, logit link, on y > 0):"
  } else {
    "Log-cover (Gaussian on log y > 0):"
  }
  cat("\n", pos_header, "\n", sep = "")
  print(.coef_table(x$beta_pos, x$se_pos))
  invisible(x)
}

#' @export
summary.cover_fit <- function(object, ...) {
  out <- list(
    family       = object$family,
    positive     = object$positive %||% "lognormal",
    n_total      = object$n_total,
    n_positive   = object$n_positive,
    sigma_pos    = object$sigma_pos,
    phi_pos      = object$phi_pos,
    converged    = object$converged,
    occurrence   = .coef_table(object$beta_occ, object$se_occ),
    positive_arm = .coef_table(object$beta_pos, object$se_pos),
    log_marginal = object$log_marginal,
    hyperpar     = object$hyperpar
  )
  class(out) <- "summary.cover_fit"
  out
}

#' @export
print.summary.cover_fit <- function(x, ...) {
  cat("Cover hurdle fit summary\n")
  cat(sprintf("  positive part: %s\n", x$positive))
  cat(sprintf("  N total = %d, N positive = %d\n", x$n_total, x$n_positive))
  if (x$positive == "lognormal") {
    cat(sprintf("  sigma_pos = %.4f\n", x$sigma_pos))
  } else {
    cat(sprintf("  phi_pos   = %.4f\n", x$phi_pos))
  }
  cat(sprintf("  log marginal: occ = %.3f, pos = %.3f\n",
              x$log_marginal["occ"], x$log_marginal["pos"]))
  cat("\nOccurrence:\n"); print(x$occurrence)
  pos_header <- if (x$positive == "beta") "Cover (beta, logit):" else "Log-cover (Gaussian):"
  cat("\n", pos_header, "\n", sep = ""); print(x$positive_arm)
  invisible(x)
}


# ---------------------------------------------------------------------------
# Joint nested-Laplace fit (Phase 1c lognormal, Phase 1d beta)
# ---------------------------------------------------------------------------

#' Fit cover_hurdle as a joint binomial+(gaussian|beta) model with shared
#' spatial field via [tulpa::tulpa_nested_laplace_joint()].
#'
#' For `positive = "lognormal"` the positive arm runs as Gaussian on
#' `log(cover)` with `sigma_pos` estimated post-hoc as the residual standard
#' error at the joint mode (matches Phase 1a).
#'
#' For `positive = "beta"` the positive arm runs as Beta on `cover` (logit
#' link). The precision `phi` is **profiled**: a pre-fit via
#' [tulpa::tulpa_laplace_beta()] on the positive subset (no spatial) gives
#' `phi_hat`, which is then held fixed in the joint engine while the
#' spatial hyperparameters are integrated out. Full posterior integration
#' over `phi` is scheduled for Phase 3.
#'
#' @param enc Output of [encode_cover_hurdle()].
#' @param data The original (un-subsetted) data frame — required to resolve
#'   the spatial spec (group_var lookup, n_spatial_units check).
#' @param positive `"lognormal"` or `"beta"`.
#' @param control List with optional `max_iter`, `tol`, `n_threads`,
#'   `sigma_grid`, `rho_grid`, `tau_grid`, `rho_car_grid`, `alpha_grid`,
#'   `phi_init`, `phi_bounds` (the last two are forwarded to the beta
#'   pre-fit when `positive = "beta"`).
#' @return List shaped like the single-Laplace fit output but with extra
#'   `joint` field carrying the raw `tulpa_nested_laplace_joint` result.
#' @keywords internal
fit_cover_hurdle_joint_nested <- function(enc, data, positive = enc$positive,
                                          control = list()) {
  if (!positive %in% c("lognormal", "beta")) {
    stop("engine = 'nested_laplace' for cover() supports positive = ",
         "'lognormal' or 'beta'. Got '", positive, "'.", call. = FALSE)
  }
  if (is.null(enc$spatial_spec)) {
    stop("engine = 'nested_laplace' for cover() requires a spatial spec. ",
         "Pass one of `tulpa::spatial_bym2(adj)`, `tulpa::spatial_icar(adj)`, ",
         "or `tulpa::spatial_car_proper(adj)`.", call. = FALSE)
  }
  spec <- enc$spatial_spec
  if (!inherits(spec, "tulpa_spatial")) {
    stop("engine = 'nested_laplace' for cover(): `spatial` must be a ",
         "tulpa_spatial spec.", call. = FALSE)
  }
  spec_type <- tolower(spec$type)
  # tulpa::spatial_car() returns type = "car" but prior_from_spec maps it
  # to backend = "icar"; treat the two as equivalent at dispatch time.
  supported <- c("bym2", "icar", "car", "car_proper")
  if (!spec_type %in% supported) {
    stop("engine = 'nested_laplace' for cover() supports spatial types: ",
         paste(shQuote(supported), collapse = ", "),
         ". Got type = '", spec$type, "'.", call. = FALSE)
  }

  # Resolve obs -> spatial unit via tulpa's prior_from_spec. The dropped-NA
  # rows in encode_cover_hurdle (obs_keep) shrink the obs set; subset the
  # spatial_idx vector accordingly.
  data_obs <- data[enc$obs_keep, , drop = FALSE]
  prior    <- tulpa::prior_from_spec(spec, data_obs)
  spi_full <- prior$spatial_idx                 # length N (post-NA-drop)
  spi_pos  <- spi_full[enc$idx_pos]             # length N_pos

  N     <- enc$N
  N_pos <- length(enc$pos_data$y)

  arm_occ <- list(
    y           = as.numeric(enc$occ_data$y),
    n_trials    = enc$occ_data$n_trials,
    X           = enc$occ_data$X,
    spatial_idx = as.integer(spi_full),
    re_idx      = rep(0, N),
    n_re_groups = 0L,
    sigma_re    = 1.0,
    family      = "binomial",
    phi         = 1.0
  )

  # Positive-arm dispersion: gaussian phi is the noise SD — pre-fit it on
  # the positive subset, otherwise the joint integrand sees noise scale 1
  # regardless of truth, and the marginal likelihood across the alpha grid
  # is unable to discriminate alpha because the alpha-driven field variance
  # is observationally indistinguishable from residual noise of size 1.
  # (issue #4). Beta needs a real phi for the same Newton-scaling reason.
  if (positive == "lognormal") {
    pos_family <- "gaussian"
    phi_hat    <- .prefit_lognormal_sigma(enc, control)
  } else {
    phi_hat <- .prefit_beta_phi(enc, control)
    pos_family <- "beta"
  }

  arm_pos <- list(
    y           = as.numeric(enc$pos_data$y),
    n_trials    = rep(1L, N_pos),
    X           = enc$pos_data$X,
    spatial_idx = as.integer(spi_pos),
    re_idx      = rep(0, N_pos),
    n_re_groups = 0L,
    sigma_re    = 1.0,
    family      = pos_family,
    phi         = phi_hat
  )

  # Strip the per-obs spatial_idx (tulpa_nested_laplace_joint takes it per
  # arm) and the legacy rho_bounds field (joint car_proper uses rho_car_grid).
  # Forward control-grid overrides per backend.
  prior_for_joint <- prior
  prior_for_joint$spatial_idx <- NULL
  prior_for_joint$rho_bounds  <- NULL
  if (!is.null(control$sigma_grid))   prior_for_joint$sigma_grid   <- control$sigma_grid
  if (!is.null(control$rho_grid))     prior_for_joint$rho_grid     <- control$rho_grid
  if (!is.null(control$tau_grid))     prior_for_joint$tau_grid     <- control$tau_grid
  if (!is.null(control$rho_car_grid)) prior_for_joint$rho_car_grid <- control$rho_car_grid

  copy_spec <- list(
    arm        = "pos",
    alpha_grid = control$alpha_grid %||% c(0.0, 0.5, 1.0, 1.5, 2.0)
  )

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = list(occ = arm_occ, pos = arm_pos),
    prior     = prior_for_joint,
    copy      = copy_spec,
    max_iter  = control$max_iter  %||% 50L,
    tol       = control$tol       %||% 1e-6,
    n_threads = control$n_threads %||% 1L,
    store_Q   = TRUE
  )

  # Posterior-weighted mean / SE for the per-arm beta blocks.
  layout <- fit$arm_layout
  p_occ  <- layout$p[1]
  p_pos  <- layout$p[2]
  bocc_idx <- layout$beta_start[1] + seq_len(p_occ)
  bpos_idx <- layout$beta_start[2] + seq_len(p_pos)

  beta_occ <- as.numeric(crossprod(fit$weights, fit$modes[, bocc_idx, drop = FALSE]))
  beta_pos <- as.numeric(crossprod(fit$weights, fit$modes[, bpos_idx, drop = FALSE]))

  # Total posterior variance per beta = Var-of-means (across grid) +
  # Mean-of-Var (per-grid inner Laplace variance, posterior-weighted). The
  # joint-Laplace kernel returns the per-grid joint precision Q_k (lower
  # triangle, CSC) when `store_Q = TRUE`; per-grid inner variance for the
  # beta sub-block is diag(Q_k^{-1})[beta_idx], computed via sparse Cholesky.
  # See gcol33/tulpaObs#2 for the ablation that motivated this.
  var_of_means_occ <- as.numeric(crossprod(fit$weights,
    fit$modes[, bocc_idx, drop = FALSE]^2)) - beta_occ^2
  var_of_means_pos <- as.numeric(crossprod(fit$weights,
    fit$modes[, bpos_idx, drop = FALSE]^2)) - beta_pos^2

  beta_idx_all <- c(bocc_idx, bpos_idx)
  inner_var <- .joint_inner_var(fit, beta_idx_all)
  if (is.null(inner_var)) {
    # No Q stored (older tulpa or unsupported backend); fall back to the
    # var-of-means-only SE that the validation harness flagged as too narrow.
    mean_of_var_occ <- rep(0, p_occ)
    mean_of_var_pos <- rep(0, p_pos)
  } else {
    occ_cols <- seq_along(bocc_idx)
    pos_cols <- length(bocc_idx) + seq_along(bpos_idx)
    mean_of_var_occ <- as.numeric(crossprod(fit$weights,
      inner_var[, occ_cols, drop = FALSE]))
    mean_of_var_pos <- as.numeric(crossprod(fit$weights,
      inner_var[, pos_cols, drop = FALSE]))
  }

  se_occ <- sqrt(pmax(0, var_of_means_occ + mean_of_var_occ))
  se_pos <- sqrt(pmax(0, var_of_means_pos + mean_of_var_pos))

  # Dispersion summary on the positive arm at the posterior-weighted-mean
  # beta_pos. For both families we subtract the alpha-scaled posterior
  # spatial field at the positive-arm locations before fitting the
  # dispersion — without that subtraction the residual / over-dispersion
  # signal absorbs alpha^2 * field_sigma^2 (lognormal -> sigma_pos in #4,
  # beta -> phi_pos in #5).
  eta_pos      <- as.numeric(enc$pos_data$X %*% beta_pos)
  field_at_pos <- .joint_field_at_obs_copy(fit, prior_for_joint, spi_pos)
  if (positive == "lognormal") {
    resid     <- enc$pos_data$y - eta_pos - field_at_pos
    sigma_pos <- sqrt(sum(resid^2) / max(N_pos - p_pos, 1L))
    phi_pos   <- NA_real_
  } else {
    sigma_pos <- NA_real_
    phi_pos   <- .refit_beta_phi_postfield(
      y      = enc$pos_data$y,
      eta    = eta_pos + field_at_pos,
      bounds = control$phi_bounds %||% c(0.1, 1e4)
    )
  }

  m_occ <- list(mode = beta_occ, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  m_pos <- list(mode = beta_pos, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)

  list(
    m_occ      = m_occ,
    m_pos      = m_pos,
    positive   = positive,
    sigma_pos  = sigma_pos,
    phi_pos    = phi_pos,
    pos_fit_n  = N_pos,
    pos_fit_p  = p_pos,
    beta_occ   = beta_occ,
    beta_pos   = beta_pos,
    se_occ     = se_occ,
    se_pos     = se_pos,
    joint      = fit
  )
}

# Pre-fit the lognormal residual SD on the positive subset before handing
# control to the joint engine. The joint integrand reads `phi` as the noise
# SD; without a sensible pre-fit it sees scale 1 regardless of truth and the
# log-marginal across the alpha grid becomes near-flat (issue #4).
#
# Strategy: non-spatial Gaussian fit on the positive subset, residual SD as
# the point estimate. This is an upper bound on the true noise SD (it
# includes the alpha-mediated field variance), but it sits inside the same
# order of magnitude as the truth, which is enough to restore the joint
# engine's discrimination across the alpha grid. The post-hoc sigma_pos in
# `fit_cover_hurdle_joint_nested` then refines this by subtracting the
# alpha-scaled posterior field.
.prefit_lognormal_sigma <- function(enc, control) {
  y <- enc$pos_data$y
  X <- enc$pos_data$X
  n <- length(y); p <- ncol(X)
  if (n <= p) return(1.0)
  beta_init  <- tryCatch(qr.solve(X, y), error = function(e) NULL)
  if (is.null(beta_init)) return(1.0)
  resid_init <- as.numeric(y - X %*% beta_init)
  sigma_init <- sqrt(sum(resid_init^2) / max(n - p, 1L))
  if (!is.finite(sigma_init) || sigma_init <= 0) return(1.0)
  sigma_init
}

# Project the posterior-weighted alpha-scaled spatial field at a set of
# spatial-unit indices (1-based) under a `tulpa_nested_laplace_joint` fit
# with a copy arm. Backend-aware: BYM2 combines (phi, theta) sub-blocks
# with the per-grid sigma/rho factors; ICAR / CAR_proper use the phi block
# unscaled. Returns a numeric vector of the same length as `spi_obs` —
# zero when `fit` has no modes, no copy spec, or theta_grid lacks an alpha
# column.
.joint_field_at_obs_copy <- function(fit, prior, spi_obs) {
  N <- length(spi_obs)
  if (N == 0L) return(numeric(0))
  modes <- fit$modes
  if (is.null(modes) || !is.matrix(modes)) return(rep(0, N))
  layout <- fit$arm_layout
  theta_grid <- fit$theta_grid
  weights <- fit$weights
  if (is.null(theta_grid) || !"alpha" %in% colnames(theta_grid)) {
    return(rep(0, N))
  }
  alpha_k <- theta_grid[, "alpha"]
  n_grid <- length(alpha_k)
  n_s <- prior$n_spatial_units %||% 0L
  if (n_s == 0L) return(rep(0, N))
  phi_start <- layout$phi_start
  theta_start <- layout$theta_start
  if (is.null(phi_start) || n_s <= 0L) return(rep(0, N))

  type <- tolower(prior$type %||% "")
  scale_factor <- as.numeric(prior$scale_factor %||% 1.0)
  field_at_obs <- numeric(N)
  for (k in seq_len(n_grid)) {
    if (type == "bym2") {
      sigma_k <- theta_grid[k, "sigma"]
      rho_k   <- theta_grid[k, "rho"]
      d_phi   <- sigma_k * sqrt(rho_k + 1e-10) * scale_factor
      d_theta <- sigma_k * sqrt(1 - rho_k + 1e-10)
      phi_k   <- modes[k, phi_start + seq_len(n_s)]
      theta_k <- modes[k, theta_start + seq_len(n_s)]
      field_k <- d_phi * phi_k + d_theta * theta_k
    } else {
      # ICAR / CAR_proper: latent x[s] is the field directly (d_phi = 1).
      phi_k   <- modes[k, phi_start + seq_len(n_s)]
      field_k <- phi_k
    }
    field_at_obs <- field_at_obs + weights[k] * alpha_k[k] * field_k[spi_obs]
  }
  field_at_obs
}

# Post-hoc refit of the beta dispersion phi conditional on the
# posterior-weighted linear predictor (including the alpha-scaled spatial
# field). The pre-fit in .prefit_beta_phi() runs without the spatial
# component, so when alpha is non-zero the field variance leaks into phi
# as apparent over-dispersion and phi_pos comes out downward-biased
# (issue #5; mirrors the lognormal correction in #4).
#
# Beta(mu*phi, (1-mu)*phi) has no closed-form MLE for phi given mu, so we
# profile the negative log-likelihood by 1D Brent search:
#
#   -loglik(phi | y, mu) =
#       - sum_i [ lgamma(phi) - lgamma(mu_i*phi) - lgamma((1-mu_i)*phi)
#                 + (mu_i*phi - 1) log(y_i)
#                 + ((1-mu_i)*phi - 1) log(1 - y_i) ]
#
# mu is clipped away from {0, 1} so lgamma(0) doesn't fire; y likewise (the
# encoder already guards against {0, 1} but the post-hoc mu can land near
# the boundary even when y is interior).
.refit_beta_phi_postfield <- function(y, eta, bounds = c(0.1, 1e4)) {
  if (length(y) < 2L) return(NA_real_)
  eps <- 1e-6
  mu  <- pmin(pmax(plogis(eta), eps), 1 - eps)
  y_c <- pmin(pmax(y,           eps), 1 - eps)
  nloglik <- function(phi) {
    a <- mu * phi
    b <- (1 - mu) * phi
    -sum(lgamma(phi) - lgamma(a) - lgamma(b) +
         (a - 1) * log(y_c) + (b - 1) * log(1 - y_c))
  }
  res <- tryCatch(
    stats::optimize(nloglik, interval = bounds, tol = 1e-4),
    error = function(e) NULL
  )
  if (is.null(res) || !is.finite(res$minimum)) return(NA_real_)
  res$minimum
}

# Pre-fit the beta precision on the positive subset (no spatial). Profiled
# treatment of phi: held fixed at this value while the joint engine integrates
# over the spatial hyperparameters. Cheap, runs in well under a second on
# Phase 1d test sizes.
.prefit_beta_phi <- function(enc, control) {
  fit <- tulpa::tulpa_laplace_beta(
    y          = enc$pos_data$y,
    X          = enc$pos_data$X,
    spatial    = NULL,
    max_iter   = control$max_iter   %||% 100L,
    tol        = control$tol        %||% 1e-6,
    n_threads  = control$n_threads  %||% 1L,
    phi_init   = control$phi_init,
    phi_bounds = control$phi_bounds %||% c(0.1, 1e4)
  )
  phi <- fit$phi
  if (!is.finite(phi) || phi <= 0) {
    stop("Pre-fit of beta precision returned a non-positive phi (",
         phi, "). Pass `control$phi_init` or `control$phi_bounds`.",
         call. = FALSE)
  }
  phi
}

#' Decode a joint-nested-Laplace cover-hurdle fit into a `cover_fit`.
#'
#' Lighter-weight than the single-Laplace decode: the joint engine has
#' already produced posterior moments for beta and the spatial hyperparameters,
#' so we just shape them into the existing `cover_fit` structure.
#'
#' @keywords internal
decode_cover_hurdle_joint <- function(fits, enc, family) {
  beta_occ <- fits$beta_occ
  beta_pos <- fits$beta_pos
  names(beta_occ) <- colnames(enc$occ_data$X)
  names(beta_pos) <- colnames(enc$pos_data$X)

  se_occ <- fits$se_occ
  se_pos <- fits$se_pos
  if (length(se_occ)) names(se_occ) <- names(beta_occ)
  if (length(se_pos)) names(se_pos) <- names(beta_pos)

  hyperpar <- list(
    spatial = fits$joint$theta_mean,
    engine  = "nested_laplace"
  )
  if (fits$positive == "lognormal") {
    hyperpar$sigma_pos <- fits$sigma_pos
  } else {
    hyperpar$phi_pos <- fits$phi_pos
  }

  out <- structure(
    list(
      occ          = fits$m_occ,
      pos          = fits$m_pos,
      beta_occ     = beta_occ,
      beta_pos     = beta_pos,
      se_occ       = se_occ,
      se_pos       = se_pos,
      positive     = fits$positive,
      sigma_pos    = fits$sigma_pos,
      phi_pos      = fits$phi_pos,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = length(enc$idx_pos),
      converged    = TRUE,
      log_marginal = c(joint = max(fits$joint$log_marginal)),
      joint        = fits$joint
    ),
    class = c("cover_fit", "tobs_fit", "tulpa_fit")
  )
  out
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Per-grid inner posterior variance for selected latent coordinates,
# applying sum-to-zero constraints on the BYM2/ICAR/CAR_proper spatial
# blocks (phi, theta) so the fixed-effect intercept is data-identified
# rather than prior-bounded.
#
# The joint-Laplace precision Q_k is near-singular along the
# (intercept, mean(phi)) direction whenever the prior on the spatial
# block has a sum-to-zero soft-null (ICAR is rank-deficient, BYM2's
# phi block likewise). The unconstrained inverse maps that direction
# onto the weak beta prior (1e-4 in the kernel = sd 100), producing
# meaningless intercept SEs. The fix is the standard INLA constraint
# correction:
#
#   Sigma_c = Q^{-1} - Q^{-1} A^T (A Q^{-1} A^T)^{-1} A Q^{-1}
#
# where A picks the per-block sums of phi (and theta, for BYM2). With
# A = 0 (no spatial block) this reduces to the unconstrained inverse,
# which is the right behaviour for the SPDE case (Q non-singular).
#
# Returns an `n_grid x length(beta_idx)` matrix of constrained
# Var(beta_j | data, theta_k), or NULL when no Q matrices were stored.
.joint_inner_var <- function(fit, beta_idx) {
  Qp <- fit$Q_csc_p_per_grid
  Qi <- fit$Q_csc_i_per_grid
  Qx <- fit$Q_csc_x_per_grid
  n_x <- fit$Q_csc_n
  if (is.null(Qp) || is.null(Qi) || is.null(Qx) || is.null(n_x)) return(NULL)

  layout <- fit$arm_layout
  n_s <- layout$n_x - max(layout$phi_start %||% layout$n_x,
                          layout$theta_start %||% layout$n_x)
  # Build constraint matrix A (k x n_x): one row of all-ones on each
  # structured spatial block. layout offsets are 0-based.
  A_cols <- list()
  if (!is.null(layout$phi_start)) {
    n_s_phi <- (layout$theta_start %||% layout$n_x) - layout$phi_start
    A_cols[[length(A_cols) + 1L]] <- layout$phi_start + seq_len(n_s_phi)
  }
  if (!is.null(layout$theta_start)) {
    n_s_theta <- layout$n_x - layout$theta_start
    A_cols[[length(A_cols) + 1L]] <- layout$theta_start + seq_len(n_s_theta)
  }
  k_constr <- length(A_cols)

  n_grid <- length(Qp)
  p <- length(beta_idx)
  out <- matrix(NA_real_, n_grid, p)

  E <- Matrix::sparseMatrix(
    i = beta_idx, j = seq_len(p), x = 1,
    dims = c(n_x, p)
  )
  A_t <- if (k_constr > 0L) {
    ii <- unlist(A_cols)
    jj <- rep(seq_len(k_constr), vapply(A_cols, length, integer(1)))
    Matrix::sparseMatrix(i = ii, j = jj, x = 1,
                         dims = c(n_x, k_constr))
  } else NULL

  for (k in seq_len(n_grid)) {
    if (is.null(Qp[[k]]) || length(Qx[[k]]) == 0L) next
    Qk_lt <- Matrix::sparseMatrix(
      i = as.integer(Qi[[k]]) + 1L,
      p = as.integer(Qp[[k]]),
      x = as.numeric(Qx[[k]]),
      dims = c(n_x, n_x),
      symmetric = FALSE,
      index1 = TRUE
    )
    Qk <- Matrix::forceSymmetric(Qk_lt, uplo = "L")
    V <- tryCatch(Matrix::solve(Qk, E), error = function(e) NULL)
    if (is.null(V)) next
    var_uncon <- vapply(seq_len(p),
      function(j) as.numeric(V[beta_idx[j], j]), numeric(1))

    if (k_constr > 0L) {
      W <- tryCatch(Matrix::solve(Qk, A_t), error = function(e) NULL)
      if (!is.null(W)) {
        AV <- as.matrix(Matrix::crossprod(A_t, V))     # k_constr x p
        M  <- as.matrix(Matrix::crossprod(A_t, W))     # k_constr x k_constr
        corr <- vapply(seq_len(p), function(j) {
          v <- AV[, j]
          as.numeric(crossprod(v, solve(M, v)))
        }, numeric(1))
        out[k, ] <- pmax(var_uncon - corr, 0)
        next
      }
    }
    out[k, ] <- pmax(var_uncon, 0)
  }
  out
}

.se_from_hessian <- function(H, scale = 1) {
  if (is.null(H)) return(numeric(0))
  cov <- tryCatch(scale * solve(H), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, nrow(H)))
  sqrt(pmax(diag(cov), 0))
}

.coef_table <- function(beta, se) {
  if (length(se) != length(beta)) se <- rep(NA_real_, length(beta))
  z <- beta / se
  data.frame(
    estimate = beta,
    std.err  = se,
    z.value  = z,
    row.names = names(beta)
  )
}

.extract_spatial_hyperpar <- function(fit, spec) {
  if (is.null(spec)) return(NULL)
  out <- list()
  for (nm in c("range", "sigma", "sigma2_gp", "phi_gp", "tau_spatial",
               "sigma_spatial", "rho")) {
    if (!is.null(fit[[nm]])) out[[nm]] <- fit[[nm]]
  }
  if (length(out) == 0) NULL else out
}
