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
                            visit_data, spatial, temporal, re, engine,
                            priors, control,
                            approx = "gaussian_laplace", ...) {
  positive <- family$params$positive
  if (!positive %in% c("lognormal", "beta")) {
    stop("cover(positive = '", positive, "') is not supported. ",
         "Use 'lognormal' or 'beta'.", call. = FALSE)
  }
  has_multi <- !is.null(temporal) || (!is.null(re) && length(re) > 0L)
  if (has_multi && !identical(engine, "nested_laplace")) {
    stop("`temporal = ` and `re = ` for cover() require ",
         "engine = 'nested_laplace' (got '", engine, "'). ",
         "The single-Laplace path is fixed-effects + spatial only.",
         call. = FALSE)
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
      fit_cover_hurdle_joint_nested(enc, data, positive, control,
                                    temporal = temporal, re = re),
      enc, family, approx = approx
    ))
  }

  fits <- fit_cover_hurdle(enc, positive, engine, priors, control)
  decode_cover_hurdle(fits, enc, family, approx = approx)
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
#' When `approx = "simplified_laplace"`, the cover-hurdle SLA gamma is
#' computed via [`.sla_compute_cover_hurdle()`]: a per-arm 5-point FD of
#' the *original* Bernoulli / Beta / Lognormal log-likelihood against the
#' arm's `solve(H_beta)` Sigma (raw Hessian — no Louis correction needed
#' here because both arms run real likelihoods at the mode, not the
#' pseudo-binomial M-step encoding). Per-arm pseudo-draws are then
#' resampled from skew-normals fit by moment-matching `(beta_arm,
#' se_arm, gamma_arm)`.
#'
#' @keywords internal
decode_cover_hurdle <- function(fits, enc, family,
                                approx = "gaussian_laplace") {
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

  # Simplified-Laplace gamma + skew-normal pseudo-draws per arm.
  skew_occ <- NULL
  skew_pos <- NULL
  draws_occ <- NULL
  draws_pos <- NULL
  sla_status <- "off"
  if (identical(approx, "simplified_laplace")) {
    sla_res <- .sla_compute_cover_hurdle(fits, enc, fits$positive)
    sla_draws <- .sla_build_cover_hurdle_draws(
      beta_occ, se_occ, beta_pos, se_pos, sla_res
    )
    draws_occ <- sla_draws$draws_occ
    draws_pos <- sla_draws$draws_pos
    sla_status <- sla_draws$sla_status
    if (isTRUE(sla_res$valid)) {
      skew_occ <- sla_res$gamma_occ
      skew_pos <- sla_res$gamma_pos
    }
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
      sigma_pos_sd = NA_real_,
      phi_pos      = if (fits$positive == "beta")      fits$phi_pos    else NA_real_,
      phi_pos_sd   = NA_real_,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = length(enc$idx_pos),
      converged    = isTRUE(fits$m_occ$converged) && isTRUE(fits$m_pos$converged),
      log_marginal = c(occ = fits$m_occ$log_marginal,
                       pos = fits$m_pos$log_marginal),
      skew_occ     = skew_occ,
      skew_pos     = skew_pos,
      draws_occ    = draws_occ,
      draws_pos    = draws_pos,
      sla_status   = sla_status
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
  if (!is.null(x$sla_status) && !identical(x$sla_status, "off")) {
    cat(sprintf("  marginals    : %s\n", x$sla_status))
  }
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
#' For both positive parts the dispersion scalar is integrated on the outer
#' joint hyperparameter grid (per-arm `phi_pos` axis):
#'
#' * `positive = "lognormal"` (gaussian arm): the residual SD is the per-grid
#'   phi. The default 7-point log-spaced grid is centred on the non-spatial
#'   prefit from [.prefit_lognormal_sigma()] and spans `[sigma_hat / 3,
#'   sigma_hat * 3]`. The posterior mean and SD across that axis are
#'   surfaced as `sigma_pos` / `sigma_pos_sd` on the returned `cover_fit`.
#' * `positive = "beta"`: the beta precision is the per-grid phi. The
#'   default 7-point log-spaced grid spans `[2, 300]`; posterior mean and
#'   SD are surfaced as `phi_pos` / `phi_pos_sd`.
#'
#' Override the per-arm phi grid via `control$phi_grid`.
#'
#' @param enc Output of [encode_cover_hurdle()].
#' @param data The original (un-subsetted) data frame — required to resolve
#'   the spatial spec (group_var lookup, n_spatial_units check).
#' @param positive `"lognormal"` or `"beta"`.
#' @param control List with optional `max_iter`, `tol`, `n_threads`,
#'   `sigma_grid`, `rho_grid`, `rho_car_grid`, `sigma_pos_grid`,
#'   `phi_init`, `phi_bounds` (the last two are forwarded to the beta
#'   pre-fit when `positive = "beta"`). For ICAR / CAR_proper backends
#'   `tau_grid` is also accepted and translated to `sigma_grid` as
#'   `sigma = 1 / sqrt(tau)` — see gcol33/tulpa#18 for the rationale
#'   behind moving the cover-arm field amplitude onto its own
#'   `sigma_pos_grid` axis. Regularizing hyperpriors on the spatial
#'   field amplitudes can be set via `prior_sigma_occ` / `prior_sigma_pos`
#'   (each a length-2 list `list(family, params)` matching tulpa's
#'   `prior_sigma_*` args — see gcol33/tulpa#22 for the donor-vs-skewed
#'   small-n_pos motivation).
#' @return List shaped like the single-Laplace fit output but with extra
#'   `joint` field carrying the raw `tulpa_nested_laplace_joint` result.
#' @keywords internal
fit_cover_hurdle_joint_nested <- function(enc, data, positive = enc$positive,
                                          control = list(),
                                          temporal = NULL, re = NULL) {
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

  has_multi <- !is.null(temporal) || (!is.null(re) && length(re) > 0L)

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

  # Positive-arm dispersion. Both regimes integrate phi on the outer joint
  # hyperparameter grid; `arm_pos$phi` is a placeholder overridden per grid
  # point by the joint engine.
  #   * lognormal: gaussian phi is the noise SD. The non-spatial residual SD
  #     from `.prefit_lognormal_sigma()` is an upper bound on the truth (it
  #     absorbs the alpha-mediated field variance), so we use it as the
  #     *centre* of a 7-point log-spaced grid spanning [sigma_hat/3,
  #     sigma_hat*3]. Neighbour-ratio ~1.44 keeps the inner Laplace
  #     warm-starts close enough that adaptive densification rarely fires.
  #     Override via `control$phi_grid`.
  #   * beta:      phi is integrated on the outer joint hyperparameter grid
  #     (tulpaObs#7). 7 log-spaced points span 2..300 (neighbour-ratio
  #     ~2.4); the joint engine's mode-tracked interior densification
  #     (gcol33/tulpa#19 follow-up) adds 1-2 midpoint cells around the peak
  #     when adjacent grid levels carry density above the edge threshold,
  #     so the *effective* phi resolution near the peak matches the
  #     previous 13-point default while the baseline cell count drops ~46%.
  #     History: fixed 13 was set when adaptive_grid was off-by-default;
  #     5 points gave ~18% mean bias, 9 points ~12% under the static grid.
  #     Refinement now closes that gap dynamically.
  if (positive == "lognormal") {
    pos_family   <- "gaussian"
    sigma_hat    <- .prefit_lognormal_sigma(enc, control)
    phi_hat      <- sigma_hat
    phi_grid_pos <- control$phi_grid %||%
      exp(seq(log(sigma_hat / 3), log(sigma_hat * 3), length.out = 7))
  } else {
    pos_family   <- "beta"
    phi_hat      <- 1.0
    phi_grid_pos <- control$phi_grid %||%
      exp(seq(log(2), log(300), length.out = 7))
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
  #
  # gcol33/tulpa#18: the engine now parameterizes the joint outer grid as
  # (sigma_occ, sigma_pos) instead of (sigma, alpha). For ICAR / CAR_proper
  # the spatial field is unit-precision and the donor-arm amplitude lives
  # on `prior$sigma_grid` (legacy `tau_grid` from prior_from_spec is
  # translated below). The cover-arm field amplitude lives on
  # `copy$sigma_pos_grid`. alpha is recovered post-hoc as
  # sigma_pos / sigma_occ from the joint posterior.
  prior_for_joint <- prior
  prior_for_joint$spatial_idx <- NULL
  prior_for_joint$rho_bounds  <- NULL
  # Translate any legacy tau_grid that prior_from_spec attached to an
  # ICAR / CAR_proper prior — the joint engine takes sigma = 1/sqrt(tau)
  # as its donor-amplitude axis under the new parameterization.
  if (!is.null(prior_for_joint$tau_grid) &&
      is.null(prior_for_joint$sigma_grid)) {
    prior_for_joint$sigma_grid <- 1.0 / sqrt(as.numeric(prior_for_joint$tau_grid))
    prior_for_joint$tau_grid   <- NULL
  } else if (!is.null(prior_for_joint$tau_grid)) {
    prior_for_joint$tau_grid <- NULL
  }
  if (!is.null(control$sigma_grid))   prior_for_joint$sigma_grid   <- control$sigma_grid
  if (!is.null(control$rho_grid))     prior_for_joint$rho_grid     <- control$rho_grid
  if (!is.null(control$tau_grid)) {
    prior_for_joint$sigma_grid <- 1.0 / sqrt(as.numeric(control$tau_grid))
  }
  if (!is.null(control$rho_car_grid)) prior_for_joint$rho_car_grid <- control$rho_car_grid

  if (!is.null(control$sigma_pos_grid)) {
    sigma_pos_grid <- as.numeric(control$sigma_pos_grid)
  } else {
    # Default: mirror the donor sigma grid (gives alpha = sigma_pos /
    # sigma_occ posterior centered on 1.0 under flat per-axis priors).
    sigma_donor <- prior_for_joint$sigma_grid %||%
      exp(seq(log(0.1), log(3), length.out = 5))
    sigma_pos_grid <- as.numeric(sigma_donor)
  }

  copy_spec <- list(
    arm            = "pos",
    sigma_pos_grid = sigma_pos_grid
  )

  # ---- Multi-block path (Phase J-D) -----------------------------------
  # When `temporal` or `re` components are supplied, stack the spatial
  # block with AR1/RW/IID blocks and dispatch through the multi-block
  # joint engine. Copy semantics remain on the spatial block (sigma_occ /
  # sigma_pos), other blocks are shared identically across the two arms
  # (no per-arm scale).
  if (has_multi) {
    multi <- .cover_build_multi_prior(
      prior_spatial = prior_for_joint,
      spi_full      = spi_full,
      spi_pos       = spi_pos,
      data_obs      = data_obs,
      idx_pos       = enc$idx_pos,
      temporal      = temporal,
      re            = re,
      control       = control,
      sigma_pos_grid = sigma_pos_grid
    )
    # Strip spatial_idx from the arms — it lives inside the spatial
    # block's per-arm spatial_idx list in the multi-block prior.
    arm_occ$spatial_idx <- NULL
    arm_pos$spatial_idx <- NULL
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = multi$prior,
      copy      = multi$copy,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma_occ = control$prior_sigma_occ,
      prior_sigma_pos = control$prior_sigma_pos,
      max_iter  = control$max_iter  %||% 50L,
      tol       = control$tol       %||% 1e-6,
      n_threads = control$n_threads %||% 1L,
      store_Q   = TRUE,
      adaptive_grid             = control$adaptive_grid             %||% TRUE,
      adaptive_grid_edge_thresh = control$adaptive_grid_edge_thresh %||% 0.02,
      adaptive_grid_max_passes  = control$adaptive_grid_max_passes  %||% 1L
    )
  } else {
    # Adaptive grid forwarding. Defaults match the joint engine's defaults
    # (`adaptive_grid = TRUE`, threshold 0.02, one pass) and triggered the
    # under-coverage fix in INLAabun D3 — see gcol33/tulpaObs#8. Pass
    # `control$adaptive_grid = FALSE` to recover the legacy fixed-grid
    # behaviour for reproducibility checks.
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = prior_for_joint,
      copy      = copy_spec,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma_occ = control$prior_sigma_occ,
      prior_sigma_pos = control$prior_sigma_pos,
      max_iter  = control$max_iter  %||% 50L,
      tol       = control$tol       %||% 1e-6,
      n_threads = control$n_threads %||% 1L,
      store_Q   = TRUE,
      adaptive_grid             = control$adaptive_grid             %||% TRUE,
      adaptive_grid_edge_thresh = control$adaptive_grid_edge_thresh %||% 0.02,
      adaptive_grid_max_passes  = control$adaptive_grid_max_passes  %||% 1L
    )
  }

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

  # Dispersion summary on the positive arm. Both regimes integrate the
  # dispersion scalar on the outer joint hyperparameter grid; read the
  # posterior mean and SD from the engine's `theta_mean` / `theta_sd`. Those
  # are computed against the phi-axis marginal (foreign-axis slice cells
  # filtered out by `.joint_recalibrate_axis_moments`) with Laplace-at-mode
  # SD at the modal cell (gcol33/tulpa#20), so they are grid-spacing-
  # independent. Hand-rolling `sum(weights * theta_grid^2) - mean^2` against
  # `theta_grid[, "phi_pos"]` underestimates SD on sharply peaked axes and
  # additionally collapses on slice cells that pin phi at the modal value
  # while varying other axes -- that's the legacy pattern tulpa#20/#21 were
  # added to replace.
  #
  # The phi axis carries the gaussian residual SD for lognormal and the
  # beta precision for beta; surface under the respective slot names.
  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive == "lognormal") {
    sigma_pos    <- phi_mu
    sigma_pos_sd <- phi_sd
    phi_pos      <- NA_real_
    phi_pos_sd   <- NA_real_
  } else {
    sigma_pos    <- NA_real_
    sigma_pos_sd <- NA_real_
    phi_pos      <- phi_mu
    phi_pos_sd   <- phi_sd
  }

  m_occ <- list(mode = beta_occ, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  m_pos <- list(mode = beta_pos, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)

  # Stash the field-decomposition scale_factor (BYM2 Riebler scaling) on the
  # joint fit so the SLA path can reconstruct per-grid field amplitude
  # without re-deriving it. Dispersion is always integrated on `phi_pos`
  # (both lognormal and beta regimes), so the SLA path reads it directly
  # from `fit$theta_grid[k, "phi_pos"]` and needs no attr fallback.
  if (has_multi) {
    sf_attr <- as.numeric(multi$prior[[1L]]$scale_factor %||% 1.0)
  } else {
    sf_attr <- as.numeric(prior_for_joint$scale_factor %||% 1.0)
  }
  attr(fit, "scale_factor") <- sf_attr

  list(
    m_occ        = m_occ,
    m_pos        = m_pos,
    positive     = positive,
    sigma_pos    = sigma_pos,
    sigma_pos_sd = sigma_pos_sd,
    phi_pos      = phi_pos,
    phi_pos_sd   = phi_pos_sd,
    pos_fit_n    = N_pos,
    pos_fit_p    = p_pos,
    beta_occ     = beta_occ,
    beta_pos     = beta_pos,
    se_occ       = se_occ,
    se_pos       = se_pos,
    spi_full     = as.integer(spi_full),
    spi_pos      = as.integer(spi_pos),
    joint        = fit
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

#' Decode a joint-nested-Laplace cover-hurdle fit into a `cover_fit`.
#'
#' Lighter-weight than the single-Laplace decode: the joint engine has
#' already produced posterior moments for beta and the spatial hyperparameters,
#' so we just shape them into the existing `cover_fit` structure.
#'
#' When `approx = "simplified_laplace"`, the per-arm marginal skewness is
#' computed via [`.sla_compute_cover_hurdle_joint()`] (mixture third-moment
#' over the outer grid; per-grid FD of the joint inner log-lik along the
#' constraint-corrected Sigma columns), and per-arm pseudo-draws are
#' resampled from moment-matched skew-normals via
#' [`.sla_build_cover_hurdle_draws()`].
#'
#' @keywords internal
decode_cover_hurdle_joint <- function(fits, enc, family,
                                      approx = "gaussian_laplace") {
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
    hyperpar$sigma_pos    <- fits$sigma_pos
    hyperpar$sigma_pos_sd <- fits$sigma_pos_sd
  } else {
    hyperpar$phi_pos    <- fits$phi_pos
    hyperpar$phi_pos_sd <- fits$phi_pos_sd
  }

  # Simplified-Laplace marginal correction on the joint path. The SLA
  # orchestrator wants `enc$..spi_full` / `enc$..spi_pos` for per-arm field
  # gather inside the inner log-lik evaluator; we attach them here from
  # the fits list (computed once inside `fit_cover_hurdle_joint_nested`).
  skew_occ <- NULL
  skew_pos <- NULL
  draws_occ <- NULL
  draws_pos <- NULL
  sla_status <- "off"
  if (identical(approx, "simplified_laplace")) {
    enc_sla <- enc
    enc_sla$..spi_full <- as.integer(fits$spi_full %||% integer(0))
    enc_sla$..spi_pos  <- as.integer(fits$spi_pos  %||% integer(0))
    sla_res <- .sla_compute_cover_hurdle_joint(fits$joint, enc_sla,
                                               fits$positive)
    sla_draws <- .sla_build_cover_hurdle_draws(
      beta_occ, se_occ, beta_pos, se_pos, sla_res
    )
    draws_occ <- sla_draws$draws_occ
    draws_pos <- sla_draws$draws_pos
    sla_status <- sla_draws$sla_status
    if (isTRUE(sla_res$valid)) {
      skew_occ <- sla_res$gamma_occ
      skew_pos <- sla_res$gamma_pos
    } else {
      # The orchestrator may still return numeric (possibly non-finite)
      # gamma vectors alongside `valid = FALSE`; surface them only when
      # they are finite so downstream consumers can inspect them.
      if (!is.null(sla_res$gamma_occ) && all(is.finite(sla_res$gamma_occ))) {
        skew_occ <- sla_res$gamma_occ
      }
      if (!is.null(sla_res$gamma_pos) && all(is.finite(sla_res$gamma_pos))) {
        skew_pos <- sla_res$gamma_pos
      }
    }
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
      sigma_pos_sd = fits$sigma_pos_sd,
      phi_pos      = fits$phi_pos,
      phi_pos_sd   = fits$phi_pos_sd,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = length(enc$idx_pos),
      converged    = TRUE,
      log_marginal = c(joint = max(fits$joint$log_marginal)),
      joint        = fits$joint,
      skew_occ     = skew_occ,
      skew_pos     = skew_pos,
      draws_occ    = draws_occ,
      draws_pos    = draws_pos,
      sla_status   = sla_status
    ),
    class = c("cover_fit", "tobs_fit", "tulpa_fit")
  )
  out
}


# ---------------------------------------------------------------------------
# Multi-block prior assembly (Phase J-D)
# ---------------------------------------------------------------------------

# Build a multi-block joint prior (spatial + optional temporal + optional
# IID RE blocks) for tulpa::tulpa_nested_laplace_joint() under
# cover-hurdle copy semantics (copy on the spatial block).
#
# Non-spatial blocks are shared identically across the two arms — no
# per-arm scaling (INLA convention). This matches the typical cover
# hurdle use case: a year RE that influences both occurrence and cover
# magnitude in the same way, an observer RE that introduces a shared
# offset on both arms.
.cover_build_multi_prior <- function(prior_spatial, spi_full, spi_pos,
                                     data_obs, idx_pos, temporal, re,
                                     control, sigma_pos_grid) {
  # Spatial block — fill missing grids with defaults and attach per-arm
  # spatial_idx vectors. (Single-block path stores spi inside the arms;
  # multi-block puts it in the block.)
  sp <- prior_spatial
  if (is.null(sp$sigma_grid)) {
    sp$sigma_grid <- exp(seq(log(0.1), log(3), length.out = 5))
  }
  if (tolower(sp$type) == "bym2" && is.null(sp$rho_grid)) {
    sp$rho_grid <- c(0.25, 0.5, 0.75)
  }
  sp$spatial_idx <- list(as.integer(spi_full), as.integer(spi_pos))

  blocks <- list(sp)

  if (!is.null(temporal)) {
    blocks[[length(blocks) + 1L]] <- .cover_temporal_block(
      temporal, data_obs, idx_pos, control
    )
  }

  if (!is.null(re) && length(re) > 0L) {
    for (re_i in re) {
      blocks[[length(blocks) + 1L]] <- .cover_re_block(
        re_i, data_obs, idx_pos, control
      )
    }
  }

  list(
    prior = blocks,
    copy  = list(block = 1L, arm = "pos",
                 sigma_pos_grid = as.numeric(sigma_pos_grid))
  )
}

.cover_temporal_block <- function(temporal, data_obs, idx_pos, control) {
  if (!inherits(temporal, "tobs_temporal")) {
    stop("`temporal` must be a tobs_temporal() object.", call. = FALSE)
  }
  # tobs_temporal()'s `shared = c(TRUE, FALSE)` default was designed for
  # occupancy + detection (state vs. observation). cover() has two
  # likelihood arms (occurrence + cover magnitude) and the temporal term
  # enters both identically -- the `shared` field is ignored here.
  t_full <- .cover_resolve_idx(temporal$time, data_obs, "tobs_temporal$time")
  t_pos  <- t_full[idx_pos]
  n_t <- max(t_full)
  type <- temporal$type

  if (type == "ar1") {
    list(
      type         = "ar1",
      n_times      = as.integer(n_t),
      tau_grid     = as.numeric(control$tau_temporal_grid %||%
                                 c(1, 4, 16)),
      rho_grid     = as.numeric(control$rho_temporal_grid %||%
                                 c(0.3, 0.7)),
      temporal_idx = list(as.integer(t_full), as.integer(t_pos))
    )
  } else if (type == "iid") {
    list(
      type       = "iid",
      n_units    = as.integer(n_t),
      sigma_grid = as.numeric(control$sigma_temporal_grid %||%
                               exp(seq(log(0.1), log(1), length.out = 3))),
      obs_idx    = list(as.integer(t_full), as.integer(t_pos))
    )
  } else if (type %in% c("rw1", "rw2")) {
    list(
      type         = type,
      n_times      = as.integer(n_t),
      tau_grid     = as.numeric(control$tau_temporal_grid %||%
                                 c(1, 4, 16)),
      temporal_idx = list(as.integer(t_full), as.integer(t_pos))
    )
  } else {
    stop(sprintf("Unsupported tobs_temporal$type: '%s'", type),
         call. = FALSE)
  }
}

.cover_re_block <- function(re_i, data_obs, idx_pos, control) {
  if (!inherits(re_i, "tobs_re")) {
    stop("`re` elements must be tobs_re() objects.", call. = FALSE)
  }
  if (!identical(re_i$type, "intercept") && !identical(re_i$type, "iid")) {
    stop("cover() multi-block: tobs_re(type = 'intercept' | 'iid') is the ",
         "only supported config. Random slopes land in a later phase.",
         call. = FALSE)
  }
  if (!identical(re_i$model, "iid")) {
    stop("cover() multi-block: tobs_re(model = 'iid') is the only ",
         "supported temporal structure on RE blocks. AR1/RW1/RW2 on RE ",
         "land in a later phase.", call. = FALSE)
  }
  # Same as in .cover_temporal_block: `shared` is ignored in cover-hurdle
  # context. The RE term enters both arms identically.
  g_full <- .cover_resolve_idx(re_i$group, data_obs, "tobs_re$group")
  g_pos  <- g_full[idx_pos]
  n_g <- max(g_full)
  list(
    type       = "iid",
    n_units    = as.integer(n_g),
    sigma_grid = as.numeric(control$sigma_re_grid %||%
                             exp(seq(log(0.1), log(1.5), length.out = 3))),
    obs_idx    = list(as.integer(g_full), as.integer(g_pos))
  )
}

# Resolve a tobs component's `time` / `group` reference to a 1-based
# integer index vector of length nrow(data).
.cover_resolve_idx <- function(x, data, what) {
  if (is.character(x) && length(x) == 1L) {
    if (!x %in% names(data)) {
      stop(sprintf("`%s` = '%s' not found in data.", what, x),
           call. = FALSE)
    }
    v <- data[[x]]
    if (is.factor(v))         return(as.integer(v))
    if (is.character(v))      return(as.integer(as.factor(v)))
    if (is.numeric(v) || is.integer(v)) {
      iv <- as.integer(v)
      if (any(!is.finite(iv) | iv < 1L)) {
        stop(sprintf("`%s` = '%s' must be 1-based integer indices.",
                     what, x), call. = FALSE)
      }
      return(iv)
    }
    stop(sprintf("`%s` = '%s' must be a factor, character, or 1-based ",
                 "integer column.", what, x), call. = FALSE)
  }
  if (is.integer(x) || is.numeric(x)) {
    if (length(x) != nrow(data)) {
      stop(sprintf("length(%s) (%d) must equal nrow(data) (%d).",
                   what, length(x), nrow(data)), call. = FALSE)
    }
    iv <- as.integer(x)
    if (any(!is.finite(iv) | iv < 1L)) {
      stop(sprintf("`%s` must be 1-based integer indices.", what),
           call. = FALSE)
    }
    return(iv)
  }
  stop(sprintf("`%s` must be a column name or 1-based integer vector.",
               what), call. = FALSE)
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
# IMPORTANT: the same constraint is applied to BOTH BYM2 sub-blocks
# (phi: rank-deficient ICAR; theta: proper IID). This is a modelling
# choice rather than a mathematical necessity for theta — the IID
# prior identifies mean(theta) at precision n_s — but matching INLA's
# `f(..., model = "bym2", constr = TRUE)` default keeps the reported
# intercept comparable. Any simulator generating BYM2 data for the
# joint engine must demean both phi_f and theta_f before scaling,
# otherwise the constrained-intercept estimator targets
# `beta_pos_0_truth + alpha * mean(w_s_sim)` rather than the
# population truth and coverage of the population truth collapses
# with alpha (see `simulate_cover_joint()` for a demeaned simulator;
# diagnosed in INLAabun `example/validation/SUMMARY.md` Demo 3).
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
