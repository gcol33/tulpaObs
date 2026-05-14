# =============================================================================
# family_cover_hurdle.R — Vegetation cover hurdle on the tulpa Laplace backend
#
# Phase 1a: lognormal-positive variant via two independent tulpa_laplace()
# calls (binomial on occur, gaussian on log(cover[occur==1])). The joint
# shared-field model with proper hyperparameter integration is Phase 1c;
# the beta-positive variant is Phase 1d.
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
#' @param engine Currently `"laplace"` is the only supported value.
#' @param priors Currently ignored — passed through for forward compat.
#' @param control List with optional `max_iter`, `tol`, `n_threads`.
#' @return List with `m_occ`, `m_pos`, `positive`, `pos_fit_n`, `pos_fit_p`,
#'   plus one of `sigma_pos` (lognormal) or `phi_pos` (beta).
#' @keywords internal
fit_cover_hurdle <- function(enc, positive = enc$positive,
                             engine = "laplace",
                             priors = NULL, control = list()) {
  if (!engine %in% c("laplace", "auto")) {
    stop("cover() currently supports only engine = 'laplace' ",
         "(got '", engine, "'). nested_laplace / nuts land in later phases.",
         call. = FALSE)
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
# Helpers
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

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
