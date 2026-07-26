# =============================================================================
# sla_cover_hurdle.R — Simplified-Laplace skewness correction for cover hurdle
#
# Phase 3.5 extension. The cover hurdle has two independent arms (no shared
# spatial field in the single-Laplace path):
#
#   1. Occurrence arm: y_indicator_i ~ Bernoulli(plogis(X_occ_i beta_occ))
#   2. Positive arm:   y_i | y_i > 0 ~ Beta(mu_i phi, (1-mu_i) phi)
#                                       or Lognormal(mu_i, sigma_pos^2)
#
# Both arms run the *real* likelihood at the mode (no M-step pseudo-binomial
# encoding), so the Sigma used for SLA gamma is the raw `solve(H_beta)`
# returned by tulpa_laplace / tulpa_laplace_beta — no Louis correction needed.
#
# The two arms are independent given y_indicator, so the joint gamma is
# (gamma_occ, gamma_pos) concatenated. We compute each block separately,
# using `.sla_gamma_fd()` from `R/simplified_laplace.R` against the
# observation log-likelihood for that arm. phi (beta arm) / sigma_pos
# (lognormal arm) are held fixed at the fitted values: SLA gamma is only
# for the fixed-effect coefficients beta_occ and beta_pos.
# =============================================================================


# ---------------------------------------------------------------------------
# Arm-wise log-likelihood evaluators
# ---------------------------------------------------------------------------

#' Occurrence-arm log-likelihood (Bernoulli)
#'
#' Computes the per-arm log-likelihood for the cover-hurdle occurrence
#' indicator at a candidate `beta_occ`, with no prior or pseudo-binomial
#' encoding. Used by [`.sla_gamma_fd()`].
#'
#' @param beta_occ Length-`p_occ` numeric coefficient vector.
#' @param enc Encoded data list from [`encode_cover_hurdle()`] (uses
#'   `enc$occ_data$y` and `enc$occ_data$X`).
#' @return Scalar log-likelihood.
#' @keywords internal
.loglik_cover_occ <- function(beta_occ, enc) {
  y <- enc$occ_data$y
  X <- enc$occ_data$X
  eta <- as.numeric(X %*% beta_occ)
  # Bernoulli log-likelihood using log1p for numerical stability.
  #   y log(p) + (1-y) log(1-p)
  # = y * eta - log1p(exp(eta))
  ll <- sum(y * eta - log1p(exp(eta)))
  ll
}

#' Beta-arm log-likelihood (positive subset only)
#'
#' Beta log-density at the encoded positive responses with logit link,
#' summed over the positive subset. See `dev_notes/` for the closed form;
#' phi is treated as a fixed nuisance parameter (held at the fitted value)
#' and SLA gamma is computed for beta_pos only.
#'
#' @param beta_pos Length-`p_pos` numeric coefficient vector.
#' @param phi Beta precision (positive scalar; held fixed at the fitted value).
#' @param enc Encoded data list from `encode_cover_hurdle()` (uses
#'   `enc$pos_data$y` and `enc$pos_data$X`).
#' @return Scalar log-likelihood.
#' @keywords internal
.loglik_cover_pos_beta <- function(beta_pos, phi, enc) {
  y <- enc$pos_data$y
  X <- enc$pos_data$X
  if (length(y) == 0L) return(0)
  mu  <- plogis(as.numeric(X %*% beta_pos))
  a   <- mu * phi
  b   <- (1 - mu) * phi
  ll <- sum(lgamma(phi) - lgamma(a) - lgamma(b) +
              (a - 1) * log(y) + (b - 1) * log1p(-y))
  ll
}

#' Lognormal-arm log-likelihood (positive subset only)
#'
#' Note: the encoded `enc$pos_data$y` is already on the log scale (i.e.
#' `log(cover[occur == 1])`; see `encode_cover_hurdle()`). On the log
#' scale the lognormal density reduces to a Gaussian on `log y` with
#' mean `X beta_pos` and SD `sigma_pos`. The Jacobian `- log y` is
#' constant in `beta_pos` and so does not affect the third derivative;
#' we drop it. (For a finite-difference d3 along the beta_pos direction,
#' beta-independent terms cancel.)
#'
#' @param beta_pos Length-`p_pos` numeric coefficient vector.
#' @param sigma_pos Positive scalar (held fixed at the fitted value).
#' @param enc Encoded data list from `encode_cover_hurdle()`.
#' @return Scalar log-likelihood.
#' @keywords internal
.loglik_cover_pos_lognormal <- function(beta_pos, sigma_pos, enc) {
  z <- enc$pos_data$y   # already log(cover)
  X <- enc$pos_data$X
  if (length(z) == 0L) return(0)
  mu_z <- as.numeric(X %*% beta_pos)
  # Drop the constant Jacobian term -log(cover); it's beta-independent and
  # cancels under finite differences along beta_pos.
  ll <- sum(dnorm(z, mean = mu_z, sd = sigma_pos, log = TRUE))
  ll
}


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

#' Compute SLA gamma for cover hurdle arms
#'
#' Returns per-coefficient SLA gamma for each of the two arms of a cover
#' hurdle fit. Each arm's gamma is computed via [`.sla_gamma_fd()`] using
#' the arm's `solve(H_beta)` as Sigma (no Louis correction — both arms are
#' real likelihoods at the mode under the single-Laplace path).
#'
#' phi (beta arm) / sigma_pos (lognormal arm) are treated as fixed nuisance
#' parameters: SLA gamma is computed only for the fixed-effect coefficients.
#'
#' @param fits The list returned by [`fit_cover_hurdle()`] (with `m_occ`,
#'   `m_pos`, `positive`, and one of `phi_pos` / `sigma_pos`).
#' @param enc The encoded data from [`encode_cover_hurdle()`].
#' @param positive One of `"beta"`, `"lognormal"`.
#' @return List(gamma_occ, gamma_pos, valid, reason).
#' @keywords internal
.sla_compute_cover_hurdle <- function(fits, enc, positive) {
  bail  <- .sla_bailer(c("gamma_occ", "gamma_pos"))
  p_occ <- ncol(enc$occ_data$X)
  p_pos <- ncol(enc$pos_data$X)
  if (is.null(fits$m_occ) || is.null(fits$m_pos)) {
    return(bail("missing occ or pos arm fit"))
  }
  H_occ <- fits$m_occ$H_beta
  H_pos <- fits$m_pos$H_beta
  if (is.null(H_occ) || is.null(H_pos)) {
    return(bail("H_beta missing on occ or pos arm"))
  }

  beta_occ <- as.numeric(fits$m_occ$mode)[seq_len(p_occ)]
  beta_pos <- as.numeric(fits$m_pos$mode)[seq_len(p_pos)]

  Sigma_occ <- .sla_solve(H_occ)
  if (is.null(Sigma_occ)) return(bail("H_beta (occ) not invertible"))
  # The Gaussian-arm Hessian is computed under phi = 1; scale it back up
  # to the true noise scale so Sigma_pos = sigma_pos^2 * solve(H_beta) is
  # the right marginal precision. For the beta arm tulpa_laplace_beta()
  # already weights H_beta by phi (Fisher info), so no scaling needed.
  scale_pos <- if (identical(positive, "lognormal")) {
    sig <- fits$sigma_pos
    if (is.null(sig) || !is.finite(sig) || sig <= 0) 1 else sig^2
  } else 1
  Sigma_pos <- .sla_solve(H_pos)
  if (is.null(Sigma_pos)) return(bail("H_beta (pos) not invertible"))
  Sigma_pos <- scale_pos * Sigma_pos

  # Occurrence arm: Bernoulli FD gamma
  gamma_occ <- tryCatch(
    .sla_gamma_fd(beta_occ, Sigma_occ,
                  function(b) .loglik_cover_occ(b, enc)),
    error = function(e) NULL
  )
  if (is.null(gamma_occ) || !all(is.finite(gamma_occ))) {
    return(bail("FD gamma on occ arm produced non-finite values"))
  }
  names(gamma_occ) <- colnames(enc$occ_data$X)

  # Positive arm: Beta or Lognormal FD gamma
  gamma_pos <- if (identical(positive, "beta")) {
    phi <- fits$phi_pos
    if (is.null(phi) || !is.finite(phi) || phi <= 0) {
      return(bail("phi_pos missing or non-positive", gamma_occ = gamma_occ))
    }
    tryCatch(
      .sla_gamma_fd(beta_pos, Sigma_pos,
                    function(b) .loglik_cover_pos_beta(b, phi, enc)),
      error = function(e) NULL
    )
  } else if (identical(positive, "lognormal")) {
    sig <- fits$sigma_pos
    if (is.null(sig) || !is.finite(sig) || sig <= 0) {
      return(bail("sigma_pos missing or non-positive", gamma_occ = gamma_occ))
    }
    tryCatch(
      .sla_gamma_fd(beta_pos, Sigma_pos,
                    function(b) .loglik_cover_pos_lognormal(b, sig, enc)),
      error = function(e) NULL
    )
  } else {
    return(bail(sprintf("unsupported positive family '%s'", positive),
                gamma_occ = gamma_occ))
  }
  if (is.null(gamma_pos) || !all(is.finite(gamma_pos))) {
    return(bail("FD gamma on pos arm produced non-finite values",
                gamma_occ = gamma_occ))
  }
  names(gamma_pos) <- colnames(enc$pos_data$X)

  list(gamma_occ = gamma_occ, gamma_pos = gamma_pos, valid = TRUE,
       reason = "ok")
}


# ---------------------------------------------------------------------------
# Resample draws on both arms
#
# For each cover_fit we keep two parallel sets of pseudo-draws (occ + pos)
# in `fit$draws_occ` / `fit$draws_pos`. Joint correlations across arms are
# not modelled (the two arms are independent under the single-Laplace path);
# joint correlations within each arm are not preserved either (matches the
# single-season occu SLA behaviour).
# ---------------------------------------------------------------------------

#' Build SN-resampled per-arm pseudo-draws for a cover hurdle fit
#'
#' Returns `list(draws_occ, draws_pos, sla_status)`. If `sla_res$valid`
#' is `FALSE`, draws fall back to Gaussian samples around the per-arm
#' MAP +- SE.
#'
#' @keywords internal
.sla_build_cover_hurdle_draws <- function(beta_occ, se_occ,
                                          beta_pos, se_pos,
                                          sla_res,
                                          V_occ = NULL, V_pos = NULL,
                                          n_pseudo = 1000L) {
  # Draw each arm from its full per-arm covariance (V_occ / V_pos) so the
  # coefficient correlation propagates to predicted cover; fall back to the
  # diagonal of the reported SEs only when the covariance is unavailable
  # (gcol33/tulpaObs#45). The subsequent .sla_replace_draws() preserves the
  # joint rank-correlation while imposing the skew-normal marginals.
  draw_arm <- function(beta, se, V) {
    p <- length(beta)
    if (!is.null(V) && all(dim(as.matrix(V)) == p)) {
      d <- .rmvn(n_pseudo, beta, (as.matrix(V) + t(as.matrix(V))) / 2)
    } else {
      d <- matrix(NA_real_, n_pseudo, p)
      for (j in seq_len(p)) d[, j] <- rnorm(n_pseudo, beta[j], max(se[j], 1e-4))
    }
    colnames(d) <- names(beta)
    d
  }
  draws_occ <- draw_arm(beta_occ, se_occ, V_occ)
  draws_pos <- draw_arm(beta_pos, se_pos, V_pos)

  if (isTRUE(sla_res$valid)) {
    draws_occ <- .sla_replace_draws(draws_occ, beta_occ, se_occ,
                                    sla_res$gamma_occ)
    draws_pos <- .sla_replace_draws(draws_pos, beta_pos, se_pos,
                                    sla_res$gamma_pos)
    sla_status <- "simplified_laplace"
  } else {
    sla_status <- paste0("fallback_gaussian (", sla_res$reason, ")")
  }
  list(draws_occ = draws_occ, draws_pos = draws_pos,
       sla_status = sla_status)
}
