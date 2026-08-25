# =============================================================================
# sla_cover_hurdle.R — Simplified-Laplace skewness correction for cover hurdle
#
# The cover hurdle has two independent arms (no shared spatial field in the
# single-Laplace path):
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
  eta <- .tobs_clamp_eta(as.numeric(X %*% beta_occ))
  # Bernoulli log-likelihood:
  #   y log(p) + (1-y) log(1-p) = y * eta - log(1 + exp(eta))
  # `.sla_gamma_fd()` steps beta by +-2h along a Sigma column, so eta is
  # evaluated at displaced coefficients and can saturate the predictor on data
  # the fitted mode does not; the clamp and the stable log-normalizer keep the
  # displaced value finite so the finite difference is defined.
  ll <- sum(y * eta - .tobs_log1pexp(eta))
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
  mu  <- plogis(.tobs_clamp_eta(as.numeric(X %*% beta_pos)))
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
