# ============================================================================
# simplified_laplace.R — Simplified Laplace skewness correction for tobs fits
#
# Implements the third-cumulant correction of Rue, Martino & Chopin (2009 §3.2),
# adapted for the EM-Laplace M-step structure used by tulpaObs.
#
# Math:  see dev_notes/simplified_laplace_derivation.md
# Spec:  see dev_notes/upstream_tulpa_sla_spec.md
#
# Self-contained: needs no tulpa engine changes for fixed-effects-coefficient
# marginals. Uses H_beta and X from tulpa_laplace()'s return list directly.
#
# The skew-normal utilities here are temporary — they live in tulpaObs until
# tulpa absorbs them (see upstream spec §2.1).
# ============================================================================


# ---------------------------------------------------------------------------
# Per-site third derivatives, by family
#
# All are dl^3 / d eta^3 evaluated at eta_i, where eta_i is the linear predictor
# at site i under the *original* observation model (not the EM-encoded pseudo-
# binomial). The pseudo-binomial encoding scales l''' by M and is NOT what SLA
# needs — we evaluate against the original likelihood at the EM-converged
# parameter values.
# ---------------------------------------------------------------------------

#' Third derivative of binomial log-likelihood wrt logit(p)
#'
#' For y_i ~ Binomial(n_i, plogis(eta_i)),
#'   d^3 log p / d eta^3 = -n_i * p_i * (1 - p_i) * (1 - 2 * p_i)
#'
#' @param eta Linear predictor at sites (numeric vector).
#' @param n_trials Per-site trial counts (integer or numeric vector, recycled
#'   to length(eta) if scalar).
#' @return Vector of third derivatives (length = length(eta)).
#' @keywords internal
.l3_binomial_logit <- function(eta, n_trials = 1L) {
  p <- plogis(eta)
  -as.numeric(n_trials) * p * (1 - p) * (1 - 2 * p)
}

#' Third derivative of Poisson log-likelihood wrt log(lambda)
#'
#' For y_i ~ Poisson(exp(eta_i)),
#'   d^3 log p / d eta^3 = -lambda_i = -exp(eta_i)
#'
#' @keywords internal
.l3_poisson_log <- function(eta) {
  -exp(eta)
}

#' Third derivative of Gaussian log-likelihood (identity link)
#'
#' Always zero — the SLA correction reduces to the Laplace approximation
#' for Gaussian likelihoods. Required regression test.
#'
#' @keywords internal
.l3_gaussian_identity <- function(eta) {
  rep(0, length(eta))
}


# ---------------------------------------------------------------------------
# SLA assembly: gamma_j = sigma_j^{-3} * sum_i l3_i * (X Sigma)_{i,j}^3
#
# Diagonal-likelihood case (single linear predictor per site). For families
# with multiple linear predictors per site (e.g. occu with separate occ + det),
# this is applied once per process block.
# ---------------------------------------------------------------------------

#' Assemble per-coefficient skewness from per-site third derivatives
#'
#' Implements equation (2.2) of `dev_notes/simplified_laplace_derivation.md`:
#'   gamma_j = sigma_j^{-3} * sum_i l3_i * v_{i,j}^3,
#'   v_{i,j} = (X Sigma)_{i,j}
#'
#' For *diagonal-in-eta* likelihoods only (binomial, Poisson, Gaussian).
#'
#' @param l3 Per-site third derivative (length n).
#' @param X Design matrix (n x p).
#' @param Sigma Posterior covariance of the fixed-effect coefficients
#'   (p x p, must be symmetric PD).
#' @return Vector of length p with one skewness per coefficient.
#' @keywords internal
.sla_gamma_diag <- function(l3, X, Sigma) {
  stopifnot(length(l3) == nrow(X))
  stopifnot(nrow(Sigma) == ncol(X), ncol(Sigma) == ncol(X))

  sigma_j <- sqrt(diag(Sigma))
  XSig <- X %*% Sigma                                # n x p
  v3   <- XSig^3                                     # n x p (element-wise)

  # kappa_3_j = sum_i l3_i * v_{i,j}^3
  kappa3 <- as.numeric(crossprod(v3, l3))            # p-vector
  kappa3 / sigma_j^3
}


# ---------------------------------------------------------------------------
# Family-specific orchestrators
#
# Each takes the converged EM result + model metadata and returns a list:
#   list(gamma = named numeric, l3_diag = named numeric (debug),
#        valid = logical (FALSE if any process produced NaN/Inf))
#
# Process blocks ("occ", "det", ...) are concatenated in process_info order.
# ---------------------------------------------------------------------------

#' Compute SLA skewness coefficients for a single-season occu fit
#'
#' Evaluates l_i''' at the *original* single-season occupancy log-likelihood
#' at the EM-converged beta, NOT at the M-step pseudo-binomial encoding.
#'
#' Why this matters: the EM M-step encodes the occupancy block as a pseudo-
#' binomial with M = 1000 pseudo-trials per site. That encoding's Hessian
#' is M-inflated and would give Sigma_pseudo ~ Sigma_true / M, so any SLA
#' formula using it would be wrong by factor M^{3/2}. The correct posterior
#' precision for the occ block is the Louis observed information
#'   I_obs = X_occ' diag(psi(1-psi) - w(1-w)) X_occ  (+ prior penalty)
#' which the package already computes via `.louis_info_psi_single()`. We
#' use Sigma_occ = solve(I_obs).
#'
#' For the detection block: the EM M-step encodes a real weighted binomial
#' (y = n_det, n_trials = n_valid, weights = w_i), no pseudo-trial inflation,
#' so `fit_det$H_beta` is approximately the correct detection-posterior
#' precision and `Sigma_det = solve(fit_det$H_beta)` is used directly.
#'
#' For the occ process:
#'   - "any-detection" sites (z = 1 forced): l_i''' wrt logit(psi) is
#'     the Bernoulli l''' = -psi(1-psi)(1-2*psi).
#'   - "no-detection" sites: l_i = log(psi*q + (1-psi)) with q = prod(1-p).
#'     The closed-form third derivative is derived from
#'     l(eta) = log(1 - sigma(eta) * (1 - q_i)).
#'
#' @param model A `tobs_model` (model_type = "single").
#' @param em_result The EM-Laplace return list (with `$fits$occ`, `$fits$det`,
#'   `$weights`).
#' @param spatial Optional `tobs_spatial` (skewness disabled when set,
#'   pending Phase 3.5 / spatial Sigma exposure).
#' @param prior_spec Optional prior spec; passed to `.louis_info_psi_single()`
#'   so the penalty is included in I_obs.
#' @return List(gamma, valid, reason).
#' @keywords internal
.sla_compute_occu_single <- function(model, em_result, spatial = NULL,
                                     prior_spec = NULL) {
  if (!is.null(spatial)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "SLA on spatial Sigma not yet implemented (Phase 3.5)"))
  }
  if (is.null(em_result$weights)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "em_result$weights missing — needed for Louis Sigma_occ"))
  }

  fit_occ <- em_result$fits$occ
  fit_det <- em_result$fits$det
  if (is.null(fit_occ) || is.null(fit_det)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "EM fits missing for occ or det block"))
  }
  if (is.null(fit_det$H_beta)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "fit_det$H_beta missing — was return_hessian = TRUE?"))
  }

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  y     <- model$y

  pi_list <- model$process_info
  p_occ <- pi_list[[1]]$p
  p_det <- pi_list[[2]]$p

  beta_occ <- as.numeric(fit_occ$mode)[seq_len(p_occ)]
  beta_det <- as.numeric(fit_det$mode)[seq_len(p_det)]

  eta_occ <- as.numeric(X_occ %*% beta_occ)
  eta_det <- as.numeric(X_det %*% beta_det)
  psi <- plogis(eta_occ)
  p   <- plogis(eta_det)

  # Per-site detection structure
  n_sites <- model$n_sites
  n_valid <- integer(n_sites); n_det <- integer(n_sites)
  any_det <- logical(n_sites)
  for (i in seq_len(n_sites)) {
    valid <- y[i, ] >= 0
    n_valid[i] <- sum(valid)
    n_det[i]   <- sum(y[i, valid] == 1)
    any_det[i] <- n_det[i] > 0
  }

  # ---- Occupancy process l''' (at original likelihood) -------------------
  q <- (1 - p)^n_valid    # vectorised; uniform site-level p
  l3_occ <- numeric(n_sites)
  for (i in seq_len(n_sites)) {
    if (any_det[i]) {
      l3_occ[i] <- -psi[i] * (1 - psi[i]) * (1 - 2 * psi[i])
    } else {
      one_m_q <- 1 - q[i]
      pp      <- psi[i] * (1 - psi[i])
      up      <- pp * one_m_q
      upp     <- pp * (1 - 2 * psi[i]) * one_m_q
      uppp    <- pp * (1 - 6 * pp) * one_m_q
      v       <- 1 - psi[i] * one_m_q
      vp      <- -up;  vpp <- -upp;  vppp <- -uppp
      l3_occ[i] <- vppp / v - 3 * vpp * vp / v^2 + 2 * (vp / v)^3
    }
  }

  # Sigma_occ via Louis observed information (the M-inflation-free precision)
  I_obs <- .louis_info_psi_single(
    X_occ       = X_occ,
    beta_psi    = beta_occ,
    weights     = em_result$weights,
    spatial     = spatial,
    spatial_fit = fit_occ,
    prior_spec  = prior_spec,
    coef_names  = pi_list[[1]]$coef_names
  )
  Sigma_occ <- tryCatch(solve(I_obs), error = function(e) NULL)
  if (is.null(Sigma_occ)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "Louis I_obs (occ) not invertible"))
  }
  gamma_occ <- .sla_gamma_diag(l3 = l3_occ, X = X_occ, Sigma = Sigma_occ)

  # ---- Detection process l''' (any-det sites only) -----------------------
  # Detection likelihood: Binomial(n_valid, p) with observed n_det; only
  # any-det sites contribute (E-step weight on z=0 sites zeroes detection).
  l3_det_per_site <- -n_valid * p * (1 - p) * (1 - 2 * p)
  l3_det_per_site[!any_det] <- 0
  Sigma_det <- tryCatch(solve(fit_det$H_beta), error = function(e) NULL)
  if (is.null(Sigma_det)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "fit_det$H_beta not invertible"))
  }
  X_det_rows <- X_det
  if (nrow(X_det_rows) > n_sites) {
    X_det_rows <- X_det[seq_len(n_sites), , drop = FALSE]
  }
  gamma_det <- .sla_gamma_diag(l3 = l3_det_per_site, X = X_det_rows,
                               Sigma = Sigma_det)

  # ---- Concatenate in process_info order ---------------------------------
  gamma <- c(gamma_occ, gamma_det)
  occ_names <- paste0(pi_list[[1]]$name, "_", pi_list[[1]]$coef_names)
  det_names <- paste0(pi_list[[2]]$name, "_", pi_list[[2]]$coef_names)
  names(gamma) <- c(occ_names, det_names)

  list(gamma = gamma, valid = all(is.finite(gamma)),
       reason = if (all(is.finite(gamma))) "ok" else "non-finite gamma")
}


# ---------------------------------------------------------------------------
# Skew-normal utilities (temporary home — to be moved to tulpa per upstream
# spec §2.1)
# ---------------------------------------------------------------------------

# Skew-normal ceiling: max |skewness| representable by SN.
# = (4 - pi)/2 * (2/pi)^(3/2) / (1 - 2/pi)^(3/2)
.SN_GAMMA_MAX <- 0.9952717

#' Match cumulants (mu, sigma, gamma) to skew-normal (xi, omega, alpha)
#'
#' Inverse of the skew-normal moment formulas. Returns NULL when
#' |gamma| >= the skew-normal ceiling (~0.995). Caller falls back to
#' numerical-quadrature quantiles in that case.
#'
#' @param mu Mean (scalar).
#' @param sigma Standard deviation (positive scalar).
#' @param gamma Skewness coefficient (scalar).
#' @return List with elements xi, omega, alpha; or NULL.
#' @keywords internal
.sn_match <- function(mu, sigma, gamma) {
  if (!is.finite(gamma) || abs(gamma) >= .SN_GAMMA_MAX) return(NULL)
  c1 <- ((4 - pi) / 2)^(2 / 3)
  ag <- abs(gamma)^(2 / 3)
  delta_sq <- (pi / 2) * ag / (ag + c1)
  delta <- sign(gamma) * sqrt(delta_sq)
  one_m_2d2_over_pi <- 1 - 2 * delta^2 / pi
  if (one_m_2d2_over_pi <= 0) return(NULL)
  omega <- sigma / sqrt(one_m_2d2_over_pi)
  xi    <- mu - omega * delta * sqrt(2 / pi)
  alpha <- delta / sqrt(max(1 - delta^2, .Machine$double.eps))
  list(xi = xi, omega = omega, alpha = alpha)
}

#' Sample from a skew-normal distribution
#'
#' Two-component method (Azzalini 1985): Z1, Z2 ~ N(0,1) i.i.d., then
#'   Z = delta * |Z1| + sqrt(1 - delta^2) * Z2
#'   X = xi + omega * Z
#'
#' @param n Number of samples.
#' @param sn Skew-normal parameter list from `.sn_match()`.
#' @return Numeric vector of length n.
#' @keywords internal
.sn_sample <- function(n, sn) {
  if (is.null(sn)) return(rep(NA_real_, n))
  delta <- sn$alpha / sqrt(1 + sn$alpha^2)
  z1 <- abs(rnorm(n))
  z2 <- rnorm(n)
  z <- delta * z1 + sqrt(1 - delta^2) * z2
  sn$xi + sn$omega * z
}

#' Owen's T function via numerical quadrature
#'
#' T(h, a) = (1 / (2*pi)) * integral_0^a exp(-h^2 (1 + x^2) / 2) / (1 + x^2) dx
#'
#' Used in the skew-normal CDF:
#'   Phi_SN(x; xi, omega, alpha) = Phi(z) - 2 * T(z, alpha), z = (x - xi) / omega.
#'
#' Pure quadrature implementation. Closed-form series (Patefield-Tandy)
#' is faster but more error-prone; we only call this from `.sn_cdf()` which
#' is itself only called from `.sn_quantile()` (rare path, used for CIs
#' when SN representation kicks in). Speed is not the bottleneck.
#'
#' Identities used to handle endpoint cases:
#'   T(0, a) = atan(a) / (2*pi)
#'   T(h, 0) = 0
#'   T(h, a) = sign(a) * T(h, |a|)   (antisymmetric in a)
#'
#' @keywords internal
.owens_t <- function(h, a) {
  if (a == 0) return(0)
  if (h == 0) return(atan(a) / (2 * pi))
  sgn <- sign(a)
  a_abs <- abs(a)
  if (is.infinite(a_abs)) {
    return(sgn * 0.5 * (1 - pnorm(abs(h))))
  }
  integrand <- function(x) exp(-h^2 * (1 + x^2) / 2) / (1 + x^2)
  val <- tryCatch(stats::integrate(integrand, 0, a_abs, rel.tol = 1e-8)$value,
                  error = function(e) NA_real_)
  sgn * val / (2 * pi)
}

#' Skew-normal CDF
#' @keywords internal
.sn_cdf <- function(x, sn) {
  z <- (x - sn$xi) / sn$omega
  vapply(z, function(zi) pnorm(zi) - 2 * .owens_t(zi, sn$alpha), numeric(1))
}

#' Skew-normal quantile (Newton iteration)
#' @keywords internal
.sn_quantile <- function(p, sn) {
  if (is.null(sn)) return(rep(NA_real_, length(p)))
  delta <- sn$alpha / sqrt(1 + sn$alpha^2)
  # Start from Gaussian quantile centred at SN mean
  mean_sn <- sn$xi + sn$omega * delta * sqrt(2 / pi)
  sd_sn   <- sn$omega * sqrt(1 - 2 * delta^2 / pi)
  q0 <- qnorm(p, mean_sn, sd_sn)
  vapply(seq_along(p), function(k) {
    q <- q0[k]
    for (it in seq_len(40)) {
      f <- .sn_cdf(q, sn) - p[k]
      # SN density: 2 * phi(z) * Phi(alpha z) / omega
      z <- (q - sn$xi) / sn$omega
      dens <- 2 * dnorm(z) * pnorm(sn$alpha * z) / sn$omega
      if (!is.finite(dens) || dens < 1e-12) break
      step <- f / dens
      q <- q - step
      if (abs(step) < 1e-8) break
    }
    q
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# Replace per-coefficient draws with skew-normal samples
# ---------------------------------------------------------------------------

#' Replace Gaussian pseudo-draws with SN-sampled draws per coefficient
#'
#' For each parameter j, fit skew-normal (xi_j, omega_j, alpha_j) by moment-
#' matching (mu_j, sigma_j, gamma_j), and resample the column of `draws`.
#' Joint correlations are NOT preserved — only marginals are SLA-corrected,
#' matching INLA's behaviour (the SLA paper is about marginal quality).
#'
#' Behaviour at large |gamma|:
#'   - |gamma| < `cap`: SN draws via moment match.
#'   - |gamma| in [cap, SN_GAMMA_MAX): SN draws with gamma clipped to
#'     `sign(gamma) * cap`. The default `cap = 0.5` aligns with the
#'     validity envelope identified in dev_notes/simplified_laplace_derivation.md
#'     §2.6 — the cumulant expansion saturates above |gamma| ~ 0.5 and
#'     overstates magnitude. Clipping there avoids over-correcting CIs
#'     in the high-skew regime where SLA itself is unreliable.
#'   - |gamma| >= SN_GAMMA_MAX or NaN: Gaussian draws (no SLA correction).
#'
#' Returns the modified draws matrix with attributes:
#'   attr "sla_clipped"  — character vector of param names whose gamma was capped
#'   attr "sla_fallback" — character vector of param names that stayed Gaussian
#'
#' @keywords internal
.sla_replace_draws <- function(draws, means, sds, gamma, cap = 0.5) {
  n_pseudo <- nrow(draws)
  n_params <- ncol(draws)
  clipped <- character(0)
  fallback <- character(0)
  cap <- min(cap, .SN_GAMMA_MAX - 1e-3)
  for (j in seq_len(n_params)) {
    g <- gamma[j]
    if (!is.finite(g) || abs(g) < 1e-6) next                 # SLA is a no-op
    g_use <- g
    if (abs(g) > cap) {
      g_use <- sign(g) * cap
      clipped <- c(clipped, colnames(draws)[j])
    }
    sn <- .sn_match(means[j], sds[j], g_use)
    if (is.null(sn)) {
      fallback <- c(fallback, colnames(draws)[j])
      next
    }
    draws[, j] <- .sn_sample(n_pseudo, sn)
  }
  if (length(clipped))  attr(draws, "sla_clipped")  <- clipped
  if (length(fallback)) attr(draws, "sla_fallback") <- fallback
  draws
}
