# ============================================================================
# simplified_laplace.R -- Simplified Laplace skewness correction for tobs fits
#
# Implements the third-cumulant correction of Rue, Martino & Chopin (2009 sec.3.2),
# adapted for the EM-Laplace M-step structure used by tulpaObs.
#
# Math:  see dev_notes/simplified_laplace_derivation.md
# Spec:  see dev_notes/upstream_tulpa_sla_spec.md
#
# Self-contained: needs no tulpa engine changes for fixed-effects-coefficient
# marginals. Uses H_beta and X from tulpa_laplace()'s return list directly.
#
# Cumulant -> skew-normal conversion (sn_match, sn_cdf, sn_quantile) lives
# upstream in tulpa; only the four-line Azzalini sampler stays local.
# ============================================================================


# ---------------------------------------------------------------------------
# Per-site third derivatives, by family
#
# All are dl^3 / d eta^3 evaluated at eta_i, where eta_i is the linear predictor
# at site i under the *original* observation model (not the EM-encoded pseudo-
# binomial). The pseudo-binomial encoding scales l''' by M and is NOT what SLA
# needs -- we evaluate against the original likelihood at the EM-converged
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
#' Always zero -- the SLA correction reduces to the Laplace approximation
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
#'   gamma_j = sigma_j^(-3) * sum_i l3_i * v_(i,j)^3,
#'   v_(i,j) = (X Sigma)_(i,j)
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
# Generic finite-difference d3 for non-diagonal families
#
# For families whose log-likelihood does not decompose as Sum_i l_i(eta_i)
# with eta_i = X_i beta (HMM dyn_occu, integrated shared-psi, cover hurdle),
# the analytical formula above does not apply: off-diagonal third derivatives
# in eta-space are non-zero. We compute kappa_3[beta_j] directly via finite
# difference along the direction v_j = Sigma[, j] in beta-space:
#
#   kappa_3[beta_j] = d^3/dh^3 l(beta_hat + h v_j) |_{h=0}
#
# Proof: l(beta_hat + h v_j) Taylor-expanded in beta gives, by chain rule,
#   d^3/dh^3 l = Sum_{a,b,c} d^3 l / (d beta_a d beta_b d beta_c)
#                * v_{a,j} v_{b,j} v_{c,j}
# = Sum_{a,b,c} d^3 l / (d beta_a d beta_b d beta_c) Sigma_{a,j} Sigma_{b,j} Sigma_{c,j}
# = kappa_3[beta_j] under the Laplace-approximated posterior (see
#   dev_notes/simplified_laplace_derivation.md eq (2.1)).
#
# Then gamma_j = kappa_3[beta_j] / sigma_j^3.
#
# Cost: O(p) likelihood evaluations (5-point central difference), all at
# the mode + small displacement. With a fast R or compiled log-likelihood
# evaluator this is feasible for the p <= 30 fixed-effect dimensions
# typical of tulpaObs fits.
# ---------------------------------------------------------------------------

#' Generic finite-difference SLA gamma for non-diagonal families
#'
#' Computes per-coefficient SLA gamma by 5-point central finite difference
#' of the original observation log-likelihood along the direction
#' `Sigma[, j]` in beta-space. Use this when the family's third derivative
#' in eta-space does not decompose into a per-site sum (e.g. HMM forward
#' likelihood, integrated shared-process, hurdle joints).
#'
#' Step size defaults to `eps^{1/5} * sigma_j / ||v_j||`, which keeps
#' the displacement in beta-space on the natural scale of the j-th
#' marginal posterior. Override `h` (scalar or length-p vector) to inspect
#' truncation/cancellation behaviour.
#'
#' @param beta_hat Mode of the joint posterior (length p).
#' @param Sigma Posterior covariance at the mode (p x p, symmetric PD).
#' @param log_lik_fn Callable: takes a length-p beta vector, returns
#'   scalar log-likelihood at that beta (no prior contribution -- SLA
#'   gamma uses the observation likelihood only, per RMC 2009 sec.3.2).
#' @param h Optional step size override (scalar or length-p vector).
#' @return Numeric vector of length p with per-coefficient gamma.
#' @keywords internal
.sla_gamma_fd <- function(beta_hat, Sigma, log_lik_fn, h = NULL) {
  p <- length(beta_hat)
  stopifnot(nrow(Sigma) == p, ncol(Sigma) == p)

  sigma_vec <- sqrt(diag(Sigma))
  eps_h <- .Machine$double.eps^(1 / 5)
  if (!is.null(h)) {
    h <- rep_len(as.numeric(h), p)
  }

  gamma <- numeric(p)
  for (j in seq_len(p)) {
    vj  <- Sigma[, j]
    nvj <- sqrt(sum(vj^2))
    hj  <- if (is.null(h)) {
      eps_h * sigma_vec[j] / max(nvj, .Machine$double.eps)
    } else h[j]

    L_p2 <- log_lik_fn(beta_hat + 2 * hj * vj)
    L_p1 <- log_lik_fn(beta_hat +     hj * vj)
    L_m1 <- log_lik_fn(beta_hat -     hj * vj)
    L_m2 <- log_lik_fn(beta_hat - 2 * hj * vj)

    d3 <- (L_p2 - 2 * L_p1 + 2 * L_m1 - L_m2) / (2 * hj^3)
    gamma[j] <- d3 / sigma_vec[j]^3
  }
  gamma
}


# ---------------------------------------------------------------------------
# R-side log-likelihood evaluators (used by FD path and as ground-truth
# cross-checks against the analytical .sla_gamma_diag in tests).
#
# All evaluators take beta = c(beta_occ, beta_det, ...) in process_info
# order and return the *original* observation log-likelihood (no priors,
# no pseudo-binomial encoding).
# ---------------------------------------------------------------------------

#' R-side single-season occupancy log-likelihood
#'
#' Computes the marginal observation log-likelihood
#'   `log P(y | beta) = Sum_i log[ psi_i * P(y_i | z=1) + (1 - psi_i) I(no detections) ]`
#' with psi_i = plogis(X_occ_i beta_occ), p_i = plogis(X_det_i beta_det).
#'
#' Used to cross-check `.sla_gamma_fd()` against the analytical
#' `.sla_gamma_diag()` on single-season fits.
#'
#' @keywords internal
.loglik_occu_single <- function(beta, model) {
  pi_list <- model$process_info
  p_occ <- pi_list[[1]]$p
  p_det <- pi_list[[2]]$p
  beta_occ <- beta[seq_len(p_occ)]
  beta_det <- beta[p_occ + seq_len(p_det)]

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  y     <- model$y
  n_sites <- model$n_sites

  psi <- plogis(as.numeric(X_occ %*% beta_occ))
  p   <- plogis(as.numeric(X_det %*% beta_det))

  ll <- 0
  for (i in seq_len(n_sites)) {
    valid <- y[i, ] >= 0
    n_v <- sum(valid)
    if (n_v == 0L) next
    n_d <- sum(y[i, valid] == 1L)
    if (n_d > 0L) {
      # any-detection: z=1 forced
      ll <- ll + log(psi[i]) + n_d * log(p[i]) + (n_v - n_d) * log1p(-p[i])
    } else {
      # all zeros: marginalise over z
      q_i <- (1 - p[i])^n_v
      ll <- ll + log(psi[i] * q_i + (1 - psi[i]))
    }
  }
  ll
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
#' formula using it would be wrong by factor M^(3/2). The correct posterior
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
#' @param spatial Optional `tobs_spatial`. When set, the skewness correction is
#'   intentionally NOT applied and the fit keeps Gaussian marginals -- this is the
#'   correct conservative behaviour, not a stub: the
#'   third-cumulant correction (Rue, Martino & Chopin 2009 sec.3.2) captures the
#'   skewness of the fixed-effect (hyperparameter-free) marginals, but for a
#'   spatial latent field the dominant marginal skewness comes from integrating
#'   over the field-precision hyperparameter, which a correction conditioned on a
#'   single hyperparameter value does not capture. Validated against NUTS: every
#'   simplified-Laplace construction (modal-hyper, grid-mixture, mixture of
#'   skew-normals) disagreed with the NUTS posterior skewness in sign and/or
#'   magnitude, so applying one would be worse than the Gaussian fallback.
#' @param prior_spec Optional prior spec; passed to `.louis_info_psi_single()`
#'   so the penalty is included in I_obs.
#' @return List(gamma, valid, reason).
#' @keywords internal
.sla_compute_occu_single <- function(model, em_result, spatial = NULL,
                                     prior_spec = NULL) {
  if (!is.null(spatial)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = paste0(
                  "Gaussian marginals retained for spatial Sigma by design: the ",
                  "simplified-Laplace third-cumulant correction is valid for ",
                  "hyperparameter-free fixed-effect marginals only; for a spatial ",
                  "field the dominant skewness is hyperparameter-marginalisation, ",
                  "which it does not capture (validated against NUTS, ",
                  "gcol33/tulpaObs#55).")))
  }
  if (is.null(em_result$weights)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "em_result$weights missing -- needed for Louis Sigma_occ"))
  }

  fit_occ <- em_result$fits$occ
  fit_det <- em_result$fits$det
  if (is.null(fit_occ) || is.null(fit_det)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "EM fits missing for occ or det block"))
  }
  if (is.null(fit_det$H_beta)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "fit_det$H_beta missing -- was return_hessian = TRUE?"))
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
# Skew-normal sampling
#
# Cumulant -> skew-normal matching and quantile/CDF evaluation live upstream
# in tulpa (sn_match, sn_cdf, sn_quantile, exported since the upstream SLA
# spec in dev_notes/upstream_tulpa_sla_spec.md sec.2.1 landed). Sampling is not
# upstream: it has no use outside marginal resampling and is a four-line
# implementation of Azzalini's (1985) two-component construction.
# ---------------------------------------------------------------------------

#' Sample from a skew-normal distribution
#'
#' Two-component method (Azzalini 1985): Z1, Z2 ~ N(0,1) i.i.d., then
#'   Z = delta * |Z1| + sqrt(1 - delta^2) * Z2
#'   X = xi + omega * Z
#'
#' @param n Number of samples.
#' @param sn Skew-normal parameter list from [tulpa::sn_match()] (elements
#'   `xi`, `omega`, `alpha`).
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


# ---------------------------------------------------------------------------
# Replace per-coefficient draws with skew-normal samples
# ---------------------------------------------------------------------------

#' Replace Gaussian pseudo-draws with SN-sampled draws per coefficient
#'
#' For each parameter j, fit skew-normal (xi_j, omega_j, alpha_j) by moment-
#' matching (mu_j, sigma_j, gamma_j), and replace the column of `draws` with
#' skew-normal samples arranged in the same rank order as the incoming column.
#' The incoming draws are the correlated Gaussian pseudo-draws, so the
#' reordering is a Gaussian-copula transform: the skew-normal marginals are
#' exact and the joint rank-correlation of the Laplace covariance is preserved
#' (a derived quantity such as predicted psi keeps the coefficient dependence).
#' A column whose gamma is a no-op keeps its correlated
#' Gaussian draw unchanged.
#'
#' Behaviour at large |gamma|:
#'   - |gamma| < `cap`: SN draws via moment match.
#'   - |gamma| in [cap, SN_GAMMA_MAX): SN draws with gamma clipped to
#'     `sign(gamma) * cap`. The default `cap = 0.5` aligns with the
#'     validity envelope identified in dev_notes/simplified_laplace_derivation.md
#'     sec.2.6 -- the cumulant expansion saturates above |gamma| ~ 0.5 and
#'     overstates magnitude. Clipping there avoids over-correcting CIs
#'     in the high-skew regime where SLA itself is unreliable.
#'   - |gamma| >= SN_GAMMA_MAX or NaN: Gaussian draws (no SLA correction).
#'
#' Returns the modified draws matrix with attributes:
#'   attr "sla_clipped"  -- character vector of param names whose gamma was capped
#'   attr "sla_fallback" -- character vector of param names that stayed Gaussian
#'
#' @keywords internal
.sla_replace_draws <- function(draws, means, sds, gamma, cap = 0.5) {
  n_pseudo <- nrow(draws)
  n_params <- ncol(draws)
  clipped <- character(0)
  fallback <- character(0)
  # cap stays well below the SN ceiling (~0.995), so tulpa::sn_match's
  # ceiling-warning path is never triggered after clipping. We still
  # suppress warnings defensively so a degenerate sigma slip doesn't surface
  # as a noisy fit-side warning.
  for (j in seq_len(n_params)) {
    g <- gamma[j]
    if (!is.finite(g) || abs(g) < 1e-6) next                 # SLA is a no-op
    g_use <- g
    if (abs(g) > cap) {
      g_use <- sign(g) * cap
      clipped <- c(clipped, colnames(draws)[j])
    }
    sn <- tryCatch(
      suppressWarnings(tulpa::sn_match(means[j], sds[j], g_use)),
      error = function(e) NULL
    )
    if (is.null(sn)) {
      fallback <- c(fallback, colnames(draws)[j])
      next
    }
    # Rank-reorder skew-normal samples onto the column's existing (correlated)
    # ordering: exact SN marginal, joint rank-correlation preserved.
    sn_samp <- .sn_sample(n_pseudo, sn)
    draws[order(draws[, j]), j] <- sort(sn_samp)
  }
  if (length(clipped))  attr(draws, "sla_clipped")  <- clipped
  if (length(fallback)) attr(draws, "sla_fallback") <- fallback
  draws
}
