# ============================================================================
# sla_int_occu.R — Simplified-Laplace skewness for integrated multi-source
# occupancy (shared psi across data sources).
#
# Math
# ----
# Integrated occupancy couples each site i across S detection sources via a
# shared latent z_i with prior P(z_i = 1) = psi_i. Source s contributes
# detections y_{i,s,t} ~ Bernoulli(p_{i,s}) when z_i = 1, all zeros when z_i = 0.
# The marginal log-likelihood at site i is:
#
#   any source detected at i:
#     l_i = log psi_i + Sum_s [ y_{i,s,.} log p_{i,s} + (n_v - y_{i,s,.}) log(1 - p_{i,s}) ]
#   no source detected anywhere at i:
#     l_i = log( psi_i * Prod_s (1 - p_{i,s})^{n_v(i,s)} + (1 - psi_i) )
#
# The "no detection" branch is NOT diagonal in beta-space: the d^3 cross-
# couples beta_psi with every beta_p,s because they all enter the same scalar
# A_i = psi_i * Prod_s q_{i,s} + (1 - psi_i). The analytical diagonal SLA
# formula `.sla_gamma_diag` does not apply. We therefore route int_occu
# through `.sla_gamma_fd` (5-point central FD along Sigma[, j]) using the
# original-likelihood evaluator below.
#
# Sigma assembly
# --------------
# Block-diagonal joint Sigma:
#   occ block: Louis-corrected observed information
#     I_obs(beta_psi) = X' diag(psi(1-psi) - w(1-w)) X (+ prior penalty),
#     re-using the single-season helper. For int_occu, `weights` is the
#     converged E-step posterior P(z_i = 1 | y, theta_hat), provided by
#     `build_integrated_callbacks$e_step` and stored in `em_result$weights`.
#   det_s blocks: `solve(em_result$fits$det<s>$H_beta)`. The integrated
#     M-step encoding for each source is a weighted binomial (y = n_det_s,
#     n_trials = n_valid_s, weights = w_i) with NO M-inflation, so the
#     inner-fit Hessian is approximately the correct posterior precision
#     for that detection block.
#
# Off-diagonal blocks (cross-block covariance) are dropped — this matches
# what `.sla_compute_occu_single` does, and is consistent with the
# block-coordinate structure of the M-step.
# ============================================================================


#' Integrated multi-source occupancy log-likelihood (R-side)
#'
#' Computes the marginal observation log-likelihood
#' \verb{log P(y | beta) = Sum_i log( psi_i * Prod_s P(y_is | z=1, p_is)
#' + (1 - psi_i) * I(no detections at any source) )}
#' where \verb{psi_i = plogis(X_occ[i,] beta_psi)},
#' \verb{p_is = plogis(X_det_s[i,] beta_det_s)}.
#'
#' Reads `model$y_sources` (list of integer matrices, source-local row indexing,
#' NA encoded as -1), `model$site_maps` (0-indexed global site per source row),
#' `model$X_processes` (occ design + one det design per source, all n_sites x p),
#' `model$n_sources`, `model$n_sites`. No priors, no pseudo-binomial encoding.
#'
#' @param beta Joint parameter vector in process_info order:
#'   `c(beta_psi, beta_det_1, beta_det_2, …, beta_det_S)`.
#' @param model A `tobs_model` with `model_type = "integrated"`.
#' @return Scalar log-likelihood.
#' @keywords internal
.loglik_int_occu <- function(beta, model) {
  pi_list   <- model$process_info
  n_sources <- model$n_sources
  n_sites   <- model$n_sites

  p_occ    <- pi_list[[1]]$p
  beta_psi <- beta[seq_len(p_occ)]

  X_occ <- model$X_processes[[1]]
  eta_psi <- as.numeric(X_occ %*% beta_psi)
  # Clamp on logit scale to avoid 0/1 saturation (the C++ likelihood applies
  # the same kind of guard implicitly via log_inv_logit).
  eta_psi <- pmin(pmax(eta_psi, -30), 30)
  log_psi   <- plogis(eta_psi, log.p = TRUE)
  log_1m_psi <- plogis(-eta_psi, log.p = TRUE)

  # Pre-compute per-source linear predictors at the global-site design rows.
  # X_det_s in X_processes is `n_sites x p_det_s` with rows filled at
  # `site_maps[[s]] + 1`; we use only those rows for sites that this source
  # actually observed.
  off <- p_occ
  beta_det_list <- vector("list", n_sources)
  eta_det_full  <- vector("list", n_sources)
  for (s in seq_len(n_sources)) {
    p_det_s <- pi_list[[1 + s]]$p
    beta_det_list[[s]] <- beta[off + seq_len(p_det_s)]
    off <- off + p_det_s
    X_det_s <- model$X_processes[[1 + s]]
    eta_s <- as.numeric(X_det_s %*% beta_det_list[[s]])
    eta_s <- pmin(pmax(eta_s, -30), 30)
    eta_det_full[[s]] <- eta_s
  }

  # Per-site contribution: marginalise z over {0, 1}. We build, for each
  # site i, (a) `any_det_i` flag across all sources, (b) sum_s log P(y_is | z=1)
  # for any-det sites, (c) sum_s log P(y_is = all 0 | z=1) for no-det sites.
  log_det_given_occ_anysite <- numeric(n_sites)   # used when any_det
  log_p_y_z1_allzero        <- numeric(n_sites)   # used when no_det
  any_det                   <- logical(n_sites)

  for (s in seq_len(n_sources)) {
    ys      <- model$y_sources[[s]]
    smap    <- model$site_maps[[s]] + 1L   # 1-indexed global sites
    eta_s   <- eta_det_full[[s]]           # length n_sites; only smap rows used
    n_local <- nrow(ys)
    max_v   <- ncol(ys)
    for (j in seq_len(n_local)) {
      i <- smap[j]
      yij <- ys[j, ]
      valid <- yij >= 0
      if (!any(valid)) next
      log_p_i   <- plogis(eta_s[i],  log.p = TRUE)
      log_1mp_i <- plogis(-eta_s[i], log.p = TRUE)
      n_det_ij  <- sum(yij[valid] == 1L)
      n_zero_ij <- sum(valid) - n_det_ij

      # Conditional on z=1: log Π_t p^y (1-p)^{1-y}
      log_p_y_z1 <- n_det_ij * log_p_i + n_zero_ij * log_1mp_i
      log_det_given_occ_anysite[i] <- log_det_given_occ_anysite[i] + log_p_y_z1

      if (n_det_ij > 0L) {
        any_det[i] <- TRUE
      } else {
        # All zeros at this source: this also enters the no-det branch as
        # Π_t (1 - p)^valid. Accumulate; will be reused below if site i ends
        # up classified as no_det globally.
        log_p_y_z1_allzero[i] <- log_p_y_z1_allzero[i] + log_p_y_z1
      }
    }
  }

  ll <- 0
  for (i in seq_len(n_sites)) {
    if (any_det[i]) {
      ll <- ll + log_psi[i] + log_det_given_occ_anysite[i]
    } else {
      # log( psi * Prod_s q_s + (1 - psi) ); the all-zero conditional
      # log-prob already sums to log(Prod_s q_s).
      log_term1 <- log_psi[i] + log_p_y_z1_allzero[i]
      log_term2 <- log_1m_psi[i]
      mx <- max(log_term1, log_term2)
      ll <- ll + mx + log(exp(log_term1 - mx) + exp(log_term2 - mx))
    }
  }
  ll
}


#' Compute SLA skewness coefficients for an integrated occu fit
#'
#' Builds joint block-diagonal Sigma (Louis I_obs on the occ block, raw
#' H_beta on each det_s block), then calls the generic FD primitive
#' against the original-likelihood evaluator `.loglik_int_occu`.
#'
#' @param model A `tobs_model` (model_type = "integrated").
#' @param em_result EM-Laplace return list (with `$fits$occ`,
#'   `$fits$det1`,…,`$fits$det<S>`, `$weights`).
#' @param spatial Optional `tobs_spatial` (skewness disabled when set —
#'   integrated + spatial Laplace is not yet plumbed in `.tobs_laplace`).
#' @param prior_spec Optional prior spec; passed to `.louis_info_psi_single()`
#'   so the penalty enters Louis I_obs.
#' @return `list(gamma, valid, reason)`.
#' @keywords internal
.sla_compute_int_occu <- function(model, em_result, spatial = NULL,
                                  prior_spec = NULL) {
  if (!is.null(spatial)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "SLA on spatial Sigma not yet implemented for int_occu"))
  }
  if (is.null(em_result$weights)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "em_result$weights missing — needed for Louis Sigma_occ"))
  }

  pi_list   <- model$process_info
  n_sources <- model$n_sources
  n_sites   <- model$n_sites

  fit_occ <- em_result$fits$occ
  if (is.null(fit_occ)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "EM fit missing for occ block"))
  }

  p_occ    <- pi_list[[1]]$p
  beta_psi <- extract_beta(fit_occ, p_occ)

  # Occ block: Louis-corrected observed information (re-use the single-season
  # helper; it computes X' diag(psi(1-psi) - w(1-w)) X with optional prior).
  X_occ <- model$X_processes[[1]]
  I_obs_occ <- .louis_info_psi_single(
    X_occ       = X_occ,
    beta_psi    = beta_psi,
    weights     = em_result$weights,
    spatial     = spatial,
    spatial_fit = fit_occ,
    prior_spec  = prior_spec,
    coef_names  = pi_list[[1]]$coef_names
  )
  Sigma_occ <- tryCatch(solve(I_obs_occ), error = function(e) NULL)
  if (is.null(Sigma_occ)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "Louis I_obs (occ) not invertible"))
  }

  # Per-source detection blocks: invert each fit_det_s$H_beta.
  beta_det_list <- vector("list", n_sources)
  Sigma_det_list <- vector("list", n_sources)
  p_det_vec <- integer(n_sources)
  for (s in seq_len(n_sources)) {
    nm  <- paste0("det", s)
    fi  <- em_result$fits[[nm]]
    if (is.null(fi)) {
      return(list(gamma = NULL, valid = FALSE,
                  reason = sprintf("EM fit missing for %s block", nm)))
    }
    p_det_s <- pi_list[[1 + s]]$p
    p_det_vec[s] <- p_det_s
    beta_det_list[[s]] <- extract_beta(fi, p_det_s)
    H <- fi$H_beta
    if (is.null(H)) {
      return(list(gamma = NULL, valid = FALSE,
                  reason = sprintf("fit_%s$H_beta missing — was return_hessian = TRUE?", nm)))
    }
    Sigma_s <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(Sigma_s)) {
      return(list(gamma = NULL, valid = FALSE,
                  reason = sprintf("fit_%s$H_beta not invertible", nm)))
    }
    # H from the M-step weighted binomial may have rows for the full set of
    # det design columns plus auxiliary blocks (e.g. spatial). Take the
    # leading p_det_s x p_det_s block.
    if (nrow(Sigma_s) > p_det_s) {
      Sigma_s <- Sigma_s[seq_len(p_det_s), seq_len(p_det_s), drop = FALSE]
    }
    Sigma_det_list[[s]] <- Sigma_s
  }

  # Assemble joint block-diagonal Sigma in process_info order:
  # (occ | det1 | det2 | … | detS).
  p_total <- p_occ + sum(p_det_vec)
  Sigma <- matrix(0, p_total, p_total)
  Sigma[seq_len(p_occ), seq_len(p_occ)] <- Sigma_occ
  off <- p_occ
  for (s in seq_len(n_sources)) {
    ps <- p_det_vec[s]
    idx <- off + seq_len(ps)
    Sigma[idx, idx] <- Sigma_det_list[[s]]
    off <- off + ps
  }
  beta_hat <- c(beta_psi, unlist(beta_det_list))

  # FD path against the original-likelihood evaluator.
  log_lik_fn <- function(b) .loglik_int_occu(b, model)
  gamma <- tryCatch(
    .sla_gamma_fd(beta_hat, Sigma, log_lik_fn),
    error = function(e) NULL
  )
  if (is.null(gamma)) {
    return(list(gamma = NULL, valid = FALSE,
                reason = "FD gamma evaluation failed"))
  }

  # Process-info order names: occ_<coef>, <source_name>_<coef>, …
  nms <- paste0(pi_list[[1]]$name, "_", pi_list[[1]]$coef_names)
  for (s in seq_len(n_sources)) {
    nms <- c(nms, paste0(pi_list[[1 + s]]$name, "_", pi_list[[1 + s]]$coef_names))
  }
  names(gamma) <- nms

  list(gamma = gamma, valid = all(is.finite(gamma)),
       reason = if (all(is.finite(gamma))) "ok" else "non-finite gamma")
}
