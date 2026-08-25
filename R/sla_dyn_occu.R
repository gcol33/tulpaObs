# ============================================================================
# sla_dyn_occu.R — Simplified Laplace for the dynamic (HMM) occupancy family
#
# The single-season occu likelihood decomposes as sum_i l_i(eta_i) with
# eta_i = X_i beta — the analytical SLA path in
# `.sla_compute_occu_single()` rests on that. The HMM forward does *not*
# decompose that way: alpha_t depends recursively on alpha_{t-1}, so the
# third derivative in eta-space has non-trivial off-diagonal terms across the
# four process blocks (psi1, p, gamma, epsilon). We therefore drive the
# generic finite-difference d3 primitive `.sla_gamma_fd()` against an
# R-side HMM log-likelihood evaluator `.loglik_dyn_occu()` that mirrors
# the C++ forward in `src/dyn_occ_likelihood.h`.
#
# The joint posterior covariance Sigma used by the FD path is assembled
# block-diagonally:
#
#   psi1 block (init): Louis-corrected observed Fisher info
#                      I_obs = X_psi' diag(psi1(1-psi1) - w1(1-w1)) X_psi
#                      (the M-step pseudo-binomial M-inflation is removed
#                       analytically, matching `.louis_info_psi_single()`).
#
#   p block (det):     fit_det$H_beta directly — the M-step encoding for
#                      detection is a *weighted* binomial (no pseudo trials),
#                      so the inner Hessian already equals the observed info.
#
#   gamma block (col): fit_col$H_beta / M.  The col M-step is pseudo-binomial
#                      with M = 1000 trials per averaged transition; dividing
#                      by M strips the artefactual inflation and leaves the
#                      effective-sample-size Fisher info. This is a *complete-
#                      data* approximation (the missing-z correction would
#                      add cross-time covariance terms that are O(1/T) for
#                      moderate T and not worth the algebra here).
#
#   epsilon block (ext): same as col — fit_ext$H_beta / M.
#
# Cross-block off-diagonal terms in Sigma are zero by construction. The
# diagonal-Sigma assumption is consistent with how SLA marginals are used
# elsewhere in tulpaObs: only per-coefficient marginals are corrected, joint
# correlations are left to Gaussian Laplace.
# ============================================================================


# Inner-M-step pseudo-binomial trial count for the col/ext blocks (must match
# `build_dynamic_callbacks()::m_step_encode`).
.SLA_DYN_M_PSEUDO <- 1000L


# ---------------------------------------------------------------------------
# R-side HMM forward log-likelihood for dyn_occu
# ---------------------------------------------------------------------------

#' R-side multi-season HMM log-likelihood for dyn_occu
#'
#' Computes the marginal observation log-likelihood for the dynamic occupancy
#' model under the C++ forward-algorithm semantics in
#' `src/dyn_occ_likelihood.h`. Used by `.sla_gamma_fd()` to evaluate the
#' third-cumulant correction; not exported.
#'
#' Joint coefficient vector layout (matches `model$process_info`):
#'   beta = c(beta_psi1, beta_p, beta_gamma, beta_epsilon)
#'
#' Site-level or season-varying detection: the `dyn_occu()` detection design is
#' site-level, optionally season-indexed (a `[n_sites x n_seasons]` covariate,
#' #124); there is no visit-level (within-season) detection covariate.
#'
#' No priors. No pseudo-binomial encoding. Returns the original-data
#' log-likelihood, matching the convention of `.loglik_occu_single()`.
#'
#' @param beta Numeric joint coefficient vector, length p_psi+p_p+p_col+p_ext.
#' @param model A `tobs_model` of `model_type = "dynamic"`.
#' @return Scalar log-likelihood.
#' @keywords internal
.loglik_dyn_occu <- function(beta, model) {
  pi_list <- model$process_info
  p_psi <- pi_list[[1]]$p
  p_p   <- pi_list[[2]]$p
  p_col <- pi_list[[3]]$p
  p_ext <- pi_list[[4]]$p

  beta_psi <- beta[seq_len(p_psi)]
  beta_p   <- beta[p_psi + seq_len(p_p)]
  beta_col <- beta[p_psi + p_p + seq_len(p_col)]
  beta_ext <- beta[p_psi + p_p + p_col + seq_len(p_ext)]

  X_psi <- model$X_processes[[1]]
  X_p   <- model$X_processes[[2]]
  X_col <- model$X_processes[[3]]
  X_ext <- model$X_processes[[4]]

  n_sites    <- model$n_sites
  n_seasons  <- model$n_seasons
  max_visits <- model$max_visits

  # Read the 3D y array directly (layout-agnostic). This is the canonical
  # user-supplied form, stored as model$y by .tobs_build_dynamic. The flat
  # y_flat layout (site-major as of the indexing fix) is for the C++ HMM
  # and the M-step encoder; we don't touch it here.
  y_arr <- model$y
  nv    <- model$n_visits         # length n_sites * n_seasons
  ad    <- model$any_detected     # length n_sites * n_seasons

  psi1 <- plogis(as.numeric(X_psi %*% beta_psi))
  p_v  <- plogis(as.numeric(X_p   %*% beta_p))
  gam  <- plogis(as.numeric(X_col %*% beta_col))
  eps  <- plogis(as.numeric(X_ext %*% beta_ext))

  # Pre-compute log(p) and log(1-p) per site (uniform across visits & seasons)
  log_p   <- log(pmax(p_v,     .Machine$double.eps))
  log1m_p <- log(pmax(1 - p_v, .Machine$double.eps))

  # For numerical safety, work in log space and mirror the forward algorithm.
  NEG_INF <- -1e10
  log_sum_exp2 <- function(a, b) {
    m <- max(a, b)
    if (!is.finite(m)) return(NEG_INF)
    m + log(exp(a - m) + exp(b - m))
  }

  log_lik <- 0
  for (i in seq_len(n_sites)) {
    log_a_occ   <- log(max(psi1[i],     .Machine$double.eps))
    log_a_unocc <- log(max(1 - psi1[i], .Machine$double.eps))

    log_gam_i      <- log(max(gam[i],     .Machine$double.eps))
    log1m_gam_i    <- log(max(1 - gam[i], .Machine$double.eps))
    log_eps_i      <- log(max(eps[i],     .Machine$double.eps))
    log1m_eps_i    <- log(max(1 - eps[i], .Machine$double.eps))

    for (t in seq_len(n_seasons)) {
      idx <- (i - 1) * n_seasons + t
      nv_it <- nv[idx]
      det_this <- ad[idx]

      if (nv_it > 0L) {
        log_p_data_occ <- 0
        sum_log1m_p <- 0
        for (j in seq_len(max_visits)) {
          y_ij <- y_arr[i, j, t]
          if (is.na(y_ij) || y_ij < 0L) next
          if (y_ij == 1L) log_p_data_occ <- log_p_data_occ + log_p[i]
          else            log_p_data_occ <- log_p_data_occ + log1m_p[i]
          sum_log1m_p <- sum_log1m_p + log1m_p[i]
        }

        if (det_this) {
          # P(data | unocc) = 0 — only the occupied path contributes
          log_lik <- log_lik + log_a_occ + log_p_data_occ
          # After detection, z_t = 1 is known: reset alpha
          log_a_occ   <- 0
          log_a_unocc <- NEG_INF
        } else {
          # P(data | occ) = prod(1-p), P(data | unocc) = 1
          term1 <- log_a_occ   + sum_log1m_p
          term2 <- log_a_unocc
          log_norm <- log_sum_exp2(term1, term2)
          log_lik <- log_lik + log_norm
          # Renormalise to a proper conditional alpha so next-season
          # forward step compounds correctly (matches C++).
          log_a_occ   <- term1 - log_norm
          log_a_unocc <- term2 - log_norm
        }
      }
      # nv_it == 0 -> no data, alpha unchanged, no lik contribution

      # Transition to t+1
      if (t < n_seasons) {
        new_log_occ <- log_sum_exp2(log_a_occ   + log1m_eps_i,
                                    log_a_unocc + log_gam_i)
        new_log_unocc <- log_sum_exp2(log_a_occ   + log_eps_i,
                                      log_a_unocc + log1m_gam_i)
        log_a_occ   <- new_log_occ
        log_a_unocc <- new_log_unocc
      }
    }
  }

  log_lik
}


# ---------------------------------------------------------------------------
# Per-block Sigma assembly
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

#' Compute SLA skewness coefficients for a dyn_occu (HMM) fit
#'
#' Drives the generic finite-difference d3 primitive against the R-side HMM
#' forward log-likelihood. The joint Sigma is assembled block-diagonally
#' with per-block precision matrices (see file-header note for the
#' justification).
#'
#' @param model A `tobs_model` (model_type = "dynamic").
#' @param em_result The EM-Laplace return list with `$fits$occ/det/col/ext`
#'   and `$weights` (n_sites x n_seasons matrix).
#' @param spatial Optional `tobs_spatial`. A state-arm field reaches this
#'   evaluator (`.validate_spatial_laplace()` blocks only a detection-arm SPDE
#'   for the dynamic model), and the correction is intentionally NOT applied
#'   under one -- see `.sla_spatial_reason()`, the same decision the single-
#'   season and integrated paths take. Returns `valid = FALSE`.
#' @param prior_spec Optional prior spec; passed to the init-block Louis
#'   info so the penalty is included.
#' @return List(gamma, valid, reason).
#' @keywords internal
.sla_compute_dyn_occu <- function(model, em_result, spatial = NULL,
                                  prior_spec = NULL) {
  bail <- .sla_bailer("gamma")
  if (!is.null(spatial)) {
    return(bail(.sla_spatial_reason("dyn_occu")))
  }
  if (is.null(em_result$weights)) {
    return(bail("em_result$weights missing (need n_sites x n_seasons matrix)"))
  }

  fits <- em_result$fits
  fit_occ <- fits$occ
  fit_det <- fits$det
  fit_col <- fits$col
  fit_ext <- fits$ext
  if (is.null(fit_occ) || is.null(fit_det) ||
      is.null(fit_col) || is.null(fit_ext)) {
    return(bail("EM fits missing for one of occ/det/col/ext blocks"))
  }
  for (nm in c("det", "col", "ext")) {
    if (is.null(fits[[nm]]$H_beta)) {
      return(bail(sprintf("fit_%s$H_beta missing (was return_hessian = TRUE?)", nm)))
    }
  }

  pi_list <- model$process_info
  p_psi <- pi_list[[1]]$p
  p_p   <- pi_list[[2]]$p
  p_col <- pi_list[[3]]$p
  p_ext <- pi_list[[4]]$p
  p_total <- p_psi + p_p + p_col + p_ext

  X_psi <- model$X_processes[[1]]

  beta_psi <- as.numeric(fit_occ$mode)[seq_len(p_psi)]
  beta_p   <- as.numeric(fit_det$mode)[seq_len(p_p)]
  beta_col <- as.numeric(fit_col$mode)[seq_len(p_col)]
  beta_ext <- as.numeric(fit_ext$mode)[seq_len(p_ext)]
  beta_hat <- c(beta_psi, beta_p, beta_col, beta_ext)

  # ---- psi1 block: Louis observed info ------------------------------------
  w <- em_result$weights
  if (!is.matrix(w) || nrow(w) != model$n_sites || ncol(w) != model$n_seasons) {
    return(bail("weights shape != n_sites x n_seasons"))
  }
  # The season-1 E-step posterior P(z_{i,1} = 1 | y) is the state weight the
  # Louis identity takes; the initial-occupancy prior keys on "psi1".
  I_psi <- .louis_info_psi_single(
    X_occ      = X_psi,
    beta_psi   = beta_psi,
    weights    = w[, 1],
    prior_spec = prior_spec,
    coef_names = pi_list[[1]]$coef_names,
    prior_arms = c("psi1", "psi")
  )
  Sigma_psi <- .sla_solve(I_psi)
  if (is.null(Sigma_psi)) return(bail("Louis I_obs (psi1) not invertible"))

  # ---- det block: raw H_beta (weighted binomial, no M inflation) ----------
  Sigma_p <- .sla_solve(as.matrix(fit_det$H_beta))
  if (is.null(Sigma_p)) return(bail("fit_det$H_beta not invertible"))

  # ---- col/ext blocks: H_beta / M (strip pseudo-binomial inflation) -------
  M <- .SLA_DYN_M_PSEUDO
  Sigma_col <- .sla_solve(as.matrix(fit_col$H_beta) / M)
  if (is.null(Sigma_col)) return(bail("fit_col$H_beta/M not invertible"))
  Sigma_ext <- .sla_solve(as.matrix(fit_ext$H_beta) / M)
  if (is.null(Sigma_ext)) return(bail("fit_ext$H_beta/M not invertible"))

  # ---- Block-diagonal joint Sigma -----------------------------------------
  Sigma <- matrix(0, p_total, p_total)
  i1 <- seq_len(p_psi)
  i2 <- p_psi + seq_len(p_p)
  i3 <- p_psi + p_p + seq_len(p_col)
  i4 <- p_psi + p_p + p_col + seq_len(p_ext)
  Sigma[i1, i1] <- Sigma_psi
  Sigma[i2, i2] <- Sigma_p
  Sigma[i3, i3] <- Sigma_col
  Sigma[i4, i4] <- Sigma_ext

  # ---- Finite-difference gamma against the HMM log-likelihood --------------
  log_lik_fn <- function(b) .loglik_dyn_occu(b, model)
  gamma <- .sla_gamma_fd(beta_hat = beta_hat, Sigma = Sigma,
                         log_lik_fn = log_lik_fn)

  nms <- c(
    paste0(pi_list[[1]]$name, "_", pi_list[[1]]$coef_names),
    paste0(pi_list[[2]]$name, "_", pi_list[[2]]$coef_names),
    paste0(pi_list[[3]]$name, "_", pi_list[[3]]$coef_names),
    paste0(pi_list[[4]]$name, "_", pi_list[[4]]$coef_names)
  )
  names(gamma) <- nms

  list(gamma = gamma, valid = all(is.finite(gamma)),
       reason = if (all(is.finite(gamma))) "ok" else "non-finite gamma")
}
