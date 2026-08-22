# =============================================================================
# Exact-marginal refinement for dynamic (multi-season HMM) occupancy.
#
# Like the single-season occupancy marginal (R/occu_marginal.R), the dynamic
# occupancy marginal is exact and cheap: the latent occupancy sequence
# integrates out by the HMM forward recursion -- `.tobs_ploglik_dynamic()`, the
# same kernel logLik() / WAIC / LOO score. The default deterministic engine fits
# this with a forward-backward E-step and a per-block pseudo-binomial Laplace
# M-step (the colonization / extinction blocks average over M = 1000 pseudo
# transitions). That encoding leaves a small EM-discretisation residual below
# the exact marginal MLE and reads the wrong curvature for the marginal, so the
# block standard errors are mis-scaled.
#
# Per tulpa's "nested approximation + debias" design, this is the debias step: a
# Newton refinement of the EM mode on the exact HMM-forward marginal
# log-posterior, with standard errors read from its Hessian. It moves the fixed
# effects onto the exact marginal MLE (matching unmarked::colext) and restores
# calibrated, near-nominal-coverage SEs. The refinement shares the
# optimise / vcov / draw-update / log-likelihood-refresh scaffold with the
# single-season path (`.tobs_marginal_refine_apply`); only the marginal
# objective differs (HMM forward vs the closed-form two-state mixture).
#
# Spatial (nested-Laplace) dynamic fits keep the EM field -- they carry a latent
# field, not a closed-form coefficient marginal -- so the refine is gated to the
# plain (non-spatial) Laplace path in `.tobs_laplace()`. Dynamic detection is
# site-level (the dyn_occu builder has no visit-level detection design), matching
# the HMM-forward kernel.
# =============================================================================


# Exact-marginal negative log-posterior for dynamic occupancy at the packed
# fixed-effect vector par = c(beta_psi1, beta_p, beta_gamma, beta_epsilon), in
# the process-block order `.tobs_ploglik_dynamic()` indexes. `pmean` / `pprec`
# are the Gaussian-prior mean and precision (1/sd^2, 0 = flat) aligned with par,
# so the refine maximises the same penalised marginal the EM M-steps target.
.tobs_dyn_occu_marginal_nlp <- function(par, model, pmean, pprec) {
  draw <- matrix(par, nrow = 1L)
  ll   <- .tobs_ploglik_dynamic(model, draw)   # [1 x n_sites]
  val  <- sum(ll[1L, ])
  if (!is.finite(val)) return(1e10)
  -val + 0.5 * sum(pprec * (par - pmean)^2)
}


# Per-(site, season) detection sufficient statistics for the R HMM-forward
# marginal: n_valid (visits with an observed y), n_det (sum of detections), and
# any-detected. Cached on first use; a season with no valid visit contributes a
# neutral emission (both states 1).
.tobs_dyn_occu_emit_stats <- function(model) {
  if (!is.null(model$.emit_stats)) return(model$.emit_stats)
  y <- model$y                                  # [n_sites x mv x T]
  n_sites <- dim(y)[1]; T_s <- dim(y)[3]
  nvalid <- matrix(0L, n_sites, T_s); ndet <- matrix(0L, n_sites, T_s)
  for (t in seq_len(T_s)) {
    yt <- y[, , t]
    v  <- !is.na(yt)
    nvalid[, t] <- rowSums(v)
    yt0 <- yt; yt0[!v] <- 0L
    ndet[, t] <- rowSums(yt0)
  }
  list(nvalid = nvalid, ndet = ndet)
}


# Exact HMM-forward marginal log-likelihood for dynamic occupancy with
# SEASON-VARYING colonization / extinction. The compiled cpp_occu_dynamic_ploglik
# reads one gamma / epsilon per site, so the interval-indexed refine needs its
# own forward pass. Per site the latent occupancy sequence integrates out by the
# 2-state forward recursion whose transition at interval t uses interval-t rates;
# the emission is the per-season detection likelihood (occupied: prod p^y
# (1-p)^(1-y); empty: 1 if no detection, 0 otherwise). Returns the total
# penalised negative log-posterior. `par` is the packed c(beta_psi1, beta_p,
# beta_gamma, beta_epsilon).
.tobs_dyn_occu_marginal_nlp_sv <- function(par, model, pmean, pprec, p_sizes) {
  n_sites <- model$n_sites; T_s <- model$n_seasons; n_int <- T_s - 1L
  o <- 0L
  b_psi <- par[(o + 1L):(o + p_sizes[1])]; o <- o + p_sizes[1]
  b_p   <- par[(o + 1L):(o + p_sizes[2])]; o <- o + p_sizes[2]
  b_gam <- par[(o + 1L):(o + p_sizes[3])]; o <- o + p_sizes[3]
  b_eps <- par[(o + 1L):(o + p_sizes[4])]

  eta_psi <- as.numeric(model$X_processes[[1]] %*% b_psi)     # [n_sites]
  # Detection eta [n_sites x T]: season-varying reads the long-form design
  # site-major season-minor; a constant arm's per-site eta broadcasts across the
  # seasons (byte-identical to the pre-#124 path).
  ep <- as.numeric(model$X_processes[[2]] %*% b_p)
  eta_p_mat <- if (isTRUE(model$det_season_varying))
    matrix(ep, n_sites, T_s, byrow = TRUE) else matrix(ep, n_sites, T_s)
  # gamma / epsilon per interval: long-form design -> [n_sites x n_int] byrow;
  # a constant arm's site-level eta recycles across the intervals.
  eg <- as.numeric(model$X_processes[[3]] %*% b_gam)
  ee <- as.numeric(model$X_processes[[4]] %*% b_eps)
  gam_mat <- if (isTRUE(model$col_season_varying))
    matrix(eg, n_sites, n_int, byrow = TRUE) else matrix(eg, n_sites, n_int)
  eps_mat <- if (isTRUE(model$ext_season_varying))
    matrix(ee, n_sites, n_int, byrow = TRUE) else matrix(ee, n_sites, n_int)

  lg  <- plogis(eta_psi, log.p = TRUE)          # log psi1
  l1g <- plogis(-eta_psi, log.p = TRUE)         # log(1 - psi1)
  lp_mat  <- plogis(eta_p_mat, log.p = TRUE)    # log p        [n_sites x T]
  l1p_mat <- plogis(-eta_p_mat, log.p = TRUE)   # log(1 - p)   [n_sites x T]
  # Transition log-probs per interval [n_sites x n_int]
  lgam  <- plogis(gam_mat, log.p = TRUE);  l1gam <- plogis(-gam_mat, log.p = TRUE)
  leps  <- plogis(eps_mat, log.p = TRUE);  l1eps <- plogis(-eps_mat, log.p = TRUE)

  st <- .tobs_dyn_occu_emit_stats(model)
  nvalid <- st$nvalid; ndet <- st$ndet
  lse2 <- function(a, b) {                       # log_sum_exp of two vectors
    m <- pmax(a, b); m + log(exp(a - m) + exp(b - m))
  }
  # Emission log-prob per (site, season): occupied uses the per-visit detection
  # likelihood; empty is 0 (log 1) when no visit detected, -Inf otherwise.
  emit1 <- ndet * lp_mat + (nvalid - ndet) * l1p_mat            # [n_sites x T]
  emit0 <- ifelse(ndet > 0, -Inf, 0)                            # [n_sites x T]

  # Forward recursion, vectorised over sites.
  la1 <- lg  + emit1[, 1]
  la0 <- l1g + emit0[, 1]
  for (t in 2:T_s) {
    iv <- t - 1L
    n1 <- lse2(la1 + l1eps[, iv], la0 + lgam[, iv]) + emit1[, t]
    n0 <- lse2(la1 + leps[, iv],  la0 + l1gam[, iv]) + emit0[, t]
    la1 <- n1; la0 <- n0
  }
  ll_site <- lse2(la1, la0)
  val <- sum(ll_site[is.finite(ll_site)])
  if (!is.finite(val)) return(1e10)
  -val + 0.5 * sum(pprec * (par - pmean)^2)
}


# Refine a dynamic-occupancy `tobs_fit` on the exact HMM-forward marginal and
# overwrite its fixed-effect estimates + SEs + pseudo-draws with the debiased,
# Hessian-calibrated values. Falls back to the unmodified fit on any failure
# (never worse than the EM result).
.tobs_dyn_occu_marginal_refine <- function(fit, model, prior_spec = NULL) {
  tryCatch({
    pi_list <- model$process_info
    if (length(pi_list) < 4L) return(fit)

    sv_dyn <- isTRUE(model$col_season_varying) ||
              isTRUE(model$ext_season_varying) ||
              isTRUE(model$det_season_varying)

    # Fixed-effect names in process-block order (psi1, p, gamma, epsilon), the
    # layout fit$means / the draw matrix already carry.
    all_nm <- unlist(lapply(pi_list, function(p)
      paste0(p$name, "_", p$coef_names)), use.names = FALSE)

    # Gaussian-prior mean / precision per process block, aligned with par. The
    # bucket keying (psi1 / p / gamma / epsilon) matches `.attach_priors_to_blocks`
    # so the refine penalises exactly as the EM did; flat when prior_spec is NULL
    # (priors = FALSE), giving the exact marginal MLE.
    pr    <- lapply(pi_list, function(p)
      .prior_for_submodel(prior_spec, p$name, p$coef_names))
    pmean <- unlist(lapply(pr, `[[`, "mean"), use.names = FALSE)
    psd   <- unlist(lapply(pr, `[[`, "sd"),   use.names = FALSE)
    pprec <- ifelse(is.finite(psd) & psd > 0, 1 / psd^2, 0)

    # Season-varying transitions need the R HMM-forward marginal (the compiled
    # cpp kernel reads one gamma / epsilon per site); the constant-rate path uses
    # the fast cpp marginal. Either way the exact-marginal Newton refine moves the
    # fixed effects onto the marginal MLE (escaping any EM local optimum) and
    # reads calibrated SEs from the Hessian.
    if (sv_dyn) {
      p_sizes <- vapply(pi_list, function(p) as.integer(p$p), integer(1))
      model$.emit_stats <- .tobs_dyn_occu_emit_stats(model)  # cache once
      nlp <- function(par)
        .tobs_dyn_occu_marginal_nlp_sv(par, model, pmean, pprec, p_sizes)
    } else {
      nlp <- function(par) .tobs_dyn_occu_marginal_nlp(par, model, pmean, pprec)
    }
    .tobs_marginal_refine_apply(fit, model, all_nm, nlp)
  }, error = function(e) fit)
}
