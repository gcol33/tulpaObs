# ============================================================================
# Occupancy-specific diagnostics
# Generic diagnostics (moran_i, durbin_watson, variogram, compare_models,
# modelAverage) are in tulpa — inherited via tulpa_fit class.
# ============================================================================

#' Model criteria for occupancy / cover models
#'
#' `tobs_waic()`, `tobs_dic()`, and `tobs_cpo()` build the family pointwise
#' log-likelihood matrix once and hand it to the engine's single criteria layer
#' [tulpa::tulpa_criteria()], which derives WAIC / DIC / CPO / LPML / PSIS-LOO.
#' DIC additionally evaluates the deviance at the posterior mean of the
#' parameters, supplied by the family-specific [.tobs_loglik_at_mean()].
#'
#' @param object A `tobs_fit` object.
#' @param n.draws Posterior draws used to build the pointwise log-likelihood
#'   (the cover / occu_cover paths sample this many; the draw-matrix families use
#'   the first `n.draws` stored draws). Default 1000.
#' @param ... Forwarded to [tulpa::tulpa_criteria()] (e.g. `chunk_size`).
#' @return A `tulpa_criteria` object. `tobs_waic()` also carries `$elpd`
#'   (an alias for `elpd_waic`) for back-compatibility.
#' @seealso [tulpa::tulpa_criteria()]
#' @export
tobs_waic <- function(object, n.draws = 1000L, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws)
  cr <- tulpa::tulpa_criteria(ll_mat, criteria = "waic", ...)
  cr$elpd <- cr$elpd_waic
  cr
}

#' @rdname tobs_waic
#' @export
tobs_dic <- function(object, n.draws = 1000L, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws)
  lam <- .tobs_loglik_at_mean(object, n.draws = n.draws)
  tulpa::tulpa_criteria(ll_mat, criteria = "dic", loglik_at_mean = lam, ...)
}

#' @rdname tobs_waic
#' @export
tobs_cpo <- function(object, n.draws = 1000L, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws)
  tulpa::tulpa_criteria(ll_mat, criteria = c("loo", "cpo", "lpml"),
                        pointwise = TRUE, ...)
}

# Pointwise log-likelihood matrix [n_draws x n_obs], marginalized over the
# latent state, for every family. This is the input WAIC, PSIS-LOO, and LOO
# stacking all consume, so it lives in one place (tobs_stack reuses it). The
# per-family marginal likelihoods mirror the C++ kernels in src/*_likelihood.h.
#
# Fidelity note: the linear predictors are evaluated from the *process fixed-
# effect* coefficient draws (the leading columns of `draws`), exactly as
# `fitted()` and the historical WAIC do. Structured-term contributions
# (spatial / temporal / random-effect fields) are not added to the predictor,
# and dynamic visit-level detection covariates are not folded in. For models
# with those components the score is conditional on the fixed-effect predictor.
.tobs_pointwise_loglik <- function(object, n.draws = NULL) {
  nd <- n.draws %||% 1000L
  if (inherits(object, "cover_fit")) return(.tobs_ploglik_cover(object, nd))
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ploglik_occu_cover(object, nd))
  }

  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("Pointwise log-likelihood needs a posterior draw matrix; ",
         "`object$draws` is missing or not a matrix.", call. = FALSE)
  }
  if (!is.null(n.draws) && n.draws < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  .tobs_ploglik_from_draws(object$model, draws)
}

# Per-family pointwise log-likelihood given an explicit [n_draws x p] draw
# matrix. Split out from the dispatcher so the posterior-mean evaluation
# (.tobs_loglik_at_mean) drives the same per-family kernels with a one-row mean.
.tobs_ploglik_from_draws <- function(model, draws) {
  mt <- model$model_type %||% "NULL"
  switch(
    mt,
    single     = .tobs_ploglik_replicated(model, draws),
    dynamic    = .tobs_ploglik_dynamic(model, draws),
    integrated = .tobs_ploglik_integrated(model, draws),
    jsdm       = .tobs_ploglik_jsdm(model, draws),
    stop("Pointwise log-likelihood is not implemented for model_type = '",
         mt, "'.", call. = FALSE)
  )
}

# Pointwise log-likelihood at the posterior mean of the parameters, the plug-in
# DIC needs (length n_obs). The draw-matrix families evaluate their per-family
# kernel at the column-mean draw; the cover / occu_cover families plug in the
# posterior-mean linear predictors (and mean dispersion) via family helpers.
.tobs_loglik_at_mean <- function(object, n.draws = 1000L) {
  if (inherits(object, "cover_fit")) {
    return(.tobs_cover_loglik_at_mean(object, n.draws))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_occu_cover_loglik_at_mean(object, n.draws))
  }
  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("DIC needs a posterior draw matrix to evaluate the deviance at the ",
         "posterior mean; `object$draws` is missing or not a matrix.",
         call. = FALSE)
  }
  mean_draw <- matrix(colMeans(draws), nrow = 1L)
  as.numeric(.tobs_ploglik_from_draws(object$model, mean_draw))
}

# --- shared numerics --------------------------------------------------------

# log(inv_logit(eta)) = log(p) and log(1 - inv_logit(eta)) = log(1 - p),
# computed stably via plogis(log.p) -- the R analogues of the C++
# log_inv_logit / log1m_inv_logit.
.tobs_log_p   <- function(eta) stats::plogis(eta,  log.p = TRUE)
.tobs_log_1mp <- function(eta) stats::plogis(-eta, log.p = TRUE)

# Numerically stable elementwise log(exp(a) + exp(b)).
.tobs_logaddexp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  both_inf <- is.infinite(m) & (a == b)
  out[both_inf] <- m[both_inf]
  out
}

# Cumulative column offset of process k's beta block within the draw matrix
# (processes are laid out in order, fixed-effect betas first).
.tobs_proc_offset <- function(model, k) {
  if (k == 1L) return(0L)
  sum(vapply(model$process_info[seq_len(k - 1L)],
             function(p) p$p, integer(1)))
}

# Linear-predictor draws for process k: [n_draws x nrow(X_k)].
.tobs_eta_draws <- function(model, draws, k) {
  p_k    <- model$process_info[[k]]$p
  off    <- .tobs_proc_offset(model, k)
  beta_k <- draws[, off + seq_len(p_k), drop = FALSE]
  beta_k %*% t(model$X_processes[[k]])
}

# --- per-family marginal likelihoods ---------------------------------------

# Single-season occupancy: per replicate row, marginalized over z.
.tobs_ploglik_replicated <- function(model, draws) {
  eta_psi <- .tobs_eta_draws(model, draws, 1L)   # [S x n_obs]
  eta_p   <- .tobs_eta_draws(model, draws, 2L)
  y       <- model$y                             # [n_obs x max_visits], <0 = NA
  S <- nrow(eta_psi); n_obs <- nrow(y)
  ll <- matrix(0, S, n_obs)
  for (i in seq_len(n_obs)) {
    yi <- y[i, ]; valid <- yi >= 0; nv <- sum(valid)
    if (nv == 0L) next                           # no data -> 0 contribution
    log_p   <- .tobs_log_p(eta_p[, i]);   log_1mp <- .tobs_log_1mp(eta_p[, i])
    log_psi <- .tobs_log_p(eta_psi[, i]); log1m_psi <- .tobs_log_1mp(eta_psi[, i])
    k1 <- sum(yi[valid] == 1); k0 <- nv - k1
    if (k1 > 0L) {
      ll[, i] <- log_psi + k1 * log_p + k0 * log_1mp
    } else {
      ll[, i] <- .tobs_logaddexp(log_psi + nv * log_1mp, log1m_psi)
    }
  }
  ll
}

# JSDM: per site x species, y ~ Bernoulli(psi). No detection, no marginalization.
.tobs_ploglik_jsdm <- function(model, draws) {
  eta <- .tobs_eta_draws(model, draws, 1L)       # [S x N]
  y   <- model$y_jsdm
  Y   <- matrix(y, nrow(eta), length(y), byrow = TRUE)
  Y * .tobs_log_p(eta) + (1 - Y) * .tobs_log_1mp(eta)
}

# Integrated multi-source: per site, shared psi, detection summed over the
# sources that observed it (src/integrated_occ_likelihood.h).
.tobs_ploglik_integrated <- function(model, draws) {
  eta_psi <- .tobs_eta_draws(model, draws, 1L)   # [S x n_sites]
  S <- nrow(eta_psi); n_sites <- model$n_sites
  n_sources <- model$n_sources
  log_psi   <- .tobs_log_p(eta_psi); log1m_psi <- .tobs_log_1mp(eta_psi)
  # per-source log(p) / log(1-p) at every site (valid only where observed)
  log_p_src   <- vector("list", n_sources)
  log_1mp_src <- vector("list", n_sources)
  for (s in seq_len(n_sources)) {
    eta_s <- .tobs_eta_draws(model, draws, 1L + s)
    log_p_src[[s]]   <- .tobs_log_p(eta_s)
    log_1mp_src[[s]] <- .tobs_log_1mp(eta_s)
  }

  ll <- matrix(0, S, n_sites)
  for (i in seq_len(n_sites)) {
    log_det_occ <- numeric(S)   # sum_s log P(y_is | occupied)
    any_det <- FALSE
    for (s in seq_len(n_sources)) {
      local <- which(model$site_maps[[s]] + 1L == i)
      if (!length(local)) next
      yvec <- model$y_sources[[s]][local[1L], ]
      valid <- yvec >= 0; nv <- sum(valid)
      if (nv == 0L) next
      k1 <- sum(yvec[valid] == 1); k0 <- nv - k1
      log_det_occ <- log_det_occ + k1 * log_p_src[[s]][, i] +
                                   k0 * log_1mp_src[[s]][, i]
      if (k1 > 0L) any_det <- TRUE
    }
    if (any_det) {
      ll[, i] <- log_psi[, i] + log_det_occ
    } else {
      # log_det_occ here is sum_s nv_s * log(1-p_s) (all-zero across sources)
      ll[, i] <- .tobs_logaddexp(log_psi[, i] + log_det_occ, log1m_psi[, i])
    }
  }
  ll
}

# Dynamic (multi-season HMM): per-site forward recursion in log space,
# mirroring src/dyn_occ_likelihood.h. Site-level detection only.
.tobs_ploglik_dynamic <- function(model, draws) {
  eta_psi1 <- .tobs_eta_draws(model, draws, 1L)
  eta_p    <- .tobs_eta_draws(model, draws, 2L)
  eta_gam  <- .tobs_eta_draws(model, draws, 3L)
  eta_eps  <- .tobs_eta_draws(model, draws, 4L)
  S <- nrow(eta_psi1); n_sites <- model$n_sites; Tn <- model$n_seasons
  n_visits <- model$n_visits; any_det <- model$any_detected
  NEG <- -1e10

  out <- matrix(0, S, n_sites)
  for (i in seq_len(n_sites)) {
    lp    <- .tobs_log_p(eta_p[, i]);    l1mp   <- .tobs_log_1mp(eta_p[, i])
    lgam  <- .tobs_log_p(eta_gam[, i]);  l1mgam <- .tobs_log_1mp(eta_gam[, i])
    leps  <- .tobs_log_p(eta_eps[, i]);  l1meps <- .tobs_log_1mp(eta_eps[, i])
    a_occ <- .tobs_log_p(eta_psi1[, i]); a_un   <- .tobs_log_1mp(eta_psi1[, i])
    site_ll <- numeric(S)
    for (t in seq_len(Tn)) {
      idx <- (i - 1L) * Tn + (t - 1L) + 1L
      nv  <- n_visits[idx]
      if (nv > 0L) {
        yvec <- model$y[i, , t]; yvec[is.na(yvec)] <- -1L
        valid <- yvec >= 0
        k1 <- sum(yvec[valid] == 1); k0 <- sum(yvec[valid] == 0)
        if (any_det[idx]) {
          site_ll <- site_ll + a_occ + (k1 * lp + k0 * l1mp)
          a_occ <- numeric(S); a_un <- rep(NEG, S)   # z_t = 1 known
        } else {
          term1 <- a_occ + nv * l1mp                 # occupied, all non-detections
          term2 <- a_un                              # unoccupied
          lnorm <- .tobs_logaddexp(term1, term2)
          site_ll <- site_ll + lnorm
          a_occ <- term1 - lnorm
          a_un  <- term2 - lnorm
        }
      }
      if (t < Tn) {
        new_occ <- .tobs_logaddexp(a_occ + l1meps, a_un + lgam)
        new_un  <- .tobs_logaddexp(a_occ + leps,   a_un + l1mgam)
        a_occ <- new_occ; a_un <- new_un
      }
    }
    out[, i] <- site_ll
  }
  out
}

#' Posterior predictive check
#' @param object A `tobs_fit` object.
#' @param fit.stat `"freeman-tukey"` (default) or `"chi-squared"`.
#' @param n.samples Number of posterior samples (default 500).
#' @return A list with `fit.y`, `fit.y.rep`, and `bayesian.p`.
#' @export
tobs_ppc <- function(object, fit.stat = c("freeman-tukey", "chi-squared"),
                     n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  if (inherits(object, "cover_fit")) {
    return(.tobs_ppc_cover(object, fit.stat, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ppc_occu_cover(object, fit.stat, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("tobs_ppc supports single-season occupancy, cover(), and occu_cover() ",
         "fits.", call. = FALSE)
  }

  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y; n_sites <- model$n_sites; max_visits <- model$max_visits
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n.samples <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n.samples)

  stat_fn <- if (fit.stat == "freeman-tukey") {
    function(obs, exp) sum((sqrt(obs) - sqrt(exp))^2, na.rm = TRUE)
  } else {
    function(obs, exp) sum((obs - exp)^2 / (exp + 1e-10), na.rm = TRUE)
  }

  # Per-site valid mask, visit count, and whether the species was ever
  # detected. The latent z is sampled from its full conditional given this
  # detection history (the spOccupancy ppcOcc construction), not from the
  # prior psi: a site with any detection is occupied with probability 1, an
  # all-zero history occupied with probability
  #   psi (1-p)^J / [psi (1-p)^J + (1-psi)].
  valid_mat <- y >= 0
  n_valid   <- rowSums(valid_mat)
  any_det   <- rowSums(y * valid_mat) > 0

  fit_y <- fit_y_rep <- numeric(n.samples)
  for (s in seq_len(n.samples)) {
    idx <- draw_idx[s]
    psi <- plogis(as.vector(X_occ %*% draws[idx, seq_len(p_occ)]))
    p   <- plogis(as.vector(X_det %*% draws[idx, p_occ + seq_len(p_det)]))
    z_prob <- ifelse(
      any_det, 1,
      psi * (1 - p)^n_valid / (psi * (1 - p)^n_valid + (1 - psi))
    )
    z_prob[n_valid == 0L] <- psi[n_valid == 0L]   # no data -> prior
    z <- rbinom(n_sites, 1, z_prob)
    expected <- y_rep <- matrix(NA, n_sites, max_visits)
    for (i in seq_len(n_sites)) for (j in seq_len(max_visits)) if (y[i,j] >= 0) {
      expected[i,j] <- z[i] * p[i]; y_rep[i,j] <- rbinom(1, 1, z[i] * p[i])
    }
    valid <- !is.na(expected)
    fit_y[s] <- stat_fn(y[valid], expected[valid])
    fit_y_rep[s] <- stat_fn(y_rep[valid], expected[valid])
  }
  list(fit.y = fit_y, fit.y.rep = fit_y_rep, bayesian.p = mean(fit_y_rep > fit_y))
}

#' PIT residuals
#' @param object A `tobs_fit` object.
#' @param n.samples Number of posterior samples (default 250).
#' @return Numeric vector of PIT residuals.
#' @export
tobs_pit_residuals <- function(object, n.samples = 250) {
  if (inherits(object, "cover_fit")) {
    return(.tobs_pit_cover(object, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_pit_occu_cover(object, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("tobs_pit_residuals supports single-season occupancy, cover(), and ",
         "occu_cover() fits.", call. = FALSE)
  }
  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y; n_sites <- model$n_sites
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n_draws <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n_draws)

  pit <- numeric(n_sites)
  for (i in seq_len(n_sites)) {
    yi <- y[i, ]; valid <- yi >= 0; n_det <- sum(yi[valid] == 1); n_valid <- sum(valid)
    if (n_valid == 0) { pit[i] <- runif(1); next }
    cdf_vals <- numeric(n_draws)
    for (s in seq_len(n_draws)) {
      psi <- plogis(sum(X_occ[i, ] * draws[draw_idx[s], seq_len(p_occ)]))
      p <- plogis(sum(X_det[i, ] * draws[draw_idx[s], p_occ + seq_len(p_det)]))
      cdf_vals[s] <- if (n_det > 0) 1 else psi * (1 - p)^n_valid + (1 - psi)
    }
    pit[i] <- min(1, max(0, mean(cdf_vals) + runif(1, 0, 1/n_draws)))
  }
  pit
}

#' @export
tobs_test_uniformity <- function(pit) {
  # PIT residuals can have ties when the response is discrete (counts /
  # detections). The asymptotic KS p-value remains valid; only the
  # ties-warning is noisy, so suppress just that warning class.
  withCallingHandlers(
    ks.test(pit, "punif"),
    warning = function(w) {
      if (grepl("ties should not be present", conditionMessage(w),
                fixed = TRUE)) invokeRestart("muffleWarning")
    }
  )
}

#' @export
tobs_test_dispersion <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  obs_var <- var(rowSums(y_obs * (y_obs >= 0), na.rm = TRUE))
  sim_vars <- vapply(sims, function(ys) var(rowSums(ys * (ys >= 0), na.rm = TRUE)), double(1))
  list(observed = obs_var, expected = mean(sim_vars),
       ratio = obs_var / mean(sim_vars), p.value = mean(sim_vars >= obs_var))
}

#' @export
tobs_test_zero_inflation <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  count_zeros <- function(y) sum(apply(y, 1, function(r) { v <- r >= 0; all(r[v] == 0) }))
  obs <- count_zeros(y_obs); sim <- vapply(sims, count_zeros, integer(1))
  list(observed = obs, expected = mean(sim), ratio = obs / max(mean(sim), 1),
       p.value = mean(sim >= obs))
}

#' @export
tobs_test_outliers <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  n_sites <- nrow(y_obs)
  obs_det <- apply(y_obs, 1, function(r) sum(r[r >= 0] == 1))
  sim_det <- vapply(sims, function(ys) apply(ys, 1, function(r) sum(r[r >= 0] == 1)), double(n_sites))
  lower <- apply(sim_det, 1, quantile, 0.025); upper <- apply(sim_det, 1, quantile, 0.975)
  n_out <- sum(obs_det < lower | obs_det > upper)
  sim_out <- vapply(seq_len(n.samples), function(s) sum(sim_det[,s] < lower | sim_det[,s] > upper), integer(1))
  list(n_outliers = n_out, expected = mean(sim_out), p.value = mean(sim_out >= n_out))
}

#' Comprehensive model checking
#' @param object A `tobs_fit` object.
#' @param coords Optional coordinates for spatial diagnostics.
#' @param n.samples Posterior samples for simulation tests.
#' @return Invisibly, diagnostic results.
#' @export
tobs_check <- function(object, coords = NULL, n.samples = 250) {
  cat("=== tobs Model Diagnostics ===\n\n")
  if (identical(object$method, "nuts")) {
    cat(sprintf("Sampler: %d samples, %d divergent, mean accept = %.3f\n",
                object$n_samples, sum(object$divergent), mean(object$accept_prob)))
  } else {
    cat(sprintf("Fit: %s, %d posterior draws (no NUTS sampler diagnostics)\n",
                object$method %||% "laplace", object$n_samples))
  }

  w <- tryCatch(tobs_waic(object), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("\nWAIC: %.1f (p_waic = %.1f)\n", w$waic, w$p_waic))

  dic <- tryCatch(tobs_dic(object), error = function(e) NULL)
  if (!is.null(dic) && is.finite(dic$dic))
    cat(sprintf("DIC:  %.1f (p_DIC = %.1f)\n", dic$dic, dic$p_dic))

  cpo <- tryCatch(tobs_cpo(object), error = function(e) NULL)
  if (!is.null(cpo)) cat(sprintf("LPML: %.1f (elpd_loo = %.1f)\n",
                                 cpo$lpml, cpo$elpd_loo))

  ppc <- tryCatch(tobs_ppc(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(ppc)) {
    cat(sprintf("\nPPC: Bayesian p = %.3f\n", ppc$bayesian.p))
    if (ppc$bayesian.p < 0.05 || ppc$bayesian.p > 0.95) cat("  WARNING: poor fit\n")
  }

  zi <- tryCatch(tobs_test_zero_inflation(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(zi)) cat(sprintf("\nZero-inflation: obs=%d, exp=%.1f, p=%.3f\n", zi$observed, zi$expected, zi$p.value))

  if (!is.null(coords)) {
    mi <- tryCatch(tulpa::moran_i(residuals(object)$occ, coords), error = function(e) NULL)
    if (!is.null(mi)) {
      cat(sprintf("\nMoran's I: %.3f (p = %.3f)\n", unname(mi$statistic), mi$p.value))
      if (mi$p.value < 0.05) cat("  WARNING: spatial autocorrelation\n")
    }
  }
  cat("\n")
  invisible(list(waic = w, dic = dic, cpo = cpo, ppc = ppc,
                 zero_inflation = zi))
}
