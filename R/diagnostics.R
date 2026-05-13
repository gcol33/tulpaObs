# ============================================================================
# Occupancy-specific diagnostics
# Generic diagnostics (moranI, durbinWatson, variogram, compare_models,
# modelAverage) are in tulpa — inherited via tulpa_fit class.
# ============================================================================

#' Compute WAIC for occupancy models
#'
#' @param object A `tobs_fit` object.
#' @param ... Ignored.
#' @return A list with `waic`, `elpd`, `p_waic`, and pointwise values.
#' @export
waicOccu <- function(object, ...) {
  model <- object$model
  draws <- object$draws
  pi_list <- model$process_info
  n_draws <- nrow(draws)

  if (model$model_type != "single") {
    stop("waicOccu currently supports single-season models only")
  }

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  y <- model$y
  n_sites <- model$n_sites
  p_occ <- pi_list[[1]]$p
  p_det <- pi_list[[2]]$p

  ll_mat <- matrix(NA_real_, n_draws, n_sites)
  for (s in seq_len(n_draws)) {
    beta_occ <- draws[s, seq_len(p_occ)]
    beta_det <- draws[s, p_occ + seq_len(p_det)]
    psi <- plogis(as.vector(X_occ %*% beta_occ))
    p <- plogis(as.vector(X_det %*% beta_det))
    for (i in seq_len(n_sites)) {
      yi <- y[i, ]; valid <- yi >= 0; n_valid <- sum(valid)
      if (n_valid == 0) { ll_mat[s, i] <- 0; next }
      if (any(yi[valid] == 1)) {
        ll_mat[s, i] <- log(psi[i]) + sum(yi[valid] * log(p[i]) + (1 - yi[valid]) * log(1 - p[i]))
      } else {
        ll_mat[s, i] <- log(psi[i] * (1 - p[i])^n_valid + (1 - psi[i]))
      }
    }
  }

  lppd <- sum(log(colMeans(exp(ll_mat))))
  p_waic <- sum(apply(ll_mat, 2, var))
  list(waic = -2 * (lppd - p_waic), elpd = lppd - p_waic,
       p_waic = p_waic, lppd = lppd)
}

#' Posterior predictive check
#' @param object A `tobs_fit` object.
#' @param fit.stat `"freeman-tukey"` (default) or `"chi-squared"`.
#' @param n.samples Number of posterior samples (default 500).
#' @return A list with `fit.y`, `fit.y.rep`, and `bayesian.p`.
#' @export
ppcOccu <- function(object, fit.stat = c("freeman-tukey", "chi-squared"),
                    n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  model <- object$model
  if (model$model_type != "single") stop("ppcOccu supports single-season only")

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

  fit_y <- fit_y_rep <- numeric(n.samples)
  for (s in seq_len(n.samples)) {
    idx <- draw_idx[s]
    psi <- plogis(as.vector(X_occ %*% draws[idx, seq_len(p_occ)]))
    p <- plogis(as.vector(X_det %*% draws[idx, p_occ + seq_len(p_det)]))
    z <- rbinom(n_sites, 1, psi)
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
pitResiduals <- function(object, n.samples = 250) {
  model <- object$model
  if (model$model_type != "single") stop("pitResiduals supports single-season only")
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
testUniformity <- function(pit) {
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
testDispersion <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  obs_var <- var(rowSums(y_obs * (y_obs >= 0), na.rm = TRUE))
  sim_vars <- vapply(sims, function(ys) var(rowSums(ys * (ys >= 0), na.rm = TRUE)), double(1))
  list(observed = obs_var, expected = mean(sim_vars),
       ratio = obs_var / mean(sim_vars), p.value = mean(sim_vars >= obs_var))
}

#' @export
testZeroInflation <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  count_zeros <- function(y) sum(apply(y, 1, function(r) { v <- r >= 0; all(r[v] == 0) }))
  obs <- count_zeros(y_obs); sim <- vapply(sims, count_zeros, integer(1))
  list(observed = obs, expected = mean(sim), ratio = obs / max(mean(sim), 1),
       p.value = mean(sim >= obs))
}

#' @export
testOutliers <- function(object, n.samples = 250) {
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
checkModel <- function(object, coords = NULL, n.samples = 250) {
  cat("=== tobs Model Diagnostics ===\n\n")
  cat(sprintf("Sampler: %d samples, %d divergent, mean accept = %.3f\n",
              object$n_samples, sum(object$divergent), mean(object$accept_prob)))

  w <- tryCatch(waicOccu(object), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("\nWAIC: %.1f (p_waic = %.1f)\n", w$waic, w$p_waic))

  ppc <- tryCatch(ppcOccu(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(ppc)) {
    cat(sprintf("\nPPC: Bayesian p = %.3f\n", ppc$bayesian.p))
    if (ppc$bayesian.p < 0.05 || ppc$bayesian.p > 0.95) cat("  WARNING: poor fit\n")
  }

  zi <- tryCatch(testZeroInflation(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(zi)) cat(sprintf("\nZero-inflation: obs=%d, exp=%.1f, p=%.3f\n", zi$observed, zi$expected, zi$p.value))

  if (!is.null(coords)) {
    mi <- tryCatch(tulpa::moranI(residuals(object)$occ, coords), error = function(e) NULL)
    if (!is.null(mi)) {
      cat(sprintf("\nMoran's I: %.3f (p = %.3f)\n", mi$I, mi$p.value))
      if (mi$p.value < 0.05) cat("  WARNING: spatial autocorrelation\n")
    }
  }
  cat("\n")
  invisible(list(waic = w, ppc = ppc, zero_inflation = zi))
}
