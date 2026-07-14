# =============================================================================
# count_spatial.R - areal-spatial count / relative-abundance GLMM
# (count() + a plain areal field; gcol33/tulpaObs#117). The spAbundance spAbund
# analogue.
#
# The count response is observed directly (no detection, no latent state), so a
# count areal fit is one GLMM block plus a shared ICAR / proper-CAR field:
#
#   log mu_i = X_i . beta + f_{u(i)},   y_i ~ Poisson(mu_i)
#
# There is nothing to marginalise out, so this is a single
# tulpa::tulpa_nested_laplace() call over the count block with the areal GMRF
# as its latent prior -- NOT the occupancy EM (whose E-step integrates the
# binary state). The field hyperparameters are integrated on the nested outer
# grid; the fixed effects and their covariance are grid-integrated here (law of
# total covariance over the grid, using the per-cell FE Hessians the engine
# returns under keep_grid_hessians). One field node per site (identity map).
#
# Poisson only: with one field node per site a negbin size / gaussian residual
# variance is not jointly identified with the field under the fixed-dispersion
# nested loop (the areal-count gate in .dispatch_count enforces this).
# =============================================================================


# Grid-integrated fixed-effect mean + covariance from a nested-Laplace fit run
# with keep_grid_hessians = TRUE. `res$grid_modes[[k]]` is the per-cell FE mode
# and `res$grid_hessians[[k]]` the per-cell FE precision (the field Schur-folded
# out); the marginal covariance is the law of total covariance over the outer
# grid: V = sum_k w_k [C_k + (m_k - beta)(m_k - beta)'], C_k = solve(H_k).
.count_spatial_fe_moments <- function(res, p) {
  w  <- tulpa:::.nl_normalise_weights_safe(res$log_marginal)
  ok <- is.finite(w) & w > 0
  if (!any(ok)) {
    stop("Areal count: every grid point produced a non-finite log-marginal. ",
         "Check the adjacency graph / the field hyperparameter grid.",
         call. = FALSE)
  }
  w[!ok] <- 0; w <- w / sum(w)

  modes <- res$grid_modes
  hess  <- res$grid_hessians
  if (is.null(modes) || is.null(hess)) {
    stop("Areal count: the nested fit did not return grid Hessians. This is a ",
         "tulpaObs bug (keep_grid_hessians was not honoured).", call. = FALSE)
  }

  beta <- numeric(p)
  for (k in which(ok)) beta <- beta + w[k] * as.numeric(modes[[k]])

  V <- matrix(0, p, p)
  for (k in which(ok)) {
    Hk <- as.matrix(hess[[k]])
    Ck <- tryCatch(solve(Hk), error = function(e) NULL)
    if (is.null(Ck) || anyNA(Ck)) next
    dk <- as.numeric(modes[[k]]) - beta
    V <- V + w[k] * (Ck + tcrossprod(dk))
  }
  V <- (V + t(V)) / 2
  list(beta = beta, vcov = V, weights = w)
}


# Fit an areal-spatial count model. `model` is the (autoscaled) count tobs_model;
# `spatial` is the resolved tobs_spatial areal term on the abundance formula.
# Returns a `tobs_fit`; the caller (`.tobs_fit_model`) transforms the per-process
# betas / SEs / draws back to natural scale via `.unscale_fit_per_process`.
.tobs_fit_count_spatial <- function(model, spatial, max_iter = 50L, tol = 1e-6,
                                    verbose = FALSE, ...) {
  if (!identical(model$model_type, "count")) {
    stop("`.tobs_fit_count_spatial` expects a count model.", call. = FALSE)
  }
  prior <- .tobs_to_multi_block_prior(spatial = spatial, model = model)
  if (is.null(prior)) {
    stop("Areal count needs an areal field block (icar / car_proper); none was ",
         "resolved from the spatial term.", call. = FALSE)
  }

  X <- model$X_processes[[1]]
  p <- ncol(X)
  y <- as.numeric(model$y_count)
  N <- length(y)
  fam <- switch(model$response %||% "poisson",
    poisson  = "poisson",
    negbin   = "neg_binomial_2",
    gaussian = "gaussian",
    stop(sprintf("Areal count: unsupported response '%s'.", model$response),
         call. = FALSE))

  res <- tulpa::tulpa_nested_laplace(
    y = y, n_trials = rep(1L, N), X = X, prior = prior,
    family = fam, phi = as.numeric(model$count_phi %||% 1.0),
    control = list(max_iter = as.integer(max_iter), tol = as.numeric(tol),
                   keep_grid_hessians = TRUE, verbose = isTRUE(verbose))
  )

  fe <- .count_spatial_fe_moments(res, p)
  pi_list <- model$process_info
  nms <- paste0(pi_list[[1L]]$name, "_", pi_list[[1L]]$coef_names)

  means <- fe$beta; names(means) <- nms
  V <- fe$vcov; dimnames(V) <- list(nms, nms)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nms

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- nms

  fit <- structure(c(list(
    draws = draws, means = means, sds = sds,
    skew = NULL, sla_status = "off",
    n_samples = n_draws, n_params = length(means),
    log_prob = rep(NA_real_, n_draws)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names = nms, param_names = nms,
    intercepts = compute_intercepts(model, means),
    model = model,
    spatial = spatial,
    process_info = pi_list,
    method = "nested_laplace",
    nested_laplace = list(multi_prior = prior, occ_fit = res),
    convergence = list(converged = TRUE, n_iter = as.integer(res$n_iter %||% 1L)),
    correction = "none"
  )), class = c("tobs_fit", "tulpa_fit"))

  # Shared field summary: the areal block's SD (sigma = 1/sqrt(tau) marginalized
  # over the outer grid) and the demeaned per-cell field, read off `res`. Reuses
  # the occu nested field summary so `fit$spatial_field` / `fit$field_table` /
  # `sigma` match what occu() / occu_cover() expose.
  .tobs_nested_attach_field_summary(fit, model, res, prior)
}
