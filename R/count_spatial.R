# =============================================================================
# count_spatial.R
# - areal-spatial count / relative-abundance GLMM (count() + a plain areal
# field). The spAbundance spAbund analogue.
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

  mode_mat <- do.call(rbind, lapply(modes, as.numeric))
  beta     <- as.numeric(crossprod(w, mode_mat))

  # The kernel returns per-cell FE precisions; the law of total covariance takes
  # covariances, so invert each cell first (a singular cell drops out).
  cov_k <- lapply(hess, function(H) tryCatch(solve(as.matrix(H)),
                                             error = function(e) NULL))
  V <- .tobs_grid_vcov(mode_mat, w, cov_k, center = beta, symmetrize = TRUE)
  list(beta = beta, vcov = V, weights = w)
}


# Fit an areal-spatial count model. `model` is the (autoscaled) count tobs_model;
# `spatial` is the resolved tobs_spatial areal term on the abundance formula.
# Returns a `tobs_fit`; the caller (`.tobs_fit_model`) transforms the per-process
# betas / SEs / draws back to natural scale via `.unscale_fit_per_process`.
.tobs_fit_count_spatial <- function(model, spatial, max_iter = 50L, tol = 1e-6,
                                    sigma.grid = NULL, rho.grid = NULL,
                                    tau.grid = NULL, range.grid = NULL,
                                    verbose = FALSE, ...) {
  if (!identical(model$model_type, "count")) {
    stop("`.tobs_fit_count_spatial` expects a count model.", call. = FALSE)
  }
  # A continuous NNGP Gaussian-process field on the abundance arm (gcol33/
  # follow-up). tulpa's nested-Laplace hosts a single-block `nngp` kernel (its
  # own cpp_fn, integrated over the GP marginal variance and range), so this
  # routes around the multi-block areal builder to tulpa directly.
  if (identical(spatial$type, "gp"))
    return(.tobs_fit_count_gp(model, spatial, max_iter = max_iter, tol = tol,
                              verbose = verbose))
  if (identical(spatial$type, "multiscale_gp")) {
    stop("Areal count: the two-scale multiscale_gp() field is not hosted by the ",
         "nested-Laplace engine (a single-scale continuous field). Use gp() for ",
         "a one-scale NNGP field, or spde() for a mesh-based continuous Matern ",
         "field with a reconstructed per-cell map.", call. = FALSE)
  }
  prior <- .tobs_to_multi_block_prior(
    spatial = spatial, model = model,
    grids = .tobs_outer_grids(sigma.grid, rho.grid, tau.grid, range.grid))
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
    binomial = "binomial",
    stop(sprintf("Areal count: unsupported response '%s'.", model$response),
         call. = FALSE))
  # Binomial carries a per-site trial count; every other areal response is one
  # trial per site (the identity used by the Poisson / Gaussian field fit).
  n_trials <- if (identical(fam, "binomial"))
                as.integer(model$n_trials %||% rep(1L, N)) else rep(1L, N)

  res <- tulpa::tulpa_nested_laplace(
    y = y, n_trials = n_trials, X = X, prior = prior,
    family = fam, phi = as.numeric(model$count_phi %||% 1.0),
    # `verbose` is not a tulpa_nested_laplace() knob (it validates its control
    # names); the driver's own reporting is the `progress` option tobs() scopes.
    control = list(max_iter = as.integer(max_iter), tol = as.numeric(tol),
                   keep_grid_hessians = TRUE)
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

  # Shared field summary: each areal block's SD (sigma = 1/sqrt(tau) marginalized
  # over the outer grid) and its demeaned per-cell field, read off `res`. Reuses
  # the occu nested field summary so `fit$spatial_field` / `fit$trend_fields` /
  # `fit$field_table` / `sigma` match what occu() / occu_cover() expose. A
  # varying-coefficient bar contributes an intercept field plus one field per
  # covariate (svcAbund), all summarized by the same loop.
  # The per-site field offset (sum_k W[i,k] f_k[i]) is attached by the caller
  # (.tobs_fit_model), which swaps in the unscaled model afterwards -- attaching
  # it here would be dropped by that swap.
  .tobs_nested_attach_field_summary(fit, model, res, prior)
}


# ---------------------------------------------------------------------------
# Continuous NNGP Gaussian-process field on the count arm (gp())
# ---------------------------------------------------------------------------

# Build the tulpa `nngp` nested-Laplace prior from a resolved gp() spatial term.
# The tobs gp() term caches the coordinate matrix (row-major, `spatial$coords`)
# and the covariance settings; tulpa's `spatial_gp()` + `prior_from_spec()` are
# the single source of truth for the NNGP neighbour structure and the nested
# prior block, so the coordinates are handed back to them rather than re-deriving
# the block here.
.count_gp_prior <- function(spatial, n_sites) {
  co <- matrix(as.numeric(spatial$coords), ncol = 2L, byrow = TRUE)
  if (nrow(co) != n_sites) {
    stop(sprintf(paste0("Areal count gp(): the term carries %d coordinates but ",
         "the model has %d sites."), nrow(co), n_sites), call. = FALSE)
  }
  tmp <- data.frame(.gp_lon = co[, 1L], .gp_lat = co[, 2L])
  gspec <- tulpa::spatial_gp(
    coords = c(".gp_lon", ".gp_lat"),
    cov    = spatial$cov_type %||% "exponential",
    nu     = spatial$nu %||% 1.5,
    nn     = as.integer(spatial$nn %||% 15L))
  tulpa::prior_from_spec(gspec, tmp)
}

# Fit a count / relative-abundance GLMM with a continuous NNGP GP field on the
# abundance arm. tulpa's nested-Laplace `nngp` kernel integrates the GP marginal
# variance and range on its own outer grid and Schur-folds the field out, so this
# returns grid-integrated fixed effects (via the shared `.count_spatial_fe_moments`)
# plus the GP hyperparameter posterior; the per-cell field itself is integrated out
# (use spde() for a reconstructed continuous field map). Poisson / binomial only,
# as for the areal path (a negbin size / gaussian residual variance is not jointly
# identified with a per-node field under the fixed-dispersion nested loop).
.tobs_fit_count_gp <- function(model, spatial, max_iter = 50L, tol = 1e-6,
                               verbose = FALSE, ...) {
  X <- model$X_processes[[1]]
  p <- ncol(X)
  y <- as.numeric(model$y_count)
  N <- length(y)
  resp <- model$response %||% "poisson"
  if (!resp %in% c("poisson", "binomial")) {
    stop(sprintf(paste0("Areal count gp(): response '%s' is not identified ",
         "against a continuous per-node GP field (the size / residual variance ",
         "and the field both absorb overdispersion). Use 'poisson' / 'binomial', ",
         "or spde() / an areal field."), resp), call. = FALSE)
  }
  fam <- if (identical(resp, "binomial")) "binomial" else "poisson"
  n_trials <- if (identical(fam, "binomial"))
                as.integer(model$n_trials %||% rep(1L, N)) else rep(1L, N)

  prior <- .count_gp_prior(spatial, N)
  res <- tulpa::tulpa_nested_laplace(
    y = y, n_trials = n_trials, X = X, prior = prior,
    family = fam, phi = as.numeric(model$count_phi %||% 1.0),
    control = list(max_iter = as.integer(max_iter), tol = as.numeric(tol),
                   keep_grid_hessians = TRUE))

  fe <- .count_spatial_fe_moments(res, p)
  pi_list <- model$process_info
  nms <- paste0(pi_list[[1L]]$name, "_", pi_list[[1L]]$coef_names)

  means <- fe$beta; names(means) <- nms
  V <- fe$vcov; dimnames(V) <- list(nms, nms)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nms

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- nms

  # GP hyperparameter posterior (marginal SD sqrt(sigma2) and the range phi_gp),
  # grid-integrated by tulpa; surfaced on fit$spatial for the user.
  gp_hyper <- list(theta_names = res$theta_names,
                   mean = res$theta_mean, sd = res$theta_sd,
                   median = res$theta_median,
                   ci_lo = res$theta_ci_lo, ci_hi = res$theta_ci_hi)

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
    spatial_field = NULL,          # GP field integrated out (no per-cell map)
    gp_hyper = gp_hyper,
    process_info = pi_list,
    method = "nested_laplace",
    nested_laplace = list(prior = prior, occ_fit = res),
    convergence = list(converged = TRUE, n_iter = as.integer(res$n_iter %||% 1L)),
    correction = "none"
  )), class = c("tobs_fit", "tulpa_fit"))
  fit
}
