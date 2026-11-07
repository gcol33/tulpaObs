# removal_nuts.R - NUTS target density for the removal-sampling family (removal()).
#
# Same flat coefficient block as the N-mixture NUTS path (abun_nuts.R) -- there
# is no latent field or community covariance, only
# theta = (beta_lambda, beta_p[, log_r]) plus the optional single-intercept RE
# tail -- so the layout (.tobs_abun_nuts_layout) and the shared sampler front
# door (.tobs_count_nuts_front_door, R/nuts_chains.R) are reused. The only
# difference is the per-site marginal: the removal depleting-binomial product
# (cpp_removal_total_log_lik) instead of the N-mixture product. The C++ FullGradFn
# (src/removal_nuts.cpp) mirrors this R target and is cross-checked against it.

# Per-site removal marginal closure (the NUTS oracle's data + eval_beta), mirroring
# nmix_site_marginal() but over the depleting-binomial removal likelihood.
.tobs_removal_nuts_marginal <- function(model, mixture = "P", K_max = NULL) {
  y        <- as.integer(model$y_long)
  site_idx <- as.integer(model$site_idx)
  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  n_sites  <- model$n_sites
  if (is.null(K_max)) {
    site_tot <- tapply(y, factor(site_idx, levels = seq_len(n_sites)), sum)
    site_tot[is.na(site_tot)] <- 0L
    K_max <- as.integer(max(as.integer(site_tot)) + 100L)
  }
  K_max <- as.integer(K_max)
  resolve_r <- function(r) {
    if (identical(mixture, "P")) return(Inf)
    if (is.null(r) || !is.finite(r) || r <= 0) {
      stop("NB removal marginal requires a finite positive `r`.", call. = FALSE)
    }
    as.numeric(r)
  }
  eval_beta <- function(beta_lambda, beta_p, r = Inf) {
    eta_lambda <- as.numeric(X_lambda %*% beta_lambda)
    eta_p      <- as.numeric(X_p %*% beta_p)
    out <- cpp_removal_total_log_lik(y, site_idx, eta_p, eta_lambda, K_max,
                                     r = resolve_r(r))
    out
  }
  list(X_lambda = X_lambda, X_p = X_p, K_max = K_max, mixture = mixture,
       eval_beta = eval_beta)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the removal family
# ---------------------------------------------------------------------------

# Sample the exact coefficient posterior of a non-spatial removal model via
# tulpa's NUTS engine and the in-tree C++ FullGradFn (cpp_removal_nuts),
# warm-started at the Laplace mode with a diagonal Laplace metric, then package
# the draws into the build_nmix_fit shape. Operates on the (autoscaled) `model`;
# the caller .tobs_fit_model() unscales the means / draws / vcov back to natural.
.tobs_fit_removal_nuts <- function(model, mixture = "poisson", K_max = NULL,
                                   sigma.beta = NULL, sigma.logr = NULL, re = NULL,
                                   n.iter = NULL, n.warmup = NULL, n.chains = NULL, n.thin = NULL,
                                   n.threads = NULL,
                                   max.treedepth = NULL, adapt.delta = NULL,
                                   seed = NULL, verbose = FALSE) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  is_nb    <- identical(mixture, "negbin")
  mix_code <- if (is_nb) "NB" else "P"
  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  y_long   <- model$y_long
  site_idx <- model$site_idx
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  n_sites <- model$n_sites
  if (is.null(K_max)) {
    site_tot <- tapply(y_long, factor(site_idx, levels = seq_len(n_sites)), sum)
    site_tot[is.na(site_tot)] <- 0L
    K_max <- max(as.integer(site_tot)) + 100L
  }
  K_max <- as.integer(K_max)

  # Single intercept RE on one arm, via the shared count-NUTS RE helpers (same
  # machinery as abun()).
  re_info <- .tobs_count_nuts_re_info(re, model)
  n_re_groups <- if (!is.null(re_info)) re_info$n_groups else 0L
  lay <- .tobs_abun_nuts_layout(p_lam, p_p, is_nb, re_groups = n_re_groups)

  warm <- removal_laplace(y = y_long, site_idx = site_idx, X_lambda = X_lambda,
                          X_p = X_p, mixture = mix_code, K_max = K_max,
                          max_iter = 100L, verbose = FALSE)
  run <- .tobs_count_nuts_front_door(
    cpp_removal_nuts,
    spec = list(y = as.integer(y_long), site_idx = as.integer(site_idx),
                X_lambda = X_lambda, X_p = X_p,
                n_sites = n_sites, K_max = K_max, is_nb = is_nb),
    priors = list(sigma_beta = sigma.beta, sigma_logr = sigma.logr),
    theta0 = c(warm$beta_lambda, warm$beta_p, if (is_nb) warm$log_r),
    V = warm$vcov,
    nms = c(paste0("lambda_", model$process_info[[1]]$coef_names),
            paste0("p_",      model$process_info[[2]]$coef_names),
            if (is_nb) "log_r"),
    re_info = re_info, sigma.logr = sigma.logr,
    sampler = .tobs_sampler_control_snapshot(environment()), verbose = verbose)
  par <- run$par; cov <- run$cov

  marg <- .tobs_removal_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$p],
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik

  raw <- list(
    mixture = mix_code,
    beta_lambda = unname(par[lay$lambda]),
    beta_p      = unname(par[lay$p]),
    log_r       = if (is_nb) unname(par[lay$log_r]) else NA_real_,
    r           = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
    vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max,
    mean_N = warm$mean_N, var_N = warm$var_N,
    boundary_weight = warm$boundary_weight)
  fit <- build_nmix_fit(raw, model, spatial = NULL)

  .tobs_count_nuts_attach(
    fit, run, ll_mean,
    extra = list(is_nb = is_nb, K_max = K_max,
                 re_arm = if (!is.null(re_info)) re_info$arm else -1L,
                 sigma_beta = sigma.beta, sigma_logr = sigma.logr))
}
