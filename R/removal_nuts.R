# removal_nuts.R - NUTS target density for the removal-sampling family (removal()).
#
# Same flat coefficient block as the N-mixture NUTS path (abun_nuts.R) -- there
# is no latent field, random effect, or community covariance, only
# theta = (beta_lambda, beta_p[, log_r]) -- so the layout / warm-start helpers
# (.tobs_abun_nuts_layout / .tobs_abun_nuts_pack_init) are reused. The only
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

# Joint log-posterior + gradient of the removal coefficient vector. `marg` is a
# `.tobs_removal_nuts_marginal()` object; `lay` the shared abun layout. Weak
# Gaussian priors keep the ridge off the boundary. Returns list(lp, grad). This
# R version is the oracle the C++ FullGradFn (cpp_removal_nuts_joint_logpost) is
# cross-checked against.
.tobs_removal_nuts_logpost <- function(theta, marg, lay,
                                       sigma.beta = 10, sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_p      <- theta[lay$p]
  r <- if (lay$is_nb) exp(theta[lay$log_r]) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_p, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_p,      ev$grad_eta_p)))
  if (lay$is_nb) grad <- c(grad, sum(ev$grad_theta))

  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * (sum(beta_lambda^2) + sum(beta_p^2))
  grad[lay$lambda] <- grad[lay$lambda] - ib2 * beta_lambda
  grad[lay$p]      <- grad[lay$p]      - ib2 * beta_p
  if (lay$is_nb) {
    lr  <- theta[lay$log_r]; ilr2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * ilr2 * lr^2
    grad[lay$log_r] <- grad[lay$log_r] - ilr2 * lr
  }
  list(lp = lp, grad = grad)
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
                                   sigma.beta = 10, sigma.logr = 1.5,
                                   n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
                                   max.treedepth = 10L, adapt.delta = 0.9,
                                   seed = 1L, verbose = FALSE) {
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
  lay <- .tobs_abun_nuts_layout(p_lam, p_p, is_nb)

  warm <- removal_laplace(y = y_long, site_idx = site_idx, X_lambda = X_lambda,
                          X_p = X_p, mixture = mix_code, K_max = K_max,
                          max_iter = 100L, verbose = FALSE)
  init <- .tobs_abun_nuts_pack_init(warm, lay)

  spec <- list(y = as.integer(y_long), site_idx = as.integer(site_idx),
               X_lambda = X_lambda, X_p = X_p,
               n_sites = n_sites, K_max = K_max, is_nb = is_nb)

  run_chain <- function(ch) {
    cpp_removal_nuts(spec, theta0 = init$theta0,
                     sigma_beta = sigma.beta, sigma_logr = sigma.logr,
                     inv_metric = init$inv_metric,
                     n_iter = as.integer(n.iter + n.warmup),
                     n_warmup = as.integer(n.warmup),
                     max_treedepth = as.integer(max.treedepth),
                     adapt_delta = adapt.delta,
                     seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           if (is_nb) "log_r")
  colnames(draws) <- nms
  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- chains[[1L]]$epsilon

  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws)
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

  n_draws <- nrow(draws)
  fit$draws       <- draws
  fit$n_samples   <- n_draws
  fit$log_prob    <- rep(ll_mean, n_draws)
  fit$accept_prob <- accept
  fit$divergent   <- divergent
  fit$treedepth   <- treedepth
  fit$epsilon     <- epsilon
  fit$method      <- "nuts"
  fit$nuts <- list(accept_prob = accept, divergent = divergent,
                   treedepth = treedepth, epsilon = epsilon,
                   n_chains = as.integer(n.chains),
                   divergent_total = sum(divergent),
                   is_nb = is_nb, K_max = K_max,
                   sigma_beta = sigma.beta, sigma_logr = sigma.logr)
  fit
}
