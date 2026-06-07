# dyn_abun_nuts.R - NUTS target for the Dail-Madsen open N-mixture family.
#
# Flat coefficient vector theta = (beta_lambda, beta_p, beta_omega, beta_gamma);
# the joint log-posterior is the forward marginal (cpp_dyn_abun_total_log_lik)
# plus weak Gaussian priors. The C++ FullGradFn (src/dyn_abun_nuts.cpp) mirrors
# this R target and is cross-checked against it.

# `use_nb` appends a single trailing log r coordinate (NB initial abundance).
.tobs_dyn_abun_nuts_layout <- function(p_lam, p_p, p_om, p_gm, use_nb = FALSE) {
  idx <- .tobs_nuts_arm_idx(c("lambda", "p", "omega", "gamma"),
                            c(p_lam, p_p, p_om, p_gm))
  base <- p_lam + p_p + p_om + p_gm
  out <- c(list(p_lam = p_lam, p_p = p_p, p_om = p_om, p_gm = p_gm),
           idx, list(use_nb = use_nb))
  if (use_nb) out <- c(out, list(logr = base + 1L, total = base + 1L))
  else        out <- c(out, list(logr = NA_integer_, total = base))
  out
}

.tobs_dyn_abun_nuts_marginal <- function(model) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_omega  <- model$X_processes[[3]]; X_gamma <- model$X_processes[[4]]
  y_flat <- as.integer(model$y_flat)
  n_sites <- model$n_sites; T <- model$n_seasons; J <- model$max_visits
  K <- model$K_max
  use_nb <- identical(model$mixture %||% "poisson", "negbin")
  # eta_logr defaults to 0; the caller passes the current log r under NB.
  eval_beta <- function(beta_lambda, beta_p, beta_omega, beta_gamma, eta_logr = 0) {
    cpp_dyn_abun_total_log_lik(
      y_flat, n_sites, T, J, K,
      as.numeric(X_lambda %*% beta_lambda), as.numeric(X_p %*% beta_p),
      as.numeric(X_omega %*% beta_omega), as.numeric(X_gamma %*% beta_gamma),
      use_nb = use_nb, eta_logr = as.numeric(eta_logr))
  }
  list(X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
       use_nb = use_nb, eval_beta = eval_beta)
}

.tobs_dyn_abun_nuts_logpost <- function(theta, marg, lay, sigma.beta = 10) {
  use_nb <- isTRUE(lay$use_nb)
  eta_logr <- if (use_nb) theta[lay$logr] else 0
  ev <- marg$eval_beta(theta[lay$lambda], theta[lay$p], theta[lay$omega],
                       theta[lay$gamma], eta_logr)
  arms <- list(
    list(idx = lay$lambda, X = marg$X_lambda, grad = "grad_eta_lambda"),
    list(idx = lay$p,      X = marg$X_p,      grad = "grad_eta_p"),
    list(idx = lay$omega,  X = marg$X_omega,  grad = "grad_eta_omega"),
    list(idx = lay$gamma,  X = marg$X_gamma,  grad = "grad_eta_gamma"))
  # The dispersion log r has no design: its gradient is the scalar
  # grad_eta_logr (already summed over sites), folded in as a 1x1 arm.
  if (use_nb) {
    ev$grad_eta_logr <- as.numeric(ev$grad_eta_logr)
    arms <- c(arms, list(list(idx = lay$logr, X = matrix(1, 1, 1),
                              grad = "grad_eta_logr")))
  }
  .tobs_nuts_logpost_k(theta, ev, arms, lay$total, sigma.beta)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the open N-mixture family
# ---------------------------------------------------------------------------

.tobs_fit_dyn_abun_nuts <- function(model, sigma.beta = 10,
                                    n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
                                    max.treedepth = 10L, adapt.delta = 0.9,
                                    seed = 1L, verbose = FALSE) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_omega  <- model$X_processes[[3]]; X_gamma <- model$X_processes[[4]]
  use_nb <- identical(model$mixture %||% "poisson", "negbin")
  lay <- .tobs_dyn_abun_nuts_layout(ncol(X_lambda), ncol(X_p), ncol(X_omega),
                                    ncol(X_gamma), use_nb = use_nb)

  warm <- dyn_abun_laplace(
    y_flat = model$y_flat, n_sites = model$n_sites, T = model$n_seasons,
    J = model$max_visits, K_max = model$K_max,
    X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
    mixture = model$mixture %||% "poisson", verbose = FALSE)
  theta0 <- warm$means
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == lay$total) && all(is.finite(diag(V))))
                  pmax(diag(V), 1e-6) else rep(1, lay$total)

  spec <- list(y = as.integer(model$y_flat), n_sites = model$n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
               use_nb = use_nb)

  run_chain <- function(ch) {
    cpp_dyn_abun_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      inv_metric = inv_metric,
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
           paste0("omega_",  model$process_info[[3]]$coef_names),
           paste0("gamma_",  model$process_info[[4]]$coef_names))
  if (use_nb) nms <- c(nms, "log_r")
  colnames(draws) <- nms
  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- chains[[1L]]$epsilon

  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws)
  marg <- .tobs_dyn_abun_nuts_marginal(model)
  log_r <- if (use_nb) as.numeric(par[lay$logr]) else NA_real_
  ev_mean <- marg$eval_beta(par[lay$lambda], par[lay$p], par[lay$omega],
                            par[lay$gamma], if (use_nb) log_r else 0)

  raw <- list(means = unname(par), vcov = cov, coef_names = nms,
              log_lik = ev_mean$log_lik, log_lik_site = ev_mean$log_lik_site,
              mean_N1 = ev_mean$mean_N1, K_max = model$K_max, converged = TRUE,
              mixture = model$mixture %||% "poisson",
              log_r = log_r, r = if (use_nb) exp(log_r) else NA_real_)
  fit <- build_dyn_abun_fit(raw, model)

  n_draws <- nrow(draws)
  fit$draws       <- draws
  fit$n_samples   <- n_draws
  fit$log_prob    <- rep(ev_mean$log_lik, n_draws)
  fit$accept_prob <- accept
  fit$divergent   <- divergent
  fit$treedepth   <- treedepth
  fit$epsilon     <- epsilon
  fit$method      <- "nuts"
  fit$nuts <- list(accept_prob = accept, divergent = divergent,
                   treedepth = treedepth, epsilon = epsilon,
                   n_chains = as.integer(n.chains),
                   divergent_total = sum(divergent), sigma_beta = sigma.beta,
                   K_max = model$K_max)
  fit
}
