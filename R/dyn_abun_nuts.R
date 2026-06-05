# dyn_abun_nuts.R - NUTS target for the Dail-Madsen open N-mixture family.
#
# Flat coefficient vector theta = (beta_lambda, beta_p, beta_omega, beta_gamma);
# the joint log-posterior is the forward marginal (cpp_dyn_abun_total_log_lik)
# plus weak Gaussian priors. The C++ FullGradFn (src/dyn_abun_nuts.cpp) mirrors
# this R target and is cross-checked against it.

.tobs_dyn_abun_nuts_layout <- function(p_lam, p_p, p_om, p_gm) {
  off <- cumsum(c(0L, p_lam, p_p, p_om))
  list(p_lam = p_lam, p_p = p_p, p_om = p_om, p_gm = p_gm,
       lambda = off[1] + seq_len(p_lam), p = off[2] + seq_len(p_p),
       omega = off[3] + seq_len(p_om), gamma = off[4] + seq_len(p_gm),
       total = p_lam + p_p + p_om + p_gm)
}

.tobs_dyn_abun_nuts_marginal <- function(model) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_omega  <- model$X_processes[[3]]; X_gamma <- model$X_processes[[4]]
  y_flat <- as.integer(model$y_flat)
  n_sites <- model$n_sites; T <- model$n_seasons; J <- model$max_visits
  K <- model$K_max
  eval_beta <- function(beta_lambda, beta_p, beta_omega, beta_gamma) {
    cpp_dyn_abun_total_log_lik(
      y_flat, n_sites, T, J, K,
      as.numeric(X_lambda %*% beta_lambda), as.numeric(X_p %*% beta_p),
      as.numeric(X_omega %*% beta_omega), as.numeric(X_gamma %*% beta_gamma))
  }
  list(X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
       eval_beta = eval_beta)
}

.tobs_dyn_abun_nuts_logpost <- function(theta, marg, lay, sigma.beta = 10) {
  ev <- marg$eval_beta(theta[lay$lambda], theta[lay$p], theta[lay$omega], theta[lay$gamma])
  lp <- ev$log_lik
  grad <- numeric(lay$total)
  grad[lay$lambda] <- as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda))
  grad[lay$p]      <- as.numeric(crossprod(marg$X_p,      ev$grad_eta_p))
  grad[lay$omega]  <- as.numeric(crossprod(marg$X_omega,  ev$grad_eta_omega))
  grad[lay$gamma]  <- as.numeric(crossprod(marg$X_gamma,  ev$grad_eta_gamma))
  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * sum(theta^2)
  grad <- grad - ib2 * theta
  list(lp = lp, grad = grad)
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
  lay <- .tobs_dyn_abun_nuts_layout(ncol(X_lambda), ncol(X_p), ncol(X_omega), ncol(X_gamma))

  warm <- dyn_abun_laplace(
    y_flat = model$y_flat, n_sites = model$n_sites, T = model$n_seasons,
    J = model$max_visits, K_max = model$K_max,
    X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
    verbose = FALSE)
  theta0 <- warm$means
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == lay$total) && all(is.finite(diag(V))))
                  pmax(diag(V), 1e-6) else rep(1, lay$total)

  spec <- list(y = as.integer(model$y_flat), n_sites = model$n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma)

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
  colnames(draws) <- nms
  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- chains[[1L]]$epsilon

  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws)
  marg <- .tobs_dyn_abun_nuts_marginal(model)
  ev_mean <- marg$eval_beta(par[lay$lambda], par[lay$p], par[lay$omega], par[lay$gamma])

  raw <- list(means = unname(par), vcov = cov, coef_names = nms,
              log_lik = ev_mean$log_lik, log_lik_site = ev_mean$log_lik_site,
              mean_N1 = ev_mean$mean_N1, K_max = model$K_max, converged = TRUE)
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
