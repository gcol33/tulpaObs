# fp_occu_nuts.R - NUTS target for the multistate false-positive occupancy family.
#
# Flat coefficient vector theta = (beta_psi, beta_p11, beta_p10, beta_b); the
# joint log-posterior is the false-positive occupancy marginal
# (cpp_fp_occu_total_log_lik) plus weak Gaussian priors. The C++ FullGradFn
# (src/fp_occu_nuts.cpp) mirrors this R target and is cross-checked against it.

.tobs_fp_occu_nuts_layout <- function(p_psi, p_p11, p_p10, p_b) {
  idx <- .tobs_nuts_arm_idx(c("psi", "p11", "p10", "b"),
                            c(p_psi, p_p11, p_p10, p_b))
  c(list(p_psi = p_psi, p_p11 = p_p11, p_p10 = p_p10, p_b = p_b),
    idx, list(total = p_psi + p_p11 + p_p10 + p_b))
}

.tobs_fp_occu_nuts_marginal <- function(model) {
  y        <- as.integer(model$y_long)
  site_idx <- as.integer(model$site_idx)
  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]
  eval_beta <- function(beta_psi, beta_p11, beta_p10, beta_b) {
    cpp_fp_occu_total_log_lik(
      y, site_idx,
      as.numeric(X_psi %*% beta_psi), as.numeric(X_p11 %*% beta_p11),
      as.numeric(X_p10 %*% beta_p10), as.numeric(X_b %*% beta_b))
  }
  list(X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
       eval_beta = eval_beta)
}

.tobs_fp_occu_nuts_logpost <- function(theta, marg, lay, sigma.beta = 10) {
  ev <- marg$eval_beta(theta[lay$psi], theta[lay$p11], theta[lay$p10], theta[lay$b])
  arms <- list(
    list(idx = lay$psi, X = marg$X_psi, grad = "grad_eta_psi"),
    list(idx = lay$p11, X = marg$X_p11, grad = "grad_eta_p11"),
    list(idx = lay$p10, X = marg$X_p10, grad = "grad_eta_p10"),
    list(idx = lay$b,   X = marg$X_b,   grad = "grad_eta_b"))
  .tobs_nuts_logpost_k(theta, ev, arms, lay$total, sigma.beta)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the false-positive occupancy family
# ---------------------------------------------------------------------------

.tobs_fit_fp_occu_nuts <- function(model, sigma.beta = 10,
                                   n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
                                   max.treedepth = 10L, adapt.delta = 0.9,
                                   seed = 1L, verbose = FALSE) {
  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]
  lay <- .tobs_fp_occu_nuts_layout(ncol(X_psi), ncol(X_p11), ncol(X_p10), ncol(X_b))

  warm <- fp_occu_laplace(y = model$y_long, site_idx = model$site_idx,
                          X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
                          verbose = FALSE)
  theta0 <- warm$means
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == lay$total) && all(is.finite(diag(V))))
                  pmax(diag(V), 1e-6) else rep(1, lay$total)

  spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
               X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
               n_sites = model$n_sites)

  run_chain <- function(ch) {
    cpp_fp_occu_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                     inv_metric = inv_metric,
                     n_iter = as.integer(n.iter + n.warmup),
                     n_warmup = as.integer(n.warmup),
                     max_treedepth = as.integer(max.treedepth),
                     adapt_delta = adapt.delta,
                     seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  nms <- c(paste0("psi_", model$process_info[[1]]$coef_names),
           paste0("p11_", model$process_info[[2]]$coef_names),
           paste0("p10_", model$process_info[[3]]$coef_names),
           paste0("b_",   model$process_info[[4]]$coef_names))
  colnames(draws) <- nms
  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- chains[[1L]]$epsilon

  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws)
  marg <- .tobs_fp_occu_nuts_marginal(model)
  ev_mean <- marg$eval_beta(par[lay$psi], par[lay$p11], par[lay$p10], par[lay$b])

  raw <- list(means = unname(par), vcov = cov, coef_names = nms,
              log_lik = ev_mean$log_lik, log_lik_site = ev_mean$log_lik_site,
              w1 = ev_mean$w1, converged = TRUE)
  fit <- build_fp_occu_fit(raw, model)

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
                   divergent_total = sum(divergent), sigma_beta = sigma.beta)
  fit
}
