# fp_occu_nuts.R - NUTS target for the multistate false-positive occupancy family.
#
# Flat coefficient vector theta = (beta_psi, beta_p11, beta_p10, beta_b); the
# joint log-posterior is the false-positive occupancy marginal
# (cpp_fp_occu_total_log_lik) plus weak Gaussian priors. The C++ FullGradFn
# (src/fp_occu_nuts.cpp) mirrors this R target and is cross-checked against it.

# `re_groups` > 0 appends a trailing [z_1..z_G, log_sigma_re] block (a single
# intercept RE on the occupancy (psi) arm, tulpaObs#51).
.tobs_fp_occu_nuts_layout <- function(p_psi, p_p11, p_p10, p_b, re_groups = 0L) {
  idx <- .tobs_nuts_arm_idx(c("psi", "p11", "p10", "b"),
                            c(p_psi, p_p11, p_p10, p_b))
  base <- p_psi + p_p11 + p_p10 + p_b
  out <- c(list(p_psi = p_psi, p_p11 = p_p11, p_p10 = p_p10, p_b = p_b),
           idx, list(re_groups = as.integer(re_groups)))
  if (re_groups > 0L) {
    out$z <- base + seq_len(re_groups); out$log_sigma <- base + re_groups + 1L
    out$total <- base + re_groups + 1L
  } else {
    out$z <- integer(0); out$log_sigma <- integer(0); out$total <- base
  }
  out
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

.tobs_fit_fp_occu_nuts <- function(model, sigma.beta = 10, re = NULL,
                                   n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
                                   max.treedepth = 10L, adapt.delta = 0.9,
                                   seed = 1L, verbose = FALSE) {
  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]

  # Single intercept RE on the occupancy (psi) arm (tulpaObs#51), via the shared
  # count-NUTS RE helpers. The false-positive arms (p11 / p10 / b) keep fixed
  # effects; a non-psi RE is rejected with a pointer.
  re_info <- .tobs_count_nuts_re_info(re, model, arms = c("psi", "p11"))
  if (!is.null(re_info) && re_info$arm != 0L)
    stop("fp_occu() NUTS supports a random effect on the occupancy (psi) arm ",
         "only; put the RE on the state formula.", call. = FALSE)
  n_re_groups <- if (!is.null(re_info)) re_info$n_groups else 0L
  lay <- .tobs_fp_occu_nuts_layout(ncol(X_psi), ncol(X_p11), ncol(X_p10),
                                   ncol(X_b), re_groups = n_re_groups)

  warm <- fp_occu_laplace(y = model$y_long, site_idx = model$site_idx,
                          X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
                          verbose = FALSE)
  n_base <- length(warm$means)
  V <- as.matrix(warm$vcov)
  base_metric <- if (!is.null(V) && nrow(V) == n_base && all(is.finite(diag(V))))
                   pmax(diag(V), 1e-6) else rep(1, n_base)
  init <- .tobs_count_nuts_re_init(
    list(theta0 = as.numeric(warm$means), inv_metric = base_metric), lay, re_info)
  theta0 <- init$theta0; inv_metric <- init$inv_metric

  spec <- .tobs_count_nuts_re_spec(
    list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
         X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
         n_sites = model$n_sites),
    re_info, 1.5)

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
           paste0("b_",   model$process_info[[4]]$coef_names),
           .tobs_count_nuts_re_names(re_info))
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
  fit <- .tobs_count_nuts_re_finish(fit, draws, par, cov, nms, re_info)
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
                   re_arm = if (!is.null(re_info)) re_info$arm else -1L)
  fit
}
