# fp_occu_nuts.R - NUTS target for the multistate false-positive occupancy family.
#
# Flat coefficient vector theta = (beta_psi, beta_p11, beta_p10, beta_b); the
# joint log-posterior is the false-positive occupancy marginal
# (cpp_fp_occu_total_log_lik) plus weak Gaussian priors. The C++ FullGradFn
# (src/fp_occu_nuts.cpp) mirrors this R target and is cross-checked against it.

# `re_groups` > 0 appends a trailing [z_1..z_G, log_sigma_re] block (a single
# intercept RE on the occupancy (psi) arm).
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

# `sigma.logr` is the prior SD on the RE log-SD (the shared count-NUTS RE
# block).
.tobs_fit_fp_occu_nuts <- function(model, sigma.beta = NULL, sigma.logr = NULL,
                                   re = NULL,
                                   n.iter = NULL, n.warmup = NULL, n.chains = NULL, n.thin = NULL,
                                   n.threads = NULL,
                                   max.treedepth = NULL, adapt.delta = NULL,
                                   seed = NULL, verbose = FALSE) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]

  # Single intercept RE on the occupancy (psi) arm, via the shared count-NUTS RE
  # helpers. The false-positive arms (p11 / p10 / b) keep fixed effects; a
  # non-psi RE is rejected with a pointer.
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
  run <- .tobs_count_nuts_front_door(
    cpp_fp_occu_nuts,
    spec = list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
                X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
                n_sites = model$n_sites),
    priors = list(sigma_beta = sigma.beta),
    theta0 = warm$means, V = warm$vcov,
    nms = c(paste0("psi_", model$process_info[[1]]$coef_names),
            paste0("p11_", model$process_info[[2]]$coef_names),
            paste0("p10_", model$process_info[[3]]$coef_names),
            paste0("b_",   model$process_info[[4]]$coef_names)),
    re_info = re_info, sigma.logr = sigma.logr,
    sampler = .tobs_sampler_control_snapshot(environment()), verbose = verbose)
  par <- run$par; cov <- run$cov; nms <- run$nms

  marg <- .tobs_fp_occu_nuts_marginal(model)
  ev_mean <- marg$eval_beta(par[lay$psi], par[lay$p11], par[lay$p10], par[lay$b])

  raw <- list(means = unname(par), vcov = cov, coef_names = nms,
              log_lik = ev_mean$log_lik, log_lik_site = ev_mean$log_lik_site,
              w1 = ev_mean$w1, converged = TRUE)
  fit <- build_fp_occu_fit(raw, model)

  .tobs_count_nuts_attach(
    fit, run, ev_mean$log_lik,
    extra = list(sigma_beta = sigma.beta, sigma_logr = sigma.logr,
                 re_arm = if (!is.null(re_info)) re_info$arm else -1L))
}
