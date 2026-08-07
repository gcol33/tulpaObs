# dyn_abun_nuts.R - NUTS target for the Dail-Madsen open N-mixture family.
#
# Flat coefficient vector theta = (beta_lambda, beta_p, beta_omega, beta_gamma);
# the joint log-posterior is the forward marginal (cpp_dyn_abun_total_log_lik)
# plus weak Gaussian priors. The C++ FullGradFn (src/dyn_abun_nuts.cpp) mirrors
# this R target and is cross-checked against it.

# `use_nb` appends a single trailing log r coordinate (NB initial abundance);
# `re_groups` > 0 then appends [z_1..z_G, log_sigma_re] (a single intercept RE on
# the initial-abundance arm, tulpaObs#51).
.tobs_dyn_abun_nuts_layout <- function(p_lam, p_p, p_om, p_gm, use_nb = FALSE,
                                       re_groups = 0L) {
  idx <- .tobs_nuts_arm_idx(c("lambda", "p", "omega", "gamma"),
                            c(p_lam, p_p, p_om, p_gm))
  base <- p_lam + p_p + p_om + p_gm
  out <- c(list(p_lam = p_lam, p_p = p_p, p_om = p_om, p_gm = p_gm),
           idx, list(use_nb = use_nb))
  if (use_nb) { out <- c(out, list(logr = base + 1L)); base <- base + 1L }
  else          out <- c(out, list(logr = NA_integer_))
  if (re_groups > 0L) {
    out <- c(out, list(re_groups = as.integer(re_groups),
                       z = base + seq_len(re_groups),
                       log_sigma = base + re_groups + 1L,
                       total = base + re_groups + 1L))
  } else {
    out <- c(out, list(re_groups = 0L, z = integer(0), log_sigma = integer(0),
                       total = base))
  }
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

.tobs_fit_dyn_abun_nuts <- function(model, sigma.beta = NULL, re = NULL,
                                    n.iter = NULL, n.warmup = NULL, n.chains = NULL,
                                    max.treedepth = NULL, adapt.delta = NULL,
                                    seed = NULL, verbose = FALSE) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_omega  <- model$X_processes[[3]]; X_gamma <- model$X_processes[[4]]
  use_nb <- identical(model$mixture %||% "poisson", "negbin")

  # Single intercept RE on the initial-abundance (lambda, tulpaObs#51) OR the
  # detection (p, tulpaObs#82) arm, via the shared count-NUTS RE helpers. The
  # offset is non-centered and routed to eta_lambda or eta_p in the C++ eval; the
  # survival / recruitment arms never carry structured terms (rejected upstream).
  re_info <- .tobs_count_nuts_re_info(re, model)
  n_re_groups <- if (!is.null(re_info)) re_info$n_groups else 0L
  lay <- .tobs_dyn_abun_nuts_layout(ncol(X_lambda), ncol(X_p), ncol(X_omega),
                                    ncol(X_gamma), use_nb = use_nb,
                                    re_groups = n_re_groups)

  warm <- dyn_abun_laplace(
    y_flat = model$y_flat, n_sites = model$n_sites, T = model$n_seasons,
    J = model$max_visits, K_max = model$K_max,
    X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
    mixture = model$mixture %||% "poisson", verbose = FALSE)
  V <- as.matrix(warm$vcov)
  n_base <- length(warm$means)
  base_metric <- if (!is.null(V) && nrow(V) == n_base && all(is.finite(diag(V))))
                   pmax(diag(V), 1e-6) else rep(1, n_base)
  init <- .tobs_count_nuts_re_init(
    list(theta0 = as.numeric(warm$means), inv_metric = base_metric), lay, re_info)
  theta0 <- init$theta0; inv_metric <- init$inv_metric

  spec <- .tobs_count_nuts_re_spec(
    list(y = as.integer(model$y_flat), n_sites = model$n_sites,
         T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
         X_lambda = X_lambda, X_p = X_p, X_omega = X_omega, X_gamma = X_gamma,
         use_nb = use_nb),
    re_info, 1.5)

  run_chain <- function(ch) {
    cpp_dyn_abun_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      inv_metric = inv_metric,
                      n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta,
                      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           paste0("omega_",  model$process_info[[3]]$coef_names),
           paste0("gamma_",  model$process_info[[4]]$coef_names))
  if (use_nb) nms <- c(nms, "log_r")
  nms <- c(nms, .tobs_count_nuts_re_names(re_info))
  run <- .tobs_count_nuts_run(run_chain, n.chains, nms)
  par <- run$par; cov <- run$cov

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

  .tobs_count_nuts_attach(
    fit, run, ev_mean$log_lik, n.chains, re_info,
    extra = list(sigma_beta = sigma.beta, K_max = model$K_max))
}
