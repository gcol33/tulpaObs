# distance_nuts.R - NUTS target density for the binned distance-sampling family.
#
# The flat coefficient vector is
#   theta = (beta_lambda [p_lam], beta_sigma [p_sig],
#            log_shape [hazard-rate only], log_r [NB only])
# and the joint log-posterior is the distance marginal (cpp_distance_total_log_lik)
# plus weak Gaussian priors. The C++ FullGradFn (src/distance_nuts.cpp) mirrors
# this R target and is cross-checked against it; this R version is the oracle.

# Parameter layout.
.tobs_distance_nuts_layout <- function(p_lam, p_sig, hazard, is_nb) {
  d <- p_lam + p_sig + (if (hazard) 1L else 0L) + (if (is_nb) 1L else 0L)
  off <- p_lam + p_sig
  log_shape <- if (hazard) { off <- off + 1L; off } else integer(0)
  log_r     <- if (is_nb)  off + 1L else integer(0)
  list(p_lam = p_lam, p_sig = p_sig, hazard = isTRUE(hazard), is_nb = isTRUE(is_nb),
       lambda = seq_len(p_lam), sigma = p_lam + seq_len(p_sig),
       log_shape = log_shape, log_r = log_r, total = d)
}

# Per-site distance marginal closure (the NUTS oracle's data + eval_beta).
.tobs_distance_nuts_marginal <- function(model, mixture = "P", K_max = NULL) {
  y         <- model$y
  X_lambda  <- model$X_processes[[1]]
  X_sigma   <- model$X_processes[[2]]
  cutpoints <- model$cutpoints
  key_code  <- .dist_key_code(model$key)
  trans_code <- .dist_transect_code(model$transect)
  quad_order <- model$quad_order
  if (is.null(K_max)) K_max <- 3L * max(rowSums(y)) + 100L
  K_max <- as.integer(K_max)
  resolve_r <- function(r) {
    if (identical(mixture, "P")) return(Inf)
    if (is.null(r) || !is.finite(r) || r <= 0)
      stop("NB distance marginal requires a finite positive `r`.", call. = FALSE)
    as.numeric(r)
  }
  eval_beta <- function(beta_lambda, beta_sigma, eta_b = 0, r = Inf) {
    eta_lambda <- as.numeric(X_lambda %*% beta_lambda)
    eta_sigma  <- as.numeric(X_sigma  %*% beta_sigma)
    cpp_distance_total_log_lik(y, eta_lambda, eta_sigma, as.numeric(eta_b),
                               cutpoints, trans_code, key_code, K_max,
                               resolve_r(r), quad_order)
  }
  list(X_lambda = X_lambda, X_sigma = X_sigma, K_max = K_max, mixture = mixture,
       eval_beta = eval_beta)
}

# Joint log-posterior + gradient of the distance coefficient vector (R oracle).
.tobs_distance_nuts_logpost <- function(theta, marg, lay,
                                        sigma.beta = 10, sigma.shape = 1.5,
                                        sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_sigma  <- theta[lay$sigma]
  eta_b <- if (lay$hazard) theta[lay$log_shape] else 0
  r     <- if (lay$is_nb)  exp(theta[lay$log_r]) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_sigma, eta_b = eta_b, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_sigma,  ev$grad_eta_sigma)))
  if (lay$hazard) grad <- c(grad, ev$grad_eta_b)
  if (lay$is_nb)  grad <- c(grad, sum(ev$grad_theta))

  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * (sum(beta_lambda^2) + sum(beta_sigma^2))
  grad[lay$lambda] <- grad[lay$lambda] - ib2 * beta_lambda
  grad[lay$sigma]  <- grad[lay$sigma]  - ib2 * beta_sigma
  if (lay$hazard) {
    is2 <- 1 / sigma.shape^2
    lp  <- lp - 0.5 * is2 * eta_b^2
    grad[lay$log_shape] <- grad[lay$log_shape] - is2 * eta_b
  }
  if (lay$is_nb) {
    lr  <- theta[lay$log_r]; ilr2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * ilr2 * lr^2
    grad[lay$log_r] <- grad[lay$log_r] - ilr2 * lr
  }
  list(lp = lp, grad = grad)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the distance family
# ---------------------------------------------------------------------------

.tobs_fit_distance_nuts <- function(model, mixture = "poisson", K_max = NULL,
                                    sigma.beta = 10, sigma.shape = 1.5,
                                    sigma.logr = 1.5,
                                    n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
                                    max.treedepth = 10L, adapt.delta = 0.9,
                                    seed = 1L, verbose = FALSE) {
  is_nb    <- identical(mixture, "negbin")
  mix_code <- if (is_nb) "NB" else "P"
  hazard   <- identical(model$key, "hazard")
  X_lambda <- model$X_processes[[1]]
  X_sigma  <- model$X_processes[[2]]
  y        <- model$y
  p_lam <- ncol(X_lambda); p_sig <- ncol(X_sigma)
  if (is.null(K_max)) K_max <- 3L * max(rowSums(y)) + 100L
  K_max <- as.integer(K_max)
  lay <- .tobs_distance_nuts_layout(p_lam, p_sig, hazard, is_nb)

  warm <- distance_laplace(y = y, X_lambda = X_lambda, X_sigma = X_sigma,
                           cutpoints = model$cutpoints, key = model$key,
                           transect = model$transect, mixture = mix_code,
                           K_max = K_max, quad_order = model$quad_order,
                           max_iter = 100L, verbose = FALSE)
  theta0 <- c(as.numeric(warm$beta_lambda), as.numeric(warm$beta_sigma))
  if (hazard) theta0 <- c(theta0, as.numeric(warm$eta_b))
  if (is_nb)  theta0 <- c(theta0, as.numeric(warm$log_r))
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == lay$total)) pmax(diag(V), 1e-6)
                else rep(1, lay$total)

  spec <- list(y = y, X_lambda = X_lambda, X_sigma = X_sigma,
               cutpoints = as.numeric(model$cutpoints),
               transect = .dist_transect_code(model$transect),
               key = .dist_key_code(model$key), K_max = K_max,
               is_nb = is_nb, quad_order = as.integer(model$quad_order))

  run_chain <- function(ch) {
    cpp_distance_nuts(spec, theta0 = theta0,
                      sigma_beta = sigma.beta, sigma_shape = sigma.shape,
                      sigma_logr = sigma.logr, inv_metric = inv_metric,
                      n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta,
                      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("sigma_",  model$process_info[[2]]$coef_names),
           if (hazard) "log_shape", if (is_nb) "log_r")
  colnames(draws) <- nms
  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- chains[[1L]]$epsilon

  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws)
  marg <- .tobs_distance_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$sigma],
                            eta_b = if (hazard) par[lay$log_shape] else 0,
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik

  raw <- list(
    mixture = mix_code, key = model$key, transect = model$transect,
    hazard = hazard, nb = is_nb,
    beta_lambda = unname(par[lay$lambda]),
    beta_sigma  = unname(par[lay$sigma]),
    eta_b  = if (hazard) unname(par[lay$log_shape]) else NA_real_,
    shape  = if (hazard) exp(unname(par[lay$log_shape])) else NA_real_,
    log_r  = if (is_nb) unname(par[lay$log_r]) else NA_real_,
    r      = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
    vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max,
    mean_N = warm$mean_N, var_N = warm$var_N, p_det = warm$p_det,
    boundary_weight = warm$boundary_weight)
  fit <- build_distance_fit(raw, model)

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
                   is_nb = is_nb, hazard = hazard, K_max = K_max,
                   sigma_beta = sigma.beta, sigma_shape = sigma.shape,
                   sigma_logr = sigma.logr)
  fit
}
