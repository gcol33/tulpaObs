# abun_nuts.R - NUTS target density for the single-species N-mixture (abun()).
#
# The Laplace fit (abun.R -> nmix_laplace) sums the latent abundance N out in
# closed form and returns a Gaussian observed-Fisher posterior over the
# coefficients c(beta_lambda, beta_p[, log_r]). NUTS instead samples the exact
# marginal posterior of those coefficients, which gives calibrated (non-Gaussian)
# intervals and the per-draw pointwise likelihood WAIC / LOO need.
#
# There is no latent field, no random effect, and no community covariance here,
# so -- unlike the spatial-factor community target (ms_occu_cover_spatial_nuts.R)
# -- the parameter vector is just the flat coefficient block. The joint
# log-posterior is
#
#   log p(theta | y) = sum_i log m_i(theta)              # per-site N-mixture marginal
#                      - 0.5 ||beta||^2 / sigma.beta^2    # weak Gaussian coef priors
#                      - 0.5  log_r^2  / sigma.logr^2     # (NB only)
#
# where m_i is the Royle (2004) per-site marginal exposed by nmix_site_marginal()
# (the same kernel the Laplace fit and the AGHQ RE path use); it already
# differentiates through both arms, returning grad_eta_lambda, grad_eta_p, and
# (NB) grad_theta = d log m_i / d log_r. The coefficient gradient is the design-
# sandwiched eta-gradient, so the whole target reuses the existing kernel with no
# new likelihood math. The C++ FullGradFn (src/abun_nuts.cpp) mirrors this R
# target byte-for-byte and is cross-checked against it before driving tulpa's NUTS
# engine; this R version is the oracle.

# Parameter layout: (beta_lambda [p_lam], beta_p [p_p], [log_r under NB],
# [z_1..z_G, log_sigma_re] when a single intercept RE is present). `re_groups`
# is 0 for no RE.
.tobs_abun_nuts_layout <- function(p_lam, p_p, is_nb, re_groups = 0L) {
  base <- p_lam + p_p + (if (is_nb) 1L else 0L)
  out <- list(
    p_lam   = p_lam,
    p_p     = p_p,
    is_nb   = isTRUE(is_nb),
    lambda  = seq_len(p_lam),
    p       = p_lam + seq_len(p_p),
    log_r   = if (is_nb) p_lam + p_p + 1L else integer(0),
    re_groups = as.integer(re_groups)
  )
  if (re_groups > 0L) {
    out$z         <- base + seq_len(re_groups)
    out$log_sigma <- base + re_groups + 1L
    out$total     <- base + re_groups + 1L
  } else {
    out$z <- integer(0); out$log_sigma <- integer(0); out$total <- base
  }
  out
}

# Joint log-posterior + gradient of the N-mixture coefficient vector. `marg` is a
# `nmix_marginal` object (nmix_site_marginal()); `lay` the layout above. Weak
# Gaussian priors N(0, sigma.beta^2) on every coefficient and N(0, sigma.logr^2)
# on log_r keep the ridge off the boundary without materially shifting the
# data-dominated optimum, matching the other Laplace paths. Returns
# list(lp, grad) over the packed coordinates.
.tobs_abun_nuts_logpost <- function(theta, marg, lay,
                                    sigma.beta = 10, sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_p      <- theta[lay$p]
  r <- if (lay$is_nb) exp(theta[lay$log_r]) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_p, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_p,      ev$grad_eta_p)))
  if (lay$is_nb) grad <- c(grad, sum(ev$grad_theta))

  # Weak Gaussian coefficient priors (data-dominated ridge).
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

# Build the per-site marginal kernel for a bound nmix model, ready to feed the
# NUTS target. Mirrors the long form .tobs_build_abun() already produced
# (y_long / site_idx, the autoscaled designs), so the NUTS target sees exactly
# the data the Laplace fit did.
.tobs_abun_nuts_marginal <- function(model, mixture = "P", K_max = NULL) {
  nmix_site_marginal(
    y        = model$y_long,
    site_idx = model$site_idx,
    X_lambda = model$X_processes[[1]],
    X_p      = model$X_processes[[2]],
    mixture  = mixture,
    K_max    = K_max)
}

# Warm-start the sampler at the Laplace mode and a diagonal inverse-metric from
# the Laplace curvature (the #67 recipe: init = mode, inv_metric = 1 / diag of
# the observed-information Hessian, floored). `raw` is a nmix_laplace() fit.
.tobs_abun_nuts_pack_init <- function(raw, lay) {
  theta0 <- c(as.numeric(raw$beta_lambda), as.numeric(raw$beta_p))
  if (lay$is_nb) theta0 <- c(theta0, as.numeric(raw$log_r))
  V <- as.matrix(raw$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == lay$total)) {
    pmax(diag(V), 1e-6)              # metric = posterior-scale (Laplace vcov)
  } else rep(1, lay$total)
  list(theta0 = theta0, inv_metric = inv_metric)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the single-species N-mixture
# ---------------------------------------------------------------------------

# Sample the exact coefficient posterior of a non-spatial N-mixture via tulpa's
# NUTS engine and the in-tree C++ FullGradFn (cpp_abun_nuts), warm-started at the
# Laplace mode with a diagonal Laplace metric, then package the draws into the
# build_nmix_fit shape so coef / vcov / confint / predict / WAIC read the NUTS
# posterior. Operates on the (autoscaled) `model`; the caller .tobs_fit_model()
# unscales the means / draws / vcov back to the natural coefficient scale, so the
# real NUTS draws land in `fit$draws` on the natural scale (WAIC reads them).
# `fit$nuts` carries the sampler diagnostics only (no draws -- fit$draws is the
# single, unscaled posterior copy).
.tobs_fit_abun_nuts <- function(model, mixture = "poisson", K_max = NULL,
                                sigma.beta = 10, sigma.logr = 1.5, re = NULL,
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
  if (is.null(K_max)) K_max <- max(y_long) + 100L
  K_max <- as.integer(K_max)

  # Single intercept random effect on one arm (tulpaObs#51). Reuse the Laplace
  # path's arm split; v1 NUTS supports one grouping factor, intercept only, on
  # the abundance OR detection arm. Slopes / correlated / multi-term / both-arm
  # RE stay on the AGHQ Laplace path (richer than the single-block NUTS target).
  re_arm <- -1L; re_group <- integer(0); n_re_groups <- 0L; re_label <- NULL
  if (!is.null(re) && length(re) > 0L) {
    arms <- .tobs_nmix_re_split_arms(re, model)
    if (length(arms$lambda) && length(arms$p))
      stop("method = \"nuts\" for abun() with a random effect supports the RE on ",
           "ONE arm; put it on lambda OR p, or use method = \"laplace\".",
           call. = FALSE)
    design <- if (length(arms$lambda)) arms$lambda else arms$p
    re_arm <- if (length(arms$lambda)) 0L else 1L
    if (length(design) != 1L || design[[1L]]$n_coefs != 1L ||
        !isTRUE(design[[1L]]$has_intercept))
      stop("method = \"nuts\" for abun() supports a single intercept random ",
           "effect (1|g) on one arm; random slopes / multiple grouping factors ",
           "fit under method = \"laplace\" (AGHQ).", call. = FALSE)
    d1          <- design[[1L]]
    re_group    <- as.integer(d1$idx)
    n_re_groups <- as.integer(d1$n_groups)
    re_label    <- d1$group_label %||% "g1"
  }
  has_re <- re_arm >= 0L
  lay <- .tobs_abun_nuts_layout(p_lam, p_p, is_nb, re_groups = n_re_groups)

  # Warm start at the Laplace mode (+ diagonal Laplace metric from its vcov).
  warm <- nmix_laplace(y = y_long, site_idx = site_idx, X_lambda = X_lambda,
                       X_p = X_p, mixture = mix_code, K_max = K_max,
                       max_iter = 100L, verbose = FALSE)
  init <- .tobs_abun_nuts_pack_init(warm, lay)
  if (has_re) {
    # z warm-started at 0 (no group deviation), log_sigma_re at log(0.5); the
    # metric is unit on z (standard-normal non-centred prior scale) and 0.25 on
    # log_sigma_re.
    init$theta0     <- c(init$theta0, rep(0, n_re_groups), log(0.5))
    init$inv_metric <- c(init$inv_metric[seq_len(lay$total - n_re_groups - 1L)],
                         rep(1, n_re_groups), 0.25)
  }

  spec <- list(y = as.integer(y_long), site_idx = as.integer(site_idx),
               X_lambda = X_lambda, X_p = X_p,
               n_sites = model$n_sites, K_max = K_max, is_nb = is_nb,
               re_arm = re_arm)
  if (has_re) {
    spec$re_group     <- re_group
    spec$n_re_groups  <- n_re_groups
    spec$sigma_re_lsd <- sigma.logr
  }

  run_chain <- function(ch) {
    cpp_abun_nuts(spec, theta0 = init$theta0,
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
  arm_tag <- if (re_arm == 1L) "p" else "lambda"
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           if (is_nb) "log_r",
           if (has_re) c(paste0("re_", re_label, "_z", seq_len(n_re_groups)),
                         paste0("log_sigma_", arm_tag, "_", re_label)))
  colnames(draws) <- nms
  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- chains[[1L]]$epsilon

  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws)
  # Data log-likelihood at the posterior mean (scale-invariant), so logLik() on
  # the NUTS fit matches the laplace-path convention (mean(log_prob)).
  marg <- .tobs_abun_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$p],
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik

  raw <- list(
    mixture     = mix_code,
    beta_lambda = unname(par[lay$lambda]),
    beta_p      = unname(par[lay$p]),
    log_r       = if (is_nb) unname(par[lay$log_r]) else NA_real_,
    r           = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
    vcov        = cov,
    log_lik     = ll_mean,
    converged   = TRUE,
    K_max       = K_max,
    mean_N      = warm$mean_N, var_N = warm$var_N,
    boundary_weight = warm$boundary_weight)
  fit <- build_nmix_fit(raw, model, spatial = NULL)

  # Replace the moment-matched MVN draws + NA sampler diagnostics with the actual
  # NUTS posterior and real diagnostics. fit$draws is the single posterior copy;
  # .tobs_fit_model() unscales it to the natural coefficient scale. With an RE
  # block present, means/sds/vcov carry the trailing RE coordinates too (column
  # order [lambda, p, (log_r), z_1..z_G, log_sigma_re] matches the draws, so the
  # per-process unscaler leaves the RE tail untouched, like log_r).
  n_draws <- nrow(draws)
  fit$draws       <- draws
  if (has_re) {
    fit$means <- par
    fit$sds   <- sqrt(pmax(diag(cov), 0)); names(fit$sds) <- nms
    fit$vcov  <- cov
    # RE summary on the natural scale: sigma_re = exp(log_sigma_re), per-group
    # BLUP b_g = sigma_re * z_g (posterior means over draws).
    ls_col  <- paste0("log_sigma_", arm_tag, "_", re_label)
    z_cols  <- paste0("re_", re_label, "_z", seq_len(n_re_groups))
    sig_dr  <- exp(draws[, ls_col])
    blup_dr <- sig_dr * draws[, z_cols, drop = FALSE]
    fit$re <- list(
      arm = arm_tag, group_label = re_label, n_groups = n_re_groups,
      sigma = mean(sig_dr), sigma_sd = stats::sd(sig_dr),
      blup = colMeans(blup_dr),
      blup_sd = apply(blup_dr, 2L, stats::sd))
  }
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
                   is_nb = is_nb, K_max = K_max, re_arm = re_arm,
                   sigma_beta = sigma.beta, sigma_logr = sigma.logr)
  fit
}
