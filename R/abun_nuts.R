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
# Shared single-intercept RE wiring for the count-marginal NUTS families
# (abun, removal): one grouping factor, intercept only, on the abundance OR
# detection arm. The C++ side (marginal_count_nuts.h) carries the non-centered
# per-site offset; these helpers resolve the RE term, thread it into the spec /
# warm-start / draw names, and surface fit$re. Slopes / correlated / multi-term
# / both-arm RE stay on the AGHQ Laplace path.
# ---------------------------------------------------------------------------

# Resolve a single intercept RE from the formula `re` terms. NULL = no RE.
# `arms` names the two predictor blocks the term may sit on (state arm first,
# detection arm second) so the count families (abun/removal: lambda/p), the
# open N-mixture and distance (lambda only) and the false-positive occupancy
# family (psi/p11) share one resolver; the returned `arm` is 0 for the first,
# 1 for the second.
.tobs_count_nuts_re_info <- function(re, model, arms = c("lambda", "p")) {
  if (is.null(re) || length(re) == 0L) return(NULL)
  if (inherits(re, "tobs_re")) re <- list(re)
  split <- .tobs_re_split_two_arms(
    re, model, arms[1L], arms[2L],
    sprintf(paste0("method = \"nuts\" with a random effect supports the RE on ",
                   "ONE arm; put it on %s OR %s, or use method = \"laplace\"."),
            arms[1L], arms[2L]))
  a1 <- split[[arms[1L]]]; a2 <- split[[arms[2L]]]
  if (length(a1) && length(a2))
    stop(sprintf(paste0("method = \"nuts\" with a random effect supports the RE ",
                        "on ONE arm; put it on %s OR %s, or use ",
                        "method = \"laplace\"."), arms[1L], arms[2L]),
         call. = FALSE)
  design <- if (length(a1)) a1 else a2
  if (length(design) != 1L || design[[1L]]$n_coefs != 1L ||
      !isTRUE(design[[1L]]$has_intercept))
    stop("method = \"nuts\" supports a single intercept random effect (1|g) on ",
         "one arm; random slopes / multiple grouping factors fit under ",
         "method = \"laplace\" (AGHQ).", call. = FALSE)
  list(arm = if (length(a1)) 0L else 1L,
       arm_tag = if (length(a1)) arms[1L] else arms[2L],
       group = as.integer(design[[1L]]$idx),
       n_groups = as.integer(design[[1L]]$n_groups),
       label = design[[1L]]$group_label %||% "g1")
}

# Merge the RE block fields into the NUTS spec list (re_arm = -1 -> no RE).
.tobs_count_nuts_re_spec <- function(spec, re_info, sigma.logr) {
  if (is.null(re_info)) { spec$re_arm <- -1L; return(spec) }
  spec$re_arm       <- re_info$arm
  spec$re_group     <- re_info$group
  spec$n_re_groups  <- re_info$n_groups
  spec$sigma_re_lsd <- sigma.logr
  spec
}

# Append z (warm-started at 0, unit metric) + log_sigma_re (log(0.5), 0.25
# metric) to the warm-start init.
.tobs_count_nuts_re_init <- function(init, lay, re_info) {
  if (is.null(re_info)) return(init)
  G <- re_info$n_groups
  init$theta0     <- c(init$theta0, rep(0, G), log(0.5))
  init$inv_metric <- c(init$inv_metric[seq_len(lay$total - G - 1L)],
                       rep(1, G), 0.25)
  init
}

# Draw-column names for the RE block (z_1..z_G + log_sigma_<arm>_<label>).
.tobs_count_nuts_re_names <- function(re_info) {
  if (is.null(re_info)) return(character(0))
  c(paste0("re_", re_info$label, "_z", seq_len(re_info$n_groups)),
    paste0("log_sigma_", re_info$arm_tag, "_", re_info$label))
}

# With an RE present, set means/sds/vcov to the full coordinate set (so they
# align with the RE draws columns and the per-process unscaler leaves the RE
# tail untouched) and attach fit$re (sigma + per-group BLUPs on the natural
# scale, b_g = sigma_re * z_g).
.tobs_count_nuts_re_finish <- function(fit, draws, par, cov, nms, re_info) {
  if (is.null(re_info)) return(fit)
  fit$means <- par
  fit$sds   <- sqrt(pmax(diag(cov), 0)); names(fit$sds) <- nms
  fit$vcov  <- cov
  ls_col  <- paste0("log_sigma_", re_info$arm_tag, "_", re_info$label)
  z_cols  <- paste0("re_", re_info$label, "_z", seq_len(re_info$n_groups))
  sig_dr  <- exp(draws[, ls_col])
  blup_dr <- sig_dr * draws[, z_cols, drop = FALSE]
  fit$re <- list(arm = re_info$arm_tag, group_label = re_info$label,
                 n_groups = re_info$n_groups,
                 sigma = mean(sig_dr), sigma_sd = stats::sd(sig_dr),
                 blup = colMeans(blup_dr), blup_sd = apply(blup_dr, 2L, stats::sd))
  fit
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

  # Single intercept random effect on one arm (tulpaObs#51), via the shared
  # count-NUTS RE helpers.
  re_info <- .tobs_count_nuts_re_info(re, model)
  has_re  <- !is.null(re_info)
  n_re_groups <- if (has_re) re_info$n_groups else 0L
  lay <- .tobs_abun_nuts_layout(p_lam, p_p, is_nb, re_groups = n_re_groups)

  # Warm start at the Laplace mode (+ diagonal Laplace metric from its vcov).
  warm <- nmix_laplace(y = y_long, site_idx = site_idx, X_lambda = X_lambda,
                       X_p = X_p, mixture = mix_code, K_max = K_max,
                       max_iter = 100L, verbose = FALSE)
  init <- .tobs_count_nuts_re_init(.tobs_abun_nuts_pack_init(warm, lay), lay, re_info)

  spec <- .tobs_count_nuts_re_spec(
    list(y = as.integer(y_long), site_idx = as.integer(site_idx),
         X_lambda = X_lambda, X_p = X_p,
         n_sites = model$n_sites, K_max = K_max, is_nb = is_nb),
    re_info, sigma.logr)

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
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           if (is_nb) "log_r",
           .tobs_count_nuts_re_names(re_info))
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
  fit <- .tobs_count_nuts_re_finish(fit, draws, par, cov, nms, re_info)
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
                   re_arm = if (has_re) re_info$arm else -1L,
                   sigma_beta = sigma.beta, sigma_logr = sigma.logr)
  fit
}


# Areal-spatial N-mixture via NUTS (gcol33/tulpaObs#51, #71): a FIXED-HYPER non-
# centered field on the abundance arm. The field precision (tau, rho) is fixed at
# the nested-Laplace posterior mean; the whitened field raw ~ N(0, I) with
# z = L %*% raw is sampled jointly with the coefficients (the tulpa#87 fixed-hyper
# pattern -- avoids the field-hyperparameter funnel and the log|Q(rho)| gradient).
# car_proper uses the square inverse Cholesky of tau Q(rho); the intrinsic
# icar / bym2 fields use the SUM-TO-ZERO reparameterisation (#71): L drops the
# precision null-space (constant) direction, so z is automatically centred and
# the geometry is well conditioned (no flat field-mean direction maxing the tree
# depth). bym2 combines the structured (centred eigen-loading, scaled by
# sigma sqrt(rho / scale)) and unstructured (iid, sigma sqrt(1 - rho)) blocks into
# one loading. Spatial XOR RE. Poisson or NB.
.tobs_fit_abun_nuts_spatial <- function(model, spatial, mixture = "poisson",
                                        K_max = NULL, sigma.beta = 10, sigma.logr = 1.5,
                                        n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
                                        max.treedepth = 10L, adapt.delta = 0.9,
                                        seed = 1L, verbose = FALSE) {
  .tobs_reject_weighted_spatial(spatial, "abun NUTS abundance spatial")
  if (!spatial$type %in% c("icar", "car_proper", "bym2"))
    stop(sprintf(paste0("abun() NUTS + areal spatial supports icar() / car_proper() / ",
                        "bym2() on the abundance arm; got '%s'. (tulpaObs#71)"),
                 spatial$type), call. = FALSE)
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites) {
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for abun NUTS."),
                 spatial$n_units, n_sites), call. = FALSE)
  }
  is_nb <- identical(mixture, "negbin")
  mix_code <- if (is_nb) "NB" else "P"
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  y_long <- model$y_long; site_idx <- model$site_idx
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  if (is.null(K_max)) K_max <- max(y_long) + 100L
  K_max <- as.integer(K_max)
  adj <- as.matrix(spatial$graph)
  csr <- if (!is.null(spatial$adj_row_ptr))
    list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
         n_neighbors = spatial$n_neighbors) else adjacency_to_csr(spatial$graph)

  # Fixed hyper + warm coefficients from the nested-Laplace areal fit, and the
  # whitened-field loading L (square for car_proper, sum-to-zero for icar / bym2).
  common <- list(y = y_long, site_idx = site_idx, map_site_to_unit = seq_len(n_sites),
                 X_lambda = X_lambda, X_p = X_p, adj_row_ptr = csr$row_ptr,
                 adj_col_idx = csr$col_idx, n_neighbors = csr$n_neighbors,
                 n_spatial = n_sites, mixture = mix_code, K_max = K_max,
                 max_iter = 100L, tol = 1e-6, verbose = FALSE)
  nl <- switch(spatial$type,
    icar       = do.call(nmix_laplace_icar, common),
    car_proper = do.call(nmix_laplace_car_proper, common),
    bym2       = do.call(nmix_laplace_bym2, c(common,
                   list(scale_factor = spatial$scale_factor %||%
                          compute_bym2_scale(spatial$graph)))))

  # Whitened-field loading + fixed hyper (shared single source of truth, #71/#113).
  fl <- .tobs_nuts_field_loading(adj, spatial$type, n_sites,
                                 tau = nl$tau_mean, rho = nl$rho_mean,
                                 sigma = nl$sigma_mean,
                                 scale_factor = spatial$scale_factor)
  field_load <- fl$field_load; n_raw <- fl$n_raw

  spec <- list(y = as.integer(y_long), site_idx = as.integer(site_idx),
               X_lambda = X_lambda, X_p = X_p, n_sites = n_sites, K_max = K_max,
               is_nb = is_nb, n_field_units = n_sites,
               field_map = seq_len(n_sites), field_load = field_load)

  n_base <- p_lam + p_p + if (is_nb) 1L else 0L
  beta0 <- c(nl$beta_lambda_mean, nl$beta_p_mean)
  if (is_nb && is.finite(nl$r_mean %||% NA_real_)) beta0 <- c(beta0, log(nl$r_mean))
  else if (is_nb) beta0 <- c(beta0, log(2))
  theta0 <- c(beta0, numeric(n_raw))
  inv_metric <- c(rep(0.1, n_base), rep(1, n_raw))

  run_chain <- function(ch)
    cpp_abun_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logr = sigma.logr,
                  inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                  n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
                  adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                  verbose = isTRUE(verbose))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           if (is_nb) "log_r", paste0("raw_", seq_len(n_raw)))
  run <- .tobs_nuts_field_draws(run_chain, n.chains, nms, n_base, n_raw, field_load)
  par <- run$par; cov <- run$cov

  marg <- .tobs_abun_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  lay  <- .tobs_abun_nuts_layout(p_lam, p_p, is_nb)
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$p],
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik
  raw_fit <- list(mixture = mix_code,
                  beta_lambda = unname(par[lay$lambda]), beta_p = unname(par[lay$p]),
                  log_r = if (is_nb) unname(par[lay$log_r]) else NA_real_,
                  r = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
                  vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max)
  fit <- build_nmix_fit(raw_fit, model, spatial = spatial)
  .tobs_nuts_field_attach(fit, run, ll_mean, n.chains,
                          prior_type = spatial$type, fl = fl)
}
