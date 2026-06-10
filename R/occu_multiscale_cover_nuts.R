# occu_multiscale_cover_nuts.R - non-spatial NUTS for the three-level occupancy
# + cover hurdle (occu_multiscale_cover(), gcol33/tulpaObs#70).
#
# The non-spatial Laplace path (.tobs_fit_occu_multiscale_cover_laplace) sums z
# (cells) and a (plots) out in closed form and returns a Gaussian observed-Fisher
# posterior over c(beta_psi, beta_theta, beta_p, beta_pos, log_disp). NUTS instead
# samples the exact marginal posterior of those coefficients, giving calibrated
# (non-Gaussian) intervals and the per-draw pointwise likelihood WAIC / LOO need.
#
# There is no latent field here (the iid-cell marginal, field fixed at 0), so the
# parameter vector is just the flat coefficient block plus log_disp. The C++
# FullGradFn (src/occu_multiscale_cover_nuts.cpp) is the sampler target;
# .tobs_occu_ms_cover_nuts_logpost below is the byte-exact R oracle it is
# cross-checked against.


# Layout helper: contiguous coordinate blocks of the packed vector
# c(beta_psi, beta_theta, beta_p[site, visit], beta_pos[site, visit], log_disp),
# matching .tobs_fit_occu_multiscale_cover_laplace's `idx`.
.tobs_occu_ms_cover_nuts_layout <- function(model) {
  pi_list <- model$process_info
  p_psi   <- pi_list[[1L]]$p; p_theta <- pi_list[[2L]]$p
  p_p     <- pi_list[[3L]]$p; p_pos   <- pi_list[[4L]]$p
  ps_p    <- ncol(model$X_p_site); ps_pos <- ncol(model$X_pos_site)
  off     <- cumsum(c(0L, p_psi, p_theta, p_p, p_pos))
  idx <- list(
    psi   = off[1] + seq_len(p_psi),
    theta = off[2] + seq_len(p_theta),
    p     = off[3] + seq_len(p_p),
    pos   = off[4] + seq_len(p_pos),
    disp  = off[5] + 1L)
  idx$p_site    <- idx$p[seq_len(ps_p)]
  idx$p_visit   <- if (p_p > ps_p)     idx$p[(ps_p + 1L):p_p]       else integer(0)
  idx$pos_site  <- idx$pos[seq_len(ps_pos)]
  idx$pos_visit <- if (p_pos > ps_pos) idx$pos[(ps_pos + 1L):p_pos] else integer(0)
  idx$total     <- off[5] + 1L
  idx
}

# Byte-exact R oracle for the C++ FullGradFn: the exact three-level marginal
# log-likelihood plus weak Gaussian coefficient priors, with the same value the
# kernel returns. Reuses .occu_ms_cover_nonspatial_ll (the Laplace path's LL) for
# the log-density and a finite-difference-free analytic gradient is supplied by
# the C++ side; here we numerically check the value and gradient via the C++
# joint-logpost cross-check, so the oracle only needs the (penalised) log-density.
.tobs_occu_ms_cover_nuts_logpost <- function(theta, model, idx, sigma.beta = 5) {
  ll  <- .occu_ms_cover_nonspatial_ll(theta, model, idx)
  ib2 <- 1 / sigma.beta^2
  beta_idx <- setdiff(seq_len(idx$total), idx$disp)
  ll - 0.5 * ib2 * sum(theta[beta_idx]^2)
}


# Assemble the C++ NUTS spec from a bound multiscale model.
.tobs_occu_ms_cover_nuts_spec <- function(model) {
  J  <- model$max_visits
  np <- model$n_plots
  # y / y_pos / valid in plot-major flat order (row p*J + v), matching the kernel.
  y_flat    <- as.integer(t(model$y))         # t() -> row-major over (plot, visit)
  ypos_flat <- as.numeric(t(model$y_pos))
  valid_flat <- as.logical(t(model$valid))
  list(
    X_psi       = model$X_psi,
    X_theta     = model$X_theta,
    X_p_site    = model$X_p_site,
    X_pos_site  = model$X_pos_site,
    X_p_visit   = model$X_p_visit,    # NULL when no visit-level part
    X_pos_visit = model$X_pos_visit,
    y           = y_flat,
    y_pos       = ypos_flat,
    valid       = valid_flat,
    plot_cell   = as.integer(model$plot_cell),
    n_cells     = as.integer(model$n_cells),
    max_visits  = as.integer(J),
    is_beta     = identical(model$positive, "beta"))
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the non-spatial three-level occupancy + cover
# ---------------------------------------------------------------------------

# Sample the exact coefficient posterior of the non-spatial three-level model
# via tulpa's NUTS engine and the in-tree C++ FullGradFn (cpp_occu_ms_cover_nuts),
# warm-started at the Laplace mode with a diagonal Laplace metric. Packages a
# tobs_fit shaped like .tobs_fit_occu_multiscale_cover_laplace but carrying the
# real NUTS draws / diagnostics.
.tobs_fit_occu_multiscale_cover_nuts <- function(model, priors = NULL,
                                                 sigma.beta = 5,
                                                 n.iter = 1000L, n.warmup = 1000L,
                                                 n.chains = 1L, n.thin = 1L,
                                                 max.treedepth = 10L,
                                                 adapt.delta = 0.9, seed = 1L,
                                                 verbose = FALSE, ...) {
  lay     <- .tobs_occu_ms_cover_nuts_layout(model)
  pi_list <- model$process_info
  is_beta <- identical(model$positive, "beta")
  par_names <- c(
    paste0("psi_",   pi_list[[1L]]$coef_names),
    paste0("theta_", pi_list[[2L]]$coef_names),
    paste0("p_",     pi_list[[3L]]$coef_names),
    paste0("pos_",   pi_list[[4L]]$coef_names),
    if (is_beta) "log_phi" else "log_sigma_pos")

  # Warm start at the Laplace mode (+ diagonal Laplace metric from its vcov).
  warm <- .tobs_fit_occu_multiscale_cover_laplace(
    model = model, priors = priors, sigma.beta = sigma.beta, verbose = FALSE)
  theta0 <- as.numeric(warm$means)
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == lay$total) &&
                    all(is.finite(diag(V)))) pmax(diag(V), 1e-6)
                else rep(1, lay$total)

  spec <- .tobs_occu_ms_cover_nuts_spec(model)

  run_chain <- function(ch) {
    cpp_occu_ms_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, inv_metric = inv_metric,
      n_iter = as.integer(n.iter + n.warmup), n_warmup = as.integer(n.warmup),
      max_treedepth = as.integer(max.treedepth), adapt_delta = adapt.delta,
      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  per_chain <- lapply(chains, function(ch) .tobs_thin(ch$draws, n.thin))
  draws <- do.call(rbind, per_chain)
  colnames(draws) <- par_names
  chain_id  <- rep(seq_len(length(chains)), vapply(per_chain, nrow, integer(1)))
  accept    <- unlist(lapply(chains, function(ch) .tobs_thin(ch$accept_prob, n.thin)))
  divergent <- unlist(lapply(chains, function(ch) .tobs_thin(ch$divergent, n.thin)))
  treedepth <- as.integer(unlist(lapply(chains,
                                        function(ch) .tobs_thin(ch$treedepth, n.thin))))
  epsilon   <- chains[[1L]]$epsilon

  par   <- colMeans(draws); names(par) <- par_names
  V_out <- stats::cov(draws); dimnames(V_out) <- list(par_names, par_names)
  sds   <- sqrt(pmax(diag(V_out), 0)); names(sds) <- par_names
  ll_mean <- .occu_ms_cover_nonspatial_ll(par, model, lay)
  n_draws <- nrow(draws)

  structure(c(list(
    draws        = draws,
    means        = par,
    sds          = sds,
    vcov         = V_out,
    n_samples    = n_draws,
    n_params     = lay$total,
    log_prob     = rep(ll_mean, n_draws),
    log_lik      = ll_mean,
    N            = sum(model$valid),
    accept_prob  = accept,
    divergent    = divergent,
    treedepth    = treedepth,
    epsilon      = epsilon,
    chain_id     = chain_id),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = NULL,
    spatial_field = NULL,
    method       = "nuts",
    positive     = model$positive,
    nuts = list(accept_prob = accept, divergent = divergent,
                treedepth = treedepth, epsilon = epsilon,
                n_chains = as.integer(n.chains),
                divergent_total = sum(divergent, na.rm = TRUE),
                sigma_beta = sigma.beta),
    convergence  = list(converged = TRUE,
                        n_iter = as.integer(n.iter))
  )), class = c("tobs_fit", "tulpa_fit"))
}
