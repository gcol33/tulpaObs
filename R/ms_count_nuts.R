# =============================================================================
# ms_count_nuts.R - community relative-abundance (count) NUTS (msAbund NUTS).
# Samples the exact joint posterior of the non-spatial community GLMM (community
# means mu, per-species coefficient deviations b_s, community covariance Sigma)
# via the in-tree C++ FullGradFn (src/ms_count_nuts.cpp), warm-started at the
# community Laplace-EM mode. NON-CENTERED: b_{s,arm} = C_arm z_{s,arm},
# z_s ~ N(0, I), so Sigma (log-Cholesky C) enters only the data term.
#
# Three response families share one target, mirroring the ms_count Laplace-EM:
#   Poisson   - one beta arm.
#   NegBin    - beta arm + a per-species dispersion RE log_r_s ~
#               N(mu_log_r, sigma_log_r^2), the second community arm (the same
#               structure as ms_abun / the ms_count negbin Laplace-EM).
#   Gaussian  - beta arm + S FREE per-species residual variances log_phi_s (no
#               community prior; matches the Laplace outer loop that estimates
#               each phi_s independently), each with a weakly-informative prior.
# The reduced counterpart of ms_abun_nuts.R (no detection, no latent N).
# =============================================================================


# Packed layout: mu (P), z species-major (S*P), chol_beta (q_beta),
#   [chol_logr (1) if negbin], [log_phi (S) if gaussian]. P = p_beta (+1 negbin).
.tobs_ms_count_nuts_layout <- function(p_beta, n_species, family = "poisson") {
  is_nb    <- identical(family, "negbin")
  is_gauss <- identical(family, "gaussian")
  arms <- list(list(name = "beta", width = p_beta))
  # The per-species log-size is a one-dimensional community arm.
  if (is_nb) arms <- c(arms, list(list(name = "logr", width = 1L)))
  # The gaussian per-species log-variances are FREE (no community prior), so
  # they trail the chol section as a per-species vector rather than an arm.
  lay <- .ms_ocs_layout(
    arms, n_species,
    trailing = list(list(name = "logphi",
                         size = if (is_gauss) n_species else 0L)))
  c(.ms_ocs_layout_base(lay),
    list(family = family, is_nb = is_nb, is_gauss = is_gauss,
         p_beta = p_beta, q_beta = lay$q[["beta"]],
         beta = lay$idx$beta,
         logr = if (is_nb) lay$idx$logr else integer(0),
         chol_beta = lay$chol$beta,
         chol_logr = if (is_nb) lay$chol$logr else integer(0),
         logphi    = lay$trailing$logphi))
}

# Log-Cholesky hyperprior + weakly-informative gaussian log_phi prior scalars.
.tobs_ms_count_nuts_priors <- function() {
  list(chol_logdiag_mean = log(0.5), chol_logdiag_sd = 1.5, chol_offdiag_sd = 1.0,
       logphi_mean = 0, logphi_sd = 2)
}


# Pack a community Laplace-EM warm start into theta: the means, each arm's
# covariance as log-Cholesky, the whitened per-species deviations z = C^{-1} b,
# and (gaussian) the per-species log residual variances.
.tobs_ms_count_nuts_pack_init <- function(lay, mu_beta, Sigma_beta, B_beta,
                                          mu_logr = NA_real_, sigma_logr = NA_real_,
                                          b_logr = NULL, logphi = NULL) {
  theta <- numeric(lay$total)
  P <- lay$P; S <- lay$n_species
  mu <- as.numeric(mu_beta)
  if (lay$is_nb) mu <- c(mu, mu_logr)
  theta[lay$mu] <- mu

  Sig <- (Sigma_beta + t(Sigma_beta)) / 2
  if (inherits(try(chol(Sig), silent = TRUE), "try-error"))
    Sig <- Sig + diag(max(1e-6, 1e-6 * mean(diag(Sig))), ncol(Sig))
  C_beta <- t(chol(Sig))
  theta[lay$chol_beta] <- .ms_ocs_chol_pack(C_beta)
  if (lay$is_nb)    theta[lay$chol_logr] <- log(sigma_logr)   # 1x1 log-chol = log(sd)
  if (lay$is_gauss) theta[lay$logphi]    <- as.numeric(logphi)

  for (s in seq_len(S)) {
    z_s <- numeric(P)
    z_s[lay$beta] <- forwardsolve(C_beta, B_beta[s, ])   # C z = b, C lower-tri
    if (lay$is_nb) z_s[lay$logr] <- b_logr[s] / sigma_logr
    theta[.ms_ocs_b_idx(lay, s)] <- z_s
  }
  theta
}


# Fit the community count model by NUTS. Warm-starts at the community Laplace-EM
# mode, then samples the exact joint posterior. Poisson / negbin / gaussian.
.tobs_fit_ms_count_nuts <- function(model, sigma.beta = NULL, sigma.logr = NULL,
                                    n.iter = NULL, n.warmup = NULL,
                                    n.chains = NULL, n.thin = NULL,
                                    n.threads = NULL,
                                    n.threads.grad = NULL, max.treedepth = NULL,
                                    adapt.delta = NULL, seed = NULL,
                                    verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "nuts", family = "ms_count")

  response <- model$response %||% "poisson"
  is_nb    <- identical(response, "negbin")
  is_gauss <- identical(response, "gaussian")

  # Warm start: the community Laplace-EM (raw output, so the negbin log_r arm is
  # available to pack).
  em     <- .tobs_ms_count_run_em(model, verbose = FALSE)
  fit    <- em$fit
  P_beta <- em$P_beta; S <- model$n_species
  Y <- matrix(as.numeric(model$y), model$n_sites, S)
  if (!is_gauss) Y <- round(Y)
  X <- model$X

  lay    <- .tobs_ms_count_nuts_layout(P_beta, S, response)
  priors <- .tobs_ms_count_nuts_priors()
  if (is_gauss) {
    # Weakly-informative log_phi prior, centred at the pooled log-variance so it
    # is on the data's scale but wide enough to be swamped by the residuals.
    yv <- model$y[model$valid]
    priors$logphi_mean <- log(max(stats::var(yv), 1e-4))
    priors$logphi_sd   <- 2
  }

  # Warm-start pieces from the raw EM fit.
  mu_beta    <- fit$mu[seq_len(P_beta)]
  Sigma_beta <- as.matrix(fit$Sigma$mu)
  B_beta     <- do.call(rbind, lapply(fit$b_list, function(b) b[seq_len(P_beta)]))
  mu_logr    <- if (is_nb) fit$mu[P_beta + 1L] else NA_real_
  sigma_logr_w <- if (is_nb) max(sqrt(fit$Sigma$disp[1, 1]), 1e-3) else NA_real_
  b_logr     <- if (is_nb) vapply(fit$b_list, function(b) b[P_beta + 1L],
                                  numeric(1)) else NULL
  logphi_w   <- if (is_gauss) log(pmax(em$disp$variance, 1e-8)) else NULL

  theta0 <- .tobs_ms_count_nuts_pack_init(
    lay, mu_beta = mu_beta, Sigma_beta = Sigma_beta, B_beta = B_beta,
    mu_logr = mu_logr, sigma_logr = sigma_logr_w, b_logr = b_logr,
    logphi = logphi_w)

  spec <- list(X = X, y = Y, family = response,
               n_threads = as.integer(n.threads.grad))
  if (is_gauss) {
    spec$logphi_mean <- priors$logphi_mean
    spec$logphi_sd   <- priors$logphi_sd
  }

  run_chains <- function(n_chains, n_threads) {
    # tulpa's NUTS engine takes n_iter as the TOTAL (warmup + sampling) count and
    # returns n_iter - n_warmup post-warmup draws.
    cpp_ms_count_nuts(spec, theta0, priors, sigma.beta, sigma.logr, NULL,
                      as.integer(n.iter + n.warmup), as.integer(n.warmup),
                      as.integer(max.treedepth), adapt.delta,
                      as.integer(seed), isTRUE(verbose), n_chains, n_threads)
  }
  rc <- .ms_ocs_run_chains(run_chains, n.chains, n.thin = n.thin,
                           n.threads = n.threads)
  chains    <- rc$chains
  draws_all <- rc$draws                             # (n_chains*n_sample) x total
  accept    <- rc$accept
  divergent <- rc$divergent
  treedepth <- rc$treedepth
  epsilon   <- rc$epsilon

  # Community means + covariance from the beta-coefficient mu draws.
  mu_draws <- draws_all[, lay$beta, drop = FALSE]
  cn  <- model$process_info[[1L]]$coef_names
  nms <- .tobs_community_mean_names(cn)
  means <- colMeans(mu_draws); names(means) <- nms
  V <- stats::cov(mu_draws); dimnames(V) <- list(nms, nms)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nms
  colnames(mu_draws) <- nms

  # Posterior-mean per-species deviations, community covariance, dispersion.
  nd <- nrow(draws_all)
  mom <- .tobs_cem_moments_new(S, P_beta); Sig_acc <- matrix(0, P_beta, P_beta)
  logr_acc <- numeric(S); sigma_logr_acc <- 0; phi_acc <- numeric(S)
  for (r in seq_len(nd)) {
    th <- draws_all[r, ]
    Bd <- .ms_ocs_b_from_z(th, lay)
    mom <- .tobs_cem_moments_add(mom, Bd[, lay$beta, drop = FALSE])
    C <- .ms_ocs_chol_unpack(th[lay$chol_beta], P_beta); Sig_acc <- Sig_acc + tcrossprod(C)
    if (is_nb) {
      logr_acc <- logr_acc + Bd[, lay$logr]
      sigma_logr_acc <- sigma_logr_acc + exp(th[lay$chol_logr])
    }
    if (is_gauss) phi_acc <- phi_acc + exp(th[lay$logphi])
  }
  blup <- .tobs_cem_moments_mean(mom); coef <- sweep(blup, 2L, means, "+")
  rownames(blup) <- rownames(coef) <- model$species_names
  colnames(blup) <- colnames(coef) <- cn
  Sigma_mu <- Sig_acc / nd; dimnames(Sigma_mu) <- list(cn, cn)

  disp <- NULL
  if (is_nb) {
    mu_log_r_hat <- mean(draws_all[, lay$logr])
    disp <- list(response = "negbin", mu_log_r = mu_log_r_hat,
                 sigma_log_r = sigma_logr_acc / nd,
                 r_s = exp(mu_log_r_hat + logr_acc / nd))
  } else if (is_gauss) {
    disp <- list(response = "gaussian", variance = phi_acc / nd)
  }

  fit <- structure(c(list(
    draws = mu_draws, means = means, sds = sds, vcov = V,
    n_samples = nd, n_params = length(means),
    log_prob = rep(NA_real_, nd), N = sum(model$valid),
    accept_prob = accept, divergent = divergent, treedepth = treedepth,
    epsilon = epsilon),
    list(
    col_names = nms, param_names = nms, n_fixed = length(means),
    fixed_names = nms, process_info = model$process_info,
    model = model, spatial = NULL, method = "nuts",
    ms_community = list(Sigma_mu = Sigma_mu,
                        sd_mu = sqrt(pmax(diag(Sigma_mu), 0)),
                        coef_mu = coef, blup_mu = blup,
                        Cinv = .tobs_cem_moments_cov(mom),
                        blup_idx = list(blup_mu = seq_len(P_beta))),
    ms_dispersion = disp,
    convergence = list(converged = NA, n_iter = as.integer(n.iter))
  )), class = c("tobs_fit", "tulpa_fit"))
  # The reported coefficients are the community beta means; the sampler also
  # carries the per-species z blocks, the packed community Cholesky and (negbin)
  # the community mean log_r, which no summary row names.
  fit <- .tobs_nuts_attach_convergence(
    fit, chains, par_names = nms, cols = lay$beta,
    n_iter = as.integer(n.iter),
    sampler_control = .tobs_sampler_control_snapshot(environment()))
  # Top-level aliases for the per-parameter record.
  fit$rhat <- fit$convergence$rhat
  fit$ess  <- fit$convergence$ess_bulk
  fit
}
