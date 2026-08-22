# ms_occu_nuts.R - NUTS target density for the community / multispecies
# single-season occupancy (ms_occu()).
#
# The Laplace-EM fit (ms_occu.R -> .tobs_community_em) profiles the per-species
# coefficient deviations and the two INDEPENDENT per-arm community covariances
# (Sigma_psi, Sigma_p) out by an arrowhead Newton + closed-form covariance
# M-step, then reports a Gaussian community-mean posterior c(mu_psi, mu_p). NUTS
# instead samples the EXACT joint posterior -- the community means, the
# per-species deviations {b_s}, AND the community covariances -- which removes
# the Gaussian approximation, gives calibrated (non-Gaussian) community
# intervals, and removes the Laplace small-cluster attenuation of the variance
# components (means are unbiased under both; the per-arm SD is what the sampler
# de-attenuates).
#
# The target is the occupancy analogue of the community N-mixture NUTS target
# (R/ms_abun_nuts.R): the data term is the occupancy two-state per-site marginal
# (the single-source case of .ms_int_occu_sp_ll, site-level detection), summed
# over (species, site) with the per-species coefficient eta_s = X . (mu + b_s),
# and a hierarchical Gaussian community prior b_{s,arm} ~ N(0, Sigma_arm). The
# two community covariances are INDEPENDENT (psi and p each carry their own),
# matching .tobs_community_em and spOccupancy msPGOcc; this independence is what
# a correct sampler needs (a single shared block would tie the two arms). The
# covariances are carried by their log-Cholesky factors (the same parametrisation
# and helpers .ms_ocs_chol_* as the community N-mixture / spatial-factor
# targets). The joint log-posterior is
#
#   log p = sum_{s,i} log L_{s,i}(theta)               # per-species-site marginal
#         - 0.5 ||mu_coef||^2 / sigma.beta^2           # community-mean priors
#         - 0.5 sum_s ||z_s||^2                        # whitened RE prior (N(0,I))
#         + log p(Sigma log-Cholesky coords)           # weakly-informative hyperpriors
#
# under the NON-CENTERED map b_{s,arm} = C_arm z_{s,arm}; the data + chol
# coordinate gradients mirror the C++ FullGradFn (src/ms_occu_nuts.cpp)
# byte-for-byte, this R version being the oracle the C++ port is cross-checked
# against.


# ---------------------------------------------------------------------------
# Parameter layout
# ---------------------------------------------------------------------------

# Packed NUTS coordinate layout for the community occupancy model:
#   theta = ( mu [P], {b_s} species-major [S*P], chol_psi [q_psi], chol_p [q_p] )
# with P = p_psi + p_p, mu = (mu_psi, mu_p), b_s = (b_psi_s, b_p_s). `psi` / `p`
# are the within-arm coordinate indices (used to slice both mu and each b_s).
.tobs_ms_occu_nuts_layout <- function(p_psi, p_p, n_species) {
  lay <- .ms_ocs_layout(list(list(name = "psi", width = p_psi),
                             list(name = "p",   width = p_p)), n_species)
  c(.ms_ocs_layout_base(lay),
    list(p_psi = p_psi, p_p = p_p,
         q_psi = lay$q[["psi"]], q_p = lay$q[["p"]],
         psi = lay$idx$psi, p = lay$idx$p,
         chol_psi = lay$chol$psi, chol_p = lay$chol$p))
}


# ---------------------------------------------------------------------------
# Hyperprior specification
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Joint log-posterior + gradient (the NUTS target density / oracle)
# ---------------------------------------------------------------------------

# Full-vector joint log-posterior and its gradient for the community occupancy
# model. `X_psi` / `X_p` are the site-level designs; `summaries` the per-species
# .ms_int_occu_sp_summary list (single-source). `lay` the layout. Returns
# list(lp, grad) over the packed coordinates. Mirrors the C++ ms_occu_nuts_eval
# (src/ms_occu_nuts.cpp).
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm} (C_arm the log-Cholesky factor of Sigma_arm), so the
# community covariance leaves the b-prior (z ~ N(0, I)) and enters ONLY the data
# term through b = C z.
.tobs_ms_occu_nuts_logpost <- function(theta, X_psi, X_p, summaries, lay, priors,
                                       sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0

  C_psi <- .ms_ocs_chol_unpack(theta[lay$chol_psi], lay$p_psi)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$p_p)

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_psi <- matrix(0, lay$p_psi, lay$p_psi)
  A_p   <- matrix(0, lay$p_p,   lay$p_p)

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zpsi <- z_s[lay$psi]; zp <- z_s[lay$p]
    bpsi <- mu[lay$psi] + as.numeric(C_psi %*% zpsi)
    bp   <- mu[lay$p]   + as.numeric(C_p   %*% zp)
    eta_psi <- as.numeric(X_psi %*% bpsi)
    eta_p   <- as.numeric(X_p   %*% bp)
    lp <- lp + .ms_int_occu_sp_ll(eta_psi, list(eta_p), summaries[[s]])
    if (grad) {
      gvec <- .ms_int_occu_sp_grad(eta_psi, list(eta_p), summaries[[s]],
                                   X_psi, list(X_p))
      gpsi <- gvec[lay$psi]              # grad_b psi
      gp   <- gvec[lay$p]               # grad_b p
      g_mu[lay$psi] <- g_mu[lay$psi] + gpsi
      g_mu[lay$p]   <- g_mu[lay$p]   + gp
      # z gradient (data part) = C' grad_b.
      g[bidx[lay$psi]] <- g[bidx[lay$psi]] + as.numeric(crossprod(C_psi, gpsi))
      g[bidx[lay$p]]   <- g[bidx[lay$p]]   + as.numeric(crossprod(C_p,   gp))
      A_psi <- A_psi + outer(gpsi, zpsi)
      A_p   <- A_p   + outer(gp,   zp)
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior ----
  arms <- list(list(chol = lay$chol_psi, A = A_psi, C = C_psi, Pa = lay$p_psi),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$p_p))
  for (arm in arms) {
    Pa <- arm$Pa
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  coef_idx <- c(lay$psi, lay$p)
  lp <- lp - 0.5 * ib2 * sum(mu[coef_idx]^2)
  g_mu[coef_idx] <- g_mu[coef_idx] - ib2 * mu[coef_idx]
  g[lay$mu] <- g[lay$mu] + g_mu

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}


# ---------------------------------------------------------------------------
# Warm-start packing + inverse-mass metric
# ---------------------------------------------------------------------------


# Pack a community Laplace-EM fit (.tobs_community_em output) into the full NUTS
# coordinate vector: the community means, the two community covariances as
# log-Cholesky coordinates, and -- under the non-centered map -- the whitened
# per-species deviations z_s = C_arm^{-1} b_s (so reconstructing b = C z returns
# the warm BLUPs exactly). This is the NUTS initial position and the FD-gradient
# reference point.
.tobs_ms_occu_nuts_pack_init <- function(em, lay, arm_idx) {
  theta <- numeric(lay$total)
  theta[lay$mu] <- as.numeric(em$mu)

  C_psi <- t(chol(.ms_ocs_pd(as.matrix(em$Sigma$psi))))
  C_p   <- t(chol(.ms_ocs_pd(as.matrix(em$Sigma$p))))
  theta[lay$chol_psi] <- .ms_ocs_chol_pack(C_psi)
  theta[lay$chol_p]   <- .ms_ocs_chol_pack(C_p)

  B <- do.call(rbind, em$b_list)              # S x P
  for (s in seq_len(lay$n_species)) {
    z_s <- numeric(lay$P)
    z_s[lay$psi] <- forwardsolve(C_psi, B[s, arm_idx$psi])
    z_s[lay$p]   <- forwardsolve(C_p,   B[s, arm_idx$p])
    theta[.ms_ocs_b_idx(lay, s)] <- z_s
  }
  theta
}

# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the community occupancy model
# ---------------------------------------------------------------------------

# Build the per-species single-source detection summaries + the warm-start sp_ll
# / sp_grad closures the community EM uses, so the NUTS warm start runs the same
# EM the laplace path does (single source of truth).
.tobs_ms_occu_nuts_pieces <- function(model) {
  pi_list <- model$process_info
  P_psi <- pi_list[[1L]]$p
  P_p   <- pi_list[[2L]]$p
  P     <- P_psi + P_p
  S     <- model$n_species
  occ_idx <- seq_len(P_psi)
  p_idx   <- P_psi + seq_len(P_p)
  arm_idx <- list(psi = occ_idx, p = p_idx)
  X_psi <- model$X_occ
  X_p   <- model$X_det
  summaries <- model$summaries

  sp_ll <- function(s, theta, global) {
    eta_psi <- as.numeric(X_psi %*% theta[occ_idx])
    eta_p   <- as.numeric(X_p   %*% theta[p_idx])
    .ms_int_occu_sp_ll(eta_psi, list(eta_p), summaries[[s]])
  }
  sp_grad <- function(s, theta, global) {
    eta_psi <- as.numeric(X_psi %*% theta[occ_idx])
    eta_p   <- as.numeric(X_p   %*% theta[p_idx])
    .ms_int_occu_sp_grad(eta_psi, list(eta_p), summaries[[s]], X_psi, list(X_p))
  }

  clp <- function(q) min(max(q, 1e-3), 1 - 1e-3)
  any_det_prop <- mean(vapply(summaries, function(z) mean(z$any_det), numeric(1)))
  det_n <- sum(vapply(summaries, function(z) sum(z$n_det[, 1L]), numeric(1)))
  val_n <- sum(vapply(summaries, function(z) sum(z$n_valid[z$any_det, 1L]),
                      numeric(1)))
  rate  <- if (val_n > 0) det_n / val_n else 0.5
  mu <- numeric(P)
  mu[occ_idx][1L] <- stats::qlogis(clp(any_det_prop))
  mu[p_idx][1L]   <- stats::qlogis(clp(rate))

  list(P_psi = P_psi, P_p = P_p, P = P, S = S, arm_idx = arm_idx,
       X_psi = X_psi, X_p = X_p, summaries = summaries,
       sp_ll = sp_ll, sp_grad = sp_grad, mu0 = mu)
}

# Sample the exact joint posterior of a non-spatial community single-season
# occupancy model via tulpa's NUTS engine and the in-tree C++ FullGradFn
# (cpp_ms_occu_nuts), warm-started at the community Laplace-EM mode with a
# diagonal Laplace metric, then reconstruct the EM-shaped `fit` from the draws
# and package it through build_ms_occu_fit so coef / vcov / confint / ranef /
# fitted / simulate / richness read the NUTS posterior. The full per-draw
# parameter vector is kept under `fit$nuts$draws` (with the layout).
.tobs_fit_ms_occu_nuts <- function(model,
                                   sigma.beta = NULL,
                                   n.iter = NULL, n.warmup = NULL,
                                   n.chains = NULL, n.thin = NULL,
                                   n.threads = NULL,
                                   n.threads.grad = NULL, max.treedepth = NULL,
                                   adapt.delta = NULL, seed = NULL,
                                   max.iter = 100L, tol = 1e-4,
                                   newton.max = 30L, verbose = FALSE,
                                   ...) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "nuts")

  pieces <- .tobs_ms_occu_nuts_pieces(model)
  lay <- .tobs_ms_occu_nuts_layout(pieces$P_psi, pieces$P_p, pieces$S)
  pri <- .ms_ocs_nuts_priors()

  # Warm start at the community Laplace-EM mode (the same engine the laplace path
  # uses; verbose forced off here).
  em <- .tobs_community_em(
    S = pieces$S, P = pieces$P, arm_idx = pieces$arm_idx,
    sp_ll = pieces$sp_ll, sp_grad = pieces$sp_grad,
    init_mu = pieces$mu0, init_global = numeric(0),
    penalize_global = FALSE, sigma_beta = sigma.beta, priors = NULL,
    sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = as.integer(newton.max), verbose = FALSE)

  theta0 <- .tobs_ms_occu_nuts_pack_init(em, lay, pieces$arm_idx)
  mats   <- .ms_occu_spatial_count_mats(pieces$summaries, model$n_sites,
                                          pieces$S)
  spec <- list(X_psi = pieces$X_psi, X_p = pieces$X_p,
               n_sites = model$n_sites, n_species = pieces$S,
               n_valid = mats$n_valid, n_det = mats$n_det,
               n_threads = as.integer(n.threads.grad))
  inv_metric <- .ms_ocs_fd_metric(
    function(th) cpp_ms_occu_nuts_joint_logpost(spec, th, pri, sigma.beta)$grad,
    theta0)

  run_chain <- function(ch) {
    cpp_ms_occu_nuts(
      spec, theta0 = theta0, pri = pri, sigma_beta = sigma.beta,
      inv_metric = inv_metric,
      n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup),
      max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta,
      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }

  rc <- .ms_ocs_run_chains(run_chain, n.chains, n.thin = n.thin,
                           n.threads = n.threads)
  draws <- rc$draws

  # ---- reconstruct the .tobs_community_em `fit` shape from the draws ----
  par    <- colMeans(draws)
  mu_hat <- par[lay$mu]
  vcov_mu <- stats::cov(draws[, lay$mu, drop = FALSE])

  Sigma_psi <- .ms_ocs_sig_mean(draws, lay$chol_psi, lay$p_psi)
  Sigma_p   <- .ms_ocs_sig_mean(draws, lay$chol_p,   lay$p_p)

  # Per-species BLUPs = posterior mean of the reconstructed deviation b = C z.
  B_bar <- matrix(0, pieces$S, lay$P)
  for (i in seq_len(nrow(draws)))
    B_bar <- B_bar + .ms_ocs_b_from_z(draws[i, ], lay)
  B_bar <- B_bar / nrow(draws)
  b_list <- lapply(seq_len(pieces$S), function(s) B_bar[s, ])

  # Data-only log-lik at the posterior-mean coefficients (over the reconstructed
  # per-species deviations), reusing the per-species marginal.
  ll_mean <- 0
  for (s in seq_len(pieces$S)) {
    eta_psi <- as.numeric(pieces$X_psi %*% (mu_hat[lay$psi] + B_bar[s, lay$psi]))
    eta_p   <- as.numeric(pieces$X_p   %*% (mu_hat[lay$p]   + B_bar[s, lay$p]))
    ll_mean <- ll_mean + .ms_int_occu_sp_ll(eta_psi, list(eta_p),
                                            pieces$summaries[[s]])
  }

  fit_em <- list(
    mu = unname(mu_hat), global = numeric(0), b_list = b_list,
    Sigma = list(psi = Sigma_psi, p = Sigma_p),
    Vf = vcov_mu, logML = ll_mean,
    converged = TRUE, n_iter = em$n_iter)

  fit <- build_ms_occu_fit(model, fit_em, pieces$arm_idx)
  fit$method <- "nuts"
  fit$log_prob <- rep(ll_mean, nrow(draws))

  .ms_ocs_finalize_nuts_fit(fit, rc, lay, n.chains, sigma_beta = sigma.beta)
}
