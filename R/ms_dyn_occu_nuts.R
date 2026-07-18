# ms_dyn_occu_nuts.R - NUTS target density for the community / multispecies
# DYNAMIC (multi-season, HMM) occupancy family (ms_dyn_occu()).
#
# The Laplace-EM fit (ms_dyn_occu.R -> .tobs_community_em) profiles the
# per-species first-season occupancy / detection deviations and the two
# INDEPENDENT per-arm community covariances (Sigma_psi1, Sigma_p) out by an
# arrowhead Newton + closed-form covariance M-step, and estimates the SHARED
# community colonisation / extinction (gamma / eps) transition coefficients as
# penalised globals, then reports a Gaussian community-mean posterior. NUTS
# instead samples the EXACT joint posterior -- the community means, the
# per-species deviations {b_s}, the two community covariances, AND the shared
# transition globals -- which removes the Gaussian approximation and the Laplace
# small-cluster attenuation of the variance components.
#
# The target is the dynamic analogue of the community single-season occupancy
# NUTS target (R/ms_occu_nuts.R): the same non-centered per-species blocks
# b_{s,arm} = C_arm z_{s,arm} with the two log-Cholesky community covariances,
# but the data term is the HMM forward marginal (.ms_dyn_occu_fwd_ll_vec) rather
# than the single-season two-state marginal, and there are two extra SHARED
# transition arms (gamma, eps) that carry no per-species random effect. The
# per-arm eta gradients come from one vectorised forward-backward smoothing per
# species (.ms_dyn_occu_fb_vec, shared verbatim with the stMsPGOcc field
# fitter): the Fisher-identity scores
#
#   d ll / d eta_psi1 = w1 - psi1                 (season-1 smoothed posterior)
#   d ll / d eta_p    = sum_t w_t (ndet - nvalid p)
#   d ll / d eta_gam  = col_y - gamma col_n       (colonisation suff-stats)
#   d ll / d eta_eps  = ext_y - eps   ext_n       (extinction  suff-stats)
#
# sandwiched through the site-level designs. The gamma / eps globals accumulate
# their gradient across ALL species (they are shared). The joint log-posterior is
#
#   log p = sum_{s} log L_s(theta)                 # per-species HMM forward
#         - 0.5 ||mu_coef||^2 / sigma.beta^2       # community-mean priors
#         - 0.5 ||global||^2 / sigma.beta^2        # shared-transition priors
#         - 0.5 sum_s ||z_s||^2                    # whitened RE prior (N(0,I))
#         + log p(Sigma log-Cholesky coords)       # weakly-informative hyperpriors
#
# under the NON-CENTERED map b_{s,arm} = C_arm z_{s,arm}. This R version is the
# oracle the C++ FullGradFn port (src/ms_dyn_occu_nuts.cpp) is cross-checked
# against byte-for-byte.


# ---------------------------------------------------------------------------
# Parameter layout
# ---------------------------------------------------------------------------

# Packed NUTS coordinate layout for the community dynamic occupancy model:
#   theta = ( mu [P], {z_s} species-major [S*P], chol_psi1 [q_psi1],
#             chol_p [q_p], global [G] )
# with P = p_psi1 + p_p, mu = (mu_psi1, mu_p), z_s = (z_psi1_s, z_p_s), and
# global = (beta_gamma, beta_eps) the two SHARED transition coefficient blocks.
# `psi1` / `p` are the within-arm coordinate indices (used to slice both mu and
# each z_s); `gam` / `eps` index within `global`.
.tobs_ms_dyn_occu_nuts_layout <- function(p_psi1, p_p, p_gam, p_eps, n_species) {
  P      <- p_psi1 + p_p
  q_psi1 <- as.integer(p_psi1 * (p_psi1 + 1L) / 2L)
  q_p    <- as.integer(p_p    * (p_p    + 1L) / 2L)
  G      <- p_gam + p_eps
  b_off         <- P
  chol_psi1_off <- P + n_species * P
  chol_p_off    <- chol_psi1_off + q_psi1
  global_off    <- chol_p_off + q_p
  total <- global_off + G
  list(
    P = P, p_psi1 = p_psi1, p_p = p_p, p_gam = p_gam, p_eps = p_eps,
    n_species = n_species, q_psi1 = q_psi1, q_p = q_p, G = G,
    psi1 = seq_len(p_psi1),
    p    = p_psi1 + seq_len(p_p),
    mu   = seq_len(P),
    b_off = b_off,
    chol_psi1 = chol_psi1_off + seq_len(q_psi1),
    chol_p    = chol_p_off    + seq_len(q_p),
    global    = global_off + seq_len(G),
    gam = seq_len(p_gam),
    eps = p_gam + seq_len(p_eps),
    global_off = global_off,
    total = total)
}




# ---------------------------------------------------------------------------
# Joint log-posterior + gradient (the NUTS target density / oracle)
# ---------------------------------------------------------------------------

# Full-vector joint log-posterior and its gradient for the community dynamic
# occupancy model. `X_psi1` / `X_p` / `X_gamma` / `X_eps` are the site-level
# designs; `em_stats` a per-species list of precomputed (nvalid, ndet)
# [n_sites x n_seasons] detection sufficient statistics (from
# .ms_dyn_occu_emit_stats). `lay` the layout. Returns list(lp, grad) over the
# packed coordinates. Mirrors the C++ ms_dyn_occu_nuts_eval.
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm}, so the community covariance enters ONLY the data
# term through b = C z. The gamma / eps globals are shared across species and so
# accumulate their gradient over the full species loop.
.tobs_ms_dyn_occu_nuts_logpost <- function(theta, X_psi1, X_p, X_gamma, X_eps,
                                           em_stats, n_sites, n_seasons, lay,
                                           priors, sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu     <- theta[lay$mu]
  global <- theta[lay$global]
  beta_gam <- global[lay$gam]
  beta_eps <- global[lay$eps]
  eta_gam0 <- as.numeric(X_gamma %*% beta_gam)
  eta_eps0 <- as.numeric(X_eps   %*% beta_eps)
  gamma <- stats::plogis(eta_gam0)
  eps   <- stats::plogis(eta_eps0)

  g      <- numeric(lay$total)
  g_mu   <- numeric(P)
  g_gam  <- numeric(n_sites)         # accumulated d ll / d eta_gam over species
  g_eps  <- numeric(n_sites)
  lp     <- 0

  C_psi1 <- .ms_ocs_chol_unpack(theta[lay$chol_psi1], lay$p_psi1)
  C_p    <- .ms_ocs_chol_unpack(theta[lay$chol_p],    lay$p_p)

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_psi1 <- matrix(0, lay$p_psi1, lay$p_psi1)
  A_p    <- matrix(0, lay$p_p,    lay$p_p)

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zpsi1 <- z_s[lay$psi1]; zp <- z_s[lay$p]
    bpsi1 <- mu[lay$psi1] + as.numeric(C_psi1 %*% zpsi1)
    bp    <- mu[lay$p]    + as.numeric(C_p    %*% zp)
    eta_psi1 <- as.numeric(X_psi1 %*% bpsi1)
    eta_p    <- as.numeric(X_p    %*% bp)
    psi1 <- stats::plogis(eta_psi1)
    p    <- stats::plogis(eta_p)

    em <- .ms_dyn_occu_emissions(p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
    if (grad) {
      fb <- .ms_dyn_occu_fb_vec(psi1, p, gamma, eps, em,
                                em_stats[[s]]$nvalid, em_stats[[s]]$ndet,
                                n_sites, n_seasons)
      lp <- lp + fb$ll
      # per-arm eta scores (Fisher identity).
      s_psi1 <- fb$w1 - psi1
      s_p    <- fb$p_score
      # design-sandwiched b gradients (psi1, p arms).
      gpsi1 <- as.numeric(crossprod(X_psi1, s_psi1))
      gp    <- as.numeric(crossprod(X_p,    s_p))
      g_mu[lay$psi1] <- g_mu[lay$psi1] + gpsi1
      g_mu[lay$p]    <- g_mu[lay$p]    + gp
      # z gradient (data part) = C' grad_b.
      g[bidx[lay$psi1]] <- g[bidx[lay$psi1]] + as.numeric(crossprod(C_psi1, gpsi1))
      g[bidx[lay$p]]    <- g[bidx[lay$p]]    + as.numeric(crossprod(C_p,    gp))
      A_psi1 <- A_psi1 + outer(gpsi1, zpsi1)
      A_p    <- A_p    + outer(gp,    zp)
      # shared transition eta scores, accumulated across species.
      g_gam <- g_gam + (fb$col_y - gamma * fb$col_n)
      g_eps <- g_eps + (fb$ext_y - eps   * fb$ext_n)
    } else {
      lp <- lp + .ms_dyn_occu_fwd_ll_vec(psi1, p, gamma, eps, em,
                                         n_sites, n_seasons)
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior ----
  arms <- list(
    list(chol = lay$chol_psi1, A = A_psi1, C = C_psi1, Pa = lay$p_psi1),
    list(chol = lay$chol_p,    A = A_p,    C = C_p,    Pa = lay$p_p))
  for (arm in arms) {
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], arm$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, arm$Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  lp <- lp - 0.5 * ib2 * sum(mu^2)
  g_mu <- g_mu - ib2 * mu
  if (grad) g[lay$mu] <- g[lay$mu] + g_mu

  # ---- shared-transition (global) priors + design-sandwiched gradient ----
  lp <- lp - 0.5 * ib2 * sum(global^2)
  if (grad) {
    g_glob <- numeric(lay$G)
    g_glob[lay$gam] <- as.numeric(crossprod(X_gamma, g_gam)) - ib2 * beta_gam
    g_glob[lay$eps] <- as.numeric(crossprod(X_eps,   g_eps)) - ib2 * beta_eps
    g[lay$global] <- g_glob
  }

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Reconstruct the per-species deviation matrix b (S x P) from a packed coordinate
# vector under the non-centered map b_{s,arm} = C_arm z_{s,arm}.
.tobs_ms_dyn_occu_nuts_b_from_z <- function(theta, lay) {
  C_psi1 <- .ms_ocs_chol_unpack(theta[lay$chol_psi1], lay$p_psi1)
  C_p    <- .ms_ocs_chol_unpack(theta[lay$chol_p],    lay$p_p)
  B <- matrix(0, lay$n_species, lay$P)
  for (s in seq_len(lay$n_species)) {
    z <- theta[.ms_ocs_b_idx(lay, s)]
    B[s, lay$psi1] <- as.numeric(C_psi1 %*% z[lay$psi1])
    B[s, lay$p]    <- as.numeric(C_p    %*% z[lay$p])
  }
  B
}


# ---------------------------------------------------------------------------
# Warm-start pieces (mirrors .tobs_fit_ms_dyn_occu + .tobs_ms_occu_nuts_pieces)
# ---------------------------------------------------------------------------

# Build the per-species emission sufficient statistics, the C++-spec (n_valid /
# n_det) matrix lists, and the warm-start sp_ll / sp_grad closures (the same
# analytic Fisher-identity forward-backward gradient the stMsPGOcc fitter uses),
# so the NUTS warm start runs the same community EM the laplace path does.
.tobs_ms_dyn_occu_nuts_pieces <- function(model) {
  pi_list <- model$process_info
  P_psi1 <- pi_list[[1L]]$p; P_p <- pi_list[[2L]]$p
  P_gam  <- pi_list[[3L]]$p; P_eps <- pi_list[[4L]]$p
  P <- P_psi1 + P_p; G <- P_gam + P_eps
  S <- model$n_species; Ns <- model$n_sites; T <- model$n_seasons

  psi1_idx <- seq_len(P_psi1); p_idx <- P_psi1 + seq_len(P_p)
  gam_idx  <- seq_len(P_gam);  eps_idx <- P_gam + seq_len(P_eps)
  arm_idx  <- list(psi1 = psi1_idx, p = p_idx)

  X_psi1  <- model$X_psi1; X_p <- model$X_p
  X_gamma <- model$X_gamma; X_eps <- model$X_eps
  ys_list <- lapply(seq_len(S), function(s) model$y[, , , s])
  vs_list <- lapply(seq_len(S), function(s) model$valid[, , , s])
  em_stats <- lapply(seq_len(S), function(s)
    .ms_dyn_occu_emit_stats(ys_list[[s]], vs_list[[s]], Ns, T))
  nv_list <- lapply(em_stats, function(e) matrix(as.integer(e$nvalid), Ns, T))
  nd_list <- lapply(em_stats, function(e) matrix(as.integer(e$ndet),   Ns, T))

  eta_of <- function(theta, global) {
    list(psi1 = stats::plogis(as.numeric(X_psi1 %*% theta[psi1_idx])),
         p    = stats::plogis(as.numeric(X_p    %*% theta[p_idx])),
         gam  = stats::plogis(as.numeric(X_gamma %*% global[gam_idx])),
         eps  = stats::plogis(as.numeric(X_eps   %*% global[eps_idx])))
  }
  sp_ll <- function(s, theta, global) {
    e  <- eta_of(theta, global)
    em <- .ms_dyn_occu_emissions(e$p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
    .ms_dyn_occu_fwd_ll_vec(e$psi1, e$p, e$gam, e$eps, em, Ns, T)
  }
  sp_grad <- function(s, theta, global) {
    e  <- eta_of(theta, global)
    em <- .ms_dyn_occu_emissions(e$p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
    fb <- .ms_dyn_occu_fb_vec(e$psi1, e$p, e$gam, e$eps, em,
                              em_stats[[s]]$nvalid, em_stats[[s]]$ndet, Ns, T)
    c(as.numeric(crossprod(X_psi1, fb$w1 - e$psi1)),
      as.numeric(crossprod(X_p,    fb$p_score)),
      as.numeric(crossprod(X_gamma, fb$col_y - e$gam * fb$col_n)),
      as.numeric(crossprod(X_eps,   fb$ext_y - e$eps * fb$ext_n)))
  }

  # ---- warm-start inits (mirror .tobs_fit_ms_dyn_occu) ----
  clamp01 <- function(q) min(max(q, 1e-3), 1 - 1e-3)
  occ_props <- numeric(S); det_rates <- numeric(S)
  for (s in seq_len(S)) {
    v <- vs_list[[s]]; yy <- ys_list[[s]]
    site_det1 <- vapply(seq_len(Ns), function(i)
      any(yy[i, , 1L][v[i, , 1L]] == 1L), logical(1))
    occ_props[s] <- mean(site_det1)
    detected <- yy[v]
    det_rates[s] <- if (length(detected)) mean(detected == 1L) else NA_real_
  }
  init_mu <- numeric(P)
  init_mu[psi1_idx][1L] <- stats::qlogis(clamp01(mean(occ_props)))
  dr <- mean(det_rates[is.finite(det_rates)]); if (!is.finite(dr)) dr <- 0.3
  init_mu[p_idx][1L] <- stats::qlogis(clamp01(dr))
  init_global <- numeric(G)
  init_global[gam_idx][1L] <- stats::qlogis(0.15)
  init_global[eps_idx][1L] <- stats::qlogis(0.10)

  list(P_psi1 = P_psi1, P_p = P_p, P = P, P_gam = P_gam, P_eps = P_eps, G = G,
       S = S, Ns = Ns, T = T, arm_idx = arm_idx, gam_idx = gam_idx,
       eps_idx = eps_idx, X_psi1 = X_psi1, X_p = X_p, X_gamma = X_gamma,
       X_eps = X_eps, em_stats = em_stats, nv_list = nv_list, nd_list = nd_list,
       sp_ll = sp_ll, sp_grad = sp_grad, init_mu = init_mu,
       init_global = init_global)
}

# Pack a community Laplace-EM fit into the full NUTS coordinate vector: community
# means, the two community covariances as log-Cholesky coordinates, the whitened
# per-species deviations z_s = C_arm^{-1} b_s, and the shared transition globals.
.tobs_ms_dyn_occu_nuts_pack_init <- function(em, lay, pieces) {
  theta <- numeric(lay$total)
  theta[lay$mu]     <- as.numeric(em$mu)
  theta[lay$global] <- as.numeric(em$global)
  C_psi1 <- t(chol(.ms_ocs_pd(as.matrix(em$Sigma$psi1))))
  C_p    <- t(chol(.ms_ocs_pd(as.matrix(em$Sigma$p))))
  theta[lay$chol_psi1] <- .ms_ocs_chol_pack(C_psi1)
  theta[lay$chol_p]    <- .ms_ocs_chol_pack(C_p)
  B <- do.call(rbind, em$b_list)                     # S x P
  for (s in seq_len(lay$n_species)) {
    z_s <- numeric(lay$P)
    z_s[lay$psi1] <- forwardsolve(C_psi1, B[s, pieces$arm_idx$psi1])
    z_s[lay$p]    <- forwardsolve(C_p,    B[s, pieces$arm_idx$p])
    theta[.ms_ocs_b_idx(lay, s)] <- z_s
  }
  theta
}



# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the community dynamic occupancy model
# ---------------------------------------------------------------------------

# Sample the exact joint posterior of a non-spatial community dynamic occupancy
# model via tulpa's NUTS engine and the in-tree C++ FullGradFn
# (cpp_ms_dyn_occu_nuts), warm-started at the community Laplace-EM mode with a
# diagonal Laplace metric, then reconstruct the EM-shaped `res` from the draws and
# package it through build_ms_dyn_occu_fit so coef / vcov / confint / ranef /
# fitted / simulate read the NUTS posterior. The full per-draw parameter vector is
# kept under `fit$nuts$draws`.
.tobs_fit_ms_dyn_occu_nuts <- function(model,
                                       priors = NULL,
                                       sigma.beta = 5,
                                       n.iter = 1000L, n.warmup = 1000L,
                                       n.chains = 1L, max.treedepth = 10L,
                                       adapt.delta = 0.9, seed = 1L,
                                       max.iter = 200L, tol = 1e-4,
                                       verbose = FALSE, ...) {
  dots <- list(...)
  newton.max <- as.integer(dots$newton.max %||% 30L)
  pieces <- .tobs_ms_dyn_occu_nuts_pieces(model)
  lay <- .tobs_ms_dyn_occu_nuts_layout(pieces$P_psi1, pieces$P_p,
                                       pieces$P_gam, pieces$P_eps, pieces$S)
  pri <- .ms_ocs_nuts_priors()

  # Warm start at the community Laplace-EM mode (same engine as the laplace path).
  em <- .tobs_community_em(
    S = pieces$S, P = pieces$P, arm_idx = pieces$arm_idx,
    sp_ll = pieces$sp_ll, sp_grad = pieces$sp_grad,
    init_mu = pieces$init_mu, init_global = pieces$init_global,
    penalize_global = TRUE, sigma_beta = sigma.beta, priors = priors,
    sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = newton.max, verbose = FALSE)

  theta0 <- .tobs_ms_dyn_occu_nuts_pack_init(em, lay, pieces)
  spec <- list(X_psi1 = pieces$X_psi1, X_p = pieces$X_p,
               X_gamma = pieces$X_gamma, X_eps = pieces$X_eps,
               n_sites = pieces$Ns, n_seasons = pieces$T, n_species = pieces$S,
               n_valid = pieces$nv_list, n_det = pieces$nd_list)
  inv_metric <- .ms_ocs_fd_metric(
    function(th) cpp_ms_dyn_occu_nuts_joint_logpost(spec, th, pri, sigma.beta)$grad,
    theta0)

  run_chain <- function(ch) {
    cpp_ms_dyn_occu_nuts(
      spec, theta0 = theta0, pri = pri, sigma_beta = sigma.beta,
      inv_metric = inv_metric,
      n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup),
      max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta,
      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }

  rc <- .ms_ocs_run_chains(run_chain, n.chains)
  draws <- rc$draws

  # ---- reconstruct the .tobs_community_em `res` shape from the draws ----
  par     <- colMeans(draws)
  mu_hat  <- par[lay$mu]
  glob_hat <- par[lay$global]
  fe_idx  <- c(lay$mu, lay$global)
  Vf      <- stats::cov(draws[, fe_idx, drop = FALSE])

  Sigma_psi1 <- .ms_ocs_sig_mean(draws, lay$chol_psi1, lay$p_psi1)
  Sigma_p    <- .ms_ocs_sig_mean(draws, lay$chol_p,    lay$p_p)

  B_bar <- matrix(0, pieces$S, lay$P)
  for (i in seq_len(nrow(draws)))
    B_bar <- B_bar + .tobs_ms_dyn_occu_nuts_b_from_z(draws[i, ], lay)
  B_bar <- B_bar / nrow(draws)
  b_list <- lapply(seq_len(pieces$S), function(s) B_bar[s, ])

  # Data-only marginal log-lik at the posterior mean (over reconstructed b_s).
  ll_mean <- 0
  gam <- stats::plogis(as.numeric(pieces$X_gamma %*% glob_hat[pieces$gam_idx]))
  eps <- stats::plogis(as.numeric(pieces$X_eps   %*% glob_hat[pieces$eps_idx]))
  for (s in seq_len(pieces$S)) {
    psi1 <- stats::plogis(as.numeric(pieces$X_psi1 %*% (mu_hat[lay$psi1] + B_bar[s, lay$psi1])))
    p    <- stats::plogis(as.numeric(pieces$X_p    %*% (mu_hat[lay$p]    + B_bar[s, lay$p])))
    em <- .ms_dyn_occu_emissions(p, pieces$em_stats[[s]]$nvalid,
                                 pieces$em_stats[[s]]$ndet)
    ll_mean <- ll_mean + .ms_dyn_occu_fwd_ll_vec(psi1, p, gam, eps, em,
                                                 pieces$Ns, pieces$T)
  }

  res_em <- list(
    mu = unname(mu_hat), global = unname(glob_hat), b_list = b_list,
    Sigma = list(psi1 = Sigma_psi1, p = Sigma_p),
    Vf = Vf, logML = ll_mean, converged = TRUE, n_iter = em$n_iter)

  fit <- build_ms_dyn_occu_fit(model, res_em, pieces$arm_idx,
                               pieces$gam_idx, pieces$eps_idx)
  fit$method <- "nuts"
  fit$log_prob <- rep(ll_mean, nrow(draws))
  .ms_ocs_finalize_nuts_fit(fit, rc, lay, n.chains, sigma_beta = sigma.beta)
}
