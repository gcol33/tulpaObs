# ms_int_occu_nuts.R - NUTS target density (R oracle) for the community /
# multispecies INTEGRATED occupancy family (ms_int_occu()).
#
# The Laplace-EM fit (ms_int_occu.R -> .tobs_community_em) profiles the
# per-species occupancy / per-source detection deviations and the D + 1
# INDEPENDENT per-arm community covariances (Sigma_psi, Sigma_p1, ..., Sigma_pD)
# out by an arrowhead Newton + closed-form covariance M-step, then reports a
# Gaussian community-mean posterior. NUTS instead samples the EXACT joint
# posterior -- the community means, the per-species deviations {b_s}, AND the
# D + 1 community covariances -- which removes the Gaussian approximation and the
# Laplace small-cluster attenuation of the variance components.
#
# The target is the multi-source generalisation of the community single-season
# occupancy NUTS target (R/ms_occu_nuts.R): the SAME non-centered per-species
# blocks b_{s,arm} = C_arm z_{s,arm} with a log-Cholesky community covariance per
# arm, and the data term is the multi-source two-state per-(species, site)
# marginal (.ms_int_occu_sp_ll, D detection sources) rather than the single-source
# case. There are NO shared globals (unlike the dynamic family's gamma / eps): a
# source's detection arm is a plain per-species random effect. The joint
# log-posterior is
#
#   log p = sum_{s,i} log L_{s,i}(theta)               # per-species-site marginal
#         - 0.5 ||mu_coef||^2 / sigma.beta^2           # community-mean priors
#         - 0.5 sum_s ||z_s||^2                        # whitened RE prior (N(0,I))
#         + log p(Sigma log-Cholesky coords)           # weakly-informative hyperpriors
#
# under the NON-CENTERED map b_{s,arm} = C_arm z_{s,arm}. This R version is the
# oracle a C++ FullGradFn port (src/ms_int_occu_nuts.cpp) will be cross-checked
# against byte-for-byte, mirroring the ms_occu / ms_dyn_occu recipe.


# ---------------------------------------------------------------------------
# Parameter layout
# ---------------------------------------------------------------------------

# Packed NUTS coordinate layout for the community integrated occupancy model:
#   theta = ( mu [P], {z_s} species-major [S*P], chol_psi [q_psi],
#             chol_p1 [q_p1], ..., chol_pD [q_pD] )
# with P = p_psi + sum_d p_pd, mu = (mu_psi, mu_p1, ..., mu_pD), z_s stacked the
# same way. `psi` is the within-arm occupancy coordinate slice; `p` a list of D
# within-arm detection slices (each used to slice both mu and each z_s). `P_p` is
# the per-source detection coefficient count.
.tobs_ms_int_occu_nuts_layout <- function(P_psi, P_p, n_species) {
  D <- length(P_p)
  P <- P_psi + sum(P_p)
  psi <- seq_len(P_psi)
  p_slices <- vector("list", D); off <- P_psi
  for (d in seq_len(D)) { p_slices[[d]] <- off + seq_len(P_p[d]); off <- off + P_p[d] }

  b_off <- P
  q_psi <- .ms_ocs_chol_dim(P_psi)
  q_p   <- vapply(P_p, .ms_ocs_chol_dim, integer(1))
  coff  <- P + n_species * P
  chol_psi <- coff + seq_len(q_psi); coff <- coff + q_psi
  chol_p   <- vector("list", D)
  for (d in seq_len(D)) { chol_p[[d]] <- coff + seq_len(q_p[d]); coff <- coff + q_p[d] }

  list(P = P, P_psi = P_psi, P_p = P_p, D = D, n_species = n_species,
       psi = psi, p = p_slices, mu = seq_len(P), b_off = b_off,
       q_psi = q_psi, q_p = q_p, chol_psi = chol_psi, chol_p = chol_p,
       total = coff)
}


# ---------------------------------------------------------------------------
# Joint log-posterior + gradient (the NUTS target density / oracle)
# ---------------------------------------------------------------------------

# Full-vector joint log-posterior and its gradient for the community integrated
# occupancy model. `X_psi` [n_sites x P_psi]; `X_p` a list of D site-level
# detection designs; `summaries` the per-species .ms_int_occu_sp_summary list;
# `lay` the layout. Returns list(lp, grad) over the packed coordinates.
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm}, so each community covariance leaves the b-prior
# (z ~ N(0, I)) and enters ONLY the data term through b = C z.
.tobs_ms_int_occu_nuts_logpost <- function(theta, X_psi, X_p, summaries, lay,
                                           priors, sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species; D <- lay$D
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0

  # Unpack the per-arm Cholesky factors (psi + D detection).
  C_psi <- .ms_ocs_chol_unpack(theta[lay$chol_psi], lay$P_psi)
  C_p   <- lapply(seq_len(D),
                  function(d) .ms_ocs_chol_unpack(theta[lay$chol_p[[d]]], lay$P_p[d]))

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_psi <- matrix(0, lay$P_psi, lay$P_psi)
  A_p   <- lapply(seq_len(D), function(d) matrix(0, lay$P_p[d], lay$P_p[d]))

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zpsi <- z_s[lay$psi]
    bpsi <- mu[lay$psi] + as.numeric(C_psi %*% zpsi)
    eta_psi <- as.numeric(X_psi %*% bpsi)
    zp   <- vector("list", D); eta_p <- vector("list", D)
    for (d in seq_len(D)) {
      zp[[d]]    <- z_s[lay$p[[d]]]
      bpd        <- mu[lay$p[[d]]] + as.numeric(C_p[[d]] %*% zp[[d]])
      eta_p[[d]] <- as.numeric(X_p[[d]] %*% bpd)
    }
    lp <- lp + .ms_int_occu_sp_ll(eta_psi, eta_p, summaries[[s]])
    if (grad) {
      gvec <- .ms_int_occu_sp_grad(eta_psi, eta_p, summaries[[s]], X_psi, X_p)
      gpsi <- gvec[lay$psi]                       # grad_b psi (coefficient space)
      g_mu[lay$psi] <- g_mu[lay$psi] + gpsi
      g[bidx[lay$psi]] <- g[bidx[lay$psi]] + as.numeric(crossprod(C_psi, gpsi))
      A_psi <- A_psi + outer(gpsi, zpsi)
      for (d in seq_len(D)) {
        gpd <- gvec[lay$p[[d]]]
        g_mu[lay$p[[d]]] <- g_mu[lay$p[[d]]] + gpd
        g[bidx[lay$p[[d]]]] <- g[bidx[lay$p[[d]]]] +
          as.numeric(crossprod(C_p[[d]], gpd))
        A_p[[d]] <- A_p[[d]] + outer(gpd, zp[[d]])
      }
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior, per arm ----
  arms <- c(list(list(chol = lay$chol_psi, A = A_psi, C = C_psi, Pa = lay$P_psi)),
            lapply(seq_len(D), function(d)
              list(chol = lay$chol_p[[d]], A = A_p[[d]], C = C_p[[d]], Pa = lay$P_p[d])))
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
  g[lay$mu] <- g[lay$mu] + g_mu

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Reconstruct the per-species deviation matrix b (S x P) from a packed coordinate
# vector under the non-centered map b_{s,arm} = C_arm z_{s,arm}.
.tobs_ms_int_occu_nuts_b_from_z <- function(theta, lay) {
  D <- lay$D
  C_psi <- .ms_ocs_chol_unpack(theta[lay$chol_psi], lay$P_psi)
  C_p   <- lapply(seq_len(D),
                  function(d) .ms_ocs_chol_unpack(theta[lay$chol_p[[d]]], lay$P_p[d]))
  B <- matrix(0, lay$n_species, lay$P)
  for (s in seq_len(lay$n_species)) {
    z <- theta[.ms_ocs_b_idx(lay, s)]
    B[s, lay$psi] <- as.numeric(C_psi %*% z[lay$psi])
    for (d in seq_len(D))
      B[s, lay$p[[d]]] <- as.numeric(C_p[[d]] %*% z[lay$p[[d]]])
  }
  B
}

# Pack a community Laplace-EM fit into the full NUTS coordinate vector: the
# community means, the D + 1 community covariances as log-Cholesky coordinates,
# and the whitened per-species deviations z_s = C_arm^{-1} b_s. `arm_idx` is the
# EM arm layout (list(psi=, p1=, ..., pD=)) mapping each arm to its coefficient
# columns of the EM b matrix.
.tobs_ms_int_occu_nuts_pack_init <- function(em, lay, arm_idx) {
  D <- lay$D
  theta <- numeric(lay$total)
  theta[lay$mu] <- as.numeric(em$mu)

  C_psi <- t(chol(.ms_ocs_pd(as.matrix(em$Sigma$psi))))
  theta[lay$chol_psi] <- .ms_ocs_chol_pack(C_psi)
  C_p <- vector("list", D)
  for (d in seq_len(D)) {
    C_p[[d]] <- t(chol(.ms_ocs_pd(as.matrix(em$Sigma[[names(arm_idx)[d + 1L]]]))))
    theta[lay$chol_p[[d]]] <- .ms_ocs_chol_pack(C_p[[d]])
  }

  B <- do.call(rbind, em$b_list)                 # S x P
  for (s in seq_len(lay$n_species)) {
    z_s <- numeric(lay$P)
    z_s[lay$psi] <- forwardsolve(C_psi, B[s, arm_idx$psi])
    for (d in seq_len(D))
      z_s[lay$p[[d]]] <- forwardsolve(C_p[[d]], B[s, arm_idx[[d + 1L]]])
    theta[.ms_ocs_b_idx(lay, s)] <- z_s
  }
  theta
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the community integrated occupancy model
# ---------------------------------------------------------------------------

# Warm-start the community Laplace-EM (the same engine the laplace path runs),
# pack the mode into the NUTS coordinate vector, and sample the exact joint
# posterior via the in-tree C++ FullGradFn (cpp_ms_int_occu_nuts). The community
# means, per-species deviations, and D + 1 community covariances are sampled
# jointly; the returned fit reuses build_ms_int_occu_fit so the tobs_fit surface
# matches the laplace path. Mirrors .tobs_fit_ms_dyn_occu_nuts (no shared globals).
.tobs_fit_ms_int_occu_nuts <- function(model,
                                       sigma.beta = 5,
                                       n.iter = 1000L, n.warmup = 1000L,
                                       n.chains = 1L, max.treedepth = 10L,
                                       adapt.delta = 0.9, seed = 1L,
                                       max.iter = 100L, tol = 1e-4,
                                       newton.max = 30L, verbose = FALSE,
                                       ...) {
  pi_list <- model$process_info
  D       <- model$n_sources
  P_psi   <- pi_list[[1L]]$p
  P_p     <- vapply(seq_len(D), function(d) pi_list[[d + 1L]]$p, integer(1))
  P       <- P_psi + sum(P_p)
  S       <- model$n_species

  psi_idx <- seq_len(P_psi)
  p_idx   <- vector("list", D); off <- P_psi
  for (d in seq_len(D)) { p_idx[[d]] <- off + seq_len(P_p[d]); off <- off + P_p[d] }
  arm_idx <- c(list(psi = psi_idx),
               stats::setNames(p_idx, model$process_names))

  X_psi <- model$X_psi
  X_p   <- model$X_p
  summaries <- model$summaries

  eta_from_theta <- function(theta) {
    eta_psi <- as.numeric(X_psi %*% theta[psi_idx])
    eta_p   <- lapply(seq_len(D), function(d) as.numeric(X_p[[d]] %*% theta[p_idx[[d]]]))
    list(psi = eta_psi, p = eta_p)
  }
  sp_ll <- function(s, theta, global) {
    e <- eta_from_theta(theta); .ms_int_occu_sp_ll(e$psi, e$p, summaries[[s]])
  }
  sp_grad <- function(s, theta, global) {
    e <- eta_from_theta(theta)
    .ms_int_occu_sp_grad(e$psi, e$p, summaries[[s]], X_psi, X_p)
  }

  clp <- function(x) min(max(x, 1e-3), 1 - 1e-3)
  any_det_prop <- mean(vapply(summaries, function(z) mean(z$any_det), numeric(1)))
  mu0 <- numeric(P)
  mu0[psi_idx][1L] <- stats::qlogis(clp(any_det_prop))
  for (d in seq_len(D)) {
    det_sites <- vapply(summaries, function(z) sum(z$n_det[, d]), numeric(1))
    val_sites <- vapply(summaries, function(z) sum(z$n_valid[z$any_det, d]), numeric(1))
    rate <- if (sum(val_sites) > 0) sum(det_sites) / sum(val_sites) else 0.5
    mu0[p_idx[[d]]][1L] <- stats::qlogis(clp(rate))
  }

  em <- .tobs_community_em(
    S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
    init_mu = mu0, init_global = numeric(0), penalize_global = FALSE,
    sigma_beta = sigma.beta, priors = NULL, sigma_init = 0.3,
    max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = as.integer(newton.max), verbose = FALSE)

  lay    <- .tobs_ms_int_occu_nuts_layout(P_psi, P_p, S)
  pri    <- .ms_ocs_nuts_priors()
  theta0 <- .tobs_ms_int_occu_nuts_pack_init(em, lay, arm_idx)

  # C++ spec: n_valid / n_det as D integer matrices [n_sites x n_species].
  n_sites <- model$n_sites
  nv_list <- lapply(seq_len(D), function(d) {
    m <- vapply(summaries, function(z) as.integer(z$n_valid[, d]), integer(n_sites))
    matrix(as.integer(m), n_sites, S)
  })
  nd_list <- lapply(seq_len(D), function(d) {
    m <- vapply(summaries, function(z) as.integer(z$n_det[, d]), integer(n_sites))
    matrix(as.integer(m), n_sites, S)
  })
  spec <- list(X_psi = X_psi, X_p = X_p, n_valid = nv_list, n_det = nd_list,
               n_sites = n_sites, n_species = S, D = D)

  inv_metric <- .ms_ocs_fd_metric(
    function(th) cpp_ms_int_occu_nuts_joint_logpost(spec, th, pri, sigma.beta)$grad,
    theta0)

  run_chain <- function(ch) {
    cpp_ms_int_occu_nuts(
      spec, theta0 = theta0, pri = pri, sigma_beta = sigma.beta,
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose))
  }
  rc <- .ms_ocs_run_chains(run_chain, n.chains)
  draws     <- rc$draws
  accept    <- rc$accept
  divergent <- rc$divergent
  treedepth <- rc$treedepth
  epsilon   <- rc$epsilon
  rhat_ess  <- rc$rhat_ess

  # ---- reconstruct the .tobs_community_em `fit` shape from the draws ----
  par    <- colMeans(draws)
  mu_hat <- par[lay$mu]
  Vf     <- stats::cov(draws[, lay$mu, drop = FALSE])

  Sigma <- list(psi = .ms_ocs_sig_mean(draws, lay$chol_psi, lay$P_psi))
  for (d in seq_len(D)) {
    Sigma[[names(arm_idx)[d + 1L]]] <-
      .ms_ocs_sig_mean(draws, lay$chol_p[[d]], lay$P_p[d])
  }

  B_bar <- matrix(0, S, lay$P)
  for (i in seq_len(nrow(draws)))
    B_bar <- B_bar + .tobs_ms_int_occu_nuts_b_from_z(draws[i, ], lay)
  B_bar <- B_bar / nrow(draws)
  b_list <- lapply(seq_len(S), function(s) B_bar[s, ])

  # Data-only marginal log-lik at the posterior mean over reconstructed b_s.
  ll_mean <- 0
  for (s in seq_len(S)) {
    beta_s <- mu_hat + B_bar[s, ]
    e <- eta_from_theta(beta_s)
    ll_mean <- ll_mean + .ms_int_occu_sp_ll(e$psi, e$p, summaries[[s]])
  }

  res_em <- list(mu = unname(mu_hat), b_list = b_list, Sigma = Sigma,
                 Vf = Vf, logML = ll_mean, converged = TRUE, n_iter = em$n_iter)

  fit <- build_ms_int_occu_fit(model, res_em, arm_idx)
  fit$method   <- "nuts"
  fit$log_prob <- rep(ll_mean, nrow(draws))
  fit$nuts <- list(
    draws = draws, layout = lay, accept_prob = accept, divergent = divergent,
    treedepth = treedepth, epsilon = epsilon, n_chains = as.integer(n.chains),
    divergent_total = sum(divergent), sigma_beta = sigma.beta)
  if (!is.null(rhat_ess)) {
    fit$nuts$rhat     <- rhat_ess$rhat
    fit$nuts$ess      <- rhat_ess$ess
    fit$nuts$max_rhat <- max(rhat_ess$rhat, na.rm = TRUE)
    fit$nuts$min_ess  <- min(rhat_ess$ess,  na.rm = TRUE)
  }
  fit
}
