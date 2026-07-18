# ms_occu_cover_nuts.R - NUTS target density (R oracle) for the community /
# multispecies joint occupancy-detection + cover family (ms_occu_cover(), the
# NON-spatial per-arm-covariance community model; #115 part B7).
#
# The Laplace-EM fit (ms_occu_cover.R -> .tobs_fit_ms_occu_cover) profiles the
# per-species occupancy / detection / cover deviations and the THREE independent
# per-arm community covariances (Sigma_occ, Sigma_p, Sigma_pos) out by an
# arrowhead joint Newton + closed-form covariance M-step, sharing one community
# log-dispersion for the positive-cover arm, then reports a Gaussian community-mean
# posterior. Above the AGHQ dimension cap the community VARIANCE carries the
# documented Laplace small-cluster attenuation. NUTS instead samples the EXACT
# joint posterior -- the community means, the per-species deviations {b_s}, the
# three community covariances, AND the shared log-dispersion -- which removes the
# attenuation.
#
# The target is the joint-cover analogue of the community occupancy NUTS targets
# (R/ms_occu_nuts.R, R/ms_int_occu_nuts.R): the SAME non-centered per-species
# blocks b_{s,arm} = C_arm z_{s,arm} with a log-Cholesky community covariance per
# arm (occ + p + pos), plus ONE shared community log-dispersion scalar (the
# positive-arm beta precision / lognormal-or-gaussian residual SD, on the log
# scale) that carries no per-species random effect -- the analogue of the dynamic
# family's shared gamma/eps globals, here a single coordinate. The data term is
# the per-(species, cell) two-state occu_cover marginal (.occu_cover_sp_ll, reused
# verbatim from the laplace path); its coefficient + log-dispersion gradient come
# from .occu_cover_sp_grad. The joint log-posterior is
#
#   log p = sum_{s} log L_s(theta)                     # per-species occu_cover marginal
#         - 0.5 ||mu_coef||^2 / sigma.beta^2           # community-mean priors
#         - 0.5 sum_s ||z_s||^2                        # whitened RE prior (N(0,I))
#         + sum_arm log p(Sigma_arm coords)            # log-Cholesky hyperpriors
#         + log p(log_disp)                            # weakly-informative dispersion prior
#
# under the NON-CENTERED map b_{s,arm} = C_arm z_{s,arm}. This R version is the
# oracle a C++ FullGradFn port (src/ms_occu_cover_nuts.cpp) will be cross-checked
# against, mirroring the ms_occu / ms_dyn_occu / ms_int_occu recipe.


# ---------------------------------------------------------------------------
# Parameter layout
# ---------------------------------------------------------------------------

# Packed NUTS coordinate layout for the community joint occu+cover model:
#   theta = ( mu [P], {z_s} species-major [S*P], chol_occ [q_occ], chol_p [q_p],
#             chol_pos [q_pos], log_disp [1] )
# with P = p_occ + p_p + p_pos, mu = (mu_occ, mu_p, mu_pos), z_s stacked the same
# way, and log_disp the single SHARED community log-dispersion coordinate.
.tobs_ms_occu_cover_nuts_layout <- function(P_occ, P_p, P_pos, n_species) {
  P <- P_occ + P_p + P_pos
  occ <- seq_len(P_occ)
  p   <- P_occ + seq_len(P_p)
  pos <- P_occ + P_p + seq_len(P_pos)

  b_off <- P
  q_occ <- .ms_ocs_chol_dim(P_occ)
  q_p   <- .ms_ocs_chol_dim(P_p)
  q_pos <- .ms_ocs_chol_dim(P_pos)
  coff  <- P + n_species * P
  chol_occ <- coff + seq_len(q_occ); coff <- coff + q_occ
  chol_p   <- coff + seq_len(q_p);   coff <- coff + q_p
  chol_pos <- coff + seq_len(q_pos); coff <- coff + q_pos
  log_disp <- coff + 1L; coff <- coff + 1L

  list(P = P, P_occ = P_occ, P_p = P_p, P_pos = P_pos, n_species = n_species,
       occ = occ, p = p, pos = pos, mu = seq_len(P), b_off = b_off,
       q_occ = q_occ, q_p = q_p, q_pos = q_pos,
       chol_occ = chol_occ, chol_p = chol_p, chol_pos = chol_pos,
       log_disp = log_disp, total = coff)
}


# ---------------------------------------------------------------------------
# Joint log-posterior + gradient (the NUTS target density / oracle)
# ---------------------------------------------------------------------------

# Full-vector joint log-posterior and its gradient. `views` is the per-species
# .ms_occu_cover_species_view list (the single-species occu_cover model views the
# laplace path builds); `lay` the layout. Returns list(lp, grad) over the packed
# coordinates.
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm}, so each community covariance leaves the b-prior and
# enters ONLY the data term. The shared log_disp carries a weakly-informative
# Normal prior on the log scale and accumulates its gradient over all species.
.tobs_ms_occu_cover_nuts_logpost <- function(theta, views, lay, priors,
                                             sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu <- theta[lay$mu]
  log_disp <- theta[lay$log_disp]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  g_ld <- 0
  lp   <- 0

  C_occ <- .ms_ocs_chol_unpack(theta[lay$chol_occ], lay$P_occ)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$P_p)
  C_pos <- .ms_ocs_chol_unpack(theta[lay$chol_pos], lay$P_pos)

  A_occ <- matrix(0, lay$P_occ, lay$P_occ)
  A_p   <- matrix(0, lay$P_p,   lay$P_p)
  A_pos <- matrix(0, lay$P_pos, lay$P_pos)

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zocc <- z_s[lay$occ]; zp <- z_s[lay$p]; zpos <- z_s[lay$pos]
    bocc <- mu[lay$occ] + as.numeric(C_occ %*% zocc)
    bp   <- mu[lay$p]   + as.numeric(C_p   %*% zp)
    bpos <- mu[lay$pos] + as.numeric(C_pos %*% zpos)

    lp <- lp + .occu_cover_sp_ll(views[[s]], bocc, bp, bpos, log_disp)
    if (grad) {
      gvec <- .occu_cover_sp_grad(views[[s]], bocc, bp, bpos, log_disp)
      gocc <- gvec[lay$occ]; gp <- gvec[lay$p]; gpos <- gvec[lay$pos]
      g_ld <- g_ld + gvec[P + 1L]
      g_mu[lay$occ] <- g_mu[lay$occ] + gocc
      g_mu[lay$p]   <- g_mu[lay$p]   + gp
      g_mu[lay$pos] <- g_mu[lay$pos] + gpos
      g[bidx[lay$occ]] <- g[bidx[lay$occ]] + as.numeric(crossprod(C_occ, gocc))
      g[bidx[lay$p]]   <- g[bidx[lay$p]]   + as.numeric(crossprod(C_p,   gp))
      g[bidx[lay$pos]] <- g[bidx[lay$pos]] + as.numeric(crossprod(C_pos, gpos))
      A_occ <- A_occ + outer(gocc, zocc)
      A_p   <- A_p   + outer(gp,   zp)
      A_pos <- A_pos + outer(gpos, zpos)
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior, per arm ----
  arms <- list(list(chol = lay$chol_occ, A = A_occ, C = C_occ, Pa = lay$P_occ),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$P_p),
               list(chol = lay$chol_pos, A = A_pos, C = C_pos, Pa = lay$P_pos))
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

  # ---- shared log-dispersion prior (weakly-informative Normal on log scale) ----
  ld_mean <- priors$log_disp_mean %||% log(0.5)
  ld_sd   <- priors$log_disp_sd   %||% 2.0
  lp <- lp - 0.5 * ((log_disp - ld_mean) / ld_sd)^2
  g_ld <- g_ld - (log_disp - ld_mean) / ld_sd^2
  if (grad) g[lay$log_disp] <- g_ld

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Reconstruct the per-species deviation matrix b (S x P) from a packed coordinate
# vector under the non-centered map b_{s,arm} = C_arm z_{s,arm}.
.tobs_ms_occu_cover_nuts_b_from_z <- function(theta, lay) {
  C_occ <- .ms_ocs_chol_unpack(theta[lay$chol_occ], lay$P_occ)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$P_p)
  C_pos <- .ms_ocs_chol_unpack(theta[lay$chol_pos], lay$P_pos)
  B <- matrix(0, lay$n_species, lay$P)
  for (s in seq_len(lay$n_species)) {
    z <- theta[.ms_ocs_b_idx(lay, s)]
    B[s, lay$occ] <- as.numeric(C_occ %*% z[lay$occ])
    B[s, lay$p]   <- as.numeric(C_p   %*% z[lay$p])
    B[s, lay$pos] <- as.numeric(C_pos %*% z[lay$pos])
  }
  B
}


# Pack a community occu_cover Laplace-EM mode into the full NUTS coordinate
# vector: community means, the shared log-dispersion, the three community
# covariances as log-Cholesky coordinates, and the whitened per-species
# deviations z_s = C_arm^{-1} b_s. `Sigma` = list(occ=, p=, pos=); `b_list` a
# per-species list of full P-vectors (occ, p, pos coefficients concatenated).
.tobs_ms_occu_cover_nuts_pack_init <- function(mu, ld, Sigma, b_list, lay) {
  theta <- numeric(lay$total)
  theta[lay$mu]       <- as.numeric(mu)
  theta[lay$log_disp] <- as.numeric(ld)
  C_occ <- t(chol(.ms_ocs_pd(as.matrix(Sigma$occ))))
  C_p   <- t(chol(.ms_ocs_pd(as.matrix(Sigma$p))))
  C_pos <- t(chol(.ms_ocs_pd(as.matrix(Sigma$pos))))
  theta[lay$chol_occ] <- .ms_ocs_chol_pack(C_occ)
  theta[lay$chol_p]   <- .ms_ocs_chol_pack(C_p)
  theta[lay$chol_pos] <- .ms_ocs_chol_pack(C_pos)
  B <- do.call(rbind, b_list)                      # S x P
  for (s in seq_len(lay$n_species)) {
    z_s <- numeric(lay$P)
    z_s[lay$occ] <- forwardsolve(C_occ, B[s, lay$occ])
    z_s[lay$p]   <- forwardsolve(C_p,   B[s, lay$p])
    z_s[lay$pos] <- forwardsolve(C_pos, B[s, lay$pos])
    theta[.ms_ocs_b_idx(lay, s)] <- z_s
  }
  theta
}

# Marshal the C++ NUTS spec: the shared occ / detection / cover designs (visit
# blocks are 0-column matrices when absent) + per-species y / y_pos / valid
# [n_sites x max_visits] matrices + the positive-arm code.
.tobs_ms_occu_cover_nuts_spec <- function(model) {
  N <- model$n_sites; J <- model$max_visits; S <- model$n_species
  emptyv <- matrix(0, N * J, 0)
  Xdv <- if (is.null(model$X_det_visit)) emptyv else model$X_det_visit
  Xpv <- if (is.null(model$X_pos_visit)) emptyv else model$X_pos_visit
  intm <- function(m) { storage.mode(m) <- "integer"; m }
  list(n_sites = N, max_visits = J, n_species = S,
       pos_code = .occu_cover_pos_code(model$positive),
       X_occ = model$X_occ, X_det_site = model$X_det_site, X_det_visit = Xdv,
       X_pos_site = model$X_pos_site, X_pos_visit = Xpv,
       y     = lapply(seq_len(S), function(s) intm(model$y[, , s])),
       y_pos = lapply(seq_len(S), function(s) model$y_pos[, , s]),
       valid = lapply(seq_len(S), function(s) intm(model$valid[, , s])))
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the community joint occu+cover model
# ---------------------------------------------------------------------------

# Warm-start the community occu_cover Laplace-EM (the same engine the laplace path
# runs), pack the mode into the NUTS coordinate vector, and sample the exact joint
# posterior via the in-tree C++ FullGradFn (cpp_ms_occu_cover_nuts). Samples the
# community means, per-species deviations, the three community covariances, AND the
# shared log-dispersion jointly -> removes the Laplace variance attenuation. The
# returned fit reuses build_ms_occu_cover_fit so the tobs_fit surface matches the
# laplace path. Mirrors .tobs_fit_ms_int_occu_nuts (3 arms + a shared scalar).
.tobs_fit_ms_occu_cover_nuts <- function(model,
                                         sigma.beta = 5,
                                         n.iter = 1000L, n.warmup = 1000L,
                                         n.chains = 1L, max.treedepth = 10L,
                                         adapt.delta = 0.9, seed = 1L,
                                         max.iter = 200L, tol = 1e-4,
                                         verbose = FALSE, ...) {
  pil   <- model$process_info
  P_occ <- pil[[1L]]$p; P_p <- pil[[2L]]$p; P_pos <- pil[[3L]]$p
  P     <- P_occ + P_p + P_pos
  S     <- model$n_species
  arm_idx <- list(occ = seq_len(P_occ), p = P_occ + seq_len(P_p),
                  pos = P_occ + P_p + seq_len(P_pos))

  # Warm start at the community Laplace-EM mode (single source of truth).
  warm <- .tobs_fit_ms_occu_cover(model, sigma.beta = sigma.beta,
                                  max.iter = as.integer(max.iter),
                                  tol = as.numeric(tol), verbose = FALSE)
  ms  <- warm$ms_community
  mu  <- unname(warm$means[seq_len(P)])
  ld  <- unname(warm$means[P + 1L])
  Sigma_w <- list(occ = ms$Sigma_occ, p = ms$Sigma_p, pos = ms$Sigma_pos)
  b_list  <- lapply(seq_len(S), function(s)
    c(ms$blup_occ[s, ], ms$blup_p[s, ], ms$blup_pos[s, ]))

  lay    <- .tobs_ms_occu_cover_nuts_layout(P_occ, P_p, P_pos, S)
  pri    <- .ms_ocs_nuts_priors()
  theta0 <- .tobs_ms_occu_cover_nuts_pack_init(mu, ld, Sigma_w, b_list, lay)
  spec   <- .tobs_ms_occu_cover_nuts_spec(model)

  inv_metric <- .ms_ocs_fd_metric(
    function(th) cpp_ms_occu_cover_nuts_joint_logpost(spec, th, pri, sigma.beta)$grad,
    theta0)

  run_chain <- function(ch) {
    cpp_ms_occu_cover_nuts(
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

  # ---- reconstruct the EM-shaped outputs from the draws ----
  par     <- colMeans(draws)
  mu_hat  <- par[lay$mu]
  ld_hat  <- par[lay$log_disp]
  Vf      <- stats::cov(draws[, c(lay$mu, lay$log_disp), drop = FALSE])
  Sigma_hat <- list(
    occ = .ms_ocs_sig_mean(draws, lay$chol_occ, P_occ),
    p   = .ms_ocs_sig_mean(draws, lay$chol_p,   P_p),
    pos = .ms_ocs_sig_mean(draws, lay$chol_pos, P_pos))

  B_bar <- matrix(0, S, lay$P)
  for (i in seq_len(nrow(draws)))
    B_bar <- B_bar + .tobs_ms_occu_cover_nuts_b_from_z(draws[i, ], lay)
  B_bar <- B_bar / nrow(draws)
  b_list_hat <- lapply(seq_len(S), function(s) B_bar[s, ])

  # Data-only marginal log-lik at the posterior mean over reconstructed b_s.
  views <- lapply(seq_len(S), function(s) .ms_occu_cover_species_view(model, s))
  ll_mean <- 0
  for (s in seq_len(S)) {
    bs <- mu_hat + B_bar[s, ]
    ll_mean <- ll_mean + .occu_cover_sp_ll(views[[s]], bs[lay$occ], bs[lay$p],
                                           bs[lay$pos], ld_hat)
  }

  fit <- build_ms_occu_cover_fit(model, unname(mu_hat), unname(ld_hat), b_list_hat,
                                 Sigma_hat, Cinv_list = NULL, Vf, arm_idx,
                                 F_val = ll_mean, converged = TRUE,
                                 n_iter = warm$convergence$n_iter %||% NA_integer_,
                                 debias_method = "none")
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


# ---------------------------------------------------------------------------
# Per-species dispersion random effect (#115 B7 follow-up)
# ---------------------------------------------------------------------------
#
# The shared-dispersion target above carries ONE community log-dispersion scalar.
# The dispersion-RE variant instead gives each species its own log-dispersion
#   log_disp_s = mu_ld + sigma_ld * z_ld_s,   z_ld_s ~ N(0, 1)
# i.e. a FOURTH community arm of dimension 1 (the ms_abun log_r_s analogue), with
# community mean mu_ld and community SD sigma_ld. It reuses the same non-centered
# machinery: the 1x1 log-Cholesky "factor" is C_ld = exp(chol_ld) = sigma_ld, and
# the per-species dispersion score (.occu_cover_sp_grad's g_ld entry) chains to
# mu_ld / z_ld_s / the 1x1 covariance exactly like a coefficient arm. mu_ld carries
# the weakly-informative log-dispersion prior (not the coefficient ridge).

# Packed layout for the dispersion-RE variant:
#   theta = ( mu_coef [P_coef], mu_ld [1], {z_s} [S*(P_coef+1)],
#             chol_occ, chol_p, chol_pos, chol_ld [1] )
# with z_s = (z_occ, z_p, z_pos, z_ld). `mu_ld` is the community mean log-dispersion
# and `chol_ld` = log(sigma_ld).
.tobs_ms_occu_cover_re_disp_layout <- function(P_occ, P_p, P_pos, n_species) {
  P_coef <- P_occ + P_p + P_pos
  Pz     <- P_coef + 1L                          # per-species block width (+ z_ld)
  occ <- seq_len(P_occ)
  p   <- P_occ + seq_len(P_p)
  pos <- P_occ + P_p + seq_len(P_pos)
  ld  <- P_coef + 1L                             # z_ld position within a block
  mu_coef <- seq_len(P_coef)
  mu_ld   <- P_coef + 1L
  b_off <- P_coef + 1L                           # after (mu_coef, mu_ld)
  q_occ <- .ms_ocs_chol_dim(P_occ)
  q_p   <- .ms_ocs_chol_dim(P_p)
  q_pos <- .ms_ocs_chol_dim(P_pos)
  coff  <- (P_coef + 1L) + n_species * Pz
  chol_occ <- coff + seq_len(q_occ); coff <- coff + q_occ
  chol_p   <- coff + seq_len(q_p);   coff <- coff + q_p
  chol_pos <- coff + seq_len(q_pos); coff <- coff + q_pos
  chol_ld  <- coff + 1L; coff <- coff + 1L
  list(P_coef = P_coef, Pz = Pz, P_occ = P_occ, P_p = P_p, P_pos = P_pos,
       n_species = n_species, occ = occ, p = p, pos = pos, ld = ld,
       mu_coef = mu_coef, mu_ld = mu_ld, b_off = b_off,
       q_occ = q_occ, q_p = q_p, q_pos = q_pos,
       chol_occ = chol_occ, chol_p = chol_p, chol_pos = chol_pos,
       chol_ld = chol_ld, total = coff)
}

# Per-species z-block indices for the dispersion-RE layout (block width Pz).
.tobs_ms_occu_cover_re_disp_b_idx <- function(lay, s) {
  lay$b_off + (s - 1L) * lay$Pz + seq_len(lay$Pz)
}

# Joint log-posterior + gradient for the dispersion-RE variant. `views` the
# per-species occu_cover model views; `lay` the RE layout. Non-centered on all
# four arms (occ/p/pos coefficients + the 1-D log-dispersion).
.tobs_ms_occu_cover_re_disp_logpost <- function(theta, views, lay, priors,
                                                sigma.beta = 5, grad = TRUE) {
  S <- lay$n_species
  mu_coef <- theta[lay$mu_coef]
  mu_ld   <- theta[lay$mu_ld]
  g    <- numeric(lay$total)
  g_mc <- numeric(lay$P_coef)
  g_mld <- 0
  lp   <- 0

  C_occ <- .ms_ocs_chol_unpack(theta[lay$chol_occ], lay$P_occ)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$P_p)
  C_pos <- .ms_ocs_chol_unpack(theta[lay$chol_pos], lay$P_pos)
  sigma_ld <- exp(theta[lay$chol_ld])            # 1x1 Cholesky factor

  A_occ <- matrix(0, lay$P_occ, lay$P_occ)
  A_p   <- matrix(0, lay$P_p,   lay$P_p)
  A_pos <- matrix(0, lay$P_pos, lay$P_pos)
  A_ld  <- 0

  for (s in seq_len(S)) {
    bidx <- .tobs_ms_occu_cover_re_disp_b_idx(lay, s)
    z_s  <- theta[bidx]
    zocc <- z_s[lay$occ]; zp <- z_s[lay$p]; zpos <- z_s[lay$pos]; zld <- z_s[lay$ld]
    bocc <- mu_coef[lay$occ] + as.numeric(C_occ %*% zocc)
    bp   <- mu_coef[lay$p]   + as.numeric(C_p   %*% zp)
    bpos <- mu_coef[lay$pos] + as.numeric(C_pos %*% zpos)
    log_disp_s <- mu_ld + sigma_ld * zld

    lp <- lp + .occu_cover_sp_ll(views[[s]], bocc, bp, bpos, log_disp_s)
    if (grad) {
      gvec <- .occu_cover_sp_grad(views[[s]], bocc, bp, bpos, log_disp_s)
      gocc <- gvec[lay$occ]; gp <- gvec[lay$p]; gpos <- gvec[lay$pos]
      g_ld_s <- gvec[lay$P_coef + 1L]
      g_mc[lay$occ] <- g_mc[lay$occ] + gocc
      g_mc[lay$p]   <- g_mc[lay$p]   + gp
      g_mc[lay$pos] <- g_mc[lay$pos] + gpos
      g_mld <- g_mld + g_ld_s
      g[bidx[lay$occ]] <- g[bidx[lay$occ]] + as.numeric(crossprod(C_occ, gocc))
      g[bidx[lay$p]]   <- g[bidx[lay$p]]   + as.numeric(crossprod(C_p,   gp))
      g[bidx[lay$pos]] <- g[bidx[lay$pos]] + as.numeric(crossprod(C_pos, gpos))
      g[bidx[lay$ld]]  <- g[bidx[lay$ld]]  + sigma_ld * g_ld_s
      A_occ <- A_occ + outer(gocc, zocc)
      A_p   <- A_p   + outer(gp,   zp)
      A_pos <- A_pos + outer(gpos, zpos)
      A_ld  <- A_ld  + g_ld_s * zld
    }
  }

  # ---- z prior: standard normal over the whole per-species block ----
  z_idx <- lay$b_off + seq_len(S * lay$Pz)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior, per arm ----
  arms <- list(list(chol = lay$chol_occ, A = A_occ, C = C_occ, Pa = lay$P_occ),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$P_p),
               list(chol = lay$chol_pos, A = A_pos, C = C_pos, Pa = lay$P_pos),
               list(chol = lay$chol_ld,  A = matrix(A_ld, 1, 1),
                    C = matrix(sigma_ld, 1, 1), Pa = 1L))
  for (arm in arms) {
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], arm$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, arm$Pa) +
        pr$grad
  }

  # ---- community-mean priors: coefficients ridge + log-dispersion prior ----
  ib2 <- 1 / sigma.beta^2
  lp <- lp - 0.5 * ib2 * sum(mu_coef^2)
  g_mc <- g_mc - ib2 * mu_coef
  ld_mean <- priors$log_disp_mean %||% log(0.5)
  ld_sd   <- priors$log_disp_sd   %||% 2.0
  lp <- lp - 0.5 * ((mu_ld - ld_mean) / ld_sd)^2
  g_mld <- g_mld - (mu_ld - ld_mean) / ld_sd^2
  if (grad) {
    g[lay$mu_coef] <- g[lay$mu_coef] + g_mc
    g[lay$mu_ld]   <- g_mld
  }

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}
