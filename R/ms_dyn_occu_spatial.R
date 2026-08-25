# =============================================================================
# ms_dyn_occu_spatial.R - community DYNAMIC occupancy with a shared areal field
# on the first-season occupancy (psi1) arm (spOccupancy stMsPGOcc /
# svcTMsPGOcc).
#
#   logit psi1_{s,i} = X_i . (mu_psi1 + b_psi1_s) + sum_k W[i,k] F[u(i),k]
#   logit p_{s,i}    = Xdet_i . (mu_p + b_p_s)
#   logit gamma_i / eps_i = shared community transition coefficients
#
# The field is shared across species and enters season-1 occupancy only. Because
# psi1 sets ONLY the initial mixing weight of the per-(species, site) HMM, the
# marginal factorises as a two-component mixture in psi1:
#
#   L_{s,i}(psi1) = (1 - psi1) A_{s,i} + psi1 B_{s,i}
#
# where A_{s,i} = P(data | z_1 = 0) and B_{s,i} = P(data | z_1 = 1) are the HMM
# forward likelihoods conditional on the first-season state (independent of
# psi1). This is the SAME mixture the single-season community occupancy field
# oracle uses (R/ms_occu_field.R), so the Louis score / curvature wrt a psi1
# offset carry over verbatim -- only A / B come from the HMM forward instead of
# the single-season emission. The block-coordinate driver (R/community_latent.R)
# then alternates the community EM (field as a psi1 offset) with the areal field
# Newton, exactly as every other community field family. icar only.
# =============================================================================


# Per-(site, season) detection sufficient statistics for one species, from its
# [n_sites x max_visits x n_seasons] detection (ys) / validity (vs) arrays:
# nvalid[i, t] = observed visits, ndet[i, t] = detections. Precomputed once per
# species (they do not depend on any parameter).
.ms_dyn_occu_emit_stats <- function(ys, vs, n_sites, n_seasons) {
  nvalid <- matrix(0, n_sites, n_seasons); ndet <- matrix(0, n_sites, n_seasons)
  for (t in seq_len(n_seasons)) {
    vt <- vs[, , t]; yt <- ys[, , t]; yt[!vt] <- 0L
    nvalid[, t] <- rowSums(vt); ndet[, t] <- rowSums(yt)
  }
  list(nvalid = nvalid, ndet = ndet)
}

# Vectorised (over sites) per-season emissions for one species given a per-site
# detection prob p: em_occ[i, t] = p_i^ndet (1 - p_i)^(nvalid - ndet), and
# em_unocc[i, t] = 1 if the season has no detection else 0.
.ms_dyn_occu_emissions <- function(p, nvalid, ndet) {
  lp <- log(pmin(pmax(p, 1e-12), 1 - 1e-12)); l1p <- log(1 - pmin(pmax(p, 1e-12), 1 - 1e-12))
  em_occ   <- exp(ndet * lp + (nvalid - ndet) * l1p)
  em_unocc <- ifelse(ndet > 0, 0, 1)
  list(occ = em_occ, unocc = em_unocc)
}

# HMM forward likelihoods conditional on the first-season state, VECTORISED over
# sites: logA[i] = log P(data_i | z_{i,1} = 0), logB[i] = log P(data_i |
# z_{i,1} = 1). Independent of psi1 (the oracle mixes them). A season whose
# conditional mass hits zero (e.g. a detection under z_1 = 0 with no
# colonisation route) carries -Inf, correctly zeroing that branch.
.ms_dyn_occu_condAB_vec <- function(p, gamma, eps, em, n_sites, n_seasons) {
  gamma <- pmin(pmax(gamma, 1e-12), 1 - 1e-12)
  eps   <- pmin(pmax(eps,   1e-12), 1 - 1e-12)
  em_occ <- em$occ; em_unocc <- em$unocc
  fwd_cond <- function(a0, a1) {                       # seeded season-1 states
    v0 <- a0 * em_unocc[, 1]; v1 <- a1 * em_occ[, 1]
    ct <- v0 + v1
    ll <- ifelse(ct > 0, log(ct), -Inf)
    a0 <- ifelse(ct > 0, v0 / ct, 0); a1 <- ifelse(ct > 0, v1 / ct, 0)
    for (t in 2:n_seasons) {
      pr1 <- a1 * (1 - eps) + a0 * gamma
      pr0 <- a1 * eps       + a0 * (1 - gamma)
      v0 <- pr0 * em_unocc[, t]; v1 <- pr1 * em_occ[, t]
      ct <- v0 + v1
      ll <- ll + ifelse(ct > 0, log(ct), -Inf)
      a0 <- ifelse(ct > 0, v0 / ct, 0); a1 <- ifelse(ct > 0, v1 / ct, 0)
    }
    ll
  }
  one <- rep(1, n_sites); zero <- rep(0, n_sites)
  list(logA = fwd_cond(one, zero), logB = fwd_cond(zero, one))
}


# The per-species HMM forward marginal log-likelihood, summed over sites: the
# one kernel every ms_dyn_occu route evaluates (laplace, NUTS, spatial, SBC).
# The latent occupancy path z integrates out by a scaled forward filter,
# vectorised over sites. Season 1 starts from the prior c(1 - psi1, psi1); each
# season multiplies by the emission (unoccupied: 1 if the season carries no
# detection else 0; occupied: p^ndet (1 - p)^(nvalid - ndet), so a season with
# no valid visits is uninformative at 1); between seasons the per-site
# gamma_i / eps_i transition applies. The log of each season's forward
# normalizer accumulates into the log-likelihood.
#
# psi1 / gamma / eps are per-site vectors (length n_sites); the detection
# probability enters only through `em`, the emission pair from
# .ms_dyn_occu_emissions().
#
# A site whose forward mass underflows to zero contributes -Inf and takes the
# species with it. The likelihood is over every site it was handed: a step that
# kills a site is one the EM line search must reject, not one it accepts on the
# surviving subset (gcol33/tulpaObs#259).
.ms_dyn_occu_fwd_ll_vec <- function(psi1, gamma, eps, em, n_sites, n_seasons) {
  psi1  <- pmin(pmax(psi1,  1e-12), 1 - 1e-12)
  gamma <- pmin(pmax(gamma, 1e-12), 1 - 1e-12)
  eps   <- pmin(pmax(eps,   1e-12), 1 - 1e-12)
  em_occ <- em$occ; em_unocc <- em$unocc
  v0 <- (1 - psi1) * em_unocc[, 1]; v1 <- psi1 * em_occ[, 1]
  ct <- v0 + v1
  ll <- ifelse(ct > 0, log(ct), -Inf)
  a0 <- ifelse(ct > 0, v0 / ct, 0); a1 <- ifelse(ct > 0, v1 / ct, 0)
  for (t in seq_len(n_seasons - 1L) + 1L) {
    pr1 <- a1 * (1 - eps) + a0 * gamma
    pr0 <- a1 * eps       + a0 * (1 - gamma)
    v0 <- pr0 * em_unocc[, t]; v1 <- pr1 * em_occ[, t]
    ct <- v0 + v1
    ll <- ll + ifelse(ct > 0, log(ct), -Inf)
    a0 <- ifelse(ct > 0, v0 / ct, 0); a1 <- ifelse(ct > 0, v1 / ct, 0)
  }
  sum(ll)
}


# Vectorised forward-backward smoothing for one species. Returns the per-site
# marginal log-likelihood, the season-1 smoothed posterior w1[i] = P(z_1 = 1 |
# y), the detection score sum_t w_it (ndet_it - nvalid_it p_i), and the
# aggregated colonisation / extinction sufficient statistics (col_y = sum_t
# xi01_t etc.). These are the pieces of the Fisher-identity gradient of the HMM
# marginal wrt the four site-level linear predictors.
.ms_dyn_occu_fb_vec <- function(psi1, p, gamma, eps, em, nvalid, ndet,
                                n_sites, n_seasons) {
  psi1  <- pmin(pmax(psi1,  1e-12), 1 - 1e-12)
  p     <- pmin(pmax(p,     1e-12), 1 - 1e-12)
  gamma <- pmin(pmax(gamma, 1e-12), 1 - 1e-12)
  eps   <- pmin(pmax(eps,   1e-12), 1 - 1e-12)
  em_occ <- em$occ; em_unocc <- em$unocc
  A0 <- matrix(0, n_sites, n_seasons); A1 <- matrix(0, n_sites, n_seasons)
  cs <- matrix(1, n_sites, n_seasons)
  # forward (scaled), storing the filtered a and normalisers.
  v0 <- (1 - psi1) * em_unocc[, 1]; v1 <- psi1 * em_occ[, 1]
  ct <- v0 + v1; cs[, 1] <- ct
  a0 <- ifelse(ct > 0, v0 / ct, 0); a1 <- ifelse(ct > 0, v1 / ct, 0)
  A0[, 1] <- a0; A1[, 1] <- a1
  for (t in 2:n_seasons) {
    pr1 <- a1 * (1 - eps) + a0 * gamma
    pr0 <- a1 * eps       + a0 * (1 - gamma)
    v0 <- pr0 * em_unocc[, t]; v1 <- pr1 * em_occ[, t]
    ct <- v0 + v1; cs[, t] <- ct
    a0 <- ifelse(ct > 0, v0 / ct, 0); a1 <- ifelse(ct > 0, v1 / ct, 0)
    A0[, t] <- a0; A1[, t] <- a1
  }
  ll_site <- rowSums(ifelse(cs > 0, log(cs), -Inf))
  # backward (scaled) + smoothed marginals w and pairwise joints xi.
  w  <- matrix(0, n_sites, n_seasons); w[, n_seasons] <- A1[, n_seasons]
  bw0 <- rep(1, n_sites); bw1 <- rep(1, n_sites)
  col_y <- numeric(n_sites); col_n <- numeric(n_sites)
  ext_y <- numeric(n_sites); ext_n <- numeric(n_sites)
  for (t in (n_seasons - 1):1) {
    bb0 <- em_unocc[, t + 1] * bw0; bb1 <- em_occ[, t + 1] * bw1
    invc <- ifelse(cs[, t + 1] > 0, 1 / cs[, t + 1], 0)
    xi01 <- A0[, t] * gamma       * bb1 * invc          # colonisation
    xi00 <- A0[, t] * (1 - gamma) * bb0 * invc
    xi10 <- A1[, t] * eps         * bb0 * invc          # extinction
    xi11 <- A1[, t] * (1 - eps)   * bb1 * invc
    col_y <- col_y + xi01; col_n <- col_n + (xi00 + xi01)
    ext_y <- ext_y + xi10; ext_n <- ext_n + (xi10 + xi11)
    nb0 <- ((1 - gamma) * bb0 + gamma       * bb1) * invc
    nb1 <- (eps         * bb0 + (1 - eps)   * bb1) * invc
    bw0 <- nb0; bw1 <- nb1
    w[, t] <- A1[, t] * bw1
  }
  # detection score: sum_t P(z_it = 1 | y) (ndet_it - nvalid_it p_i).
  p_score <- rowSums(w * (ndet - nvalid * p))
  list(ll = sum(ll_site[is.finite(ll_site)]), w1 = w[, 1L], p_score = p_score,
       col_y = col_y, col_n = col_n, ext_y = ext_y, ext_n = ext_n)
}


# Community dynamic-occupancy working oracle for the psi1 arm. `logA` / `logB`
# are [n_sites x n_species] conditional HMM log-likelihoods (held fixed across a
# field update, since psi1 only re-mixes them). The marginal is
# L = (1 - psi1) A + psi1 B, and integrating the latent season-1 state gives the
# Louis score / observed information wrt an additive offset on logit psi1:
#
#   r     = psi1 B / L            (posterior P(z_1 = 1 | data))
#   score = r - psi1             (= psi1 (1 - psi1)(B - A) / L)
#   curv  = psi1 (1 - psi1) - r (1 - r)   (>= 0, the missing-information subtract)
#
# identical to the single-season two-state oracle with A = 1 / B = q. Everything
# is computed in log space for numerical stability.
.tobs_ms_dyn_occu_oracle <- function(logA, logB) {
  Ns <- nrow(logA); S <- ncol(logA)
  # log L and log r via log-sum-exp of the two mixture branches.
  branch <- function(eta) {
    lpsi  <- stats::plogis(eta,  log.p = TRUE)      # log psi1
    l1psi <- stats::plogis(-eta, log.p = TRUE)      # log(1 - psi1)
    lo <- l1psi + logA                              # log((1-psi1) A)
    hi <- lpsi  + logB                              # log(psi1 B)
    m  <- pmax(lo, hi)
    logL <- m + log(exp(lo - m) + exp(hi - m))      # log L
    r    <- exp(hi - logL)                          # psi1 B / L in (0, 1)
    list(logL = logL, r = pmin(pmax(r, 0), 1), psi = stats::plogis(eta))
  }
  list(
    n_sites   = Ns,
    n_species = S,
    working = function(eta) {
      b <- branch(eta)
      list(score = b$r - b$psi,
           curv  = pmax(b$psi * (1 - b$psi) - b$r * (1 - b$r), 1e-8))
    },
    # Undetected sites drive the observed curvature toward zero; the complete-data
    # Fisher curvature psi (1 - psi) keeps the tau M-step covariance conditioned.
    cov_curv = function(eta) {
      psi <- stats::plogis(eta); psi * (1 - psi)
    },
    # Per-(site, species) marginal log-likelihood. The joint site marginal that
    # sets the factor magnitude integrates zeta with the species at a site kept
    # together, so it needs the cells rather than their sum.
    ll_cell = function(eta) branch(eta)$logL,
    data_ll = function(eta) sum(branch(eta)$logL))
}


# Fit the community dynamic occupancy model with a shared areal field on the
# first-season occupancy arm. Mirrors .tobs_fit_ms_occu_field: the community EM
# (with the field as a psi1 offset) alternates with the areal field Newton via
# the shared block-coordinate driver. icar only.
.tobs_fit_ms_dyn_occu_field <- function(model, spatial,
                                        max.iter = 200L, tol = 1e-4,
                                        sigma.beta = 5, priors = NULL,
                                        max.outer = NULL, verbose = FALSE, ...) {
  pi_list <- model$process_info
  P_psi1 <- pi_list[[1L]]$p; P_p <- pi_list[[2L]]$p
  P_gam  <- pi_list[[3L]]$p; P_eps <- pi_list[[4L]]$p
  P <- P_psi1 + P_p; G <- P_gam + P_eps
  S <- model$n_species; Ns <- model$n_sites; n_seasons <- model$n_seasons

  psi1_idx <- seq_len(P_psi1); p_idx <- P_psi1 + seq_len(P_p)
  gam_idx  <- seq_len(P_gam);  eps_idx <- P_gam + seq_len(P_eps)
  arm_idx  <- list(psi1 = psi1_idx, p = p_idx)

  X_psi1 <- model$X_psi1; X_p <- model$X_p
  X_gamma <- model$X_gamma; X_eps <- model$X_eps
  ys_list <- lapply(seq_len(S), function(s) model$y[, , , s])
  vs_list <- lapply(seq_len(S), function(s) model$valid[, , , s])
  # Per-species per-(site, season) detection sufficient stats (fixed).
  em_stats <- lapply(seq_len(S), function(s)
    .ms_dyn_occu_emit_stats(ys_list[[s]], vs_list[[s]], Ns, n_seasons))

  # ---- warm start (mirrors .tobs_fit_ms_dyn_occu) ----
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

  em_fit <- function(site_off, fac_off, em_prev) {
    # Marginal log-lik and its analytic Fisher-identity gradient share one
    # vectorised forward-backward per species, so the community EM skips its
    # O(U^2) finite-difference Hessian (it finite-differences this cheap
    # analytic gradient instead of the expensive marginal).
    eta_of <- function(s, theta, global) {
      list(psi1 = stats::plogis(as.numeric(X_psi1 %*% theta[psi1_idx]) + site_off),
           p    = stats::plogis(as.numeric(X_p    %*% theta[p_idx])),
           gam  = stats::plogis(as.numeric(X_gamma %*% global[gam_idx])),
           eps  = stats::plogis(as.numeric(X_eps   %*% global[eps_idx])))
    }
    sp_ll <- function(s, theta, global) {
      e  <- eta_of(s, theta, global)
      em <- .ms_dyn_occu_emissions(e$p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
      .ms_dyn_occu_fwd_ll_vec(e$psi1, e$gam, e$eps, em, Ns, n_seasons)
    }
    sp_grad <- function(s, theta, global) {
      e  <- eta_of(s, theta, global)
      em <- .ms_dyn_occu_emissions(e$p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
      fb <- .ms_dyn_occu_fb_vec(e$psi1, e$p, e$gam, e$eps, em,
                                em_stats[[s]]$nvalid, em_stats[[s]]$ndet,
                                Ns, n_seasons)
      c(as.numeric(crossprod(X_psi1, fb$w1 - e$psi1)),
        as.numeric(crossprod(X_p,    fb$p_score)),
        as.numeric(crossprod(X_gamma, fb$col_y - e$gam * fb$col_n)),
        as.numeric(crossprod(X_eps,   fb$ext_y - e$eps * fb$ext_n)))
    }
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em_prev)) init_mu else em_prev$mu,
      init_b = em_prev$b_list, init_Sigma = em_prev$Sigma,
      init_global = if (is.null(em_prev)) init_global else em_prev$global,
      penalize_global = TRUE, sigma_beta = sigma.beta, priors = priors,
      sigma_init = 0.3, max_iter = min(as.integer(max.iter), 30L),
      tol = as.numeric(tol), newton_max = 20L, verbose = FALSE)
  }
  offset_of <- function(em) {
    vapply(seq_len(S), function(s)
      as.numeric(X_psi1 %*% (em$mu[psi1_idx] + em$b_list[[s]][psi1_idx])),
      numeric(Ns))
  }
  make_oracle <- function(em) {
    gam <- stats::plogis(as.numeric(X_gamma %*% em$global[gam_idx]))
    eps <- stats::plogis(as.numeric(X_eps   %*% em$global[eps_idx]))
    logA <- matrix(0, Ns, S); logB <- matrix(0, Ns, S)
    for (s in seq_len(S)) {
      p_s <- stats::plogis(as.numeric(X_p %*% (em$mu[p_idx] + em$b_list[[s]][p_idx])))
      em_s <- .ms_dyn_occu_emissions(p_s, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
      ab  <- .ms_dyn_occu_condAB_vec(p_s, gam, eps, em_s, Ns, n_seasons)
      logA[, s] <- ab$logA; logB[, s] <- ab$logB
    }
    .tobs_ms_dyn_occu_oracle(logA, logB)
  }

  res <- .tobs_community_latent_ascent(
    spatial = spatial, latent = NULL, model = model, what = "ms_dyn_occu()",
    make_oracle = make_oracle, em_fit = em_fit, offset_of = offset_of,
    allow = "icar", tol = tol, max.outer = max.outer, verbose = verbose)

  fit <- build_ms_dyn_occu_fit(model, res$em, arm_idx, gam_idx, eps_idx)
  fit$method <- "laplace"
  fit <- .tobs_latent_attach_field(fit, res, spatial, "psi1_field_offset")
  fit
}
