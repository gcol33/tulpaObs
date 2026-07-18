# ms_occu_pg_gibbs.R - Polya-Gamma Gibbs for COMMUNITY single-season occupancy
# (spOccupancy msPGOcc; gcol33/tulpaObs#115, #126). The hierarchical extension of
# the single-species PGOcc engine (R/occu_pg_gibbs.R): per-species occupancy /
# detection coefficients with Gaussian community hyperpriors,
#
#   logit psi_{s,i} = X_occ_i . beta_psi_s,   beta_psi_s ~ N(mu_psi, diag(tau_psi^2))
#   logit p_{s,i}   = X_det_i . beta_p_s,      beta_p_s   ~ N(mu_p,   diag(tau_p^2))
#
# with community means mu ~ N(0, sigma.mu^2) and independent community variances
# tau^2 ~ Inverse-Gamma(a, b) (the diagonal community covariance spOccupancy's
# msPGOcc uses). Conditional on the PG auxiliaries every coefficient update is
# exactly conjugate Gaussian, so a Gibbs sweep is: per species sample z, then the
# PG-augmented conjugate beta_s; then the conjugate community mean mu and the
# Inverse-Gamma community variance tau^2 per coordinate. This gives a calibrated
# community-VARIANCE posterior -- the Laplace-EM (R/community_em.R) leaves those
# components with small-cluster attenuation (a documented lower bound); the Gibbs
# does not. PG draws use tulpa's Polson-Scott-Windle sampler (tulpa:::cpp_rpg).
#
# v1: single-season community occupancy, site-level detection, no structured
# terms / spatial field (the shared-field sfMsPGOcc is a follow-up).

.tobs_fit_ms_occu_pg_gibbs <- function(model, priors = NULL, sigma.beta = 2.5,
                                       n.iter = 3000L, n.warmup = 1500L,
                                       n.chains = 2L, n.thin = 1L, seed = 1L,
                                       verbose = FALSE, ...) {
  rpg    <- get("cpp_rpg", envir = asNamespace("tulpa"))
  X_psi  <- model$X_occ
  X_p    <- model$X_det
  y      <- model$y                              # [n x mv x S]
  valid  <- model$valid
  n      <- nrow(X_psi); S <- model$n_species
  p_psi  <- ncol(X_psi); p_p <- ncol(X_p)

  # Per-species site sufficient statistics.
  kdet <- nvis <- matrix(0L, n, S); anydet <- matrix(FALSE, n, S)
  for (s in seq_len(S)) {
    ys <- y[, , s]; vs <- valid[, , s]
    nvis[, s]   <- rowSums(vs)
    kdet[, s]   <- rowSums(ys == 1L & vs)
    anydet[, s] <- kdet[, s] > 0L
  }

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  par_names <- c(paste0("psi_", model$process_info[[1L]]$coef_names),
                 paste0("p_",   model$process_info[[2L]]$coef_names))

  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    b_psi <- matrix(stats::rnorm(S * p_psi, 0, 0.2), S, p_psi)
    b_p   <- matrix(stats::rnorm(S * p_p,   0, 0.2), S, p_p)
    mu_psi <- rep(0, p_psi); tau2_psi <- rep(1, p_psi)
    mu_p   <- rep(0, p_p);   tau2_p   <- rep(1, p_p)
    mu_draws  <- matrix(NA_real_, n_keep, p_psi + p_p)
    tau_draws <- matrix(NA_real_, n_keep, p_psi + p_p)
    b_psi_sum <- matrix(0, S, p_psi); b_p_sum <- matrix(0, S, p_p); nsum <- 0L
    ki <- 0L
    for (it in seq_len(n.iter)) {
      for (s in seq_len(S)) {
        eta_psi <- as.vector(X_psi %*% b_psi[s, ]); psi <- stats::plogis(eta_psi)
        eta_p   <- as.vector(X_p   %*% b_p[s, ]);   pd  <- stats::plogis(eta_p)
        # latent occupancy
        z <- integer(n); z[anydet[, s]] <- 1L
        und <- !anydet[, s]
        if (any(und)) {
          l1 <- psi[und] * (1 - pd[und])^nvis[und, s]; l0 <- 1 - psi[und]
          z[und] <- stats::rbinom(sum(und), 1L, l1 / (l1 + l0))
        }
        # occupancy coefficients
        om <- rpg(rep(1, n), eta_psi)
        b_psi[s, ] <- .tobs_pg_draw_beta(X_psi, om, z - 0.5,
                                         1 / tau2_psi, mu_psi / tau2_psi)
        # detection coefficients at occupied sites
        occ <- which(z == 1L & nvis[, s] > 0L)
        if (length(occ) >= p_p) {
          Xo <- X_p[occ, , drop = FALSE]
          om_p <- rpg(nvis[occ, s], as.vector(Xo %*% b_p[s, ]))
          b_p[s, ] <- .tobs_pg_draw_beta(Xo, om_p, kdet[occ, s] - nvis[occ, s] / 2,
                                         1 / tau2_p, mu_p / tau2_p)
        }
      }
      # community means + variances (per coordinate, both arms)
      cu <- .tobs_pg_community_update(b_psi, mu_psi, tau2_psi, S)
      mu_psi <- cu$mu; tau2_psi <- cu$tau2
      cu <- .tobs_pg_community_update(b_p, mu_p, tau2_p, S)
      mu_p <- cu$mu; tau2_p <- cu$tau2
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L
        mu_draws[ki, ]  <- c(mu_psi, mu_p)
        tau_draws[ki, ] <- sqrt(c(tau2_psi, tau2_p))
        b_psi_sum <- b_psi_sum + b_psi; b_p_sum <- b_p_sum + b_p; nsum <- nsum + 1L
      }
    }
    list(mu = mu_draws, tau = tau_draws,
         b_psi = b_psi_sum / nsum, b_p = b_p_sum / nsum)
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  summ <- .tobs_pg_summarize(lapply(chains, `[[`, "mu"), par_names)
  means <- summ$means

  # Community SD (a derived quantity = sqrt(tau^2)): report the posterior MEDIAN,
  # robust to the right-skew of the variance-component posterior at moderate S
  # (the mean over-states it; marginalize-derived-quantities rule).
  tau_all <- do.call(rbind, lapply(chains, function(c) c$tau))
  tau_med <- apply(tau_all, 2L, stats::median)
  sd_psi <- tau_med[seq_len(p_psi)]
  sd_p   <- tau_med[p_psi + seq_len(p_p)]
  names(sd_psi) <- model$process_info[[1L]]$coef_names
  names(sd_p)   <- model$process_info[[2L]]$coef_names

  # Per-species coefficients = chain-averaged posterior-mean betas.
  coef_psi <- Reduce(`+`, lapply(chains, `[[`, "b_psi")) / n.chains
  coef_p   <- Reduce(`+`, lapply(chains, `[[`, "b_p"))   / n.chains
  rownames(coef_psi) <- rownames(coef_p) <- model$species_names
  colnames(coef_psi) <- model$process_info[[1L]]$coef_names
  colnames(coef_p)   <- model$process_info[[2L]]$coef_names
  blup_psi <- sweep(coef_psi, 2L, means[seq_len(p_psi)], "-")
  blup_p   <- sweep(coef_p,   2L, means[p_psi + seq_len(p_p)], "-")

  .tobs_pg_finalize_fit(
    summ, par_names, model, model$process_info, N = sum(model$valid),
    n.iter = n.iter, n.chains = n.chains,
    extra = list(ms_community = list(
      Sigma_psi = diag(sd_psi^2, p_psi), Sigma_p = diag(sd_p^2, p_p),
      sd_psi = sd_psi, sd_p = sd_p,
      coef_psi = coef_psi, coef_p = coef_p,
      blup_psi = blup_psi, blup_p = blup_p)))
}
