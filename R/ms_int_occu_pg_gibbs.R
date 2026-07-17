# ms_int_occu_pg_gibbs.R - Polya-Gamma Gibbs for COMMUNITY multi-source
# integrated occupancy (spOccupancy-style; gcol33/tulpaObs#115, #126). The
# community integrated extension of msPGOcc (R/ms_occu_pg_gibbs.R): a single
# latent occupancy state per (species, site) observed by D detection sources,
#
#   logit psi_{s,i}   = X_psi_i . beta_psi_s,   beta_psi_s ~ N(mu_psi, diag(tau_psi^2))
#   logit p_{s,i,d}   = X_pd_i  . beta_pd_s,     beta_pd_s  ~ N(mu_pd,  diag(tau_pd^2))
#
# with community means mu ~ N(0, sigma.mu^2) and diagonal community variances
# tau^2 ~ Inverse-Gamma(a, b). Conditional on the Polya-Gamma auxiliaries every
# coefficient update is exactly conjugate Gaussian: per species sample the single
# latent z (occupied if any source detects, else Bernoulli on the pooled
# occupied-undetected mass), then the PG-augmented conjugate beta_psi_s and the D
# per-source beta_pd_s (each at that species' occupied sites the source covers);
# then the conjugate community mean + Inverse-Gamma variance per coordinate per
# arm. This gives a calibrated community-VARIANCE posterior -- the community
# Laplace-EM leaves those components with small-cluster attenuation. PG draws use
# tulpa's Polson-Scott-Windle sampler (tulpa:::cpp_rpg).
#
# v1: full/partial source coverage (site_maps folded into the per-species n_valid
# summaries -- an uncovered (site, source) has n_valid = 0 and drops out), site-
# level per-source detection, constant/covariate arms, no spatial field.

.tobs_fit_ms_int_occu_pg_gibbs <- function(model, priors = NULL, sigma.beta = 2.5,
                                           n.iter = 3000L, n.warmup = 1500L,
                                           n.chains = 2L, n.thin = 1L, seed = 1L,
                                           verbose = FALSE, ...) {
  rpg   <- get("cpp_rpg", envir = asNamespace("tulpa"))
  X_psi <- model$X_psi; X_p <- model$X_p            # X_p: list of D designs
  D     <- model$n_sources; S <- model$n_species
  n     <- model$n_sites; p_psi <- ncol(X_psi)
  p_pd  <- vapply(X_p, ncol, integer(1))            # per-source detection dims
  summ  <- model$summaries                          # per-species suff. stats
  proc  <- model$process_info
  arm_names <- c("psi", model$process_names)        # psi, p1, ..., pD

  sigma.mu2 <- 100; ig_a <- 0.1; ig_b <- 0.1

  draw_beta <- function(X, omega, kappa, mu, tau2) {
    XtOX <- crossprod(X, X * omega); diag(XtOX) <- diag(XtOX) + 1 / tau2
    V <- chol2inv(chol(XtOX))
    m <- V %*% (crossprod(X, kappa) + mu / tau2)
    as.vector(m + t(chol(V)) %*% stats::rnorm(ncol(X)))
  }

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  par_names <- c(paste0("psi_", proc[[1L]]$coef_names),
                 unlist(lapply(seq_len(D), function(d)
                   paste0(proc[[1L + d]]$name, "_", proc[[1L + d]]$coef_names))))
  p_tot <- p_psi + sum(p_pd)

  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    b_psi <- matrix(stats::rnorm(S * p_psi, 0, 0.2), S, p_psi)
    b_p   <- lapply(seq_len(D), function(d)
      matrix(stats::rnorm(S * p_pd[d], 0, 0.2), S, p_pd[d]))
    mu_psi <- rep(0, p_psi); tau2_psi <- rep(1, p_psi)
    mu_p   <- lapply(p_pd, function(pd) rep(0, pd))
    tau2_p <- lapply(p_pd, function(pd) rep(1, pd))
    mu_draws  <- matrix(NA_real_, n_keep, p_tot)
    tau_draws <- matrix(NA_real_, n_keep, p_tot)
    b_psi_sum <- matrix(0, S, p_psi)
    b_p_sum   <- lapply(p_pd, function(pd) matrix(0, S, pd))
    nsum <- 0L; ki <- 0L
    for (it in seq_len(n.iter)) {
      for (s in seq_len(S)) {
        sm <- summ[[s]]; nv <- sm$n_valid; nd_mat <- sm$n_det; anyd <- sm$any_det
        eta_psi <- as.vector(X_psi %*% b_psi[s, ]); psi <- stats::plogis(eta_psi)
        pd_list <- lapply(seq_len(D), function(d)
          stats::plogis(as.vector(X_p[[d]] %*% b_p[[d]][s, ])))
        # Latent occupancy: a detection at any source forces z = 1; else the
        # occupied-undetected mass pools (1 - p_d)^{n_valid_d} across sources.
        z <- integer(n); z[anyd] <- 1L; und <- !anyd
        if (any(und)) {
          logmass <- log(pmax(psi[und], 1e-12))
          for (d in seq_len(D))
            logmass <- logmass + nv[und, d] * log1p(-pmin(pd_list[[d]][und], 1 - 1e-12))
          l1 <- exp(logmass); l0 <- 1 - psi[und]
          z[und] <- stats::rbinom(sum(und), 1L, l1 / (l1 + l0))
        }
        # Occupancy coefficients (all sites).
        om <- rpg(rep(1, n), eta_psi)
        b_psi[s, ] <- draw_beta(X_psi, om, z - 0.5, mu_psi, tau2_psi)
        # Per-source detection coefficients at this species' occupied, covered
        # sites (n_valid_d > 0).
        for (d in seq_len(D)) {
          occ <- which(z == 1L & nv[, d] > 0L)
          if (length(occ) >= p_pd[d]) {
            Xo <- X_p[[d]][occ, , drop = FALSE]
            om_p <- rpg(nv[occ, d], as.vector(Xo %*% b_p[[d]][s, ]))
            b_p[[d]][s, ] <- draw_beta(Xo, om_p, nd_mat[occ, d] - nv[occ, d] / 2,
                                       mu_p[[d]], tau2_p[[d]])
          }
        }
      }
      # Community means + Inverse-Gamma variances, per coordinate, per arm.
      for (j in seq_len(p_psi)) {
        vj <- 1 / (S / tau2_psi[j] + 1 / sigma.mu2)
        mu_psi[j] <- stats::rnorm(1, vj * sum(b_psi[, j]) / tau2_psi[j], sqrt(vj))
        tau2_psi[j] <- 1 / stats::rgamma(1, ig_a + S / 2,
                                         ig_b + 0.5 * sum((b_psi[, j] - mu_psi[j])^2))
      }
      for (d in seq_len(D)) for (j in seq_len(p_pd[d])) {
        vj <- 1 / (S / tau2_p[[d]][j] + 1 / sigma.mu2)
        mu_p[[d]][j] <- stats::rnorm(1, vj * sum(b_p[[d]][, j]) / tau2_p[[d]][j], sqrt(vj))
        tau2_p[[d]][j] <- 1 / stats::rgamma(1, ig_a + S / 2,
                                            ig_b + 0.5 * sum((b_p[[d]][, j] - mu_p[[d]][j])^2))
      }
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L
        mu_draws[ki, ]  <- c(mu_psi, unlist(mu_p))
        tau_draws[ki, ] <- sqrt(c(tau2_psi, unlist(tau2_p)))
        b_psi_sum <- b_psi_sum + b_psi
        for (d in seq_len(D)) b_p_sum[[d]] <- b_p_sum[[d]] + b_p[[d]]
        nsum <- nsum + 1L
      }
    }
    list(mu = mu_draws, tau = tau_draws, b_psi = b_psi_sum / nsum,
         b_p = lapply(b_p_sum, function(m) m / nsum))
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  mu_chains <- lapply(chains, function(c) { colnames(c$mu) <- par_names; c$mu })

  draws <- do.call(rbind, mu_chains)
  means <- colMeans(draws); names(means) <- par_names
  V <- stats::cov(draws); dimnames(V) <- list(par_names, par_names)
  sds <- apply(draws, 2L, stats::sd); names(sds) <- par_names

  # Community SD (derived): posterior MEDIAN, robust to the variance-component
  # right-skew at moderate S (marginalize-derived-quantities rule).
  tau_all <- do.call(rbind, lapply(chains, function(c) c$tau))
  tau_med <- apply(tau_all, 2L, stats::median)

  # Per-arm community structure to match the Laplace fit's ms_community layout.
  off <- cumsum(c(0L, p_psi, p_pd))
  coef_by_arm <- c(list(Reduce(`+`, lapply(chains, `[[`, "b_psi")) / n.chains),
                   lapply(seq_len(D), function(d)
                     Reduce(`+`, lapply(chains, function(c) c$b_p[[d]])) / n.chains))
  Sigma_list <- list(); sd_list <- list(); coef_list <- list(); blup_list <- list()
  for (k in seq_along(arm_names)) {
    arm <- arm_names[k]; cn <- proc[[k]]$coef_names
    idx <- off[k] + seq_along(cn)
    coef <- coef_by_arm[[k]]
    rownames(coef) <- model$species_names; colnames(coef) <- cn
    blup <- sweep(coef, 2L, means[idx], "-")
    sdk  <- tau_med[idx]; names(sdk) <- cn
    Sig  <- diag(sdk^2, length(cn)); dimnames(Sig) <- list(cn, cn)
    Sigma_list[[paste0("Sigma_", arm)]] <- Sig
    sd_list[[paste0("sd_", arm)]]       <- sdk
    coef_list[[paste0("coef_", arm)]]   <- coef
    blup_list[[paste0("blup_", arm)]]   <- blup
  }
  ms_community <- c(Sigma_list, sd_list, coef_list, blup_list)

  re <- .tobs_nuts_rhat_ess(mu_chains)
  rhat <- re$rhat; ess <- re$ess; names(rhat) <- names(ess) <- par_names

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = V,
    n_samples = nrow(draws), n_params = length(means),
    log_prob = rep(NA_real_, nrow(draws)), log_lik = NA_real_,
    N = sum(vapply(model$valid, sum, integer(1))), rhat = rhat, ess = ess),
    list(
    col_names = par_names, param_names = par_names,
    n_fixed = length(means), fixed_names = par_names,
    process_info = proc, model = model, spatial = NULL,
    method = "pg_gibbs", n_chains = n.chains,
    ms_community = ms_community,
    convergence = list(converged = any(is.finite(rhat)) &&
                         max(rhat, na.rm = TRUE) < 1.1, n_iter = n.iter)
  )), class = c("tobs_fit", "tulpa_fit"))
}
