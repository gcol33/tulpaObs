# ms_dyn_occu_pg_gibbs.R - Polya-Gamma Gibbs for COMMUNITY multi-season (dynamic)
# occupancy (spOccupancy tMsPGOcc; gcol33/tulpaObs#115, #126). Combines the
# community PG machinery (msPGOcc, R/ms_occu_pg_gibbs.R) with a 2-state HMM
# forward-filter backward-sample (FFBS) latent-state step:
#
#   z_{s,i,1}          ~ Bernoulli(psi1_{s,i})      per-species season-1 occupancy
#   z_{s,i,t}          : colonization gamma_i (0->1), survival 1 - eps_i (1->1)
#   y_{s,i,t,j}|z=1    ~ Bernoulli(p_{s,i})          per-species detection
#
# The season-1 occupancy psi1 and detection p carry per-species Gaussian
# community hyperpriors (like msPGOcc); the colonization gamma and extinction eps
# are SHARED community-level coefficients (as ms_dyn_occu's Laplace-EM has them).
# Each Gibbs sweep: FFBS the occupancy path z per species, then the PG-augmented
# conjugate coefficient updates -- per-species beta_psi1 / beta_p, the SHARED
# beta_gamma / beta_eps from the aggregated 0->/1-> transitions, and the
# conjugate community mean + Inverse-Gamma community variance for the psi1 / p
# arms. PG draws via tulpa's Polson-Scott-Windle sampler. Gives a calibrated
# community-variance posterior (the community Laplace-EM attenuates it).
#
# v1: constant transitions, site-level detection, no structured terms.

.tobs_fit_ms_dyn_occu_pg_gibbs <- function(model, priors = NULL, sigma.beta = NULL,
                                           n.iter = NULL, n.warmup = NULL,
                                           n.chains = NULL, n.thin = NULL, seed = NULL,
                                           verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "pg_gibbs")

  rpg   <- get("cpp_rpg", envir = asNamespace("tulpa"))
  Xp1   <- model$X_psi1; Xp <- model$X_p; Xg <- model$X_gamma; Xe <- model$X_eps
  n     <- model$n_sites; T_s <- model$n_seasons; S <- model$n_species
  p_p1  <- ncol(Xp1); p_pd <- ncol(Xp); p_g <- ncol(Xg); p_e <- ncol(Xe)

  # Per-(species) [n x T] detection sufficient statistics.
  kmat <- nmat <- vector("list", S)
  for (s in seq_len(S)) {
    ys <- model$y[, , , s]; vs <- model$valid[, , , s]
    k <- nv <- matrix(0L, n, T_s)
    for (t in seq_len(T_s)) {
      yt <- ys[, , t]; vt <- vs[, , t]
      nv[, t] <- rowSums(vt); k[, t] <- rowSums(yt == 1L & vt)
    }
    kmat[[s]] <- k; nmat[[s]] <- nv
  }

  # 2-state FFBS for one species: returns the sampled z [n x T].
  ffbs <- function(psi1, p, gamma, eps, k, nv) {
    e1 <- p^k * (1 - p)^(nv - k)                    # occupied emission [n x T]
    e0 <- ifelse(k > 0L, 0, 1)                      # empty: 0 if a detection
    a1 <- a0 <- matrix(0, n, T_s)
    a1[, 1] <- psi1 * e1[, 1]; a0[, 1] <- (1 - psi1) * e0[, 1]
    ct <- a1[, 1] + a0[, 1]; ct[ct <= 0] <- 1
    a1[, 1] <- a1[, 1] / ct; a0[, 1] <- a0[, 1] / ct
    for (t in 2:T_s) {
      pr1 <- a1[, t - 1] * (1 - eps) + a0[, t - 1] * gamma
      pr0 <- a1[, t - 1] * eps       + a0[, t - 1] * (1 - gamma)
      a1[, t] <- pr1 * e1[, t]; a0[, t] <- pr0 * e0[, t]
      ct <- a1[, t] + a0[, t]; ct[ct <= 0] <- 1
      a1[, t] <- a1[, t] / ct; a0[, t] <- a0[, t] / ct
    }
    z <- matrix(0L, n, T_s)
    z[, T_s] <- stats::rbinom(n, 1L, a1[, T_s])
    for (t in (T_s - 1):1) {
      nxt1 <- z[, t + 1] == 1L
      w1 <- ifelse(nxt1, a1[, t] * (1 - eps), a1[, t] * eps)
      w0 <- ifelse(nxt1, a0[, t] * gamma,     a0[, t] * (1 - gamma))
      wt <- w1 + w0; wt[wt <= 0] <- 1
      z[, t] <- stats::rbinom(n, 1L, w1 / wt)
    }
    z
  }

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  par_names <- c(paste0("psi1_",  model$process_info[[1L]]$coef_names),
                 paste0("p_",     model$process_info[[2L]]$coef_names),
                 paste0("gamma_", model$process_info[[3L]]$coef_names),
                 paste0("eps_",   model$process_info[[4L]]$coef_names))

  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    b_p1 <- matrix(stats::rnorm(S * p_p1, 0, 0.2), S, p_p1)
    b_pd <- matrix(stats::rnorm(S * p_pd, 0, 0.2), S, p_pd)
    b_g  <- rep(0, p_g); b_e <- rep(0, p_e)
    mu_p1 <- rep(0, p_p1); tau2_p1 <- rep(1, p_p1)
    mu_pd <- rep(0, p_pd); tau2_pd <- rep(1, p_pd)
    out <- matrix(NA_real_, n_keep, p_p1 + p_pd + p_g + p_e)
    tau_out <- matrix(NA_real_, n_keep, p_p1 + p_pd)
    b_p1_sum <- matrix(0, S, p_p1); b_pd_sum <- matrix(0, S, p_pd); nsum <- 0L
    ki <- 0L
    for (it in seq_len(n.iter)) {
      gamma <- stats::plogis(as.vector(Xg %*% b_g))
      eps   <- stats::plogis(as.vector(Xe %*% b_e))
      # transition aggregation accumulators (shared arms)
      n0 <- k0 <- n1 <- k1 <- numeric(n)
      for (s in seq_len(S)) {
        psi1_s <- stats::plogis(as.vector(Xp1 %*% b_p1[s, ]))
        p_s    <- stats::plogis(as.vector(Xp  %*% b_pd[s, ]))
        z <- ffbs(psi1_s, p_s, gamma, eps, kmat[[s]], nmat[[s]])
        # season-1 occupancy -> beta_psi1_s
        om1 <- rpg(rep(1, n), as.vector(Xp1 %*% b_p1[s, ]))
        b_p1[s, ] <- .tobs_pg_draw_beta(Xp1, om1, z[, 1] - 0.5,
                                        1 / tau2_p1, mu_p1 / tau2_p1)
        # detection at occupied site-seasons -> beta_p_s (aggregated per site)
        occ <- z == 1L
        nocc <- rowSums(nmat[[s]] * occ); kocc <- rowSums(kmat[[s]] * occ)
        si <- which(nocc > 0)
        if (length(si) >= p_pd) {
          Xo <- Xp[si, , drop = FALSE]
          omp <- rpg(nocc[si], as.vector(Xo %*% b_pd[s, ]))
          b_pd[s, ] <- .tobs_pg_draw_beta(Xo, omp, kocc[si] - nocc[si] / 2,
                                          1 / tau2_pd, mu_pd / tau2_pd)
        }
        # transitions for the shared gamma / eps
        if (T_s >= 2L) {
          prev <- z[, 1:(T_s - 1), drop = FALSE]; cur <- z[, 2:T_s, drop = FALSE]
          from0 <- prev == 0L; from1 <- prev == 1L
          n0 <- n0 + rowSums(from0);            k0 <- k0 + rowSums(from0 & cur == 1L)
          n1 <- n1 + rowSums(from1);            k1 <- k1 + rowSums(from1 & cur == 0L)
        }
      }
      # shared colonization gamma (0 -> 1) and extinction eps (1 -> 0)
      s0 <- which(n0 > 0)
      if (length(s0) >= p_g) {
        Xg0 <- Xg[s0, , drop = FALSE]
        omg <- rpg(n0[s0], as.vector(Xg0 %*% b_g))
        b_g <- .tobs_pg_draw_beta(Xg0, omg, k0[s0] - n0[s0] / 2, 1 / sigma.beta^2)
      }
      s1 <- which(n1 > 0)
      if (length(s1) >= p_e) {
        Xe1 <- Xe[s1, , drop = FALSE]
        ome <- rpg(n1[s1], as.vector(Xe1 %*% b_e))
        b_e <- .tobs_pg_draw_beta(Xe1, ome, k1[s1] - n1[s1] / 2, 1 / sigma.beta^2)
      }
      # community mean + Inverse-Gamma variance for the per-species arms
      cu <- .tobs_pg_community_update(b_p1, mu_p1, tau2_p1, S)
      mu_p1 <- cu$mu; tau2_p1 <- cu$tau2
      cu <- .tobs_pg_community_update(b_pd, mu_pd, tau2_pd, S)
      mu_pd <- cu$mu; tau2_pd <- cu$tau2
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L
        out[ki, ]     <- c(mu_p1, mu_pd, b_g, b_e)
        tau_out[ki, ] <- sqrt(c(tau2_p1, tau2_pd))
        b_p1_sum <- b_p1_sum + b_p1; b_pd_sum <- b_pd_sum + b_pd; nsum <- nsum + 1L
      }
    }
    list(mu = out, tau = tau_out, b_p1 = b_p1_sum / nsum, b_pd = b_pd_sum / nsum)
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  summ <- .tobs_pg_summarize(lapply(chains, `[[`, "mu"), par_names)
  means <- summ$means

  tau_all <- do.call(rbind, lapply(chains, `[[`, "tau"))
  tau_med <- apply(tau_all, 2L, stats::median)
  sd_psi1 <- tau_med[seq_len(p_p1)]; sd_p <- tau_med[p_p1 + seq_len(p_pd)]
  names(sd_psi1) <- model$process_info[[1L]]$coef_names
  names(sd_p)    <- model$process_info[[2L]]$coef_names

  coef_psi1 <- Reduce(`+`, lapply(chains, `[[`, "b_p1")) / n.chains
  coef_p    <- Reduce(`+`, lapply(chains, `[[`, "b_pd")) / n.chains
  rownames(coef_psi1) <- rownames(coef_p) <- model$species_names
  colnames(coef_psi1) <- model$process_info[[1L]]$coef_names
  colnames(coef_p)    <- model$process_info[[2L]]$coef_names

  .tobs_pg_finalize_fit(
    summ, par_names, model, model$process_info, N = sum(model$valid),
    n.iter = n.iter, n.chains = n.chains,
    extra = list(ms_community = list(
      Sigma_psi1 = diag(sd_psi1^2, p_p1), Sigma_p = diag(sd_p^2, p_pd),
      sd_psi1 = sd_psi1, sd_p = sd_p, coef_psi1 = coef_psi1, coef_p = coef_p,
      blup_psi1 = sweep(coef_psi1, 2L, means[seq_len(p_p1)], "-"),
      blup_p    = sweep(coef_p,    2L, means[p_p1 + seq_len(p_pd)], "-"))))
}
