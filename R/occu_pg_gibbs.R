# occu_pg_gibbs.R - Polya-Gamma Gibbs sampler for single-season occupancy
# (spOccupancy PGOcc; gcol33/tulpaObs#126). This is a REAL Gibbs chain over the
# exact posterior -- distinct from `method = "laplace_gibbs"`, which is a
# stochastic-EM variance correction (Rubin-pooled mode-finds, no PG augmentation,
# stationary distribution not the posterior). Conditional on the Polya-Gamma
# auxiliary variables omega, both the occupancy and detection logistic
# coefficient updates are exactly conjugate Gaussian (Polson, Scott & Windle
# 2013), so each sweep is:
#
#   1. z_i | .      : occupied at any detection; else Bernoulli with weight
#                     psi_i (1 - p_i)^{n_i} / (psi_i (1 - p_i)^{n_i} + (1 - psi_i))
#   2. beta_psi | . : omega_psi_i ~ PG(1, eta_psi_i); kappa = z - 1/2;
#                     beta_psi ~ N(V X'kappa, V), V = (X' Omega X + B0^-1)^-1
#   3. beta_p | .   : at occupied sites, omega_p_i ~ PG(n_i, eta_p_i);
#                     kappa = k_i - n_i/2; beta_p ~ N(V X'kappa, V) over those sites
#
# with weakly-informative N(0, sigma.beta^2) coefficient priors. The PG draws use
# tulpa's tested Polson-Scott-Windle sampler (`tulpa:::cpp_rpg`). v1: single
# season, site-level detection, no random effects / spatial field (those are the
# PG-spatial extensions -- pg_binomial_{icar,bym2,...} exist in tulpa and are the
# documented follow-up).

.tobs_fit_occu_pg_gibbs <- function(model, priors = NULL, sigma.beta = NULL,
                                    n.iter = NULL, n.warmup = NULL,
                                    n.chains = NULL, n.thin = NULL, seed = NULL,
                                    verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "pg_gibbs", single_species = TRUE)

  if (!is.null(model$X_det_visit))
    stop("occu() method = \"pg_gibbs\" supports site-level detection only in v1 ",
         "(visit-level detection covariates are a follow-up).", call. = FALSE)
  rpg <- get("cpp_rpg", envir = asNamespace("tulpa"))

  y      <- model$y                              # [n x mv], -1 = NA
  X_psi  <- model$X_processes[[1L]]
  X_p    <- model$X_processes[[2L]]
  n      <- nrow(y); p_psi <- ncol(X_psi); p_p <- ncol(X_p)
  valid  <- y >= 0L
  nvis   <- rowSums(valid)
  kdet   <- rowSums(y == 1L & valid)
  anydet <- kdet > 0L

  prec <- 1 / sigma.beta^2          # weakly-informative N(0, sigma.beta^2) prior

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    bpsi <- stats::rnorm(p_psi, 0, 0.1); bp <- stats::rnorm(p_p, 0, 0.1)
    out <- matrix(NA_real_, n_keep, p_psi + p_p); ki <- 0L
    for (it in seq_len(n.iter)) {
      eta_psi <- as.vector(X_psi %*% bpsi); psi  <- stats::plogis(eta_psi)
      eta_p   <- as.vector(X_p   %*% bp);   pdet <- stats::plogis(eta_p)
      # 1. latent occupancy
      z <- .tobs_pg_draw_z(psi, (1 - pdet)^nvis, anydet)
      # 2. occupancy coefficients (all sites)
      om_psi <- rpg(rep(1, n), eta_psi)
      bpsi   <- .tobs_pg_draw_beta(X_psi, om_psi, z - 0.5, prec)
      # 3. detection coefficients (occupied sites with >= 1 visit)
      occ <- which(z == 1L & nvis > 0L)
      if (length(occ) >= p_p) {
        Xo     <- X_p[occ, , drop = FALSE]
        eta_po <- as.vector(Xo %*% bp)
        om_p   <- rpg(nvis[occ], eta_po)
        bp     <- .tobs_pg_draw_beta(Xo, om_p, kdet[occ] - nvis[occ] / 2, prec)
      }
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L; out[ki, ] <- c(bpsi, bp)
      }
    }
    out
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  par_names <- c(paste0("psi_", model$process_info[[1L]]$coef_names),
                 paste0("p_",   model$process_info[[2L]]$coef_names))
  summ <- .tobs_pg_summarize(chains, par_names)
  means <- summ$means

  intercepts <- list(
    psi = stats::setNames(means[1L], par_names[1L]),
    p   = stats::setNames(means[p_psi + 1L], par_names[p_psi + 1L]))

  .tobs_pg_finalize_fit(
    summ, par_names, model, model$process_info, N = n,
    n.iter = n.iter, n.chains = n.chains,
    extra = list(intercepts = intercepts))
}


# Spatial PG Gibbs for single-season occupancy with an intrinsic areal (ICAR)
# field on the occupancy logit (spOccupancy spPGOcc; gcol33/tulpaObs#126). The
# field f (one node per site) enters psi linearly: logit psi_i = X_i beta + f_i,
# f ~ ICAR(tau). Conditional on the Polya-Gamma auxiliaries the joint (beta, f)
# update is a Gaussian Markov random field draw -- the coefficient prior on beta
# and the intrinsic tau Q prior on f -- and tau has a conjugate Gamma full
# conditional. The field is centred (sum-to-zero) each sweep for identifiability
# against the intercept. icar only (bym2 adds the iid block; a follow-up).
.tobs_fit_occu_pg_gibbs_spatial <- function(model, spatial, priors = NULL,
                                            sigma.beta = NULL,
                                            n.iter = NULL, n.warmup = NULL,
                                            n.chains = NULL, n.thin = NULL, seed = NULL,
                                            verbose = FALSE) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "pg_gibbs", single_species = TRUE)

  if (!identical(spatial$type, "icar"))
    stop("occu() method = \"pg_gibbs\" + a spatial field supports icar() only in ",
         "v1 (bym2 / car_proper are follow-ups).", call. = FALSE)
  if (!is.null(model$X_det_visit))
    stop("occu() pg_gibbs supports site-level detection only.", call. = FALSE)
  rpg <- get("cpp_rpg", envir = asNamespace("tulpa"))
  adj <- as.matrix(spatial$graph)
  n   <- model$n_sites
  if (nrow(adj) != n)
    stop(sprintf(paste0("the icar() graph has %d nodes but the model has %d ",
                        "sites; one field node per site is required."),
                 nrow(adj), n), call. = FALSE)
  Q   <- .occu_cover_icar_Q(adj)

  y      <- model$y
  X_psi  <- model$X_processes[[1L]]; X_p <- model$X_processes[[2L]]
  p_psi  <- ncol(X_psi); p_p <- ncol(X_p)
  valid  <- y >= 0L
  nvis   <- rowSums(valid); kdet <- rowSums(y == 1L & valid); anydet <- kdet > 0L
  B0inv_psi <- diag(1 / sigma.beta^2, p_psi)   # joins the joint (beta_psi, f) GMRF
  prec      <- 1 / sigma.beta^2                 # detection-arm coefficient prior

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  par_names <- c(paste0("psi_", model$process_info[[1L]]$coef_names),
                 paste0("p_",   model$process_info[[2L]]$coef_names), "log_tau")

  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    bpsi <- stats::rnorm(p_psi, 0, 0.1); bp <- stats::rnorm(p_p, 0, 0.1)
    f <- rep(0, n); tau <- 1
    out <- matrix(NA_real_, n_keep, p_psi + p_p + 1L); ki <- 0L
    f_sum <- numeric(n); nsum <- 0L
    for (it in seq_len(n.iter)) {
      eta_psi <- as.vector(X_psi %*% bpsi) + f; psi <- stats::plogis(eta_psi)
      eta_p   <- as.vector(X_p %*% bp);        pdet <- stats::plogis(eta_p)
      # 1. latent occupancy
      z <- .tobs_pg_draw_z(psi, (1 - pdet)^nvis, anydet)
      # 2. joint (beta_psi, f) GMRF update given omega_psi
      om  <- rpg(rep(1, n), eta_psi); kap <- z - 0.5
      Pbb <- crossprod(X_psi, X_psi * om) + B0inv_psi          # p x p
      Pbf <- t(X_psi * om)                                     # p x n
      Pff <- diag(om, n) + tau * Q                             # n x n
      Prec <- rbind(cbind(Pbb, Pbf), cbind(t(Pbf), Pff))
      rhs  <- c(crossprod(X_psi, kap), kap)
      L    <- chol(Prec + diag(1e-8, p_psi + n))
      mean <- backsolve(L, forwardsolve(t(L), rhs))
      samp <- mean + backsolve(L, stats::rnorm(p_psi + n))
      bpsi <- samp[seq_len(p_psi)]; f <- samp[p_psi + seq_len(n)]
      # Sum-to-zero the field, moving its level into the intercept so eta = X beta
      # + f is preserved (the field mean is confounded with the intercept; column
      # 1 of X_psi is the all-ones intercept).
      mf <- mean(f); f <- f - mf; bpsi[1L] <- bpsi[1L] + mf
      # 3. field precision tau (ICAR rank n - 1)
      tau <- stats::rgamma(1, 0.5 + (n - 1) / 2, 0.5 + 0.5 * as.numeric(t(f) %*% Q %*% f))
      # 4. detection coefficients at occupied sites
      occ <- which(z == 1L & nvis > 0L)
      if (length(occ) >= p_p) {
        Xo <- X_p[occ, , drop = FALSE]
        om_p <- rpg(nvis[occ], as.vector(Xo %*% bp))
        bp <- .tobs_pg_draw_beta(Xo, om_p, kdet[occ] - nvis[occ] / 2, prec)
      }
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L; out[ki, ] <- c(bpsi, bp, log(tau))
        f_sum <- f_sum + f; nsum <- nsum + 1L
      }
    }
    list(draws = out, field = f_sum / nsum)
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  summ <- .tobs_pg_summarize(lapply(chains, `[[`, "draws"), par_names)
  spatial_field <- Reduce(`+`, lapply(chains, `[[`, "field")) / n.chains

  .tobs_pg_finalize_fit(
    summ, par_names, model, model$process_info, N = n,
    n.iter = n.iter, n.chains = n.chains, spatial = spatial,
    extra = list(spatial_field = spatial_field))
}
