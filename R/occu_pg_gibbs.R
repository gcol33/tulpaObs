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

.tobs_fit_occu_pg_gibbs <- function(model, priors = NULL, sigma.beta = 2.5,
                                    n.iter = 2000L, n.warmup = 1000L,
                                    n.chains = 2L, n.thin = 1L, seed = 1L,
                                    verbose = FALSE, ...) {
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

  B0inv_psi <- diag(1 / sigma.beta^2, p_psi)
  B0inv_p   <- diag(1 / sigma.beta^2, p_p)

  # One conjugate Gaussian PG update: beta ~ N(V X'kappa, V),
  # V = (X' diag(omega) X + B0inv)^-1.
  draw_beta <- function(X, omega, kappa, B0inv) {
    XtOX <- crossprod(X, X * omega) + B0inv
    V    <- chol2inv(chol(XtOX))
    m    <- V %*% crossprod(X, kappa)
    as.vector(m + t(chol(V)) %*% stats::rnorm(ncol(X)))
  }

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    bpsi <- stats::rnorm(p_psi, 0, 0.1); bp <- stats::rnorm(p_p, 0, 0.1)
    out <- matrix(NA_real_, n_keep, p_psi + p_p); ki <- 0L
    for (it in seq_len(n.iter)) {
      eta_psi <- as.vector(X_psi %*% bpsi); psi  <- stats::plogis(eta_psi)
      eta_p   <- as.vector(X_p   %*% bp);   pdet <- stats::plogis(eta_p)
      # 1. latent occupancy
      z <- integer(n); z[anydet] <- 1L
      und <- !anydet
      if (any(und)) {
        l1 <- psi[und] * (1 - pdet[und])^nvis[und]
        l0 <- 1 - psi[und]
        z[und] <- stats::rbinom(sum(und), 1L, l1 / (l1 + l0))
      }
      # 2. occupancy coefficients (all sites)
      om_psi <- rpg(rep(1, n), eta_psi)
      bpsi   <- draw_beta(X_psi, om_psi, z - 0.5, B0inv_psi)
      # 3. detection coefficients (occupied sites with >= 1 visit)
      occ <- which(z == 1L & nvis > 0L)
      if (length(occ) >= p_p) {
        Xo     <- X_p[occ, , drop = FALSE]
        eta_po <- as.vector(Xo %*% bp)
        om_p   <- rpg(nvis[occ], eta_po)
        bp     <- draw_beta(Xo, om_p, kdet[occ] - nvis[occ] / 2, B0inv_p)
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
  for (c in seq_along(chains)) colnames(chains[[c]]) <- par_names

  draws <- do.call(rbind, chains)
  means <- colMeans(draws); names(means) <- par_names
  V <- stats::cov(draws); dimnames(V) <- list(par_names, par_names)
  sds <- apply(draws, 2L, stats::sd); names(sds) <- par_names

  # Split-Rhat + bulk-ESS across chains (shared NUTS diagnostic helper).
  re <- .tobs_nuts_rhat_ess(chains)
  rhat <- re$rhat; ess <- re$ess
  names(rhat) <- names(ess) <- par_names

  intercepts <- list(
    psi = stats::setNames(means[1L], par_names[1L]),
    p   = stats::setNames(means[p_psi + 1L], par_names[p_psi + 1L]))

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = nrow(draws),
    n_params     = length(means),
    log_prob     = rep(NA_real_, nrow(draws)),
    log_lik      = NA_real_,
    N            = n,
    rhat         = rhat,
    ess          = ess),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    intercepts   = intercepts,
    process_info = model$process_info,
    model        = model,
    spatial      = NULL,
    method       = "pg_gibbs",
    n_chains     = n.chains,
    convergence  = list(converged = any(is.finite(rhat)) &&
                          (max(rhat, na.rm = TRUE) < 1.1),
                        n_iter = n.iter)
  )), class = c("tobs_fit", "tulpa_fit"))
}
