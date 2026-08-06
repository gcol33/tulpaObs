# =============================================================================
# ms_count_nuts.R - community relative-abundance (count) NUTS (msAbund NUTS).
# Samples the exact joint posterior of the non-spatial community GLMM (community
# means mu, per-species coefficient deviations b_s, community covariance Sigma)
# via the in-tree C++ FullGradFn (src/ms_count_nuts.cpp), warm-started at the
# community Laplace-EM mode. NON-CENTERED: b_{s,arm} = C_arm z_{s,arm},
# z_s ~ N(0, I), so Sigma (log-Cholesky C) enters only the data term.
#
# Three response families share one target, mirroring the ms_count Laplace-EM:
#   Poisson   - one beta arm.
#   NegBin    - beta arm + a per-species dispersion RE log_r_s ~
#               N(mu_log_r, sigma_log_r^2), the second community arm (the same
#               structure as ms_abun / the ms_count negbin Laplace-EM).
#   Gaussian  - beta arm + S FREE per-species residual variances log_phi_s (no
#               community prior; matches the Laplace outer loop that estimates
#               each phi_s independently), each with a weakly-informative prior.
# The reduced counterpart of ms_abun_nuts.R (no detection, no latent N).
# =============================================================================


# Packed layout: mu (P), z species-major (S*P), chol_beta (q_beta),
#   [chol_logr (1) if negbin], [log_phi (S) if gaussian]. P = p_beta (+1 negbin).
# `lay$chol` aliases `chol_beta` (the Poisson single-arm cross-check reads it).
.tobs_ms_count_nuts_layout <- function(p_beta, n_species, family = "poisson") {
  is_nb    <- identical(family, "negbin")
  is_gauss <- identical(family, "gaussian")
  q_beta   <- as.integer(p_beta * (p_beta + 1L) / 2L)
  P        <- p_beta + (if (is_nb) 1L else 0L)
  b_off         <- P
  chol_beta_off <- P + n_species * P
  chol_logr_off <- chol_beta_off + q_beta
  logphi_off    <- chol_beta_off + q_beta
  q_logr <- if (is_nb) 1L else 0L
  n_phi  <- if (is_gauss) n_species else 0L
  chol_beta <- chol_beta_off + seq_len(q_beta)
  list(
    family = family, is_nb = is_nb, is_gauss = is_gauss,
    p_beta = p_beta, P = P, n_species = n_species, q_beta = q_beta,
    beta  = seq_len(p_beta),
    logr  = if (is_nb) p_beta + 1L else integer(0),
    mu    = seq_len(P), b_off = b_off,
    chol_beta = chol_beta,
    chol      = chol_beta,                          # backward-compat alias
    chol_logr = if (is_nb) chol_logr_off + seq_len(q_logr) else integer(0),
    logphi    = if (is_gauss) logphi_off + seq_len(n_phi) else integer(0),
    total = chol_beta_off + q_beta + q_logr + n_phi)
}

# Log-Cholesky hyperprior + weakly-informative gaussian log_phi prior scalars.
.tobs_ms_count_nuts_priors <- function() {
  list(chol_logdiag_mean = log(0.5), chol_logdiag_sd = 1.5, chol_offdiag_sd = 1.0,
       logphi_mean = 0, logphi_sd = 2)
}

# Full-vector joint log-posterior + gradient (the NUTS target / oracle). `Y` is
# the n_sites x n_species response matrix; `X` the site-level design. An NA entry
# Y[i, s] is a missing site x species observation: that (species, site) is dropped
# from the data term (matching the Laplace-EM per-species `valid` subsets), so a
# species keeps only its observed sites. Mirrors the C++ ms_count_nuts_eval
# (src/ms_count_nuts.cpp) exactly, including the NA skip.
.tobs_ms_count_nuts_logpost <- function(theta, X, Y, lay, priors,
                                        sigma.beta = 10, sigma.logr = 1.5,
                                        grad = TRUE) {
  pb <- lay$p_beta; S <- lay$n_species; P <- lay$P
  is_nb <- lay$is_nb; is_gauss <- lay$is_gauss
  is_bern <- identical(lay$family, "bernoulli")
  mu     <- theta[lay$mu]
  C_beta <- .ms_ocs_chol_unpack(theta[lay$chol_beta], pb)
  C_lr   <- if (is_nb) exp(theta[lay$chol_logr]) else NULL
  logphi <- if (is_gauss) theta[lay$logphi] else NULL

  g <- numeric(lay$total); g_mu <- numeric(P)
  A_beta <- matrix(0, pb, pb); A_lr <- 0; lp <- 0

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]; zb <- z_s[lay$beta]
    b_beta <- mu[lay$beta] + as.numeric(C_beta %*% zb)
    # Drop missing (NA) site x species observations for this species.
    ys_all <- Y[, s]
    obs    <- !is.na(ys_all)
    Xs     <- X[obs, , drop = FALSE]
    ys     <- ys_all[obs]
    eta    <- as.numeric(Xs %*% b_beta)
    if (is_nb) {
      zr  <- z_s[lay$logr]
      r   <- exp(min(mu[lay$logr] + C_lr * zr, 30))
      muv <- pmax(exp(pmin(eta, 700)), 1e-10)
      lp  <- lp + sum(stats::dnbinom(ys, size = r, mu = muv, log = TRUE))
      if (grad) {
        gl     <- as.numeric(crossprod(Xs, r * (ys - muv) / (r + muv)))
        dLL_dr <- digamma(ys + r) - digamma(r) + log(r / (r + muv)) + 1 -
                  (ys + r) / (r + muv)
        g_logr <- r * sum(dLL_dr)
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g_mu[lay$logr] <- g_mu[lay$logr] + g_logr
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        g[bidx[lay$logr]] <- g[bidx[lay$logr]] + C_lr * g_logr
        A_beta <- A_beta + outer(gl, zb)
        A_lr   <- A_lr + g_logr * zr
      }
    } else if (is_gauss) {
      phi_s <- exp(min(logphi[s], 30))
      lp    <- lp + sum(stats::dnorm(ys, mean = eta, sd = sqrt(phi_s), log = TRUE))
      if (grad) {
        gl <- as.numeric(crossprod(Xs, (ys - eta) / phi_s))
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        g[lay$logphi[s]]  <- g[lay$logphi[s]] + sum(-0.5 + (ys - eta)^2 / (2 * phi_s))
        A_beta <- A_beta + outer(gl, zb)
      }
    } else if (is_bern) {
      # jsdm(): log p = y eta - log(1 + exp(eta)); d log p / d eta = y - plogis(eta).
      # Written via plogis(log.p) to match the C++ stable log-sum-exp branch.
      lp <- lp + sum(ifelse(ys > 0, stats::plogis(eta,  log.p = TRUE),
                                    stats::plogis(-eta, log.p = TRUE)))
      if (grad) {
        gl <- as.numeric(crossprod(Xs, ys - stats::plogis(eta)))
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        A_beta <- A_beta + outer(gl, zb)
      }
    } else {
      lam <- exp(pmin(eta, 700))
      lp  <- lp + sum(ys * eta - lam - lgamma(ys + 1))
      if (grad) {
        gl <- as.numeric(crossprod(Xs, ys - lam))
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        A_beta <- A_beta + outer(gl, zb)
      }
    }
  }
  # z prior N(0, I)
  z_idx <- lay$b_off + seq_len(S * P); z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all
  # chol_beta coords: data gradient (b = C z) + log-Cholesky hyperprior
  pr <- .ms_ocs_chol_logprior(theta[lay$chol_beta], pb, priors)
  lp <- lp + pr$lp
  if (grad) g[lay$chol_beta] <- .ms_abun_nuts_chol_data_grad(A_beta, C_beta, pb) +
      pr$grad
  if (is_nb) {                                # scalar log-dispersion covariance
    pr_r <- .ms_ocs_chol_logprior(theta[lay$chol_logr], 1L, priors)
    lp   <- lp + pr_r$lp
    if (grad) g[lay$chol_logr] <-
      .ms_abun_nuts_chol_data_grad(matrix(A_lr, 1, 1), matrix(C_lr, 1, 1), 1L) +
      pr_r$grad
  }
  if (is_gauss) {                             # free log_phi prior N(mean, sd^2)
    lpm <- priors$logphi_mean; lps <- priors$logphi_sd
    lp  <- lp - 0.5 * sum(((logphi - lpm) / lps)^2)
    if (grad) g[lay$logphi] <- g[lay$logphi] - (logphi - lpm) / lps^2
  }
  # community-mean priors
  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * sum(mu[lay$beta]^2)
  g_mu[lay$beta] <- g_mu[lay$beta] - ib2 * mu[lay$beta]
  if (is_nb) {
    il2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * il2 * mu[lay$logr]^2
    g_mu[lay$logr] <- g_mu[lay$logr] - il2 * mu[lay$logr]
  }
  g[lay$mu] <- g[lay$mu] + g_mu
  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Per-species deviation matrix B (S x P) from a packed vector, b = C z per arm.
.tobs_ms_count_nuts_b_from_z <- function(theta, lay) {
  C_beta <- .ms_ocs_chol_unpack(theta[lay$chol_beta], lay$p_beta)
  C_lr   <- if (lay$is_nb) exp(theta[lay$chol_logr]) else NULL
  B <- matrix(0, lay$n_species, lay$P)
  for (s in seq_len(lay$n_species)) {
    z <- theta[.ms_ocs_b_idx(lay, s)]
    B[s, lay$beta] <- as.numeric(C_beta %*% z[lay$beta])
    if (lay$is_nb) B[s, lay$logr] <- C_lr * z[lay$logr]
  }
  B
}

# Pack a community Laplace-EM warm start into theta: the means, each arm's
# covariance as log-Cholesky, the whitened per-species deviations z = C^{-1} b,
# and (gaussian) the per-species log residual variances.
.tobs_ms_count_nuts_pack_init <- function(lay, mu_beta, Sigma_beta, B_beta,
                                          mu_logr = NA_real_, sigma_logr = NA_real_,
                                          b_logr = NULL, logphi = NULL) {
  theta <- numeric(lay$total)
  P <- lay$P; S <- lay$n_species
  mu <- as.numeric(mu_beta)
  if (lay$is_nb) mu <- c(mu, mu_logr)
  theta[lay$mu] <- mu

  Sig <- (Sigma_beta + t(Sigma_beta)) / 2
  if (inherits(try(chol(Sig), silent = TRUE), "try-error"))
    Sig <- Sig + diag(max(1e-6, 1e-6 * mean(diag(Sig))), ncol(Sig))
  C_beta <- t(chol(Sig))
  theta[lay$chol_beta] <- .ms_ocs_chol_pack(C_beta)
  if (lay$is_nb)    theta[lay$chol_logr] <- log(sigma_logr)   # 1x1 log-chol = log(sd)
  if (lay$is_gauss) theta[lay$logphi]    <- as.numeric(logphi)

  for (s in seq_len(S)) {
    z_s <- numeric(P)
    z_s[lay$beta] <- forwardsolve(C_beta, B_beta[s, ])   # C z = b, C lower-tri
    if (lay$is_nb) z_s[lay$logr] <- b_logr[s] / sigma_logr
    theta[.ms_ocs_b_idx(lay, s)] <- z_s
  }
  theta
}


# Fit the community count model by NUTS. Warm-starts at the community Laplace-EM
# mode, then samples the exact joint posterior. Poisson / negbin / gaussian.
.tobs_fit_ms_count_nuts <- function(model, sigma.beta = 10, sigma.logr = 1.5,
                                    n.iter = 1000L, n.warmup = 1000L,
                                    n.chains = 1L, max.treedepth = 10L,
                                    adapt.delta = 0.9, seed = 1L,
                                    verbose = FALSE, ...) {
  response <- model$response %||% "poisson"
  is_nb    <- identical(response, "negbin")
  is_gauss <- identical(response, "gaussian")

  # Warm start: the community Laplace-EM (raw output, so the negbin log_r arm is
  # available to pack).
  em     <- .tobs_ms_count_run_em(model, verbose = FALSE)
  fit    <- em$fit
  P_beta <- em$P_beta; S <- model$n_species
  Y <- matrix(as.numeric(model$y), model$n_sites, S)
  X <- model$X

  lay    <- .tobs_ms_count_nuts_layout(P_beta, S, response)
  priors <- .tobs_ms_count_nuts_priors()
  if (is_gauss) {
    # Weakly-informative log_phi prior, centred at the pooled log-variance so it
    # is on the data's scale but wide enough to be swamped by the residuals.
    yv <- model$y[model$valid]
    priors$logphi_mean <- log(max(stats::var(yv), 1e-4))
    priors$logphi_sd   <- 2
  }

  # Warm-start pieces from the raw EM fit.
  mu_beta    <- fit$mu[seq_len(P_beta)]
  Sigma_beta <- as.matrix(fit$Sigma$mu)
  B_beta     <- do.call(rbind, lapply(fit$b_list, function(b) b[seq_len(P_beta)]))
  mu_logr    <- if (is_nb) fit$mu[P_beta + 1L] else NA_real_
  sigma_logr_w <- if (is_nb) max(sqrt(fit$Sigma$disp[1, 1]), 1e-3) else NA_real_
  b_logr     <- if (is_nb) vapply(fit$b_list, function(b) b[P_beta + 1L],
                                  numeric(1)) else NULL
  logphi_w   <- if (is_gauss) log(pmax(em$disp$variance, 1e-8)) else NULL

  theta0 <- .tobs_ms_count_nuts_pack_init(
    lay, mu_beta = mu_beta, Sigma_beta = Sigma_beta, B_beta = B_beta,
    mu_logr = mu_logr, sigma_logr = sigma_logr_w, b_logr = b_logr,
    logphi = logphi_w)

  spec <- list(X = X, y = Y, family = response)
  if (is_gauss) {
    spec$logphi_mean <- priors$logphi_mean
    spec$logphi_sd   <- priors$logphi_sd
  }

  run_chain <- function(ch) {
    # tulpa's NUTS engine takes n_iter as the TOTAL (warmup + sampling) count and
    # returns n_iter - n_warmup post-warmup draws.
    cpp_ms_count_nuts(spec, theta0, priors, sigma.beta, sigma.logr, NULL,
                      as.integer(n.iter + n.warmup), as.integer(n.warmup),
                      as.integer(max.treedepth), adapt.delta,
                      as.integer(seed + ch - 1L), isTRUE(verbose))
  }
  rc <- .ms_ocs_run_chains(run_chain, n.chains)
  chains    <- rc$chains
  draws_all <- rc$draws                             # (n_chains*n_sample) x total
  accept    <- rc$accept
  divergent <- rc$divergent
  treedepth <- rc$treedepth
  epsilon   <- rc$epsilon

  # Community means + covariance from the beta-coefficient mu draws.
  mu_draws <- draws_all[, lay$beta, drop = FALSE]
  cn  <- model$process_info[[1L]]$coef_names
  nms <- paste0("mu_", cn)
  means <- colMeans(mu_draws); names(means) <- nms
  V <- stats::cov(mu_draws); dimnames(V) <- list(nms, nms)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nms
  colnames(mu_draws) <- nms

  # Posterior-mean per-species deviations, community covariance, dispersion.
  nd <- nrow(draws_all)
  B_acc <- matrix(0, S, P_beta); Sig_acc <- matrix(0, P_beta, P_beta)
  logr_acc <- numeric(S); sigma_logr_acc <- 0; phi_acc <- numeric(S)
  for (r in seq_len(nd)) {
    th <- draws_all[r, ]
    Bd <- .tobs_ms_count_nuts_b_from_z(th, lay)
    B_acc <- B_acc + Bd[, lay$beta, drop = FALSE]
    C <- .ms_ocs_chol_unpack(th[lay$chol_beta], P_beta); Sig_acc <- Sig_acc + tcrossprod(C)
    if (is_nb) {
      logr_acc <- logr_acc + Bd[, lay$logr]
      sigma_logr_acc <- sigma_logr_acc + exp(th[lay$chol_logr])
    }
    if (is_gauss) phi_acc <- phi_acc + exp(th[lay$logphi])
  }
  blup <- B_acc / nd; coef <- sweep(blup, 2L, means, "+")
  rownames(blup) <- rownames(coef) <- model$species_names
  colnames(blup) <- colnames(coef) <- cn
  Sigma_mu <- Sig_acc / nd; dimnames(Sigma_mu) <- list(cn, cn)

  disp <- NULL
  if (is_nb) {
    mu_log_r_hat <- mean(draws_all[, lay$logr])
    disp <- list(response = "negbin", mu_log_r = mu_log_r_hat,
                 sigma_log_r = sigma_logr_acc / nd,
                 r_s = exp(mu_log_r_hat + logr_acc / nd))
  } else if (is_gauss) {
    disp <- list(response = "gaussian", variance = phi_acc / nd)
  }

  fit <- structure(c(list(
    draws = mu_draws, means = means, sds = sds, vcov = V,
    n_samples = nd, n_params = length(means),
    log_prob = rep(NA_real_, nd), N = sum(model$valid),
    accept_prob = accept, divergent = divergent, treedepth = treedepth,
    epsilon = epsilon),
    list(
    col_names = nms, param_names = nms, n_fixed = length(means),
    fixed_names = nms, process_info = model$process_info,
    model = model, spatial = NULL, method = "nuts",
    ms_community = list(Sigma_mu = Sigma_mu,
                        sd_mu = sqrt(pmax(diag(Sigma_mu), 0)),
                        coef_mu = coef, blup_mu = blup),
    ms_dispersion = disp,
    convergence = list(converged = NA, n_iter = as.integer(n.iter))
  )), class = c("tobs_fit", "tulpa_fit"))
  # The reported coefficients are the community beta means; the sampler also
  # carries the per-species z blocks, the packed community Cholesky and (negbin)
  # the community mean log_r, which no summary row names.
  fit <- .tobs_nuts_attach_convergence(fit, chains, par_names = nms,
                                       cols = lay$beta,
                                       n_iter = as.integer(n.iter))
  # Top-level aliases for the per-parameter record.
  fit$rhat <- fit$convergence$rhat
  fit$ess  <- fit$convergence$ess_bulk
  fit
}
