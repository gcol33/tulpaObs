# =============================================================================
# ms_count_nuts.R - community relative-abundance (count) NUTS (msAbund NUTS).
# Samples the exact joint posterior of the non-spatial community Poisson GLMM
# (community means mu, per-species coefficient deviations b_s, community
# covariance Sigma) via the in-tree C++ FullGradFn (src/ms_count_nuts.cpp),
# warm-started at the community Laplace-EM mode. NON-CENTERED: b_s = C z_s,
# z_s ~ N(0, I), so Sigma (log-Cholesky C) enters only the data term. Poisson.
# The reduced counterpart of ms_abun_nuts.R (no detection, no latent N).
# =============================================================================


# Packed layout: mu (P), z species-major (S*P), chol (P*(P+1)/2).
.tobs_ms_count_nuts_layout <- function(P, n_species) {
  q <- as.integer(P * (P + 1L) / 2L)
  list(P = P, n_species = n_species, q = q,
       mu = seq_len(P), b_off = P,
       chol = P + n_species * P + seq_len(q),
       total = P + n_species * P + q)
}

.tobs_ms_count_nuts_b_idx <- function(lay, s) {
  lay$b_off + (s - 1L) * lay$P + seq_len(lay$P)
}

.tobs_ms_count_nuts_priors <- function() {
  list(chol_logdiag_mean = log(0.5), chol_logdiag_sd = 1.5, chol_offdiag_sd = 1.0)
}

# Full-vector joint log-posterior + gradient (the NUTS target / oracle). `Y` is
# the n_sites x n_species count matrix; `X` the site-level design. Mirrors the
# C++ ms_count_nuts_eval (src/ms_count_nuts.cpp) exactly.
.tobs_ms_count_nuts_logpost <- function(theta, X, Y, lay, priors,
                                        sigma.beta = 10, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu <- theta[lay$mu]
  C  <- .ms_ocs_chol_unpack(theta[lay$chol], P)
  g  <- numeric(lay$total); g_mu <- numeric(P); A <- matrix(0, P, P); lp <- 0

  for (s in seq_len(S)) {
    bidx <- .tobs_ms_count_nuts_b_idx(lay, s)
    z_s  <- theta[bidx]
    b    <- mu + as.numeric(C %*% z_s)
    eta  <- as.numeric(X %*% b)
    lam  <- exp(pmin(eta, 700))
    lp   <- lp + sum(Y[, s] * eta - lam - lgamma(Y[, s] + 1))
    if (grad) {
      gl <- as.numeric(crossprod(X, Y[, s] - lam))     # grad_b (data)
      g_mu     <- g_mu + gl
      g[bidx]  <- g[bidx] + as.numeric(crossprod(C, gl))  # z grad = C' grad_b
      A <- A + outer(gl, z_s)
    }
  }
  # z prior N(0, I)
  z_idx <- lay$b_off + seq_len(S * P); z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all
  # chol coords: data gradient (b = C z) + log-Cholesky hyperprior
  pr <- .ms_ocs_chol_logprior(theta[lay$chol], P, priors)
  lp <- lp + pr$lp
  if (grad) g[lay$chol] <- .ms_abun_nuts_chol_data_grad(A, C, P) + pr$grad
  # community-mean prior N(0, sigma.beta^2)
  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * sum(mu^2)
  g_mu <- g_mu - ib2 * mu
  g[lay$mu] <- g[lay$mu] + g_mu
  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Per-species deviation matrix B (S x P) from a packed vector, b_s = C z_s.
.tobs_ms_count_nuts_b_from_z <- function(theta, lay) {
  C <- .ms_ocs_chol_unpack(theta[lay$chol], lay$P)
  B <- matrix(0, lay$n_species, lay$P)
  for (s in seq_len(lay$n_species))
    B[s, ] <- as.numeric(C %*% theta[.tobs_ms_count_nuts_b_idx(lay, s)])
  B
}

# Pack the community Laplace-EM warm (mu, Sigma, per-species B) into theta: the
# means, Sigma as log-Cholesky, and the whitened deviations z_s = C^{-1} b_s.
.tobs_ms_count_nuts_pack_init <- function(mu, Sigma, B, lay) {
  theta <- numeric(lay$total)
  theta[lay$mu] <- as.numeric(mu)
  Sig <- (Sigma + t(Sigma)) / 2
  if (inherits(try(chol(Sig), silent = TRUE), "try-error"))
    Sig <- Sig + diag(max(1e-6, 1e-6 * mean(diag(Sig))), ncol(Sig))
  C <- t(chol(Sig))
  theta[lay$chol] <- .ms_ocs_chol_pack(C)
  for (s in seq_len(lay$n_species))
    theta[.tobs_ms_count_nuts_b_idx(lay, s)] <- forwardsolve(C, B[s, ])
  theta
}


# Fit the community count model by NUTS. Warm-starts at the Poisson community
# Laplace-EM mode, then samples the exact joint posterior.
.tobs_fit_ms_count_nuts <- function(model, sigma.beta = 10,
                                    n.iter = 1000L, n.warmup = 1000L,
                                    n.chains = 1L, max.treedepth = 10L,
                                    adapt.delta = 0.9, seed = 1L,
                                    verbose = FALSE, ...) {
  if (!identical(model$response %||% "poisson", "poisson")) {
    stop("ms_count() NUTS is Poisson-only in this release; negbin / gaussian ",
         "community NUTS are follow-ups (gcol33/tulpaObs#117).", call. = FALSE)
  }
  # Warm start: the Poisson community Laplace-EM.
  warm <- .tobs_fit_ms_count(model, verbose = FALSE)
  cm   <- warm$ms_community
  P    <- model$process_info[[1L]]$p; S <- model$n_species
  mu   <- unname(as.numeric(warm$means[seq_len(P)]))
  Sigma <- as.matrix(cm$Sigma_mu); B <- as.matrix(cm$blup_mu)
  Y <- matrix(as.numeric(model$y), model$n_sites, S)
  X <- model$X

  lay    <- .tobs_ms_count_nuts_layout(P, S)
  priors <- .tobs_ms_count_nuts_priors()
  theta0 <- .tobs_ms_count_nuts_pack_init(mu, Sigma, B, lay)
  spec   <- list(X = X, y = Y)

  chains <- lapply(seq_len(n.chains), function(ch) {
    res <- cpp_ms_count_nuts(spec, theta0, priors, sigma.beta, NULL,
                             as.integer(n.iter), as.integer(n.warmup),
                             as.integer(max.treedepth), adapt.delta,
                             as.integer(seed + ch - 1L), isTRUE(verbose))
    res$draws
  })
  draws_all <- do.call(rbind, chains)               # (n_chains*n_sample) x total

  # Community means + covariance from the mu draws.
  mu_draws <- draws_all[, lay$mu, drop = FALSE]
  cn  <- model$process_info[[1L]]$coef_names
  nms <- paste0("mu_", cn)
  means <- colMeans(mu_draws); names(means) <- nms
  V <- stats::cov(mu_draws); dimnames(V) <- list(nms, nms)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nms
  colnames(mu_draws) <- nms

  # Posterior-mean per-species deviations + community covariance.
  nd <- nrow(draws_all)
  B_acc <- matrix(0, S, P); Sig_acc <- matrix(0, P, P)
  for (r in seq_len(nd)) {
    th <- draws_all[r, ]
    B_acc <- B_acc + .tobs_ms_count_nuts_b_from_z(th, lay)
    C <- .ms_ocs_chol_unpack(th[lay$chol], P); Sig_acc <- Sig_acc + tcrossprod(C)
  }
  blup <- B_acc / nd; coef <- sweep(blup, 2L, means, "+")
  rownames(blup) <- rownames(coef) <- model$species_names
  colnames(blup) <- colnames(coef) <- cn
  Sigma_mu <- Sig_acc / nd; dimnames(Sigma_mu) <- list(cn, cn)

  rhat_ess <- .tobs_nuts_rhat_ess(chains, n.chains)

  structure(c(list(
    draws = mu_draws, means = means, sds = sds, vcov = V,
    n_samples = nd, n_params = length(means),
    log_prob = rep(NA_real_, nd), N = sum(model$valid)),
    .tobs_na_nuts_diagnostics(nd),
    list(
    col_names = nms, param_names = nms, n_fixed = length(means),
    fixed_names = nms, process_info = model$process_info,
    model = model, spatial = NULL, method = "nuts",
    ms_community = list(Sigma_mu = Sigma_mu,
                        sd_mu = sqrt(pmax(diag(Sigma_mu), 0)),
                        coef_mu = coef, blup_mu = blup),
    ms_dispersion = NULL,
    rhat = rhat_ess$rhat, ess = rhat_ess$ess,
    convergence = list(converged = TRUE, n_iter = as.integer(n.iter))
  )), class = c("tobs_fit", "tulpa_fit"))
}
