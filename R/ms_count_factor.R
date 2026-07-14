# =============================================================================
# ms_count_factor.R - community latent-factor relative-abundance GLMM
# (ms_count() + latent(); the spAbundance lfMsAbund analogue, Poisson;
# gcol33/tulpaObs#117). Residual species co-occurrence via Q per-site latent
# factors with per-species loadings:
#
#   log mu_{s,i} = X_i . (mu + b_s) + sum_q lambda_{s,q} eta_{q,i}
#   eta_{q,i} ~ N(0, 1),  lambda_{s,q} = per-species factor loadings
#
# Fit by block coordinate ascent, reusing the pure-R community Laplace-EM
# unchanged: the factor term sum_q lambda_{s,q} eta_{q,i} enters each species'
# log-likelihood as a fixed per-observation offset (captured in the sp_ll
# closure), so a coefficient update is an ordinary community EM; the factor
# update given the coefficients is a Poisson factor analysis on the residuals
# (alternate eta | lambda and lambda | eta, with a unit-variance anchor on the
# factor columns). The loadings / factors are identified only up to rotation, so
# the recoverable target is the residual species covariance Sigma_res = lambda
# lambda' (reported alongside the implied correlation).
# =============================================================================


# Poisson factor update: fit the per-site factors eta [n_sites x Q] and per-
# species loadings lambda [n_species x Q] to the residuals over the fixed
# coefficient offsets `offset_mat` [n_sites x n_species]. Alternating Newton
# (eta_i | lambda; lambda_s | eta) with an N(0, I) factor prior and a weak
# loading ridge; the factor columns are centred + scaled to unit variance each
# pass (the scale folds into lambda), a standard factor-analysis anchor.
.ms_count_factor_update <- function(offset_mat, y, eta, lambda, inner = 15L) {
  Ns <- nrow(y); S <- ncol(y); Qk <- ncol(eta)
  for (it in seq_len(inner)) {
    for (i in seq_len(Ns)) {
      mu <- exp(pmin(offset_mat[i, ] + as.numeric(lambda %*% eta[i, ]), 700))
      g  <- as.numeric(crossprod(lambda, y[i, ] - mu)) - eta[i, ]
      H  <- crossprod(lambda * sqrt(mu)) + diag(Qk)
      eta[i, ] <- eta[i, ] + solve(H, g)
    }
    for (s in seq_len(S)) {
      mu <- exp(pmin(offset_mat[, s] + as.numeric(eta %*% lambda[s, ]), 700))
      g  <- as.numeric(crossprod(eta, y[, s] - mu)) - 1e-3 * lambda[s, ]
      H  <- crossprod(eta * sqrt(mu)) + 1e-3 * diag(Qk)
      lambda[s, ] <- lambda[s, ] + solve(H, g)
    }
    for (q in seq_len(Qk)) {
      eta[, q] <- eta[, q] - mean(eta[, q])
      sdq <- stats::sd(eta[, q])
      if (sdq > 1e-6) { eta[, q] <- eta[, q] / sdq; lambda[, q] <- lambda[, q] * sdq }
    }
  }
  list(eta = eta, lambda = lambda)
}


# Fit the community latent-factor count model. `model` is the ms_count model;
# `latent` the resolved tobs_latent term (n_factors). Block coordinate ascent
# between the community EM (factor offset) and the Poisson factor update.
.tobs_fit_ms_count_factor <- function(model, latent,
                                      max.iter = 200L, tol = 1e-4,
                                      sigma.beta = 5, priors = NULL,
                                      max.outer = 25L, verbose = FALSE, ...) {
  if (!identical(model$response %||% "poisson", "poisson")) {
    stop("Community latent-factor count (ms_count + latent()) is Poisson-only ",
         "in this release (gcol33/tulpaObs#117).", call. = FALSE)
  }
  Qk <- as.integer(latent$n_factors %||% 1L)
  if (Qk < 1L) stop("latent(): n_factors must be >= 1.", call. = FALSE)
  X  <- model$X; P <- ncol(X); S <- model$n_species; Ns <- model$n_sites
  su <- model$summaries
  if (any(!model$valid)) {
    stop("ms_count() latent factors need a complete y (no NA species-site ",
         "cells) (gcol33/tulpaObs#117).", call. = FALSE)
  }
  if (Qk > S - 1L) {
    stop(sprintf("latent(): n_factors (%d) must be < n_species (%d).", Qk, S),
         call. = FALSE)
  }
  y_mat <- matrix(as.numeric(model$y), Ns, S)

  arm_idx <- list(mu = seq_len(P))
  # Deterministic init (Date/Random unavailable in some contexts; a fixed small
  # grid pattern is enough to break the factor symmetry).
  eta    <- matrix(0, Ns, Qk)
  for (q in seq_len(Qk)) eta[, q] <- scale(cos(seq_len(Ns) * q))[, 1]
  lambda <- matrix(0.1, S, Qk)
  mu0 <- numeric(P); mu0[1L] <- log(max(mean(rowMeans(y_mat)), 0.1))
  em  <- NULL

  for (outer in seq_len(max.outer)) {
    foff <- tcrossprod(eta, lambda)                    # Ns x S factor offset
    sp_ll <- function(s, theta, global) {
      eta_s <- as.numeric(su[[s]]$X %*% theta) + foff[, s]
      sum(stats::dpois(su[[s]]$y, exp(pmin(eta_s, 700)), log = TRUE))
    }
    sp_grad <- function(s, theta, global) {
      mu_s <- exp(pmin(as.numeric(su[[s]]$X %*% theta) + foff[, s], 700))
      as.numeric(crossprod(su[[s]]$X, su[[s]]$y - mu_s))
    }
    em <- .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em)) mu0 else em$mu, init_global = numeric(0),
      penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
      sigma_init = 0.3, max_iter = min(as.integer(max.iter), 50L),
      tol = as.numeric(tol), newton_max = 30L, verbose = FALSE)

    offset_mat <- vapply(seq_len(S),
                         function(s) as.numeric(X %*% (em$mu + em$b_list[[s]])),
                         numeric(Ns))
    fu <- .ms_count_factor_update(offset_mat, y_mat, eta, lambda)
    delta <- max(abs(tcrossprod(fu$eta, fu$lambda) - foff))
    eta <- fu$eta; lambda <- fu$lambda
    if (isTRUE(verbose)) {
      message(sprintf("[ms_count factor %d] offset delta=%.2e", outer, delta))
    }
    if (outer > 2L && delta < tol) break
  }

  fit <- build_ms_count_fit(model, em, arm_idx, disp = NULL)
  Sigma_res <- tcrossprod(lambda)
  dimnames(Sigma_res) <- list(model$species_names, model$species_names)
  cor_res <- stats::cov2cor(Sigma_res + diag(1e-10, S))
  rownames(lambda) <- model$species_names
  colnames(lambda) <- paste0("factor", seq_len(Qk))

  fit$method  <- "laplace"
  fit$latent  <- latent
  fit$ms_factor <- list(
    n_factors = Qk, loadings = lambda, factors = eta,
    residual_cov = Sigma_res, residual_cor = cor_res)
  # fitted() / WAIC add the per-(species, site) factor offset.
  fit$model$count_factor_offset <- tcrossprod(eta, lambda)
  fit
}
