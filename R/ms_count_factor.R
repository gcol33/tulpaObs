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
.ms_count_factor_update <- function(offset_mat, y, eta, lambda, inner = 15L,
                                    center = FALSE) {
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
    # When composed with a shared areal field, centre the loadings across species
    # so the field owns the shared spatial mean and the factors own the between-
    # species differences (otherwise a near-constant factor trades off with the
    # field). No-op for the factor-only model.
    if (isTRUE(center)) lambda <- sweep(lambda, 2L, colMeans(lambda), "-")
  }
  list(eta = eta, lambda = lambda)
}


# Fit the community latent-factor count model -> the unified latent fitter
# (R/ms_count_spatial.R) with no shared field. Kept as a named entry for the
# dispatcher; the block coordinate ascent + assembly live in one place.
.tobs_fit_ms_count_factor <- function(model, latent, ...) {
  .tobs_fit_ms_count_latent(model, spatial = NULL, latent = latent, ...)
}
