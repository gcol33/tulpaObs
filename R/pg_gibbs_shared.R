# pg_gibbs_shared.R - single source of truth for the Polya-Gamma Gibbs conjugate
# machinery shared by every PG-Gibbs fitter (occu_pg_gibbs, t_occu,
# ms_occu_pg_gibbs, ms_int_occu_pg_gibbs, ms_dyn_occu_pg_gibbs,
# ms_count_pg_gibbs; gcol33/tulpaObs#135). Conditional on the Polya-Gamma
# auxiliaries omega, a logistic-coefficient update is exactly conjugate Gaussian
# and every fitter re-derived the same three pieces: the conjugate draw, the
# community mean + Inverse-Gamma variance update, and the post-sampling
# summary / fit assembly. Extracting them keeps the sampler math in one place;
# the per-family fitters keep only their latent-state step and design layout.

# Community hyperprior defaults: the community-mean prior variance and the
# near-Jeffreys Inverse-Gamma(a, b) shape/rate on each community variance
# component. Shared by every community PG-Gibbs fitter.
.tobs_pg_hyper <- list(sigma.mu2 = 100, ig_a = 0.1, ig_b = 0.1)

# One conjugate Gaussian Polya-Gamma coefficient update:
#   beta ~ N(V (X'kappa + lin), V),  V = (X' diag(omega) X + diag(prec_diag))^-1
# `prec_diag` = the diagonal prior precision, either a scalar (recycled) or a
# length-ncol(X) vector -- 1 / tau^2 for a community-mean prior, or
# 1 / sigma.beta^2 for a fixed weakly-informative prior. `lin` = the
# prior-precision-weighted prior mean diag(prec_diag) %*% prior_mean, supplied
# only when the prior mean is non-zero (community arms: mu / tau^2); NULL leaves
# it out. One rnorm(ncol(X)) call, so the caller's RNG stream is unchanged when
# an inlined draw is routed through here.
.tobs_pg_draw_beta <- function(X, omega, kappa, prec_diag, lin = NULL) {
  XtOX <- crossprod(X, X * omega)
  diag(XtOX) <- diag(XtOX) + prec_diag
  V   <- chol2inv(chol(XtOX))
  rhs <- crossprod(X, kappa)
  if (!is.null(lin)) rhs <- rhs + lin
  as.vector(V %*% rhs + t(chol(V)) %*% stats::rnorm(ncol(X)))
}

# Per-coordinate community update for one arm: the conjugate Gaussian community
# mean and the Inverse-Gamma community variance, given the S per-species
# coefficient rows `b` (S x p). Returns the updated `mu` (length p) and `tau2`
# (length p). One rnorm + one rgamma per coordinate, in coordinate order, so the
# RNG stream matches the inlined loop it replaces.
.tobs_pg_community_update <- function(b, mu, tau2, S, hyper = .tobs_pg_hyper) {
  sigma.mu2 <- hyper$sigma.mu2; ig_a <- hyper$ig_a; ig_b <- hyper$ig_b
  for (j in seq_len(ncol(b))) {
    vj      <- 1 / (S / tau2[j] + 1 / sigma.mu2)
    mu[j]   <- stats::rnorm(1, vj * sum(b[, j]) / tau2[j], sqrt(vj))
    tau2[j] <- 1 / stats::rgamma(1, ig_a + S / 2,
                                 ig_b + 0.5 * sum((b[, j] - mu[j])^2))
  }
  list(mu = mu, tau2 = tau2)
}

# Post-sampling summary shared by every PG-Gibbs fitter: pool the per-chain
# posterior-draw matrices, name by `par_names`, and compute the posterior mean /
# SD / covariance plus split-Rhat + bulk-ESS (the shared NUTS diagnostic). The
# chains carry a parameter per column in `par_names` order.
.tobs_pg_summarize <- function(chains, par_names) {
  draws <- do.call(rbind, chains); colnames(draws) <- par_names
  means <- colMeans(draws); names(means) <- par_names
  V   <- stats::cov(draws); dimnames(V) <- list(par_names, par_names)
  sds <- apply(draws, 2L, stats::sd); names(sds) <- par_names
  re  <- .tobs_nuts_rhat_ess(chains)
  rhat <- re$rhat; ess <- re$ess; names(rhat) <- names(ess) <- par_names
  list(draws = draws, means = means, vcov = V, sds = sds, rhat = rhat, ess = ess)
}

# Assemble the common `tobs_fit` skeleton from a `.tobs_pg_summarize` result.
# `extra` carries the family-specific tail (ms_community, intercepts,
# temporal_field, spatial_field); `spatial` overrides the default NULL slot.
.tobs_pg_finalize_fit <- function(summ, par_names, model, process_info, N,
                                  n.iter, n.chains, spatial = NULL,
                                  extra = list()) {
  rhat <- summ$rhat
  base <- list(
    draws        = summ$draws, means = summ$means, sds = summ$sds, vcov = summ$vcov,
    n_samples    = nrow(summ$draws), n_params = length(summ$means),
    log_prob     = rep(NA_real_, nrow(summ$draws)), log_lik = NA_real_,
    N            = N, rhat = rhat, ess = summ$ess,
    col_names    = par_names, param_names = par_names,
    n_fixed      = length(summ$means), fixed_names = par_names,
    process_info = process_info, model = model, spatial = spatial,
    method       = "pg_gibbs", n_chains = n.chains)
  structure(c(base, extra, list(
    convergence = list(converged = any(is.finite(rhat)) &&
                         max(rhat, na.rm = TRUE) < 1.1, n_iter = n.iter))),
    class = c("tobs_fit", "tulpa_fit"))
}
