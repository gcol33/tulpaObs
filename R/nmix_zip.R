# nmix_zip.R - Zero-inflated N-mixture (ZIP / ZINB) for abun(mixture = "zip" /
# "zinb"). Structural-zero on the latent abundance:
#   N_i = 0 with probability omega (structural), else N_i ~ Pois/NB(lambda_i),
#   y_ij | N_i ~ Binom(N_i, p_ij).
# The observed per-site marginal is a two-component mixture,
#   L_i = omega * 1{all y_i = 0} + (1 - omega) * L_royle_i,
# with L_royle_i the (Poisson / NB) Royle marginal that nmix_site_marginal()
# already exposes with its eta-level derivatives. This is a PURE-R additive layer
# over the C++ per-site Royle pieces -- no marginal-kernel change, so the plain
# abun() / removal() / distance() paths are untouched. omega is an intercept-only
# structural-zero probability (logit). Mode found by BFGS on the exact ZIP
# marginal; vcov = inverse of the numeric Hessian at the mode.
#
# Scope (v1, gcol33/tulpaObs#116): non-spatial laplace only, intercept-only
# omega. A zero-inflation covariate design, an areal field, and a NUTS path are
# follow-ups (the marginal + its gradient are already the additive layer they
# would share).

.tobs_fit_nmix_zip <- function(model, mixture = "zip", K_max = NULL,
                               max_iter = 300L, verbose = TRUE, ...) {
  is_nb    <- identical(mixture, "zinb")
  X_lambda <- model$X_processes[[1L]]
  X_p      <- model$X_processes[[2L]]
  y_long   <- as.integer(model$y_long)
  site_idx <- as.integer(model$site_idx)
  n_sites  <- nrow(X_lambda); p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  if (is.null(K_max)) K_max <- max(y_long) + 100L
  K_max <- as.integer(K_max)

  # Per-site all-zero indicator (no detection at any visit of the site).
  az  <- rep(TRUE, n_sites)
  agg <- tapply(y_long, site_idx, function(v) all(v == 0L))
  az[as.integer(names(agg))] <- as.logical(agg)

  mrg <- nmix_site_marginal(y = y_long, site_idx = site_idx,
                            X_lambda = X_lambda, X_p = X_p,
                            mixture = if (is_nb) "NB" else "P", K_max = K_max)

  # Exact ZIP negative log-likelihood at theta = [beta_lambda | beta_p |
  # logit_omega | (log_r if ZINB)].
  neg_ll <- function(theta) {
    bl <- theta[seq_len(p_lam)]
    bp <- theta[p_lam + seq_len(p_p)]
    lo <- theta[p_lam + p_p + 1L]
    r  <- if (is_nb) exp(theta[p_lam + p_p + 2L]) else Inf
    om <- stats::plogis(lo)
    ev <- tryCatch(mrg$eval_beta(bl, bp, r), error = function(e) NULL)
    if (is.null(ev)) return(1e10)
    llr <- ev$log_lik_site
    if (any(!is.finite(llr))) return(1e10)
    log1m <- log1p(-om)                       # log(1 - omega)
    ll <- numeric(n_sites)
    ll[!az] <- log1m + llr[!az]               # a detection rules out N = 0
    a  <- log1m + llr[az]; b <- log(om)       # all-zero: mix in the structural 0
    mx <- pmax(a, b)
    ll[az] <- mx + log(exp(a - mx) + exp(b - mx))
    val <- -sum(ll)
    if (is.finite(val)) val else 1e10
  }

  # Warm start: the no-ZI Royle betas, plus a structural-zero logit seeded from
  # a modest share of the observed all-zero sites.
  warm <- tryCatch(
    nmix_laplace(y = y_long, site_idx = site_idx, X_lambda = X_lambda,
                 X_p = X_p, mixture = if (is_nb) "NB" else "P",
                 K_max = K_max, verbose = FALSE),
    error = function(e) NULL)
  bl0 <- if (!is.null(warm)) as.numeric(warm$beta_lambda)
         else c(log(mean(y_long) + 0.5), rep(0, p_lam - 1L))
  bp0 <- if (!is.null(warm)) as.numeric(warm$beta_p) else rep(0, p_p)
  lo0 <- stats::qlogis(min(max(mean(az) * 0.3, 0.02), 0.6))
  theta0 <- c(bl0, bp0, lo0)
  if (is_nb) {
    lr0 <- if (!is.null(warm) && is.finite(warm$log_r %||% NA_real_)) warm$log_r
           else log(2)
    theta0 <- c(theta0, lr0)
  }

  # Box constraints only on the pathological corners: the structural-zero logit
  # and (ZINB) log_r can run away because a huge NB overdispersion (r -> 0) mimics
  # structural zeros -- the well-known ZINB zero-source confounding. Betas stay
  # unbounded; the bounds do not bias an interior (identified) fit, they only stop
  # a degenerate seed from diverging. L-BFGS-B still yields the numeric Hessian at
  # the mode for the observed-information vcov.
  lo_b <- c(rep(-Inf, p_lam + p_p), -8)
  hi_b <- c(rep( Inf, p_lam + p_p),  8)
  if (is_nb) { lo_b <- c(lo_b, log(1e-2)); hi_b <- c(hi_b, log(1e3)) }
  opt <- stats::optim(theta0, neg_ll, method = "L-BFGS-B",
                      lower = lo_b, upper = hi_b,
                      control = list(maxit = max_iter, factr = 1e7),
                      hessian = TRUE)

  p_tot <- length(theta0)
  vcov  <- tryCatch(solve(opt$hessian), error = function(e) {
    d <- diag(opt$hessian); d[d <= 0] <- NA_real_
    diag(1 / d, p_tot)
  })

  est <- opt$par
  raw <- list(
    beta_lambda = est[seq_len(p_lam)],
    beta_p      = est[p_lam + seq_len(p_p)],
    logit_omega = est[p_lam + p_p + 1L],
    mixture     = if (is_nb) "NB" else "P",
    vcov        = vcov,
    log_lik     = -opt$value,
    K_max       = K_max,
    converged   = opt$convergence == 0L,
    n_iter      = NA_integer_)
  if (is_nb) { raw$log_r <- est[p_lam + p_p + 2L]; raw$r <- exp(raw$log_r) }

  build_nmix_fit(raw, model, spatial = NULL)
}
