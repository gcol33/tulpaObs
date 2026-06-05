# ms_abun_nuts.R - NUTS target density for the community / multispecies
# N-mixture (ms_abun()).
#
# The Laplace-EM fit (ms_abun.R -> nmix_laplace_re) profiles the per-species
# coefficient deviations and the community covariances out by an arrowhead
# Newton + closed-form covariance M-step, then reports a Gaussian community-mean
# posterior c(mu_lambda, mu_p[, mu_log_r]). NUTS instead samples the EXACT joint
# posterior -- the community means, the per-species deviations {b_s}, AND the
# community covariances (Sigma_lambda, Sigma_p[, sigma_log_r]) -- which removes
# the Gaussian approximation, gives calibrated (non-Gaussian) community
# intervals, and yields the per-(species, site) pointwise likelihood WAIC / LOO
# need.
#
# The target is the community generalisation of the single-species abun() NUTS
# target (R/abun_nuts.R): the data term is the same Royle (2004) per-site
# marginal (nmix_site_marginal()), but now summed over (species, site) with the
# per-species coefficient eta_{s} = X . (mu + b_s), and a hierarchical Gaussian
# community prior b_{s,arm} ~ N(0, Sigma_arm) replaces the single flat coefficient
# ridge. The community covariances are carried by their log-Cholesky factors (the
# same parametrisation and helpers as the spatial-factor community occu_cover
# target, R/ms_occu_cover_spatial_nuts.R). The joint log-posterior is
#
#   log p = sum_{s,i} log m_{s,i}(theta)               # per-species-site marginal
#         - 0.5 ||mu_coef||^2 / sigma.beta^2           # community-mean priors
#         [ - 0.5 mu_log_r^2 / sigma.logr^2 ]          # (NB)
#         - 0.5 sum_s b_{s,arm}' Sigma_arm^{-1} b_{s,arm}   (per arm)
#         - 0.5 S sum_arm log|Sigma_arm|               # MVN normalisers
#         + log p(Sigma log-Cholesky coords)           # weakly-informative hyperpriors
#
# Block (data + b-quadratic + mu prior + log-det) and the chol coordinate
# gradients all mirror the C++ FullGradFn (src/ms_abun_nuts.cpp) byte-for-byte;
# this R version is the oracle the C++ port is cross-checked against. The
# log-Cholesky packing / unpacking / block-gradient helpers (.ms_ocs_chol_*) are
# shared with the spatial-factor community target.


# ---------------------------------------------------------------------------
# Parameter layout
# ---------------------------------------------------------------------------

# Packed NUTS coordinate layout for the community N-mixture:
#   theta = ( mu [P], {b_s} species-major [S*P], chol_lambda [q_lam],
#             chol_p [q_p] [, chol_logr [1]] )
# with P = p_lam + p_p [+ 1 under NB], mu = (mu_lambda, mu_p[, mu_log_r]),
# b_s = (b_lambda_s, b_p_s[, b_logr_s]). The trailing 1x1 chol_logr block is the
# log-Cholesky of the scalar log-dispersion community covariance sigma_log_r^2
# (its single packed coordinate is log(sigma_log_r)). `lambda` / `p` / `logr` are
# the within-arm coordinate indices (used to slice both mu and each b_s).
.tobs_ms_abun_nuts_layout <- function(p_lam, p_p, n_species, is_nb) {
  P     <- p_lam + p_p + (if (is_nb) 1L else 0L)
  q_lam <- as.integer(p_lam * (p_lam + 1L) / 2L)
  q_p   <- as.integer(p_p   * (p_p   + 1L) / 2L)
  q_logr <- if (is_nb) 1L else 0L
  b_off        <- P
  chol_lam_off <- P + n_species * P
  chol_p_off   <- chol_lam_off + q_lam
  chol_logr_off <- chol_p_off + q_p
  total <- chol_logr_off + q_logr
  list(
    P = P, p_lam = p_lam, p_p = p_p, n_species = n_species, is_nb = isTRUE(is_nb),
    q_lam = q_lam, q_p = q_p, q_logr = q_logr,
    lambda = seq_len(p_lam),
    p      = p_lam + seq_len(p_p),
    logr   = if (is_nb) p_lam + p_p + 1L else integer(0),
    mu     = seq_len(P),
    b_off  = b_off,
    chol_lam  = chol_lam_off  + seq_len(q_lam),
    chol_p    = chol_p_off    + seq_len(q_p),
    chol_logr = if (is_nb) chol_logr_off + seq_len(q_logr) else integer(0),
    total = total)
}

# Coordinate indices of species s's b_s sub-vector within the packed theta.
.tobs_ms_abun_nuts_b_idx <- function(lay, s) {
  lay$b_off + (s - 1L) * lay$P + seq_len(lay$P)
}


# ---------------------------------------------------------------------------
# Hyperprior specification
# ---------------------------------------------------------------------------

# Weakly-informative priors on the sampled log-Cholesky coordinates of the
# community covariances: the log-diagonal carries a Normal centred at log(0.5)
# (community SDs of order 0.5 on the link scale) with a wide SD; the off-diagonals
# a mean-zero Normal shrinking community correlations toward independence. Shared
# in spirit with the spatial-factor target's chol priors; the field / dispersion
# axes there are absent here.
.tobs_ms_abun_nuts_priors <- function() {
  list(chol_logdiag_mean = log(0.5), chol_logdiag_sd = 1.5,
       chol_offdiag_sd   = 1.0)
}


# ---------------------------------------------------------------------------
# Per-species marginal kernels
# ---------------------------------------------------------------------------

# Build one nmix_site_marginal() per species from the stacked long form. Each
# species' marginal carries the FULL site-level abundance design X_lambda (so its
# grad_eta_lambda is length n_sites) and that species' detection-design slice;
# eval_beta() then differentiates through both arms exactly as the single-species
# abun() NUTS target does, once per species.
.tobs_ms_abun_nuts_marginals <- function(lf, X_lambda, n_sites, mix_code, K_max) {
  S <- max(lf$species_idx)
  lapply(seq_len(S), function(s) {
    idx <- which(lf$species_idx == s)
    nmix_site_marginal(
      y        = lf$y[idx],
      site_idx = lf$site_idx[idx],
      X_lambda = X_lambda,
      X_p      = lf$X_p[idx, , drop = FALSE],
      mixture  = mix_code,
      K_max    = K_max)
  })
}


# ---------------------------------------------------------------------------
# Joint log-posterior + gradient (the NUTS target density / oracle)
# ---------------------------------------------------------------------------

# Full-vector joint log-posterior and its gradient for the community N-mixture.
# `margs` is the per-species list from .tobs_ms_abun_nuts_marginals(); `lay` the
# layout. Returns list(lp, grad) over the packed coordinates. Mirrors the C++
# ms_abun_nuts_eval (src/ms_abun_nuts.cpp).
.tobs_ms_abun_nuts_logpost <- function(theta, margs, lay, priors,
                                       sigma.beta = 10, sigma.logr = 1.5,
                                       grad = TRUE) {
  P <- lay$P; S <- lay$n_species; is_nb <- lay$is_nb
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0
  b_list <- vector("list", S)

  # ---- data log-lik + inner gradient (per species) ----
  for (s in seq_len(S)) {
    bidx <- .tobs_ms_abun_nuts_b_idx(lay, s)
    b_s  <- theta[bidx]; b_list[[s]] <- b_s
    bl   <- mu[lay$lambda] + b_s[lay$lambda]
    bp   <- mu[lay$p]      + b_s[lay$p]
    r    <- if (is_nb) exp(mu[lay$logr] + b_s[lay$logr]) else Inf
    ev   <- margs[[s]]$eval_beta(bl, bp, r = r)
    lp   <- lp + ev$log_lik
    if (grad) {
      gl <- as.numeric(crossprod(margs[[s]]$X_lambda, ev$grad_eta_lambda))
      gp <- as.numeric(crossprod(margs[[s]]$X_p,      ev$grad_eta_p))
      g_mu[lay$lambda]    <- g_mu[lay$lambda] + gl
      g_mu[lay$p]         <- g_mu[lay$p]      + gp
      g[bidx[lay$lambda]] <- g[bidx[lay$lambda]] + gl
      g[bidx[lay$p]]      <- g[bidx[lay$p]]      + gp
      if (is_nb) {
        gt <- sum(ev$grad_theta)
        g_mu[lay$logr]    <- g_mu[lay$logr] + gt
        g[bidx[lay$logr]] <- g[bidx[lay$logr]] + gt
      }
    }
  }

  # ---- community covariance: per-arm b-quadratic + log-det normaliser + chol
  #      block gradient; accumulate the b-prior into g. ----
  arms <- list(list(coords = lay$lambda, chol = lay$chol_lam, Pa = lay$p_lam),
               list(coords = lay$p,      chol = lay$chol_p,   Pa = lay$p_p))
  if (is_nb) arms <- c(arms, list(list(coords = lay$logr, chol = lay$chol_logr,
                                       Pa = 1L)))
  for (arm in arms) {
    Pa <- arm$Pa; if (Pa == 0L) next
    C  <- .ms_ocs_chol_unpack(theta[arm$chol], Pa)
    Si <- chol2inv(t(C))
    logdet <- 2 * sum(log(diag(C)))
    M <- matrix(0, Pa, Pa); quad <- 0
    for (s in seq_len(S)) {
      bb  <- b_list[[s]][arm$coords]
      sib <- as.numeric(Si %*% bb)
      quad <- quad + sum(bb * sib)
      M    <- M + outer(bb, bb)
      if (grad) {
        bidx <- .tobs_ms_abun_nuts_b_idx(lay, s)[arm$coords]
        g[bidx] <- g[bidx] - sib
      }
    }
    lp <- lp - 0.5 * quad - 0.5 * S * logdet
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_ocs_chol_block_grad(C, M, S, theta[arm$chol],
                                                     Pa, priors)
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  coef_idx <- c(lay$lambda, lay$p)
  lp <- lp - 0.5 * ib2 * sum(mu[coef_idx]^2)
  g_mu[coef_idx] <- g_mu[coef_idx] - ib2 * mu[coef_idx]
  if (is_nb) {
    il2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * il2 * mu[lay$logr]^2
    g_mu[lay$logr] <- g_mu[lay$logr] - il2 * mu[lay$logr]
  }
  g[lay$mu] <- g[lay$mu] + g_mu

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Data-only log-likelihood (no priors) at a packed coordinate vector, summed over
# (species, site). Used for the fit's reported log_lik (the scale-invariant value
# logLik() / glance() read) at the posterior-mean coefficients.
.tobs_ms_abun_nuts_data_loglik <- function(theta, margs, lay) {
  mu <- theta[lay$mu]; tot <- 0
  for (s in seq_len(lay$n_species)) {
    b_s <- theta[.tobs_ms_abun_nuts_b_idx(lay, s)]
    r   <- if (lay$is_nb) exp(mu[lay$logr] + b_s[lay$logr]) else Inf
    tot <- tot + margs[[s]]$eval_beta(mu[lay$lambda] + b_s[lay$lambda],
                                      mu[lay$p]      + b_s[lay$p], r = r)$log_lik
  }
  tot
}


# ---------------------------------------------------------------------------
# Warm-start packing + inverse-mass metric
# ---------------------------------------------------------------------------

# Symmetrise + (if needed) jitter a covariance to a Cholesky-able PD matrix.
.tobs_ms_abun_pd <- function(S) {
  S <- (S + t(S)) / 2
  ok <- tryCatch({ chol(S); TRUE }, error = function(e) FALSE)
  if (ok) return(S)
  S + diag(max(1e-6, 1e-6 * mean(diag(S))), ncol(S))
}

# Pack a nmix_laplace_re() community fit into the full NUTS coordinate vector: the
# community means, the per-species BLUP deviations, and the community covariances
# as log-Cholesky coordinates (the 1x1 dispersion covariance as log(sigma_log_r)).
# This is the NUTS initial position and the FD-gradient reference point.
.tobs_ms_abun_nuts_pack_init <- function(warm, lay) {
  theta <- numeric(lay$total)
  mu <- c(as.numeric(warm$mu_lambda), as.numeric(warm$mu_p))
  if (lay$is_nb) mu <- c(mu, as.numeric(warm$mu_log_r))
  theta[lay$mu] <- mu

  bl <- as.matrix(warm$b_lambda); bp <- as.matrix(warm$b_p)
  blogr <- if (lay$is_nb) as.numeric(warm$b_logr) else NULL
  for (s in seq_len(lay$n_species)) {
    b_s <- c(bl[s, ], bp[s, ])
    if (lay$is_nb) b_s <- c(b_s, blogr[s])
    theta[.tobs_ms_abun_nuts_b_idx(lay, s)] <- b_s
  }

  chol_coords <- function(Sig) .ms_ocs_chol_pack(t(chol(.tobs_ms_abun_pd(Sig))))
  theta[lay$chol_lam] <- chol_coords(as.matrix(warm$Sigma_lambda))
  theta[lay$chol_p]   <- chol_coords(as.matrix(warm$Sigma_p))
  if (lay$is_nb) {
    s_lr <- max(as.numeric(warm$sigma_log_r), 1e-3)
    theta[lay$chol_logr] <- log(s_lr)        # 1x1 log-Cholesky = log(sd)
  }
  theta
}

# Inverse-mass diagonal for the NUTS warm start: the posterior variance per
# coordinate from the finite-difference diagonal of the joint log-posterior
# Hessian at the mode (the Laplace metric), via the fast C++ gradient.
.tobs_ms_abun_nuts_metric <- function(spec, theta, pri, sigma.beta, sigma.logr,
                                      h = 1e-4) {
  np <- length(theta); md <- numeric(np)
  g <- function(th) cpp_ms_abun_nuts_joint_logpost(spec, th, pri,
                                                   sigma.beta, sigma.logr)$grad
  for (j in seq_len(np)) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    md[j] <- -(g(tp)[j] - g(tm)[j]) / (2 * h)
  }
  1 / pmax(md, 1e-3)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the community N-mixture
# ---------------------------------------------------------------------------

# Sample the exact joint posterior of a non-spatial community N-mixture via
# tulpa's NUTS engine and the in-tree C++ FullGradFn (cpp_ms_abun_nuts),
# warm-started at the Laplace-EM mode with a diagonal Laplace metric, then package
# the draws into the build_ms_nmix_fit shape so coef / vcov / confint / ranef /
# simulate read the NUTS posterior. The community means / their covariance come
# from the posterior of the mu block; the community covariances Sigma_arm from the
# posterior mean of C C' over the chol draws; the per-species BLUPs from the
# posterior mean of {b_s}. The full per-draw parameter vector is kept under
# `fit$nuts$draws` (with the layout) so the calibrated per-(species, site) WAIC /
# LOO can be scored from it.
.tobs_fit_ms_abun_nuts <- function(model, mixture = "poisson", K_max = NULL,
                                   sigma.beta = 10, sigma.logr = 1.5,
                                   n.iter = 1000L, n.warmup = 1000L,
                                   n.chains = 1L, max.treedepth = 10L,
                                   adapt.delta = 0.9, seed = 1L,
                                   lkj_eta = 1.5, n.quad = 5L, max.iter = 100L,
                                   verbose = FALSE) {
  if (!identical(mixture, "poisson") && !identical(mixture, "negbin")) {
    stop("Community N-mixture NUTS supports mixture = \"poisson\" or \"negbin\" ",
         "(got \"", mixture, "\").", call. = FALSE)
  }
  is_nb    <- identical(mixture, "negbin")
  mix_code <- if (is_nb) "NB" else "P"

  X_lambda <- model$X_processes[[1]]
  p_lam    <- ncol(X_lambda)
  p_p      <- model$process_info[[2]]$p
  n_species <- model$n_species
  lf <- .tobs_ms_nmix_longform(model)
  if (is.null(K_max)) K_max <- max(lf$y) + 100L
  K_max <- as.integer(K_max)
  lay  <- .tobs_ms_abun_nuts_layout(p_lam, p_p, n_species, is_nb)
  pri  <- .tobs_ms_abun_nuts_priors()

  # Warm start at the Laplace-EM mode (Poisson: closed-form EM; NB: the joint_grad
  # / AGHQ path, matching the .tobs_fit_ms_nmix NB defaults so mu_log_r /
  # sigma_log_r / b_logr are returned).
  warm <- nmix_laplace_re(
    y = lf$y, site_idx = lf$site_idx, species_idx = lf$species_idx,
    X_lambda = X_lambda, X_p = lf$X_p, n_sites = model$n_sites,
    n_species = n_species, K_max = K_max, max_iter = as.integer(max.iter),
    mixture = mix_code,
    optimizer = if (is_nb) "joint_grad" else "em",
    n_quad = if (is_nb) as.integer(n.quad) else 1L,
    lkj_eta = lkj_eta, verbose = FALSE)

  theta0 <- .tobs_ms_abun_nuts_pack_init(warm, lay)
  spec <- list(y = as.integer(lf$y), site_idx = as.integer(lf$site_idx),
               species_idx = as.integer(lf$species_idx),
               X_lambda = X_lambda, X_p = lf$X_p,
               n_sites = model$n_sites, n_species = n_species,
               K_max = K_max, is_nb = is_nb)
  inv_metric <- .tobs_ms_abun_nuts_metric(spec, theta0, pri, sigma.beta,
                                          sigma.logr)

  run_chain <- function(ch) {
    cpp_ms_abun_nuts(
      spec, theta0 = theta0, pri = pri,
      sigma_beta = sigma.beta, sigma_logr = sigma.logr,
      inv_metric = inv_metric,
      n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup),
      max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta,
      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }

  n_chains <- as.integer(n.chains)
  rhat_ess <- NULL
  if (n_chains > 1L) {
    rcs    <- lapply(seq_len(n_chains), run_chain)
    chains <- lapply(rcs, `[[`, "draws")
    rhat_ess <- .ms_ocs_rhat_ess(chains)
    draws  <- do.call(rbind, chains)
    accept    <- unlist(lapply(rcs, `[[`, "accept_prob"))
    divergent <- unlist(lapply(rcs, `[[`, "divergent"))
    treedepth <- as.integer(unlist(lapply(rcs, `[[`, "treedepth")))
    epsilon   <- mean(vapply(rcs, `[[`, 0, "epsilon"))
  } else {
    res <- run_chain(1L)
    draws     <- res$draws
    accept    <- res$accept_prob
    divergent <- res$divergent
    treedepth <- as.integer(res$treedepth)
    epsilon   <- res$epsilon
  }

  # ---- reconstruct the build_ms_nmix_fit `raw` shape from the draws ----
  par    <- colMeans(draws)
  mu_hat <- par[lay$mu]
  vcov_mu <- stats::cov(draws[, lay$mu, drop = FALSE])

  sig_mean <- function(cols, Pa) {
    acc <- matrix(0, Pa, Pa)
    for (i in seq_len(nrow(draws)))
      acc <- acc + tcrossprod(.ms_ocs_chol_unpack(draws[i, cols], Pa))
    acc <- acc / nrow(draws); (acc + t(acc)) / 2
  }
  Sigma_lambda <- sig_mean(lay$chol_lam, lay$p_lam)
  Sigma_p      <- sig_mean(lay$chol_p,   lay$p_p)

  b_lambda <- matrix(0, n_species, lay$p_lam)
  b_p      <- matrix(0, n_species, lay$p_p)
  for (s in seq_len(n_species)) {
    b_s <- par[.tobs_ms_abun_nuts_b_idx(lay, s)]
    b_lambda[s, ] <- b_s[lay$lambda]
    b_p[s, ]      <- b_s[lay$p]
  }

  margs   <- .tobs_ms_abun_nuts_marginals(lf, X_lambda, model$n_sites, mix_code,
                                          K_max)
  ll_mean <- .tobs_ms_abun_nuts_data_loglik(par, margs, lay)

  raw <- list(
    mu_lambda = mu_hat[lay$lambda], mu_p = mu_hat[lay$p],
    vcov = vcov_mu, Sigma_lambda = Sigma_lambda, Sigma_p = Sigma_p,
    b_lambda = b_lambda, b_p = b_p,
    log_lik = ll_mean, converged = TRUE, n_iter = NA_integer_,
    optimizer = "nuts", n_quad = if (is_nb) as.integer(n.quad) else 1L,
    lkj_eta = lkj_eta)
  if (is_nb) {
    mu_log_r    <- unname(mu_hat[lay$logr])
    sigma_log_r <- mean(exp(draws[, lay$chol_logr]))
    b_logr <- vapply(seq_len(n_species),
                     function(s) par[.tobs_ms_abun_nuts_b_idx(lay, s)][lay$logr],
                     0)
    raw$mu_log_r    <- mu_log_r
    raw$sigma_log_r <- sigma_log_r
    raw$b_logr      <- b_logr
    raw$r_s         <- exp(mu_log_r + b_logr)
  }

  fit <- build_ms_nmix_fit(raw, model, mixture = mixture, spatial = NULL)
  fit$method <- "nuts"
  fit$log_prob <- rep(ll_mean, nrow(draws))

  fit$nuts <- list(
    draws = draws, layout = lay,
    accept_prob = accept, divergent = divergent, treedepth = treedepth,
    epsilon = epsilon, n_chains = n_chains,
    divergent_total = sum(divergent),
    is_nb = is_nb, K_max = K_max,
    sigma_beta = sigma.beta, sigma_logr = sigma.logr)
  if (!is.null(rhat_ess)) {
    fit$nuts$rhat     <- rhat_ess$rhat
    fit$nuts$ess      <- rhat_ess$ess
    fit$nuts$max_rhat <- max(rhat_ess$rhat, na.rm = TRUE)
    fit$nuts$min_ess  <- min(rhat_ess$ess,  na.rm = TRUE)
  }
  fit
}
