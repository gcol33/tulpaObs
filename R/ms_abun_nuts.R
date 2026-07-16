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
#
# NON-CENTERED parameterisation: the per-species block holds standard-normal
# z_s, and the deviation is reconstructed per arm as b_{s,arm} = C_arm z_{s,arm}
# (C_arm the log-Cholesky factor of Sigma_arm). The community covariance then
# leaves the b-prior entirely -- z ~ N(0, I), no b' Sigma^{-1} b quadratic and no
# -0.5 S log|Sigma| normaliser -- and enters ONLY the data term through b = C z.
# This breaks the centered b<->Sigma funnel that otherwise saturates the NUTS
# treedepth (mean ~7, eps ~0.025 measured for the centered map); the chol
# gradient now flows from the data term via b = C z plus the hyperprior.
.tobs_ms_abun_nuts_logpost <- function(theta, margs, lay, priors,
                                       sigma.beta = 10, sigma.logr = 1.5,
                                       grad = TRUE) {
  P <- lay$P; S <- lay$n_species; is_nb <- lay$is_nb
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0

  # Cholesky factors per arm (Sigma_arm = C_arm C_arm'); logr arm is the 1x1
  # scalar SD = exp(chol_logr).
  C_lam <- .ms_ocs_chol_unpack(theta[lay$chol_lam], lay$p_lam)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$p_p)
  C_lr  <- if (is_nb) exp(theta[lay$chol_logr]) else NULL

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_lam <- matrix(0, lay$p_lam, lay$p_lam)
  A_p   <- matrix(0, lay$p_p,   lay$p_p)
  A_lr  <- 0

  # ---- data log-lik + inner gradient (per species), non-centered b = C z ----
  for (s in seq_len(S)) {
    bidx <- .tobs_ms_abun_nuts_b_idx(lay, s)
    z_s  <- theta[bidx]
    zl <- z_s[lay$lambda]; zp <- z_s[lay$p]
    bl   <- mu[lay$lambda] + as.numeric(C_lam %*% zl)
    bp   <- mu[lay$p]      + as.numeric(C_p   %*% zp)
    if (is_nb) { zr <- z_s[lay$logr]; r <- exp(mu[lay$logr] + C_lr * zr) }
    else        r <- Inf
    ev   <- margs[[s]]$eval_beta(bl, bp, r = r)
    lp   <- lp + ev$log_lik
    if (grad) {
      gl <- as.numeric(crossprod(margs[[s]]$X_lambda, ev$grad_eta_lambda)) # grad_b lambda
      gp <- as.numeric(crossprod(margs[[s]]$X_p,      ev$grad_eta_p))      # grad_b p
      g_mu[lay$lambda]    <- g_mu[lay$lambda] + gl
      g_mu[lay$p]         <- g_mu[lay$p]      + gp
      # z gradient (data part) = C' grad_b; the -z prior is added below.
      g[bidx[lay$lambda]] <- g[bidx[lay$lambda]] + as.numeric(crossprod(C_lam, gl))
      g[bidx[lay$p]]      <- g[bidx[lay$p]]      + as.numeric(crossprod(C_p,   gp))
      A_lam <- A_lam + outer(gl, zl)
      A_p   <- A_p   + outer(gp, zp)
      if (is_nb) {
        gt <- sum(ev$grad_theta)                                          # grad_b logr
        g_mu[lay$logr]    <- g_mu[lay$logr] + gt
        g[bidx[lay$logr]] <- g[bidx[lay$logr]] + C_lr * gt
        A_lr <- A_lr + gt * zr
      }
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior ----
  arms <- list(list(chol = lay$chol_lam, A = A_lam, C = C_lam, Pa = lay$p_lam),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$p_p))
  if (is_nb) arms <- c(arms, list(list(chol = lay$chol_logr,
                                       A = matrix(A_lr, 1, 1),
                                       C = matrix(C_lr, 1, 1), Pa = 1L)))
  for (arm in arms) {
    Pa <- arm$Pa; if (Pa == 0L) next
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, Pa) +
        pr$grad
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

# Data-term gradient w.r.t. one arm's packed log-Cholesky coordinates under the
# non-centered map b = C z. With A[i,j] = sum_s grad_b_{s,i} z_{s,j}, the packed
# gradient is A's lower triangle laid column-major (matching .ms_ocs_chol_pack),
# the diagonal scaled by C[j,j] (chain rule for the log-diagonal coordinate
# l_j = log C[j,j], since d b_j / d l_j = C[j,j] z_j).
.ms_abun_nuts_chol_data_grad <- function(A, C, Pa) {
  out <- numeric(.ms_ocs_chol_dim(Pa)); pos <- 0L
  for (j in seq_len(Pa)) {
    out[pos + 1L] <- A[j, j] * C[j, j]; pos <- pos + 1L
    if (j < Pa) {
      ni <- Pa - j
      out[pos + seq_len(ni)] <- A[(j + 1L):Pa, j]; pos <- pos + ni
    }
  }
  out
}

# Reconstruct the per-species deviation matrix b (S x P) from a packed coordinate
# vector under the non-centered map b_{s,arm} = C_arm z_{s,arm}.
.tobs_ms_abun_nuts_b_from_z <- function(theta, lay) {
  C_lam <- .ms_ocs_chol_unpack(theta[lay$chol_lam], lay$p_lam)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$p_p)
  C_lr  <- if (lay$is_nb) exp(theta[lay$chol_logr]) else NULL
  B <- matrix(0, lay$n_species, lay$P)
  for (s in seq_len(lay$n_species)) {
    z <- theta[.tobs_ms_abun_nuts_b_idx(lay, s)]
    B[s, lay$lambda] <- as.numeric(C_lam %*% z[lay$lambda])
    B[s, lay$p]      <- as.numeric(C_p   %*% z[lay$p])
    if (lay$is_nb) B[s, lay$logr] <- C_lr * z[lay$logr]
  }
  B
}

# Data-only log-likelihood (no priors) at a packed coordinate vector, summed over
# (species, site). Used for the fit's reported log_lik (the scale-invariant value
# logLik() / glance() read) at the posterior-mean coefficients.
.tobs_ms_abun_nuts_data_loglik <- function(theta, margs, lay) {
  mu <- theta[lay$mu]; tot <- 0
  B  <- .tobs_ms_abun_nuts_b_from_z(theta, lay)
  for (s in seq_len(lay$n_species)) {
    b_s <- B[s, ]
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
# community means, the community covariances as log-Cholesky coordinates (the 1x1
# dispersion covariance as log(sigma_log_r)), and -- under the non-centered map --
# the whitened per-species deviations z_s = C_arm^{-1} b_s (so reconstructing
# b = C z returns the warm BLUPs exactly). This is the NUTS initial position and
# the FD-gradient reference point.
.tobs_ms_abun_nuts_pack_init <- function(warm, lay) {
  theta <- numeric(lay$total)
  mu <- c(as.numeric(warm$mu_lambda), as.numeric(warm$mu_p))
  if (lay$is_nb) mu <- c(mu, as.numeric(warm$mu_log_r))
  theta[lay$mu] <- mu

  C_lam <- t(chol(.tobs_ms_abun_pd(as.matrix(warm$Sigma_lambda))))
  C_p   <- t(chol(.tobs_ms_abun_pd(as.matrix(warm$Sigma_p))))
  theta[lay$chol_lam] <- .ms_ocs_chol_pack(C_lam)
  theta[lay$chol_p]   <- .ms_ocs_chol_pack(C_p)
  s_lr <- if (lay$is_nb) max(as.numeric(warm$sigma_log_r), 1e-3) else NA_real_
  if (lay$is_nb) theta[lay$chol_logr] <- log(s_lr)   # 1x1 log-Cholesky = log(sd)

  bl <- as.matrix(warm$b_lambda); bp <- as.matrix(warm$b_p)
  blogr <- if (lay$is_nb) as.numeric(warm$b_logr) else NULL
  for (s in seq_len(lay$n_species)) {
    z_s <- numeric(lay$P)
    z_s[lay$lambda] <- forwardsolve(C_lam, bl[s, ])   # C z = b, C lower-tri
    z_s[lay$p]      <- forwardsolve(C_p,   bp[s, ])
    if (lay$is_nb) z_s[lay$logr] <- blogr[s] / s_lr
    theta[.tobs_ms_abun_nuts_b_idx(lay, s)] <- z_s
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
  user_K  <- !is.null(K_max)
  K_warm  <- if (user_K) as.integer(K_max) else max(lf$y) + 100L
  lay  <- .tobs_ms_abun_nuts_layout(p_lam, p_p, n_species, is_nb)
  pri  <- .tobs_ms_abun_nuts_priors()

  # Warm start at the Laplace-EM mode (Poisson: closed-form EM; NB: the joint_grad
  # / AGHQ path, matching the .tobs_fit_ms_nmix NB defaults so mu_log_r /
  # sigma_log_r / b_logr are returned).
  warm <- nmix_laplace_re(
    y = lf$y, site_idx = lf$site_idx, species_idx = lf$species_idx,
    X_lambda = X_lambda, X_p = lf$X_p, n_sites = model$n_sites,
    n_species = n_species, K_max = K_warm, max_iter = as.integer(max.iter),
    mixture = mix_code,
    optimizer = if (is_nb) "joint_grad" else "em",
    n_quad = if (is_nb) as.integer(n.quad) else 1L,
    lkj_eta = lkj_eta, verbose = FALSE)

  # Data-driven K_max for the (expensive) NUTS marginal: the latent N is summed
  # over [max(y), K_max] on EVERY leapfrog step, so a too-generous cap wastes the
  # inner loop. Cap at the abundance scale's 10-sigma upper tail (Poisson, or NB
  # variance lambda + lambda^2/r), where the dropped mass is < 1e-12 -- the
  # marginal is unchanged to machine precision while the per-step cost drops with
  # the latent range. A user-supplied K_max is respected verbatim.
  if (user_K) {
    K_max <- K_warm
  } else {
    bl <- as.matrix(warm$b_lambda)
    lam_max <- 0
    for (s in seq_len(n_species))
      lam_max <- max(lam_max, exp(as.numeric(X_lambda %*% (warm$mu_lambda + bl[s, ]))))
    var_max <- if (is_nb) {
      r_min <- max(min(as.numeric(warm$r_s)), 1e-6)
      lam_max + lam_max^2 / r_min
    } else lam_max
    K_tail <- as.integer(ceiling(lam_max + 10 * sqrt(var_max))) + 10L
    K_max  <- min(K_warm, max(max(lf$y) + 5L, K_tail))
  }
  K_max <- as.integer(K_max)

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

  # Per-species BLUPs = posterior mean of the RECONSTRUCTED deviation b = C z
  # (non-centered: the stored coordinates are the whitened z, so b is rebuilt per
  # draw with that draw's Cholesky factor before averaging).
  B_bar <- matrix(0, n_species, lay$P)
  for (i in seq_len(nrow(draws)))
    B_bar <- B_bar + .tobs_ms_abun_nuts_b_from_z(draws[i, ], lay)
  B_bar <- B_bar / nrow(draws)
  b_lambda <- B_bar[, lay$lambda, drop = FALSE]
  b_p      <- B_bar[, lay$p,      drop = FALSE]

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
    b_logr      <- B_bar[, lay$logr]   # reconstructed b = sd * z, posterior mean
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


# ---------------------------------------------------------------------------
# Spatial community N-mixture NUTS (shared fixed-hyper areal field; tulpaObs#73)
# ---------------------------------------------------------------------------

# Sample the exact joint posterior of a spatial community N-mixture: the
# community means, per-species deviations, community covariances, AND a SHARED
# fixed-hyper non-centered proper-CAR field on the abundance arm, jointly. The
# field precision tau Q(rho) is fixed at the nested-Laplace (#12 sfMsNMix)
# posterior mean; the whitened raw ~ N(0, I) with f = Linv %*% raw is sampled by
# the in-tree C++ FullGradFn (the shared-field block in src/ms_abun_nuts.cpp).
# proper-CAR only -- its full-rank precision gives a well-conditioned non-centered
# geometry; intrinsic icar/bym2 have a flat field-mean direction needing a
# sum-to-zero reparameterisation, so they stay on nested_laplace (the same gate
# the single-species abun() NUTS+areal path uses).
.tobs_fit_ms_abun_nuts_spatial <- function(model, spatial, mixture = "poisson",
                                           K_max = NULL, sigma.beta = 10,
                                           sigma.logr = 1.5, n.iter = 1000L,
                                           n.warmup = 1000L, n.chains = 1L,
                                           max.treedepth = 10L, adapt.delta = 0.9,
                                           seed = 1L, max.iter = 100L,
                                           verbose = FALSE) {
  .tobs_reject_weighted_spatial(spatial, "ms_abun NUTS abundance spatial")
  if (!identical(spatial$type, "car_proper")) {
    stop(sprintf(paste0(
      "ms_abun() NUTS + areal spatial supports the proper-CAR field ",
      "car_proper() (full-rank precision -> well-conditioned non-centered ",
      "geometry); the intrinsic '%s' field has a flat field-mean direction ",
      "needing a sum-to-zero reparameterisation for NUTS -- use ",
      "method = \"nested_laplace\" for the icar()/bym2() areal community fit. ",
      "(gcol33/tulpaObs#73)"), spatial$type), call. = FALSE)
  }
  if (!identical(mixture, "poisson")) {
    stop("ms_abun() NUTS + areal spatial is Poisson-only; negative-binomial ",
         "areal community N-mixture uses method = \"nested_laplace\".",
         call. = FALSE)
  }
  n_sites   <- model$n_sites
  n_species <- model$n_species
  if ((spatial$n_units %||% n_sites) != n_sites) {
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for ms_abun NUTS."),
                 spatial$n_units, n_sites), call. = FALSE)
  }
  X_lambda <- model$X_processes[[1L]]
  p_lam <- ncol(X_lambda); p_p <- model$process_info[[2L]]$p
  lf <- .tobs_ms_nmix_longform(model)
  K_warm <- if (is.null(K_max)) max(lf$y) + 100L else as.integer(K_max)
  lay  <- .tobs_ms_abun_nuts_layout(p_lam, p_p, n_species, FALSE)
  pri  <- .tobs_ms_abun_nuts_priors()
  csr  <- .nmix_spatial_csr(spatial)

  # Warm the field precision (tau, rho) + community means / covariances / field
  # from the nested-Laplace community-spatial (sfMsNMix) fit (#12).
  nl <- nmix_community_laplace_car_proper(
    lf = lf, X_lambda = X_lambda, n_sites = n_sites, n_species = n_species,
    csr = csr, n_spatial = n_sites, graph = spatial$graph, mixture = "P",
    K_max = K_warm, max_iter = as.integer(max.iter), verbose = FALSE)
  tau <- max(nl$hyper$tau[["mean"]], 1e-3)
  rho <- min(max(nl$hyper$rho[["mean"]], 0.01), 0.99)

  # Fixed field precision tau Q(rho) -> Linv = L^{-1} (f = Linv %*% raw).
  Q  <- .areal_Q(as.matrix(spatial$graph), rho)
  Qr <- tau * Q + diag(1e-4 * tau, n_sites)
  L  <- tryCatch(chol(Qr), error = function(e) NULL)
  if (is.null(L)) stop("ms_abun NUTS spatial: field precision not PD.",
                       call. = FALSE)
  Linv <- backsolve(L, diag(n_sites))

  # K_max for the NUTS marginal (data-driven, as the non-spatial path).
  if (!is.null(K_max)) {
    K_max <- K_warm
  } else {
    bl <- as.matrix(nl$b_lambda)
    lam_max <- 0
    for (s in seq_len(n_species))
      lam_max <- max(lam_max, exp(as.numeric(
        X_lambda %*% (nl$mu_lambda + bl[s, ])) + max(nl$spatial_field)))
    K_tail <- as.integer(ceiling(lam_max + 10 * sqrt(lam_max))) + 10L
    K_max  <- min(K_warm, max(max(lf$y) + 5L, K_tail))
  }
  K_max <- as.integer(K_max)

  # Warm the per-species block from the nested-Laplace community fit; pack raw0 so
  # f = Linv %*% raw0 returns the warm field (raw0 = L %*% f_warm).
  warm <- list(mu_lambda = nl$mu_lambda, mu_p = nl$mu_p,
               Sigma_lambda = nl$Sigma_lambda, Sigma_p = nl$Sigma_p,
               b_lambda = nl$b_lambda, b_p = nl$b_p)
  theta0_base <- .tobs_ms_abun_nuts_pack_init(warm, lay)
  raw0 <- as.numeric(L %*% nl$spatial_field)
  theta0 <- c(theta0_base, raw0)

  spec <- list(y = as.integer(lf$y), site_idx = as.integer(lf$site_idx),
               species_idx = as.integer(lf$species_idx),
               X_lambda = X_lambda, X_p = lf$X_p,
               n_sites = n_sites, n_species = n_species, K_max = K_max,
               is_nb = FALSE, n_field_units = n_sites,
               field_map = seq_len(n_sites), field_Linv = Linv)
  inv_metric <- .tobs_ms_abun_nuts_metric(spec, theta0, pri, sigma.beta,
                                          sigma.logr)

  run_chain <- function(ch) cpp_ms_abun_nuts(
    spec, theta0 = theta0, pri = pri, sigma_beta = sigma.beta,
    sigma_logr = sigma.logr, inv_metric = inv_metric,
    n_iter = as.integer(n.iter + n.warmup), n_warmup = as.integer(n.warmup),
    max_treedepth = as.integer(max.treedepth), adapt_delta = adapt.delta,
    seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))

  n_chains <- as.integer(n.chains); rhat_ess <- NULL
  if (n_chains > 1L) {
    rcs    <- lapply(seq_len(n_chains), run_chain)
    chains <- lapply(rcs, `[[`, "draws")
    rhat_ess <- .ms_ocs_rhat_ess(chains)
    draws  <- do.call(rbind, chains)
    accept <- unlist(lapply(rcs, `[[`, "accept_prob"))
    divergent <- unlist(lapply(rcs, `[[`, "divergent"))
    treedepth <- as.integer(unlist(lapply(rcs, `[[`, "treedepth")))
    epsilon   <- mean(vapply(rcs, `[[`, 0, "epsilon"))
  } else {
    res <- run_chain(1L)
    draws <- res$draws; accept <- res$accept_prob
    divergent <- res$divergent; treedepth <- as.integer(res$treedepth)
    epsilon <- res$epsilon
  }

  # Reconstruct the community block (drop the trailing raw columns).
  par     <- colMeans(draws)
  mu_hat  <- par[lay$mu]
  vcov_mu <- stats::cov(draws[, lay$mu, drop = FALSE])
  sig_mean <- function(cols, Pa) {
    acc <- matrix(0, Pa, Pa)
    for (i in seq_len(nrow(draws)))
      acc <- acc + tcrossprod(.ms_ocs_chol_unpack(draws[i, cols], Pa))
    acc <- acc / nrow(draws); (acc + t(acc)) / 2
  }
  Sigma_lambda <- sig_mean(lay$chol_lam, lay$p_lam)
  Sigma_p      <- sig_mean(lay$chol_p,   lay$p_p)
  B_bar <- matrix(0, n_species, lay$P)
  for (i in seq_len(nrow(draws)))
    B_bar <- B_bar + .tobs_ms_abun_nuts_b_from_z(draws[i, ], lay)
  B_bar <- B_bar / nrow(draws)

  # Field posterior mean f = Linv %*% mean(raw).
  raw_idx  <- lay$total + seq_len(n_sites)
  field_mean <- as.numeric(Linv %*% colMeans(draws[, raw_idx, drop = FALSE]))

  margs   <- .tobs_ms_abun_nuts_marginals(lf, X_lambda, n_sites, "P", K_max)
  ll_mean <- .tobs_ms_abun_nuts_data_loglik(par[seq_len(lay$total)], margs, lay)

  raw <- list(
    mu_lambda = mu_hat[lay$lambda], mu_p = mu_hat[lay$p],
    vcov = vcov_mu, Sigma_lambda = Sigma_lambda, Sigma_p = Sigma_p,
    b_lambda = B_bar[, lay$lambda, drop = FALSE],
    b_p = B_bar[, lay$p, drop = FALSE],
    log_lik = ll_mean, converged = TRUE, n_iter = NA_integer_,
    optimizer = "nuts", n_quad = 1L, lkj_eta = 1.5)

  fit <- build_ms_nmix_fit(raw, model, mixture = "poisson", spatial = spatial)
  fit$method <- "nuts"
  fit$log_prob <- rep(ll_mean, nrow(draws))
  fit$spatial_field <- field_mean
  fit$nuts <- list(
    draws = draws, layout = lay, n_field_units = n_sites,
    accept_prob = accept, divergent = divergent, treedepth = treedepth,
    epsilon = epsilon, n_chains = n_chains, divergent_total = sum(divergent),
    is_nb = FALSE, K_max = K_max, field_tau = tau, field_rho = rho,
    sigma_beta = sigma.beta, sigma_logr = sigma.logr)
  if (!is.null(rhat_ess)) {
    fit$nuts$rhat     <- rhat_ess$rhat
    fit$nuts$ess      <- rhat_ess$ess
    fit$nuts$max_rhat <- max(rhat_ess$rhat, na.rm = TRUE)
    fit$nuts$min_ess  <- min(rhat_ess$ess,  na.rm = TRUE)
  }
  fit
}
