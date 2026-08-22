# cover_nuts.R - NUTS target for the non-spatial standalone cover hurdle family
# (cover()).
#
# The cover hurdle is occu_cover() minus the occupancy / detection latent
# mixture: there is no z to marginalise, and the two arms are conditionally
# independent given the data. The Laplace fit (fit_cover_hurdle) runs two
# independent tulpa_laplace() calls and returns a Gaussian observed-Fisher
# posterior per arm. NUTS instead samples the exact joint posterior of the
# packed coefficient vector
#   theta = c(beta_presence, beta_positive, log_dispersion),
# giving calibrated (non-Gaussian) intervals and the per-draw pointwise
# likelihood WAIC / LOO need. There is no latent field, random effect, or
# community covariance on the non-spatial path, so the parameter vector is just
# the two flat coefficient blocks plus one log-dispersion scalar.
#
# .tobs_cover_nuts_logpost is the R oracle: it recomputes the joint
# log-posterior and gradient exactly as the C++ FullGradFn (src/cover_nuts.cpp,
# cpp_cover_nuts) does, and a byte-exact test cross-checks the two before the
# sampler is trusted. The positive arm reuses the same beta / lognormal closed
# forms the Laplace path uses (BetaPositive / LognormalPositive in
# occu_coupling_shared.h).


# Joint log-posterior + gradient of the non-spatial cover-hurdle coefficient
# vector theta = c(beta_presence, beta_pos, log_disp). Weak Gaussian priors
# N(0, sigma.beta^2) on every coefficient and a broad N(0, sigma.logdisp^2) on
# log_disp keep the optimum proper without materially shifting the data-dominated
# mode. `enc` is a natural-scale cover encoding (encode_cover_hurdle with
# autoscale = FALSE): enc$occ_data carries the presence Bernoulli arm
# (`y` = 1{cover > 0}, `X`), enc$pos_data the positive arm (`y` = the arm
# response the policy consumes, `X`). The positive-arm response is
# exp(enc$pos_data$y) for lognormal (the encoder stores log(cover); the policy
# takes the log internally) and enc$pos_data$y for beta. Returns list(lp, grad)
# over the packed coordinates. This is the oracle the C++ FullGradFn mirrors.
.tobs_cover_nuts_logpost <- function(theta, enc, sigma.beta = 5,
                                     sigma.logdisp = 5) {
  X_pres <- enc$occ_data$X
  X_pos  <- enc$pos_data$X
  p_pres <- ncol(X_pres)
  p_pos  <- ncol(X_pos)
  total  <- p_pres + p_pos + 1L
  pos_code <- .occu_cover_pos_code(enc$positive)

  b_pres   <- theta[seq_len(p_pres)]
  b_pos    <- theta[p_pres + seq_len(p_pos)]
  log_disp <- theta[total]
  disp     <- exp(log_disp)

  present <- enc$occ_data$y
  y_pos   <- .tobs_cover_pos_response(enc)
  sgm     <- function(e) 1 / (1 + exp(-e))

  # Presence arm: plain Bernoulli on 1{cover > 0}.
  eta_pres <- as.numeric(X_pres %*% b_pres)
  pr       <- sgm(eta_pres)
  lp       <- sum(ifelse(present == 1L, log(pr), log(1 - pr)))
  g_eta_pres <- ifelse(present == 1L, 1 - pr, -pr)
  grad_pres  <- as.numeric(crossprod(X_pres, g_eta_pres))

  # Positive arm: beta / lognormal density at the present rows.
  eta_pos    <- as.numeric(X_pos %*% b_pos)
  g_eta_pos  <- numeric(length(eta_pos))
  g_logdisp  <- 0
  if (pos_code == 3L) {              # beta
    mu   <- sgm(eta_pos)
    a    <- mu * disp; b <- (1 - mu) * disp
    ly   <- log(y_pos); l1my <- log(1 - y_pos)
    lp   <- lp + sum(lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my)
    g_eta_pos <- disp * mu * (1 - mu) * (-digamma(a) + digamma(b) + ly - l1my)
    g_logdisp <- sum(disp * (digamma(disp) - mu * digamma(a) -
                             (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my))
  } else if (pos_code == 4L) {       # identity-Gaussian (#112): raw response
    sig <- disp; r <- (y_pos - eta_pos) / sig
    lp  <- lp + sum(-log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r)
    g_eta_pos <- r / sig
    g_logdisp <- sum(r * r - 1)
  } else {                           # lognormal: Gaussian on log(cover)
    sig <- disp; r <- (log(y_pos) - eta_pos) / sig
    lp  <- lp + sum(-log(y_pos) - log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r)
    g_eta_pos <- r / sig
    g_logdisp <- sum(r * r - 1)
  }
  grad_pos <- as.numeric(crossprod(X_pos, g_eta_pos))

  grad <- c(grad_pres, grad_pos, g_logdisp)

  # Weak Gaussian priors.
  ib2 <- 1 / sigma.beta^2
  nb  <- total - 1L
  bv  <- theta[seq_len(nb)]
  lp  <- lp - 0.5 * ib2 * sum(bv^2)
  grad[seq_len(nb)] <- grad[seq_len(nb)] - ib2 * bv
  ild2 <- 1 / sigma.logdisp^2
  lp   <- lp - 0.5 * ild2 * log_disp^2
  grad[total] <- grad[total] - ild2 * log_disp
  list(lp = lp, grad = grad)
}

# Positive-arm response the beta / lognormal policy consumes. The cover encoder
# stores log(cover) on the lognormal arm and the raw cover in (0, 1) on the beta
# arm; the policy takes the log internally for lognormal, so exponentiate it back
# to the natural cover scale here. Single source for both the R oracle and the
# C++ spec builder.
.tobs_cover_pos_response <- function(enc) {
  # beta stores the raw cover in (0, 1); gaussian stores the raw unbounded
  # response (#112); lognormal stores log(cover), so exponentiate it back.
  if (enc$positive %in% c("beta", "gaussian")) enc$pos_data$y
  else exp(enc$pos_data$y)
}

# Build the C++ NUTS spec list from a natural-scale cover encoding. The presence
# and positive designs are passed straight through (raw natural scale), so the
# draws land on the natural coefficient scale.
.tobs_cover_nuts_spec <- function(enc) {
  list(
    pos_code = .occu_cover_pos_code(enc$positive),
    present = as.integer(enc$occ_data$y),
    y_pos   = as.numeric(.tobs_cover_pos_response(enc)),
    X_pres  = enc$occ_data$X,
    X_pos   = enc$pos_data$X
  )
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the non-spatial cover hurdle
# ---------------------------------------------------------------------------

# Front-door from .dispatch_cover: splat the sampler control keys (n.iter /
# n.warmup / n.chains / adapt.delta / max.treedepth / seed / verbose / progress*)
# into the NUTS fitter. Matches the occu_cover NUTS dispatch (control names are
# dotted and align with the fitter's formals; unused keys fall into `...`).
.tobs_fit_cover_nuts_dispatch <- function(formula, data, y, positive, family,
                                          priors, control) {
  do.call(.tobs_fit_cover_nuts,
          c(list(formula = formula, data = data, y = y, positive = positive,
                 family = family, priors = priors), control))
}

# Sample the exact non-spatial cover-hurdle coefficient posterior via tulpa's
# NUTS engine and the in-tree C++ FullGradFn (cpp_cover_nuts), warm-started at
# the Laplace mode with a diagonal Laplace metric, then package the draws into
# the same cover_fit shape fit_cover_hurdle / decode_cover_hurdle return so
# print / summary / predict / WAIC read the NUTS posterior. The cover hurdle is
# fit on a natural-scale encoding here (autoscale = FALSE) so the warm-start
# mode, the sampler draws, and the returned betas are all on the natural
# coefficient scale with no unscale round-trip. `sigma.logdisp` is an internal
# weak-prior width (no control knob, like abun's sigma.logr). `...` absorbs
# unused sampler controls (n.thin / n.threads / progress.*).
.tobs_fit_cover_nuts <- function(formula, data, y, positive, family,
                                 priors = NULL,
                                 sigma.beta = NULL, sigma.logdisp = 5,
                                 n.iter = NULL, n.warmup = NULL,
                                 n.chains = NULL, max.treedepth = NULL,
                                 adapt.delta = NULL, seed = NULL,
                                 verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts")

  # Natural-scale encoding: the NUTS spec and the warm-start mode share the same
  # (unscaled) design, so the draws need no unscale transform.
  enc <- encode_cover_hurdle(formula, data, y, positive = positive,
                             autoscale = FALSE)
  if (!is.null(enc$spatial_spec) || !is.null(enc$trend) || !is.null(enc$mcar) ||
      !is.null(enc$armspec) || !is.null(enc$temporal) || !is.null(enc$re)) {
    stop("cover() NUTS is the non-spatial sampler: a spatial / temporal / re ",
         "term in the formula is not wired for method = 'nuts'. Use ",
         "method = 'nested_laplace' for structured terms.", call. = FALSE)
  }

  p_pres <- ncol(enc$occ_data$X)
  p_pos  <- ncol(enc$pos_data$X)
  n_par  <- p_pres + p_pos + 1L

  ld_name <- if (identical(positive, "beta")) "log_phi" else "log_sigma_pos"
  par_names <- c(
    paste0("presence_", colnames(enc$occ_data$X)),
    paste0("positive_", colnames(enc$pos_data$X)),
    ld_name
  )

  # Warm start at the Laplace mode + a diagonal Laplace metric from its vcov.
  # The natural-scale fit returns the per-arm mode and Hessian directly.
  warm <- fit_cover_hurdle(enc, positive = positive, engine = "laplace",
                           priors = priors,
                           control = list(max.iter = 200L))
  beta_pres0 <- warm$m_occ$mode[seq_len(p_pres)]
  beta_pos0  <- warm$m_pos$mode[seq_len(p_pos)]
  log_disp0  <- if (identical(positive, "beta")) log(warm$phi_pos)
                else log(warm$sigma_pos)
  theta0 <- c(beta_pres0, beta_pos0, log_disp0)

  V_pres <- tryCatch(solve(warm$m_occ$H_beta), error = function(e) NULL)
  pos_vscale <- if (identical(positive, "lognormal")) warm$sigma_pos^2 else 1
  V_pos  <- tryCatch(pos_vscale * solve(warm$m_pos$H_beta), error = function(e) NULL)
  diag_pres <- if (!is.null(V_pres)) pmax(diag(as.matrix(V_pres)), 1e-6)
               else rep(1, p_pres)
  diag_pos  <- if (!is.null(V_pos))  pmax(diag(as.matrix(V_pos)), 1e-6)
               else rep(1, p_pos)
  inv_metric <- c(diag_pres, diag_pos, 0.01)

  spec <- .tobs_cover_nuts_spec(enc)

  run_chain <- function(ch) {
    cpp_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta,
      sigma_logdisp = sigma.logdisp, inv_metric = inv_metric,
      n_iter = as.integer(n.iter + n.warmup), n_warmup = as.integer(n.warmup),
      max_treedepth = as.integer(max.treedepth), adapt_delta = adapt.delta,
      seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose) && ch == 1L)
  }
  n_chains <- max(1L, as.integer(n.chains))
  chains   <- lapply(seq_len(n_chains), run_chain)
  per_chain_draws <- lapply(chains, `[[`, "draws")
  draws    <- do.call(rbind, per_chain_draws)
  colnames(draws) <- par_names
  n_draws  <- nrow(draws)

  means  <- colMeans(draws); names(means) <- par_names
  V_post <- stats::cov(draws); dimnames(V_post) <- list(par_names, par_names)
  sds    <- sqrt(pmax(diag(V_post), 0)); names(sds) <- par_names

  # Per-arm posterior summaries on the natural coefficient scale, named to match
  # the design columns (so print / summary / predict.cover_fit read them).
  beta_occ <- means[seq_len(p_pres)]
  se_occ   <- sds[seq_len(p_pres)]
  beta_pos <- means[p_pres + seq_len(p_pos)]
  se_pos   <- sds[p_pres + seq_len(p_pos)]
  names(beta_occ) <- names(se_occ) <- colnames(enc$occ_data$X)
  names(beta_pos) <- names(se_pos) <- colnames(enc$pos_data$X)

  log_disp_draws <- draws[, n_par]
  if (identical(positive, "beta")) {
    phi_pos   <- mean(exp(log_disp_draws))
    sigma_pos <- NA_real_
    disp_sd   <- stats::sd(exp(log_disp_draws))
  } else {
    sigma_pos <- mean(exp(log_disp_draws))
    phi_pos   <- NA_real_
    disp_sd   <- stats::sd(exp(log_disp_draws))
  }

  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- as.integer(unlist(lapply(chains, `[[`, "divergent")))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_,
                           numeric(1)), na.rm = TRUE)

  nuts <- list(accept_prob = accept, divergent = divergent,
               treedepth = treedepth, epsilon = epsilon, n_chains = n_chains,
               divergent_total = sum(divergent), sigma_beta = sigma.beta,
               sigma_logdisp = sigma.logdisp)

  hyperpar <- list(occ = NULL, pos = NULL)
  if (identical(positive, "lognormal")) hyperpar$sigma_pos <- sigma_pos
  else                                  hyperpar$phi_pos   <- phi_pos

  # Data log-likelihood at the posterior mean (scale-invariant), so logLik() on
  # the NUTS fit matches the laplace-path convention.
  ll_mean <- sum(.tobs_cover_loglik_at_mean_nuts(enc, beta_occ, beta_pos,
                                                 if (identical(positive, "beta"))
                                                   phi_pos else sigma_pos))

  fit <- structure(list(
    occ          = warm$m_occ,
    pos          = warm$m_pos,
    beta_occ     = beta_occ,
    beta_pos     = beta_pos,
    se_occ       = se_occ,
    se_pos       = se_pos,
    positive     = positive,
    sigma_pos    = sigma_pos,
    sigma_pos_sd = if (identical(positive, "lognormal")) disp_sd else NA_real_,
    phi_pos      = phi_pos,
    phi_pos_sd   = if (identical(positive, "beta")) disp_sd else NA_real_,
    hyperpar     = hyperpar,
    encoding     = enc,
    family       = family,
    n_total      = enc$N,
    n_positive   = length(enc$idx_pos),
    converged    = TRUE,
    log_marginal = c(occ = ll_mean, pos = NA_real_),
    skew_occ     = NULL,
    skew_pos     = NULL,
    draws_occ    = NULL,
    draws_pos    = NULL,
    sla_status   = "off",
    # NUTS posterior + generic inference surface (coef / vcov / confint / logLik).
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V_post,
    n_samples    = n_draws,
    n_params     = n_par,
    # logLik.tulpa_fit reads log_prob (the per-draw vector); the data
    # log-likelihood at the posterior mean is scale-invariant, so a constant
    # vector matches the laplace-path logLik convention.
    log_prob     = rep(ll_mean, n_draws),
    log_lik      = ll_mean,
    N            = enc$N,
    method       = "nuts",
    nuts         = nuts,
    convergence  = list(converged = NA, n_iter = as.integer(n.iter))
  ), class = c("cover_fit", "tobs_multiarm_fit", "tobs_fit", "tulpa_fit"))

  # Per-parameter split-R-hat / bulk + tail ESS, through the writer every sampled
  # path shares, so summary.cover_fit surfaces them per parameter; `fit$nuts`
  # carries the same two vectors alongside the sampler diagnostics.
  fit <- .tobs_nuts_attach_convergence(fit, per_chain_draws,
                                       par_names = par_names)
  fit$nuts$rhat <- fit$convergence$rhat
  fit$nuts$ess  <- fit$convergence$ess_bulk
  fit
}


# Per-parameter posterior summary for a cover NUTS fit: posterior mean / sd and
# 2.5% / 50% / 97.5% quantiles from the draws, plus the cross-chain Rhat / ESS
# the convergence list carries (NA on a single chain). Mirrors the generic
# summary.tobs_fit table so the sampler diagnostics are visible for the cover_fit
# class (which has its own bespoke list-style summary on the Laplace path).
.tobs_cover_nuts_summary <- function(object) {
  draws <- object$draws
  nm    <- colnames(draws)
  q     <- t(apply(draws, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975),
                   names = FALSE))
  cv    <- object$convergence
  idx   <- match(nm, cv$parameter)
  data.frame(
    parameter = nm,
    mean      = colMeans(draws),
    sd        = apply(draws, 2L, stats::sd),
    `2.5%`    = q[, 1L],
    `50%`     = q[, 2L],
    `97.5%`   = q[, 3L],
    rhat      = cv$rhat[idx],
    ess_bulk  = cv$ess_bulk[idx],
    ess_tail  = cv$ess_tail[idx],
    row.names = nm,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Pointwise data log-likelihood (length N) at the posterior-mean coefficients
# for a cover NUTS fit: presence Bernoulli at every row plus, at present rows,
# the positive-part density. Mirrors .tobs_cover_hurdle_ll with one draw, but
# reads the natural-scale NUTS encoding directly so it does not depend on the
# per-arm Laplace mode + Hessian the separate-Laplace eta-draw path needs.
.tobs_cover_loglik_at_mean_nuts <- function(enc, beta_occ, beta_pos, disp) {
  X_pres <- enc$occ_data$X
  X_pos  <- enc$pos_data$X
  present <- enc$occ_data$y
  idx_pos <- enc$idx_pos
  positive <- enc$positive

  eta_pres <- as.numeric(X_pres %*% beta_occ)
  ll <- ifelse(present == 1L,
               stats::plogis(eta_pres, log.p = TRUE),
               stats::plogis(-eta_pres, log.p = TRUE))

  eta_pos <- as.numeric(X_pos %*% beta_pos)
  if (identical(positive, "beta")) {
    y_pos <- enc$pos_data$y
    mu    <- stats::plogis(eta_pos)
    dens  <- stats::dbeta(y_pos, mu * disp, (1 - mu) * disp, log = TRUE)
  } else {
    # enc$pos_data$y = log(cover); natural-scale density adds the Jacobian -log y.
    log_y <- enc$pos_data$y
    dens  <- stats::dnorm(log_y, mean = eta_pos, sd = disp, log = TRUE) - log_y
  }
  ll[idx_pos] <- ll[idx_pos] + dens
  ll
}
