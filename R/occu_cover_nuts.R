# occu_cover_nuts.R - NUTS target for the non-spatial joint occupancy +
# cover-hurdle family (occu_cover()).
#
# The Laplace fit (.tobs_fit_occu_cover) sums the latent occupancy state z out in
# closed form (two states per cell) and returns a Gaussian observed-Fisher
# posterior over the packed coefficient vector
#   theta = c(beta_psi, beta_p, beta_pos, log_dispersion).
# NUTS instead samples the exact marginal posterior of that vector, giving
# calibrated (non-Gaussian) intervals and the per-draw pointwise likelihood that
# WAIC / LOO need. There is no latent field, random effect, or community
# covariance on the non-spatial path, so the parameter vector is just the flat
# three-arm coefficient block plus one log-dispersion scalar.
#
# .tobs_occu_cover_nuts_logpost is the R oracle: it recomputes the joint
# log-posterior and gradient exactly as the C++ FullGradFn (src/occu_cover_nuts.cpp,
# cpp_occu_cover_nuts) does, and a byte-exact test cross-checks the two before the
# sampler is trusted. The per-cell marginal mirrors .occu_cover_site_ll; the cover
# arm and the no-detection mixture reuse the same closed forms the Laplace path uses.


# Joint log-posterior + gradient of the non-spatial occu_cover coefficient vector
# theta = c(beta_psi, beta_p, beta_pos, log_disp) (beta_p / beta_pos each packed as
# the site-level block then the optional visit-level block, exactly as the fitter
# stacks them). Weak Gaussian priors N(0, sigma.beta^2) on every coefficient and a
# broad N(0, sigma.logdisp^2) on log_disp keep the ridge / dispersion proper
# without materially shifting the data-dominated optimum. Returns list(lp, grad)
# over the packed coordinates. This is the oracle the C++ FullGradFn mirrors.
.tobs_occu_cover_nuts_logpost <- function(theta, model, sigma.beta = 5,
                                          sigma.logdisp = 5) {
  pin         <- model$process_info
  p_occ       <- pin[[1L]]$p
  p_det_site  <- ncol(model$X_det_site)
  p_det_visit <- if (!is.null(model$X_det_visit)) ncol(model$X_det_visit) else 0L
  p_pos_site  <- ncol(model$X_pos_site)
  p_pos_visit <- if (!is.null(model$X_pos_visit)) ncol(model$X_pos_visit) else 0L
  p_p   <- p_det_site + p_det_visit
  p_pos <- p_pos_site + p_pos_visit
  total <- p_occ + p_p + p_pos + 1L

  bo         <- theta[seq_len(p_occ)]
  bp_site    <- theta[p_occ + seq_len(p_det_site)]
  bp_visit   <- if (p_det_visit > 0L) theta[p_occ + p_det_site + seq_len(p_det_visit)] else numeric(0)
  bpos_site  <- theta[p_occ + p_p + seq_len(p_pos_site)]
  bpos_visit <- if (p_pos_visit > 0L) theta[p_occ + p_p + p_pos_site + seq_len(p_pos_visit)] else numeric(0)
  log_disp   <- theta[total]
  disp       <- exp(log_disp)
  is_beta    <- identical(model$positive, "beta")

  N <- model$n_sites; J <- model$max_visits
  sgm <- function(e) 1 / (1 + exp(-e))

  eta_psi <- as.numeric(model$X_occ %*% bo)
  psi     <- sgm(eta_psi)

  eta_p <- matrix(as.numeric(model$X_det_site %*% bp_site), N, J)   # site block broadcast
  if (p_det_visit > 0L) {
    eta_p <- eta_p + matrix(as.numeric(model$X_det_visit %*% bp_visit), N, J, byrow = TRUE)
  }
  eta_pos <- matrix(as.numeric(model$X_pos_site %*% bpos_site), N, J)
  if (p_pos_visit > 0L) {
    eta_pos <- eta_pos + matrix(as.numeric(model$X_pos_visit %*% bpos_visit), N, J, byrow = TRUE)
  }

  valid <- model$valid; y <- model$y; y_pos <- model$y_pos

  g_eta_psi <- numeric(N)
  g_eta_p   <- matrix(0, N, J)
  g_eta_pos <- matrix(0, N, J)
  g_logdisp <- 0
  lp_data   <- 0

  for (i in seq_len(N)) {
    vv <- which(valid[i, ])
    if (length(vv) == 0L) next
    any_det <- any(y[i, vv] == 1L)
    if (any_det) {
      lp_data <- lp_data + log(psi[i])
      g_eta_psi[i] <- 1 - psi[i]
      for (v in vv) {
        pv <- sgm(eta_p[i, v])
        if (y[i, v] == 1L) { lp_data <- lp_data + log(pv);     g_eta_p[i, v] <- 1 - pv }
        else               { lp_data <- lp_data + log(1 - pv); g_eta_p[i, v] <- -pv }
      }
      for (v in vv[y[i, vv] == 1L]) {
        ev <- eta_pos[i, v]; yy <- y_pos[i, v]
        if (is_beta) {
          mu <- sgm(ev); a <- mu * disp; b <- (1 - mu) * disp
          ly <- log(yy); l1my <- log(1 - yy)
          lp_data <- lp_data + lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my
          g_eta_pos[i, v] <- disp * mu * (1 - mu) *
                             (-digamma(a) + digamma(b) + ly - l1my)
          g_logdisp <- g_logdisp + disp * (digamma(disp) - mu * digamma(a) -
                       (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my)
        } else {
          sig <- disp; r <- (log(yy) - ev) / sig
          lp_data <- lp_data - log(yy) - log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        }
      }
    } else {
      P0 <- exp(sum(log(1 - sgm(eta_p[i, vv]))))
      L  <- psi[i] * P0 + (1 - psi[i])
      lp_data <- lp_data + log(L)
      g_eta_psi[i] <- -psi[i] * (1 - psi[i]) * (1 - P0) / L
      for (v in vv) g_eta_p[i, v] <- -psi[i] * P0 * sgm(eta_p[i, v]) / L
    }
  }

  # Design-sandwich the eta-gradients onto the coefficient blocks. The visit-level
  # long vector follows site-major order (row i visit v -> (i-1)*J + v), matching
  # the byrow = TRUE design broadcast and the C++ row map i * J + v.
  grad_bo  <- as.numeric(crossprod(model$X_occ, g_eta_psi))
  grad_bp  <- as.numeric(crossprod(model$X_det_site, rowSums(g_eta_p)))
  if (p_det_visit > 0L)
    grad_bp <- c(grad_bp, as.numeric(crossprod(model$X_det_visit, as.numeric(t(g_eta_p)))))
  grad_bpos <- as.numeric(crossprod(model$X_pos_site, rowSums(g_eta_pos)))
  if (p_pos_visit > 0L)
    grad_bpos <- c(grad_bpos, as.numeric(crossprod(model$X_pos_visit, as.numeric(t(g_eta_pos)))))
  grad <- c(grad_bo, grad_bp, grad_bpos, g_logdisp)

  lp  <- lp_data
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

# Build the C++ NUTS spec list from a bound non-spatial occu_cover model. NULL
# visit-level designs become explicit 0-column matrices so the C++ side reads a
# uniform layout. The matrices are passed straight through (raw design scale; the
# occu_cover path is not autoscaled), so the draws land on the natural scale.
.tobs_occu_cover_nuts_spec <- function(model) {
  N <- model$n_sites; J <- model$max_visits
  empty_visit <- function(X) if (is.null(X)) matrix(0, N * J, 0L) else X
  list(
    n_sites     = as.integer(N),
    max_visits  = as.integer(J),
    is_beta     = identical(model$positive, "beta"),
    y           = matrix(as.integer(model$y), N, J),
    y_pos       = matrix(as.numeric(model$y_pos), N, J),
    valid       = matrix(as.integer(model$valid), N, J),
    X_occ       = model$X_occ,
    X_det_site  = model$X_det_site,
    X_det_visit = empty_visit(model$X_det_visit),
    X_pos_site  = model$X_pos_site,
    X_pos_visit = empty_visit(model$X_pos_visit)
  )
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the non-spatial joint occupancy + cover hurdle
# ---------------------------------------------------------------------------

# Sample the exact non-spatial occu_cover coefficient posterior via tulpa's NUTS
# engine and the in-tree C++ FullGradFn (cpp_occu_cover_nuts), warm-started at the
# Laplace mode with a diagonal Laplace metric, then package the draws into the
# same tobs_fit shape .tobs_fit_occu_cover returns so coef / vcov / confint /
# predict / WAIC read the NUTS posterior. The occu_cover path is not autoscaled,
# so the returned draws / means / vcov are already on the natural coefficient
# scale. `sigma.logdisp` is an internal weak-prior width (no control knob, like
# abun's sigma.logr). `...` absorbs unused sampler controls (n.thin / n.threads /
# progress.*).
.tobs_fit_occu_cover_nuts <- function(model, priors = NULL,
                                      sigma.beta = 5, sigma.logdisp = 5,
                                      n.iter = 2000L, n.warmup = 1000L,
                                      n.chains = 1L, max.treedepth = 10L,
                                      adapt.delta = 0.9, seed = 1L,
                                      verbose = FALSE, ...) {
  pin   <- model$process_info
  p_occ <- pin[[1L]]$p; p_p <- pin[[2L]]$p; p_pos <- pin[[3L]]$p
  n_par <- p_occ + p_p + p_pos + 1L

  par_names <- c(
    paste0("psi_", pin[[1L]]$coef_names),
    paste0("p_",   pin[[2L]]$coef_names),
    paste0("pos_", pin[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos"
  )

  # Warm start at the Laplace mode + a diagonal Laplace metric from its vcov.
  warm <- .tobs_fit_occu_cover(model, method = "laplace", priors = priors,
                               sigma.beta = sigma.beta, verbose = FALSE)
  theta0 <- as.numeric(warm$means)
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == n_par) && all(is.finite(diag(V))))
                  pmax(diag(V), 1e-6) else rep(1, n_par)

  spec <- .tobs_occu_cover_nuts_spec(model)

  run_chain <- function(ch) {
    cpp_occu_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
      inv_metric = inv_metric, n_iter = as.integer(n.iter),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
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

  # Data log-likelihood at the posterior mean (scale-invariant), so logLik() on
  # the NUTS fit matches the laplace-path convention. Reuses the shared marginal.
  bo   <- means[seq_len(p_occ)]
  bp   <- means[p_occ + seq_len(p_p)]
  bpos <- means[p_occ + p_p + seq_len(p_pos)]
  eta  <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  ll_mean <- sum(.occu_cover_site_ll(model, eta$psi, eta$p_mat, eta$ep_mat, means[n_par]))

  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- as.integer(unlist(lapply(chains, `[[`, "divergent")))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_, numeric(1)),
                    na.rm = TRUE)

  nuts <- list(accept_prob = accept, divergent = divergent, treedepth = treedepth,
               epsilon = epsilon, n_chains = n_chains, divergent_total = sum(divergent),
               sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp)
  # Split-R-hat / bulk-ESS when more than one chain (the shared diagnostic). The
  # convergence list mirrors the cpp_occu_fit NUTS shape (parameter / rhat /
  # ess_bulk / ess_tail) so summary.tobs_fit surfaces them per parameter.
  rhat <- ess <- rep(NA_real_, n_par)
  if (n_chains > 1L) {
    re <- .tobs_nuts_rhat_ess(per_chain_draws)
    rhat <- re$rhat; ess <- re$ess
    names(rhat) <- par_names; names(ess) <- par_names
    nuts$rhat <- rhat; nuts$ess <- ess
  }

  structure(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V_post,
    n_samples    = n_draws,
    n_params     = n_par,
    log_prob     = rep(ll_mean, n_draws),
    log_lik      = ll_mean,
    N            = sum(model$valid),
    accept_prob  = accept,
    divergent    = divergent,
    treedepth    = treedepth,
    epsilon      = epsilon,
    col_names    = par_names,
    param_names  = par_names,
    process_info = pin,
    model        = model,
    spatial      = NULL,
    method       = "nuts",
    positive     = model$positive,
    nuts         = nuts,
    convergence  = list(converged = TRUE, n_iter = as.integer(n.iter),
                        parameter = par_names, rhat = rhat,
                        ess_bulk = ess, ess_tail = rep(NA_real_, n_par))
  ), class = c("tobs_fit", "tulpa_fit"))
}
