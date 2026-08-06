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
                                          sigma.logdisp = 5, field = NULL) {
  pin         <- model$process_info
  p_occ       <- pin[[1L]]$p
  p_det_site  <- ncol(model$X_det_site)
  p_det_visit <- if (!is.null(model$X_det_visit)) ncol(model$X_det_visit) else 0L
  p_pos_site  <- ncol(model$X_pos_site)
  p_pos_visit <- if (!is.null(model$X_pos_visit)) ncol(model$X_pos_visit) else 0L
  p_p   <- p_det_site + p_det_visit
  p_pos <- p_pos_site + p_pos_visit
  n_coef <- p_occ + p_p + p_pos
  total  <- n_coef + 1L

  # Optional fixed-hyper coupled field (gcol33/tulpaObs#74): f = Linv %*% raw,
  # raw the trailing block of theta. f enters psi additively and the cover arm
  # scaled by alpha (both fixed at the nested-Laplace estimate). The R oracle
  # mirrors the C++ FullGradFn byte-for-byte; absent, the non-spatial target.
  has_field <- !is.null(field)
  if (has_field) {
    n_field  <- field$n_field_units
    raw      <- theta[total + seq_len(n_field)]
    f_field  <- as.numeric(field$Linv %*% raw)
    f_site   <- f_field[field$field_map]          # per-site field value
    alpha    <- field$alpha
  } else {
    f_site <- numeric(model$n_sites)
    alpha  <- 0
  }

  bo         <- theta[seq_len(p_occ)]
  bp_site    <- theta[p_occ + seq_len(p_det_site)]
  bp_visit   <- if (p_det_visit > 0L) theta[p_occ + p_det_site + seq_len(p_det_visit)] else numeric(0)
  bpos_site  <- theta[p_occ + p_p + seq_len(p_pos_site)]
  bpos_visit <- if (p_pos_visit > 0L) theta[p_occ + p_p + p_pos_site + seq_len(p_pos_visit)] else numeric(0)
  log_disp   <- theta[total]
  disp       <- exp(log_disp)
  pos_code   <- .tobs_cover_pos_code(model$positive)

  N <- model$n_sites; J <- model$max_visits
  sgm <- function(e) 1 / (1 + exp(-e))

  eta_psi <- as.numeric(model$X_occ %*% bo) + f_site
  psi     <- sgm(eta_psi)

  eta_p <- matrix(as.numeric(model$X_det_site %*% bp_site), N, J)   # site block broadcast
  if (p_det_visit > 0L) {
    eta_p <- eta_p + matrix(as.numeric(model$X_det_visit %*% bp_visit), N, J, byrow = TRUE)
  }
  eta_pos <- matrix(as.numeric(model$X_pos_site %*% bpos_site) + alpha * f_site, N, J)
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
      # Cover factor at detected visits with an observed cover; a missing (NA)
      # cover drops out (missing-at-random cover), the detection loop above still
      # counts the visit.
      for (v in vv[y[i, vv] == 1L & is.finite(y_pos[i, vv])]) {
        ev <- eta_pos[i, v]; yy <- y_pos[i, v]
        if (pos_code == 3L) {            # beta
          mu <- sgm(ev); a <- mu * disp; b <- (1 - mu) * disp
          ly <- log(yy); l1my <- log(1 - yy)
          lp_data <- lp_data + lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my
          g_eta_pos[i, v] <- disp * mu * (1 - mu) *
                             (-digamma(a) + digamma(b) + ly - l1my)
          g_logdisp <- g_logdisp + disp * (digamma(disp) - mu * digamma(a) -
                       (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my)
        } else if (pos_code == 4L) {     # identity-Gaussian (#112): raw response
          sig <- disp; r <- (yy - ev) / sig
          lp_data <- lp_data - log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        } else {                         # lognormal: Gaussian on log(cover)
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

  # Field block (gcol33/tulpaObs#74): accumulate the per-cell field score (psi
  # arm + alpha * cover arm), push it back to raw via Linv^T, and add the
  # whitened N(0, I) prior. The field gradient appends to the coefficient block.
  if (has_field) {
    g_f_site <- g_eta_psi + alpha * rowSums(g_eta_pos)
    g_field  <- numeric(n_field)
    for (s in seq_len(N)) {
      cell <- field$field_map[s]
      g_field[cell] <- g_field[cell] + g_f_site[s]
    }
    g_raw <- as.numeric(crossprod(field$Linv, g_field))   # Linv^T g_field
    grad  <- c(grad, g_raw)
  }

  lp  <- lp_data
  ib2 <- 1 / sigma.beta^2
  nb  <- n_coef
  bv  <- theta[seq_len(nb)]
  lp  <- lp - 0.5 * ib2 * sum(bv^2)
  grad[seq_len(nb)] <- grad[seq_len(nb)] - ib2 * bv
  ild2 <- 1 / sigma.logdisp^2
  lp   <- lp - 0.5 * ild2 * log_disp^2
  grad[total] <- grad[total] - ild2 * log_disp
  if (has_field) {
    raw_v <- theta[total + seq_len(n_field)]
    lp    <- lp - 0.5 * sum(raw_v^2)
    grad[total + seq_len(n_field)] <- grad[total + seq_len(n_field)] - raw_v
  }
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
    pos_code    = .tobs_cover_pos_code(model$positive),
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
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
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

  fit <- structure(list(
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
    convergence  = list(converged = NA, n_iter = as.integer(n.iter))
  ), class = c("tobs_fit", "tulpa_fit"))

  # Per-parameter split-R-hat / bulk + tail ESS, through the writer every sampled
  # path shares, so summary.tobs_fit surfaces them per parameter; `fit$nuts`
  # carries the same two vectors alongside the sampler diagnostics.
  fit <- .tobs_nuts_attach_convergence(fit, per_chain_draws,
                                       par_names = par_names)
  fit$nuts$rhat <- fit$convergence$rhat
  fit$nuts$ess  <- fit$convergence$ess_bulk
  fit
}


# ---------------------------------------------------------------------------
# Spatial occu_cover NUTS: fixed-hyper non-centered coupled proper-CAR field
# (gcol33/tulpaObs#74)
# ---------------------------------------------------------------------------

# Resolve the single spatial term + fixed-effect psi formula for the NUTS
# spatial path. Unlike .occu_cover_spatial_fields (which gates to icar/bym2 for
# the grid-integrated nested-Laplace engine), the NUTS path accepts car_proper()
# (the full-rank precision the fixed-hyper non-centered field needs) and rejects
# icar/bym2/SVC/trend/temporal/RE/correlated bars with a pointer to the
# nested-Laplace route. Returns NULL when the psi formula carries no spatial term
# (the non-spatial NUTS sampler), or list(fe, spatial, group_var).
.occu_cover_nuts_spatial_term <- function(formula, data) {
  bind <- .tobs_bind_formulas(list(psi = formula), data)
  if (length(bind$terms) == 0L) return(NULL)
  spatial <- Filter(function(t) inherits(t$spec, "tobs_spatial"), bind$terms)
  bars    <- Filter(function(t) isTRUE(t$spec$is_bar), spatial)
  re_terms <- Filter(function(t) inherits(t$spec, "tobs_re"), bind$terms)
  other   <- Filter(function(t) !inherits(t$spec, "tobs_spatial") &&
                                !inherits(t$spec, "tobs_re"), bind$terms)
  if (length(spatial) == 0L) return(NULL)
  if (length(spatial) > 1L || length(bars) > 0L || length(re_terms) > 0L ||
      length(other) > 0L) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples a SINGLE shared coupled field ",
      "(one car_proper() term on the psi formula); SVC / trend / correlated bars ",
      "/ per-group RE / temporal terms compose only on the grid-integrated ",
      "method = \"nested_laplace\" path."), call. = FALSE)
  }
  spec <- spatial[[1L]]$spec
  if (!is.null(spec$weight) || isTRUE(spec$is_multifield)) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples a single unweighted shared ",
      "field; a weighted SVC field needs method = \"nested_laplace\"."),
      call. = FALSE)
  }
  list(fe = bind$fe$psi, spatial = spec, group_var = spec$group_var)
}


# Run the nested-Laplace joint engine once with a proper-CAR copy block
# to obtain the FIXED field hyperparameters (sigma, rho_car) and the copy
# amplitude (alpha), plus the grid-weighted posterior-mean coupled field. This is
# the occu_cover analogue of the warm fit .tobs_fit_abun_nuts_spatial reads from
# nmix_laplace_car_proper: the same proper-CAR precision tau Q(rho) the NUTS field
# block then fixes, here estimated through the shared two-state cell-coupling
# marginal (the occu_cover_{lognormal,beta} spec) so the field hyper sits at the
# nested-Laplace estimate of THIS model (the marginalize-then-fix recipe). Returns
# the betas, log-dispersion, (sigma, rho, alpha), and field f at its grid-weighted
# posterior mean.
.tobs_occu_cover_nuts_carproper_warm <- function(model, adj, priors,
                                                 max.iter = 200L, tol = 1e-6,
                                                 type = "car_proper",
                                                 sigma.grid = NULL,
                                                 rho.car.grid = NULL,
                                                 alpha.grid = NULL) {
  is_beta   <- identical(model$positive, "beta")
  spec_name <- switch(model$positive,
                      beta     = "occu_cover_beta",
                      gaussian = "occu_cover_gaussian",
                      "occu_cover_lognormal")
  n_sites   <- model$n_sites
  site_cell <- model$site_cell %||% seq_len(n_sites)
  n_cells   <- nrow(adj)

  # Pre-fit the pos-arm dispersion at the empirical cover spread (matching the
  # joint non-latent path); it rides the spec's phi slot, fixed here.
  # Observed covers only (a detected visit may carry a missing cover).
  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  pos_vals <- pos_vals[is.finite(pos_vals)]
  sigma_pos_init <- if (is_beta) {
    if (length(pos_vals) >= 2L) {
      mu_hat <- mean(pos_vals); var_hat <- max(stats::var(pos_vals), 1e-6)
      max((mu_hat * (1 - mu_hat)) / var_hat - 1, 1)
    } else 10
  } else if (identical(model$positive, "gaussian")) {
    # Residual SD on the raw response scale (no log; the response may be negative).
    if (length(pos_vals) >= 2L) max(stats::sd(pos_vals), 0.05) else 1
  } else {
    if (length(pos_vals) > 0L) max(stats::sd(log(pos_vals)), 0.05) + 0.05 else 0.4
  }

  alpha_grid <- alpha.grid %||% .tobs_default_alpha_grid()
  sigma_grid <- sigma.grid %||% .tobs_default_sigma_grid()
  rho_car_grid <- rho.car.grid %||% c(0.5, 0.8, 0.95, 0.99)

  # Single-block (multi = FALSE): the pos arm carries
  # field_coef = list(name = "alpha", grid = alpha_grid), so the copy alpha axis
  # rides the one shared field. The single-block joint path takes the copy
  # coefficient on the arm, not a top-level `copy` block.
  arms_out <- .occu_cover_build_joint_arms(
    model = model, sigma_pos_init = sigma_pos_init, alpha_grid = alpha_grid,
    positive = model$positive, multi = FALSE, n_cells = n_cells,
    site_cell = site_cell, cover_aggregate = "none")
  responses <- arms_out$responses

  arm_priors <- .occu_cover_coupled_arm_priors(priors, responses)
  for (nm in c("psi", "p", "pos")) {
    ap <- arm_priors[[nm]]
    if (!is.null(ap)) {
      responses[[nm]]$beta_prior_mean <- ap$mean
      responses[[nm]]$beta_prior_prec <- ap$prec
    }
  }

  csr <- .occu_cover_adj_to_csr(adj)
  prior_arg <- list(
    type            = type,
    n_spatial_units = csr$n_spatial_units,
    adj_row_ptr     = csr$adj_row_ptr,
    adj_col_idx     = csr$adj_col_idx,
    n_neighbors     = csr$n_neighbors,
    sigma_grid      = sigma_grid,
    spatial_idx     = lapply(responses, function(a) as.integer(a$spatial_idx)))
  # Only the proper-CAR field carries a spatial-correlation rho grid; the
  # intrinsic icar / bym2 fields fix rho (icar rho = 1, bym2 mixing gridded by the
  # engine's own bym2 axis).
  if (identical(type, "car_proper")) prior_arg$rho_car_grid <- rho_car_grid

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = responses, prior = prior_arg,
    cell_coupling = spec_name,
    control = list(max_iter = as.integer(max.iter), tol = as.numeric(tol),
                   n_threads = 1L, store_Q = FALSE, adaptive_grid = FALSE,
                   var_of_means_consistency = FALSE, diagnose_k = FALSE,
                   progress = FALSE))

  ok <- which(is.finite(fit$log_marginal))
  if (length(ok) == 0L)
    stop("occu_cover NUTS spatial: warm car_proper fit failed at every grid ",
         "cell. Bump control$max.iter or tighten control$tol.", call. = FALSE)
  w_raw <- exp(fit$log_marginal[ok] - max(fit$log_marginal[ok]))
  w     <- w_raw / sum(w_raw)

  layout <- fit$arm_layout
  p_psi  <- layout$p[1L]; p_p <- layout$p[2L]; p_pos <- layout$p[3L]
  bpsi_idx <- layout$beta_start[1L] + seq_len(p_psi)
  bp_idx   <- layout$beta_start[2L] + seq_len(p_p)
  bpos_idx <- layout$beta_start[3L] + seq_len(p_pos)
  f0       <- (layout$field_starts %||% layout$phi_start)
  field_idx <- f0[[1L]] + seq_len(n_cells)
  modes <- fit$modes[ok, , drop = FALSE]

  beta_psi <- as.numeric(crossprod(w, modes[, bpsi_idx, drop = FALSE]))
  beta_p   <- as.numeric(crossprod(w, modes[, bp_idx,   drop = FALSE]))
  beta_pos <- as.numeric(crossprod(w, modes[, bpos_idx, drop = FALSE]))
  field    <- as.numeric(crossprod(w, modes[, field_idx, drop = FALSE]))
  field    <- field - mean(field)            # sum-to-zero convention

  tg     <- fit$theta_grid[ok, , drop = FALSE]
  pick   <- function(nm) {
    j <- match(nm, colnames(fit$theta_grid)); if (is.na(j)) return(NA_real_)
    sum(w * as.numeric(tg[, j]))
  }
  sigma <- pick("sigma"); alpha <- pick("alpha")
  rho   <- if (identical(type, "car_proper")) pick("rho_car")
           else if (identical(type, "bym2")) pick("rho") else 1.0

  list(beta_psi = beta_psi, beta_p = beta_p, beta_pos = beta_pos,
       log_disp = log(sigma_pos_init), field = field, type = type,
       sigma = sigma, rho = rho, alpha = alpha, joint_fit = fit)
}


# Sample the exact occu_cover coefficient posterior jointly with a FIXED-HYPER
# non-centered coupled PROPER-CAR field on the latent state z (gcol33/tulpaObs#74):
# the psi-arm field f (one value per cell) enters psi linearly and is copied to the
# cover (positive) arm with the FIXED scaling alpha; the field precision tau Q(rho)
# and alpha are fixed at the nested-Laplace joint estimate, and the
# whitened raw ~ N(0, I) (f = Linv %*% raw) is sampled alongside the coefficient
# marginal. Parameter vector: c(beta_psi, beta_p, beta_pos, log_disp, raw_field).
# This is the occu_cover analogue of .tobs_fit_abun_nuts_spatial (tulpa#87): a
# full-rank proper-CAR precision gives a well-conditioned non-centered geometry, so
# the field block reuses the abun#51 field machinery (the optional field block in
# src/occu_cover_nuts.cpp, byte-exact vs the R oracle's field branch). Intrinsic
# icar / bym2 fields also sample here via the sum-to-zero eigen-loading that drops
# the precision null-space (constant) direction (gcol33/tulpaObs#71/#113), the same
# reparam abun NUTS+areal uses. `...` absorbs unused sampler controls.
.tobs_fit_occu_cover_nuts_spatial <- function(model, spatial, priors = NULL,
                                              sigma.beta = 5, sigma.logdisp = 5,
                                              n.iter = 2000L, n.warmup = 1000L,
                                              n.chains = 1L, max.treedepth = 10L,
                                              adapt.delta = 0.9, seed = 1L,
                                              max.iter = 200L, tol = 1e-6,
                                              verbose = FALSE,
                                              sigma.grid = NULL, rho.car.grid = NULL,
                                              alpha.grid = NULL, ...) {
  .tobs_reject_weighted_spatial(spatial, "occu_cover NUTS psi spatial")
  if (!spatial$type %in% c("icar", "car_proper", "bym2")) {
    stop(sprintf(paste0(
      "occu_cover() NUTS + areal spatial supports icar() / car_proper() / ",
      "bym2() on the psi arm (coupled to the cover arm); got '%s'. Use ",
      "method = \"nested_laplace\" for other field kinds. (gcol33/tulpaObs#74, #113)"),
      spatial$type), call. = FALSE)
  }
  pin   <- model$process_info
  p_occ <- pin[[1L]]$p; p_p <- pin[[2L]]$p; p_pos <- pin[[3L]]$p
  n_coef <- p_occ + p_p + p_pos
  n_base <- n_coef + 1L                    # + log_dispersion

  n_sites   <- model$n_sites
  site_cell <- model$site_cell %||% seq_len(n_sites)
  adj       <- as.matrix(spatial$graph)
  n_cells   <- nrow(adj)
  if (length(site_cell) != n_sites || max(site_cell) > n_cells ||
      min(site_cell) < 1L) {
    stop(sprintf(paste0(
      "occu_cover NUTS spatial: site_cell must map %d sites into 1..%d ",
      "graph nodes."), n_sites, n_cells), call. = FALSE)
  }

  # Fixed hyper (sigma -> tau, rho_car) + copy amplitude alpha + warm betas / field
  # from the nested-Laplace joint proper-CAR fit.
  warm <- .tobs_occu_cover_nuts_carproper_warm(
    model, adj, priors, type = spatial$type, max.iter = max.iter, tol = tol,
    sigma.grid = sigma.grid, rho.car.grid = rho.car.grid, alpha.grid = alpha.grid)
  alpha <- warm$alpha

  # Whitened coupled-field loading L (f = L %*% raw): the square inverse Cholesky
  # of the fixed precision for car_proper, the sum-to-zero eigen-loading for the
  # intrinsic icar / bym2 fields (gcol33/tulpaObs#71/#113). The joint warm fit
  # reports the field marginal SD sigma; tau = 1 / sigma^2 pins the icar/car
  # precision.
  fl <- .tobs_nuts_field_loading(adj, spatial$type, n_cells,
                                 tau = 1 / max(warm$sigma, 1e-3)^2, rho = warm$rho,
                                 sigma = warm$sigma,
                                 scale_factor = spatial$scale_factor)
  field_load <- fl$field_load; n_raw <- fl$n_raw

  # car_proper warm-starts raw near the integrated field (raw0 = L %*% f_warm);
  # the non-square icar / bym2 loadings start raw at 0.
  raw0 <- if (identical(spatial$type, "car_proper")) {
    Q <- .areal_Q(adj, fl$rho)
    L <- tryCatch(chol(fl$tau * Q + diag(1e-4 * fl$tau, n_cells)),
                  error = function(e) NULL)
    if (is.null(L)) stop("occu_cover NUTS spatial: field precision not PD.",
                         call. = FALSE)
    as.numeric(L %*% warm$field)
  } else numeric(n_raw)
  theta0 <- c(warm$beta_psi, warm$beta_p, warm$beta_pos, warm$log_disp, raw0)

  spec <- .tobs_occu_cover_nuts_spec(model)
  spec$n_field_units <- n_cells
  spec$field_map     <- as.integer(site_cell)
  spec$field_load    <- field_load
  spec$field_alpha   <- as.numeric(alpha)

  inv_metric <- c(rep(0.1, n_base), rep(1, n_cells))

  par_names <- c(
    paste0("psi_", pin[[1L]]$coef_names),
    paste0("p_",   pin[[2L]]$coef_names),
    paste0("pos_", pin[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos")
  raw_names <- paste0("raw_", seq_len(n_raw))
  all_names <- c(par_names, raw_names)

  run_chain <- function(ch) {
    cpp_occu_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose) && ch == 1L)
  }
  n_chains <- max(1L, as.integer(n.chains))
  chains   <- lapply(seq_len(n_chains), run_chain)
  per_chain_draws <- lapply(chains, `[[`, "draws")
  draws    <- do.call(rbind, per_chain_draws)
  colnames(draws) <- all_names
  n_draws  <- nrow(draws)

  b_idx <- seq_len(n_base)
  means <- colMeans(draws)
  V_post <- stats::cov(draws[, b_idx, drop = FALSE])
  dimnames(V_post) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V_post), 0)); names(sds) <- par_names
  par_means <- means[b_idx]; names(par_means) <- par_names

  # Posterior-mean coupled field f = L %*% mean(raw) (n_raw whitened coords).
  raw_idx <- n_base + seq_len(n_raw)
  field_mean <- as.numeric(field_load %*% colMeans(draws[, raw_idx, drop = FALSE]))

  # Data log-likelihood at the posterior mean (field on psi + alpha*field on
  # cover), so logLik() matches the laplace convention.
  bo   <- par_means[seq_len(p_occ)]
  bp   <- par_means[p_occ + seq_len(p_p)]
  bpos <- par_means[p_occ + p_p + seq_len(p_pos)]
  eta  <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  f_site <- field_mean[site_cell]
  psi_f  <- stats::plogis(.tobs_clamp_eta(stats::qlogis(eta$psi) + f_site))
  ep_f   <- eta$ep_mat + alpha * f_site
  ll_mean <- sum(.occu_cover_site_ll(model, psi_f, eta$p_mat, ep_f, par_means[n_base]))

  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- as.integer(unlist(lapply(chains, `[[`, "divergent")))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_, numeric(1)),
                    na.rm = TRUE)

  nuts <- list(accept_prob = accept, divergent = divergent, treedepth = treedepth,
               epsilon = epsilon, n_chains = n_chains,
               divergent_total = sum(divergent),
               sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
               field_tau = fl$tau, field_rho = fl$rho, field_alpha = alpha,
               field_sigma = warm$sigma, fixed_hyper = TRUE,
               n_field_units = n_cells)

  fit <- structure(list(
    draws        = draws[, b_idx, drop = FALSE],
    means        = par_means,
    sds          = sds,
    vcov         = V_post,
    n_samples    = n_draws,
    n_params     = n_base,
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
    spatial      = list(type = spatial$type, graph = adj,
                        sigma_mean = warm$sigma, alpha_mean = alpha,
                        rho_mean = fl$rho),
    spatial_field = field_mean,
    method       = "nuts",
    positive     = model$positive,
    nuts         = nuts,
    convergence  = list(converged = NA, n_iter = as.integer(n.iter))
  ), class = c("tobs_fit", "tulpa_fit"))

  # Diagnostics over the coefficient block (the coordinates the fit reports);
  # the whitened field `raw` coordinates carry no named parameter.
  fit <- .tobs_nuts_attach_convergence(fit, per_chain_draws, par_names = par_names,
                                       cols = b_idx)
  fit$nuts$rhat <- fit$convergence$rhat
  fit$nuts$ess  <- fit$convergence$ess_bulk
  fit
}
