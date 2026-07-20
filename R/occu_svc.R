# =============================================================================
# occu_svc.R - continuous NNGP spatially-varying coefficients on the Laplace
# backends of single-season occupancy (gcol33/tulpaObs#143).
#
# `spatial(lon, lat, model = "svc", coefficients = )` -- equivalently the direct
# `svc(lon, lat, coefficients = )` -- declares K continuous coefficient surfaces:
#
#   logit psi_i = X_i . beta + sum_k X[i, col_k] * z_k(s_i)
#   z_k ~ NNGP(sigma_k^2, phi_k)      (Vecchia / Datta et al. 2016)
#
# The NUTS path samples this in the compiled tulpa engine. On the Laplace
# backends it is a latent field block like any other, so it rides the SHARED
# areal-BFGS nested-Laplace driver (`.tobs_areal_bfgs_fit`, R/areal_bfgs.R):
# per outer-grid cell, BFGS over (fixed coefficients, field values) against the
# family's analytic gradient plus the field log-prior, with the Laplace marginal
# from a finite-difference Hessian at the mode. The only new pieces are
#
#   * `.tobs_svc_nngp_field()` -- the field spec (a continuous NNGP block with an
#     optional per-site design weight), the continuous sibling of
#     `.areal_field_car(weight = )`; and
#   * `.tobs_occu_svc_marginal()` -- the single-season occupancy eval closure
#     (exact two-state marginal, Fisher-identity gradients).
#
# The Gaussian-Markov precision of the Vecchia approximation is assembled here
# from the SAME neighbour structure, kernels, jitter and conditional-variance
# floor the compiled kernel uses (tulpa src/hmc_svc.h + src/nngp_cond.h), so the
# density this path integrates is the density the sampler samples;
# `test-occu-svc-laplace-recovery.R` asserts the two agree numerically.
#
# The field hyperparameters (marginal SD, range) are integrated on the outer
# grid, so `method = "laplace"` and `method = "nested_laplace"` land on the same
# fitter; `fit$svc_hyper` reports each surface's marginalised SD and range.
# =============================================================================


# Resolve an svc() term's coefficient selector to design column POSITIONS in
# `X_arm`, the design of the arm the surfaces load on. Names (`coefficients = `)
# are matched against the design column names; positions (`indices = `) are
# range-checked. Both consumers -- the Laplace field builder
# (`.tobs_svc_field_blocks`) and the NUTS spec packer (R/occu_fit.R) -- go
# through here, so the term reads its columns against the same design the
# coefficients are and neither path can drift (gcol33/tulpaObs#146).
#
# `family` and `arm` name the model in the error, so a mistyped coefficient
# reports the design it was matched against rather than a bare index failure.
.tobs_svc_columns <- function(svc, X_arm, family, arm = NULL) {
  if (is.null(arm)) arm <- .tobs_svc_arm_label(family)
  nms <- colnames(X_arm)
  if (!is.null(svc$coefficients)) {
    if (is.null(nms)) {
      stop(sprintf(paste0(
        "svc(coefficients = ): the %s() %s design has no column names to match ",
        "against; select the columns by position with `indices = `."),
        family, arm), call. = FALSE)
    }
    hit <- match(svc$coefficients, nms)
    if (anyNA(hit)) {
      stop(sprintf(paste0(
        "svc(coefficients = ): no %s() %s design column named %s.\n",
        "  Available: %s."),
        family, arm,
        paste0("`", svc$coefficients[is.na(hit)], "`", collapse = ", "),
        paste0("`", nms, "`", collapse = ", ")), call. = FALSE)
    }
    return(as.integer(hit))
  }
  idx <- as.integer(svc$indices)
  if (any(idx < 1L) || any(idx > ncol(X_arm))) {
    stop(sprintf(paste0(
      "svc(indices = ): index out of range -- the %s() %s design has %d ",
      "column(s)%s."),
      family, arm, ncol(X_arm),
      if (is.null(nms)) "" else
        paste0(" (", paste0("`", nms, "`", collapse = ", "), ")")),
      call. = FALSE)
  }
  idx
}


# Covariance kernels of the NNGP field, keyed by tulpa's CovType code
# (0 exponential, 1 Matern-3/2, 2 Gaussian; see tulpa inst/include/tulpa/types.h
# and the compute_cov dispatch in tulpa src/hmc_svc.h). `d` may be a scalar,
# vector or matrix; the result keeps its shape. Every kernel parameterises
# distance as d / phi, so phi IS the range -- which is what makes the PC range
# prior apply to it directly.
.tobs_nngp_cov <- function(d, sigma2, phi, cov_code) {
  if (cov_code == 1L) {                       # Matern, nu = 3/2
    r <- sqrt(3) * d / phi
    return(sigma2 * (1 + r) * exp(-r))
  }
  if (cov_code == 2L) {                       # Gaussian (squared exponential)
    r <- d / phi
    return(sigma2 * exp(-r * r))
  }
  sigma2 * exp(-d / phi)                      # exponential (default)
}

# Map the `cov` string an svc() term carries to tulpa's CovType code.
.tobs_nngp_cov_code <- function(cov_type) {
  switch(as.character(cov_type %||% "exponential"),
         exponential = 0L, matern = 1L, gaussian = 2L,
         stop(sprintf("svc(): unsupported covariance '%s'.", cov_type),
              call. = FALSE))
}

# Gaussian-Markov precision of the NNGP (Vecchia) approximation at one
# (sigma2, phi), plus its log-determinant.
#
# In the neighbour ordering the field factorises as w_(1) ~ N(0, sigma2) and
# w_(i) | w_(N(i)) ~ N(a_i' w_(N(i)), d_i) with kriging weights a_i =
# C(N(i), N(i))^-1 c(i, N(i)) and conditional variance d_i = sigma2 - c' a. That
# is exactly Q = (I - A)' D^-1 (I - A) with log|Q| = -sum log d_i, so the
# quadratic form below reproduces the compiled kernel's sum of conditional
# log-densities term for term -- including its jitter on the neighbour
# covariance diagonal and its blended conditional-variance floor, which are
# arguments there rather than defaults precisely so a second consumer can match
# them (tulpa src/nngp_cond.h).
#
# `co` is the n x 2 coordinate matrix in the ORIGINAL row order; `nn_idx` /
# `nn_dist` are the n x nn neighbour structure indexed by ORDERING POSITION
# (entries of `nn_idx` are 1-based ordering positions of earlier locations, 0
# where a location has fewer than nn predecessors); `ord` maps ordering position
# to original row (`compute_nngp_neighbors()`, R/spatial.R). The returned
# precision is in the original row order. NULL when a neighbour covariance could
# not be factorised at this hyperparameter cell.
.tobs_nngp_precision <- function(co, nn_idx, nn_dist, ord, sigma2, phi, cov_code,
                                 jitter = 1e-4, var_floor = 1e-4) {
  n <- nrow(co)
  co_ord <- co[ord, , drop = FALSE]
  A <- matrix(0, n, n)
  D <- rep(sigma2, n)
  for (i in seq_len(n)[-1L]) {
    nb <- nn_idx[i, ]
    nb <- nb[nb > 0L]
    m  <- length(nb)
    if (m == 0L) next
    cvec <- .tobs_nngp_cov(nn_dist[i, seq_len(m)], sigma2, phi, cov_code)
    C <- if (m == 1L) matrix(sigma2, 1L, 1L) else
      .tobs_nngp_cov(as.matrix(stats::dist(co_ord[nb, , drop = FALSE])),
                     sigma2, phi, cov_code)
    diag(C) <- sigma2 + jitter
    a <- tryCatch(solve(C, cvec), error = function(e) NULL)
    if (is.null(a) || any(!is.finite(a))) return(NULL)
    v <- sigma2 - sum(cvec * a)
    # The compiled kernel BLENDS at the floor rather than clamping, so a mode
    # that grazes it keeps a little gradient instead of hitting a flat spot.
    if (v <= var_floor) v <- 0.99 * var_floor + 0.01 * v
    D[i] <- v
    A[i, nb] <- a
  }
  if (any(D <= 0) || any(!is.finite(D))) return(NULL)
  M <- (diag(1, n) - A) / sqrt(D)             # row-scaled (I - A)
  Q_ord <- crossprod(M)
  Q <- matrix(0, n, n)
  Q[ord, ord] <- Q_ord                        # ordering -> original rows
  list(Q = (Q + t(Q)) / 2, log_det = -sum(log(D)))
}

# Log NNGP density of a field realisation, from the precision above. Exposed as
# its own helper so the equivalence test can score it against the compiled
# kernel without going through a fit.
.tobs_nngp_log_density <- function(z, Q, log_det) {
  0.5 * log_det - 0.5 * length(z) * log(2 * pi) -
    0.5 * as.numeric(crossprod(z, Q %*% z))
}

# Log hyperprior of one NNGP field at (sigma, phi), on the log-log grid measure.
#
# The same two densities the sampler carries (tulpa src/tulpa_priors_svc.h): a
# PC prior on the range (Fuglstad et al. 2019, d = 2, calibrated to
# P(range < U) = alpha from the term's `prior_range`) and a half-Cauchy on the
# marginal SD (scale from `sigma2_prior_scale`). The outer grid is log-spaced in
# both coordinates, so each density is taken to the log scale by its Jacobian --
# omitting that would weight the grid by an implicit uniform on the natural
# scale instead of the prior the term declares.
.tobs_svc_log_hyperprior <- function(sigma, phi, prior_range, sigma_scale) {
  lambda <- -log(prior_range[2L]) * prior_range[1L]
  lp_phi <- log(lambda) - 2 * log(phi) - lambda / phi + log(phi)
  lp_sig <- -log1p(sigma^2 / sigma_scale^2) + log(sigma)
  lp_phi + lp_sig
}

# Continuous NNGP latent-field block for the areal-BFGS driver.
#
# The continuous sibling of `.areal_field_car()`: one field value per site, the
# same closure interface (offset / scatter / prior_logp / prior_grad / center /
# constrain / to_phi / cells / to_hyper / valid). An unweighted surface loads
# eta += z; a varying-coefficient surface carries the design column as a
# per-site `weight` and loads eta += weight * z, with z reported unweighted --
# the coefficient surface itself.
#
# The prior is proper (a mean-zero GP), so there is no sum-to-zero constraint;
# the field level is anchored the way the sampler anchors it, by the soft
# mean penalty mean(z) ~ N(0, 1 / lambda_mean) that keeps the surface level from
# drifting against the intercept (tulpa src/hmc_svc.h svc_sum_to_zero_penalty).
#
# The precision is assembled densely, one matrix per outer-grid cell (as the
# areal blocks assemble theirs), so the block holds n^2 x n_cells doubles: this
# path targets the moderate-n regime the recovery tests measure, not an EVA-scale
# site count.
.tobs_svc_nngp_field <- function(co, nn_idx, nn_dist, ord, cov_code,
                                 prior_range, sigma_scale = 1,
                                 weight = NULL, sigma_grid = NULL,
                                 phi_grid = NULL, lambda_mean = 1) {
  n <- nrow(co)
  wtd <- !is.null(weight)
  if (is.null(sigma_grid)) sigma_grid <- exp(seq(log(0.3), log(3), length.out = 4L))
  if (is.null(phi_grid)) {
    U <- prior_range[1L]
    phi_grid <- exp(seq(log(U * 0.5), log(U * 10), length.out = 4L))
  }
  make_cell <- function(theta) {
    sigma <- theta[1L]; phi <- theta[2L]
    pr <- .tobs_nngp_precision(co, nn_idx, nn_dist, ord, sigma^2, phi, cov_code)
    list(sigma = sigma, phi = phi,
         Q = if (is.null(pr)) NULL else pr$Q,
         ldQ = if (is.null(pr)) -Inf else pr$log_det,
         lp_hyper = .tobs_svc_log_hyperprior(sigma, phi, prior_range, sigma_scale))
  }
  cells <- list()
  for (phi in phi_grid) for (sg in sigma_grid)
    cells[[length(cells) + 1L]] <- make_cell(c(sg, phi))

  list(
    n_field = n, n_sp = n, cells = cells, axes = NULL, make_cell = make_cell,
    valid = function(cell) !is.null(cell$Q) && is.finite(cell$ldQ),
    offset = if (wtd) function(fp, cell) weight * fp
             else function(fp, cell) fp,
    scatter = if (wtd) function(grad_eta) weight * grad_eta
              else function(grad_eta) grad_eta,
    prior_logp = function(fp, cell) {
      .tobs_nngp_log_density(fp, cell$Q, cell$ldQ) -
        0.5 * lambda_mean * n * mean(fp)^2 + cell$lp_hyper
    },
    prior_grad = function(fp, cell)                    # d(-log p) / d fp
      as.numeric(cell$Q %*% fp) + lambda_mean * mean(fp),
    center = function(fp) fp,
    constrain = rep(FALSE, n),
    to_phi = function(fp, cell) fp,
    type = "svc_nngp",
    to_hyper = function(cell) c(sigma = cell$sigma, phi = cell$phi)
  )
}


# Resolve an svc() term into a LIST of continuous NNGP field blocks for the
# areal-BFGS driver -- one coefficient surface per `indices` entry, loaded by that
# design column of the arm the term sits on (a column of ones is the unweighted
# intercept surface). Shared by every family that rides the driver, so the term's
# validation and hyperparameter grid live in one place rather than once per family
# (gcol33/tulpaObs#143, #144).
#
# `X_arm` is the design matrix of the arm the surfaces load on -- the occupancy
# logit for occu()/fp_occu(), log lambda for the count families -- so `indices` is
# read against the same design the coefficients are.
#
# Each surface is its own block and the driver's outer grid is the product of the
# blocks' grids, so the per-surface grid is coarsened when several coefficients
# vary.
.tobs_svc_field_blocks <- function(svc, X_arm, n_sites, family,
                                   sigma_grid = NULL, phi_grid = NULL) {
  arm <- .tobs_svc_arm_label(family)
  if (isTRUE(svc$shared[2L])) {
    stop(sprintf(paste0("svc() is wired on the %s predictor of %s(); a ",
                        "detection-arm varying coefficient is not yet fit on the ",
                        "Laplace backends. Move the term to the %s formula%s. ",
                        "(tulpaObs#143, #144)"),
                 arm, family, arm,
                 if (identical(family, "occu")) ", or use method = \"nuts\"" else ""),
         call. = FALSE)
  }
  if (as.integer(svc$n_obs) != n_sites) {
    stop(sprintf(paste0("svc() carries %d coordinates but the model has %d ",
                        "sites; one coordinate pair per site is required."),
                 as.integer(svc$n_obs), n_sites), call. = FALSE)
  }
  idx <- .tobs_svc_columns(svc, X_arm, family, arm)
  co      <- matrix(as.numeric(svc$coords), ncol = 2L, byrow = TRUE)
  nn_idx  <- matrix(as.integer(svc$nn_idx),  nrow = n_sites, byrow = TRUE)
  nn_dist <- matrix(as.numeric(svc$nn_dist), nrow = n_sites, byrow = TRUE)
  ord     <- as.integer(svc$nn_order)
  cov_code <- .tobs_nngp_cov_code(svc$cov_type)

  n_g <- if (length(idx) > 1L) 3L else 4L
  if (is.null(sigma_grid))
    sigma_grid <- exp(seq(log(0.3), log(3), length.out = n_g))
  if (is.null(phi_grid)) {
    U <- svc$prior_range[1L]
    phi_grid <- exp(seq(log(U * 0.5), log(U * 10), length.out = n_g))
  }
  lapply(idx, function(j) {
    w <- as.numeric(X_arm[, j])
    .tobs_svc_nngp_field(
      co, nn_idx, nn_dist, ord, cov_code,
      prior_range = as.numeric(svc$prior_range),
      sigma_scale = as.numeric(svc$sigma2_prior_scale %||% 1),
      weight = if (isTRUE(all(w == 1))) NULL else w,
      sigma_grid = sigma_grid, phi_grid = phi_grid)
  })
}

# Name of the arm an svc() surface loads on, per family, for error messages.
.tobs_svc_arm_label <- function(family) {
  switch(family,
         occu = , fp_occu = "occupancy",
         "abundance")
}

# An areal-BFGS fit exposes ONE per-observation eta gradient to its field blocks,
# on whichever arm its spatial term selected. A detection-arm areal field therefore
# hands the varying-coefficient blocks the detection gradient, which would fit the
# surfaces against the wrong arm; reject the combination rather than fit it.
.tobs_check_svc_arm <- function(svc, det_arm, family) {
  if (!is.null(svc) && isTRUE(det_arm))
    stop(sprintf(paste0("%s() fits svc() on the %s arm, so it cannot be combined ",
                        "with a detection-arm areal field in the same fit. Move ",
                        "the areal term to the %s formula, or drop svc(). ",
                        "(tulpaObs#144)"),
                 family, .tobs_svc_arm_label(family),
                 .tobs_svc_arm_label(family)), call. = FALSE)
  invisible(TRUE)
}


# Exact single-season occupancy marginal with analytic gradients, shaped as the
# areal-BFGS driver's family callback.
#
# The latent occupancy z integrates out in closed form (R/occu_marginal.R), and
# the score follows from the Fisher identity with w_i = P(z_i = 1 | y_i):
#
#   d logL / d eta_psi_i  = w_i - psi_i
#   d logL / d eta_p_ij   = w_i (y_ij - p_ij)     (valid visits)
#
# `offset` is the per-site field contribution to the occupancy logit. Returns
# `eval` (the driver's callback, also usable at offset = 0 for the field-free
# warm start), the packed parameter layout, and the weakly-informative
# fixed-effect prior the deterministic path applies.
.tobs_occu_svc_marginal <- function(model, prior_spec = NULL) {
  X_occ <- model$X_processes[[1L]]
  X_det <- model$X_processes[[2L]]
  X_dv  <- model$X_det_visit
  y <- model$y
  n_sites <- nrow(X_occ); max_visits <- ncol(y)
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)
  p_vis <- if (is.null(X_dv)) 0L else ncol(X_dv)

  valid <- y >= 0
  if (!any(valid)) stop("occu() + svc(): no valid visits in the data.", call. = FALSE)
  Y <- y; Y[!valid] <- 0
  any_det <- rowSums(Y * valid) > 0

  i_occ <- seq_len(p_occ)
  i_det <- p_occ + seq_len(p_det)
  i_vis <- if (p_vis > 0L) p_occ + p_det + seq_len(p_vis) else integer(0)
  n_fixed <- p_occ + p_det + p_vis

  pi_list <- model$process_info
  pr_o <- .prior_for_submodel(prior_spec, "psi", pi_list[[1L]]$coef_names)
  pr_d <- .prior_for_submodel(prior_spec, "p",
                              c(pi_list[[2L]]$coef_names, model$det_visit_names))
  pmean <- c(pr_o$mean, pr_d$mean)
  psd   <- c(pr_o$sd,   pr_d$sd)
  pprec <- ifelse(is.finite(psd) & psd > 0, 1 / psd^2, 0)
  if (length(pmean) != n_fixed) { pmean <- numeric(n_fixed); pprec <- numeric(n_fixed) }

  cl <- .tobs_clamp_eta
  eval <- function(theta_fix, offset = 0) {
    eta_psi <- cl(as.numeric(X_occ %*% theta_fix[i_occ]) + offset)
    psi <- stats::plogis(eta_psi)
    eta_site <- as.numeric(X_det %*% theta_fix[i_det])
    logit_p <- matrix(eta_site, n_sites, max_visits)
    if (p_vis > 0L) {
      logit_p <- logit_p + matrix(as.numeric(X_dv %*% theta_fix[i_vis]),
                                  n_sites, max_visits, byrow = TRUE)
    }
    p <- stats::plogis(cl(logit_p))
    logp   <- ifelse(valid, log(p),     0)
    log1mp <- ifelse(valid, log(1 - p), 0)
    S      <- rowSums(Y * logp + (1 - Y) * log1mp)
    prod0  <- exp(rowSums(log1mp))
    Lu     <- psi * prod0 + (1 - psi)
    ll <- sum(ifelse(any_det, log(psi) + S, log(Lu)))

    w <- ifelse(any_det, 1, psi * prod0 / Lu)
    g_psi <- w - psi
    resid <- (Y - p) * valid * w                       # per-site w recycles by row

    g <- numeric(n_fixed)
    g[i_occ] <- as.numeric(crossprod(X_occ, g_psi))
    g[i_det] <- as.numeric(crossprod(X_det, rowSums(resid)))
    if (p_vis > 0L) g[i_vis] <- as.numeric(crossprod(X_dv, as.numeric(t(resid))))

    d <- theta_fix - pmean
    list(log_lik  = ll - 0.5 * sum(pprec * d^2),
         grad_fixed = g - pprec * d,
         grad_eta = g_psi,
         weights  = w)
  }

  list(eval = eval, n_fixed = n_fixed, p_occ = p_occ, p_det = p_det,
       p_vis = p_vis, any_det = any_det, n_sites = n_sites)
}


# Fit single-season occupancy with continuous NNGP varying coefficients on the
# Laplace backends (gcol33/tulpaObs#143). `model` is the autoscaled tobs_model;
# `svc` the resolved `tobs_svc` term. Returns a `tobs_fit` whose `svc_field`
# matches the NUTS path's naming (a bare vector for one varying coefficient, an
# n_sites x n_svc matrix otherwise) and whose `svc_hyper` carries each surface's
# marginal SD and range, marginalised over the outer grid.
#
# `max_iter` is the per-cell BFGS iteration budget over (coefficients, surfaces),
# so it scales with the number of field values rather than with the EM iteration
# count `control$max.iter` carries elsewhere; the caller's value raises it but
# does not lower it below the floor.
.tobs_fit_occu_svc <- function(model, svc, priors = NULL, max_iter = 300L,
                               tol = 1e-8, verbose = TRUE, ...) {
  if (!identical(model$model_type, "single")) {
    stop("occu() + svc() on the Laplace backends fits a single-season model.",
         call. = FALSE)
  }
  n_sites <- model$n_sites
  X_occ <- model$X_processes[[1L]]
  idx <- .tobs_svc_columns(svc, X_occ, "occu")

  dots <- list(...)
  blocks <- .tobs_svc_field_blocks(svc, X_occ, n_sites, "occu",
                                   sigma_grid = dots$sigma.grid,
                                   phi_grid = dots$phi.grid)

  prior_spec <- .resolve_occu_priors(priors)
  marg <- .tobs_occu_svc_marginal(model, prior_spec)
  n_fixed <- marg$n_fixed

  # Field-free warm start on the same marginal (offset = 0), so the per-cell
  # BFGS starts from the non-spatial mode rather than the origin.
  th0 <- numeric(n_fixed)
  th0[1L] <- stats::qlogis(min(max(mean(marg$any_det), 0.05), 0.95))
  warm <- tryCatch(stats::optim(
    th0, function(t) -marg$eval(t, 0)$log_lik,
    function(t) -marg$eval(t, 0)$grad_fixed, method = "BFGS",
    control = list(maxit = 300L)), error = function(e) NULL)
  theta0_fix <- if (!is.null(warm) && all(is.finite(warm$par))) warm$par else th0

  bfgs_maxit <- max(as.integer(max_iter),
                    500L, 5L * (n_fixed + length(idx) * n_sites))
  res <- .tobs_areal_bfgs_fit(marg$eval, n_fixed, blocks, theta0_fix,
                              max_iter = bfgs_maxit, tol = tol,
                              label = "occu-svc-nngp", integration = "grid")

  pi_list <- model$process_info
  nm <- c(paste0(pi_list[[1L]]$name, "_", pi_list[[1L]]$coef_names),
          paste0(pi_list[[2L]]$name, "_", pi_list[[2L]]$coef_names))
  if (marg$p_vis > 0L) nm <- c(nm, paste0("p_visit_", model$det_visit_names))

  means <- res$beta_mean; names(means) <- nm
  V <- res$vcov; dimnames(V) <- list(nm, nm)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nm
  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- nm

  # Posterior-mean coefficient surfaces, in `indices` order, plus the per-site
  # occupancy weights P(z = 1 | y) at the integrated mode and its field.
  fields <- res$field_means
  surf <- matrix(unlist(fields), nrow = n_sites, ncol = length(fields))
  offset <- numeric(n_sites)
  for (k in seq_along(fields)) {
    w <- as.numeric(X_occ[, idx[k]])
    offset <- offset + w * as.numeric(fields[[k]])
  }
  ev <- marg$eval(means, offset)

  svc_field <- if (length(fields) == 1L) as.numeric(surf[, 1L]) else surf
  hyper <- res$hyper_means
  names(hyper) <- paste0("svc", seq_along(hyper))

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = V,
    skew = NULL, sla_status = "off",
    n_samples = n_draws, n_params = length(means),
    log_prob = rep(res$log_lik, n_draws)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names = nm, param_names = nm,
    process_info = pi_list,
    model = model,
    weights = ev$weights,
    svc_eta_offset = offset,
    log_lik = res$log_lik,
    N = sum(model$y >= 0),
    svc = svc,
    svc_field = svc_field,
    svc_hyper = hyper,
    svc_indices = idx,
    svc_coefficients = svc$coefficients,
    method = "nested_laplace",
    spatial_integration = res$integration,
    convergence = list(converged = TRUE, n_iter = NA_integer_),
    correction = "none"
  )), class = c("tobs_fit", "tulpa_fit"))
}
