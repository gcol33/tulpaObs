# =============================================================================
# jsdm_spatial.R - areal-spatial joint species distribution model
# (jsdm() + shared field; tulpaObs#76).
#
# The JSDM observes presence/absence directly (no detection process), with shared
# fixed-effect coefficients, a scalar per-species random intercept, and one shared
# ICAR / BYM2 / proper-CAR field on the latent occupancy:
#
#   logit psi_{s,i} = X_i . beta + b_s + f_{u(i)}
#   b_s ~ N(0, sigma_re^2),   one f shared across species
#
# There is no latent state to integrate out (y is observed), so the per-(species,
# site) contribution is a plain Bernoulli (jsdm_spatial_kernel.h). Per outer grid
# point an in-tree C++ Laplace-EM (cpp_jsdm_spatial_*) iterates the joint
# (beta, f, {b_s}) mode-find (the scalar b_s Schur-folded) + a closed-form
# sigma_re^2 M-step, and R grid-integrates the fixed effects / their covariance /
# the field over the field-hyperparameter posterior (law of total covariance).
# One spatial unit per site. Thin R wrappers over the C++ drivers, mirroring
# ms_occu_spatial.R for the community-occupancy field.
# =============================================================================


# ---------------------------------------------------------------------------
# Shared grid post-processing
# ---------------------------------------------------------------------------

# Grid-integrate the per-grid EM results into the fixed-effect means / covariance
# / field / sigma_re / BLUP summaries, plus the field-hyperparameter moments.
.jsdm_spatial_post <- function(fit, p_occ, n_spatial, prior_type, theta_cols,
                               phi_loadings = NULL) {
  weights <- tulpa:::.nl_normalise_weights_safe(fit$log_marginal)
  ok <- is.finite(weights) & weights > 0
  if (!any(ok)) {
    stop("Spatial JSDM: every grid point produced a non-finite log-marginal. ",
         "Check the adjacency graph / the hyperparameter grids.", call. = FALSE)
  }
  w <- weights; w[!ok] <- 0; w <- w / sum(w)

  modes    <- fit$modes
  beta_mat <- modes[, seq_len(p_occ), drop = FALSE]
  beta     <- as.numeric(crossprod(w, beta_mat))

  # Fixed-effect covariance: law of total covariance over the grid.
  V <- matrix(0, p_occ, p_occ)
  for (k in which(ok)) {
    Ck <- fit$vcov_beta[[k]]
    if (is.null(Ck) || anyNA(Ck)) next
    dk <- beta_mat[k, ] - beta
    V <- V + w[k] * (as.matrix(Ck) + tcrossprod(dk))
  }
  V <- (V + t(V)) / 2

  # sigma_re^2 + per-species BLUPs: weighted means across the grid.
  S_n  <- length(as.numeric(fit$blup[[1L]]))
  sig2 <- sum(w * fit$sigma_re2)
  blup <- numeric(S_n)
  for (k in which(ok)) blup <- blup + w[k] * as.numeric(fit$blup[[k]])

  # Field posterior mean. ICAR / CAR: f (n_spatial). BYM2: phi = a v + b w.
  field_cols <- modes[, p_occ + seq_len(ncol(modes) - p_occ), drop = FALSE]
  if (identical(prior_type, "bym2")) {
    phi_grid <- matrix(0, nrow(modes), n_spatial)
    for (k in seq_len(nrow(modes))) {
      ab <- phi_loadings(fit$theta_grid[k, ])
      v  <- field_cols[k, seq_len(n_spatial)]
      wv <- field_cols[k, n_spatial + seq_len(n_spatial)]
      phi_grid[k, ] <- ab$a * v + ab$b * wv
    }
    field_mean <- as.numeric(crossprod(w, phi_grid))
  } else {
    field_mean <- as.numeric(crossprod(w, field_cols))
  }

  hyper <- list()
  for (nm in theta_cols) {
    if (!(nm %in% colnames(fit$theta_grid))) next
    vec <- fit$theta_grid[, nm]
    m   <- sum(w * vec); s <- sqrt(max(0, sum(w * vec^2) - m^2))
    hyper[[nm]] <- c(mean = m, sd = s)
  }

  list(
    beta = beta, vcov = V, sigma_re = sqrt(max(0, sig2)), blup = blup,
    log_lik = sum(w * fit$log_lik),
    converged = any(as.logical(fit$converged)),
    n_iter = max(as.integer(fit$n_iter)),
    spatial_field = field_mean, prior_type = prior_type,
    weights = weights, hyper = hyper)
}

# CSR adjacency from a spatial term (precomputed arrays if present, else from the
# dense graph). Mirrors .ms_occu_spatial_csr.
.jsdm_spatial_csr <- function(spatial) {
  if (!is.null(spatial$adj_row_ptr)) {
    list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
         n_neighbors = spatial$n_neighbors)
  } else {
    adjacency_to_csr(spatial$graph)
  }
}

# Warm start shared across grid points: intercept at the overall logit-prevalence
# of presences, slopes at 0, a small species-RE variance.
.jsdm_spatial_warm_start <- function(model, p_occ) {
  pr <- min(max(mean(model$y_mat), 1e-3), 1 - 1e-3)
  beta <- numeric(p_occ); beta[1L] <- stats::qlogis(pr)
  list(beta = beta, sigma_re2 = 0.25)
}


# ---------------------------------------------------------------------------
# Field-kind wrappers
# ---------------------------------------------------------------------------

# ICAR
.jsdm_field_icar <- function(model, csr, n_spatial, map_site_to_unit,
                             tau_grid = NULL, max_iter = 100L, verbose = FALSE) {
  X <- model$X_occ
  p_occ <- ncol(X)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 6L))
  ws <- .jsdm_spatial_warm_start(model, p_occ)
  raw <- cpp_jsdm_spatial_icar(
    X = X, y = matrix(as.integer(model$y_mat), model$n_sites, model$n_species),
    map_site_to_unit = as.integer(map_site_to_unit), n_spatial = n_spatial,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, tau_grid = as.numeric(tau_grid),
    beta_init = ws$beta, sigma_re2_init = ws$sigma_re2,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  .jsdm_spatial_post(raw, p_occ, n_spatial, "icar", c("tau"))
}

# Proper CAR
.jsdm_field_car_proper <- function(model, csr, n_spatial, map_site_to_unit,
                                   graph, tau_grid = NULL, rho_grid = NULL,
                                   max_iter = 100L, verbose = FALSE) {
  X <- model$X_occ
  p_occ <- ncol(X)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 5L))
  if (is.null(rho_grid)) rho_grid <- c(0.2, 0.5, 0.8)
  if (any(rho_grid <= 0) || any(rho_grid >= 1)) {
    stop("rho_grid values must lie strictly in (0, 1) for proper CAR.",
         call. = FALSE)
  }
  if (is.null(graph)) {
    stop("car_proper() JSDM needs the adjacency graph to precompute log|Q(rho)|.",
         call. = FALSE)
  }
  log_det_Q <- .jsdm_car_logdet_Q(graph, rho_grid)
  ws <- .jsdm_spatial_warm_start(model, p_occ)
  raw <- cpp_jsdm_spatial_car_proper(
    X = X, y = matrix(as.integer(model$y_mat), model$n_sites, model$n_species),
    map_site_to_unit = as.integer(map_site_to_unit), n_spatial = n_spatial,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, tau_grid = as.numeric(tau_grid),
    rho_grid = as.numeric(rho_grid), log_det_Q_rho = as.numeric(log_det_Q),
    beta_init = ws$beta, sigma_re2_init = ws$sigma_re2,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  .jsdm_spatial_post(raw, p_occ, n_spatial, "car_proper", c("tau", "rho"))
}

# log|Q(rho)| = log|D - rho W| per rho grid point.
.jsdm_car_logdet_Q <- function(graph, rho_grid) {
  D <- diag(rowSums(graph))
  vapply(rho_grid, function(rho) {
    Q <- D - rho * graph
    ch <- tryCatch(chol(Q), error = function(e) NULL)
    if (is.null(ch)) return(-Inf)
    2 * sum(log(diag(ch)))
  }, numeric(1))
}

# BYM2
.jsdm_field_bym2 <- function(model, csr, n_spatial, map_site_to_unit,
                             scale_factor, sigma_grid = NULL, rho_grid = NULL,
                             max_iter = 100L, verbose = FALSE) {
  X <- model$X_occ
  p_occ <- ncol(X)
  if (is.null(sigma_grid)) sigma_grid <- exp(seq(log(0.2), log(3), length.out = 4L))
  if (is.null(rho_grid))   rho_grid   <- c(0.05, 0.4, 0.7, 0.95)
  if (any(sigma_grid <= 0)) stop("sigma_grid must be strictly positive.", call. = FALSE)
  if (any(rho_grid < 0) || any(rho_grid > 1)) {
    stop("rho_grid values must lie in [0, 1].", call. = FALSE)
  }
  ws <- .jsdm_spatial_warm_start(model, p_occ)
  raw <- cpp_jsdm_spatial_bym2(
    X = X, y = matrix(as.integer(model$y_mat), model$n_sites, model$n_species),
    map_site_to_unit = as.integer(map_site_to_unit), n_spatial = n_spatial,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, sigma_grid = as.numeric(sigma_grid),
    rho_grid = as.numeric(rho_grid), scale_factor = as.numeric(scale_factor),
    beta_init = ws$beta, sigma_re2_init = ws$sigma_re2,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  phi_loadings <- function(theta_row) {
    sigma <- theta_row[["sigma"]]; rho <- theta_row[["rho"]]
    list(a = sigma * sqrt(rho / scale_factor), b = sigma * sqrt(1 - rho))
  }
  out <- .jsdm_spatial_post(raw, p_occ, n_spatial, "bym2", c("sigma", "rho"),
                            phi_loadings = phi_loadings)
  out$scale_factor <- scale_factor
  out
}

# BYM2 Riebler scale factor (geometric mean of the non-zero ICAR-Q eigenvalues).
.jsdm_bym2_scale_factor <- function(graph) {
  Q <- diag(rowSums(graph)) - graph
  eig <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  nz  <- eig[abs(eig) > 1e-10]
  exp(mean(log(nz)))
}


# ---------------------------------------------------------------------------
# Fit builder
# ---------------------------------------------------------------------------

# Build a `tobs_fit` from the grid-integrated JSDM spatial fit, matching the
# non-spatial JSDM laplace fit shape (FE means + covariance, pseudo-draws,
# spatial field, per-species random intercept).
.build_jsdm_spatial_fit <- function(model, post, ptype, n_spatial) {
  pi_list <- model$process_info
  nms <- paste0("psi_", pi_list[[1L]]$coef_names)

  means <- post$beta; names(means) <- nms
  V <- post$vcov; dimnames(V) <- list(nms, nms)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nms

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- nms

  blup <- post$blup; names(blup) <- model$species_names
  re_effects <- list(species = data.frame(
    level = model$species_names,
    `(Intercept)` = as.numeric(blup),
    check.names = FALSE, stringsAsFactors = FALSE))

  structure(c(list(
    draws = draws, means = means, sds = sds,
    skew = NULL, sla_status = "off",
    n_samples = n_draws, n_params = length(means),
    log_prob = rep(NA_real_, n_draws)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names = nms, param_names = nms,
    intercepts = compute_intercepts(model, means),
    model = model, spatial = list(type = ptype, n_spatial = n_spatial),
    spatial_field = post$spatial_field,
    spatial_field_det = NULL,
    process_info = pi_list,
    method = "nested_laplace",
    re_effects = re_effects,
    jsdm_re = list(sigma_re = post$sigma_re, blup = blup),
    jsdm_spatial = list(type = ptype, field = post$spatial_field,
                        hyper = post$hyper, n_spatial = n_spatial,
                        scale_factor = post$scale_factor),
    aghq = NULL,
    convergence = list(converged = isTRUE(post$converged), n_iter = post$n_iter),
    correction = "none"
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# Front-door fitter
# ---------------------------------------------------------------------------

# Fit the areal-spatial JSDM. `spatial` is the resolved tobs_spatial spec on the
# occupancy formula (icar / bym2 / car_proper). Returns a `tobs_fit`, with the
# shared areal field + hyperparameters + per-species random intercept attached.
.tobs_fit_jsdm_spatial <- function(model, spatial, max.iter = 100L,
                                   verbose = FALSE, ...) {
  .tobs_reject_weighted_spatial(spatial, "jsdm() spatial")
  ptype <- spatial$type %||% "icar"
  if (!ptype %in% c("icar", "bym2", "car_proper")) {
    stop(sprintf(
      "jsdm() areal field supports icar() / bym2() / car_proper(); got '%s'. ",
      ptype),
      "(car() is the improper non-intrinsic CAR; use icar() for the intrinsic ",
      "field. spde() / gp() are not wired for the JSDM areal path.)",
      call. = FALSE)
  }
  csr <- .jsdm_spatial_csr(spatial)
  n_spatial <- spatial$n_units %||%
    (if (!is.null(spatial$graph)) nrow(spatial$graph) else length(csr$n_neighbors))
  if (n_spatial != model$n_sites) {
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for the JSDM areal field.",
                 n_spatial, model$n_sites), call. = FALSE)
  }
  map_site_to_unit <- seq_len(model$n_sites)

  post <- switch(ptype,
    icar = .jsdm_field_icar(model, csr, n_spatial, map_site_to_unit,
                            max_iter = max.iter, verbose = verbose),
    car_proper = .jsdm_field_car_proper(
      model, csr, n_spatial, map_site_to_unit, graph = spatial$graph,
      max_iter = max.iter, verbose = verbose),
    bym2 = .jsdm_field_bym2(
      model, csr, n_spatial, map_site_to_unit,
      scale_factor = spatial$scale_factor %||%
        .jsdm_bym2_scale_factor(spatial$graph),
      max_iter = max.iter, verbose = verbose),
    stop(sprintf("jsdm() areal field supports icar / bym2 / car_proper; got '%s'.",
                 ptype), call. = FALSE))

  .build_jsdm_spatial_fit(model, post, ptype, n_spatial)
}
