# =============================================================================
# ms_occu_spatial.R - areal-spatial community single-season occupancy
# (ms_occu() + shared field; the occupancy analogue of sfMsNMix, tulpaObs#75).
#
# A per-species two-state occupancy model with Gaussian community hyperpriors on
# the per-species coefficients AND one shared ICAR / BYM2 / proper-CAR field on
# the OCCUPANCY arm:
#
#   logit psi_{s,i} = X_psi_i . (mu_psi + b_psi_s) + f_{u(i)}
#   logit p_{s,i}   = X_p_i   . (mu_p   + b_p_s)
#
# The latent z integrates out per species-site in closed form (the occupancy
# two-state marginal). Per outer grid point an in-tree C++ nested Laplace-EM
# (cpp_ms_occu_spatial_*) iterates the joint (mu, f, {b_s}) mode-find + a
# closed-form Sigma M-step, and R grid-integrates the community means / covariance
# / field over the field-hyperparameter posterior (law of total covariance). One
# spatial unit per site (or a many-to-one site -> unit map). Thin R wrappers over
# the C++ drivers, mirroring nmix_laplace_re_spatial.R for the count family.
# =============================================================================


# ---------------------------------------------------------------------------
# Shared grid post-processing
# ---------------------------------------------------------------------------

# Grid-integrate the per-grid EM results into the community-mean / covariance /
# field / Sigma / BLUP summaries build_ms_occu_fit() consumes, plus the spatial
# field posterior mean and the field-hyperparameter moments.
.ms_occu_spatial_post <- function(fit, p_psi, p_p, n_spatial, prior_type,
                                  theta_cols, phi_loadings = NULL) {
  d <- p_psi + p_p
  weights <- tulpa:::.nl_normalise_weights_safe(fit$log_marginal)
  ok <- is.finite(weights) & weights > 0
  if (!any(ok)) {
    stop("Spatial community occupancy: every grid point produced a non-finite ",
         "log-marginal. Check the adjacency graph / the hyperparameter grids.",
         call. = FALSE)
  }
  w <- weights; w[!ok] <- 0; w <- w / sum(w)

  modes  <- fit$modes
  mu_mat <- modes[, seq_len(d), drop = FALSE]
  mu     <- as.numeric(crossprod(w, mu_mat))

  # Community-mean covariance: law of total covariance over the grid.
  V <- .tobs_grid_vcov(mu_mat, w, fit$vcov_mu, center = mu)

  # Sigma + BLUPs: weighted means across the grid.
  S_n  <- nrow(as.matrix(fit$b_psi[[1L]]))
  Sps  <- matrix(0, p_psi, p_psi); Sp <- matrix(0, p_p, p_p)
  bpsi <- matrix(0, S_n, p_psi);   bp <- matrix(0, S_n, p_p)
  for (k in which(ok)) {
    Sps  <- Sps  + w[k] * as.matrix(fit$Sigma_psi[[k]])
    Sp   <- Sp   + w[k] * as.matrix(fit$Sigma_p[[k]])
    bpsi <- bpsi + w[k] * as.matrix(fit$b_psi[[k]])
    bp   <- bp   + w[k] * as.matrix(fit$b_p[[k]])
  }

  # Field posterior mean. ICAR / CAR: f (n_spatial). BYM2: phi = a v + b w.
  field_cols <- modes[, d + seq_len(ncol(modes) - d), drop = FALSE]
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
    hyper[[nm]] <- .tobs_weighted_moment(w, fit$theta_grid[, nm])
  }

  list(
    mu_psi = mu[seq_len(p_psi)], mu_p = mu[p_psi + seq_len(p_p)],
    vcov = V, Sigma_psi = Sps, Sigma_p = Sp,
    b_psi = bpsi, b_p = bp,
    log_lik = sum(w * fit$log_lik),
    converged = any(as.logical(fit$converged)),
    n_iter = max(as.integer(fit$n_iter)),
    spatial_field = field_mean, prior_type = prior_type,
    weights = weights, hyper = hyper)
}

# Per-species (n_valid, n_det) integer matrices [n_sites x n_species] for the C++
# spec, from the model summaries (single source: the same summaries the marginal
# uses).
.ms_occu_spatial_count_mats <- function(summaries, n_sites, n_species) {
  nv <- matrix(0L, n_sites, n_species)
  nd <- matrix(0L, n_sites, n_species)
  for (s in seq_len(n_species)) {
    nv[, s] <- as.integer(summaries[[s]]$n_valid[, 1L])
    nd[, s] <- as.integer(summaries[[s]]$n_det[, 1L])
  }
  list(n_valid = nv, n_det = nd)
}

# CSR adjacency from a spatial term (precomputed arrays if present, else from the
# dense graph). Mirrors .nmix_spatial_csr.
.ms_occu_spatial_csr <- function(spatial) {
  if (!is.null(spatial$adj_row_ptr)) {
    list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
         n_neighbors = spatial$n_neighbors)
  } else {
    adjacency_to_csr(spatial$graph)
  }
}

# Warm start shared across grid points: a crude community mean + small community
# covariances; each grid point runs its own EM to convergence.
.ms_occu_spatial_warm_start <- function(summaries, p_psi, p_p) {
  any_det_prop <- mean(vapply(summaries, function(z) mean(z$any_det), numeric(1)))
  clp <- function(q) min(max(q, 1e-3), 1 - 1e-3)
  mu <- numeric(p_psi + p_p)
  mu[1L] <- stats::qlogis(clp(any_det_prop))                       # psi intercept
  mu[p_psi + 1L] <- 0                                              # p intercept
  list(mu = mu, Sigma_psi = diag(0.25, p_psi), Sigma_p = diag(0.25, p_p))
}


# ---------------------------------------------------------------------------
# Field-kind wrappers
# ---------------------------------------------------------------------------

# ICAR
.ms_occu_community_icar <- function(model, csr, n_spatial, map_site_to_unit,
                                    tau_grid = NULL, max_iter = 100L,
                                    verbose = FALSE) {
  X_psi <- model$X_occ; X_p <- model$X_det
  p_psi <- ncol(X_psi); p_p <- ncol(X_p)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 6L))
  ws <- .ms_occu_spatial_warm_start(model$summaries, p_psi, p_p)
  mats <- .ms_occu_spatial_count_mats(model$summaries, model$n_sites,
                                      model$n_species)
  raw <- cpp_ms_occu_spatial_icar(
    X_psi = X_psi, X_p = X_p, n_valid = mats$n_valid, n_det = mats$n_det,
    map_site_to_unit = as.integer(map_site_to_unit), n_spatial = n_spatial,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, tau_grid = as.numeric(tau_grid),
    mu_init = ws$mu, Sigma_psi_init = ws$Sigma_psi, Sigma_p_init = ws$Sigma_p,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  .ms_occu_spatial_post(raw, p_psi, p_p, n_spatial, "icar", c("tau"))
}

# Proper CAR
.ms_occu_community_car_proper <- function(model, csr, n_spatial, map_site_to_unit,
                                          graph, tau_grid = NULL, rho_grid = NULL,
                                          max_iter = 100L, verbose = FALSE) {
  X_psi <- model$X_occ; X_p <- model$X_det
  p_psi <- ncol(X_psi); p_p <- ncol(X_p)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 5L))
  if (is.null(rho_grid)) rho_grid <- c(0.2, 0.5, 0.8)
  if (any(rho_grid <= 0) || any(rho_grid >= 1)) {
    stop("rho_grid values must lie strictly in (0, 1) for proper CAR.",
         call. = FALSE)
  }
  if (is.null(graph)) {
    stop("car_proper() community occupancy needs the adjacency graph to ",
         "precompute log|Q(rho)|.", call. = FALSE)
  }
  log_det_Q <- .tobs_car_logdet_Q(graph, rho_grid)
  ws <- .ms_occu_spatial_warm_start(model$summaries, p_psi, p_p)
  mats <- .ms_occu_spatial_count_mats(model$summaries, model$n_sites,
                                      model$n_species)
  raw <- cpp_ms_occu_spatial_car_proper(
    X_psi = X_psi, X_p = X_p, n_valid = mats$n_valid, n_det = mats$n_det,
    map_site_to_unit = as.integer(map_site_to_unit), n_spatial = n_spatial,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, tau_grid = as.numeric(tau_grid),
    rho_grid = as.numeric(rho_grid), log_det_Q_rho = as.numeric(log_det_Q),
    mu_init = ws$mu, Sigma_psi_init = ws$Sigma_psi, Sigma_p_init = ws$Sigma_p,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  .ms_occu_spatial_post(raw, p_psi, p_p, n_spatial, "car_proper", c("tau", "rho"))
}

# BYM2
.ms_occu_community_bym2 <- function(model, csr, n_spatial, map_site_to_unit,
                                    scale_factor, sigma_grid = NULL,
                                    rho_grid = NULL, max_iter = 100L,
                                    verbose = FALSE) {
  X_psi <- model$X_occ; X_p <- model$X_det
  p_psi <- ncol(X_psi); p_p <- ncol(X_p)
  if (is.null(sigma_grid)) sigma_grid <- exp(seq(log(0.2), log(3), length.out = 4L))
  if (is.null(rho_grid))   rho_grid   <- c(0.05, 0.4, 0.7, 0.95)
  if (any(sigma_grid <= 0)) stop("sigma_grid must be strictly positive.", call. = FALSE)
  if (any(rho_grid < 0) || any(rho_grid > 1)) {
    stop("rho_grid values must lie in [0, 1].", call. = FALSE)
  }
  ws <- .ms_occu_spatial_warm_start(model$summaries, p_psi, p_p)
  mats <- .ms_occu_spatial_count_mats(model$summaries, model$n_sites,
                                      model$n_species)
  raw <- cpp_ms_occu_spatial_bym2(
    X_psi = X_psi, X_p = X_p, n_valid = mats$n_valid, n_det = mats$n_det,
    map_site_to_unit = as.integer(map_site_to_unit), n_spatial = n_spatial,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, sigma_grid = as.numeric(sigma_grid),
    rho_grid = as.numeric(rho_grid), scale_factor = as.numeric(scale_factor),
    mu_init = ws$mu, Sigma_psi_init = ws$Sigma_psi, Sigma_p_init = ws$Sigma_p,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  phi_loadings <- function(theta_row) {
    sigma <- theta_row[["sigma"]]; rho <- theta_row[["rho"]]
    list(a = sigma * sqrt(rho / scale_factor), b = sigma * sqrt(1 - rho))
  }
  out <- .ms_occu_spatial_post(raw, p_psi, p_p, n_spatial, "bym2",
                               c("sigma", "rho"), phi_loadings = phi_loadings)
  out$scale_factor <- scale_factor
  out
}

# BYM2 Riebler scale factor (geometric mean of the non-zero ICAR-Q eigenvalues).
.ms_occu_bym2_scale_factor <- function(graph) {
  Q <- diag(rowSums(graph)) - graph
  eig <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  nz  <- eig[abs(eig) > 1e-10]
  exp(mean(log(nz)))
}


# ---------------------------------------------------------------------------
# Front-door fitter
# ---------------------------------------------------------------------------

# Fit the areal-spatial community occupancy model. `spatial` is the resolved
# tobs_spatial spec on the occupancy formula (icar / bym2 / car_proper). Returns
# a `tobs_fit` (via build_ms_occu_fit), with the spatial field + hyperparameters
# attached.
.tobs_fit_ms_occu_spatial <- function(model, spatial, max.iter = 100L,
                                      verbose = FALSE, ...) {
  .tobs_reject_weighted_spatial(spatial, "ms_occu() spatial")
  ptype <- spatial$type %||% "icar"
  if (!ptype %in% c("icar", "bym2", "car_proper")) {
    stop(sprintf(
      "ms_occu() areal field supports icar() / bym2() / car_proper(); got '%s'. ",
      ptype),
      "(car() is the improper non-intrinsic CAR; use icar() for the intrinsic ",
      "field. spde() / gp() are not wired for community occupancy.)",
      call. = FALSE)
  }
  csr   <- .ms_occu_spatial_csr(spatial)
  n_spatial <- spatial$n_units %||%
    (if (!is.null(spatial$graph)) nrow(spatial$graph) else length(csr$n_neighbors))
  if (n_spatial != model$n_sites) {
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for the community ",
                        "occupancy field."), n_spatial, model$n_sites), call. = FALSE)
  }

  # site -> spatial unit map: one unit per site (the areal occupancy field has
  # one node per occupancy unit).
  map_site_to_unit <- seq_len(model$n_sites)

  fit <- switch(ptype,
    icar = .ms_occu_community_icar(model, csr, n_spatial, map_site_to_unit,
                                   max_iter = max.iter, verbose = verbose),
    car_proper = .ms_occu_community_car_proper(
      model, csr, n_spatial, map_site_to_unit, graph = spatial$graph,
      max_iter = max.iter, verbose = verbose),
    bym2 = .ms_occu_community_bym2(
      model, csr, n_spatial, map_site_to_unit,
      scale_factor = spatial$scale_factor %||%
        .ms_occu_bym2_scale_factor(spatial$graph),
      max_iter = max.iter, verbose = verbose),
    stop(sprintf("ms_occu() areal field supports icar / bym2 / car_proper; got '%s'.",
                 ptype), call. = FALSE))

  # Reshape into the .tobs_community_em fit shape build_ms_occu_fit consumes.
  arm_idx <- list(psi = seq_len(ncol(model$X_occ)),
                  p   = ncol(model$X_occ) + seq_len(ncol(model$X_det)))
  b_list <- lapply(seq_len(model$n_species), function(s)
    c(fit$b_psi[s, ], fit$b_p[s, ]))
  fit_em <- list(
    mu = c(fit$mu_psi, fit$mu_p), global = numeric(0), b_list = b_list,
    Sigma = list(psi = fit$Sigma_psi, p = fit$Sigma_p),
    Vf = fit$vcov, logML = fit$log_lik,
    converged = isTRUE(fit$converged), n_iter = fit$n_iter)

  out <- build_ms_occu_fit(model, fit_em, arm_idx)
  out$method <- "nested_laplace"
  out$spatial <- list(type = ptype, field = fit$spatial_field,
                      hyper = fit$hyper, n_spatial = n_spatial,
                      scale_factor = fit$scale_factor)
  out$spatial_field <- fit$spatial_field
  out
}
