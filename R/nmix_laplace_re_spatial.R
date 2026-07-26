# =============================================================================
# nmix_laplace_re_spatial.R — spatial community / multispecies N-mixture
#
# Thin R wrappers over the C++ nested Laplace-EM drivers
# (src/nmix_community_spatial.cpp) for the spatial community N-mixture
# (spAbundance sfMsNMix analogue): a per-species N-mixture with Gaussian
# community hyperpriors on the per-species coefficients AND one shared ICAR /
# BYM2 / proper-CAR field on the abundance arm. One spatial unit per site.
#
# Each wrapper builds the native NMixCommunityOracle (the pre-grouped, lgamma-
# cached data container, reused from the non-spatial community fit) and the CSR
# adjacency, calls the matching C++ grid driver, and grid-integrates the
# community means / covariance / field over the hyperparameter posterior:
#   - community-mean covariance V(mu) = sum_k w_k [vcov_mu_k + (mu_k - mubar)(mu_k - mubar)']
#     (law of total covariance over the grid -- the calibrated marginal SE);
#   - field posterior mean (f, or phi = a v + b w for BYM2) as a weighted mean;
#   - Sigma / per-species BLUPs as weighted means.
#
#   nmix_community_laplace_icar() / _bym2() / _car_proper()
# =============================================================================


# ---------------------------------------------------------------------------
# Shared post-processing: grid-integrate one cpp driver result.
# ---------------------------------------------------------------------------

# Build the warm-start community means / covariances common to every grid point.
# A cold restart per grid point (the C++ default) needs a single shared start; a
# crude one suffices since each grid point runs its own EM to convergence.
.nmix_community_warm_start <- function(y, p_lam, p_p) {
  mu <- c(log(mean(y) + 1), rep(0, p_lam - 1L), rep(0, p_p))
  list(mu = as.numeric(mu),
       Sigma_lambda = diag(0.25, p_lam),
       Sigma_p      = diag(0.25, p_p))
}

# Grid-integrate the per-grid EM results into the community-mean / covariance /
# field / Sigma / BLUP summaries `build_ms_nmix_fit()` consumes, plus the
# spatial-field posterior mean and hyperparameter moments.
.nmix_community_spatial_post <- function(fit, p_lam, p_p, n_spatial,
                                         prior_type, mixture, theta_cols,
                                         phi_loadings = NULL, weights = NULL) {
  d <- p_lam + p_p
  # Default: flat design weights over the fixed tensor grid. The mode-centred CCD
  # path (gcol33/tulpaObs#60) passes its own corrected R-INLA design weights.
  if (is.null(weights))
    weights <- tulpa:::.nl_normalise_weights_safe(fit$log_marginal)
  ok <- is.finite(weights) & weights > 0
  if (!any(ok)) {
    stop("Spatial community N-mixture: every grid point produced a non-finite ",
         "log-marginal. Check K_max / the adjacency graph / the grids.",
         call. = FALSE)
  }
  w <- weights; w[!ok] <- 0; w <- w / sum(w)

  modes <- fit$modes
  mu_mat <- modes[, seq_len(d), drop = FALSE]            # n_grid x d
  mu     <- as.numeric(crossprod(w, mu_mat))             # weighted-mean community means

  # Community-mean covariance: law of total covariance over the grid.
  V <- .tobs_grid_vcov(mu_mat, w, fit$vcov_mu, center = mu)

  # Sigma + BLUPs: weighted means across the grid.
  Sl <- matrix(0, p_lam, p_lam); Sp <- matrix(0, p_p, p_p)
  bl <- matrix(0, fit_n_species(fit), p_lam); bp <- matrix(0, nrow(bl), p_p)
  for (k in which(ok)) {
    Sl <- Sl + w[k] * as.matrix(fit$Sigma_lambda[[k]])
    Sp <- Sp + w[k] * as.matrix(fit$Sigma_p[[k]])
    bl <- bl + w[k] * as.matrix(fit$b_lambda[[k]])
    bp <- bp + w[k] * as.matrix(fit$b_p[[k]])
  }

  # Field posterior mean. ICAR / CAR: f (n_spatial). BYM2: phi = a v + b w.
  field_cols <- modes[, d + seq_len(ncol(modes) - d), drop = FALSE]
  if (identical(prior_type, "bym2")) {
    phi_grid <- matrix(0, nrow(modes), n_spatial)
    for (k in seq_len(nrow(modes))) {
      ab <- phi_loadings(fit$theta_grid[k, ])
      v <- field_cols[k, seq_len(n_spatial)]
      wv <- field_cols[k, n_spatial + seq_len(n_spatial)]
      phi_grid[k, ] <- ab$a * v + ab$b * wv
    }
    field_mean <- as.numeric(crossprod(w, phi_grid))
  } else {
    field_mean <- as.numeric(crossprod(w, field_cols))
  }

  # Hyperparameter posterior moments from the theta grid.
  hyper <- list()
  for (nm in theta_cols) {
    if (!(nm %in% colnames(fit$theta_grid))) next
    hyper[[nm]] <- .tobs_weighted_moment(w, fit$theta_grid[, nm])
  }
  r_summary <- if (identical(mixture, "negbin") && "r" %in% colnames(fit$theta_grid)) {
    list(r = unname(hyper$r["mean"]), log_r = NA_real_, r_sd = unname(hyper$r["sd"]))
  } else NULL

  list(
    mu_lambda    = mu[seq_len(p_lam)],
    mu_p         = mu[p_lam + seq_len(p_p)],
    vcov         = V,
    Sigma_lambda = Sl,
    Sigma_p      = Sp,
    b_lambda     = bl,
    b_p          = bp,
    log_lik      = sum(w * fit$log_lik),
    converged    = any(as.logical(fit$converged)),
    n_iter       = max(as.integer(fit$n_iter)),
    spatial_field = field_mean,
    prior_type   = prior_type,
    weights      = weights,
    hyper        = hyper,
    dispersion   = r_summary,
    boundary_max = max(fit$boundary_max, na.rm = TRUE),
    mixture      = mixture
  )
}

# Number of species = nrow of any per-grid BLUP matrix.
fit_n_species <- function(fit) nrow(as.matrix(fit$b_lambda[[1]]))

# Build the native community oracle (caches per species-site marginal). The
# spatial driver integrates the NB size r over the grid (Poisson kernel per grid
# point), so the oracle is always built Poisson (nb = FALSE); the per-site kernel
# is evaluated with the grid r directly by the driver.
.nmix_community_oracle <- function(lf, X_lambda, n_sites, n_species, K_max) {
  K_max <- if (is.null(K_max)) as.integer(max(lf$y) + 100L) else as.integer(K_max)
  list(ptr = cpp_nmix_community_oracle(lf$y, lf$site_idx, lf$species_idx,
                                       X_lambda, lf$X_p, n_sites, n_species,
                                       K_max, nb = FALSE),
       K_max = K_max)
}

# CSR adjacency from a spatial term (precomputed arrays if present, else from the
# dense graph). Mirrors .tobs_fit_nmix_spatial.
.nmix_spatial_csr <- function(spatial) {
  if (!is.null(spatial$adj_row_ptr)) {
    list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
         n_neighbors = spatial$n_neighbors)
  } else {
    adjacency_to_csr(spatial$graph)
  }
}


# ---------------------------------------------------------------------------
# ICAR
# ---------------------------------------------------------------------------

nmix_community_laplace_icar <- function(lf, X_lambda, n_sites, n_species,
                                        csr, n_spatial, mixture = "P",
                                        tau_grid = NULL, r_grid = NULL,
                                        K_max = NULL, max_iter = 100L,
                                        verbose = FALSE) {
  p_lam <- ncol(X_lambda); p_p <- ncol(lf$X_p)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 6L))
  r_grid_use <- .nmix_community_r_grid(mixture, r_grid)
  ws <- .nmix_community_warm_start(lf$y, p_lam, p_p)
  orc <- .nmix_community_oracle(lf, X_lambda, n_sites, n_species, K_max)
  raw <- .cpp_nmix_progress(cpp_nmix_community_spatial_icar,
    oracle = orc$ptr, map_site_to_unit_R = seq_len(n_sites),
    X_lambda_R = X_lambda,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, n_spatial = n_spatial,
    tau_grid = as.numeric(tau_grid), r_grid = as.numeric(r_grid_use),
    mu_init = ws$mu, Sigma_lambda_init = ws$Sigma_lambda,
    Sigma_p_init = ws$Sigma_p,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  out <- .nmix_community_spatial_post(raw, p_lam, p_p, n_spatial, "icar",
                                      .nmix_mix_label(mixture), c("tau", "r"))
  class(out) <- c("nmix_community_spatial_fit", "list")
  out
}


# ---------------------------------------------------------------------------
# Proper CAR
# ---------------------------------------------------------------------------

# log|Q(rho)| = log|D - rho W| per rho grid point (tau- and z-independent),
# precomputed from the dense graph via a Cholesky.
.nmix_car_logdet_Q <- function(graph, rho_grid) {
  D <- diag(rowSums(graph))
  vapply(rho_grid, function(rho) {
    Q <- D - rho * graph
    ch <- tryCatch(chol(Q), error = function(e) NULL)
    if (is.null(ch)) return(-Inf)            # Q(rho) not PD -> skip grid point
    2 * sum(log(diag(ch)))
  }, numeric(1))
}

nmix_community_laplace_car_proper <- function(lf, X_lambda, n_sites, n_species,
                                              csr, n_spatial, graph,
                                              mixture = "P",
                                              tau_grid = NULL, rho_grid = NULL,
                                              r_grid = NULL, K_max = NULL,
                                              max_iter = 100L, verbose = FALSE) {
  p_lam <- ncol(X_lambda); p_p <- ncol(lf$X_p)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 5L))
  if (is.null(rho_grid)) rho_grid <- c(0.2, 0.5, 0.8)
  if (any(rho_grid <= 0) || any(rho_grid >= 1)) {
    stop("rho_grid values must lie strictly in (0, 1) for proper CAR.",
         call. = FALSE)
  }
  if (is.null(graph)) {
    stop("car_proper() community N-mixture needs the adjacency graph to ",
         "precompute log|Q(rho)|.", call. = FALSE)
  }
  r_grid_use <- .nmix_community_r_grid(mixture, r_grid)
  log_det_Q  <- .nmix_car_logdet_Q(graph, rho_grid)
  ws <- .nmix_community_warm_start(lf$y, p_lam, p_p)
  orc <- .nmix_community_oracle(lf, X_lambda, n_sites, n_species, K_max)
  raw <- .cpp_nmix_progress(cpp_nmix_community_spatial_car_proper,
    oracle = orc$ptr, map_site_to_unit_R = seq_len(n_sites),
    X_lambda_R = X_lambda,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, n_spatial = n_spatial,
    tau_grid = as.numeric(tau_grid), rho_grid = as.numeric(rho_grid),
    log_det_Q_rho = as.numeric(log_det_Q), r_grid = as.numeric(r_grid_use),
    mu_init = ws$mu, Sigma_lambda_init = ws$Sigma_lambda,
    Sigma_p_init = ws$Sigma_p,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  out <- .nmix_community_spatial_post(raw, p_lam, p_p, n_spatial, "car_proper",
                                      .nmix_mix_label(mixture), c("tau", "rho", "r"))
  class(out) <- c("nmix_community_spatial_fit", "list")
  out
}


# ---------------------------------------------------------------------------
# BYM2
# ---------------------------------------------------------------------------

nmix_community_laplace_bym2 <- function(lf, X_lambda, n_sites, n_species,
                                        csr, n_spatial, scale_factor,
                                        mixture = "P",
                                        sigma_grid = NULL, rho_grid = NULL,
                                        r_grid = NULL, K_max = NULL,
                                        max_iter = 100L, verbose = FALSE) {
  p_lam <- ncol(X_lambda); p_p <- ncol(lf$X_p)
  if (is.null(sigma_grid)) sigma_grid <- exp(seq(log(0.2), log(3), length.out = 4L))
  if (is.null(rho_grid))   rho_grid   <- c(0.05, 0.4, 0.7, 0.95)
  if (any(sigma_grid <= 0)) stop("sigma_grid must be strictly positive.", call. = FALSE)
  if (any(rho_grid < 0) || any(rho_grid > 1)) {
    stop("rho_grid values must lie in [0, 1].", call. = FALSE)
  }
  r_grid_use <- .nmix_community_r_grid(mixture, r_grid)
  ws <- .nmix_community_warm_start(lf$y, p_lam, p_p)
  orc <- .nmix_community_oracle(lf, X_lambda, n_sites, n_species, K_max)
  raw <- .cpp_nmix_progress(cpp_nmix_community_spatial_bym2,
    oracle = orc$ptr, map_site_to_unit_R = seq_len(n_sites),
    X_lambda_R = X_lambda,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, n_spatial = n_spatial,
    sigma_grid = as.numeric(sigma_grid), rho_grid = as.numeric(rho_grid),
    scale_factor = as.numeric(scale_factor), r_grid = as.numeric(r_grid_use),
    mu_init = ws$mu, Sigma_lambda_init = ws$Sigma_lambda,
    Sigma_p_init = ws$Sigma_p,
    max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
  # phi = sigma (sqrt(rho/scale) v + sqrt(1-rho) w); the field-mean reconstruction
  # needs the (a, b) loadings per (sigma, rho) grid row.
  phi_loadings <- function(theta_row) {
    sigma <- theta_row[["sigma"]]; rho <- theta_row[["rho"]]
    list(a = sigma * sqrt(rho / scale_factor), b = sigma * sqrt(1 - rho))
  }
  out <- .nmix_community_spatial_post(raw, p_lam, p_p, n_spatial, "bym2",
                                      .nmix_mix_label(mixture),
                                      c("sigma", "rho", "r"),
                                      phi_loadings = phi_loadings)
  out$scale_factor <- scale_factor
  class(out) <- c("nmix_community_spatial_fit", "list")
  out
}


# ---------------------------------------------------------------------------
# SPDE (continuous Matern field shared across species)
# ---------------------------------------------------------------------------

# Shared continuous Matern field on the abundance arm:
#   log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s) + (A u)_i
# with u ~ N(0, Q(range, sigma)^{-1}) on the FEM mesh and A (n_sites x n_mesh)
# the projection shared across species. The outer grid integrates the SPDE
# hyperparameters (range, sigma) [, NB size r]; per grid point the proper Matern
# precision Q (and log|Q|) is built once on the R side (the same FEM assembly the
# single-species SPDE and the occupancy SPDE paths use) and passed into the same
# community Laplace-EM driver as the areal fields.
nmix_community_laplace_spde <- function(lf, X_lambda, n_sites, n_species,
                                        spatial, mixture = "P",
                                        range_grid = NULL, sigma_grid = NULL,
                                        r_grid = NULL, K_max = NULL,
                                        max_iter = 100L, integration = "grid",
                                        verbose = FALSE) {
  p_lam <- ncol(X_lambda); p_p <- ncol(lf$X_p)
  ts <- spatial$tulpa_spec
  if (is.null(ts) || !identical(ts$type, "spde")) {
    stop("nmix_community_laplace_spde() requires an SPDE tulpa_spec.", call. = FALSE)
  }
  n_mesh  <- ts$n_mesh
  A_dense <- as.matrix(ts$A)
  if (nrow(A_dense) != n_sites) {
    stop(sprintf("SPDE projection A has %d rows but the model has %d sites.",
                 nrow(A_dense), n_sites), call. = FALSE)
  }

  prior_range <- ts$prior_range
  prior_sigma <- ts$prior_sigma
  # Coarser default grid than the single-species SPDE path (3 x 3): each
  # community grid point is an n_species-fold-more-expensive EM, so the outer
  # (range x sigma [x r]) product is kept small -- the same rationale the areal
  # community fitters use for their coarser sigma/rho/r grids.
  if (is.null(range_grid)) {
    r_med <- prior_range[1]
    range_grid <- exp(seq(log(r_med * 0.4), log(r_med * 2.2), length.out = 3L))
  }
  if (is.null(sigma_grid)) {
    s_scale <- prior_sigma[1]
    sigma_grid <- exp(seq(log(s_scale * 0.4), log(s_scale * 1.8), length.out = 3L))
  }
  if (any(range_grid <= 0)) stop("range_grid must be strictly positive.", call. = FALSE)
  if (any(sigma_grid <= 0)) stop("sigma_grid must be strictly positive.", call. = FALSE)
  r_grid_use <- .nmix_community_r_grid(mixture, r_grid)
  is_nb <- identical(mixture, "NB")

  build_Q <- function(range_val, sigma_val) {
    kappa    <- sqrt(8 * ts$nu) / range_val
    tau_spde <- 1 / (sqrt(4 * pi) * kappa * sigma_val)
    Q <- Matrix::forceSymmetric(tulpa:::.spde_precision_Q(ts, kappa, tau_spde))
    list(Q = as.matrix(Q), log_det = .spde_logdet_Q(Q))
  }

  ws  <- .nmix_community_warm_start(lf$y, p_lam, p_p)
  orc <- .nmix_community_oracle(lf, X_lambda, n_sites, n_species, K_max)

  # Run the batched community EM driver over an explicit (range, sigma, r) grid,
  # folding the Matern PC prior on (range, sigma) into each log-marginal so the
  # integrated posterior is proper (the SPDE PC prior lives in R, as on the
  # single-species path). Returns the raw per-grid community fit.
  run_driver <- function(theta_grid) {
    n <- nrow(theta_grid)
    Q_list <- vector("list", n); log_dets <- numeric(n); pc_lp <- numeric(n)
    cache <- list()
    for (k in seq_len(n)) {
      key <- paste0(theta_grid[k, 1L], "_", theta_grid[k, 2L])
      if (is.null(cache[[key]])) cache[[key]] <- build_Q(theta_grid[k, 1L], theta_grid[k, 2L])
      Q_list[[k]] <- cache[[key]]$Q; log_dets[k] <- cache[[key]]$log_det
      pc_lp[k] <- tulpa:::pc_prior_log_density(theta_grid[k, 1L], theta_grid[k, 2L],
                                               prior_range, prior_sigma)
    }
    raw <- .cpp_nmix_progress(cpp_nmix_community_spatial_spde,
      oracle = orc$ptr, X_lambda_R = X_lambda, A_R = A_dense,
      Q_list = Q_list, log_det_Q = log_dets,
      theta_grid_R = theta_grid, r_grid = as.numeric(theta_grid[, 3L]),
      mu_init = ws$mu, Sigma_lambda_init = ws$Sigma_lambda,
      Sigma_p_init = ws$Sigma_p,
      max_iter_em = as.integer(max_iter), verbose = isTRUE(verbose))
    raw$log_marginal <- raw$log_marginal + pc_lp
    raw
  }

  post <- function(raw, weights = NULL)
    .nmix_community_spatial_post(raw, p_lam, p_p, n_mesh, "spde",
                                 .nmix_mix_label(mixture),
                                 c("range", "sigma", "r"), weights = weights)

  # ---- outer integration over (range, sigma [, r]): opt-in mode-centred CCD over
  # the log hyperparameters, declining to the fixed tensor grid (gcol33/tulpaObs#60).
  out <- NULL; integration_used <- "grid"; pareto_k <- NA_real_
  if (identical(integration, "ccd")) {
    rm0 <- prior_range[1]; sm0 <- prior_sigma[1]
    axes <- list(
      .tobs_ccd_axis("range", "log", lower = rm0 * 0.05, upper = rm0 * 20, start = rm0),
      .tobs_ccd_axis("sigma", "log", lower = sm0 * 0.05, upper = sm0 * 20, start = sm0))
    if (is_nb)
      axes <- c(axes, list(.tobs_ccd_axis("r", "log", lower = 0.5, upper = 40,
                                          start = exp(mean(log(c(0.5, 40)))))))
    eval_logm <- function(theta) {
      r_val <- if (is_nb) theta[3L] else Inf
      lm <- run_driver(matrix(c(theta[1L], theta[2L], r_val), 1L, 3L))$log_marginal[1L]
      if (is.finite(lm)) lm else NA_real_
    }
    cc <- tryCatch(.tobs_ccd_outer_grid(eval_logm, axes), error = function(e) NULL)
    if (!is.null(cc)) {
      r_col <- if (is_nb) cc$nodes[, "r"] else rep(Inf, nrow(cc$nodes))
      node_grid <- cbind(cc$nodes[, "range"], cc$nodes[, "sigma"], r_col)
      raw <- run_driver(node_grid)
      lm  <- raw$log_marginal; okc <- is.finite(lm)
      if (any(okc)) {
        w <- cc$dnode * exp(lm - max(lm[okc])); w[!okc] <- 0; w <- w / sum(w)
        out <- post(raw, weights = w); integration_used <- "ccd"; pareto_k <- cc$pareto_k
      }
    }
  }

  if (is.null(out)) {
    grid <- expand.grid(range = range_grid, sigma = sigma_grid,
                        r = r_grid_use, KEEP.OUT.ATTRS = FALSE)
    out <- post(run_driver(as.matrix(grid[, c("range", "sigma", "r")])))
  }

  out$spatial_integration <- integration_used
  out$spatial_pareto_k <- pareto_k
  class(out) <- c("nmix_community_spatial_fit", "list")
  out
}


# tulpaObs mixture code ("P"/"NB") -> public label ("poisson"/"negbin").
.nmix_mix_label <- function(mixture) {
  if (identical(mixture, "NB")) "negbin" else "poisson"
}

# NB size grid for the community-spatial path. Coarser than the single-species
# default (.nmix_resolve_r_grid, 6 points): each community-spatial grid point is
# an n_species-fold-more-expensive EM, so the outer (tau/sigma x r) product is
# kept small. Poisson -> the single r = Inf node.
.nmix_community_r_grid <- function(mixture, r_grid) {
  if (identical(mixture, "P")) return(Inf)
  if (is.null(r_grid)) return(exp(seq(log(0.5), log(40), length.out = 4L)))
  .nmix_resolve_r_grid(mixture, r_grid)
}
