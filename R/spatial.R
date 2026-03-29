#' Spatial specifications for occupancy models
#'
#' @description
#' Functions to specify spatial random effects for occupancy models.
#' These wrap tulpa's spatial infrastructure and handle NNGP neighbor
#' computation internally.
#'
#' @name tulpaOcc_spatial
NULL

#' ICAR spatial random effect
#'
#' @param adjacency Symmetric adjacency matrix (1 = neighbors, 0 = not)
#' @param shared Logical vector of length 2: which processes get the spatial
#'   effect. Default `c(TRUE, FALSE)` = occupancy only, not detection.
#' @return A `tulpaOcc_spatial` object
#' @export
occu_icar <- function(adjacency, shared = c(TRUE, FALSE)) {
  if (!is.matrix(adjacency)) stop("adjacency must be a matrix")
  if (!isSymmetric(adjacency)) stop("adjacency must be symmetric")
  n <- nrow(adjacency)

  # Convert to CSR (compressed sparse row)
  csr <- adjacency_to_csr(adjacency)

  structure(list(
    type = "icar",
    n_units = n,
    adj_row_ptr = csr$row_ptr,
    adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors,
    shared = shared
  ), class = "tulpaOcc_spatial")
}

#' BYM2 spatial random effect
#'
#' @param adjacency Symmetric adjacency matrix
#' @param scale_factor BYM2 scaling factor (computed from graph if NULL)
#' @param shared Which processes get spatial effect (default: occupancy only)
#' @return A `tulpaOcc_spatial` object
#' @export
occu_bym2 <- function(adjacency, scale_factor = NULL, shared = c(TRUE, FALSE)) {
  if (!is.matrix(adjacency)) stop("adjacency must be a matrix")
  if (!isSymmetric(adjacency)) stop("adjacency must be symmetric")
  n <- nrow(adjacency)

  csr <- adjacency_to_csr(adjacency)

  if (is.null(scale_factor)) {
    scale_factor <- compute_bym2_scale(adjacency)
  }

  structure(list(
    type = "bym2",
    n_units = n,
    adj_row_ptr = csr$row_ptr,
    adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors,
    scale_factor = scale_factor,
    shared = shared
  ), class = "tulpaOcc_spatial")
}

#' GP (NNGP) spatial random effect
#'
#' @param coords Matrix of coordinates (n_sites x 2)
#' @param cov Covariance function: "exponential", "matern", "gaussian"
#' @param nu Matern smoothness parameter (default 1.5)
#' @param nn Number of nearest neighbors for NNGP (default 15)
#' @param shared Which processes get spatial effect (default: occupancy only)
#' @param sigma2_prior_U Upper bound for half-Cauchy prior on sigma2
#' @param sigma2_prior_alpha Shape for half-Cauchy prior on sigma2
#' @param phi_prior_lower Lower bound for uniform prior on phi (range)
#' @param phi_prior_upper Upper bound for uniform prior on phi (range)
#' @return A `tulpaOcc_spatial` object
#' @export
occu_gp <- function(coords, cov = "exponential", nu = 1.5, nn = 15,
                   shared = c(TRUE, FALSE),
                   sigma2_prior_U = 1.0, sigma2_prior_alpha = 0.01,
                   phi_prior_lower = 0.01, phi_prior_upper = 10.0) {
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be a matrix with 2 columns")
  }
  cov <- match.arg(cov, c("exponential", "matern", "gaussian", "spherical"))
  n <- nrow(coords)
  nn <- min(nn, n - 1)

  nngp <- compute_nngp_neighbors(coords, nn)

  structure(list(
    type = "gp",
    n_obs = n,
    nn = nn,
    coords = as.vector(t(coords)),
    nn_idx = as.vector(t(nngp$nn_idx)),
    nn_dist = as.vector(t(nngp$nn_dist)),
    nn_neighbor_dist = as.vector(nngp$nn_neighbor_dist),
    nn_order = nngp$nn_order,
    nn_order_inv = nngp$nn_order_inv,
    cov_type = cov,
    nu = nu,
    shared = shared,
    sigma2_prior_U = sigma2_prior_U,
    sigma2_prior_alpha = sigma2_prior_alpha,
    phi_prior_lower = phi_prior_lower,
    phi_prior_upper = phi_prior_upper
  ), class = "tulpaOcc_spatial")
}

#' Multi-scale GP spatial random effect
#'
#' @param coords Matrix of coordinates (n_sites x 2)
#' @param cov Covariance function
#' @param nu Matern smoothness
#' @param nn_local Nearest neighbors for local scale (default 15)
#' @param nn_regional Nearest neighbors for regional scale (default 15)
#' @param shared Which processes get spatial effect
#' @param range_local_lower,range_local_upper Range bounds for local scale
#' @param range_regional_lower,range_regional_upper Range bounds for regional
#' @param sigma2_local_prior_U,sigma2_local_prior_alpha Prior for local sigma2
#' @param sigma2_regional_prior_U,sigma2_regional_prior_alpha Prior for regional sigma2
#' @return A `tulpaOcc_spatial` object
#' @export
occu_multiscale_gp <- function(coords, cov = "exponential", nu = 1.5,
                              nn_local = 15, nn_regional = 15,
                              shared = c(TRUE, FALSE),
                              range_local_lower = 0.01,
                              range_local_upper = 10.0,
                              range_regional_lower = 0.01,
                              range_regional_upper = 100.0,
                              sigma2_local_prior_U = 1.0,
                              sigma2_local_prior_alpha = 0.01,
                              sigma2_regional_prior_U = 1.0,
                              sigma2_regional_prior_alpha = 0.01) {
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be a matrix with 2 columns")
  }
  cov <- match.arg(cov, c("exponential", "matern", "gaussian", "spherical"))
  n <- nrow(coords)
  nn_local <- min(nn_local, n - 1)
  nn_regional <- min(nn_regional, n - 1)

  nngp_local <- compute_nngp_neighbors(coords, nn_local)
  nngp_regional <- compute_nngp_neighbors(coords, nn_regional)

  structure(list(
    type = "multiscale_gp",
    n_obs = n,
    coords = as.vector(t(coords)),
    nn_local = nn_local,
    nn_idx_local = as.vector(t(nngp_local$nn_idx)),
    nn_dist_local = as.vector(t(nngp_local$nn_dist)),
    nn_neighbor_dist_local = as.vector(nngp_local$nn_neighbor_dist),
    nn_order_local = nngp_local$nn_order,
    nn_order_inv_local = nngp_local$nn_order_inv,
    nn_regional = nn_regional,
    nn_idx_regional = as.vector(t(nngp_regional$nn_idx)),
    nn_dist_regional = as.vector(t(nngp_regional$nn_dist)),
    nn_neighbor_dist_regional = as.vector(nngp_regional$nn_neighbor_dist),
    nn_order_regional = nngp_regional$nn_order,
    nn_order_inv_regional = nngp_regional$nn_order_inv,
    cov_type = cov,
    nu = nu,
    shared = shared,
    range_local_lower = range_local_lower,
    range_local_upper = range_local_upper,
    range_regional_lower = range_regional_lower,
    range_regional_upper = range_regional_upper,
    sigma2_local_prior_U = sigma2_local_prior_U,
    sigma2_local_prior_alpha = sigma2_local_prior_alpha,
    sigma2_regional_prior_U = sigma2_regional_prior_U,
    sigma2_regional_prior_alpha = sigma2_regional_prior_alpha
  ), class = "tulpaOcc_spatial")
}

#' SPDE Spatial Random Effect (Matérn via Triangular Mesh)
#'
#' Specifies a continuous Matérn spatial field for occupancy models using the
#' SPDE approach. Builds a triangular mesh and FEM matrices via tulpaMesh,
#' then uses tulpa's CHOLMOD-accelerated SPDE Laplace engine.
#'
#' @param coords Matrix of coordinates (n_sites x 2) or a formula `~ x + y`.
#' @param data Optional data.frame for formula evaluation.
#' @param mesh A pre-built `tulpa_mesh` object (from tulpaMesh). If NULL,
#'   built automatically from coords.
#' @param max_edge Maximum mesh edge length. Scalar or `c(inner, outer)`.
#' @param cutoff Minimum distance between mesh vertices. Default 0.
#' @param nu Matérn smoothness parameter. Default 1. Fractional values (0.5, 1.5)
#'   use rational SPDE approximation.
#' @param shared Which processes get the spatial effect.
#'   Default `c(TRUE, FALSE)` = occupancy only.
#' @param prior_range Prior for range: `c(U, alpha)` where P(range < U) = alpha.
#' @param prior_sigma Prior for sigma: `c(U, alpha)` where P(sigma > U) = alpha.
#' @return A `tulpaOcc_spatial` object
#' @export
occu_spde <- function(coords, data = NULL, mesh = NULL,
                     max_edge = NULL, cutoff = 0,
                     nu = 1, shared = c(TRUE, FALSE),
                     prior_range = c(0.5, 0.5),
                     prior_sigma = c(1, 0.5)) {
  tulpa_spec <- tulpa::spatial_spde(
    coords = coords, data = data, mesh = mesh,
    max_edge = max_edge, cutoff = cutoff, nu = nu,
    prior_range = prior_range, prior_sigma = prior_sigma
  )

  structure(list(
    type = "spde",
    tulpa_spec = tulpa_spec,
    n_units = tulpa_spec$n_mesh,
    shared = shared,
    nu = nu,
    prior_range = prior_range,
    prior_sigma = prior_sigma
  ), class = "tulpaOcc_spatial")
}

#' @export
print.tulpaOcc_spatial <- function(x, ...) {
  n <- if (!is.null(x$n_units)) x$n_units else x$n_obs
  cat(sprintf("tulpaOcc spatial: %s (%d units)\n", x$type, n))
  shared_str <- ifelse(x$shared, "yes", "no")
  cat(sprintf("  Shared: psi=%s, p=%s\n", shared_str[1], shared_str[2]))
  if (x$type == "spde") {
    cat(sprintf("  Matern nu=%g, mesh=%d nodes\n", x$nu, x$n_units))
  }
  invisible(x)
}

# ============================================================================
# Internal helpers
# ============================================================================

adjacency_to_csr <- function(adj) {
  n <- nrow(adj)
  row_ptr <- integer(n + 1)
  col_idx <- integer(0)
  n_neighbors <- integer(n)

  for (i in seq_len(n)) {
    neighbors <- which(adj[i, ] > 0)
    n_neighbors[i] <- length(neighbors)
    col_idx <- c(col_idx, neighbors - 1L)  # 0-based
    row_ptr[i + 1] <- row_ptr[i] + n_neighbors[i]
  }

  list(
    row_ptr = as.integer(row_ptr),
    col_idx = as.integer(col_idx),
    n_neighbors = as.integer(n_neighbors)
  )
}

compute_bym2_scale <- function(adj) {
  n <- nrow(adj)
  D <- diag(rowSums(adj))
  Q <- D - adj
  # Geometric mean of generalized eigenvalues (excluding zero eigenvalue)
  evals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  evals <- evals[evals > 1e-10]
  exp(mean(log(evals)))
}

compute_nngp_neighbors <- function(coords, k) {
  N <- nrow(coords)
  order_idx <- order(coords[, 1], coords[, 2])
  coords_ordered <- coords[order_idx, , drop = FALSE]

  nn_idx <- matrix(0L, nrow = N, ncol = k)
  nn_dist <- matrix(Inf, nrow = N, ncol = k)
  nn_neighbor_dist <- array(0, dim = c(N, k, k))

  for (i in 2:N) {
    n_candidates <- min(i - 1, k)
    if (n_candidates > 0) {
      dists <- sqrt(
        (coords_ordered[1:(i-1), 1] - coords_ordered[i, 1])^2 +
        (coords_ordered[1:(i-1), 2] - coords_ordered[i, 2])^2
      )

      if (length(dists) <= k) {
        nn_order <- order(dists)
        nn_idx[i, seq_len(length(dists))] <- nn_order
        nn_dist[i, seq_len(length(dists))] <- dists[nn_order]
      } else {
        nn_order <- order(dists)[1:k]
        nn_idx[i, ] <- nn_order
        nn_dist[i, ] <- dists[nn_order]
      }

      n_neighbors <- sum(nn_idx[i, ] > 0)
      if (n_neighbors > 1) {
        neighbor_indices <- nn_idx[i, 1:n_neighbors]
        neighbor_coords <- coords_ordered[neighbor_indices, , drop = FALSE]
        for (j1 in 1:n_neighbors) {
          for (j2 in 1:n_neighbors) {
            if (j1 != j2) {
              nn_neighbor_dist[i, j1, j2] <- sqrt(
                (neighbor_coords[j1, 1] - neighbor_coords[j2, 1])^2 +
                (neighbor_coords[j1, 2] - neighbor_coords[j2, 2])^2
              )
            }
          }
        }
      }
    }
  }

  list(
    nn_idx = nn_idx,
    nn_dist = nn_dist,
    nn_neighbor_dist = nn_neighbor_dist,
    nn_order = order_idx,
    nn_order_inv = order(order_idx),
    k = k
  )
}
