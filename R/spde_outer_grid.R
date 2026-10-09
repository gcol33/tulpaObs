# spde_outer_grid.R - the (range, sigma) outer grid of a continuous Matern
# (SPDE) field, and the per-grid-row FEM precision the nested-Laplace drivers
# integrate over (abun() + spde(), ms_occu() + spde(), ms_abun() + spde()). Per
# grid row the proper Matern precision Q(range, sigma) and log|Q| are assembled
# once on the R side (the same FEM assembly the occupancy SPDE path uses) and
# the PC log-prior on (range, sigma) is evaluated, so the driver's log-marginal
# can be folded into a proper integrated posterior.
#
# `ts` is the term's tulpa SPDE spec (mesh, FEM matrices, nu, prior_range,
# prior_sigma). The default grid spans `range_span` / `sigma_span` times the
# PC-prior medians on a log scale with `n_grid` nodes per axis; the community
# callers take the coarse default (each of their grid points is an
# n_species-fold-more-expensive EM), the single-species N-mixture path passes
# its wider 5-node spans.
.tobs_spde_outer_grid <- function(ts, range_grid = NULL, sigma_grid = NULL,
                                  range_span = c(0.4, 2.2),
                                  sigma_span = c(0.4, 1.8), n_grid = 3L) {
  prior_range <- ts$prior_range
  prior_sigma <- ts$prior_sigma
  if (is.null(range_grid)) {
    r_med <- prior_range[1]
    range_grid <- exp(seq(log(r_med * range_span[1]), log(r_med * range_span[2]),
                          length.out = n_grid))
  }
  if (is.null(sigma_grid)) {
    s_scale <- prior_sigma[1]
    sigma_grid <- exp(seq(log(s_scale * sigma_span[1]), log(s_scale * sigma_span[2]),
                          length.out = n_grid))
  }
  if (any(range_grid <= 0)) stop("range_grid must be strictly positive.", call. = FALSE)
  if (any(sigma_grid <= 0)) stop("sigma_grid must be strictly positive.", call. = FALSE)

  build_Q <- function(range_val, sigma_val) {
    kappa    <- sqrt(8 * ts$nu) / range_val
    tau_spde <- 1 / (sqrt(4 * pi) * kappa * sigma_val)
    Q <- Matrix::forceSymmetric(tulpa::tulpa_spde_precision_Q(ts, kappa, tau_spde))
    list(Q = as.matrix(Q), log_det = .spde_logdet_Q(Q))
  }

  # Per row of `theta_grid` (columns 1-2 = range, sigma; further axes such as an
  # NB size are the caller's and leave Q untouched): the dense precision, its
  # log-determinant and the PC log-prior. Rows sharing (range, sigma) build Q
  # once.
  precision <- function(theta_grid) {
    n <- nrow(theta_grid)
    Q_list <- vector("list", n); log_dets <- numeric(n); pc_lp <- numeric(n)
    cache <- list()
    for (k in seq_len(n)) {
      key <- paste0(theta_grid[k, 1L], "_", theta_grid[k, 2L])
      if (is.null(cache[[key]])) cache[[key]] <- build_Q(theta_grid[k, 1L], theta_grid[k, 2L])
      Q_list[[k]] <- cache[[key]]$Q; log_dets[k] <- cache[[key]]$log_det
      pc_lp[k] <- tulpa::tulpa_spde_log_hyperprior(
        theta_grid[k, 1L], theta_grid[k, 2L],
        list(prior_range = prior_range, prior_sigma = prior_sigma))
    }
    list(Q_list = Q_list, log_det = log_dets, pc_lp = pc_lp)
  }

  list(range_grid = range_grid, sigma_grid = sigma_grid,
       prior_range = prior_range, prior_sigma = prior_sigma,
       precision = precision)
}
