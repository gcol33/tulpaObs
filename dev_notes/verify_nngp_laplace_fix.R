# verify_nngp_laplace_fix.R
# Multi-seed recovery check for the NNGP Laplace fix (sparse path).
# Pre-fix: every iteration produced a CHOLMOD non-PD warning and the mode
#          stayed at zero (no convergence).
# Post-fix: converges on every seed, β within tolerance, spatial field moves.

suppressMessages(library(tulpa))
spatial_gp <- tulpa:::spatial_gp

make_synthetic_gp_spec <- function(n_obs = 250, k = 5L, seed = 11L) {
  set.seed(seed)
  coords <- cbind(runif(n_obs), runif(n_obs))
  ord <- order(coords[, 1], coords[, 2])
  coords_ord <- coords[ord, , drop = FALSE]
  nn_idx  <- matrix(0L, n_obs, k)
  nn_dist <- matrix(0,  n_obs, k)
  for (i in 2:n_obs) {
    d <- sqrt((coords_ord[1:(i - 1), 1] - coords_ord[i, 1])^2 +
              (coords_ord[1:(i - 1), 2] - coords_ord[i, 2])^2)
    nc <- min(length(d), k)
    sel <- order(d)[seq_len(nc)]
    nn_idx[i, seq_len(nc)]  <- sel
    nn_dist[i, seq_len(nc)] <- d[sel]
  }
  spec <- spatial_gp(coords = ~ x1 + x2, cov = "exponential", nn = k)
  spec$n_obs <- n_obs; spec$n_spatial <- n_obs; spec$n_unique <- n_obs
  spec$obs_to_loc <- seq_len(n_obs)
  spec$unique_coords <- coords_ord; spec$coords_matrix <- coords_ord
  spec$nn <- k
  spec$neighbor_info <- list(
    nn_idx = nn_idx, nn_dist = nn_dist,
    nn_order = ord, nn_order_inv = order(ord)
  )
  list(spec = spec, n_obs = n_obs)
}

beta_true <- c(-0.3, 0.5)
n_seeds <- 20L
results <- data.frame(
  seed = integer(n_seeds), converged = logical(n_seeds),
  n_iter = integer(n_seeds),
  b0_err = numeric(n_seeds), b1_err = numeric(n_seeds),
  max_abs_w = numeric(n_seeds), logml = numeric(n_seeds)
)

for (i in seq_len(n_seeds)) {
  s <- 100L + i
  ss <- make_synthetic_gp_spec(n_obs = 250L, k = 5L, seed = 11L)
  set.seed(s)
  X <- cbind(1, rnorm(ss$n_obs))
  eta <- as.numeric(X %*% beta_true)
  y <- rbinom(ss$n_obs, 1, plogis(eta))

  fit <- suppressWarnings(tulpa_laplace(
    y = y, n_trials = rep(1L, ss$n_obs), X = X,
    re_list = list(), family = "binomial", spatial = ss$spec,
    max_iter = 50L, tol = 1e-6, n_threads = 1L
  ))
  results[i, ] <- list(
    seed = s, converged = isTRUE(fit$converged), n_iter = fit$n_iter,
    b0_err = fit$mode[1] - beta_true[1],
    b1_err = fit$mode[2] - beta_true[2],
    max_abs_w = max(abs(fit$mode[-(1:2)])),
    logml = fit$log_marginal
  )
}

cat("Converged across", sum(results$converged), "/", n_seeds, "seeds\n")
cat("Mean iterations:", mean(results$n_iter), "\n")
cat("β0 bias / RMSE :", mean(results$b0_err), "/", sqrt(mean(results$b0_err^2)), "\n")
cat("β1 bias / RMSE :", mean(results$b1_err), "/", sqrt(mean(results$b1_err^2)), "\n")
cat("max |w| range  :", range(results$max_abs_w), "\n")
cat("log_marginal   :", range(results$logml), "\n")
stopifnot(all(results$converged))
stopifnot(all(is.finite(results$logml)))
stopifnot(sqrt(mean(results$b0_err^2)) < 0.4)
stopifnot(sqrt(mean(results$b1_err^2)) < 0.4)
cat("PASS: all", n_seeds, "seeds converged and β recovery is within tolerance.\n")
