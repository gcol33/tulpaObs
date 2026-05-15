# repro_nngp_laplace_cholmod.R
# Minimal repro: tulpa_laplace(spatial = spatial_gp(...)) with n_obs >= 200
# triggers CHOLMOD "not positive definite" warning every iteration and the
# mode stays at zero.
#
# Reuses make_synthetic_gp_spec() from tulpa's own test file, just pushes
# n_obs above SPARSE_THRESHOLD (200) so we land on the sparse Newton path.

suppressMessages({
  library(tulpa)
})

# spatial_gp() is internal in tulpa (NAMESPACE only exports tulpa_laplace*);
# pull it in directly so the repro doesn't depend on devtools::load_all().
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
  spec$n_obs        <- n_obs
  spec$n_spatial    <- n_obs
  spec$n_unique     <- n_obs
  spec$obs_to_loc   <- seq_len(n_obs)
  spec$unique_coords <- coords_ord
  spec$coords_matrix <- coords_ord
  spec$nn           <- k
  spec$neighbor_info <- list(
    nn_idx       = nn_idx,
    nn_dist      = nn_dist,
    nn_order     = ord,
    nn_order_inv = order(ord)
  )
  list(spec = spec, n_obs = n_obs, coords = coords_ord, ord = ord)
}

# n_obs = 250 → n_x = p + n_obs > SPARSE_THRESHOLD (200) → sparse path.
ss <- make_synthetic_gp_spec(n_obs = 250L, k = 5L, seed = 11L)
set.seed(101)
X  <- cbind(1, rnorm(ss$n_obs))
beta_true <- c(-0.3, 0.5)
eta <- as.numeric(X %*% beta_true)
y  <- rbinom(ss$n_obs, 1, plogis(eta))

cat("Fitting NNGP Laplace with n_obs=", ss$n_obs,
    " (sparse path: n_x=", ncol(X) + ss$spec$n_spatial, ")\n", sep = "")

fit <- tulpa_laplace(
  y = y, n_trials = rep(1L, ss$n_obs), X = X,
  re_list = list(),
  family = "binomial",
  spatial = ss$spec,
  max_iter = 50L, tol = 1e-6, n_threads = 1L,
  return_hessian = TRUE
)

cat("converged       :", isTRUE(fit$converged), "\n")
cat("n_iter          :", fit$n_iter, "\n")
cat("log_marginal    :", fit$log_marginal, "\n")
cat("beta_true       :", beta_true, "\n")
cat("beta_hat        :", fit$mode[1:2], "\n")
cat("max |w|         :", max(abs(fit$mode[-(1:2)])), "\n")
cat("max |w| > 1e-6? :", max(abs(fit$mode[-(1:2)])) > 1e-6, "\n")
