## After the tulpaMesh fix: direct SPDE Laplace via tulpa::tulpa_laplace.
## Expect: converged, field magnitudes ~ sigma_prior (0.7), correlation with
## truth > 0.5.

suppressMessages({
  devtools::load_all("../tulpaMesh", quiet = TRUE)
  devtools::load_all("../tulpa", quiet = TRUE)
})

set.seed(42)
n_sites <- 400
coords <- cbind(runif(n_sites), runif(n_sites))

u_true <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
x_cov <- rnorm(n_sites)
beta_true <- c(-0.5, 0.7)
eta <- beta_true[1] + beta_true[2] * x_cov + u_true
z <- rbinom(n_sites, 1, plogis(eta))

X <- cbind(1, x_cov)

sp <- spatial_spde(coords = coords, max_edge = c(0.3, 0.6), nu = 1,
                   prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5))
cat("n_mesh =", sp$n_mesh, " n_obs =", n_sites, "\n\n")

opts <- options(warn = 1)
res <- tulpa_laplace(y = z, n_trials = rep(1L, n_sites), X = X,
                     family = "binomial", spatial = sp,
                     max_iter = 100L, tol = 1e-6)
options(opts)

cat("converged:", res$converged, "  n_iter:", res$n_iter, "\n")
cat("beta_hat:", res$mode[1:2], "  (truth:", beta_true, ")\n\n")

u_hat <- res$mode[3:length(res$mode)]
field_at_sites <- as.numeric(sp$A %*% u_hat)
cat("field-at-sites summary:\n"); print(summary(field_at_sites))
cat("\nu_true summary:\n"); print(summary(u_true))
cat("\nCorrelation field vs truth:",
    round(cor(field_at_sites, u_true), 3), "\n")
