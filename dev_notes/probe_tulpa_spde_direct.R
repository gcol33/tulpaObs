## Probe: call tulpa::tulpa_laplace directly with an SPDE spec on clean
## Bernoulli occupancy data (no EM). If the field magnitudes are still
## absurd here, the issue is in tulpa's SPDE Laplace engine, not in
## tulpaObs's M-step encoding.

suppressMessages({
  library(tulpa)
})

set.seed(42)
n_sites <- 400
coords <- cbind(runif(n_sites), runif(n_sites))

## True latent state z: smooth spatial pattern
u_true <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
x_cov <- rnorm(n_sites)
beta_true <- c(-0.5, 0.7)
eta <- beta_true[1] + beta_true[2] * x_cov + u_true
z <- rbinom(n_sites, 1, plogis(eta))

X <- cbind(1, x_cov)

## Build SPDE spec
sp <- spatial_spde(coords = coords, max_edge = c(0.3, 0.6), nu = 1,
                   prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5))
cat("n_mesh =", sp$n_mesh, "  n_obs =", n_sites, "\n\n")

## CHECK: FEM mass matrix from tulpaMesh
cat("C0_diag summary (FEM mass diagonal):\n"); print(summary(sp$C0_diag))
cat("Number of zero C0_diag:", sum(sp$C0_diag == 0), " of", sp$n_mesh, "\n")

mesh <- sp$mesh
cat("\nmesh class:", class(mesh), "\n")
cat("mesh n vertices:", nrow(mesh$loc), "  n triangles:", nrow(mesh$tv), "\n")

## Recompute FEM directly
fem <- tulpaMesh::fem_matrices(mesh, obs_coords = coords)
cat("\nfem$C dim:", dim(fem$C), "\n")
cat("fem$C diag summary:\n"); print(summary(Matrix::diag(fem$C)))
cat("fem$G dim:", dim(fem$G), "\n")
cat("fem$G diag summary:\n"); print(summary(Matrix::diag(fem$G)))
cat("fem$A dim:", dim(fem$A), "\n\n")

## Inspect Q construction at prior modes (range=0.3, sigma=0.7, nu=1, alpha=2):
##   Q = tau^2 * (k4 C + 2 k2 G + G C^-1 G)
## Marginal variance should be sigma^2 = 0.49 -> field magnitudes ~ 0.7.
range_use <- sp$prior_range[1]; sigma_use <- sp$prior_sigma[1]
kappa <- sqrt(8 * sp$nu) / range_use
tau   <- 1.0 / (sqrt(4 * pi) * kappa * sigma_use)
cat("kappa =", kappa, "  tau =", tau, "  tau^2 =", tau^2, "\n")

C0 <- sp$C0_diag
G  <- as(sp$G, "CsparseMatrix")
A  <- sp$A
n_mesh <- sp$n_mesh
## Build Q the same way the C++ does
C0_mat <- Matrix::Diagonal(n_mesh, C0)
C0_inv <- Matrix::Diagonal(n_mesh, 1 / pmax(C0, 1e-15))
Q_unitary <- (kappa^2)^2 * C0_mat + 2 * kappa^2 * G + G %*% C0_inv %*% G
Q <- tau^2 * Q_unitary
cat("Q diag range: [", min(Matrix::diag(Q)), ",", max(Matrix::diag(Q)), "]\n")

## Marginal variance of latent field = diag(Q^-1).
## Solve Q V = I one column at a time would be slow; just look at total signal:
cat("Mean Q_diag * sigma^2 (should be ~ 1):", mean(Matrix::diag(Q)) * sigma_use^2, "\n")

## Sample from N(0, Q^-1) by Cholesky and inspect spread
L <- Matrix::Cholesky(Q, LDL = FALSE)
set.seed(7)
w_draw <- as.numeric(Matrix::solve(L, rnorm(n_mesh)))
cat("Field draw summary (one prior sample):\n"); print(summary(w_draw))

## Fit directly
res <- tulpa_laplace(y = z, n_trials = rep(1L, n_sites), X = X,
                     family = "binomial", spatial = sp,
                     max_iter = 100L, tol = 1e-6)

cat("converged:", res$converged, "  n_iter:", res$n_iter, "\n")
cat("beta_hat:", res$mode[1:2], "  (truth:", beta_true, ")\n\n")

u_hat <- res$mode[3:length(res$mode)]
cat("field summary:\n"); print(summary(u_hat))
cat("\nfield-at-sites summary (A %*% u):\n")
field_at_sites <- as.numeric(sp$A %*% u_hat)
print(summary(field_at_sites))

cat("\nCompare to u_true at sites:\n"); print(summary(u_true))

cat("\nCorrelation field-at-sites vs u_true:",
    round(cor(field_at_sites, u_true), 3), "\n")

cat("\n\n===========================================\n")
cat("Now try tulpa::fit_spde (nested over range/sigma)\n")
cat("===========================================\n")
res2 <- fit_spde(y = z, X = X, spatial = sp, family = "binomial",
                 nested_laplace = TRUE, n_grid = 4L)
cat("converged:", res2$converged, "\n")
cat("beta_hat:", res2$mode[1:2], "  (truth:", beta_true, ")\n")
u_hat2 <- res2$mode[3:length(res2$mode)]
cat("field summary:\n"); print(summary(u_hat2))
field_at_sites2 <- as.numeric(sp$A %*% u_hat2)
cat("\nCorrelation field-at-sites vs u_true:",
    round(cor(field_at_sites2, u_true), 3), "\n")
cat("range used:", res2$range, " sigma used:", res2$sigma, "\n")
