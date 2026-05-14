repo <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
suppressPackageStartupMessages(devtools::load_all(repo, quiet = TRUE))

set.seed(13)
n_s <- 25; N <- 200
spatial_idx <- sample.int(n_s, N, replace = TRUE)
phi_s <- rnorm(n_s, 0, 0.6)
x <- rnorm(N)
eta_occ <- -0.3 + 0.7 * x + phi_s[spatial_idx]
y_occ <- rbinom(N, 1, plogis(eta_occ))
eta_pos <- 0.4 - 0.5 * x + 1.0 * phi_s[spatial_idx]
y <- ifelse(y_occ == 1, exp(rnorm(N, eta_pos, 0.4)), 0)
y <- pmin(y, 1 - 1e-6)

dat <- data.frame(x = x, region = factor(spatial_idx))
nbr <- lapply(seq_len(n_s),
              function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
adj <- matrix(0L, n_s, n_s)
for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L
spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

fit <- tobs(formula = ~ x, data = dat, family = cover("lognormal"),
            y = y, spatial = spatial, engine = "nested_laplace",
            control = list(sigma_grid = c(0.4, 0.8),
                           rho_grid = c(0.5, 0.9),
                           alpha_grid = c(0.0, 1.0)))

j <- fit$joint
cat("Q_csc_n =", j$Q_csc_n, "\n")
cat("n_grid =", length(j$Q_csc_p_per_grid), "\n")
cat("layout$beta_start:", j$arm_layout$beta_start, "\n")
cat("layout$re_start:", j$arm_layout$re_start, "\n")
cat("layout$phi_start:", j$arm_layout$phi_start, "\n")
cat("layout$theta_start:", j$arm_layout$theta_start, "\n")
cat("layout$n_x:", j$arm_layout$n_x, "\n")

k <- which.max(j$weights)
cat("\n--- Heaviest grid k =", k, "weight =", round(j$weights[k], 3), "---\n")
Qp <- as.integer(j$Q_csc_p_per_grid[[k]])
Qi <- as.integer(j$Q_csc_i_per_grid[[k]])
Qx <- as.numeric(j$Q_csc_x_per_grid[[k]])
n_x <- j$Q_csc_n
cat("length(Qp) =", length(Qp), "(should be", n_x + 1, ")\n")
cat("range(Qi) =", range(Qi), "(should be in [0, n_x-1])\n")
cat("nnz =", length(Qx), "\n")
cat("Qp[1:6] =", Qp[1:6], "\n")
cat("Qi[1:10] =", Qi[1:10], "\n")
cat("Qx[1:10] =", round(Qx[1:10], 3), "\n")

# Assemble lower-triangle, then symmetrize
Qk_lt <- Matrix::sparseMatrix(
  i = Qi + 1L,
  p = Qp,
  x = Qx,
  dims = c(n_x, n_x),
  symmetric = FALSE,
  index1 = TRUE
)
Qk <- Matrix::forceSymmetric(Qk_lt, uplo = "L")

# Sanity: Q diagonal should all be > 0
Qd <- as.matrix(Qk)
cat("range(diag(Qd)):", range(diag(Qd)), "\n")
cat("isSymmetric(Qd):", isSymmetric(Qd), "\n")
cat("eigen min:", min(eigen(Qd, symmetric = TRUE,
                            only.values = TRUE)$values), "\n")
cat("Qd[1:4, 1:4]:\n"); print(round(Qd[1:4, 1:4], 4))

# Solve for first 4 columns (the betas)
beta_idx <- c(1, 2, 3, 4)  # 1-based
E <- Matrix::sparseMatrix(i = beta_idx, j = seq_along(beta_idx),
                          x = 1, dims = c(n_x, length(beta_idx)))
V <- Matrix::solve(Qk, E)
cat("\nQ^{-1}[beta, beta] diag at heaviest grid:\n")
for (jj in seq_along(beta_idx)) {
  cat("  beta_idx", beta_idx[jj], ": V[i,jj] =",
      round(V[beta_idx[jj], jj], 5), "sqrt =",
      round(sqrt(max(V[beta_idx[jj], jj], 0)), 4), "\n")
}

# Check: full dense inverse for comparison (n_x = 54 small enough)
Qinv <- solve(Qd)
cat("Dense diag(Qinv)[1:4]:", round(diag(Qinv)[1:4], 5), "\n")
cat("Dense sqrt(diag(Qinv))[1:4]:", round(sqrt(pmax(diag(Qinv)[1:4], 0)), 4), "\n")
