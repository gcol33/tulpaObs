# Speed comparison: tulpaObs (NUTS) vs inlaocc (EM+INLA)
library(tulpa)
library(tulpaObs)
library(INLAocc)

cat("=== tulpaObs vs inlaocc Speed Comparison ===\n\n")

# ============================================================================
# 1. Single-season occupancy (intercept-only)
# ============================================================================
cat("--- 1. Single-season, intercept-only (200 sites, 4 visits) ---\n")
set.seed(42)
N <- 200; J <- 4
z <- rbinom(N, 1, 0.6)
y <- matrix(rbinom(N * J, 1, z * 0.4), N, J)
d <- data.frame(x = rnorm(N))

t1 <- system.time({
  mod_t <- tulpaObs::occu(~ 1, ~ 1, d, y)
  fit_t <- tulpaObs::occu_fit(mod_t, iter = 2000, warmup = 1000, seed = 42, verbose = FALSE)
})

t2 <- system.time({
  dat_i <- INLAocc::occu_format(y, occ.covs = d)
  fit_i <- INLAocc::occu(~ 1, ~ 1, data = dat_i, verbose = 0)
})

cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: %.2fs  |  ratio: %.1fx\n",
            t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"]))
cat(sprintf("  tulpaObs psi=%.3f\n", fit_t$intercepts$psi))

# ============================================================================
# 2. Single-season with covariates
# ============================================================================
cat("\n--- 2. Single-season, 3 occ + 2 det covariates (300 sites) ---\n")
set.seed(42)
N <- 300; J <- 5
d2 <- data.frame(x1 = rnorm(N), x2 = rnorm(N), x3 = rnorm(N),
                 d1 = rnorm(N), d2 = rnorm(N))
psi <- plogis(-0.5 + 0.8 * d2$x1 - 0.3 * d2$x2 + 0.5 * d2$x3)
p <- plogis(0.2 - 0.4 * d2$d1 + 0.3 * d2$d2)
z <- rbinom(N, 1, psi)
y2 <- matrix(rbinom(N * J, 1, z * p), N, J)

t1 <- system.time({
  mod_t <- tulpaObs::occu(~ x1 + x2 + x3, ~ d1 + d2, d2, y2)
  fit_t <- tulpaObs::occu_fit(mod_t, iter = 2000, warmup = 1000, seed = 42, verbose = FALSE)
})

t2 <- system.time({
  dat_i <- INLAocc::occu_format(y2, occ.covs = d2[, 1:3],
                                 det.covs = list(d1 = d2$d1, d2 = d2$d2))
  fit_i <- INLAocc::occu(~ x1 + x2 + x3, ~ d1 + d2, data = dat_i, verbose = 0)
})

cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: %.2fs  |  ratio: %.1fx\n",
            t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"]))

# ============================================================================
# 3. Larger dataset
# ============================================================================
cat("\n--- 3. Large: 1000 sites, 6 visits ---\n")
set.seed(42)
N <- 1000; J <- 6
d3 <- data.frame(x1 = rnorm(N), x2 = rnorm(N))
psi <- plogis(0.5 + 0.5 * d3$x1)
z <- rbinom(N, 1, psi)
y3 <- matrix(rbinom(N * J, 1, z * 0.5), N, J)

t1 <- system.time({
  mod_t <- tulpaObs::occu(~ x1 + x2, ~ 1, d3, y3)
  fit_t <- tulpaObs::occu_fit(mod_t, iter = 2000, warmup = 1000, seed = 42, verbose = FALSE)
})

t2 <- system.time({
  dat_i <- INLAocc::occu_format(y3, occ.covs = d3)
  fit_i <- INLAocc::occu(~ x1 + x2, ~ 1, data = dat_i, verbose = 0)
})

cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: %.2fs  |  ratio: %.1fx\n",
            t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"]))

# ============================================================================
# 4. Community model (multi-species)
# ============================================================================
cat("\n--- 4. Community: 100 sites, 5 species ---\n")
set.seed(42)
N <- 100; J <- 3; S <- 5
d4 <- data.frame(x = rnorm(N))
y4 <- list()
for (s in 1:S) {
  z_s <- rbinom(N, 1, plogis(rnorm(1, 0, 0.5)))
  y4[[paste0("sp", s)]] <- matrix(rbinom(N * J, 1, z_s * 0.4), N, J)
}

t1 <- system.time({
  mod_t <- tulpaObs::occu(~ x, ~ 1, d4, y4, species = TRUE)
  fit_t <- tulpaObs::occu_fit(mod_t, iter = 2000, warmup = 1000, seed = 42, verbose = FALSE)
})

t2 <- system.time({
  dat_i <- INLAocc::occu_format_ms(y4, occ.covs = d4)
  fit_i <- INLAocc::occu(~ x, ~ 1, data = dat_i, multispecies = TRUE, verbose = 0)
})

cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: %.2fs  |  ratio: %.1fx\n",
            t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"]))

# ============================================================================
# 5. Spatial GP model
# ============================================================================
cat("\n--- 5. Spatial GP (NNGP): 200 sites ---\n")
set.seed(42)
N <- 200; J <- 4
coords <- cbind(runif(N, 0, 10), runif(N, 0, 10))
d5 <- data.frame(x = rnorm(N))
# Simple spatial effect via distance-decay
dists <- as.matrix(dist(coords))
Sigma <- exp(-dists / 3)
w <- as.vector(chol(Sigma + diag(1e-6, N)) %*% rnorm(N))
psi <- plogis(0.5 * d5$x + w)
z <- rbinom(N, 1, psi)
y5 <- matrix(rbinom(N * J, 1, z * 0.5), N, J)

t1 <- system.time({
  mod_t <- tulpaObs::occu(~ x, ~ 1, d5, y5)
  sp_t <- tulpaObs::occu_gp(coords, nn = 10)
  fit_t <- tulpaObs::occu_fit(mod_t, spatial = sp_t, iter = 2000, warmup = 1000,
                               seed = 42, verbose = FALSE)
})

t2 <- system.time({
  dat_i <- INLAocc::occu_format(y5, occ.covs = d5, coords = coords)
  fit_i <- INLAocc::occu(~ x, ~ 1, data = dat_i, spatial = coords, verbose = 0)
})

cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: %.2fs  |  ratio: %.1fx\n",
            t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"]))

# ============================================================================
# 6. Dynamic occupancy
# ============================================================================
cat("\n--- 6. Dynamic: 100 sites, 5 seasons ---\n")
set.seed(42)
N <- 100; J <- 3; T_s <- 5
d6 <- data.frame(x = rnorm(N))
z <- matrix(0L, N, T_s)
z[, 1] <- rbinom(N, 1, 0.6)
for (t in 2:T_s) {
  z[, t] <- z[, t-1] * (1 - rbinom(N, 1, 0.1)) + (1 - z[, t-1]) * rbinom(N, 1, 0.2)
}
y6 <- array(0L, dim = c(N, J, T_s))
for (i in 1:N) for (t in 1:T_s) if (z[i, t] == 1) y6[i, , t] <- rbinom(J, 1, 0.5)

t1 <- system.time({
  mod_t <- tulpaObs::occu(~ 1, ~ 1, d6, y6, col_formula = ~ 1, ext_formula = ~ 1)
  fit_t <- tulpaObs::occu_fit(mod_t, iter = 2000, warmup = 1000, seed = 42, verbose = FALSE)
})

# inlaocc temporal uses different interface
t2 <- tryCatch(system.time({
  # Reshape for inlaocc temporal format
  dat_i <- INLAocc::occu_format(y6[,,1], occ.covs = d6)
  fit_i <- INLAocc::occu(~ 1, ~ 1, data = dat_i, temporal = "iid", verbose = 0)
}), error = function(e) {
  cat(sprintf("  inlaocc error: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(t2)) {
  cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: %.2fs  |  ratio: %.1fx\n",
              t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"]))
} else {
  cat(sprintf("  tulpaObs: %.2fs  |  inlaocc: N/A (different temporal API)\n",
              t1["elapsed"]))
}

cat("\n=== Done ===\n")
