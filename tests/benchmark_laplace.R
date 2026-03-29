# Speed comparison: tulpaOcc Laplace vs inlaocc vs tulpaOcc NUTS
library(tulpa)
library(tulpaOcc)
library(INLAocc)

cat("=== Speed Comparison: Laplace vs INLA vs NUTS ===\n\n")

# ============================================================================
# 1. Simple model (200 sites)
# ============================================================================
cat("--- 1. Single-season, intercept-only (200 sites, 4 visits) ---\n")
set.seed(42)
N <- 200; J <- 4
z <- rbinom(N, 1, 0.6)
y <- matrix(rbinom(N * J, 1, z * 0.4), N, J)
d <- data.frame(x = rnorm(N))

t_laplace <- system.time({
  mod <- tulpaOcc::occu(~ 1, ~ 1, d, y)
  fit_l <- tulpaOcc::occu_fit(mod, method = "laplace", verbose = FALSE)
})

t_nuts <- system.time({
  fit_n <- tulpaOcc::occu_fit(mod, method = "nuts", iter = 2000, warmup = 1000,
                               seed = 42, verbose = FALSE)
})

t_inla <- system.time({
  dat_i <- INLAocc::occu_format(y, occ.covs = d)
  fit_i <- INLAocc::occu(~ 1, ~ 1, data = dat_i, verbose = 0)
})

cat(sprintf("  Laplace: %.2fs  |  NUTS: %.2fs  |  inlaocc: %.2fs\n",
            t_laplace["elapsed"], t_nuts["elapsed"], t_inla["elapsed"]))
cat(sprintf("  Laplace psi=%.3f  |  NUTS psi=%.3f\n",
            fit_l$intercepts$psi, fit_n$intercepts$psi))

# ============================================================================
# 2. With covariates (300 sites)
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

t_laplace <- system.time({
  mod2 <- tulpaOcc::occu(~ x1 + x2 + x3, ~ d1 + d2, d2, y2)
  fit_l <- tulpaOcc::occu_fit(mod2, method = "laplace", verbose = FALSE)
})

t_nuts <- system.time({
  fit_n <- tulpaOcc::occu_fit(mod2, method = "nuts", iter = 2000, warmup = 1000,
                               seed = 42, verbose = FALSE)
})

t_inla <- system.time({
  dat_i <- INLAocc::occu_format(y2, occ.covs = d2[, 1:3],
                                 det.covs = list(d1 = d2$d1, d2 = d2$d2))
  fit_i <- INLAocc::occu(~ x1 + x2 + x3, ~ d1 + d2, data = dat_i, verbose = 0)
})

cat(sprintf("  Laplace: %.2fs  |  NUTS: %.2fs  |  inlaocc: %.2fs\n",
            t_laplace["elapsed"], t_nuts["elapsed"], t_inla["elapsed"]))

# ============================================================================
# 3. Large (1000 sites)
# ============================================================================
cat("\n--- 3. Large: 1000 sites, 6 visits ---\n")
set.seed(42)
N <- 1000; J <- 6
d3 <- data.frame(x1 = rnorm(N), x2 = rnorm(N))
psi <- plogis(0.5 + 0.5 * d3$x1)
z <- rbinom(N, 1, psi)
y3 <- matrix(rbinom(N * J, 1, z * 0.5), N, J)

t_laplace <- system.time({
  mod3 <- tulpaOcc::occu(~ x1 + x2, ~ 1, d3, y3)
  fit_l <- tulpaOcc::occu_fit(mod3, method = "laplace", verbose = FALSE)
})

t_nuts <- system.time({
  fit_n <- tulpaOcc::occu_fit(mod3, method = "nuts", iter = 2000, warmup = 1000,
                               seed = 42, verbose = FALSE)
})

t_inla <- system.time({
  dat_i <- INLAocc::occu_format(y3, occ.covs = d3)
  fit_i <- INLAocc::occu(~ x1 + x2, ~ 1, data = dat_i, verbose = 0)
})

cat(sprintf("  Laplace: %.2fs  |  NUTS: %.2fs  |  inlaocc: %.2fs\n",
            t_laplace["elapsed"], t_nuts["elapsed"], t_inla["elapsed"]))

cat("\n=== Done ===\n")
