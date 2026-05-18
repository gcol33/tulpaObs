# Smoke test for the joint-Laplace SLA correction on the cover hurdle.
#
# Mirrors the d3-like simulation used in
# `tests/testthat/test-cover-hurdle-nested-joint-recovery.R`: BYM2 spatial
# block on a chain adjacency, alpha = 1.0, demeaned (phi_f, theta_f) so the
# field matches the model's sum-to-zero prior. Two fits per seed: one under
# `approx = "gaussian_laplace"`, one under `approx = "simplified_laplace"`,
# and we confirm SLA returns sensible shapes (length of skew, dim of draws).
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa", quiet = TRUE, export_all = FALSE)
  load_all(".",         quiet = TRUE, export_all = FALSE)
})

set.seed(20260518L)
N    <- 200L
n_s  <- 25L
adj  <- matrix(0L, n_s, n_s)
for (s in seq_len(n_s)) {
  for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
}
spatial_idx <- sample.int(n_s, N, replace = TRUE)
phi_f   <- rnorm(n_s); phi_f   <- phi_f   - mean(phi_f)
theta_f <- rnorm(n_s); theta_f <- theta_f - mean(theta_f)
sigma   <- 0.6; rho <- 0.7; alpha <- 1.0
w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)

x        <- rnorm(N)
beta_occ <- c(0.2, 0.7)
beta_pos <- c(-1.5, 0.3)
eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
occur   <- rbinom(N, 1L, plogis(eta_occ))
eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[spatial_idx]
log_y   <- rnorm(N, eta_pos, 0.4)
y       <- ifelse(occur == 1L, exp(log_y), 0)
y       <- pmin(y, 1 - 1e-6)

dat <- data.frame(x = x, region = factor(spatial_idx))
spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

ctrl <- list(
  sigma_grid     = c(0.3, 0.6, 0.9),
  rho_grid       = c(0.5, 0.7, 0.9),
  sigma_pos_grid = c(0.3, 0.6, 0.9)
)

cat("--- fit under gaussian_laplace ---\n")
fit_gauss <- tobs(
  formula = ~ x, data = dat, family = cover("lognormal"),
  y = y, spatial = spatial, engine = "nested_laplace",
  approx = "gaussian_laplace", control = ctrl
)
cat(sprintf("  sla_status = %s\n", fit_gauss$sla_status))
cat(sprintf("  draws_occ  = %s\n",
            if (is.null(fit_gauss$draws_occ)) "NULL"
            else paste(dim(fit_gauss$draws_occ), collapse = " x ")))

cat("\n--- fit under simplified_laplace ---\n")
old_msgs <- character(0)
withCallingHandlers(
  fit_sla <- tobs(
    formula = ~ x, data = dat, family = cover("lognormal"),
    y = y, spatial = spatial, engine = "nested_laplace",
    approx = "simplified_laplace", control = ctrl
  ),
  message = function(m) {
    old_msgs <<- c(old_msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
fallback_msgs <- grep("simplified_laplace is currently wired", old_msgs,
                      value = TRUE)
cat(sprintf("  any 'wired-only' fallback message? %s\n",
            if (length(fallback_msgs)) "YES (BUG)" else "no"))
cat(sprintf("  sla_status = %s\n", fit_sla$sla_status))
cat(sprintf("  skew_occ   = length %d, finite? %s\n",
            length(fit_sla$skew_occ %||% numeric(0)),
            isTRUE(all(is.finite(fit_sla$skew_occ)))))
cat(sprintf("  skew_pos   = length %d, finite? %s\n",
            length(fit_sla$skew_pos %||% numeric(0)),
            isTRUE(all(is.finite(fit_sla$skew_pos)))))
cat(sprintf("  draws_occ  = %s\n",
            if (is.null(fit_sla$draws_occ)) "NULL"
            else paste(dim(fit_sla$draws_occ), collapse = " x ")))
cat(sprintf("  draws_pos  = %s\n",
            if (is.null(fit_sla$draws_pos)) "NULL"
            else paste(dim(fit_sla$draws_pos), collapse = " x ")))
if (!is.null(fit_sla$skew_occ)) {
  cat("  skew_occ values:\n  ");  print(fit_sla$skew_occ)
}
if (!is.null(fit_sla$skew_pos)) {
  cat("  skew_pos values:\n  ");  print(fit_sla$skew_pos)
}

# Expected shapes:
#   ncol(X_occ) = ncol(X_pos) = 2 (intercept + x).
#   draws_*    = 1000 x 2.
stopifnot(length(fit_sla$skew_occ %||% numeric(0)) == 2L,
          length(fit_sla$skew_pos %||% numeric(0)) == 2L,
          identical(dim(fit_sla$draws_occ), c(1000L, 2L)),
          identical(dim(fit_sla$draws_pos), c(1000L, 2L)),
          startsWith(fit_sla$sla_status, "simplified_laplace") ||
            startsWith(fit_sla$sla_status, "fallback_gaussian"))

cat("\nALL CHECKS PASS\n")
