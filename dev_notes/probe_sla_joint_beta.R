# Quick smoke check that SLA on the joint cover-hurdle path also works with
# `positive = "beta"`. The beta arm theta_grid carries a `phi_pos` column;
# the lognormal arm uses the post-hoc sigma_pos as the noise SD.
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
sigma   <- 0.5; rho <- 0.7
w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)

x        <- rnorm(N)
beta_occ <- c(0.2, 0.7)
beta_pos <- c(-0.4, 0.3)
eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
occur   <- rbinom(N, 1L, plogis(eta_occ))
mu_pos  <- plogis(beta_pos[1] + beta_pos[2] * x + w_s[spatial_idx])
phi_beta <- 30
y <- numeric(N)
is_pos <- occur == 1L
if (sum(is_pos)) {
  y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi_beta,
                                  (1 - mu_pos[is_pos]) * phi_beta)
}
y <- pmin(y, 1 - 1e-6)

dat <- data.frame(x = x, region = factor(spatial_idx))
spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

ctrl <- list(
  sigma_grid     = c(0.3, 0.6),
  rho_grid       = c(0.5, 0.7),
  sigma_pos_grid = c(0.3, 0.6),
  phi_grid       = c(10, 30, 90)
)

cat("--- beta cover SLA fit ---\n")
old_msgs <- character(0)
withCallingHandlers(
  fit <- tobs(
    formula = ~ x, data = dat, family = cover("beta"),
    y = y, spatial = spatial, engine = "nested_laplace",
    approx = "simplified_laplace", control = ctrl
  ),
  message = function(m) {
    old_msgs <<- c(old_msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
cat(sprintf("  sla_status = %s\n", fit$sla_status))
cat(sprintf("  skew_occ   = length %d\n", length(fit$skew_occ %||% numeric(0))))
cat(sprintf("  skew_pos   = length %d\n", length(fit$skew_pos %||% numeric(0))))
cat(sprintf("  draws_occ  = %s\n",
            if (is.null(fit$draws_occ)) "NULL"
            else paste(dim(fit$draws_occ), collapse = " x ")))
cat(sprintf("  draws_pos  = %s\n",
            if (is.null(fit$draws_pos)) "NULL"
            else paste(dim(fit$draws_pos), collapse = " x ")))
fb <- grep("simplified_laplace is currently wired", old_msgs, value = TRUE)
cat(sprintf("  any wired-only fallback message? %s\n",
            if (length(fb)) "YES (BUG)" else "no"))
if (!is.null(fit$skew_occ)) { cat("  skew_occ:\n  "); print(fit$skew_occ) }
if (!is.null(fit$skew_pos)) { cat("  skew_pos:\n  "); print(fit$skew_pos) }

stopifnot(length(fit$skew_occ %||% numeric(0)) == 2L,
          length(fit$skew_pos %||% numeric(0)) == 2L,
          identical(dim(fit$draws_occ), c(1000L, 2L)),
          identical(dim(fit$draws_pos), c(1000L, 2L)))
cat("\nALL CHECKS PASS\n")
