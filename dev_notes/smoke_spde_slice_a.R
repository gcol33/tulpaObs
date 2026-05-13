## Smoke test for Slice A: SPDE end-to-end via tobs() + tobs_spde() on the
## post-refactor API. Validates that:
##   1. existing non-spatial occupancy still works (regression baseline)
##   2. tobs(spatial = tobs_spde(...)) fits without error
##   3. recovered beta_occ is reasonable
##   4. recovered spatial field correlates with the truth

suppressMessages({
  devtools::load_all("../tulpaMesh", quiet = TRUE)
  devtools::load_all("../tulpa",     quiet = TRUE)
  devtools::load_all(quiet = TRUE)
})

set.seed(42)

n_sites <- 400
J <- 8

coords <- cbind(runif(n_sites), runif(n_sites))

## Smooth spatial signal on the unit square (amplitude ~0.8)
u_true <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])

x_cov   <- rnorm(n_sites)
beta_occ_true <- c(-0.5, 0.7)
beta_det_true <- c(-0.3, 0.4)
det_cov <- rnorm(n_sites)

eta_occ <- beta_occ_true[1] + beta_occ_true[2] * x_cov + u_true
psi <- plogis(eta_occ)
z <- rbinom(n_sites, 1, psi)

eta_det <- beta_det_true[1] + beta_det_true[2] * det_cov
p_det <- plogis(eta_det)

y <- matrix(0L, n_sites, J)
for (i in seq_len(n_sites)) {
  if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p_det[i])
}

dat <- data.frame(occ_cov = x_cov, det_cov = det_cov,
                  x = coords[, 1], y_coord = coords[, 2])

cat("Truth: beta_occ =", beta_occ_true,
    " beta_det =", beta_det_true, "\n")
cat("Naive any-detection rate =", mean(rowSums(y) > 0), "\n\n")

## ---------------------------------------------------------------------------
## (1) Non-spatial baseline
## ---------------------------------------------------------------------------
cat("--- Non-spatial tobs(occu) ---\n")
fit_ns <- tobs(
  formula   = ~ occ_cov,
  data      = dat,
  family    = occu(),
  detection = ~ det_cov,
  y         = y,
  engine    = "laplace",
  control   = list(verbose = FALSE)
)
cat("Non-spatial means:\n"); print(fit_ns$means)

## ---------------------------------------------------------------------------
## (2) SPDE Laplace
## ---------------------------------------------------------------------------
sp <- tobs_spde(coords = coords, max_edge = c(0.3, 0.6), nu = 1,
                prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5))
cat("\nSPDE spec: n_mesh =", sp$n_units, "\n")

cat("\n--- tobs(occu) + spatial = tobs_spde(...) ---\n")
fit_sp <- tobs(
  formula   = ~ occ_cov,
  data      = dat,
  family    = occu(),
  detection = ~ det_cov,
  y         = y,
  spatial   = sp,
  engine    = "laplace",
  control   = list(verbose = TRUE)
)
cat("SPDE means:\n"); print(fit_sp$means)

cat("\nDifference (sp - ns):\n"); print(fit_sp$means - fit_ns$means)

## Field-at-sites recovery check
if (!is.null(fit_sp$spatial_field)) {
  u_hat <- fit_sp$spatial_field
  field_at_sites <- as.numeric(sp$tulpa_spec$A %*% u_hat)
  cat("\nField-at-sites summary:\n"); print(summary(field_at_sites))
  cat("u_true summary:\n"); print(summary(u_true))
  cat(sprintf("\nCorrelation field vs truth: %.3f\n",
              cor(field_at_sites, u_true)))
} else {
  cat("\n(no spatial field stored on fit — skipping field-recovery check)\n")
}

cat("\nDone.\n")
