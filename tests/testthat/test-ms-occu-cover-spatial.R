# test-ms-occu-cover-spatial.R - Stage-1 reduced-rank spatial-factor community
# occu_cover (gcol33/tulpa#67). Currently exercises the ground-truth simulator;
# the fitter recovery tests are added with the fitter increments.

.mscs_grid_adj <- function(nr, nc) {
  N <- nr * nc
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (c - 1L) * nr + r
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    k <- idx(r, c)
    if (r > 1L)  adj[k, idx(r - 1L, c)] <- 1L
    if (r < nr)  adj[k, idx(r + 1L, c)] <- 1L
    if (c > 1L)  adj[k, idx(r, c - 1L)] <- 1L
    if (c < nc)  adj[k, idx(r, c + 1L)] <- 1L
  }
  adj
}

test_that("simulate_ms_occu_cover_spatial returns well-formed K=1 community data", {
  adj <- .mscs_grid_adj(6L, 6L)        # N = 36 cells
  N <- nrow(adj); J <- 4L; S <- 8L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = J, seed = 11L)

  expect_identical(dim(sim$y),     c(N, J, S))
  expect_identical(dim(sim$y_pos), c(N, J, S))
  expect_length(sim$species, S)

  # Detections are 0/1; cover is present exactly at detected visits and positive
  # (lognormal support).
  expect_true(all(sim$y %in% c(0L, 1L)))
  expect_true(all(is.na(sim$y_pos) == (sim$y != 1L)))
  expect_true(all(sim$y_pos[!is.na(sim$y_pos)] > 0))

  # A detection implies presence: y == 1 only where the latent z == 1.
  for (s in seq_len(S)) {
    det_any <- rowSums(sim$y[, , s] == 1L) > 0
    expect_true(all(sim$truth$z[det_any, s] == 1L))
  }
})

test_that("K > 1 draws lower-triangular loadings and K unit-scale fields", {
  adj <- .mscs_grid_adj(7L, 7L)
  N <- nrow(adj); S <- 10L; K <- 3L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, sd_load = 1.0,
                                        seed = 42L)
  expect_identical(sim$truth$K, K)
  expect_identical(dim(sim$truth$L), c(S, K))   # S x K loading matrix
  expect_identical(dim(sim$truth$w), c(N, K))   # N x K field matrix

  # Lower-triangular canonical form: factor k loads on species k..S only, with a
  # positive diagonal loading.
  L <- sim$truth$L
  for (k in seq_len(K)) {
    if (k > 1L) expect_true(all(L[seq_len(k - 1L), k] == 0))
    expect_gt(L[k, k], 0)
  }

  # Each field is a centred, unit-ish-scale ICAR draw.
  for (k in seq_len(K)) {
    expect_lt(abs(mean(sim$truth$w[, k])), 1e-8)
    expect_gt(stats::sd(sim$truth$w[, k]), 0.3)
    expect_lt(stats::sd(sim$truth$w[, k]), 3)
  }

  # Rank-K structure: the per-cell psi map has more than one degree of spatial
  # freedom across species (a rank-1 / single-field model cannot produce this).
  sv <- svd(scale(sim$truth$psi, center = TRUE, scale = FALSE))$d
  expect_gt(sv[2L] / sv[1L], 0.1)               # a genuine second spatial axis
})

test_that("K = 1 simulator output is unchanged (Stage-1 shapes preserved)", {
  adj <- .mscs_grid_adj(6L, 6L)
  a <- simulate_ms_occu_cover_spatial(adj, n_species = 8L, J = 4L, seed = 11L)
  b <- simulate_ms_occu_cover_spatial(adj, n_species = 8L, J = 4L, K = 1L,
                                      seed = 11L)
  expect_null(dim(a$truth$L))                   # vector, not matrix
  expect_null(dim(a$truth$w))
  expect_equal(a$truth$L, b$truth$L)
  expect_equal(a$truth$w, b$truth$w)
  expect_identical(a$y, b$y)
})

test_that("the shared factor is unit-scaled and sign-anchored", {
  adj <- .mscs_grid_adj(7L, 7L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 6L, seed = 7L)
  w <- sim$truth$w

  expect_length(w, nrow(adj))
  expect_lt(abs(mean(w)), 1e-8)              # centred on the constrained space
  expect_gt(stats::sd(w), 0.3)               # genuine spatial amplitude
  # Sorbye-Rue unit-marginal scale: the field SD is O(1), not collapsed/blown up.
  expect_lt(stats::sd(w), 3)
  # Canonical sign anchor: reference species loading is positive.
  expect_gt(sim$truth$L[1L], 0)
})

test_that("loadings give per-species range heterogeneity (not the naive shared map)", {
  adj <- .mscs_grid_adj(8L, 8L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 12L, sd_load = 1.2,
                                        seed = 99L)
  w <- sim$truth$w; L <- sim$truth$L

  # The spatial contribution to each species' occupancy predictor is L_s * w.
  # With K = 1 these are collinear in shape but differ in sign and amplitude:
  # species with opposite-sign loadings have ANTI-correlated spatial maps, which
  # a single shared field + intercept RE (naive structure) cannot produce.
  expect_gt(diff(range(L)), 0.5)             # loadings span a real range
  expect_true(any(L > 0) && any(L < 0))      # both signs present
  pos <- which(L > 0)[1L]; neg <- which(L < 0)[1L]
  expect_lt(stats::cor(L[pos] * w, L[neg] * w), -0.99)  # opposite maps

  # Per-cell psi genuinely varies across species (maps are not identical).
  cell_sd <- apply(sim$truth$psi, 1L, stats::sd)
  expect_gt(mean(cell_sd), 0.02)
})

test_that("penalised joint gradient matches finite differences", {
  adj <- .mscs_grid_adj(4L, 4L)            # N = 16 cells (small for FD)
  S <- 3L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 3L,
                                        n_occ_covs = 1L, n_det_covs = 1L,
                                        n_pos_covs = 1L, seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj)

  d <- tulpaObs:::.ms_ocs_dims(model)
  # Pack a parameter vector near (but not at) the truth so the gradient is
  # non-trivial; FD validity is point-independent.
  tr <- sim$truth
  mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))   # species-major
  par <- c(mu, b, tr$L, tr$w, log(tr$sigma_pos))
  set.seed(1L)
  par <- par + stats::rnorm(length(par), 0, 0.05)

  # Fixed hyperparameters for the inner objective.
  Sinv <- diag(c(rep(1 / 0.4^2, d$P_occ), rep(1 / 0.4^2, d$P_p),
                 rep(1 / 0.3^2, d$P_pos)))
  Pmu  <- diag(1 / 25, d$P)
  inv_sdL2 <- 1 / 1.0^2
  tau_w <- 1.3

  out <- tulpaObs:::.ms_ocs_penll_grad(model, par, Sinv, Pmu, inv_sdL2, tau_w,
                                       grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_penll_grad(model, p, Sinv, Pmu, inv_sdL2,
                                                 tau_w, grad = FALSE)$ll
  h <- 1e-5
  gnum <- numeric(length(par))
  for (k in seq_along(par)) {
    pp <- par; pp[k] <- pp[k] + h
    pm <- par; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the full packed gradient")
})

test_that("penalised joint gradient matches finite differences at K = 2", {
  adj <- .mscs_grid_adj(4L, 4L)            # N = 16 cells
  S <- 4L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                        n_occ_covs = 1L, n_det_covs = 1L,
                                        n_pos_covs = 1L, seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K)

  d <- tulpaObs:::.ms_ocs_dims(model)
  expect_identical(d$K, K)
  tr <- sim$truth
  mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))   # species-major
  # tr$L is S x K, tr$w is N x K; column-major vec matches the packed layout.
  par <- c(mu, b, as.numeric(tr$L), as.numeric(tr$w), log(tr$sigma_pos))
  set.seed(1L)
  par <- par + stats::rnorm(length(par), 0, 0.05)

  Sinv <- diag(c(rep(1 / 0.4^2, d$P_occ), rep(1 / 0.4^2, d$P_p),
                 rep(1 / 0.3^2, d$P_pos)))
  Pmu  <- diag(1 / 25, d$P)
  inv_sdL2 <- 1 / 1.0^2
  tau_w <- c(1.3, 0.9)                      # distinct per-factor precisions

  out <- tulpaObs:::.ms_ocs_penll_grad(model, par, Sinv, Pmu, inv_sdL2, tau_w,
                                       grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_penll_grad(model, p, Sinv, Pmu, inv_sdL2,
                                                 tau_w, grad = FALSE)$ll
  h <- 1e-5
  gnum <- numeric(length(par))
  for (k in seq_along(par)) {
    pp <- par; pp[k] <- pp[k] + h
    pm <- par; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the full K=2 packed gradient")
})

test_that("inner mode-find recovers the latent field + loadings at the true hyperparameters", {
  skip_on_cran()
  adj <- .mscs_grid_adj(8L, 8L)            # N = 64 cells
  S <- 16L                                 # more species -> sharper shared factor
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 6L,
                                        sd_load = 1.2, sigma_pos = 0.4,
                                        seed = 2024L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj)

  tr <- sim$truth
  # Condition on the data-generating community covariances + a unit field
  # precision (the marginal-scale convention) and unit loading SD.
  Sigma <- list(occ = diag(tr$sd_occ^2), p = diag(tr$sd_p^2),
                pos = diag(tr$sd_pos^2))
  fit <- tulpaObs:::.ms_ocs_inner_mode(model, Sigma, sd_L = 1.2, tau_w = 1.0)

  expect_identical(fit$convergence, 0L)

  # The shared factor and the loadings are recovered up to the K=1 sign anchor
  # (already applied), which the heterogeneous loadings make identifiable.
  # Conditional-recovery milestone: a latent ICAR field recovered to cor > 0.75
  # from binary occupancy + imperfect detection at N=64, and the loadings to
  # cor > 0.8.
  expect_gt(stats::cor(fit$w, tr$w), 0.75)
  expect_gt(stats::cor(fit$L, tr$L), 0.8)

  # Community means recovered (occupancy + cover intercepts and the cover scale).
  mu_occ_hat <- fit$mu[fit$d$occ_idx]
  mu_pos_hat <- fit$mu[fit$d$pos_idx]
  expect_lt(abs(mu_occ_hat[1L] - tr$mu_occ[1L]), 0.4)
  expect_lt(abs(mu_pos_hat[1L] - tr$mu_pos[1L]), 0.3)
  expect_lt(abs(exp(fit$ld) - tr$sigma_pos), 0.15)
})

test_that("Laplace-EM recovers the factor, loadings and community scales", {
  skip_on_cran()
  skip_if_fast()
  # Unconditional fit: Sigma and the field precision tau_w are estimated by the
  # EM (M-step), not supplied. The milestone is that the full Stage-1 fitter
  # recovers the shared factor, the loadings, the community means, and the
  # community covariance scales from data alone.
  adj <- .mscs_grid_adj(8L, 8L)            # N = 64 cells
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 16L, J = 6L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        sigma_pos = 0.4, seed = 4040L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj)

  fit <- tulpaObs:::.tobs_fit_ms_occu_cover_spatial(model, sd_L = 1.2,
                                                    max.em = 25L, tol = 1e-3)
  tr <- sim$truth

  # Latent structure recovered (up to the applied K=1 sign anchor).
  expect_gt(stats::cor(fit$w, tr$w), 0.75)
  expect_gt(stats::cor(fit$L, tr$L), 0.8)

  # Community means.
  expect_lt(abs(fit$mu[fit$d$occ_idx][1L] - tr$mu_occ[1L]), 0.5)
  expect_lt(abs(fit$mu[fit$d$pos_idx][1L] - tr$mu_pos[1L]), 0.3)
  expect_lt(abs(exp(fit$ld) - tr$sigma_pos), 0.15)

  # Community RE scale on the occupancy slope recovered to the right ballpark
  # (EM M-step), and the field precision is finite + positive.
  sd_occ_slope_hat <- sqrt(fit$Sigma$occ[2L, 2L])
  expect_gt(sd_occ_slope_hat, 0.2)
  expect_lt(sd_occ_slope_hat, 1.2)
  expect_true(is.finite(fit$tau_w) && fit$tau_w > 0)
})

test_that("Laplace-EM recovers the rank-2 spatial structure (K = 2)", {
  skip_on_cran()
  skip_if_fast()
  # K = 2 unconditional fit. The factors are fit unconstrained, so each factor is
  # only identified up to rotation / sign; the recovery measure is the
  # rotation-invariant spatial occupancy contribution F = W L' (N x S), the
  # quantity that enters every species' occupancy predictor.
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L; K <- 2L
  for (seed in c(2024L, 77L)) {
    sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
                                          sd_occ = 0.5, sd_load = 1.2,
                                          sigma_pos = 0.4, seed = seed)
    model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
      occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
      data = sim$data, y = sim$y, y_pos = sim$y_pos,
      positive = "lognormal", species = sim$species, adj = adj, K = K)
    fit <- tulpaObs:::.tobs_fit_ms_occu_cover_spatial(model, sd_L = 1.2,
                                                      max.em = 30L, tol = 1e-3)

    expect_identical(dim(fit$w), c(N, K))   # N x K field matrix
    expect_identical(dim(fit$L), c(S, K))   # S x K loading matrix

    F_hat  <- fit$w %*% t(fit$L)
    F_true <- sim$truth$w %*% t(sim$truth$L)
    expect_gt(stats::cor(as.numeric(F_hat), as.numeric(F_true)), 0.6)

    # The second factor is genuinely used (not collapsed to the rank-1 fit).
    expect_gt(stats::sd(fit$L[, 2L]), 0.3)
    expect_true(all(is.finite(fit$tau_w)) && length(fit$tau_w) == K)
  }
})

test_that("tobs() front door routes icar() on the occupancy arm to the spatial fit", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mscs_grid_adj(8L, 8L)            # N = 64 cells
  S <- 16L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 6L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        sigma_pos = 0.4, seed = 5151L)
  fit <- tobs(
    ~ occ_cov1 + icar(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species,
    method = "laplace",
    control = list(sd.load = 1.2, max.iter = 25L, tol = 1e-3))
  tr <- sim$truth

  # The public object is a tobs_fit carrying the shared-factor block.
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$spatial$type, "icar")
  expect_identical(fit$spatial$K, 1L)
  expect_length(fit$spatial$field, nrow(adj))
  expect_length(fit$spatial$loadings, S)
  expect_true(is.finite(fit$spatial$tau_w) && fit$spatial$tau_w > 0)

  # Same recovery as the internal fitter (field + loadings up to the sign
  # anchor), now reached through the user-facing front door.
  expect_gt(stats::cor(fit$spatial$field, tr$w), 0.75)
  expect_gt(stats::cor(as.numeric(fit$spatial$loadings), tr$L), 0.8)

  # Generic accessors work: coef() splits into the per-arm community means.
  cf <- coef(fit)
  expect_true(is.list(cf))
  psi_int <- cf[[1L]][["(Intercept)"]]
  expect_lt(abs(unname(psi_int) - tr$mu_occ[1L]), 0.5)

  # Community RE scale is reported and finite.
  expect_true(is.finite(fit$ms_community$sd_occ[1L]))
})

test_that("tobs() front door fits K > 1 via control$n.factors", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        sigma_pos = 0.4, seed = 77L)
  fit <- tobs(
    ~ occ_cov1 + icar(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species,
    method = "laplace",
    control = list(n.factors = K, sd.load = 1.2, max.iter = 30L, tol = 1e-3))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$spatial$K, K)
  expect_identical(dim(fit$spatial$field), c(N, K))
  expect_identical(dim(fit$spatial$loadings), c(S, K))

  F_hat  <- fit$spatial$field %*% t(fit$spatial$loadings)
  F_true <- sim$truth$w %*% t(sim$truth$L)
  expect_gt(stats::cor(as.numeric(F_hat), as.numeric(F_true)), 0.6)
})

test_that("a structured term off the occupancy arm is rejected (Stage 1)", {
  adj <- .mscs_grid_adj(5L, 5L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 4L, J = 3L, seed = 3L)
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1 + icar(graph = adj), positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace"),
    "occupancy arm only")
  expect_error(
    tobs(~ occ_cov1 + bym2(graph = adj), data = sim$data,
         family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
         positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
         species = sim$species, method = "laplace"),
    "icar")
})
