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

test_that("constrained (triangular) penalised gradient matches FD at K = 2", {
  adj <- .mscs_grid_adj(4L, 4L)
  S <- 4L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                        seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K)
  d <- tulpaObs:::.ms_ocs_dims(model)
  tr <- sim$truth

  # Pack a constrained par: mu, b, lfree (triangular truth -> log-diagonal), W, ld.
  mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
  lfree <- tulpaObs:::.ms_ocs_L_to_lfree(tr$L, S, K)
  expect_length(lfree, tulpaObs:::.ms_ocs_lfree_dim(S, K))
  par_c <- c(mu, b, lfree, as.numeric(tr$w), log(tr$sigma_pos))
  set.seed(2L); par_c <- par_c + stats::rnorm(length(par_c), 0, 0.05)

  Sinv <- diag(c(rep(1 / 0.4^2, d$P_occ), rep(1 / 0.4^2, d$P_p),
                 rep(1 / 0.3^2, d$P_pos)))
  Pmu  <- diag(1 / 25, d$P); inv_sdL2 <- 1; tau_w <- c(1.3, 0.9)

  out <- tulpaObs:::.ms_ocs_penll_grad_c(model, par_c, Sinv, Pmu, inv_sdL2,
                                         tau_w, grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_penll_grad_c(model, p, Sinv, Pmu, inv_sdL2,
                                                   tau_w, grad = FALSE)$ll
  h <- 1e-5; gnum <- numeric(length(par_c))
  for (k in seq_along(par_c)) {
    pp <- par_c; pp[k] <- pp[k] + h
    pm <- par_c; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the constrained K=2 gradient")

  # The triangular map round-trips and respects the structural zeros.
  Lr <- tulpaObs:::.ms_ocs_lfree_to_L(lfree, S, K)
  expect_equal(Lr, tr$L, tolerance = 1e-8)
  expect_true(all(Lr[upper.tri(Lr)] == 0))
  expect_true(all(diag(Lr) > 0))
})

test_that("joint NUTS log-posterior gradient matches FD (K = 1 unconstrained)", {
  adj <- .mscs_grid_adj(4L, 4L)            # N = 16 cells (small for FD)
  S <- 3L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 3L, seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj)
  d <- tulpaObs:::.ms_ocs_dims(model)
  tr <- sim$truth

  # Pack the full NUTS coordinate vector from truth: inner par, then the three
  # arm Cholesky blocks (from a chosen Sigma) and log tau_w. FD validity is
  # point-independent, so a perturbed truth point exercises every block.
  mu  <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b   <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
  par_inner <- c(mu, b, tr$L, tr$w, log(tr$sigma_pos))
  chol_v <- function(Sig) tulpaObs:::.ms_ocs_chol_pack(t(chol(Sig)))
  Sig <- list(occ = diag(0.4^2, d$P_occ), p = diag(0.35^2, d$P_p),
              pos = diag(0.3^2, d$P_pos))
  theta <- c(par_inner, chol_v(Sig$occ), chol_v(Sig$p), chol_v(Sig$pos),
             log(c(1.3)))
  lay <- tulpaObs:::.ms_ocs_nuts_layout(d, FALSE)
  expect_identical(length(theta), lay$total)
  set.seed(1L); theta <- theta + stats::rnorm(length(theta), 0, 0.05)

  out <- tulpaObs:::.ms_ocs_joint_logpost(model, theta, constrain = FALSE,
                                          grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_joint_logpost(model, p, constrain = FALSE,
                                                    grad = FALSE)$lp
  expect_true(is.finite(out$lp))
  h <- 1e-5; gnum <- numeric(length(theta))
  for (k in seq_along(theta)) {
    pp <- theta; pp[k] <- pp[k] + h
    pm <- theta; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the full joint NUTS gradient")
})

test_that("joint NUTS log-posterior gradient matches FD (K = 2 constrained + cover factor)", {
  adj <- .mscs_grid_adj(4L, 4L)
  S <- 4L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                        cover_factor = TRUE, seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K,
    cover_factor = TRUE)
  d <- tulpaObs:::.ms_ocs_dims(model)
  tr <- sim$truth

  mu  <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b   <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
  lfree <- tulpaObs:::.ms_ocs_L_to_lfree(tr$L, S, K)
  par_inner <- c(mu, b, lfree, as.numeric(tr$L_pos), as.numeric(tr$w),
                 log(tr$sigma_pos))
  chol_v <- function(Sig) tulpaObs:::.ms_ocs_chol_pack(t(chol(Sig)))
  theta <- c(par_inner,
             chol_v(diag(0.4^2, d$P_occ)), chol_v(diag(0.35^2, d$P_p)),
             chol_v(diag(0.3^2, d$P_pos)), log(c(1.3, 0.9)))
  lay <- tulpaObs:::.ms_ocs_nuts_layout(d, TRUE)
  expect_identical(length(theta), lay$total)
  set.seed(2L); theta <- theta + stats::rnorm(length(theta), 0, 0.05)

  out <- tulpaObs:::.ms_ocs_joint_logpost(model, theta, constrain = TRUE,
                                          grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_joint_logpost(model, p, constrain = TRUE,
                                                    grad = FALSE)$lp
  expect_true(is.finite(out$lp))
  h <- 1e-5; gnum <- numeric(length(theta))
  for (k in seq_along(theta)) {
    pp <- theta; pp[k] <- pp[k] + h
    pm <- theta; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the constrained cover-factor joint gradient")
})

test_that("joint NUTS log-posterior gradient matches FD on the field-hyper axis (car_proper / bym2)", {
  for (cfg in list(list(field = "car_proper", h = 0.85),
                   list(field = "bym2",       h = 0.7))) {
    adj <- .mscs_grid_adj(4L, 4L)
    S <- 3L
    sim <- simulate_ms_occu_cover_spatial(
      adj, n_species = S, J = 3L, field = cfg$field,
      rho = cfg$h, phi = cfg$h, seed = 321L)
    model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
      occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
      data = sim$data, y = sim$y, y_pos = sim$y_pos,
      positive = "lognormal", species = sim$species, adj = adj,
      field_type = cfg$field)
    d <- tulpaObs:::.ms_ocs_dims(model)
    expect_true(isTRUE(model$field_spec$has_hyper))
    tr <- sim$truth
    h_true <- if (identical(cfg$field, "car_proper")) tr$rho else tr$phi

    mu  <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
    b   <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
    par_inner <- c(mu, b, tr$L, tr$w, log(tr$sigma_pos))
    chol_v <- function(Sig) tulpaObs:::.ms_ocs_chol_pack(t(chol(Sig)))
    theta <- c(par_inner, chol_v(diag(0.4^2, d$P_occ)),
               chol_v(diag(0.35^2, d$P_p)), chol_v(diag(0.3^2, d$P_pos)),
               log(1.3), stats::qlogis(h_true))
    lay <- tulpaObs:::.ms_ocs_nuts_layout(d, FALSE, TRUE)
    expect_identical(length(theta), lay$total)
    expect_length(lay$logit_h, d$K)
    set.seed(3L); theta <- theta + stats::rnorm(length(theta), 0, 0.05)

    out <- tulpaObs:::.ms_ocs_joint_logpost(model, theta, constrain = FALSE,
                                            grad = TRUE)
    f <- function(p) tulpaObs:::.ms_ocs_joint_logpost(model, p, constrain = FALSE,
                                                      grad = FALSE)$lp
    expect_true(is.finite(out$lp))
    h <- 1e-5; gnum <- numeric(length(theta))
    for (k in seq_along(theta)) {
      pp <- theta; pp[k] <- pp[k] + h
      pm <- theta; pm[k] <- pm[k] - h
      gnum[k] <- (f(pp) - f(pm)) / (2 * h)
    }
    expect_lt(max(abs(out$grad - gnum)), 1e-4,
              label = sprintf("max|analytic - FD| joint gradient (%s)", cfg$field))
  }
})

test_that("C++ joint log-posterior + gradient matches the R target", {
  for (cfg in list(list(K = 1L, constrain = FALSE),
                   list(K = 2L, constrain = TRUE))) {
    adj <- .mscs_grid_adj(4L, 4L); S <- 4L; K <- cfg$K
    sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                          seed = 321L)
    model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
      occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
      data = sim$data, y = sim$y, y_pos = sim$y_pos,
      positive = "lognormal", species = sim$species, adj = adj, K = K)
    d <- tulpaObs:::.ms_ocs_dims(model)
    tr <- sim$truth
    Lmat <- if (K == 1L) matrix(tr$L, S, 1L) else tr$L
    Lblock <- if (cfg$constrain) tulpaObs:::.ms_ocs_L_to_lfree(Lmat, S, K)
              else as.numeric(Lmat)
    mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
    b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
    par_inner <- c(mu, b, Lblock, as.numeric(tr$w), log(tr$sigma_pos))
    chol_v <- function(Sig) tulpaObs:::.ms_ocs_chol_pack(t(chol(Sig)))
    theta <- c(par_inner, chol_v(diag(0.4^2, d$P_occ)),
               chol_v(diag(0.35^2, d$P_p)), chol_v(diag(0.3^2, d$P_pos)),
               log(rep(1.3, K)))
    set.seed(1L); theta <- theta + stats::rnorm(length(theta), 0, 0.05)

    spec <- tulpaObs:::.ms_ocs_nuts_spec(model)
    pri  <- tulpaObs:::.ms_ocs_nuts_priors()
    cpp <- tulpaObs:::cpp_ms_ocs_joint_logpost(spec, theta, pri, 5, 1.0,
                                               cfg$constrain)
    rr  <- tulpaObs:::.ms_ocs_joint_logpost(model, theta, priors = pri,
                                            constrain = cfg$constrain,
                                            sigma.beta = 5, sd_L = 1.0, grad = TRUE)
    expect_lt(abs(cpp$lp - rr$lp), 1e-8)
    expect_lt(max(abs(cpp$grad - rr$grad)), 1e-7)
  }
})

test_that("NUTS samples the spatial-factor community target and recovers the means", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mscs_grid_adj(5L, 5L); S <- 5L; K <- 1L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L, K = K,
          sd_occ = 0.5, sd_load = 1.1, sigma_pos = 0.4, seed = 3L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K)
  fit <- tulpaObs:::.tobs_fit_ms_occu_cover_spatial(model, sd_L = 1.1,
                                                    max.em = 20L, constrain = FALSE)
  d <- tulpaObs:::.ms_ocs_dims(model)
  theta0 <- tulpaObs:::.ms_ocs_nuts_pack_init(fit)
  spec <- tulpaObs:::.ms_ocs_nuts_spec(model)
  pri  <- tulpaObs:::.ms_ocs_nuts_priors()

  # Laplace-metric warm-start: inverse mass = posterior variance from the FD
  # Hessian diagonal of the joint log-posterior at the mode.
  g_at <- function(th) tulpaObs:::cpp_ms_ocs_joint_logpost(spec, th, pri, 5, 1.1,
                                                           FALSE)$grad
  hh <- 1e-4; np <- length(theta0); Md <- numeric(np)
  for (j in seq_len(np)) {
    tp <- theta0; tp[j] <- tp[j] + hh; tm <- theta0; tm[j] <- tm[j] - hh
    Md[j] <- -(g_at(tp)[j] - g_at(tm)[j]) / (2 * hh)
  }
  inv_metric <- 1 / pmax(Md, 1e-3)

  res <- tulpaObs:::cpp_ms_ocs_nuts(spec, theta0, pri, 5, 1.1, inv_metric,
                                    n_iter = 500L, n_warmup = 250L,
                                    max_treedepth = 10L, adapt_delta = 0.95,
                                    seed = 42L, verbose = FALSE, constrain = FALSE)
  expect_identical(nrow(res$draws), 250L)
  expect_true(all(is.finite(res$draws)))
  mu_post <- colMeans(res$draws[, seq_len(d$P)])
  mu_true <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
  expect_gt(stats::cor(mu_post, mu_true), 0.7)
  expect_lt(sum(res$divergent) / nrow(res$draws), 0.2)
})

test_that("tobs(method = 'nuts') fits the spatial-factor community occu_cover", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mscs_grid_adj(5L, 5L); S <- 5L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L, K = 1L,
          sd_occ = 0.5, sd_load = 1.1, sigma_pos = 0.4, seed = 3L)
  fit <- tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
              species = sim$species, method = "nuts",
              control = list(n.factors = 1L, sd.load = 1.1, n.iter = 400L,
                             n.warmup = 200L, adapt.delta = 0.95, seed = 42L))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nuts")
  expect_false(is.null(fit$nuts$draws))
  expect_identical(nrow(fit$nuts$draws), 400L)
  expect_false(is.null(fit$spatial$maps))
  mu_true <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
  cm <- fit$means[seq_along(mu_true)]
  expect_gt(stats::cor(cm, mu_true), 0.7)
  expect_lt(sum(fit$divergent) / length(fit$divergent), 0.2)

  # NUTS rejects auto-K (a Laplace-evidence procedure) and the non-spatial path.
  expect_error(
    tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
         family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
         positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
         species = sim$species, method = "nuts",
         control = list(n.factors = "auto")),
    "explicit n.factors")
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1, y = sim$y,
         y_pos = sim$y_pos, species = sim$species, method = "nuts"),
    "spatial-factor")
})

test_that("multi-chain NUTS reports split-R-hat / ESS and converges on the means", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mscs_grid_adj(5L, 5L); S <- 5L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L, K = 1L,
          sd_occ = 0.5, sd_load = 1.1, sigma_pos = 0.4, seed = 3L)
  fit <- tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
              species = sim$species, method = "nuts",
              control = list(n.factors = 1L, sd.load = 1.1, n.chains = 4L,
                             n.iter = 500L, n.warmup = 300L, adapt.delta = 0.95,
                             seed = 11L))
  nd <- fit$nuts
  expect_identical(nd$n_chains, 4L)
  expect_identical(nrow(nd$draws), 4L * 500L)
  expect_length(nd$rhat, ncol(nd$draws))
  expect_length(nd$ess,  ncol(nd$draws))
  expect_true(all(is.finite(nd$rhat)))
  # The community means should mix to R-hat ~ 1; ESS stays positive (the
  # occupancy intercept is autocorrelated through the ICAR field-level
  # confounding, so ESS there is modest even when R-hat is clean).
  d <- tulpaObs:::.ms_ocs_dims(fit$model)
  expect_lt(max(nd$rhat[seq_len(d$P)]), 1.2)
  expect_gt(min(nd$ess[seq_len(d$P)]), 10)
  expect_true(all(nd$ess > 0))
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

test_that("constrained (triangular) Laplace-EM recovers the rank-2 structure", {
  skip_on_cran()
  skip_if_fast()
  # The identifiability-constrained fit (lower-triangular L, positive
  # log-diagonal) must recover the same rotation-invariant spatial contribution
  # F = W L' as the unconstrained fit, while returning an identified, triangular
  # loading matrix (the proper-posterior parameterisation for per-factor
  # uncertainty / model comparison).
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        sigma_pos = 0.4, seed = 2024L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K)
  fit <- tulpaObs:::.tobs_fit_ms_occu_cover_spatial(model, sd_L = 1.2,
                                                    max.em = 30L, tol = 1e-3,
                                                    constrain = TRUE)
  expect_true(isTRUE(fit$constrained))
  expect_identical(dim(fit$L), c(S, K))

  # Loadings are in the lower-triangular canonical form with a positive diagonal.
  expect_true(all(fit$L[upper.tri(fit$L)] == 0))
  expect_true(all(diag(fit$L) > 0))

  F_hat  <- fit$w %*% t(fit$L)
  F_true <- sim$truth$w %*% t(sim$truth$L)
  expect_gt(stats::cor(as.numeric(F_hat), as.numeric(F_true)), 0.6)
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

  # Post-hoc varimax rotation: the spatial predictor F = W L' is invariant under
  # the orthogonal rotation, and the rotated loadings are at least as "simple"
  # (higher sum of squared loading variances) as the raw ones.
  rot <- fit$spatial$rotation
  expect_false(is.null(rot))
  F_rot <- rot$field %*% t(rot$loadings)
  expect_lt(max(abs(F_rot - F_hat)), 1e-8)
  simplicity <- function(M) sum(apply(M^2, 2L, stats::var))
  expect_gte(simplicity(rot$loadings) + 1e-9, simplicity(fit$spatial$loadings))
})

test_that("Laplace marginal likelihood recovers the true number of factors K", {
  skip_on_cran()
  skip_if_fast()
  # Rank selection by the empirical-Bayes Laplace marginal likelihood log Z(K).
  # The three latent-level criteria (held-out cells, unconstrained / constrained
  # WAIC) all pick K = 1 on K = 2 truth -- they track the field's effective
  # dimension, not the rank. The integrated criterion must (a) reject a spurious
  # second factor on K = 1 truth and (b) recover the second factor on K = 2 truth
  # at the scale where the recovery suite confirms it is estimable.
  build <- function(sim, adj) tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = 1L)

  # (a) K = 1 truth -> selected K = 1 (no spurious factor).
  adj1 <- .mscs_grid_adj(8L, 8L)
  sim1 <- simulate_ms_occu_cover_spatial(adj1, n_species = 16L, K = 1L, J = 6L,
                                         sd_occ = 0.5, sd_load = 1.2,
                                         sigma_pos = 0.4, seed = 11L)
  sel1 <- tulpaObs:::.ms_ocs_select_K(build(sim1, adj1), K.max = 3L, sd_L = 1.2,
                                      max.em = 30L, tol = 1e-4)
  expect_identical(sel1$K, 1L)
  expect_gt(sel1$table$logZ[1L], sel1$table$logZ[2L])   # log Z peaks at K = 1

  # (b) K = 2 truth -> selected K = 2 (the second factor clears the Occam budget).
  adj2 <- .mscs_grid_adj(9L, 9L)
  sim2 <- simulate_ms_occu_cover_spatial(adj2, n_species = 20L, K = 2L, J = 6L,
                                         sd_occ = 0.5, sd_load = 1.2,
                                         sigma_pos = 0.4, seed = 2024L)
  sel2 <- tulpaObs:::.ms_ocs_select_K(build(sim2, adj2), K.max = 3L, sd_L = 1.2,
                                      max.em = 40L, tol = 1e-4)
  expect_identical(sel2$K, 2L)
  expect_gt(sel2$table$logZ[2L], sel2$table$logZ[1L])
  expect_gt(sel2$table$logZ[2L], sel2$table$logZ[3L])   # unimodal: peak at K = 2
})

test_that("tobs() front door selects K via control$n.factors = 'auto'", {
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
    control = list(n.factors = "auto", n.factors.max = 3L, sd.load = 1.2,
                   max.iter = 40L, tol = 1e-4))

  # The auto fit is a tobs_fit at the selected rank, carrying the per-K evidence
  # table that drove the choice.
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$spatial$K, 2L)
  expect_false(is.null(fit$spatial$K_selection))
  ks <- fit$spatial$K_selection
  expect_true(all(c("K", "logZ", "best") %in% names(ks)))
  expect_identical(ks$K[which.max(ks$logZ)], 2L)
  expect_true(ks$best[ks$K == 2L])

  # The selected-rank fit recovers the rotation-invariant spatial contribution.
  F_hat  <- fit$spatial$field %*% t(fit$spatial$loadings)
  F_true <- sim$truth$w %*% t(sim$truth$L)
  expect_gt(stats::cor(as.numeric(F_hat), as.numeric(F_true)), 0.6)
})

test_that("a structured term on the detection arm, or an unsupported field, is rejected", {
  adj <- .mscs_grid_adj(5L, 5L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 4L, J = 3L, seed = 3L)
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1 + icar(graph = adj), positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace"),
    "detection arm")
  # car() (improper CAR) is not a supported field for this family (icar /
  # car_proper / bym2 are); it errors from the dispatcher.
  expect_error(
    tobs(~ occ_cov1 + car(graph = adj), data = sim$data,
         family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
         positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
         species = sim$species, method = "laplace"),
    "icar")
})

test_that("a cover-arm factor requires a matching shared field on the occupancy arm", {
  adj  <- .mscs_grid_adj(5L, 5L)
  adj2 <- .mscs_grid_adj(25L, 1L)        # same N, different graph
  sim  <- simulate_ms_occu_cover_spatial(adj, n_species = 4L, J = 3L, seed = 3L)
  # Cover icar() with a plain occupancy arm: the field is shared, so this errors.
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1 + icar(graph = adj),
         y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace"),
    "occupancy arm")
  # icar() on both arms but naming different graphs: the field must be one graph.
  expect_error(
    tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
         family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
         positive = ~ pos_cov1 + icar(graph = adj2),
         y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace"),
    "same graph")
})

test_that("simulate_ms_occu_cover_spatial cover factor is well-formed and RNG-gated", {
  adj <- .mscs_grid_adj(6L, 6L); S <- 8L; K <- 2L
  cf <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 4L,
                                       cover_factor = TRUE, seed = 11L)
  expect_true(isTRUE(cf$truth$cover_factor))
  expect_identical(dim(cf$truth$L_pos), c(S, K))

  # The cover-factor draws are gated: with cover_factor = FALSE (the default) the
  # detection / cover data is byte-identical to a call that never requests one, so
  # every Stage 1-2 fixture is unchanged.
  a <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L, seed = 11L)
  b <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L,
                                      cover_factor = FALSE, seed = 11L)
  expect_identical(a$y, b$y)
  expect_identical(a$y_pos, b$y_pos)
  expect_null(a$truth$L_pos)
})

test_that("tobs() front door recovers a shared factor on the cover arm", {
  skip_on_cran()
  skip_if_fast()
  # icar() on BOTH the occupancy and cover formulas shares one latent field: the
  # field loads on occupancy through Locc and on cover through a free Lpos. Both
  # spatial contributions -- F = W Locc' (occupancy) and Fpos = W Lpos' (cover) --
  # must be recovered, with the cover loadings carried in the same canonical sign
  # as the field (so Fpos has the right sign, gcol33/tulpa#67 Stage 3).
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 14L
  sim <- simulate_ms_occu_cover_spatial(
    adj, n_species = S, J = 6L, cover_factor = TRUE,
    mu_occ = c(0.3, 0.6), mu_pos = c(log(5), 0.3),
    sd_load = 1.0, sd_load_pos = 0.9, sigma_pos = 0.4, seed = 2024L)

  fit <- tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1 + icar(graph = adj),
              y = sim$y, y_pos = sim$y_pos, species = sim$species,
              method = "laplace", control = list(max.iter = 25L, tol = 1e-3))

  expect_identical(fit$spatial$type, "icar+cover")
  expect_true(isTRUE(fit$spatial$cover_factor))
  expect_length(fit$spatial$loadings_cover, S)

  W_hat <- as.numeric(fit$spatial$field)
  F_hat  <- outer(W_hat, as.numeric(fit$spatial$loadings))
  Fp_hat <- outer(W_hat, as.numeric(fit$spatial$loadings_cover))
  F_t    <- outer(as.numeric(sim$truth$w), as.numeric(sim$truth$L))
  Fp_t   <- outer(as.numeric(sim$truth$w), as.numeric(sim$truth$L_pos))

  expect_gt(stats::cor(as.numeric(F_hat),  as.numeric(F_t)),  0.75)
  expect_gt(stats::cor(as.numeric(Fp_hat), as.numeric(Fp_t)), 0.80)

  # The community cover-arm intercept is recovered.
  mu_pos1 <- fit$means[["pos_(Intercept)"]]
  expect_lt(abs(mu_pos1 - sim$truth$mu_pos[1L]), 0.4)
})

# ---------------------------------------------------------------------------
# Richer fields: proper-CAR factors (gcol33/tulpa#67 Stage 3)
# ---------------------------------------------------------------------------

test_that("penalised gradient matches FD on the proper-CAR field path", {
  # The car_proper field replaces the fixed ICAR structure with R(rho) = D - rho A
  # per factor; the objective reads it from model$field_R. The analytic gradient
  # must still match finite differences (the field enters only the W-block prior).
  adj <- .mscs_grid_adj(4L, 4L)            # N = 16 cells
  S <- 4L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                        field = "car_proper", rho = 0.85,
                                        seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K,
    field_type = "car_proper")
  # Per-factor structure at distinct correlations (FD validity is rho-independent).
  model$field_R <- tulpaObs:::.ms_ocs_build_field_R(model, hyper_w = c(0.7, 0.5))
  d <- tulpaObs:::.ms_ocs_dims(model)
  tr <- sim$truth
  mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
  Sinv <- diag(c(rep(1 / 0.4^2, d$P_occ), rep(1 / 0.4^2, d$P_p),
                 rep(1 / 0.3^2, d$P_pos)))
  Pmu <- diag(1 / 25, d$P); inv_sdL2 <- 1; tau_w <- c(1.3, 0.9)

  # Unconstrained packing.
  par <- c(mu, b, as.numeric(tr$L), as.numeric(tr$w), log(tr$sigma_pos))
  set.seed(1L); par <- par + stats::rnorm(length(par), 0, 0.05)
  out <- tulpaObs:::.ms_ocs_penll_grad(model, par, Sinv, Pmu, inv_sdL2, tau_w,
                                       grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_penll_grad(model, p, Sinv, Pmu, inv_sdL2,
                                                 tau_w, grad = FALSE)$ll
  h <- 1e-5; gnum <- numeric(length(par))
  for (k in seq_along(par)) {
    pp <- par; pp[k] <- pp[k] + h; pm <- par; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the car_proper K=2 gradient")

  # Constrained (triangular) packing.
  lfree <- tulpaObs:::.ms_ocs_L_to_lfree(tr$L, S, K)
  par_c <- c(mu, b, lfree, as.numeric(tr$w), log(tr$sigma_pos))
  set.seed(2L); par_c <- par_c + stats::rnorm(length(par_c), 0, 0.05)
  outc <- tulpaObs:::.ms_ocs_penll_grad_c(model, par_c, Sinv, Pmu, inv_sdL2,
                                          tau_w, grad = TRUE)
  fc <- function(p) tulpaObs:::.ms_ocs_penll_grad_c(model, p, Sinv, Pmu, inv_sdL2,
                                                    tau_w, grad = FALSE)$ll
  gnumc <- numeric(length(par_c))
  for (k in seq_along(par_c)) {
    pp <- par_c; pp[k] <- pp[k] + h; pm <- par_c; pm[k] <- pm[k] - h
    gnumc[k] <- (fc(pp) - fc(pm)) / (2 * h)
  }
  expect_lt(max(abs(outc$grad - gnumc)), 1e-4,
            label = "max|analytic - FD| over the constrained car_proper gradient")
})

test_that("simulate_ms_occu_cover_spatial proper-CAR field is well-formed and RNG-gated", {
  adj <- .mscs_grid_adj(7L, 7L); N <- nrow(adj); S <- 8L
  cr <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L,
                                       field = "car_proper", rho = 0.9, seed = 11L)
  expect_identical(cr$truth$field, "car_proper")
  expect_identical(cr$truth$rho, 0.9)
  w <- cr$truth$w
  expect_length(w, N)
  expect_gt(stats::sd(w), 0.2)            # genuine spatial amplitude
  expect_lt(stats::sd(w), 5)             # not blown up

  # The car branch is gated: field = "icar" (the default) reproduces the Stage 1-2
  # RNG stream byte for byte, so every existing fixture is unchanged.
  a <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L, seed = 11L)
  b <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L,
                                      field = "icar", seed = 11L)
  expect_identical(a$y, b$y)
  expect_identical(a$y_pos, b$y_pos)
  expect_identical(a$truth$w, b$truth$w)
  expect_null(a$truth$rho)
})

test_that("a car_proper field must be the same term on both spatial arms", {
  adj <- .mscs_grid_adj(5L, 5L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 4L, J = 3L, seed = 3L)
  # Mixed field types across the shared arms: the field is one GMRF, so its type
  # cannot differ between occupancy and cover.
  expect_error(
    tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
         family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
         positive = ~ pos_cov1 + car_proper(graph = adj),
         y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace"),
    "same term")
})

test_that("tobs() front door recovers a proper-CAR field and its correlation", {
  skip_on_cran()
  skip_if_fast()
  # car_proper() on the occupancy arm routes to the spatial fit with a proper CAR
  # field: the per-factor correlation rho is estimated by EM alongside tau, and
  # the field shape + spatial contribution F = W L' are recovered (gcol33/tulpa#67
  # Stage 3, richer fields).
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 14L
  rho_true <- 0.9
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 6L,
                                        sd_load = 1.2, field = "car_proper",
                                        rho = rho_true, seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + car_proper(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species,
    method = "laplace",
    control = list(n.factors = 1L, sd.load = 1.2, max.iter = 25L, tol = 1e-3))

  expect_identical(fit$spatial$type, "car_proper")
  expect_identical(fit$spatial$field_type, "car_proper")
  expect_length(fit$spatial$rho_w, 1L)
  # A proper, positive-dependence correlation in the right band (an empirical-Bayes
  # point estimate of a field correlation from occupancy data carries real
  # uncertainty, so the band is generous).
  expect_gt(fit$spatial$rho_w, 0.5)
  expect_lt(fit$spatial$rho_w, 1)
  expect_lt(abs(fit$spatial$rho_w - rho_true), 0.3)

  W_hat <- as.numeric(fit$spatial$field)
  F_hat <- outer(W_hat, as.numeric(fit$spatial$loadings))
  F_t   <- outer(as.numeric(sim$truth$w), as.numeric(sim$truth$L))
  expect_gt(abs(stats::cor(W_hat, as.numeric(sim$truth$w))), 0.75)
  expect_gt(stats::cor(as.numeric(F_hat), as.numeric(F_t)), 0.75)
})

test_that("auto-K rank selection composes with a proper-CAR field", {
  skip_on_cran()
  skip_if_fast()
  # The Laplace evidence integrates the field out with the full-rank car_proper
  # normaliser (|R(rho)| at rank N, not the rank-(N-1) ICAR pseudo-determinant);
  # auto-K must run that path end to end and return a fit at the selected rank.
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        field = "car_proper", rho = 0.9,
                                        seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + car_proper(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species,
    method = "laplace",
    control = list(n.factors = "auto", n.factors.max = 3L, sd.load = 1.2,
                   max.iter = 40L, tol = 1e-4))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$spatial$field_type, "car_proper")
  expect_true(is.data.frame(fit$spatial$K_selection))
  Ksel <- fit$spatial$K
  expect_true(Ksel >= 1L && Ksel <= 3L)
  expect_length(fit$spatial$rho_w, Ksel)
  expect_true(all(is.finite(fit$spatial$rho_w)))
})

# ---------------------------------------------------------------------------
# Richer fields: BYM2 factors (gcol33/tulpa#67 Stage 3)
# ---------------------------------------------------------------------------

test_that("penalised gradient matches FD on the BYM2 field path", {
  # BYM2's combined effect has marginal precision tau * R(phi) with
  # R(phi) = [(1-phi) I + phi Sigma_u]^{-1}; the objective reads it from
  # model$field_R exactly as for car_proper, so the analytic gradient must match
  # finite differences (the field enters only the W-block prior).
  adj <- .mscs_grid_adj(4L, 4L)            # N = 16 cells
  S <- 4L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                        field = "bym2", phi = 0.7, seed = 321L)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K,
    field_type = "bym2")
  model$field_R <- tulpaObs:::.ms_ocs_build_field_R(model, hyper_w = c(0.6, 0.4))
  d <- tulpaObs:::.ms_ocs_dims(model)
  tr <- sim$truth
  mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
  Sinv <- diag(c(rep(1 / 0.4^2, d$P_occ), rep(1 / 0.4^2, d$P_p),
                 rep(1 / 0.3^2, d$P_pos)))
  Pmu <- diag(1 / 25, d$P); inv_sdL2 <- 1; tau_w <- c(1.3, 0.9)

  par <- c(mu, b, as.numeric(tr$L), as.numeric(tr$w), log(tr$sigma_pos))
  set.seed(1L); par <- par + stats::rnorm(length(par), 0, 0.05)
  out <- tulpaObs:::.ms_ocs_penll_grad(model, par, Sinv, Pmu, inv_sdL2, tau_w,
                                       grad = TRUE)
  f <- function(p) tulpaObs:::.ms_ocs_penll_grad(model, p, Sinv, Pmu, inv_sdL2,
                                                 tau_w, grad = FALSE)$ll
  h <- 1e-5; gnum <- numeric(length(par))
  for (k in seq_along(par)) {
    pp <- par; pp[k] <- pp[k] + h; pm <- par; pm[k] <- pm[k] - h
    gnum[k] <- (f(pp) - f(pm)) / (2 * h)
  }
  expect_lt(max(abs(out$grad - gnum)), 1e-4,
            label = "max|analytic - FD| over the BYM2 K=2 gradient")

  lfree <- tulpaObs:::.ms_ocs_L_to_lfree(tr$L, S, K)
  par_c <- c(mu, b, lfree, as.numeric(tr$w), log(tr$sigma_pos))
  set.seed(2L); par_c <- par_c + stats::rnorm(length(par_c), 0, 0.05)
  outc <- tulpaObs:::.ms_ocs_penll_grad_c(model, par_c, Sinv, Pmu, inv_sdL2,
                                          tau_w, grad = TRUE)
  fc <- function(p) tulpaObs:::.ms_ocs_penll_grad_c(model, p, Sinv, Pmu, inv_sdL2,
                                                    tau_w, grad = FALSE)$ll
  gnumc <- numeric(length(par_c))
  for (k in seq_along(par_c)) {
    pp <- par_c; pp[k] <- pp[k] + h; pm <- par_c; pm[k] <- pm[k] - h
    gnumc[k] <- (fc(pp) - fc(pm)) / (2 * h)
  }
  expect_lt(max(abs(outc$grad - gnumc)), 1e-4,
            label = "max|analytic - FD| over the constrained BYM2 gradient")
})

test_that("simulate_ms_occu_cover_spatial BYM2 field is well-formed and RNG-gated", {
  adj <- .mscs_grid_adj(7L, 7L); N <- nrow(adj); S <- 8L
  by <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L,
                                       field = "bym2", phi = 0.7, seed = 11L)
  expect_identical(by$truth$field, "bym2")
  expect_identical(by$truth$phi, 0.7)
  w <- by$truth$w
  expect_length(w, N)
  expect_gt(stats::sd(w), 0.2)
  expect_lt(stats::sd(w), 5)

  # The bym2 branch is gated: field = "icar" (default) reproduces the Stage 1-2
  # RNG stream byte for byte.
  a <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L, seed = 11L)
  b <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 4L,
                                      field = "icar", seed = 11L)
  expect_identical(a$y, b$y)
  expect_identical(a$truth$w, b$truth$w)
  expect_null(a$truth$phi)
})

test_that("tobs() front door recovers a BYM2 field and its variance fraction", {
  skip_on_cran()
  skip_if_fast()
  # bym2() on the occupancy arm routes to the spatial fit with a BYM2 field: the
  # combined-effect shape + F = W L' are recovered, and the spatial-variance
  # fraction phi is a positive fraction. phi is a single-field-realisation
  # variance component, so its EM point estimate is high-variance and attenuated
  # at small N (it recovers toward truth as species / visits grow); the field
  # shape -- the quantity that matters -- recovers strongly throughout.
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L
  phi_true <- 0.7
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = 8L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        field = "bym2", phi = phi_true,
                                        seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + bym2(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species,
    method = "laplace",
    control = list(n.factors = 1L, sd.load = 1.2, max.iter = 25L, tol = 1e-3))

  expect_identical(fit$spatial$type, "bym2")
  expect_identical(fit$spatial$field_type, "bym2")
  expect_null(fit$spatial$rho_w)
  expect_length(fit$spatial$phi_w, 1L)
  expect_gt(fit$spatial$phi_w, 0.2)        # a real spatial fraction
  expect_lt(fit$spatial$phi_w, 1)
  expect_lt(abs(fit$spatial$phi_w - phi_true), 0.45)

  W_hat <- as.numeric(fit$spatial$field)
  F_hat <- outer(W_hat, as.numeric(fit$spatial$loadings))
  F_t   <- outer(as.numeric(sim$truth$w), as.numeric(sim$truth$L))
  expect_gt(abs(stats::cor(W_hat, as.numeric(sim$truth$w))), 0.75)
  expect_gt(stats::cor(as.numeric(F_hat), as.numeric(F_t)), 0.75)
})

test_that("auto-K rank selection composes with a BYM2 field", {
  skip_on_cran()
  skip_if_fast()
  # The Laplace evidence integrates the field out with the BYM2 full-rank
  # normaliser (|R(phi)| at rank N); auto-K must run that path end to end.
  adj <- .mscs_grid_adj(9L, 9L); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 8L,
                                        sd_occ = 0.5, sd_load = 1.2,
                                        field = "bym2", phi = 0.7, seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + bym2(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species,
    method = "laplace",
    control = list(n.factors = "auto", n.factors.max = 3L, sd.load = 1.2,
                   max.iter = 40L, tol = 1e-4))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$spatial$field_type, "bym2")
  expect_true(is.data.frame(fit$spatial$K_selection))
  Ksel <- fit$spatial$K
  expect_true(Ksel >= 1L && Ksel <= 3L)
  expect_length(fit$spatial$phi_w, Ksel)
  expect_true(all(is.finite(fit$spatial$phi_w)))
})

test_that("tobs_associations() recovers the residual species associations", {
  skip_on_cran()
  skip_if_fast()
  # The K shared latent fields imply a residual species-association matrix (the
  # spatial-JSDM / HMSC output): with unit-variance fields the occupancy
  # association is corr(L L'), invariant to the factor rotation. With a cover-arm
  # factor the cross-arm association (standardized L_occ L_pos') links a species'
  # spatial occupancy to another's spatial cover. The matrices are marginalised
  # over the loading posterior, so an interval accompanies the estimate. Recovery
  # is checked on the off-diagonals (the diagonal is 1 by construction).
  corr_self  <- function(L) { Om <- tcrossprod(L); s <- sqrt(diag(Om)); Om / outer(s, s) }
  corr_cross <- function(Lo, Lp) {
    so <- sqrt(rowSums(Lo^2)); sp <- sqrt(rowSums(Lp^2)); (Lo %*% t(Lp)) / outer(so, sp)
  }
  adj <- .mscs_grid_adj(9L, 9L); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
           sd_occ = 0.5, sd_load = 1.3, sigma_pos = 0.4, cover_factor = TRUE,
           mean_load_pos = 0, sd_load_pos = 1.2, seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + icar(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1 + icar(graph = adj),
    y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace",
    control = list(n.factors = K, sd.load = 1.3, max.iter = 30L, tol = 1e-3))

  occ <- tobs_associations(fit, "occupancy")
  expect_identical(dim(occ$median), c(S, S))
  expect_equal(unname(diag(occ$estimate)), rep(1, S), tolerance = 1e-8)
  expect_true(isSymmetric(unname(occ$median)))
  expect_true(all(occ$median >= -1 - 1e-8 & occ$median <= 1 + 1e-8))
  expect_identical(rownames(occ$median), sim$species)
  expect_true(all(occ$lower <= occ$upper + 1e-8))

  off <- upper.tri(matrix(0, S, S))
  occ_t <- corr_self(sim$truth$L)
  expect_gt(stats::cor(occ$median[off], occ_t[off]), 0.7)
  cov_occ <- mean(occ_t[off] >= occ$lower[off] & occ_t[off] <= occ$upper[off])
  expect_gte(cov_occ, 0.8)

  # Cross-arm association (occupancy vs cover); present only with a cover factor.
  cross <- tobs_associations(fit, "cross")
  expect_identical(dim(cross$median), c(S, S))
  offall  <- row(matrix(0, S, S)) != col(matrix(0, S, S))
  cross_t <- corr_cross(sim$truth$L, sim$truth$L_pos)
  expect_gt(stats::cor(cross$median[offall], cross_t[offall]), 0.7)

  # A single matrix is returned when a summary is named.
  expect_identical(dim(tobs_associations(fit, "occupancy", "median")), c(S, S))
})

test_that("tobs_associations() errors off the spatial-factor surface", {
  no_spatial <- structure(list(spatial = NULL), class = "tobs_fit")
  expect_error(tobs_associations(no_spatial), "spatial-factor")
  no_cover <- structure(
    list(spatial = list(associations = list(
      occupancy = list(median = matrix(1, 1, 1))))),
    class = "tobs_fit")
  expect_error(tobs_associations(no_cover, "cover"), "cover-arm")
})

test_that("predict() returns calibrated per-species occupancy maps", {
  skip_on_cran()
  skip_if_fast()
  # The latent fields are tied to the cell graph, so predict() returns the
  # in-sample per-species per-cell occupancy posterior: a long table (cell x
  # species) of psi with a credible interval, marginalised over the loading +
  # field posterior. The maps recover the true psi surface and the interval is
  # calibrated -- a rare species borrows strength across the shared factors.
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
           sd_occ = 0.6, sd_load = 1.2, sigma_pos = 0.4, seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + icar(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace",
    control = list(n.factors = K, sd.load = 1.2, max.iter = 30L, tol = 1e-3))

  pr <- predict(fit)
  expect_true(is.data.frame(pr))
  expect_identical(nrow(pr), N * S)
  expect_true(all(c("cell", "species", "psi", "psi_median",
                    "psi_lower", "psi_upper") %in% names(pr)))
  expect_true(all(pr$psi >= 0 & pr$psi <= 1))
  expect_true(all(pr$psi_lower <= pr$psi_upper + 1e-8))
  expect_identical(sort(unique(pr$species)), sort(sim$species))

  psi_t <- sim$truth$psi                      # N x S
  M  <- matrix(pr$psi,       N, S)
  Lo <- matrix(pr$psi_lower, N, S); Hi <- matrix(pr$psi_upper, N, S)
  expect_gt(stats::cor(as.numeric(M), as.numeric(psi_t)), 0.7)
  expect_gte(mean(psi_t >= Lo & psi_t <= Hi), 0.8)

  # New-data prediction is unsupported (no field at an unseen cell).
  expect_error(predict(fit, X.0 = matrix(1, 2, 2)), "not supported")
})

test_that("predict() returns calibrated per-species cover maps", {
  skip_on_cran()
  skip_if_fast()
  # The cover hurdle's spatial output: predict(type = "cover_cond") is the
  # conditional cover mean E[cover | present], "cover_exp" the unconditional
  # expected cover psi * E[cover | present]. With a cover-arm factor the cover
  # has its own spatial structure. Both are derived quantities (the nonlinear
  # lognormal mean, the psi product), marginalised per draw. Truth is
  # reconstructed from the fit's own cover design so it matches exactly.
  adj <- .mscs_grid_adj(9L, 9L); N <- nrow(adj); S <- 20L; K <- 2L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 6L,
           sd_occ = 0.6, sd_load = 1.2, sigma_pos = 0.4, cover_factor = TRUE,
           mean_load_pos = 0, sd_load_pos = 1.2, seed = 2024L)
  fit <- tobs(
    ~ occ_cov1 + icar(graph = adj), data = sim$data,
    family    = ms_occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1 + icar(graph = adj),
    y = sim$y, y_pos = sim$y_pos, species = sim$species, method = "laplace",
    control = list(n.factors = K, sd.load = 1.2, max.iter = 30L, tol = 1e-3))

  pc <- predict(fit, type = "cover_cond")
  pe <- predict(fit, type = "cover_exp")
  expect_identical(nrow(pc), N * S)
  expect_true(all(c("cover_cond", "cover_cond_lower", "cover_cond_upper") %in% names(pc)))
  expect_true(all(c("cover_exp", "cover_exp_lower", "cover_exp_upper") %in% names(pe)))
  expect_true(all(pc$cover_cond >= 0) && all(pe$cover_exp >= 0))
  expect_true(all(pe$cover_exp_lower <= pe$cover_exp_upper + 1e-8))

  # Truth from the fit's own cover design + field offset.
  Xp <- fit$model$X_pos_site
  Wt <- if (K == 1L) matrix(sim$truth$w, N, 1L) else sim$truth$w
  Lp <- if (K == 1L) matrix(sim$truth$L_pos, S, 1L) else sim$truth$L_pos
  ce_t <- matrix(0, N, S)
  for (s in seq_len(S)) {
    eta <- as.numeric(Xp %*% (sim$truth$mu_pos + sim$truth$b_pos[s, ])) +
           as.numeric(Wt %*% Lp[s, ])
    ce_t[, s] <- sim$truth$psi[, s] * exp(eta + sim$truth$sigma_pos^2 / 2)
  }
  CE  <- matrix(pe$cover_exp,       N, S)
  Lo  <- matrix(pe$cover_exp_lower, N, S); Hi <- matrix(pe$cover_exp_upper, N, S)
  expect_gt(stats::cor(as.numeric(CE), as.numeric(ce_t)), 0.8)
  expect_gte(mean(ce_t >= Lo & ce_t <= Hi), 0.8)

  expect_error(predict(fit, type = "detection"), "not available")
})

test_that("fitted() returns per-species occupancy / detection / cover surfaces", {
  skip_on_cran()
  skip_if_fast()
  # fitted() gives the per-species posterior-mean surfaces: field-augmented
  # occupancy psi and conditional cover (from the map posterior), and the
  # detection probability p (no field on detection). Same shape as the
  # non-spatial community fitted(); the occupancy surface matches predict().
  adj <- .mscs_grid_adj(7L, 7L); N <- nrow(adj); S <- 12L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = 2L, J = 5L,
           sd_occ = 0.5, sd_load = 1.1, sigma_pos = 0.4, seed = 3L)
  fit <- tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
              species = sim$species, method = "laplace",
              control = list(n.factors = 2L, sd.load = 1.1, max.iter = 25L))

  ft <- fitted(fit)
  expect_named(ft, c("psi", "p", "cover"))
  for (m in ft) expect_identical(dim(m), c(N, S))
  expect_true(all(ft$psi >= 0 & ft$psi <= 1))
  expect_true(all(ft$p   >= 0 & ft$p   <= 1))
  expect_true(all(ft$cover > 0))
  expect_identical(colnames(ft$psi), sim$species)
  expect_equal(as.numeric(ft$psi), predict(fit)$psi, tolerance = 1e-8)

  # nobs() counts the valid detection observations (feeds AIC / BIC).
  expect_identical(nobs(fit), sum(!is.na(sim$y)))
})

test_that("simulate() reproduces the per-species data structure", {
  skip_on_cran()
  skip_if_fast()
  # Plug-in (posterior-mean) simulation at the observed visit pattern: the fitted
  # model reproduces the data's per-species detection prevalence and cover
  # magnitude (the basis for posterior-predictive checks).
  adj <- .mscs_grid_adj(8L, 8L); N <- nrow(adj); S <- 14L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = 2L, J = 5L,
           sd_occ = 0.5, sd_load = 1.2, sigma_pos = 0.4, seed = 5L)
  fit <- tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
              species = sim$species, method = "laplace",
              control = list(n.factors = 2L, sd.load = 1.2, max.iter = 30L))

  s1 <- simulate(fit, nsim = 1, seed = 7)
  expect_named(s1, c("y", "y_pos"))
  expect_identical(dim(s1$y), dim(sim$y))
  expect_identical(dim(s1$y_pos), dim(sim$y_pos))
  expect_length(simulate(fit, nsim = 3, seed = 7), 3L)

  det_rate <- function(a) apply(a, 3L, function(m) mean(m, na.rm = TRUE))
  obs  <- det_rate(sim$y)
  sims <- simulate(fit, nsim = 30, seed = 7)
  simrate <- rowMeans(vapply(sims, function(ss) det_rate(ss$y), numeric(S)))
  expect_gt(stats::cor(obs, simrate), 0.8)
})
