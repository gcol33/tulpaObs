# Missing-at-random cover on occu_cover() (tulpaObs 0.0.104).
#
# A detected visit (y == 1) may carry a missing cover (y_pos = NA). The detection
# and occupancy arms keep the visit; only the Beta/lognormal cover factor drops
# out. These tests cover:
#   - the mechanism: the cell-coupling density drops the missing-cover visit and
#     zeroes its cover-arm derivatives, leaving psi / detection / the other cover
#     visits byte-identical (the fit path);
#   - the builder: NA at a detected visit is admitted and carried as the sentinel,
#     detection is retained, and a no-missing build introduces no NA;
#   - recovery: dropping cover at random does not bias the cover arm, and (with no
#     shared field) leaves the occupancy / detection estimates unchanged; the
#     spatial compact path (the deliverable) fits and stays close to the full-data
#     fit with NA cover present.

# Self-contained Beta cover log-density (mirrors .occu_cover_pos_logdens beta arm
# for un-clamped eta), so the mechanism test needs no engine internals.
.beta_cover_ld <- function(y, eta, phi) {
  mu <- plogis(eta); a <- mu * phi; b <- (1 - mu) * phi
  lgamma(phi) - lgamma(a) - lgamma(b) + (a - 1) * log(y) + (b - 1) * log(1 - y)
}
.lnorm_cover_ld <- function(y, eta, sigma) {
  -log(y) - log(sigma) - 0.5 * log(2 * pi) - 0.5 * ((log(y) - eta) / sigma)^2
}


test_that("cell coupling drops a detected visit with missing cover (beta)", {
  eta_psi <- 0.2
  eta_p   <- c(-0.3, 0.4)
  eta_pos <- c(0.1, -0.2)
  phi     <- 20
  y_det   <- c(1L, 1L)

  ref <- tulpaObs:::cpp_eval_occu_cover_beta_cell(eta_psi, eta_p, eta_pos,
                                       y_det, c(0.30, 0.55), phi)
  mar <- tulpaObs:::cpp_eval_occu_cover_beta_cell(eta_psi, eta_p, eta_pos,
                                       y_det, c(0.30, NA),   phi)

  # The only change is visit 2's cover factor dropping out.
  f2 <- .beta_cover_ld(0.55, eta_pos[2], phi)
  expect_equal(mar$cell_ll, ref$cell_ll - f2, tolerance = 1e-10)

  # Visit 2's cover grad / Hessian are zeroed; every other block is untouched.
  expect_equal(mar$grad_pos[2], 0)
  expect_equal(mar$neg_hess_pos[2], 0)
  expect_equal(mar$grad_pos[1], ref$grad_pos[1], tolerance = 1e-12)
  expect_equal(mar$grad_p,   ref$grad_p,   tolerance = 1e-12)
  expect_equal(mar$grad_psi, ref$grad_psi, tolerance = 1e-12)
})


test_that("cell coupling drops a detected visit with missing cover (lognormal)", {
  eta_psi <- -0.1
  eta_p   <- c(0.2, -0.4, 0.5)
  eta_pos <- c(1.1, 0.9, 1.3)
  sigma   <- 0.35
  y_det   <- c(1L, 1L, 1L)
  y_obs   <- c(2.1, 3.0, 2.5)

  ref <- tulpaObs:::cpp_eval_occu_cover_lognormal_cell(eta_psi, eta_p, eta_pos,
                                            y_det, y_obs, sigma)
  y_mar <- y_obs; y_mar[2] <- NA
  mar <- tulpaObs:::cpp_eval_occu_cover_lognormal_cell(eta_psi, eta_p, eta_pos,
                                            y_det, y_mar, sigma)

  f2 <- .lnorm_cover_ld(y_obs[2], eta_pos[2], sigma)
  expect_equal(mar$cell_ll, ref$cell_ll - f2, tolerance = 1e-10)
  expect_equal(mar$grad_pos[2], 0)
  expect_equal(mar$neg_hess_pos[2], 0)
  expect_equal(mar$grad_pos[c(1, 3)], ref$grad_pos[c(1, 3)], tolerance = 1e-12)
  expect_equal(mar$grad_p,   ref$grad_p,   tolerance = 1e-12)
  expect_equal(mar$grad_psi, ref$grad_psi, tolerance = 1e-12)
})


test_that("dense builder admits NA cover at a detected visit and keeps detection", {
  # site 1: (detected+observed, undetected); site 2: (detected+observed,
  # detected+missing-cover).
  y     <- matrix(c(1L, 0L, 1L, 1L), 2, 2, byrow = TRUE)
  y_pos <- matrix(c(0.30, 0, 0.50, NA), 2, 2, byrow = TRUE)
  data  <- data.frame(x = c(0.1, -0.2))

  m <- tulpaObs:::.tobs_build_occu_cover(~ 1, ~ 1, ~ 1, data, y, y_pos, "beta")

  # Detection is retained at the missing-cover visit.
  expect_equal(m$y[2, 2], 1L)
  # Observed covers pass through; the missing detected cover is the NA sentinel;
  # the undetected cell is the plain zero-fill.
  expect_equal(m$y_pos[1, 1], 0.30)
  expect_equal(m$y_pos[2, 1], 0.50)
  expect_true(is.na(m$y_pos[2, 2]))
  expect_equal(m$y_pos[1, 2], 0)

  # No missing cover -> no NA introduced (the historical zero-fill).
  y_pos2 <- matrix(c(0.30, 0, 0.50, 0.70), 2, 2, byrow = TRUE)
  m2 <- tulpaObs:::.tobs_build_occu_cover(~ 1, ~ 1, ~ 1, data, y, y_pos2, "beta")
  expect_false(anyNA(m2$y_pos))
})


test_that("builder range check fires on an observed out-of-range cover, not on NA", {
  y    <- matrix(c(1L, 1L), 1, 2)
  data <- data.frame(x = 0)
  # NA at a detected visit is fine; an observed cover of 1.4 is not (beta arm).
  expect_error(
    tulpaObs:::.tobs_build_occu_cover(~ 1, ~ 1, ~ 1, data,
                                      y, matrix(c(0.3, 1.4), 1, 2), "beta"),
    "0 < y_pos < 1")
  expect_silent(
    tulpaObs:::.tobs_build_occu_cover(~ 1, ~ 1, ~ 1, data,
                                      y, matrix(c(0.3, NA), 1, 2), "beta"))
})


test_that("MAR cover is unbiased and leaves occupancy / detection unchanged (non-spatial)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 10L
  N <- 300L; J <- 5L
  beta_occ_truth <- c(stats::qlogis(0.4), 0.9)
  beta_p_truth   <- c(0.0, 0.6)
  beta_pos_truth <- c(log(0.10), -0.4)
  sigma_pos_truth <- 0.4

  pos_int_full <- pos_int_mar <- numeric(n_seeds)
  pos_x_full   <- pos_x_mar   <- numeric(n_seeds)
  sig_mar      <- numeric(n_seeds)
  psi_gap <- p_gap <- numeric(n_seeds)
  ok <- logical(n_seeds)

  fit_one <- function(od, cell_dat, y_pos) {
    tryCatch(
      tobs(formula   = ~ occ_cov1, data = cell_dat,
           family    = occu_cover("lognormal"),
           detection = ~ det_cov1, positive = ~ pos_cov1,
           y         = od$y, y_pos = y_pos, visits = od$det.covs,
           method    = "laplace",
           control   = list(verbose = FALSE, max.iter = 500L)),
      error = function(e) NULL)
  }

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, sigma_pos = sigma_pos_truth,
      positive = "lognormal", seed = 4200L + s)

    long <- data.frame(
      site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
      y = as.vector(t(sim$y)),
      det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1)
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)

    y_pos_full <- sim$y_pos; y_pos_full[is.na(y_pos_full)] <- 0
    # Knock out 30% of detected covers at random (missing-at-random).
    set.seed(50L + s)
    det_cells <- which(sim$y == 1L)
    miss <- sample(det_cells, floor(0.30 * length(det_cells)))
    y_pos_mar <- y_pos_full; y_pos_mar[miss] <- NA

    f_full <- fit_one(od, cell_dat, y_pos_full)
    f_mar  <- fit_one(od, cell_dat, y_pos_mar)
    if (is.null(f_full) || is.null(f_mar) ||
        !isTRUE(f_full$convergence$converged) ||
        !isTRUE(f_mar$convergence$converged)) next
    ok[s] <- TRUE

    pos_int_full[s] <- f_full$means["pos_(Intercept)"]
    pos_int_mar[s]  <- f_mar$means["pos_(Intercept)"]
    pos_x_full[s]   <- f_full$means["pos_pos_cov1"]
    pos_x_mar[s]    <- f_mar$means["pos_pos_cov1"]
    sig_mar[s]      <- exp(f_mar$means["log_sigma_pos"])

    # No shared field: the occupancy / detection likelihood does not involve
    # cover, so dropping cover leaves those estimates unchanged.
    psi_gap[s] <- max(abs(f_mar$means[c("psi_(Intercept)", "psi_occ_cov1")] -
                          f_full$means[c("psi_(Intercept)", "psi_occ_cov1")]))
    p_gap[s]   <- max(abs(f_mar$means[c("p_(Intercept)", "p_det_cov1")] -
                          f_full$means[c("p_(Intercept)", "p_det_cov1")]))
  }

  expect_gte(mean(ok), 0.8)

  # MAR cover arm recovers truth (unbiased).
  expect_lt(abs(mean(pos_int_mar[ok]) - beta_pos_truth[1L]), 0.20)
  expect_lt(abs(mean(pos_x_mar[ok])   - beta_pos_truth[2L]), 0.20)
  expect_lt(abs(mean(sig_mar[ok])     - sigma_pos_truth),    0.10)

  # Occupancy / detection are invariant to cover-drop (exact factorisation).
  expect_lt(max(psi_gap[ok]), 1e-3)
  expect_lt(max(p_gap[ok]),   1e-3)
})


test_that("MAR cover fits and stays close to full on the spatial compact path", {
  skip_on_cran()
  skip_if_fast()

  N <- 40L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "beta",
                             adj = adj, sigma = 0.8, alpha = 1.0, seed = 321L)

  long <- data.frame(
    site_key = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    cell_idx = rep(seq_len(N), each = J),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1)

  cover_full <- as.vector(t(sim$y_pos)); cover_full[is.na(cover_full)] <- 0
  set.seed(9L)
  det_cells <- which(long$y == 1L)
  miss <- sample(det_cells, floor(0.30 * length(det_cells)))
  cover_mar <- cover_full; cover_mar[miss] <- NA

  fit_compact <- function(cover_vec) {
    dd <- cbind(long, cover.flat = cover_vec)
    od <- tobs_data(dd, y = "y", site = "site_key", visit = "visit",
                    type = "occurrence", occ.covs = "cell_idx",
                    det.covs = c("det_cov1", "pos_cov1"), compact = TRUE)
    ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                    visit = "visit", type = "cover", compact = TRUE))
    ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
                 max.iter = 200L, sigma.grid = c(0.5, 1.5),
                 phi.grid.pos = c(2, 10), integration = "grid",
                 adaptive.grid = FALSE, verbose = FALSE, progress = FALSE)
    suppressWarnings(tobs(
      occurrence = ~ spatial(~ 1 || cell_idx, graph = adj),
      data = od$occ.covs,
      family = occu_cover(response = "beta", cover_aggregate = "none"),
      detection = ~ det_cov1, positive = ~ pos_cov1 + copy(spatial()),
      y = od$y, y_pos = ocv$y, visits = od$det.covs,
      method = "nested_laplace", control = ctrl))
  }

  f_full <- fit_compact(cover_full)
  f_mar  <- fit_compact(cover_mar)

  # Both fits complete with finite summaries even with NA cover present.
  expect_true(all(is.finite(f_mar$means)))
  expect_true(all(is.finite(f_mar$sds)))

  # The cover-arm coefficient shifts only within Monte-Carlo error of the
  # full-data fit (dropping 30% of covers widens, does not bias).
  pos_nm <- grep("^pos_", names(f_mar$means), value = TRUE)
  expect_lt(max(abs(f_mar$means[pos_nm] - f_full$means[pos_nm])), 0.5)

  # WAIC is defined with NA cover in the data (the pointwise ll drops the term).
  w <- suppressWarnings(tobs_waic(f_mar))
  expect_true(is.finite(w$waic))
})
