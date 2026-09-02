# Tests for cover(response = "lognormal_trunc") -- the upper-truncated lognormal
# (Gaussian on log-cover truncated at log(1) = 0, so cover <= 1) positive arm on the
# joint nested-Laplace path (tulpa truncated_gaussian family).

test_that("cover(response = 'lognormal_trunc') constructs", {
  fam <- cover(response = "lognormal_trunc")
  expect_s3_class(fam, "tobs_family")
  expect_equal(fam$params$positive, "lognormal_trunc")
  expect_equal(fam$observation, "binomial_plus_lognormal_trunc")
  expect_equal(fam$status, "working")
  expect_error(cover(response = "lognormal_trunc", breaks = 0.3), "only used for")
})

test_that("encode stores the log-cover truncation ceiling (log 1 = 0)", {
  dat <- data.frame(x = rnorm(8))
  y   <- c(0, 0.2, 0, 0.8, 0.05, 0.5, 0, 0.95)
  enc <- encode_cover_hurdle(~ x, dat, y = y, positive = "lognormal_trunc")
  # one ceiling per positive plot, all 0 (= log of the cover ceiling 1)
  expect_equal(enc$pos_data$trunc_upper, rep(0, sum(y > 0)))
  # positive response is log-cover, like the lognormal arm
  expect_equal(enc$pos_data$y, log(y[y > 0]))
})

test_that("lognormal_trunc requires method = 'nested_laplace'", {
  dat <- data.frame(x = rnorm(40), cell = sample.int(4, 40, replace = TRUE),
                    y = ifelse(runif(40) < 0.5, 0, runif(40)^2))
  expect_error(
    tobs(y ~ x, data = dat, family = cover(response = "lognormal_trunc"),
         method = "laplace"),
    "requires method = 'nested_laplace'"
  )
})

# Self-contained spatial simulator: occurrence + cover share an areal field; the
# positive cover is an upper-truncated lognormal -- latent log-cover ~ N(eta, sigma)
# truncated to <= 0 by inverse-CDF (so cover in (0, 1] by construction, NOT clipped).
.sim_trunc <- function(seed, N = 5000L, gx = 8L, b_pos_int = -0.7, sigma = 0.8,
                       field_scale = 1) {
  set.seed(seed)
  nc  <- gx * gx
  adj <- matrix(0L, nc, nc); id <- function(i, j) (j - 1L) * gx + i
  for (i in 1:gx) for (j in 1:gx) for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next; ii <- i + di; jj <- j + dj
    if (ii >= 1 && ii <= gx && jj >= 1 && jj <= gx) adj[id(i, j), id(ii, jj)] <- 1L
  }
  fx    <- outer(1:gx, 1:gx, function(i, j) sin(i / 2) + cos(j / 2))
  field <- field_scale * as.numeric(scale(as.numeric(fx)))
  cell  <- sample.int(nc, N, replace = TRUE); x <- rnorm(N)
  b_occ <- c(-0.1, 0.7); b_pos <- c(b_pos_int, 0.5)
  occ     <- rbinom(N, 1L, plogis(b_occ[1] + b_occ[2] * x + field[cell]))
  eta_pos <- b_pos[1] + b_pos[2] * x + field[cell]
  a       <- (0 - eta_pos) / sigma                 # standardized ceiling (u = 0)
  uu      <- runif(N) * pnorm(a)                   # inverse-CDF truncated-normal
  log_cov <- eta_pos + sigma * qnorm(uu)           # <= 0 by construction
  y       <- ifelse(occ == 1L, exp(log_cov), 0)    # cover in (0, 1]
  list(data = data.frame(y = y, x = x, cell_idx = cell), adj = adj,
       truth = list(beta_occ = b_occ, beta_pos = b_pos, sigma_pos = sigma))
}

test_that("lognormal_trunc cover recovers truth from bounded cover data", {
  skip_on_cran()
  skip_if_fast()
  s <- .sim_trunc(2026L)
  expect_true(all(s$data$y >= 0 & s$data$y <= 1))
  fit <- tobs(y ~ x + spatial(~ 1 || cell_idx, graph = s$adj) +
                share(spatial(), alpha = grid(c(0.5, 1.0, 1.5))),
              data = s$data, family = cover(response = "lognormal_trunc"),
              method = "nested_laplace", control = list(progress = FALSE))
  expect_s3_class(fit, "cover_fit")
  expect_equal(fit$positive, "lognormal_trunc")
  expect_true(is.finite(fit$sigma_pos) && fit$sigma_pos > 0)
  expect_true(is.na(fit$phi_pos))
  expect_lt(abs(fit$beta_occ[1] - s$truth$beta_occ[1]), 0.3)
  expect_lt(abs(fit$beta_occ[2] - s$truth$beta_occ[2]), 0.2)
  expect_lt(abs(fit$beta_pos[1] - s$truth$beta_pos[1]), 0.35)
  expect_lt(abs(fit$beta_pos[2] - s$truth$beta_pos[2]), 0.2)
  expect_lt(abs(fit$sigma_pos / s$truth$sigma_pos - 1), 0.25)

  # The pointwise log-likelihood is per-plot and finite (natural-scale truncated
  # lognormal density at occupied plots, Bernoulli elsewhere).
  ll <- tulpaObs:::.tobs_ploglik_cover(fit, n.draws = 100L)
  expect_equal(ncol(ll), nrow(s$data))
  expect_true(all(is.finite(ll)))
})

test_that("lognormal_trunc cover recovers betas + dispersion across seeds (#140)", {
  # Multi-seed recovery for the upper-truncated lognormal cover arm on the shared-
  # field joint path. Estimands: occurrence + cover-arm slopes (point recovery) and
  # the latent dispersion.
  #
  # Coverage caveat (measured, 8-seed diagnostic dev_notes/_diag_trunc_cov.R): unlike
  # the plain lognormal / beta / ordinal arms, the truncated-lognormal slope Wald CI
  # is anti-conservative on this path. The occurrence arm carries a mild ~7% shared-
  # field slope attenuation (mean 0.65 vs 0.7; its observed-Fisher SE is honest,
  # meanSE 0.045 >= empSD 0.036, so the miss is bias not under-dispersion), and the
  # truncated positive arm's mode-Hessian under-captures the truncation curvature
  # (empSD 0.058 vs meanSE 0.042, ~1.4x under-dispersed) with a ~12% upward slope
  # bias. Pooled 95% Wald coverage therefore runs ~0.6, not the 0.85 the untruncated
  # arms reach. This mirrors the documented shared-field-slope anti-conservatism and
  # is a truncation-under-Laplace property (the observed information at the mode omits
  # the truncation's third-order term), NOT a fitter regression -- so this test asserts
  # robust point recovery and guards only against catastrophic CI breakage below.
  skip_on_cran()
  skip_if_fast()
  n_seeds <- 10L
  covered <- logical(0)
  bo2 <- bp2 <- sg <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    s <- .sim_trunc(3000L + r, N = 3000L)
    fit <- tobs(y ~ x + spatial(~ 1 || cell_idx, graph = s$adj) +
                share(spatial(), alpha = grid(c(0.5, 1.0, 1.5))),
                data = s$data, family = cover(response = "lognormal_trunc"),
                method = "nested_laplace", control = list(progress = FALSE))
    expect_s3_class(fit, "cover_fit")
    bo2[r] <- fit$beta_occ[2]; bp2[r] <- fit$beta_pos[2]; sg[r] <- fit$sigma_pos
    covered <- c(covered,
                 abs(fit$beta_occ[2] - s$truth$beta_occ[2]) <= 1.96 * fit$se_occ[2],
                 abs(fit$beta_pos[2] - s$truth$beta_pos[2]) <= 1.96 * fit$se_pos[2])
  }
  expect_lt(abs(mean(bo2) - 0.7), 0.15)                 # occurrence slope
  expect_lt(abs(mean(bp2) - 0.5), 0.12)                 # cover-arm slope
  expect_lt(abs(mean(sg) / 0.8 - 1), 0.15)              # dispersion
  # CIs recover to within a truncation-inflated band, not the nominal 0.95 (see the
  # coverage caveat above): guard against total breakage, do not assert nominal.
  expect_gte(mean(covered), 0.45)
})

test_that("negligible truncation reduces to the lognormal fit", {
  skip_on_cran()
  skip_if_fast()
  # cover far below the ceiling everywhere (deeply negative intercept + a small
  # field, so eta stays well below 0 at every plot): almost no mass is truncated,
  # so the truncated-lognormal fit must match the plain lognormal fit -- including
  # sigma_pos -- on the same data. A shared phi-grid removes grid-quantization
  # differences between the two families. (With a LARGE field some cells push eta
  # to the ceiling and the two legitimately diverge there -- that is the truncated
  # model recovering the true latent sigma, tested in the recovery block above.)
  s <- .sim_trunc(7L, N = 3000L, b_pos_int = -4.0, sigma = 0.4, field_scale = 0.25)
  fml <- y ~ x + spatial(~ 1 || cell_idx, graph = s$adj) +
    share(spatial(), alpha = grid(c(0.5, 1.0, 1.5)))
  ctl <- list(progress = FALSE,
              phi.grid = exp(seq(log(0.15), log(1.1), length.out = 9)))
  fit_ln <- tobs(fml, data = s$data, family = cover(response = "lognormal"),
                 method = "nested_laplace", control = ctl)
  fit_tr <- tobs(fml, data = s$data, family = cover(response = "lognormal_trunc"),
                 method = "nested_laplace", control = ctl)
  expect_lt(max(abs(fit_tr$beta_occ - fit_ln$beta_occ)), 0.03)
  expect_lt(max(abs(fit_tr$beta_pos - fit_ln$beta_pos)), 0.03)
  expect_lt(abs(fit_tr$sigma_pos / fit_ln$sigma_pos - 1), 0.05)
})
