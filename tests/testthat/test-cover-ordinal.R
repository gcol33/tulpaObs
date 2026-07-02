# Tests for cover(response = "ordinal") -- the interval-censored Gaussian
# (ordered-probit with KNOWN Braun-Blanquet thresholds) positive arm on the
# joint nested-Laplace path.

test_that("cover(response = 'ordinal') constructs and validates breaks", {
  fam <- cover(response = "ordinal", breaks = c(0.05, 0.25, 0.5))
  expect_s3_class(fam, "tobs_family")
  expect_equal(fam$params$positive, "ordinal")
  expect_equal(fam$params$breaks, c(0.05, 0.25, 0.5))
  expect_equal(fam$status, "working")

  expect_error(cover(response = "ordinal"), "requires `breaks`")
  expect_error(cover(response = "ordinal", breaks = c(0.5, 0.2)), "ascending")
  expect_error(cover(response = "ordinal", breaks = c(0.2, 1.5)), "in \\(0, 1\\)")
  expect_error(cover(response = "lognormal", breaks = 0.3), "only used for")
})

test_that("encode assigns each Braun-Blanquet representative to its own class", {
  # The myscale rep for class 'a' is 0.03 = the (1.5-3] | (3-5] boundary; with
  # upper-closed bands it would fall one class too low. Lower-closed assignment
  # must put every representative (incl. 0.03) in its own class. (Regression.)
  brks <- c(0.002, 0.015, 0.03, 0.05, 0.25, 0.50, 0.75)
  mids <- c(0.001, 0.01, 0.02, 0.03, 0.155, 0.38, 0.63, 0.86)
  dat  <- data.frame(x = rnorm(length(mids)))
  enc  <- encode_cover_hurdle(~ x, dat, y = mids, positive = "ordinal", breaks = brks)
  expect_equal(enc$pos_data$class, seq_along(mids))
  # bounds are the log class edges, open outer classes
  expect_true(is.infinite(enc$pos_data$lower[1]) && enc$pos_data$lower[1] < 0)
  expect_true(is.infinite(enc$pos_data$upper[length(mids)]) &&
              enc$pos_data$upper[length(mids)] > 0)
  expect_equal(enc$pos_data$lower[4], log(0.03))     # class 'a' lower edge
  expect_equal(enc$pos_data$upper[4], log(0.05))     # class 'a' upper edge
})

test_that("ordinal requires method = 'nested_laplace'", {
  dat <- data.frame(x = rnorm(40), cell = sample.int(4, 40, replace = TRUE),
                    y = ifelse(runif(40) < 0.5, 0,
                               sample(c(0.02, 0.155, 0.38), 40, TRUE)))
  expect_error(
    tobs(y ~ x, data = dat, family = cover(response = "ordinal",
         breaks = c(0.05, 0.25, 0.5)), method = "laplace"),
    "requires method = 'nested_laplace'"
  )
})

# Self-contained spatial simulator: occurrence + cover share an areal field
# (alpha-coupled); the positive cover is the latent log-cover censored to its
# Braun-Blanquet class (observed = class midpoint), so the data carry ONLY the
# ordinal class.
.sim_ordinal <- function(seed, N = 6000L, gx = 8L) {
  set.seed(seed)
  nc <- gx * gx
  adj <- matrix(0L, nc, nc); id <- function(i, j) (j - 1L) * gx + i
  for (i in 1:gx) for (j in 1:gx) for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next; ii <- i + di; jj <- j + dj
    if (ii >= 1 && ii <= gx && jj >= 1 && jj <= gx) adj[id(i, j), id(ii, jj)] <- 1L
  }
  fx <- outer(1:gx, 1:gx, function(i, j) sin(i / 2) + cos(j / 2))
  field <- as.numeric(scale(as.numeric(fx)))
  brks <- c(0.002, 0.015, 0.03, 0.05, 0.25, 0.50, 0.75)
  mids <- c(0.001, 0.01, 0.02, 0.03, 0.155, 0.38, 0.63, 0.86)
  cell <- sample.int(nc, N, replace = TRUE); x <- rnorm(N)
  b_occ <- c(-0.1, 0.7); b_pos <- c(-2.2, 0.5); sigma <- 0.85
  occ <- rbinom(N, 1L, plogis(b_occ[1] + b_occ[2] * x + field[cell]))
  cov_fr <- pmin(exp(b_pos[1] + b_pos[2] * x + field[cell] + rnorm(N, 0, sigma)),
                 1 - 1e-6)
  y <- ifelse(occ == 1L, mids[findInterval(cov_fr, brks, left.open = TRUE) + 1L], 0)
  list(data = data.frame(y = y, x = x, cell_idx = cell), adj = adj, brks = brks,
       truth = list(beta_occ = b_occ, beta_pos = b_pos, sigma_pos = sigma))
}

test_that("ordinal cover recovers truth from censored class data", {
  skip_on_cran()
  skip_if_fast()
  s <- .sim_ordinal(2026L)
  fit <- tobs(y ~ x + spatial(~ 1 || cell_idx, graph = s$adj),
              data = s$data, family = cover(response = "ordinal", breaks = s$brks),
              method = "nested_laplace", control = list(progress = FALSE))
  expect_s3_class(fit, "cover_fit")
  expect_equal(fit$positive, "ordinal")
  expect_true(is.finite(fit$sigma_pos) && fit$sigma_pos > 0)
  expect_true(is.na(fit$phi_pos))
  expect_lt(abs(fit$beta_occ[1] - s$truth$beta_occ[1]), 0.3)
  expect_lt(abs(fit$beta_occ[2] - s$truth$beta_occ[2]), 0.2)
  expect_lt(abs(fit$beta_pos[1] - s$truth$beta_pos[1]), 0.3)
  expect_lt(abs(fit$beta_pos[2] - s$truth$beta_pos[2]), 0.2)
  expect_lt(abs(fit$sigma_pos / s$truth$sigma_pos - 1), 0.2)

  # The pointwise log-likelihood is per-plot and a genuine class probability
  # (a PMF: each occupied plot's contribution is <= 0 and exp() <= 1).
  ll <- tulpaObs:::.tobs_ploglik_cover(fit, n.draws = 100L)
  expect_equal(ncol(ll), nrow(s$data))
  expect_true(all(is.finite(ll)))
})
