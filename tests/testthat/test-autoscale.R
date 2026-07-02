# Regression tests for the internal design-matrix autoscaling
# (gcol33/tulpaObs#9). The MAP / Laplace path on a design matrix with a
# numeric column on a large additive scale (e.g. raw calendar year ~ 2000)
# used to converge to a numerically nonsensical (large intercept, tiny slope)
# stationary point because the natural-scale Hessian is near-singular along
# the (intercept, mean(column)) direction. The engine now centers and scales
# numeric design columns internally; coefficients are transformed back to
# natural scale before being returned, so user-facing predictions are
# unchanged but the optimizer reaches the true MAP reliably.

test_that(".autoscale_pickcols skips intercept, binaries, and zero-sd columns", {
  X <- cbind(
    "(Intercept)"  = rep(1, 50),
    binary         = rep(0:1, 25),
    constant       = rep(2.5, 50),
    numeric_normal = rnorm(50),
    numeric_year   = sample(1990:2010, 50, replace = TRUE)
  )
  cols <- tulpaObs:::.autoscale_pickcols(X)
  nms <- colnames(X)[cols]
  expect_setequal(nms, c("numeric_normal", "numeric_year"))
})

test_that(".autoscale_design centers + scales selected columns to mean ~ 0, sd ~ 1", {
  set.seed(1)
  X <- cbind("(Intercept)" = rep(1, 100),
             year = sample(1985:2005, 100, replace = TRUE),
             effort = runif(100, 0.5, 1.5))
  ad <- tulpaObs:::.autoscale_design(X)
  expect_equal(ad$scale$cols, c(2L, 3L))
  expect_equal(mean(ad$X[, "year"]),   0, tolerance = 1e-12)
  expect_equal(sd(ad$X[, "year"]),     1, tolerance = 1e-12)
  expect_equal(mean(ad$X[, "effort"]), 0, tolerance = 1e-12)
  expect_equal(sd(ad$X[, "effort"]),   1, tolerance = 1e-12)
  expect_equal(ad$X[, "(Intercept)"], rep(1, 100))
})

test_that("unscale / scale round-trip is exact", {
  set.seed(2)
  X <- cbind("(Intercept)" = rep(1, 100),
             year   = sample(1985:2005, 100, replace = TRUE),
             effort = runif(100, 0.5, 1.5))
  sc <- tulpaObs:::.scale_meta(X)
  beta_nat <- c(-1.5, 0.05, 0.8)
  beta_sc  <- tulpaObs:::.scale_beta_vec(beta_nat, sc)
  expect_equal(tulpaObs:::.unscale_beta_vec(beta_sc, sc), beta_nat,
               tolerance = 1e-10)
})

test_that("cover(response='beta') with year ~ 2000 column reaches the well-conditioned MAP", {
  # Issue-#9 repro: a raw calendar-year column on the additive scale ~ 2000
  # forms a near-singular pair with the intercept in the natural-scale
  # Hessian. Autoscale internally so the optimizer converges to the true
  # MAP. Compare to manually-scaled columns: predictions on the training
  # data must be identical (the two parameterizations are mathematically
  # equivalent), and slope * sd(year_btw) must match the user-side scaled
  # slope coefficient.
  skip_on_cran()
  skip_if_fast()
  suppressPackageStartupMessages({library(dplyr, quietly = TRUE)})

  set.seed(2026)
  n_groups <- 60L; visits_per_grp <- 6L
  rs_groups <- data.frame(
    group_id = paste0("G_", seq_len(n_groups)),
    base_year = sample(1985:2005, n_groups, replace = TRUE)
  )
  dat <- rs_groups[rep(seq_len(n_groups), each = visits_per_grp), ] %>%
    dplyr::mutate(visit_no = rep(seq_len(visits_per_grp), n_groups),
                  year = base_year + 3 * (visit_no - 1))
  psi_true <- 0.4; p_true <- 0.6; mu_int <- -1.2; mu_slope <- 0.02
  z_per_grp <- stats::rbinom(n_groups, 1, psi_true)
  dat$occur_true <- z_per_grp[match(dat$group_id, rs_groups$group_id)]
  dat$detect <- ifelse(dat$occur_true == 1,
                       stats::rbinom(nrow(dat), 1, p_true), 0)
  mu_cover <- plogis(mu_int + mu_slope * (dat$year - 2000))
  phi_beta <- 8
  cov_pos <- stats::rbeta(nrow(dat), mu_cover * phi_beta,
                          (1 - mu_cover) * phi_beta)
  dat$y <- ifelse(dat$detect == 1, cov_pos, 0)
  dat <- within_between(dat, group = "group_id", vars = "year")

  fit_raw <- tobs(~ year_btw + year_wtn, data = dat,
                  family = cover(response = "beta"), y = dat$y,
                  method = "laplace")

  dat_sc <- dat %>%
    dplyr::mutate(year_btw_sc = as.numeric(scale(year_btw)),
                  year_wtn_sc = as.numeric(scale(year_wtn)))
  fit_sc <- tobs(~ year_btw_sc + year_wtn_sc, data = dat_sc,
                 family = cover(response = "beta"), y = dat_sc$y,
                 method = "laplace")

  # Predictions on the training data must match between the two fits
  # (both parameterizations of the same MAP).
  p_raw <- predict(fit_raw, newdata = dat,    type = "expected")
  p_sc  <- predict(fit_sc,  newdata = dat_sc, type = "expected")
  expect_lt(max(abs(p_raw - p_sc)), 1e-6)

  # Slope cross-check: a$year_btw * sd(year_btw) == b$year_btw_sc.
  sd_ybtw <- sd(dat$year_btw)
  sd_ywtn <- sd(dat$year_wtn)
  expect_equal(unname(fit_raw$beta_occ[["year_btw"]] * sd_ybtw),
               unname(fit_sc$beta_occ[["year_btw_sc"]]), tolerance = 1e-6)
  expect_equal(unname(fit_raw$beta_pos[["year_btw"]] * sd_ybtw),
               unname(fit_sc$beta_pos[["year_btw_sc"]]),
               tolerance = 1e-6)
  expect_equal(unname(fit_raw$beta_pos[["year_wtn"]] * sd_ywtn),
               unname(fit_sc$beta_pos[["year_wtn_sc"]]), tolerance = 1e-6)
})

test_that("occu() with mean-2000 covariate matches manually-scaled fit (predictions)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(7)
  N <- 80L; J <- 4L
  year <- sample(1985:2005, N, replace = TRUE)  # ~ mean 2000, sd ~ 6
  x_det <- rnorm(N)
  truth_int <- -0.5; truth_slope <- 0.03
  psi <- plogis(truth_int + truth_slope * (year - 1995))
  z   <- rbinom(N, 1L, psi)
  p   <- plogis(0 + 0.5 * x_det)
  y   <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p[i])
  }
  dat   <- data.frame(year = year, x_det = x_det)
  dat_s <- dat; dat_s$year_sc <- as.numeric(scale(dat_s$year))

  fit_raw <- tobs(~ year, data = dat, family = occu(),
                  detection = ~ x_det, y = y, method = "laplace",
                  control = list(verbose = FALSE))
  fit_sc  <- tobs(~ year_sc, data = dat_s, family = occu(),
                  detection = ~ x_det, y = y, method = "laplace",
                  control = list(verbose = FALSE))

  # Probability-scale intercepts (psi(0) under each parameterization) differ
  # because the formula labels differ. What must match is in-sample fitted
  # psi from `fitted()`.
  f_raw <- fitted(fit_raw)
  f_sc  <- fitted(fit_sc)
  expect_lt(max(abs(f_raw$psi - f_sc$psi)), 1e-4)
  expect_lt(max(abs(f_raw$p   - f_sc$p)),   1e-4)
})
