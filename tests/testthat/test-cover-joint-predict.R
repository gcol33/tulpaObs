# predict() for the cover() hurdle on the nested-Laplace shared-field path
# and the unified joint-fit substrate it shares with occu_cover().

.cjp_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < n)  adj[s, s + 1L] <- 1L
  }
  adj
}

# Joint lognormal cover hurdle with a shared field over `n_s` regions.
.cjp_build_fit <- function(N = 220L, n_s = 24L, prior = c("bym2", "icar"),
                           seed = 41L) {
  prior <- match.arg(prior)
  set.seed(seed)
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  adj <- .cjp_chain_adj(n_s)
  phi <- rnorm(n_s); theta <- rnorm(n_s)
  w_s <- 0.6 * (sqrt(0.7) * phi + sqrt(0.3) * theta)
  x   <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.3 + 0.7 * x + w_s[spatial_idx]))
  eta_pos <- 0.4 - 0.5 * x + w_s[spatial_idx]
  y <- ifelse(occur == 1L, pmin(exp(rnorm(N, eta_pos, 0.4)), 1 - 1e-6), 0)
  dat <- data.frame(x = x, region = factor(spatial_idx))

  term <- if (prior == "bym2") quote(bym2(graph = adj, group_var = "region"))
          else quote(icar(graph = adj, group_var = "region"))
  fm <- stats::reformulate(c("x", deparse(term)))

  fit <- tobs(
    formula = fm, data = dat, family = cover("lognormal"), y = y,
    method = "nested_laplace",
    control = list(sigma.grid = c(0.4, 0.8), rho.grid = c(0.5, 0.9))
  )
  list(fit = fit, n_s = n_s, N = N)
}

test_that("the unified substrate accesses both family joint objects (#24)", {
  skip_on_cran()
  skip_if_fast()
  f <- .cjp_build_fit()
  jf <- tulpaObs:::.tobs_joint_fit(f$fit)
  expect_s3_class(jf, "tulpa_nested_laplace_joint")
  # cover() stores at $joint; the accessor finds it without the slot name.
  expect_identical(jf, f$fit$joint)

  bundle <- tulpaObs:::.tobs_joint_draws(f$fit, n = 50L)
  expect_equal(bundle$n, 50L)
  expect_equal(bundle$n_cells, f$n_s)
  expect_equal(nrow(bundle$b$occ), 50L)
  expect_null(bundle$b$det)                       # cover has no detection arm
  expect_equal(ncol(bundle$blocks[[1L]]$z), f$n_s)
  expect_length(bundle$blocks, 1L)                # single shared field
})

test_that("cover() joint predict: single-quantity tobs_prediction per cell (#23)", {
  skip_on_cran()
  skip_if_fast()
  f <- .cjp_build_fit()
  nd <- data.frame(x = 0, cell = seq_len(f$n_s))
  for (ty in c("occurrence", "cover_cond", "cover_exp")) {
    pr <- predict(f$fit, newdata = nd, type = ty, nsim = 400L)
    expect_s3_class(pr, "tobs_prediction")
    expect_true(all(c("cell", "mean", "sd", "lwr", "upr") %in% names(pr)))
    expect_equal(nrow(pr), f$n_s)
    expect_true(all(pr$lwr <= pr$upr))
    dm <- attr(pr, "draws")[[ty]]
    expect_equal(dim(dm), c(f$n_s, 400L))
  }
  pr <- predict(f$fit, newdata = nd, type = "occurrence", nsim = 200L)
  expect_true(all(pr$mean >= 0 & pr$mean <= 1))
  pr <- predict(f$fit, newdata = nd, type = "cover_exp", nsim = 200L)
  expect_true(all(pr$mean >= 0))
})

test_that("cover() joint predict: change decomposition identity holds (#23)", {
  skip_on_cran()
  skip_if_fast()
  f <- .cjp_build_fit()
  nd <- data.frame(x = 0, cell = seq_len(f$n_s))
  pr <- predict(f$fit, newdata = nd, type = "change",
                times = c(-1, 1), time_col = "x", nsim = 500L)
  expect_s3_class(pr, "tobs_prediction")
  expect_equal(nrow(pr), f$n_s)
  point_cols <- c("p_T1", "p_T2", "delta_p", "cover_cond_T1", "cover_cond_T2",
                  "delta_cover_cond", "cover_exp_T1", "cover_exp_T2",
                  "delta_cover_exp", "delta_cover_from_occ", "delta_cover_from_ab")
  expect_true(all(c("cell", point_cols) %in% names(pr)))
  dr <- attr(pr, "draws")
  expect_equal(dr$delta_cover_from_occ + dr$delta_cover_from_ab,
               dr$delta_cover_exp, tolerance = 1e-6)
  expect_equal(dr$delta_cover_exp, dr$cover_exp_T2 - dr$cover_exp_T1,
               tolerance = 1e-6)
  # x drives occupancy -> delta_p varies across cells.
  expect_gt(stats::sd(pr$delta_p), 0)
})

test_that("cover() joint predict reconstructs the ICAR field too (#24)", {
  skip_on_cran()
  skip_if_fast()
  f <- .cjp_build_fit(prior = "icar", seed = 53L)
  nd <- data.frame(x = 0, cell = seq_len(f$n_s))
  pr <- predict(f$fit, newdata = nd, type = "cover_exp", nsim = 300L)
  expect_s3_class(pr, "tobs_prediction")
  expect_equal(nrow(pr), f$n_s)
  expect_true(all(pr$mean >= 0))
})

test_that("legacy fixed-effects type names map onto the joint vocabulary (#23)", {
  skip_on_cran()
  skip_if_fast()
  f <- .cjp_build_fit()
  nd <- data.frame(x = 0, cell = seq_len(f$n_s))
  a <- predict(f$fit, newdata = nd, type = "expected", nsim = 300L)
  b <- predict(f$fit, newdata = nd, type = "cover_exp", nsim = 300L)
  expect_s3_class(a, "tobs_prediction")
  expect_equal(attr(a, "quantity"), "cover_exp")
  expect_equal(attr(b, "quantity"), "cover_exp")
})

test_that("separate-Laplace cover() predict is unchanged (numeric vector)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(7)
  N <- 150L
  x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.2 + 0.6 * x))
  y <- ifelse(occur == 1L, pmin(exp(rnorm(N, 0.3 - 0.4 * x, 0.4)), 1 - 1e-6), 0)
  dat <- data.frame(x = x)
  fit <- tobs(formula = ~ x, data = dat, family = cover("lognormal"),
              y = y, method = "laplace")
  expect_null(tulpaObs:::.tobs_joint_fit(fit))
  pred <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "expected")
  expect_true(is.numeric(pred))
  expect_length(pred, 3L)
  # default type still resolves to "expected".
  expect_equal(predict(fit, newdata = data.frame(x = 0)),
               predict(fit, newdata = data.frame(x = 0), type = "expected"))
})
