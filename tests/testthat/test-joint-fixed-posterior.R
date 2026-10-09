# =============================================================================
# test-joint-fixed-posterior.R - the fixed-effect accessors of a joint
# nested-Laplace fit read the engine's grid mixture (exact moments,
# mixture-CDF bounds), not quantiles of the Monte Carlo draws (#385).
# =============================================================================


test_that("joint occu_cover accessors read the engine's mixture posterior", {
  skip_if_fast()
  skip_on_cran()
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal",
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y.pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE)
  ))
  p  <- fit$n_fixed
  nm <- names(fit$means)[seq_len(p)]

  jf <- fit$joint_fit
  jf$fixed_names <- nm
  ref <- confint(jf)

  ci <- confint(fit)
  expect_identical(attr(ci, "interval_source"), "mixture_cdf")
  expect_equal(attr(ci, "retained_mass"), 1)
  expect_identical(unclass(ci)[, 1:2], unclass(ref)[, 1:2])

  expect_equal(coef(fit), fit$means[seq_len(p)], tolerance = 1e-10)
  expect_equal(sqrt(diag(vcov(fit))), fit$sds[seq_len(p)], tolerance = 1e-10)

  s <- summary(fit)
  expect_identical(attr(s, "interval_source"), "mixture_cdf")
  expect_equal(unname(as.matrix(s[, 3:4])), unname(ci[, 1:2]))
  td <- tidy(fit)
  expect_equal(td$conf.low, unname(ci[, 1L]))
  expect_equal(td$estimate, unname(coef(fit)))

  # A hyperparameter asked for beside a fixed effect is still read off the
  # draws; the fixed effect keeps its mixture bound.
  mixed <- confint(fit, parm = c(nm[p], "sigma"))
  expect_identical(rownames(mixed), c(nm[p], "sigma"))
  expect_equal(mixed[nm[p], ], ci[nm[p], ])
  expect_true(all(is.finite(mixed["sigma", ])))
  expect_equal(confint(fit, parm = 2:3)[, 1:2], ci[2:3, 1:2])

  # An engine fit whose moments are not the reported ones is not read.
  off <- fit; off$means[1L] <- off$means[1L] + 1
  expect_null(.tobs_fixed_posterior(off))
  # Without the engine fit the mixture the draws were sampled from is read;
  # without either, the draw quantiles are.
  nojf <- fit; nojf$joint_fit <- NULL
  expect_identical(attr(confint(nojf), "interval_source"), "mixture_cdf")
  expect_equal(confint(nojf)[, 1:2], ci[, 1:2], tolerance = 1e-8)
  bare <- nojf; attr(bare$draws, "grid_mixture") <- NULL
  expect_null(attr(confint(bare), "interval_source"))
})


test_that("an autoscaled areal count fit reads the natural-scale mixture", {
  skip_if_fast()
  skip_on_cran()
  set.seed(7L)
  A <- rook_adj(8L); N <- nrow(A)
  coord <- expand.grid(r = seq_len(8L), c = seq_len(8L))
  f <- 0.7 * scale(sin(coord$r / 8 * pi) + cos(coord$c / 8 * pi))[, 1]
  data <- data.frame(x = stats::rnorm(N, mean = 3, sd = 2))
  y <- stats::rpois(N, exp(0.5 + 0.3 * data$x + f - mean(f)))
  fit <- tobs(~ x + icar(graph = A), data = data, y = y,
              family = count(), method = "nested_laplace",
              control = list(progress = FALSE, verbose = FALSE))
  p <- fit$n_fixed %||% 2L
  ci <- confint(fit)
  expect_identical(attr(ci, "interval_source"), "mixture_cdf")
  expect_equal(unname(coef(fit)[seq_len(p)]), unname(fit$means[seq_len(p)]),
               tolerance = 1e-10)
  # The natural-scale SDs are the diagonal of the transformed covariance, so
  # they agree with vcov() and not only approximately.
  expect_equal(unname(sqrt(diag(vcov(fit)))), unname(fit$sds[seq_len(p)]),
               tolerance = 1e-10)
})
