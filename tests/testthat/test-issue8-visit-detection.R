# Regression + recovery test for gcol33/tulpaObs#8: the tobs_data() output
# (det.covs as a named list of [N, J] matrices) must compose with
# tobs(visits = ...), and the visit-level detection covariate it carries
# must be estimated correctly. Before the fix the call errored with
# "object 'effort' not found"; the visit-level-detection path - the whole
# point of occu() - had no clean public route and was unexercised by tests.

test_that("tobs_data() det.covs list composes with tobs(visits = ...)", {
  set.seed(1)
  N <- 40L; J <- 4L
  df <- data.frame(
    site_id = rep(seq_len(N), each = J),
    visit   = rep(seq_len(J), times = N),
    occur   = rbinom(N * J, 1, 0.3),
    effort  = rnorm(N * J)
  )
  od <- tobs_data(df, y = "occur", site = "site_id", visit = "visit",
                  det.covs = "effort")
  # det.covs arrives as a named list of [N, J] matrices (the shape #8 is about).
  expect_true(is.list(od$det.covs))
  expect_equal(dim(od$det.covs$effort), c(N, J))

  fit <- tobs(~ 1, data = data.frame(site_id = unique(df$site_id)),
              y = od$y, detection = ~ effort, visits = od$det.covs,
              family = occu(), method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  # The visit-level covariate is carried into the detection design.
  expect_true("p_visit_effort" %in% names(fit$means))
})

test_that("visit-level detection slope is recovered (Laplace, 20 seeds)", {
  skip_on_cran()

  n_seeds <- 20L
  N <- 300L; J <- 6L
  p0 <- 0.2; p1 <- 1.2; psi <- plogis(0.4)

  est <- se <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    set.seed(5000L + s)
    zocc <- rbinom(N, 1, psi)
    eff  <- matrix(rnorm(N * J), N, J)
    y <- matrix(0L, N, J)
    for (i in seq_len(N)) {
      y[i, ] <- rbinom(J, 1, zocc[i] * plogis(p0 + p1 * eff[i, ]))
    }
    df <- data.frame(
      site_id = rep(seq_len(N), each = J),
      visit   = rep(seq_len(J), times = N),
      occur   = as.vector(t(y)),
      effort  = as.vector(t(eff))
    )
    od <- tobs_data(df, y = "occur", site = "site_id", visit = "visit",
                    det.covs = "effort")
    fit <- tryCatch(
      tobs(~ 1, data = data.frame(site_id = unique(df$site_id)),
           y = od$y, detection = ~ effort, visits = od$det.covs,
           family = occu(), method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    est[s] <- fit$means["p_visit_effort"]
    se[s]  <- fit$sds["p_visit_effort"]
  }

  keep <- is.finite(est) & is.finite(se) & se > 0
  expect_gte(sum(keep), floor(0.8 * n_seeds))

  # Point estimate unbiased for the visit-varying detection slope.
  expect_lt(abs(mean(est[keep]) - p1), 0.15)
  # 95% Wald CI covers truth at near-nominal rate (>= 0.80 floor for seed
  # noise at n_seeds = 20, per the CLAUDE.md recovery rubric).
  covered <- abs(est[keep] - p1) < 1.96 * se[keep]
  expect_gte(mean(covered), 0.80)
})
