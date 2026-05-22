# Tests for the tobs(visit_data = ...) composition fix (issue #8).
# Verifies that tobs_data() output (list of [N, J] matrices) composes with
# tobs() without requiring the user to hand-build a long visit-level frame.

simulate_visit_effort_data <- function(seed = 42) {
  set.seed(seed)
  n_sites <- 60
  max_visits <- 4
  elev <- rnorm(n_sites)
  psi <- plogis(-0.3 + 0.8 * elev)
  z <- rbinom(n_sites, 1, psi)

  effort_mat <- matrix(runif(n_sites * max_visits, 0.2, 1.5),
                       n_sites, max_visits)
  beta_eff <- 1.5
  p_mat <- plogis(qlogis(0.4) + beta_eff * (effort_mat - mean(effort_mat)))
  y <- matrix(0L, n_sites, max_visits)
  for (i in seq_len(n_sites)) {
    if (z[i] == 1) {
      for (j in seq_len(max_visits)) {
        y[i, j] <- rbinom(1, 1, p_mat[i, j])
      }
    }
  }

  long <- data.frame(
    site_id = rep(paste0("s", seq_len(n_sites)), each = max_visits),
    year    = rep(seq_len(max_visits), times = n_sites),
    occur   = as.vector(t(y)),
    effort  = as.vector(t(effort_mat))
  )
  list(long = long, y = y, n_sites = n_sites, max_visits = max_visits,
       elev = elev, effort_mat = effort_mat)
}

test_that("tobs_data() output composes with tobs(visit_data = ...)", {
  sim <- simulate_visit_effort_data()
  od <- tobs_data(sim$long, y = "occur", site = "site_id", visit = "year",
                  det.covs = "effort")
  expect_type(od$det.covs, "list")
  expect_true(is.matrix(od$det.covs$effort))
  expect_equal(dim(od$det.covs$effort), c(sim$n_sites, sim$max_visits))

  site_df <- data.frame(site_id = paste0("s", seq_len(sim$n_sites)))

  fit <- tobs(
    formula    = ~ 1,
    data       = site_df,
    y          = od$y,
    detection  = ~ effort,
    visit_data = od$det.covs,
    family     = occu(),
    engine     = "laplace",
    control    = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(any(grepl("p_visit_effort", colnames(fit$draws))))
})

test_that("visit_data accepts a long data frame without a formula attribute", {
  sim <- simulate_visit_effort_data()
  effort_long <- data.frame(effort = as.vector(t(sim$effort_mat)))
  site_df <- data.frame(dummy = rep(0, sim$n_sites))

  fit <- tobs(
    formula    = ~ 1,
    data       = site_df,
    y          = sim$y,
    detection  = ~ effort,
    visit_data = effort_long,
    family     = occu(),
    engine     = "laplace",
    control    = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(any(grepl("p_visit_effort", colnames(fit$draws))))
})

test_that("attr(visit_data, 'formula') keeps dual site/visit detection split", {
  sim <- simulate_visit_effort_data()
  effort_long <- data.frame(effort = as.vector(t(sim$effort_mat)))
  attr(effort_long, "formula") <- ~ effort
  site_df <- data.frame(observer = sample(c("A", "B"), sim$n_sites,
                                          replace = TRUE))

  fit <- tobs(
    formula    = ~ 1,
    data       = site_df,
    y          = sim$y,
    detection  = ~ observer,
    visit_data = effort_long,
    family     = occu(),
    engine     = "laplace",
    control    = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(any(grepl("p_observer", colnames(fit$draws))))
  expect_true(any(grepl("p_visit_effort", colnames(fit$draws))))
})

test_that("visit_data shape mismatch is rejected with a clear message", {
  sim <- simulate_visit_effort_data()
  bad_list <- list(effort = sim$effort_mat[, 1:2, drop = FALSE])
  site_df <- data.frame(dummy = rep(0, sim$n_sites))
  expect_error(
    tobs(formula = ~ 1, data = site_df, y = sim$y,
         detection = ~ effort, visit_data = bad_list,
         family = occu(), engine = "laplace",
         control = list(verbose = FALSE)),
    "wrong shape"
  )

  bad_df <- data.frame(effort = runif(sim$n_sites * sim$max_visits + 3))
  expect_error(
    tobs(formula = ~ 1, data = site_df, y = sim$y,
         detection = ~ effort, visit_data = bad_df,
         family = occu(), engine = "laplace",
         control = list(verbose = FALSE)),
    "must have .* rows"
  )

  unnamed_list <- list(sim$effort_mat)
  expect_error(
    tobs(formula = ~ 1, data = site_df, y = sim$y,
         detection = ~ effort, visit_data = unnamed_list,
         family = occu(), engine = "laplace",
         control = list(verbose = FALSE)),
    "named list"
  )
})
