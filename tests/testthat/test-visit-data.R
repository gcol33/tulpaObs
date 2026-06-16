# Tests for the tobs(visits = ...) composition fix (issue #8).
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

test_that("tobs_data() output composes with tobs(visits = ...)", {
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
    visits = od$det.covs,
    family     = occu(),
    method     = "laplace",
    control    = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(any(grepl("p_visit_effort", colnames(fit$draws))))
})

test_that("visits accepts a long data frame without a formula attribute", {
  sim <- simulate_visit_effort_data()
  effort_long <- data.frame(effort = as.vector(t(sim$effort_mat)))
  site_df <- data.frame(dummy = rep(0, sim$n_sites))

  fit <- tobs(
    formula    = ~ 1,
    data       = site_df,
    y          = sim$y,
    detection  = ~ effort,
    visits = effort_long,
    family     = occu(),
    method     = "laplace",
    control    = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(any(grepl("p_visit_effort", colnames(fit$draws))))
})

test_that("attr(visits, 'formula') keeps dual site/visit detection split", {
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
    visits = effort_long,
    family     = occu(),
    method     = "laplace",
    control    = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(any(grepl("p_observer", colnames(fit$draws))))
  expect_true(any(grepl("p_visit_effort", colnames(fit$draws))))
})

# --- Categorical (factor / character) visit-level detection covariates --------
# tobs_data() preserves a factor / character det.cov as categorical end-to-end:
# the detection / positive visit design expands it to k - 1 dummies (first level
# the reference), not a single numeric coefficient.

simulate_visit_factor_data <- function(seed = 7) {
  set.seed(seed)
  n_sites <- 300L
  max_visits <- 5L
  levs <- c("low", "mid", "high")
  p_off <- c(low = 0.0, mid = 1.2, high = -1.0)

  elev <- rnorm(n_sites)
  psi  <- plogis(-0.2 + 0.9 * elev)
  z    <- rbinom(n_sites, 1, psi)

  hab <- matrix(sample(levs, n_sites * max_visits, replace = TRUE),
                n_sites, max_visits)
  p_mat <- matrix(plogis(stats::qlogis(0.45) + p_off[hab]), n_sites, max_visits)
  y <- matrix(0L, n_sites, max_visits)
  for (i in seq_len(n_sites)) {
    if (z[i] == 1L) {
      for (j in seq_len(max_visits)) y[i, j] <- rbinom(1, 1, p_mat[i, j])
    }
  }

  long <- data.frame(
    site_id = rep(seq_len(n_sites), each = max_visits),
    visit   = rep(seq_len(max_visits), times = n_sites),
    occur   = as.vector(t(y)),
    habitat = factor(as.vector(t(hab)), levels = levs),
    habitat_chr = as.vector(t(hab)),
    eff     = as.vector(t(matrix(rnorm(n_sites * max_visits),
                                 n_sites, max_visits)))
  )
  list(long = long, y = y, n_sites = n_sites, max_visits = max_visits,
       elev = elev, levs = levs, p_off = p_off)
}

test_that("tobs_data() preserves a factor det.cov as a tagged character matrix", {
  sim <- simulate_visit_factor_data()
  od <- tobs_data(sim$long, y = "occur", site = "site_id", visit = "visit",
                  det.covs = c("eff", "habitat"))

  # Numeric column keeps the existing double-matrix path.
  expect_true(is.double(od$det.covs$eff))
  expect_null(attr(od$det.covs$eff, "tobs_factor"))

  # Factor column becomes a character matrix carrying its level set.
  hm <- od$det.covs$habitat
  expect_true(is.character(hm))
  expect_true(isTRUE(attr(hm, "tobs_factor")))
  expect_identical(attr(hm, "tobs_levels"), sim$levs)
  expect_equal(dim(hm), c(sim$n_sites, sim$max_visits))
})

test_that("character det.cov derives its level set from sorted unique values", {
  sim <- simulate_visit_factor_data()
  od <- tobs_data(sim$long, y = "occur", site = "site_id", visit = "visit",
                  det.covs = "habitat_chr")
  hm <- od$det.covs$habitat_chr
  expect_true(isTRUE(attr(hm, "tobs_factor")))
  expect_identical(attr(hm, "tobs_levels"), sort(unique(sim$long$habitat_chr)))
})

test_that("factor det.cov fits as k - 1 categorical detection coefficients (occu)", {
  sim <- simulate_visit_factor_data()
  od <- tobs_data(sim$long, y = "occur", site = "site_id", visit = "visit",
                  det.covs = "habitat")
  site_df <- data.frame(site_id = seq_len(sim$n_sites), elev = sim$elev)

  fit <- tobs(
    formula   = ~ elev,
    data      = site_df,
    y         = od$y,
    detection = ~ habitat,
    visits    = od$det.covs,
    family    = occu(),
    method    = "laplace",
    control   = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")

  # A 3-level factor yields exactly two categorical detection terms (mid, high
  # against the reference level low), not a single numeric "habitat" slope.
  det_terms <- grep("^p_visit_habitat", names(fit$means), value = TRUE)
  expect_setequal(det_terms, c("p_visit_habitatmid", "p_visit_habitathigh"))
  expect_false("p_visit_habitat" %in% names(fit$means))

  # Calibrated recovery: estimated contrasts match the simulated per-level
  # offsets in sign and ordering (low = 0 reference; mid > 0 > high).
  est_mid  <- unname(fit$means["p_visit_habitatmid"])
  est_high <- unname(fit$means["p_visit_habitathigh"])
  expect_gt(est_mid, 0)
  expect_lt(est_high, 0)
  expect_gt(est_mid, est_high)
  expect_lt(abs(est_mid  - sim$p_off[["mid"]]),  0.5)
  expect_lt(abs(est_high - sim$p_off[["high"]]), 0.5)
})

test_that("factor det.cov fits categorically on both occu_cover arms (beta)", {
  skip_on_cran()
  set.seed(11)
  N <- 300L; J <- 5L
  levs <- c("A", "B", "C")
  p_off   <- c(A = 0.0, B = 1.1, C = -0.9)
  pos_off <- c(A = 0.0, B = 0.8, C = -0.6)

  elev <- rnorm(N)
  z    <- rbinom(N, 1, plogis(-0.1 + 0.8 * elev))
  hab  <- matrix(sample(levs, N * J, replace = TRUE), N, J)
  p_mat    <- matrix(plogis(stats::qlogis(0.5) + p_off[hab]), N, J)
  mu_logit <- matrix(stats::qlogis(0.25) + pos_off[hab], N, J)
  phi <- 12

  y <- matrix(0L, N, J); y_pos <- matrix(0, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) for (j in seq_len(J)) {
    y[i, j] <- rbinom(1, 1, p_mat[i, j])
    if (y[i, j] == 1L) {
      mu <- plogis(mu_logit[i, j])
      y_pos[i, j] <- rbeta(1, mu * phi, (1 - mu) * phi)
    }
  }
  y_pos[y_pos <= 0] <- 1e-4; y_pos[y_pos >= 1] <- 1 - 1e-4
  y_pos[y != 1L] <- 0

  long <- data.frame(
    site_id = rep(seq_len(N), each = J),
    visit   = rep(seq_len(J), times = N),
    occur   = as.vector(t(y)),
    habitat = factor(as.vector(t(hab)), levels = levs)
  )
  od <- tobs_data(long, y = "occur", site = "site_id", visit = "visit",
                  det.covs = "habitat")
  cell <- data.frame(site_id = seq_len(N), elev = elev)

  fit <- tobs(
    formula   = ~ elev,
    data      = cell,
    family    = occu_cover("beta"),
    detection = ~ habitat,
    positive  = ~ habitat,
    y         = od$y, y_pos = y_pos, visits = od$det.covs,
    method    = "laplace", control = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")

  det_terms <- grep("^p_habitat", names(fit$means), value = TRUE)
  pos_terms <- grep("^pos_habitat", names(fit$means), value = TRUE)
  expect_setequal(det_terms, c("p_habitatB", "p_habitatC"))
  expect_setequal(pos_terms, c("pos_habitatB", "pos_habitatC"))

  expect_gt(unname(fit$means["p_habitatB"]), 0)
  expect_lt(unname(fit$means["p_habitatC"]), 0)
  expect_gt(unname(fit$means["pos_habitatB"]), 0)
  expect_lt(unname(fit$means["pos_habitatC"]), 0)
  expect_lt(abs(fit$means[["p_habitatB"]] - p_off[["B"]]),     0.5)
  expect_lt(abs(fit$means[["pos_habitatB"]] - pos_off[["B"]]), 0.5)
})

test_that("numeric det.cov output is unchanged by categorical support", {
  sim <- simulate_visit_factor_data()

  # tobs_data() numeric det.cov matrix matches a hand-built double matrix
  # (the pre-change path: 2D index fill, no categorical tagging).
  od <- tobs_data(sim$long, y = "occur", site = "site_id", visit = "visit",
                  det.covs = "eff")
  ref <- matrix(NA_real_, sim$n_sites, sim$max_visits)
  ref[cbind(sim$long$site_id, sim$long$visit)] <- as.numeric(sim$long$eff)
  expect_identical(od$det.covs$eff, ref)

  # A numeric-only fit recovers the same effort coefficient as the existing
  # visit-data case, with no stray categorical columns.
  site_df <- data.frame(site_id = seq_len(sim$n_sites), elev = sim$elev)
  fit <- tobs(
    formula   = ~ elev,
    data      = site_df,
    y         = od$y,
    detection = ~ eff,
    visits    = od$det.covs,
    family    = occu(),
    method    = "laplace",
    control   = list(verbose = FALSE)
  )
  expect_true("p_visit_eff" %in% names(fit$means))
  expect_equal(sum(grepl("^p_visit_", names(fit$means))), 1L)
})

test_that("visits shape mismatch is rejected with a clear message", {
  sim <- simulate_visit_effort_data()
  bad_list <- list(effort = sim$effort_mat[, 1:2, drop = FALSE])
  site_df <- data.frame(dummy = rep(0, sim$n_sites))
  expect_error(
    tobs(formula = ~ 1, data = site_df, y = sim$y,
         detection = ~ effort, visits = bad_list,
         family = occu(), method = "laplace",
         control = list(verbose = FALSE)),
    "wrong shape"
  )

  bad_df <- data.frame(effort = runif(sim$n_sites * sim$max_visits + 3))
  expect_error(
    tobs(formula = ~ 1, data = site_df, y = sim$y,
         detection = ~ effort, visits = bad_df,
         family = occu(), method = "laplace",
         control = list(verbose = FALSE)),
    "must have .* rows"
  )

  unnamed_list <- list(sim$effort_mat)
  expect_error(
    tobs(formula = ~ 1, data = site_df, y = sim$y,
         detection = ~ effort, visits = unnamed_list,
         family = occu(), method = "laplace",
         control = list(verbose = FALSE)),
    "named list"
  )
})
