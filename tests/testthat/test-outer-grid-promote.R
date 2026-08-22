# Promote the joint nested-Laplace outer-grid placement record to the tobs_fit
# top level + glance(). The engine reports where the outer grid ended up --
# outer_grid_placement ("fixed" / "auto_recentered"),
# outer_grid_recenter_attempts, outer_grid_prior_added, and
# outer_grid_recenter_declined, the reason a "fixed" placement stayed fixed --
# on the object the postprocess wrappers nest at $joint_fit / $joint. A caller
# asking whether the auto grid did anything should not have to reach in there;
# an inert recenter staying invisible across a batch is what filed this.

# --------------------------------------------------------------------------- #
# Extractor: present-only, and NOT gated on the grid having moved              #
# --------------------------------------------------------------------------- #

test_that(".tobs_promote_outer_grid surfaces a recentered placement", {
  jf <- list(outer_grid_placement = "auto_recentered",
             outer_grid_recenter_attempts = 2L,
             outer_grid_prior_added = TRUE)
  og <- tulpaObs:::.tobs_promote_outer_grid(jf)
  expect_identical(og$outer_grid_placement, "auto_recentered")
  expect_identical(og$outer_grid_recenter_attempts, 2L)
  expect_true(og$outer_grid_prior_added)
  # Present-only: a recentered fit carries no decline reason, so the field is
  # absent rather than NA.
  expect_false("outer_grid_recenter_declined" %in% names(og))
})

test_that(".tobs_promote_outer_grid surfaces a FIXED placement and its reason", {
  # The case #187 exists for: the recenter was applicable and declined, and the
  # reason is the whole diagnostic. Promoting only "auto_recentered" would hide
  # exactly the fits worth looking at, so the extractor is not gated on movement.
  og <- tulpaObs:::.tobs_promote_outer_grid(
    list(outer_grid_placement = "fixed",
         outer_grid_recenter_declined = "user_pinned_axis"))
  expect_identical(og$outer_grid_placement, "fixed")
  expect_identical(og$outer_grid_recenter_declined, "user_pinned_axis")
})

test_that(".tobs_promote_outer_grid is inert without a placement record", {
  # A non-joint fit, or one from an engine predating the record.
  expect_null(tulpaObs:::.tobs_promote_outer_grid(list(foo = 1)))
  expect_null(tulpaObs:::.tobs_promote_outer_grid(NULL))
  expect_null(tulpaObs:::.tobs_promote_outer_grid("not a list"))
})

# --------------------------------------------------------------------------- #
# glance.tobs_fit: placement columns, from either location                     #
# --------------------------------------------------------------------------- #

.ogp_fit <- function(jf = NULL) {
  structure(list(
    n_fixed = 4L, n_samples = 1000L, log_prob = -120.5, converged = TRUE,
    joint_fit = jf
  ), class = c("tobs_fit", "tulpa_fit"))
}

test_that("glance.tobs_fit reports placement + decline reason", {
  g <- glance(.ogp_fit(list(outer_grid_placement = "fixed",
                            outer_grid_recenter_declined = "auto_recenter_disabled")))
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  expect_identical(g$outer_grid_placement, "fixed")
  expect_identical(g$outer_grid_recenter_declined, "auto_recenter_disabled")
})

test_that("glance.tobs_fit fills a missing decline reason with NA, not a drop", {
  # A batch summary rbind()s one row per species; the column has to exist on
  # every row or the recentered fits ragged-drop it.
  g <- glance(.ogp_fit(list(outer_grid_placement = "auto_recentered",
                            outer_grid_recenter_attempts = 1L)))
  expect_identical(g$outer_grid_placement, "auto_recentered")
  expect_true("outer_grid_recenter_declined" %in% names(g))
  expect_true(is.na(g$outer_grid_recenter_declined))
})

test_that("glance.tobs_fit reads the promoted top-level fields too", {
  fit <- .ogp_fit(NULL)
  fit$outer_grid_placement <- "auto_recentered"
  fit$outer_grid_recenter_attempts <- 2L
  g <- glance(fit)
  expect_identical(g$outer_grid_placement, "auto_recentered")
})

test_that("glance.tobs_fit adds no placement columns to a non-joint fit", {
  g <- glance(.ogp_fit(NULL))
  expect_s3_class(g, "data.frame")
  expect_false(any(grepl("^outer_grid", names(g))))
})

# --------------------------------------------------------------------------- #
# End-to-end: a real spatial occu_cover fit promotes the placement             #
# --------------------------------------------------------------------------- #

test_that("occu_cover() spatial fit surfaces outer_grid_placement at the top level", {
  skip_if_fast()
  skip_on_cran()

  N <- 40L; J <- 5L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 1, alpha = 1, seed = 187L)
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
    occurrence = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE)
  ))

  # Promoted: reading fit$outer_grid_placement works without touching $joint_fit.
  expect_true("outer_grid_placement" %in% names(fit))
  expect_true(fit$outer_grid_placement %in% c("fixed", "auto_recentered"))
  # Same value as the nested raw object it was promoted from.
  expect_identical(fit$outer_grid_placement, fit$joint_fit$outer_grid_placement)
  # A "fixed" placement carries the reason it stayed fixed; a recentered one
  # carries the attempt count instead.
  if (identical(fit$outer_grid_placement, "fixed")) {
    expect_true("outer_grid_recenter_declined" %in% names(fit))
    expect_identical(fit$outer_grid_recenter_declined,
                     fit$joint_fit$outer_grid_recenter_declined)
  } else {
    expect_true(fit$outer_grid_recenter_attempts >= 1L)
  }

  g <- glance(fit)
  expect_identical(g$outer_grid_placement, fit$outer_grid_placement)
})

test_that("cover() spatial fit surfaces outer_grid_placement at the top level", {
  skip_if_fast()
  skip_on_cran()

  # cover() nests its joint object at $joint (occu_cover() at $joint_fit), so
  # the promotion has to reach both slots.
  set.seed(186L)
  n_s <- 25L; N <- 200L
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  w_s <- 0.6 * rnorm(n_s)
  x   <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.3 + 0.7 * x + w_s[spatial_idx]))
  y <- ifelse(occur == 1L,
              pmin(exp(rnorm(N, 0.4 - 0.5 * x + w_s[spatial_idx], 0.4)), 1 - 1e-6),
              0)
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    if (s > 1L)   adj[s, s - 1L] <- 1L
    if (s < n_s)  adj[s, s + 1L] <- 1L
  }

  fit <- tobs(
    formula = ~ x + icar(graph = adj, group_var = "region"),
    data    = data.frame(x = x, region = factor(spatial_idx)),
    family  = cover("lognormal"), y = y, method = "nested_laplace",
    control = list(verbose = FALSE)
  )

  expect_true("outer_grid_placement" %in% names(fit))
  expect_true(fit$outer_grid_placement %in% c("fixed", "auto_recentered"))
  expect_identical(fit$outer_grid_placement, fit$joint$outer_grid_placement)
  expect_identical(glance(fit)$outer_grid_placement, fit$outer_grid_placement)
})
