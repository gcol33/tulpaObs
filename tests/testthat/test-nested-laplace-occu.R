## Smoke tests for the single-season occupancy nested-Laplace path. These
## exercise tobs(method = "nested_laplace") end-to-end with the multi-block
## prior built from spatial + temporal + re. Recovery thresholds for the
## underlying multi-block engine are still relaxed (per
## `dev_notes/plan_multi_block.md` Phase D), so these tests only assert that
## the pipeline runs and returns a sensibly-shaped `tobs_fit`. Recovery is
## tracked under a follow-up.

simulate_panel_occu <- function(n_sites = 20, n_visits = 4, n_times = 5,
                                psi_intercept = 0, p = 0.4, seed = 1) {
  set.seed(seed)
  # Chain graph for areal spatial; each site gets one timestamp.
  adj <- matrix(0, n_sites, n_sites)
  for (i in seq_len(n_sites - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  time_idx <- sample(seq_len(n_times), n_sites, replace = TRUE)
  obs_grp  <- sample(seq_len(max(2L, n_sites %/% 4L)), n_sites, replace = TRUE)
  z <- rbinom(n_sites, 1, plogis(psi_intercept))
  y <- matrix(0L, n_sites, n_visits)
  for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(n_visits, 1, p)
  list(
    adj = adj,
    data = data.frame(year = time_idx, obs = obs_grp,
                      x = rnorm(n_sites)),
    y = y
  )
}


test_that("tobs(engine='nested_laplace') runs with spatial only", {
  d <- simulate_panel_occu()
  adj <- d$adj
  expect_silent(
    fit <- tobs(~ x + bym2(graph = adj), data = d$data, family = occu(),
                detection = ~ 1, y = d$y,
                method = "nested_laplace",
                control = list(max.iter = 5L, verbose = FALSE))
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(!is.null(fit$nested_laplace$multi_prior))
  expect_identical(fit$nested_laplace$multi_prior[[1]]$type, "bym2")
})


test_that("tobs(engine='nested_laplace') runs with spatial + temporal + re", {
  d <- simulate_panel_occu()
  adj <- d$adj

  # The 3-block joint grid (6 BYM2 x 6 AR1 x 3 IID = 108 cells) trips
  # the engine's >50-cell soft warning. CCD integration around the joint
  # pilot mode is the follow-up; suppress here to keep test output clean.
  fit <- suppressWarnings(
    tobs(~ x + bym2(graph = adj) + temporal(year, type = "ar1") + re(obs),
         data = d$data, family = occu(),
         detection = ~ 1, y = d$y,
         method = "nested_laplace",
         control = list(max.iter = 5L, verbose = FALSE))
  )

  expect_s3_class(fit, "tobs_fit")
  mp <- fit$nested_laplace$multi_prior
  expect_true(is.list(mp) && length(mp) == 3L)
  expect_identical(vapply(mp, function(b) b$type, character(1)),
                   c("bym2", "ar1", "iid"))
})


test_that("nested_laplace requires at least one latent block", {
  d <- simulate_panel_occu()
  expect_error(
    tobs(~ x, data = d$data, family = occu(),
         detection = ~ 1, y = d$y,
         method = "nested_laplace",
         control = list(max.iter = 2L, verbose = FALSE)),
    "at least one latent block"
  )
})


test_that(".map_engine routes nested_laplace for the multi-block families, errors for jsdm", {
  for (fam in c("occu", "int_occu", "ms_occu", "dyn_occu")) {
    expect_identical(.map_engine("nested_laplace", family = fam),
                     "nested_laplace")
  }
  # jsdm has no nested-Laplace driver; the registry rejects it before dispatch,
  # so reaching .map_engine with it is an internal mis-wire (not a silent
  # downgrade to single-Laplace).
  expect_error(
    .map_engine("nested_laplace", family = "jsdm"),
    "Internal error"
  )
  expect_identical(.map_engine("laplace", family = "occu"), "laplace")
  expect_identical(.map_engine("nuts", family = "occu"), "nuts")
})


test_that(".tobs_to_multi_block_prior shapes the block list correctly", {
  d <- simulate_panel_occu()
  model <- list(
    model_type  = "single",
    n_sites     = nrow(d$y),
    # The state-block geometry is read from X_processes[[1]] (one occ row per
    # site for single-season); a 1-column placeholder suffices here.
    X_processes = list(matrix(0, nrow(d$y), 1L)),
    data        = d$data
  )
  class(model) <- "tobs_model"

  sp <- .tobs_term_bym2(d$adj)
  tm <- .tobs_term_temporal(d$data$year, type = "ar1")
  re <- .tobs_term_re(d$data$obs)

  out <- .tobs_to_multi_block_prior(spatial = sp, temporal = tm,
                                    re = list(re), model = model)
  expect_true(is.list(out) && length(out) == 3L)
  types <- vapply(out, function(b) b$type, character(1))
  expect_identical(types, c("bym2", "ar1", "iid"))
  # spatial_idx is 1:n_sites for single-season occupancy
  expect_identical(out[[1]]$spatial_idx, seq_len(nrow(d$y)))
  # temporal_idx resolves from data
  expect_identical(length(out[[2]]$temporal_idx), nrow(d$y))
})
