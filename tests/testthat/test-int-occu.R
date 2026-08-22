# each source of an integrated response surveys its own subset of the sites in
# `data`, so the site a source row measures is a join on a site key --
# `site.map`, or the response's row names against rownames(data) -- and never
# the row's position. simulate_int_occu() gives source 2 onward a non-contiguous
# site set, so reading a source in row order attaches the tail of every such
# source to the wrong latent states, silently.

test_that("each source joins to the sites it names (#241)", {
  sim <- simulate_int_occu(N_total = 40, n_data = 2, J = c(3, 2),
                           n_shared = 20, seed = 1)

  # The layout the join has to survive: source 2 covers the shared block and
  # then skips a run of sites, so its rows are not the leading 1..n_rows.
  expect_false(identical(as.integer(sim$site_maps[[2]]),
                         seq_along(sim$site_maps[[2]])))

  model <- tulpaObs:::.tobs_build_integrated(~ x, ~ x, sim$data, sim$y)

  expect_identical(model$site_maps,
                   lapply(sim$site_maps, function(m) as.integer(m) - 1L))

  # Each source's detection design carries the covariate of the site it named,
  # and is empty at every site the source did not survey.
  for (s in seq_len(2L)) {
    rows <- as.integer(sim$site_maps[[s]])
    X <- model$X_processes[[1 + s]]
    expect_equal(unname(X[rows, "x"]), sim$data$x[rows])
    expect_true(all(X[-rows, ] == 0))
  }
})

test_that("row names win over row order at equal row count (#241)", {
  d <- data.frame(x = as.numeric(1:5))
  rownames(d) <- paste0("s", 1:5)
  ys <- matrix(0L, 5, 2)
  rownames(ys) <- paste0("s", 5:1)

  model <- tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(A = ys))
  expect_identical(model$site_maps[[1]], 5:1 - 1L)
})

test_that("a site name `data` does not carry is an error naming it (#241)", {
  d <- data.frame(x = as.numeric(1:5))
  rownames(d) <- paste0("s", 1:5)
  ys <- matrix(0L, 3, 2)
  rownames(ys) <- c("s1", "s9", "s7")

  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(A = ys)),
    "'s9', 's7'", fixed = TRUE)
  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(A = ys)),
    "rownames(data) does not carry", fixed = TRUE)
})

test_that("`site.map` keys a source that carries no row names (#241)", {
  d <- data.frame(x = as.numeric(1:5))
  short <- matrix(0L, 3, 2)

  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(short)),
    "carries no site key", fixed = TRUE)
  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(matrix(0L, 7, 2))),
    "carries no site key", fixed = TRUE)

  by_index <- tulpaObs:::.tobs_build_integrated(
    ~ x, ~ 1, d, list(short), site.map = list(c(4L, 2L, 5L)))
  expect_identical(by_index$site_maps[[1]], c(4L, 2L, 5L) - 1L)

  rownames(d) <- paste0("s", 1:5)
  by_name <- tulpaObs:::.tobs_build_integrated(
    ~ x, ~ 1, d, list(short), site.map = list(c("s4", "s2", "s5")))
  expect_identical(by_name$site_maps[[1]], c(4L, 2L, 5L) - 1L)

  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(short),
                                      site.map = list(c(1L, 2L, 99L))),
    "must index sites in 1..5", fixed = TRUE)
  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(short),
                                      site.map = list(c(1L, 1L, 2L))),
    "more than one source row", fixed = TRUE)
  expect_error(
    tulpaObs:::.tobs_build_integrated(~ x, ~ 1, d, list(short, short),
                                      site.map = list(c(1L, 2L, 3L))),
    "list of 2 entries", fixed = TRUE)
})

test_that("a multi-source fit carries the simulator's site maps (#241)", {
  skip_on_cran()
  skip_if_fast()

  sim <- simulate_int_occu(N_total = 120, n_data = 2, J = c(4, 3),
                           n_shared = 40, seed = 7)
  fit <- tobs(~ x, data = sim$data, family = int_occu(), detection = ~ 1,
              y = sim$y, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))

  expect_identical(fit$model$site_maps,
                   lapply(sim$site_maps, function(m) as.integer(m) - 1L))
})

test_that("a keyed compact source fits the NA-padded full source (#241)", {
  skip_on_cran()
  skip_if_fast()

  # The same data written two ways: each source as the rows it surveyed, keyed
  # by name, and each source as a full n_sites grid with NA at the sites it did
  # not survey. Both name the same site for every observation, so the fits
  # agree; under a positional read the compact form is a different model.
  sim <- simulate_int_occu(N_total = 120, n_data = 2, J = c(4, 3),
                           n_shared = 40, seed = 3)
  padded <- lapply(seq_along(sim$y), function(s) {
    full <- matrix(NA_integer_, nrow(sim$data), ncol(sim$y[[s]]))
    full[sim$site_maps[[s]], ] <- sim$y[[s]]
    full
  })

  ctl <- list(verbose = FALSE, progress = FALSE)
  fit_keyed <- tobs(~ x, data = sim$data, family = int_occu(), detection = ~ 1,
                    y = sim$y, method = "laplace", control = ctl)
  fit_padded <- tobs(~ x, data = sim$data, family = int_occu(), detection = ~ 1,
                     y = padded, method = "laplace", control = ctl)

  expect_equal(unname(fit_keyed$means), unname(fit_padded$means),
               tolerance = 1e-4)
  expect_equal(unname(fit_keyed$sds), unname(fit_padded$sds), tolerance = 1e-4)
})
