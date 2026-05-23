test_that("ICAR spatial occupancy runs", {
  set.seed(42)
  n_sites <- 20
  adj <- matrix(0, n_sites, n_sites)
  for (i in 1:(n_sites-1)) { adj[i, i+1] <- 1; adj[i+1, i] <- 1 }

  z <- rbinom(n_sites, 1, 0.5)
  y <- matrix(0L, n_sites, 3)
  for (i in seq_len(n_sites)) if (z[i] == 1) y[i, ] <- rbinom(3, 1, 0.5)

  fit <- tobs(~ icar(graph = adj), data.frame(x = rnorm(n_sites)),
              family = occu(), detection = ~ 1, y = y,
              engine = "nuts",
              control = list(iter = 100, warmup = 50, seed = 1, verbose = FALSE))
  expect_true(fit$n_params > 2)
})

test_that("GP spatial occupancy runs", {
  set.seed(42)
  n_sites <- 20
  dat <- data.frame(lon = runif(n_sites), lat = runif(n_sites))
  z <- rbinom(n_sites, 1, 0.5)
  y <- matrix(0L, n_sites, 3)
  for (i in seq_len(n_sites)) if (z[i] == 1) y[i, ] <- rbinom(3, 1, 0.5)

  fit <- tobs(~ gp(lon, lat, nn = 5), dat, family = occu(),
              detection = ~ 1, y = y,
              engine = "nuts",
              control = list(iter = 100, warmup = 50, seed = 1, verbose = FALSE))
  expect_true(fit$n_params > 2)
})

test_that("icar() validates its graph", {
  expect_error(tulpaObs:::.tobs_term_icar(matrix(1:4, 2)), "symmetric")
})
