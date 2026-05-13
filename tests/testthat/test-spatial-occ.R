test_that("ICAR spatial occupancy runs", {
  set.seed(42)
  n_sites <- 20
  adj <- matrix(0, n_sites, n_sites)
  for (i in 1:(n_sites-1)) { adj[i, i+1] <- 1; adj[i+1, i] <- 1 }

  z <- rbinom(n_sites, 1, 0.5)
  y <- matrix(0L, n_sites, 3)
  for (i in seq_len(n_sites)) if (z[i] == 1) y[i, ] <- rbinom(3, 1, 0.5)

  sp <- tobs_icar(adj)
  expect_s3_class(sp, "tobs_spatial")

  fit <- tobs(~ 1, data.frame(x = rnorm(n_sites)), family = occu(),
              detection = ~ 1, y = y, spatial = sp,
              engine = "nuts",
              control = list(iter = 100, warmup = 50, seed = 1, verbose = FALSE))
  expect_true(fit$n_params > 2)
})

test_that("GP spatial occupancy runs", {
  set.seed(42)
  n_sites <- 20
  coords <- cbind(runif(n_sites), runif(n_sites))
  z <- rbinom(n_sites, 1, 0.5)
  y <- matrix(0L, n_sites, 3)
  for (i in seq_len(n_sites)) if (z[i] == 1) y[i, ] <- rbinom(3, 1, 0.5)

  sp <- tobs_gp(coords, nn = 5)
  expect_s3_class(sp, "tobs_spatial")

  fit <- tobs(~ 1, data.frame(x = rnorm(n_sites)), family = occu(),
              detection = ~ 1, y = y, spatial = sp,
              engine = "nuts",
              control = list(iter = 100, warmup = 50, seed = 1, verbose = FALSE))
  expect_true(fit$n_params > 2)
})

test_that("tobs_icar validates inputs", {
  expect_error(tobs_icar(matrix(1:4, 2)), "symmetric")
})
