# Joint Louis observed information for single-season (occupancy, detection)
# coefficients (#369). The Louis identity holds at any parameter value when the
# occupancy weight is the posterior at that value, so the assembled joint info
# must equal minus the Hessian of the closed-form marginal log-likelihood.

test_that("joint Louis info equals the exact marginal Hessian", {
  set.seed(11)
  n <- 150L; J <- 3L
  X_occ <- cbind(1, stats::rnorm(n))
  X_det <- cbind(1, stats::rnorm(n))
  b_psi <- c(0.2, 0.8); b_det <- c(-0.4, 0.6)
  psi <- plogis(drop(X_occ %*% b_psi)); p <- plogis(drop(X_det %*% b_det))
  z <- stats::rbinom(n, 1L, psi)
  y <- matrix(stats::rbinom(n * J, 1L, rep(z * p, J)), n, J)
  y[sample(n * J, 40L)] <- -1L                    # unsurveyed visits
  valid <- y >= 0
  Y <- y; Y[!valid] <- 0L
  nv <- rowSums(valid)
  any_det <- rowSums(Y) > 0

  q <- (1 - p)^nv
  w <- ifelse(any_det, 1, psi * q / (psi * q + 1 - psi))

  I <- .louis_info_joint_single(X_occ, b_psi, X_det, b_det, w, nv)

  nlp <- function(par) .tobs_occu_marginal_nlp(
    par, X_occ, X_det, NULL, valid, Y, any_det, n, J, 2L, 2L, 0L,
    pmean = rep(0, 4), pprec = rep(0, 4))
  par <- c(b_psi, b_det)
  h <- 1e-4
  H <- matrix(0, 4, 4)
  for (a in 1:4) for (b in 1:4) {
    ea <- replace(numeric(4), a, h); eb <- replace(numeric(4), b, h)
    H[a, b] <- (nlp(par + ea + eb) - nlp(par + ea - eb) -
                nlp(par - ea + eb) + nlp(par - ea - eb)) / (4 * h^2)
  }
  expect_equal(unname(I), H, tolerance = 1e-5)
  # The cross block carries real information, so the agreement is not vacuous.
  expect_gt(max(abs(I[1:2, 3:4])), 0.1 * min(abs(diag(I))))
})

test_that("a spatial single-season fit draws the arms as correlated", {
  skip_on_cran()
  skip_if_fast()
  set.seed(5)
  n_cells <- 30L; reps <- 5L; J <- 3L
  adj <- chain_adj(n_cells)
  cell <- rep(seq_len(n_cells), each = reps)
  n_sites <- n_cells * reps
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  f0 <- 0.8 * sin(2 * pi * seq_len(n_cells) / n_cells)
  z <- stats::rbinom(n_sites, 1L, plogis(0.3 + 0.5 * x + f0[cell]))
  y <- matrix(stats::rbinom(n_sites * J, 1L, rep(z * plogis(-0.3), J)),
              n_sites, J)
  fit <- tobs(~ x + icar(graph = adj, group_var = "cell"),
              detection = ~ 1, data = data.frame(x = x, cell = cell),
              family = occu(), y = y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))

  C <- stats::cor(fit$draws[, c("psi_(Intercept)", "p_(Intercept)")])
  # With J = 3 the occupancy and detection intercepts trade off: the joint
  # Louis covariance puts a clear negative correlation between them.
  expect_lt(C[1, 2], -0.1)
  expect_equal(unname(fit$sds[c("psi_(Intercept)", "p_(Intercept)")]),
               unname(apply(fit$draws[, c("psi_(Intercept)", "p_(Intercept)")],
                            2L, stats::sd)),
               tolerance = 0.1)
})
