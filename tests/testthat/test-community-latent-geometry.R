.chain_adj <- function(n) {
  A <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) A[i, i + 1L] <- A[i + 1L, i] <- 1
  A
}

test_that("community latent icar field is centred per connected component (#424)", {
  # Two islands of 5 nodes each, one node per site.
  A <- matrix(0, 10, 10)
  A[1:5, 1:5] <- .chain_adj(5); A[6:10, 6:10] <- .chain_adj(5)
  term <- tulpaObs:::.tobs_term_icar(graph = A)
  model <- list(n_sites = 10L, data = data.frame(x = seq_len(10)))
  geom <- tulpaObs:::.tobs_latent_field_setup(term, model, what = "test")
  expect_length(geom$components, 2L)
  expect_setequal(lapply(geom$components, sort), list(1:5, 6:10))

  # A proper CAR keeps the single whole-graph constraint.
  geom_car <- tulpaObs:::.tobs_latent_field_setup(
    tulpaObs:::.tobs_term_car_proper(graph = A), model, what = "test",
    allow = c("icar", "car_proper"))
  expect_length(geom_car$components, 1L)

  # Poisson working oracle: each island carries its own within-island gradient.
  set.seed(1)
  S <- 4L
  f_true <- c(seq(-1, 1, length.out = 5), seq(1, -1, length.out = 5))
  y <- matrix(stats::rpois(10 * S, exp(1 + f_true)), 10, S)
  oracle <- list(working = function(eta) {
    mu <- exp(eta); list(score = y - mu, curv = mu)
  }, data_ll = function(eta) sum(stats::dpois(y, exp(eta), log = TRUE)))
  Q <- Matrix::Diagonal(x = rowSums(A)) - methods::as(A, "CsparseMatrix")
  fu <- tulpaObs:::.tobs_latent_field_solve(
    oracle, matrix(1, 10, S), Q, matrix(0, 10, 1), 1, matrix(1, 10, 1),
    components = geom$components)
  expect_equal(sum(fu$Fm[1:5, 1]), 0, tolerance = 1e-10)
  expect_equal(sum(fu$Fm[6:10, 1]), 0, tolerance = 1e-10)
  expect_gt(stats::cor(fu$Fm[, 1], f_true), 0.9)

  # The BYM2 structured part is centred per component as well.
  br <- tulpaObs:::.tobs_latent_bym2_solve(
    oracle, matrix(1, 10, S), Q, numeric(10), numeric(10), a = 0.8, b = 0.3,
    components = geom$components)
  expect_equal(sum(br$v[1:5]), 0, tolerance = 1e-10)
  expect_equal(sum(br$v[6:10]), 0, tolerance = 1e-10)
})

test_that("community latent SPDE precision follows spde(nu =) (#425)", {
  set.seed(2)
  co <- cbind(stats::runif(25), stats::runif(25))
  spec1 <- tulpaObs:::.tobs_term_spde(coords = co, nu = 1)$tulpa_spec
  kappa <- 3
  C0 <- Matrix::Diagonal(x = spec1$C0_diag)
  G  <- spec1$G
  Q_closed <- kappa^4 * C0 + 2 * kappa^2 * G +
    G %*% Matrix::Diagonal(x = 1 / spec1$C0_diag) %*% G
  expect_equal(as.matrix(tulpaObs:::.tobs_latent_spde_Q(spec1, kappa)),
               as.matrix(Q_closed), tolerance = 1e-10)

  spec2 <- tulpaObs:::.tobs_term_spde(coords = co, nu = 2)$tulpa_spec
  Q2 <- tulpaObs:::.tobs_latent_spde_Q(spec2, kappa)
  expect_false(isTRUE(all.equal(as.matrix(Q2), as.matrix(Q_closed))))

  # A fractional alpha has no precision of this form and is rejected.
  spec_h <- tulpaObs:::.tobs_term_spde(coords = co, nu = 0.5)$tulpa_spec
  model <- list(n_sites = 25L)
  expect_error(
    tulpaObs:::.tobs_latent_spde_setup(list(list(tulpa_spec = spec_h)), model,
                                       "ms_count()"),
    "nu must be a whole number")
  expect_silent(
    tulpaObs:::.tobs_latent_spde_setup(list(list(tulpa_spec = spec2)), model,
                                       "ms_count()"))
})
