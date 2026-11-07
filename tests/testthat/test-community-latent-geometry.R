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

  # The geometry handed to tulpa_latent_ascent() carries those constraints and
  # the intrinsic precision D - W; the proper CAR's precision moves with rho.
  g_icar <- tulpaObs:::.tobs_latent_geometry(geom)
  expect_identical(g_icar$components, geom$components)
  L <- diag(rowSums(A)) - A
  expect_equal(as.matrix(g_icar$Q(numeric(0))), L, ignore_attr = TRUE)
  expect_length(g_icar$hyper_grid, 0L)
  g_car <- tulpaObs:::.tobs_latent_geometry(geom_car)
  expect_identical(g_car$components, geom_car$components)
  expect_equal(as.matrix(g_car$Q(c(rho = 0.5))), diag(rowSums(A)) - 0.5 * A,
               ignore_attr = TRUE)
  expect_named(g_car$hyper_grid, "rho")
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
  model <- list(n_sites = 25L)
  spde_Q <- function(spec, kappa) {
    geom <- tulpaObs:::.tobs_latent_spde_setup(list(list(tulpa_spec = spec)),
                                               model, "ms_count()")
    g <- tulpaObs:::.tobs_latent_geometry(geom)
    expect_null(g$components)
    g$Q(c(kappa = kappa))
  }
  expect_equal(as.matrix(spde_Q(spec1, kappa)), as.matrix(Q_closed),
               tolerance = 1e-10)

  spec2 <- tulpaObs:::.tobs_term_spde(coords = co, nu = 2)$tulpa_spec
  Q2 <- spde_Q(spec2, kappa)
  expect_false(isTRUE(all.equal(as.matrix(Q2), as.matrix(Q_closed))))

  # A fractional alpha has no precision of this form and is rejected.
  spec_h <- tulpaObs:::.tobs_term_spde(coords = co, nu = 0.5)$tulpa_spec
  expect_error(
    tulpaObs:::.tobs_latent_spde_setup(list(list(tulpa_spec = spec_h)), model,
                                       "ms_count()"),
    "nu must be a whole number")
  expect_silent(
    tulpaObs:::.tobs_latent_spde_setup(list(list(tulpa_spec = spec2)), model,
                                       "ms_count()"))
})
