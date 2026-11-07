# Mode-centred CCD for the outer field-hyperparameter integration of the
# areal-BFGS families (tulpa_hyper_grid(design = "ccd"); the design itself is
# tested in tulpa): an end-to-end engage/decline check on an areal family.

test_that("areal fit: default grid and opt-in CCD both recover the field + slope", {
  skip_if_fast()
  skip_on_cran()
  .adj <- function(side) {
    ng <- side*side; co <- expand.grid(x = seq_len(side), y = seq_len(side))
    a <- matrix(0L, ng, ng)
    for (i in seq_len(ng)) for (j in seq_len(ng))
      if (i != j && abs(co$x[i]-co$x[j]) + abs(co$y[i]-co$y[j]) == 1L) a[i,j] <- 1L
    a
  }
  .sim <- function(adj, cuts, seed, sd_phi = 0.6) {
    set.seed(seed); ng <- nrow(adj); nb <- length(cuts) - 1L
    phi <- as.numeric(scale(rnorm(ng)))
    for (rep in 1:3) { pn <- phi; for (i in seq_len(ng)) { nb_i <- which(adj[i,]==1L); pn[i] <- 0.5*phi[i]+0.5*mean(phi[nb_i]) }; phi <- pn }
    phi <- sd_phi*as.numeric(scale(phi)); phi <- phi-mean(phi)
    x_ab <- rnorm(ng); x_sg <- rnorm(ng)
    lam <- exp(log(50) + 0.5*x_ab + phi); sg <- exp(log(0.4) + 0.2*x_sg)
    y <- matrix(0L, ng, nb)
    for (i in seq_len(ng)) {
      pib <- tulpaObs:::.distance_pi(sg[i], cuts, "halfnorm", "line")
      probs <- c(pib, max(1-sum(pib),0)); N <- rpois(1, lam[i])
      if (N>0) { cc <- rmultinom(1,N,probs); y[i,] <- cc[seq_len(nb)] }
    }
    list(y = y, data = data.frame(abund_cov1 = x_ab, sigma_cov1 = x_sg), phi = phi)
  }
  cuts <- seq(0, 1, length.out = 6); adj <- .adj(7L)
  fam  <- distance(key = "halfnorm", transect = "line", cutpoints = cuts)
  sim  <- .sim(adj, cuts, seed = 402)
  fm   <- ~ abund_cov1 + car_proper(graph = adj)

  fg <- tobs(formula = fm, data = sim$data, family = fam, detection = ~ sigma_cov1,
             y = sim$y, method = "nested_laplace",
             control = list(progress = FALSE, verbose = FALSE))           # default grid
  fc <- tobs(formula = fm, data = sim$data, family = fam, detection = ~ sigma_cov1,
             y = sim$y, method = "nested_laplace",
             control = list(progress = FALSE, verbose = FALSE, integration = "ccd"))

  # Default path is the fixed grid; the CCD path either engages or declines to it.
  expect_identical(fg$spatial_integration, "grid")
  expect_true(fc$spatial_integration %in% c("ccd", "grid"))
  # Both recover the abundance slope (truth 0.5) and the field.
  for (f in list(fg, fc)) {
    expect_lt(abs(f$means[["lambda_abund_cov1"]] - 0.5), 0.25)
    expect_gt(cor(f$spatial_field, sim$phi), 0.5)
  }
  # CCD must not move the estimate materially versus the grid.
  expect_lt(abs(fg$means[["lambda_abund_cov1"]] - fc$means[["lambda_abund_cov1"]]), 0.1)
})
