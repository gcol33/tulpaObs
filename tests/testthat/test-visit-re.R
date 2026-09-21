# Random effects over visit rows on the occu() detection logit (#363): an
# observer who changes between visits to the same site. The Laplace route fits
# them in the random-effect EM on per-visit detection rows, the NUTS route as
# extra parameters of the occupancy likelihood with the engine's own term prior.

visit_re_data <- function(seed, N = 250L, J = 4L, n_groups = 10L, sigma = 1) {
  sim <- simulate_occu(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.3, 0.8), beta_det = c(0, 0.6),
                       n_visit_groups = n_groups, sigma_visit = sigma,
                       seed = seed)
  sim$v <- data.frame(det_cov1 = rep(sim$data$det_cov1, each = J),
                      visit_group = sim$visits$visit_group)
  sim
}

test_that("simulate_occu() draws a per-visit group effect only when asked", {
  plain <- simulate_occu(N = 30, J = 3, seed = 1)
  expect_null(plain$visits)
  expect_null(plain$truth$b_visit)

  sim <- simulate_occu(N = 30, J = 3, n_visit_groups = 4, sigma_visit = 0.5,
                       seed = 1)
  expect_equal(nrow(sim$visits), 90L)
  expect_identical(levels(sim$visits$visit_group), as.character(1:4))
  expect_length(sim$truth$b_visit, 4L)
  expect_identical(sim$truth$sigma_visit, 0.5)
})

test_that("a bar in the visit formula becomes a visit-level random effect", {
  N <- 8L; J <- 3L
  y <- matrix(rbinom(N * J, 1, 0.4), N, J)
  d <- data.frame(x = rnorm(N))
  v <- data.frame(effort = runif(N * J),
                  obs = factor(rep(c("a", "b", "c", "d"), length.out = N * J)))
  # `- 1` is how tobs() hands `detection` to the visit path.
  m <- tulpaObs:::.tobs_build_single(~ x, ~ 1, d, y,
                                     ~ effort + (1 | obs) - 1, v)
  expect_identical(colnames(m$X_det_visit), "effort")
  re <- Filter(function(t) inherits(t$spec, "tobs_re"), m$structured_terms)
  expect_length(re, 1L)
  expect_identical(re[[1]]$spec$level, "visit")
  expect_identical(re[[1]]$processes, 2L)
  expect_length(re[[1]]$spec$group_idx, N * J)
  expect_identical(m$det_visit_data, v)

  # A bar alone leaves no visit fixed effect.
  m2 <- tulpaObs:::.tobs_build_single(~ x, ~ 1, d, y, ~ (1 | obs) - 1, v)
  expect_null(m2$X_det_visit)
})

test_that("only random effects can be written over visit rows", {
  N <- 8L; J <- 3L
  y <- matrix(rbinom(N * J, 1, 0.4), N, J)
  d <- data.frame(x = rnorm(N))
  v <- data.frame(t = rep(1:J, times = N))
  expect_error(
    tulpaObs:::.tobs_build_single(~ x, ~ 1, d, y, ~ temporal(t) - 1, v),
    "only random effects can be written over visit rows")
})

test_that("a visit-level random effect is refused off the laplace / nuts routes", {
  sim <- visit_re_data(1, N = 40L)
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = occu(),
         detection = ~ det_cov1 + (1 | visit_group), y = sim$y, visits = sim$v,
         method = "nested_laplace",
         control = list(verbose = FALSE, progress = FALSE)),
    "fitted with method = \"laplace\" or \"nuts\"", fixed = TRUE)
})

test_that("a correction on a random-effect Laplace fit is refused (#366)", {
  sim <- simulate_occu(N = 60, J = 3, seed = 1)
  d <- cbind(sim$data, g = factor(rep(1:6, 10)))
  expect_error(
    tobs(~ occ_cov1 + (1 | g), data = d, family = occu(),
         detection = ~ det_cov1, y = sim$y, method = "laplace_gibbs",
         control = list(verbose = FALSE, progress = FALSE)),
    "has no MI/Gibbs correction", fixed = TRUE)
})

test_that("Laplace recovers a visit-level random intercept", {
  skip_on_cran()
  sim <- visit_re_data(1)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1 + (1 | visit_group), y = sim$y,
              visits = sim$v, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  m <- fit$means
  expect_true(all(c("p_visit_det_cov1", "sigma_p1_(Intercept)") %in% names(m)))
  expect_lt(abs(m[["sigma_p1_(Intercept)"]] - sd(sim$truth$b_visit)), 0.35)
  expect_gt(cor(fit$re_effects$p1$estimate, sim$truth$b_visit), 0.85)
  expect_lt(abs(m[["psi_occ_cov1"]] - 0.8), 0.4)
  expect_lt(abs(m[["p_visit_det_cov1"]] - 0.6), 0.4)
  # The per-group quadrature does not factor over visit rows.
  expect_identical(fit$aghq$declined, "visit_rows")
})

test_that("NUTS recovers a visit-level random intercept and names it in place", {
  skip_on_cran()
  skip_if_fast()
  sim <- visit_re_data(1)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1 + (1 | visit_group), y = sim$y,
              visits = sim$v, method = "nuts",
              control = list(n.iter = 500, n.warmup = 500, n.chains = 2,
                             seed = 1, verbose = FALSE, progress = FALSE))
  m <- fit$means
  sigma_hat <- exp(m[["log_sigma_v1_(Intercept)"]])
  expect_lt(abs(sigma_hat - sd(sim$truth$b_visit)), 0.35)
  expect_gt(cor(fit$re_effects$v1$estimate, sim$truth$b_visit), 0.85)
  expect_lt(abs(m[["p_visit_det_cov1"]] - 0.6), 0.4)
  expect_length(fit$re, 1L)
  expect_identical(fit$re[[1]]$level, "visit")
})

test_that("NUTS names a visit coefficient in place beside a site random effect (#364)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(2)
  N <- 200L; J <- 4L
  sim <- simulate_occu(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.5, 1), beta_det = c(0, 0), seed = 2)
  d <- cbind(sim$data, g = factor(rep(1:10, length.out = N)))
  vx <- rnorm(N * J)
  p <- plogis(-0.2 + 2.5 * vx)
  y <- matrix(rbinom(N * J, 1, rep(sim$truth$z, each = J) * p), N, J,
              byrow = TRUE)
  v <- data.frame(vx = vx)
  attr(v, "formula") <- ~ vx
  fit <- tobs(~ occ_cov1 + (1 | g), data = d, family = occu(), detection = ~ 1,
              y = y, visits = v, method = "nuts",
              control = list(n.iter = 300, n.warmup = 300, n.chains = 1,
                             seed = 1, verbose = FALSE, progress = FALSE))
  expect_lt(abs(fit$means[["p_visit_vx"]] - 2.5), 0.8)
  expect_identical(tail(names(fit$means), 1L), "p_visit_vx")
})

test_that("Laplace and NUTS agree on a visit-level random effect", {
  skip_on_cran()
  skip_if_fast()
  sim <- visit_re_data(3)
  ctl <- list(verbose = FALSE, progress = FALSE)
  la <- tobs(~ occ_cov1, data = sim$data, family = occu(),
             detection = ~ det_cov1 + (1 | visit_group), y = sim$y,
             visits = sim$v, method = "laplace", control = ctl)
  nu <- tobs(~ occ_cov1, data = sim$data, family = occu(),
             detection = ~ det_cov1 + (1 | visit_group), y = sim$y,
             visits = sim$v, method = "nuts",
             control = c(ctl, list(n.iter = 500, n.warmup = 500, n.chains = 2,
                                   seed = 1)))
  for (nm in c("psi_occ_cov1", "p_visit_det_cov1")) {
    expect_lt(abs(la$means[[nm]] - nu$means[[nm]]), 0.2, label = nm)
  }
  expect_gt(cor(la$re_effects$p1$estimate, nu$re_effects$v1$estimate), 0.95)
})

test_that("the visit-level SD is recovered without a one-sided bias over seeds", {
  skip_on_cran()
  skip_if_fast()
  err <- vapply(1:8, function(seed) {
    sim <- visit_re_data(seed)
    fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
                detection = ~ det_cov1 + (1 | visit_group), y = sim$y,
                visits = sim$v, method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    fit$means[["sigma_p1_(Intercept)"]] - sd(sim$truth$b_visit)
  }, numeric(1))
  expect_lt(abs(mean(err)), 0.15)
})
