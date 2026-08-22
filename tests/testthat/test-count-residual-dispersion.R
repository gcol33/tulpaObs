# =============================================================================
# test-count-residual-dispersion.R -- residuals() scores a count fit at the
# variance its own mixture implies.
#
# abun(), removal() and distance() all observe a thinned latent abundance, and
# the negative binomial is closed under binomial thinning with the SAME size, so
# a per-cell marginal of an NB(r, lambda) abundance is NB(r, lambda * pi) with
# variance mu + mu^2/r. Those three share .tobs_count_residual() with count().
#
# dyn_abun() does NOT: its estimated size describes the INITIAL abundance, and
# N_t for t > 1 is a survival-thinning plus recruitment convolution that is not
# negative binomial. Its Pearson residual uses the exact moment recursion
# instead, which is asserted here against Monte Carlo.
#
# Every family also asserts the Poisson path is unchanged to the bit, so the fix
# cannot have moved a correctly specified Poisson fit.
# =============================================================================

ctl <- list(verbose = FALSE, progress = FALSE)

# The formula all four handlers used before the fix.
old_pois <- function(y, mu, type) switch(type,
  response = y - mu,
  pearson  = (y - mu) / sqrt(mu),
  deviance = { d <- 2 * (ifelse(y > 0, y * log(y / mu), 0) - (y - mu))
               sign(y - mu) * sqrt(pmax(d, 0)) })

# Written out rather than called from the package, so the assertion is against
# the negative-binomial definition and not against the code under test.
nb_pearson  <- function(y, mu, r) (y - mu) / sqrt(mu + mu^2 / r)
nb_deviance <- function(y, mu, r) {
  term <- ifelse(y > 0, y * log(y / mu), 0)
  d <- 2 * (term - (y + r) * log((y + r) / (mu + r)))
  sign(y - mu) * sqrt(pmax(d, 0))
}

test_that("the Dail-Madsen marginal variance recursion matches Monte Carlo", {
  skip_on_cran()
  # Var[N_t] = omega^2 Var[N_{t-1}] + omega (1 - omega) E[N_{t-1}] + gamma, and
  # a count is Binomial(N_t, p) -> Var[y] = p^2 Var[N_t] + p (1 - p) E[N_t].
  # Under a Poisson initial abundance Var[N_t] == E[N_t] at every season, so the
  # whole expression collapses to mu -- which is why the Poisson path below is
  # byte-identical rather than merely close.
  set.seed(1)
  lam <- 8; om <- 0.7; gam <- 1.2; pdet <- 0.6; r <- 4; Tn <- 4L; M <- 150000L
  for (mix in c("poisson", "negbin")) {
    EN <- VN <- numeric(Tn)
    EN[1] <- lam
    VN[1] <- if (mix == "negbin") lam + lam^2 / r else lam
    for (t in seq_len(Tn - 1L)) {
      EN[t + 1] <- om * EN[t] + gam
      VN[t + 1] <- om^2 * VN[t] + om * (1 - om) * EN[t] + gam
    }
    if (identical(mix, "poisson")) expect_equal(VN, EN)   # the collapse

    cur <- if (mix == "negbin") stats::rnbinom(M, size = r, mu = lam)
           else stats::rpois(M, lam)
    emp <- numeric(Tn)
    for (t in seq_len(Tn)) {
      if (t > 1L) cur <- stats::rbinom(M, cur, om) + stats::rpois(M, gam)
      emp[t] <- stats::var(stats::rbinom(M, cur, pdet))
    }
    pred <- pdet^2 * VN + pdet * (1 - pdet) * EN
    expect_lt(max(abs(emp - pred) / pred), 0.03)   # measured <= 0.009 at M = 2e5
  }
})

test_that("abun() residuals: Poisson unchanged, negbin at its own variance", {
  skip_on_cran()
  args <- list(N = 90, J = 4, n_abund_covs = 1, n_det_covs = 1,
               beta_lambda = c(log(6), 0.4), beta_p = c(0.3, -0.2), seed = 3)
  cell_mu <- function(fit) {
    fv <- .tobs_fitted_nmix(fit)
    pmax(fv$lambda[fit$model$site_idx] * fv$p, 1e-10)
  }
  flat <- function(fit, ty) residuals(fit, type = ty)[
    cbind(fit$model$site_idx, fit$model$visit_idx)]

  sim <- do.call(simulate_abun, args)
  fp  <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
              family = abun(), method = "laplace", control = ctl)
  mu  <- cell_mu(fp); y <- as.numeric(fp$model$y_long)
  for (ty in c("response", "pearson", "deviance"))
    expect_identical(flat(fp, ty), old_pois(y, mu, ty))

  simn <- do.call(simulate_abun, c(args, list(mixture = "negbin", size = 3)))
  fn   <- tobs(~ abund_cov1, data = simn$data, detection = ~ det_cov1,
               y = simn$y, family = abun(mixture = "negbin"),
               method = "laplace", control = ctl)
  mun <- cell_mu(fn); yn <- as.numeric(fn$model$y_long)
  r   <- fn$nmix_dispersion$r
  expect_true(is.finite(r))
  expect_equal(flat(fn, "pearson"),  nb_pearson(yn, mun, r))
  expect_equal(flat(fn, "deviance"), nb_deviance(yn, mun, r))
  # Discriminates: the Poisson form is materially different at this dispersion.
  expect_gt(max(abs(flat(fn, "pearson") - old_pois(yn, mun, "pearson"))), 0.05)
})

test_that("removal() residuals: Poisson unchanged, negbin at its own variance", {
  skip_on_cran()
  cell_mu <- function(fit) {
    fv <- .tobs_fitted_nmix(fit); m <- fit$model
    pm <- matrix(NA_real_, m$n_sites, m$n_passes)
    pm[cbind(m$site_idx, m$visit_idx)] <- fv$p
    mu <- fv$lambda * t(apply(pm, 1, .removal_pi))
    pmax(mu[cbind(m$site_idx, m$visit_idx)], 1e-10)
  }
  flat <- function(fit, ty) residuals(fit, type = ty)[
    cbind(fit$model$site_idx, fit$model$visit_idx)]

  sim <- simulate_removal(N = 90, K = 3, n_abund_covs = 1, n_det_covs = 1,
                          beta_lambda = c(log(7), 0.3), beta_p = c(0.4, -0.2),
                          seed = 4)
  fp <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
             family = removal(), method = "laplace", control = ctl)
  mu <- cell_mu(fp); y <- as.numeric(fp$model$y_long)
  for (ty in c("response", "pearson", "deviance"))
    expect_identical(flat(fp, ty), old_pois(y, mu, ty))

  simn <- simulate_removal(N = 90, K = 3, n_abund_covs = 1, n_det_covs = 1,
                           beta_lambda = c(log(7), 0.3), beta_p = c(0.4, -0.2),
                           mixture = "negbin", size = 3, seed = 4)
  fn <- tobs(~ abund_cov1, data = simn$data, detection = ~ det_cov1, y = simn$y,
             family = removal(mixture = "negbin"), method = "laplace",
             control = ctl)
  mun <- cell_mu(fn); yn <- as.numeric(fn$model$y_long)
  r   <- fn$nmix_dispersion$r
  expect_true(is.finite(r))
  expect_equal(flat(fn, "pearson"), nb_pearson(yn, mun, r))
})

test_that("dyn_abun() residuals: Poisson unchanged, negbin at the recursion", {
  skip_on_cran()
  args <- list(N = 50, T = 3, J = 3, n_abund_covs = 1,
               beta_lambda = c(log(8), 0), p = 0.6, omega = 0.7, gamma = 1.2,
               seed = 4)
  exact_var <- function(fit) {
    fv <- fit$model; f <- .tobs_fitted_dyn_abun(fit)
    Tn <- fv$n_seasons
    is_nb <- (fit$mixture %||% "poisson") %in% c("negbin", "zinb")
    rr <- fit$dispersion$r %||% NA_real_
    VN <- matrix(0, fv$n_sites, Tn)
    VN[, 1] <- if (is_nb && is.finite(rr)) f$lambda + f$lambda^2 / rr else f$lambda
    for (t in seq_len(Tn - 1L))
      VN[, t + 1L] <- f$omega[, t]^2 * VN[, t] +
        f$omega[, t] * (1 - f$omega[, t]) * f$EN[, t] + f$gamma[, t]
    f$p^2 * VN + f$p * (1 - f$p) * f$EN
  }

  sim <- do.call(simulate_dyn_abun, args)
  fp  <- tobs(~ 1, data = sim$data, detection = ~ 1, y = sim$y,
              family = dyn_abun(), method = "laplace", control = ctl)
  f   <- .tobs_fitted_dyn_abun(fp)
  mu_t <- f$EN * f$p
  for (ty in c("response", "pearson", "deviance")) {
    got <- residuals(fp, type = ty)
    for (t in seq_len(fp$model$n_seasons)) {
      mu <- pmax(matrix(mu_t[, t], fp$model$n_sites, fp$model$max_visits), 1e-10)
      expect_identical(got[, , t], old_pois(sim$y[, , t], mu, ty))
    }
  }

  simn <- do.call(simulate_dyn_abun, c(args, list(mixture = "negbin", r = 3)))
  fn <- tobs(~ 1, data = simn$data, detection = ~ 1, y = simn$y,
             family = dyn_abun(mixture = "negbin"), method = "laplace",
             control = ctl)
  vt <- exact_var(fn); fvn <- .tobs_fitted_dyn_abun(fn)
  mun_t <- fvn$EN * fvn$p
  got <- residuals(fn, type = "pearson")
  nsi <- fn$model$n_sites; nvi <- fn$model$max_visits
  for (t in seq_len(fn$model$n_seasons)) {
    mu <- pmax(matrix(mun_t[, t], nsi, nvi), 1e-10)
    sd <- sqrt(pmax(matrix(vt[, t], nsi, nvi), 1e-10))
    expect_equal(got[, , t], (simn$y[, , t] - mu) / sd)
  }
  # Discriminates: season 1 is exactly NB, so its variance exceeds the mean.
  expect_gt(max(vt[, 1] / mun_t[, 1]), 1.05)
})

test_that("count() residuals are unchanged by the shared helper", {
  skip_on_cran()
  for (resp in c("poisson", "negbin", "gaussian")) {
    sim <- simulate_count(N = 150, beta = c(0.8, 0.5), response = resp, seed = 5)
    fit <- tobs(~ x, data = sim$data, y = sim$y, family = count(response = resp),
                method = "laplace", control = ctl)
    y <- as.numeric(fit$model$y_count); mu <- fitted(fit)$mu
    grab <- function(ty) { r <- residuals(fit, type = ty)
                           if (is.list(r)) r$mu else r }
    expect_identical(grab("response"), y - mu)
    if (identical(resp, "poisson")) {
      mup <- pmax(mu, 1e-8)
      expect_identical(grab("pearson"),  old_pois(y, mup, "pearson"))
      expect_identical(grab("deviance"),
                       sign(y - mu) * sqrt(pmax(2 * (ifelse(y > 0,
                         y * log(y / mup), 0) - (y - mu)), 0)))
    } else if (identical(resp, "negbin")) {
      mup <- pmax(mu, 1e-8); size <- fit$count_dispersion$phi %||% Inf
      expect_equal(grab("pearson"),  nb_pearson(y, mup, size))
      expect_equal(grab("deviance"), nb_deviance(y, mup, size))
    } else {
      phi <- fit$count_dispersion$phi %||% 1
      expect_identical(grab("pearson"),
                       (y - mu) / sqrt(rep(max(phi, 1e-8), length(mu))))
      expect_identical(grab("deviance"), y - mu)
    }
  }
})
