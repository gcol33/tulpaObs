# Analytic AGHQ gradient (tulpa:::cpp_aghq_objective_grad) vs central finite
# difference of the objective (tulpa:::cpp_aghq_objective), driving the engine
# with the native community / multispecies N-mixture oracle that lives in
# tulpaObs. The analytic gradient is the Fisher-identity gradient of the true
# marginal; it agrees with the FD of the AGHQ objective as n_quad grows (the
# omitted node-placement terms are O(AGHQ truncation)), so these checks run at
# n_quad = 9 where the two match tightly. Per-parameter-block tolerance; a sign
# flip fails loudly.

.sim_community_nmix <- function(seed, S, R, J, mu_lambda, mu_p, Sig_l, Sig_p) {
  set.seed(seed)
  p_lam <- length(mu_lambda); p_p <- length(mu_p)
  Ll <- t(chol(Sig_l)); Lp <- t(chol(Sig_p))
  X_lambda <- cbind(1, as.numeric(scale(rnorm(R))))[, seq_len(p_lam), drop = FALSE]
  y <- integer(0); site_idx <- integer(0); species_idx <- integer(0)
  xp_rows <- list()
  for (s in seq_len(S)) {
    cf_l <- mu_lambda + as.numeric(Ll %*% rnorm(p_lam))
    cf_p <- mu_p      + as.numeric(Lp %*% rnorm(p_p))
    for (i in seq_len(R)) {
      lam <- exp(sum(X_lambda[i, ] * cf_l))
      N   <- rpois(1, lam)
      for (j in seq_len(J)) {
        xp  <- c(1, rnorm(1))[seq_len(p_p)]
        pij <- stats::plogis(sum(xp * cf_p))
        y           <- c(y, stats::rbinom(1, N, pij))
        site_idx    <- c(site_idx, i)
        species_idx <- c(species_idx, s)
        xp_rows[[length(xp_rows) + 1L]] <- xp
      }
    }
  }
  list(y = y, site_idx = site_idx, species_idx = species_idx,
       X_lambda = X_lambda, X_p = do.call(rbind, xp_rows), R = R, S = S)
}

# Central FD of the exposed objective at every coordinate of par.
.aghq_fd_grad <- function(par, orc, nc, full, nq, lkj, h = 1e-4) {
  vapply(seq_along(par), function(jj) {
    pp <- par; pp[jj] <- pp[jj] + h
    pm <- par; pm[jj] <- pm[jj] - h
    (tulpa:::cpp_aghq_objective(pp, orc, nc, full, nq, lkj) -
     tulpa:::cpp_aghq_objective(pm, orc, nc, full, nq, lkj)) / (2 * h)
  }, numeric(1))
}

test_that("analytic AGHQ gradient matches central FD (Poisson community N-mixture)", {
  skip_on_cran()
  skip_if_fast()
  p_lam <- 2L; p_p <- 2L
  d <- .sim_community_nmix(
    seed = 404L, S = 5L, R = 8L, J = 3L,
    mu_lambda = c(1.0, 0.4), mu_p = c(0.2, -0.3),
    Sig_l = matrix(c(0.30, 0.10, 0.10, 0.20), 2, 2),
    Sig_p = matrix(c(0.25, -0.05, -0.05, 0.15), 2, 2))
  K_max <- max(d$y) + 25L

  orc <- cpp_nmix_community_oracle(
    d$y, d$site_idx, d$species_idx, d$X_lambda, d$X_p, d$R, d$S, K_max)
  nc <- c(p_lam, p_p); full <- c(TRUE, TRUE)
  layout <- tulpa:::.re_cov_block_layout(
    list(list(n_groups = d$S, n_coefs = p_lam, correlated = TRUE),
         list(n_groups = d$S, n_coefs = p_p,   correlated = TRUE)), NULL)
  re_par <- tulpa:::.re_cov_L_list_to_theta(
    lapply(list(matrix(c(0.30, 0.10, 0.10, 0.20), 2, 2),
                matrix(c(0.25, -0.05, -0.05, 0.15), 2, 2)), tulpa:::.re_chol_spd),
    layout)
  par0 <- c(c(1.0, 0.4) + 0.10, c(0.2, -0.3) - 0.10, re_par)
  nq <- 9L; lkj <- 1

  ana <- tulpa:::cpp_aghq_objective_grad(par0, orc, nc, full, nq, lkj)$grad
  fd  <- .aghq_fd_grad(par0, orc, nc, full, nq, lkj)

  th_idx  <- seq_len(p_lam + p_p)
  sig_idx <- setdiff(seq_along(par0), th_idx)

  big <- abs(fd) > 1e-3
  expect_true(all(sign(ana[big]) == sign(fd[big])),
              info = paste0("Poisson gradient sign disagrees with FD:\n  analytic = ",
                            paste(signif(ana, 4), collapse = ", "),
                            "\n  FD       = ", paste(signif(fd, 4), collapse = ", ")))

  expect_lt(max(abs(ana[th_idx]  - fd[th_idx])),  1e-4)
  expect_lt(max(abs(ana[sig_idx] - fd[sig_idx])), 1e-4)
})

test_that("analytic AGHQ gradient matches FD with NB dispersion (community, incl log_r)", {
  skip_on_cran()
  skip_if_fast()
  p_lam <- 2L; p_p <- 2L
  d <- .sim_community_nmix(
    seed = 505L, S = 5L, R = 8L, J = 3L,
    mu_lambda = c(1.0, 0.4), mu_p = c(0.2, -0.3),
    Sig_l = matrix(c(0.30, 0.10, 0.10, 0.20), 2, 2),
    Sig_p = matrix(c(0.25, -0.05, -0.05, 0.15), 2, 2))
  K_max <- max(d$y) + 25L

  # nb = TRUE: the dispersion is a per-species random effect log_r_s ~
  # N(mu_log_r, sigma_log_r). The oracle widens the per-species RE vector to
  # d = p_lambda + p_p + 1 (the trailing log_r_s coordinate) and carries mu_log_r
  # as the trailing fixed effect (n_theta = d). The RE-covariance layout therefore
  # has THREE blocks: the abundance and detection coefficient blocks plus a scalar
  # log_r_s block (sigma_log_r), all integrated by the AGHQ engine.
  orc <- cpp_nmix_community_oracle(
    d$y, d$site_idx, d$species_idx, d$X_lambda, d$X_p, d$R, d$S, K_max, nb = TRUE)
  nc <- c(p_lam, p_p, 1L); full <- c(TRUE, TRUE, FALSE)
  layout <- tulpa:::.re_cov_block_layout(
    list(list(n_groups = d$S, n_coefs = p_lam, correlated = TRUE),
         list(n_groups = d$S, n_coefs = p_p,   correlated = TRUE),
         list(n_groups = d$S, n_coefs = 1L,    correlated = FALSE)), NULL)
  re_par <- tulpa:::.re_cov_L_list_to_theta(
    lapply(list(matrix(c(0.30, 0.10, 0.10, 0.20), 2, 2),
                matrix(c(0.25, -0.05, -0.05, 0.15), 2, 2),
                matrix(0.5^2, 1, 1)), tulpa:::.re_chol_spd),
    layout)
  par0 <- c(c(1.0, 0.4) + 0.10, c(0.2, -0.3) - 0.10, log(8), re_par)
  nq <- 9L; lkj <- 1

  ana <- tulpa:::cpp_aghq_objective_grad(par0, orc, nc, full, nq, lkj)$grad
  fd  <- .aghq_fd_grad(par0, orc, nc, full, nq, lkj)

  th_idx  <- seq_len(p_lam + p_p + 1L)            # mu + log_r
  sig_idx <- setdiff(seq_along(par0), th_idx)

  big <- abs(fd) > 1e-3
  expect_true(all(sign(ana[big]) == sign(fd[big])),
              info = paste0("NB gradient sign disagrees with FD (incl log_r):\n  ana=",
                            paste(signif(ana, 4), collapse = ", "),
                            "\n  fd =", paste(signif(fd, 4), collapse = ", ")))
  expect_lt(max(abs(ana[th_idx]  - fd[th_idx])),  5e-4)   # mu + log_r block
  expect_lt(max(abs(ana[sig_idx] - fd[sig_idx])), 5e-4)   # log-Cholesky Sigma block
})
